/* pileup_native: single-pass htslib pileup counter, replicating the exact
 * read-admission and base/quality filters of pileup_counts.PileupCountsProvider
 * ._count_locus / .window_counts (Python), but as one streaming htslib pass
 * over an arbitrary region instead of one pysam pileup() call per 64bp window
 * plus a Python-level per-read loop.
 *
 * Equivalence contract (must match pileup_counts.py exactly):
 *   - read admission: mapq >= min_mapq; not (dup|qcfail|unmapped); not
 *     (secondary|supplementary) -- checked BEFORE depth/mapping_sum increment.
 *   - depth, mapping_sum count every admitted read, including deletions/gaps.
 *   - is_del / is_refskip -> gap_reads++, no base/quality counted.
 *   - base quality (0 if BAM has no per-base qualities, matching pysam's
 *     query_qualities None-check on byte 0 == 0xff) must be >= min_baseq to
 *     count toward count_A/C/G/T and quality_sum.
 *   - only A/C/G/T bases counted (matches config.NUCLEOTIDES).
 *   - insertion_reads++ iff pileup indel length > 0, same locus as the base.
 * Output column order (0..8) matches pileup_counts.COUNT_COLUMNS[0:9]
 * (reference_index, column 9, is filled by the Python wrapper from the FASTA,
 * as it already was in the original code -- the BAM contributes nothing to
 * that column).
 */
#include <Python.h>
#include <stdlib.h>
#include <string.h>
#include "htslib/sam.h"
#include "htslib/hts.h"

typedef struct {
    samFile *fp;
    hts_itr_t *iter;
} read_iter_state;

/* Replicates pysam's stepper="samtools" read-admission callback
 * (__advance_samtools in libcalignmentfile.pyx) for the subset of behaviour
 * pileup_counts.py actually relies on: default flag_filter
 * (FUNMAP|FSECONDARY|FQCFAIL|FDUP) and default ignore_orphans=True (drop
 * paired reads whose mate did not map properly-paired). No fastafile is
 * passed by providers.py's pileup() call, so BAQ / capped-MAPQ branches of
 * __advance_samtools never trigger there and are correctly omitted here.
 * min_mapping_quality at the iterator level defaults to 0 in pysam (the real
 * MAPQ>=20 filter happens per-read afterwards, replicated in count_region's
 * inner loop below), so it is not re-applied here either.
 */
static int read_bam_cb(void *data, bam1_t *b) {
    read_iter_state *st = (read_iter_state *)data;
    int ret;
    for (;;) {
        ret = sam_itr_next(st->fp, st->iter, b);
        if (ret < 0) break;
        uint16_t flag = b->core.flag;
        if (flag & (BAM_FUNMAP | BAM_FSECONDARY | BAM_FQCFAIL | BAM_FDUP)) continue;
        if ((flag & BAM_FPAIRED) && !(flag & BAM_FPROPER_PAIR)) continue; /* ignore_orphans */
        break;
    }
    return ret;
}

static inline int base_col(char base) {
    switch (base) {
        case 'A': return 0;
        case 'C': return 1;
        case 'G': return 2;
        case 'T': return 3;
        default: return -1;
    }
}

static PyObject *count_region(PyObject *self, PyObject *args) {
    const char *bam_path;
    const char *contig;
    long long start, end;
    int min_mapq, min_baseq;

    if (!PyArg_ParseTuple(args, "ssLLii", &bam_path, &contig, &start, &end,
                           &min_mapq, &min_baseq)) {
        return NULL;
    }
    if (end <= start) {
        return PyErr_Format(PyExc_ValueError, "end must be > start");
    }

    long long region_len = end - start;
    size_t n_doubles = (size_t)region_len * 9;
    double *out = (double *)calloc(n_doubles, sizeof(double));
    if (!out) {
        return PyErr_NoMemory();
    }
    /* column 9 (reference_index) is not written here; columns 0-8 only. */

    samFile *fp = sam_open(bam_path, "r");
    if (!fp) {
        free(out);
        return PyErr_Format(PyExc_IOError, "cannot open BAM: %s", bam_path);
    }
    sam_hdr_t *hdr = sam_hdr_read(fp);
    if (!hdr) {
        sam_close(fp);
        free(out);
        return PyErr_Format(PyExc_IOError, "cannot read header: %s", bam_path);
    }
    hts_idx_t *idx = sam_index_load(fp, bam_path);
    if (!idx) {
        sam_hdr_destroy(hdr);
        sam_close(fp);
        free(out);
        return PyErr_Format(PyExc_IOError, "cannot load BAM index: %s", bam_path);
    }
    int tid = sam_hdr_name2tid(hdr, contig);
    if (tid < 0) {
        hts_idx_destroy(idx);
        sam_hdr_destroy(hdr);
        sam_close(fp);
        free(out);
        return PyErr_Format(PyExc_ValueError, "contig not found: %s", contig);
    }

    hts_itr_t *iter = sam_itr_queryi(idx, tid, (hts_pos_t)start, (hts_pos_t)end);
    if (!iter) {
        hts_idx_destroy(idx);
        sam_hdr_destroy(hdr);
        sam_close(fp);
        free(out);
        return PyErr_Format(PyExc_IOError, "cannot create iterator");
    }

    read_iter_state st = {fp, iter};
    void *data_arr[1] = {&st};
    bam_mplp_t mplp = bam_mplp_init(1, read_bam_cb, data_arr);
    /* Match samtools/pysam default max pileup depth ceiling. */
    bam_mplp_set_maxcnt(mplp, 8000);
    /* Match pysam pileup()'s default ignore_overlaps=True: for overlapping
     * mate pairs, htslib zeroes the base quality of one mate's overlapping
     * bases so they fail the min_baseq filter below instead of being
     * double-counted. */
    bam_mplp_init_overlaps(mplp);

    int tid_out, pos, n_plp_arr[1];
    const bam_pileup1_t *pil_arr[1];
    const bam_pileup1_t *pil;
    int n_plp;

    Py_BEGIN_ALLOW_THREADS
    while (bam_mplp_auto(mplp, &tid_out, &pos, n_plp_arr, pil_arr) > 0) {
        pil = pil_arr[0];
        n_plp = n_plp_arr[0];
        if (pos < start || pos >= end) continue;
        double *row = out + (size_t)(pos - start) * 9;
        int count_A = 0, count_C = 0, count_G = 0, count_T = 0;
        int gap_reads = 0, insertion_reads = 0, depth = 0;
        double quality_sum = 0.0, mapping_sum = 0.0;

        for (int i = 0; i < n_plp; i++) {
            const bam_pileup1_t *p = pil + i;
            bam1_t *b = p->b;
            int mapq = b->core.qual;
            if (mapq < min_mapq) continue;
            uint16_t flag = b->core.flag;
            if (flag & (BAM_FDUP | BAM_FQCFAIL | BAM_FUNMAP)) continue;
            if (flag & (BAM_FSECONDARY | BAM_FSUPPLEMENTARY)) continue;

            depth++;
            mapping_sum += (double)mapq;

            if (p->is_del || p->is_refskip) {
                gap_reads++;
                continue;
            }
            int qpos = p->qpos;
            uint8_t *qual = bam_get_qual(b);
            int bq = (qual[0] == 0xff) ? 0 : (int)qual[qpos];
            if (bq < min_baseq) continue;

            uint8_t *seq = bam_get_seq(b);
            int code = bam_seqi(seq, qpos);
            char base = seq_nt16_str[code];
            int col = base_col(base);
            if (col < 0) continue;

            switch (col) {
                case 0: count_A++; break;
                case 1: count_C++; break;
                case 2: count_G++; break;
                case 3: count_T++; break;
            }
            quality_sum += (double)bq;
            if (p->indel > 0) insertion_reads++;
        }

        row[0] = count_A; row[1] = count_C; row[2] = count_G; row[3] = count_T;
        row[4] = gap_reads; row[5] = insertion_reads; row[6] = depth;
        row[7] = quality_sum; row[8] = mapping_sum;
    }
    Py_END_ALLOW_THREADS

    bam_mplp_destroy(mplp);
    hts_itr_destroy(iter);
    hts_idx_destroy(idx);
    sam_hdr_destroy(hdr);
    sam_close(fp);

    PyObject *bytes = PyBytes_FromStringAndSize((const char *)out, (Py_ssize_t)(n_doubles * sizeof(double)));
    free(out);
    return bytes;
}

static PyMethodDef PileupNativeMethods[] = {
    {"count_region", count_region, METH_VARARGS,
     "count_region(bam_path, contig, start, end, min_mapq, min_baseq) -> bytes\n"
     "Row-major float64 buffer, shape (end-start, 9), columns:\n"
     "count_A,count_C,count_G,count_T,gap_reads,insertion_reads,depth,quality_sum,mapping_sum"},
    {NULL, NULL, 0, NULL}
};

static struct PyModuleDef pileup_native_module = {
    PyModuleDef_HEAD_INIT, "pileup_native", NULL, -1, PileupNativeMethods
};

PyMODINIT_FUNC PyInit_pileup_native(void) {
    return PyModule_Create(&pileup_native_module);
}
