"""Exact equivalence of bisect-based _is_confident vs the original any() scan on real BEDs and random interval sets."""
import json, os, random, tempfile, time
from common import *
from providers import GiabAlignmentProvider, load_bed_regions
def orig(self, start, end):
    if self.high_confidence_bed is None: return True
    if self._confident_regions is None: self._confident_regions = load_bed_regions(self.high_confidence_bed, self.contig)
    return any(rs <= start and end <= re_ for rs, re_ in self._confident_regions)
def fresh(bed, contig):
    a = GiabAlignmentProvider.__new__(GiabAlignmentProvider)
    a.high_confidence_bed = bed; a.contig = contig; a._confident_regions = None; a._confident_index = None
    return a
def check_bed(bed, contig, lo, hi, n_random=40000):
    rng = random.Random(1); a = fresh(bed, contig); regions = load_bed_regions(bed, contig)
    starts = [rng.randint(lo, hi - 64) for _ in range(n_random)]
    for rs, re_ in regions: starts += [rs - 1, rs, rs + 1, re_ - 65, re_ - 64, re_ - 63]
    starts = [x for x in starts if lo <= x <= hi - 64]
    t = time.perf_counter(); new = [GiabAlignmentProvider._is_confident(a, x, x + 64) for x in starts]; t_new = time.perf_counter() - t
    t = time.perf_counter(); old = [orig(a, x, x + 64) for x in starts]; t_old = time.perf_counter() - t
    return dict(bed=bed.split("/")[-1], intervals=len(regions), queries=len(new), mismatches=sum(x != y for x, y in zip(new, old)),
                confident=sum(new), t_orig_s=round(t_old, 2), t_bisect_s=round(t_new, 3))
out = [check_bed(HG005["bed"], "chr1", 0, 25_000_000), check_bed(HG002_15X["bed"], "chr21", 30_000_000, 46_000_000),
       check_bed(CHR20["bed"], "chr20", 0, 64_444_167)]
rng = random.Random(0); bad = tot = 0
for trial in range(200):
    iv = sorted((s, s + rng.randint(1, 400)) for s in [rng.randint(0, 5000) for _ in range(rng.randint(0, 25))])
    if rng.random() < .3: rng.shuffle(iv)
    f = tempfile.NamedTemporaryFile("w", suffix=".bed", delete=False)
    for s, e in iv: f.write(f"chr1\t{s}\t{e}\n")
    f.close(); a = fresh(f.name, "chr1")
    for s in range(0, 5600, 7):
        for L in (1, 64):
            tot += 1; bad += GiabAlignmentProvider._is_confident(a, s, s + L) != orig(a, s, s + L)
    os.unlink(f.name)
out.append(dict(random_interval_sets=200, queries=tot, mismatches=bad))
for o in out: print(json.dumps(o))
