"""Tests for R3.1 mode selection (automatic and manual) and the CLI front end, using a fake engine."""
import json
import os
import stat
import subprocess
import sys
import tempfile
import textwrap
import unittest
from pathlib import Path

REPO = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(REPO))

from varion.adaptive import Thresholds, choose_config  # noqa: E402
from varion import cli  # noqa: E402

T = Thresholds(global_min_depth=40.0, local_min_depth=25, status="test", provenance="unit test")
T_NONE = Thresholds(global_min_depth=None, local_min_depth=None, status="no depth found", provenance="unit test")


class ChooseConfig(unittest.TestCase):
    def test_auto_low_depth_uses_A(self):
        d = choose_config("auto", 15.0, T)
        self.assertEqual((d.config, d.use_repeat_model, d.str_min_depth, d.source), ("A", False, 0, "auto"))

    def test_auto_at_threshold_uses_AB_with_local_gate(self):
        d = choose_config("auto", 40.0, T)
        self.assertEqual((d.config, d.use_repeat_model, d.str_min_depth), ("A+B", True, 25))

    def test_auto_just_below_threshold(self):
        self.assertEqual(choose_config("auto", 39.99, T).config, "A")

    def test_auto_high_depth(self):
        self.assertEqual(choose_config("auto", 250.0, T).config, "A+B")

    def test_auto_without_supported_threshold_never_switches(self):
        for depth in (5, 60, 500):
            self.assertEqual(choose_config("auto", depth, T_NONE).config, "A")
        self.assertIn("no depth", choose_config("auto", 60, T_NONE).reason)

    def test_auto_unmeasurable_depth_falls_back_to_A(self):
        self.assertEqual(choose_config("auto", None, T).config, "A")

    def test_manual_override_ignores_depth(self):
        self.assertEqual(choose_config("A", 500.0, T).config, "A")
        d = choose_config("A+B", 5.0, T)
        self.assertEqual((d.config, d.use_repeat_model, d.str_min_depth, d.source), ("A+B", True, 0, "manual"))

    def test_invalid_mode(self):
        with self.assertRaises(ValueError):
            choose_config("fast", 30.0, T)

    def test_thresholds_file_roundtrip(self):
        with tempfile.TemporaryDirectory() as td:
            p = Path(td) / "t.json"
            p.write_text(json.dumps({"global_min_depth": 33, "local_min_depth": 20, "status": "x", "provenance": "y"}))
            t = Thresholds.load(p)
            self.assertEqual((t.global_min_depth, t.local_min_depth), (33, 20))
        self.assertIsInstance(Thresholds.load(), Thresholds)  # the shipped file parses


FAKE_ENGINE = textwrap.dedent("""\
    #!/usr/bin/env python3
    import json, sys
    cmd, a = sys.argv[1], sys.argv[2:]
    opt = lambda k: a[a.index(k) + 1] if k in a else None
    if cmd == "contigs":
        print("#sample\\tNA12878"); print("chr20\\t64444167\\t1000"); print("chr21\\t46709983\\t500"); print("chrUn\\t100\\t0")
    elif cmd == "depth":
        d = {"chr20": 12.0, "chr21": 80.0}[opt("--region")]
        json.dump({"region": opt("--region"), "window": 10000, "windows_sampled": 400, "windows_empty": 3,
                   "median": d, "mean": d, "p10": d * .9, "p25": d, "p75": d, "p90": d * 1.1, "cv": 0.1}, open(opt("--out"), "w"))
    elif cmd == "estimate-str":
        open(opt("--out"), "w").write("#period\\tunits\\tloci\\teps\\tgap_open\\n")
    elif cmd == "call":
        used = "--str-model" in a
        open(opt("--out-vcf"), "w").write("##fileformat=VCFv4.2\\n#CHROM\\tPOS\\tID\\tREF\\tALT\\tQUAL\\tFILTER\\tINFO\\tFORMAT\\tNA12878\\n"
                                       + opt("--region") + "\\t100\\t.\\tA\\tAT\\t50\\tPASS\\tK=1;DP=30\\tGT\\t0/1\\n")
        json.dump({"genotyper": opt("--genotyper"), "repeat_model": used, "str_min_depth": int(opt("--str-min-depth") or 0) if used else -1,
                   "loci": 1, "pass": 1, "loci_with_repeat_model": 1 if used else 0}, open(opt("--summary-json"), "w"))
""")


class CliWithFakeEngine(unittest.TestCase):
    def setUp(self):
        self.td = tempfile.TemporaryDirectory()
        td = Path(self.td.name)
        self.engine = td / "fake-nindels"
        self.engine.write_text(FAKE_ENGINE)
        self.engine.chmod(self.engine.stat().st_mode | stat.S_IEXEC)
        (td / "s.bam").write_bytes(b"x"); (td / "s.bam.bai").write_bytes(b"x")
        (td / "r.fa").write_text(">chr20\nACGT\n>chr21\nACGT\n")
        (td / "r.fa.fai").write_text("chr20\t64444167\t6\t60\t61\nchr21\t46709983\t6\t60\t61\n")
        (td / "t.json").write_text(json.dumps({"global_min_depth": 40, "local_min_depth": 25, "status": "test", "provenance": "t"}))
        self.td_path = td

    def tearDown(self):
        self.td.cleanup()

    def run_cli(self, *extra):
        td = self.td_path
        argv = ["indel", "--bam", str(td / "s.bam"), "--ref", str(td / "r.fa"), "--out-dir", str(td / "out"),
                "--engine", str(self.engine), "--thresholds", str(td / "t.json"), "--quiet", "--threads", "2", *extra]
        return cli.main(argv)

    def summary(self):
        return json.loads(next((self.td_path / "out").glob("*.summary.json")).read_text())

    def test_auto_picks_per_contig_and_records_it(self):
        self.assertEqual(self.run_cli(), 0)
        s = self.summary()
        by = {c["contig"]: c for c in s["contigs"]}
        self.assertEqual(sorted(by), ["chr20", "chr21"])           # chrUn has no reads and is not a primary chromosome
        self.assertEqual(by["chr20"]["decision"]["config"], "A")    # 12x < 40x
        self.assertEqual(by["chr21"]["decision"]["config"], "A+B")  # 80x >= 40x
        self.assertFalse(by["chr20"]["call"]["repeat_model"])
        self.assertTrue(by["chr21"]["call"]["repeat_model"])
        self.assertEqual(by["chr21"]["call"]["str_min_depth"], 25)
        self.assertEqual(s["mode_requested"], "auto")
        vcf = (self.td_path / "out" / "NA12878.indels.vcf").read_text()
        self.assertIn("##varion_config_chr20=A;", vcf)
        self.assertIn("##varion_config_chr21=A+B;", vcf)
        self.assertEqual(s["sample"], "NA12878")

    def test_manual_overrides(self):
        self.assertEqual(self.run_cli("--mode", "A"), 0)
        self.assertTrue(all(c["decision"]["config"] == "A" and c["decision"]["source"] == "manual" for c in self.summary()["contigs"]))
        self.assertEqual(self.run_cli("--mode", "A+B", "--force", "--region", "chr20"), 0)
        c = self.summary()["contigs"][0]
        self.assertEqual((c["decision"]["config"], c["call"]["repeat_model"], c["call"]["str_min_depth"]), ("A+B", True, 0))

    def test_refuses_to_overwrite_without_force(self):
        self.assertEqual(self.run_cli(), 0)
        self.assertEqual(self.run_cli(), 2)

    def test_input_validation(self):
        td = self.td_path
        self.assertEqual(cli.main(["indel", "--bam", str(td / "missing.bam"), "--ref", str(td / "r.fa"), "--out-dir", str(td / "o"),
                                   "--engine", str(self.engine)]), 2)
        (td / "r.fa.fai").unlink()
        self.assertEqual(self.run_cli(), 2)
        (td / "r.fa.fai").write_text("chr20\t100\t6\t60\t61\n")
        self.assertEqual(self.run_cli("--region", "chr99"), 2)   # not in BAM
        self.assertEqual(self.run_cli("--region", "chr21"), 2)   # in BAM, not in the reference index

    def test_dry_run_runs_nothing(self):
        self.assertEqual(self.run_cli("--dry-run"), 0)
        self.assertFalse((self.td_path / "out").exists())


if __name__ == "__main__":
    unittest.main()
