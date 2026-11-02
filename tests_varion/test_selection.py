"""Unit tests for the pure parts of the threshold-selection procedure."""
import sys
import unittest
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "nindels/scripts/r31"))
import select_threshold as st  # noqa: E402
import records  # noqa: E402


class Rule(unittest.TestCase):
    def test_threshold_is_first_bin_from_which_everything_qualifies(self):
        items = [("lo", -0.1, -0.2), ("mid", 0.01, 0.002), ("hi", 0.02, 0.01)]
        self.assertEqual(st.apply_rule(items), 1)

    def test_a_later_failure_disqualifies_earlier_candidates(self):
        items = [("a", 0.01, 0.002), ("b", 0.01, -0.001), ("c", 0.02, 0.01)]
        self.assertEqual(st.apply_rule(items), 2)

    def test_positive_point_estimate_without_ci_support_is_not_enough(self):
        self.assertIsNone(st.apply_rule([("a", 0.01, -0.001), ("b", 0.02, -0.002)]))

    def test_never_qualifies(self):
        self.assertIsNone(st.apply_rule([("a", -0.01, -0.02)]))


class Bins(unittest.TestCase):
    def test_sparse_bins_merge_into_the_smaller_neighbour(self):
        edges, cnt = st.merge_sparse_bins([0, 10, 20, 30, 40], [50, 400, 500, 20], 300)
        self.assertEqual((edges, cnt), ([0, 20, 40], [450, 520]))

    def test_nothing_to_merge(self):
        self.assertEqual(st.merge_sparse_bins([0, 10, 20], [400, 500], 300), ([0, 10, 20], [400, 500]))


class Bootstrap(unittest.TestCase):
    def test_clear_improvement_gives_positive_interval(self):
        import random
        units = [("d", b) for b in range(20)]
        mk = lambda tp, fp, fn: {"TPt": tp, "FN": fn, "TPq": tp, "FP": fp}
        by = {"A": {("x", "l"): {u: mk(90, 10, 10) for u in units}}, "AB": {("x", "l"): {u: mk(95, 3, 5) for u in units}}}
        old = st.N_BOOT
        st.N_BOOT = 200
        try:
            out = st.bootstrap_delta(units, by, [("x", "l")], random.Random(1))
        finally:
            st.N_BOOT = old
        d, lo, hi, p = out[("x", "l")]
        self.assertGreater(d, 0.03)
        self.assertGreater(lo, 0)

    def test_f1_from_counts(self):
        self.assertAlmostEqual(records.f1({"TPt": 90, "FN": 10, "TPq": 90, "FP": 10}), 0.9)


if __name__ == "__main__":
    unittest.main()
