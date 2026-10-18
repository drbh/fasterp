import copy
import unittest

from benchmark_fastp import compare_reports


class ReportComparisonTests(unittest.TestCase):
    def test_integer_counters_must_match_exactly(self):
        reference = {"summary": {"before_filtering": {}, "after_filtering": {}},
                     "adapter_cutting": {"adapter_trimmed_bases": 1_000_000}}
        actual = copy.deepcopy(reference)
        actual["adapter_cutting"]["adapter_trimmed_bases"] += 1
        self.assertTrue(compare_reports(reference, actual))

    def test_upstream_float_rounding_is_allowed(self):
        reference = {"summary": {"before_filtering": {}, "after_filtering": {}},
                     "duplication": {"rate": 0.123457}}
        actual = copy.deepcopy(reference)
        actual["duplication"]["rate"] = 0.123456789
        self.assertEqual(compare_reports(reference, actual), [])


if __name__ == "__main__":
    unittest.main()
