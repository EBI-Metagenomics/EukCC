import unittest
import os
from eukcc.base import union_sets, load_SCMGs, percentage_sets, load_tax_info

TESTDATA_SCMG = os.path.join(os.path.dirname(__file__), "testfiles", "SCMG_TEST.csv")
TESTDATA_SCMG_GZ = os.path.join(os.path.dirname(__file__), "testfiles", "SCMG_TEST.csv.gz")


class test_base(unittest.TestCase):
    def test_simple_union_sets(self):
        a = set([1, 2, 3])
        b = set([1, 2, 3, 4, 5])
        c = set([1, 2, 3, 8, 9])
        x = union_sets([a, b, c])
        self.assertEqual(a, x)

    def test_complex_union_sets(self):
        a = set([1, "x", "op"])
        b = a
        c = set(["y", "b", 11])
        x = union_sets([a, b, c])
        self.assertEqual(set(), x)

    def test_percentage_sets(self):
        a = set([1, 2, 3, 4, 5, 6, 7, 8, 9])
        b = set([1, 2, 3, 4, 5, 6, 7, 8, 10])
        c = set([1, 2, 3, 4, 5, 6, 7, 9, 10])
        d = set([1, 2, 3, 4, 5, 6, 8, 9, 10])
        e = set([1, 2, 3, 4, 5, 7, 8, 9, 10])
        sets = [a, b, c, d, e]
        target = set([1, 2, 3, 4, 5, 6, 7, 8, 9, 10])
        x = percentage_sets(sets, percent=80)
        self.assertEqual(target, x)
        # union_sets should produce the same as a 100% cutoff
        y = percentage_sets(sets, percent=100)
        self.assertEqual(union_sets(sets), y)
        # max limit, reproducinility
        x1 = percentage_sets(sets, percent=80, atmost=3)
        x2 = percentage_sets(sets, percent=80, atmost=3)
        self.assertEqual(x1, x2)

        # make sure the limit works
        n = 3
        x1 = percentage_sets(sets, percent=0, atmost=n)
        self.assertEqual(len(x1), n)

        # make sure its actually a set that is returned
        self.assertIsInstance(x1, set)

    def test_percentage_sets_weird_input(self):
        a = [1, 2, 3, 4, 5, 6, 7, 8, 9]
        b = [1, 2, 3, 4, 5, 6, 7, 8, 10]
        c = [1, 2, 3, 4, 5, 6, 7, 9, 10]
        d = [1, 2, 3, 4, 5, 6, 8, 9, 10]
        e = [1, 2, 3, 4, 5, 7, 8, 9, 10]
        sets = [a, b, c, d, e]
        target = set([1, 2, 3, 4, 5, 6, 7, 8, 9, 10])
        x = percentage_sets(sets, percent=80)
        self.assertEqual(target, x)

        # make sure the limit remover
        n = 0
        x1 = percentage_sets(sets, percent=80, atmost=n)
        self.assertEqual(x, x1)

    def test_error_union_sets(self):
        a = set([1, 2, 3])
        with self.assertRaises(TypeError):
            union_sets(a)

    def test_scmg_missing(self):
        with self.assertRaises(FileNotFoundError):
            load_SCMGs("adfhjdshjskfjf")

    def test_scmg_loading_gz(self):
        expected = {"A": set(["1", "2"]), "B": set(["1", "3", "2"]), "C": set(["1"])}
        found = load_SCMGs(TESTDATA_SCMG_GZ)
        self.assertEqual(expected, found)

    def test_scmg_loading_csv(self):
        csv = load_SCMGs(TESTDATA_SCMG)
        gz = load_SCMGs(TESTDATA_SCMG_GZ)
        self.assertEqual(csv, gz)

    def test_load_tax_info(self):
        TESTDATA_TAXINFO = os.path.join(os.path.dirname(__file__), "testfiles", "tree_taxinfo.csv")
        ti = load_tax_info(TESTDATA_TAXINFO)
        self.assertEqual(ti["B"], ["1"])


if __name__ == "__main__":
    unittest.main()
