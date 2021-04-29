import unittest
import os
from eukcc.treehandler import tree_sets

TESTDATA_SCMG_GZ = os.path.join(os.path.dirname(__file__), "testfiles", "SCMG_TEST.csv.gz")
TESTDATA_SCTREE = os.path.join(os.path.dirname(__file__), "testfiles", "tree_scmg.csv")
TESTDATA_TAXINFO = os.path.join(os.path.dirname(__file__), "testfiles", "tree_taxinfo.csv")


class test_tree(unittest.TestCase):
    def test_tree_sets_fail(self):
        simple_tree = "(BB:6.0,(A:5.0,C:3.0,EE:4.0):5.0,D:11.0);"
        placement = {"placements": [{"n": "EE"}, {"n": "BB"}]}
        t = tree_sets(tree_v=simple_tree, placement=placement, setp=TESTDATA_SCMG_GZ, taxinfo=TESTDATA_TAXINFO)
        self.assertIsNone(t.marker_set)

    def test_tree_sets_set(self):
        s = "((A,(B,AAA),((C,D,(E,DDD)),(F,G))),outgroup);"
        placement = {"placements": [{"n": "AAA"}, {"n": "DDD"}]}
        t = tree_sets(
            s, placement, TESTDATA_SCTREE, set_species=4, set_size=2, set_prevalence=80, taxinfo=TESTDATA_TAXINFO
        )
        self.assertEqual(t.marker_set.profiles, set(["Marker100", "Marker80"]))
