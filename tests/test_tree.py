import unittest
import os
from eukcc.treehandler import tree_sets, tax_LCA

TESTDATA_SCMG_GZ = os.path.join(os.path.dirname(__file__), "testfiles", "SCMG_TEST.csv.gz")
TESTDATA_SCTREE = os.path.join(os.path.dirname(__file__), "testfiles", "tree_scmg.csv")
TESTDATA_TAXINFO = os.path.join(os.path.dirname(__file__), "testfiles", "tree_taxinfo.csv")
TESTDATA_TAX_LCA = os.path.join(os.path.dirname(__file__), "testfiles", "tax_lca_taxinfo.csv")


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

    def test_tax_LCA_partial_placements(self):
        tree = (
            "(((GCA_015473125.1,metaeuk_OY253681.1_29),(((GCA_006384855.1,metaeuk_OY253673.1_71),"
            "metaeuk_OY253685.1_62),((((metaeuk_OY253682.1_41,GCA_004138255.1),metaeuk_OY253672.1_198),"
            "((GCA_011766145.1,metaeuk_OY253676.1_95),(GCA_004335775.1,metaeuk_OY253680.1_173))),"
            "((GCA_010646915.1,metaeuk_OY253673.1_210),((metaeuk_OY253668.1_317,GCA_002794665.1),"
            "GCA_000690575.1))))),metaeuk_OY253675.1_91);"
        )
        lng = tax_LCA(tree, TESTDATA_TAX_LCA, placements=["metaeuk_OY253681.1_29"], add_protist_common=False)

        # root, cellular organisms, Eukaryota, Viridiplantae, Chlorophyta, Pseudoscourfieldiophyceae, Chlorophyta incertae sedis, Pycnococcaceae
        self.assertEqual(lng, ["1", "131567", "2759", "33090", "3041", "3417985", "3417986", "41878"])
