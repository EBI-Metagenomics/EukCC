import os
import tempfile
import unittest
from eukcc.eukcc import eukcc
from eukcc.treehandler import subdb_vote

# QCho: choanoflagellate (Opisthokonta) whose neighbours NCBI publishes in protozoa
# QFun: sister of a fungi clade, so the vote has to walk up one node
# QRed: red alga, sister to the plants in the tree, but Rhodophyta is not
#       Viridiplantae, so its neighbours sit in protozoa
# QPla: green plant, whose neighbours sit in plant
TREE = (
    "(((ChoA:1,ChoB:1,QCho:1):1,((FunA:1,FunB:1):1,QFun:1):1):1,"
    "((RedA:1,RedB:1,QRed:1):1,(PlaA:1,PlaB:1,QPla:1):1):1);"
)

LEAF_SUBDB = {
    "FunA": "fungi",
    "FunB": "fungi",
    "RedA": "protozoa",
    "RedB": "protozoa",
    "PlaA": "plant",
    "PlaB": "plant",
    "ChoA": "protozoa",
    "ChoB": "protozoa",
}


class test_subdb_vote(unittest.TestCase):
    def test_opisthokont_protist(self):
        self.assertEqual(subdb_vote(TREE, ["QCho"], LEAF_SUBDB), ("protozoa", 1))

    def test_fungi(self):
        self.assertEqual(subdb_vote(TREE, ["QFun"], LEAF_SUBDB), ("fungi", 1))

    def test_red_alga(self):
        self.assertEqual(subdb_vote(TREE, ["QRed"], LEAF_SUBDB), ("protozoa", 1))

    def test_plant(self):
        self.assertEqual(subdb_vote(TREE, ["QPla"], LEAF_SUBDB), ("plant", 1))

    def test_split_votes(self):
        clade, frac = subdb_vote(TREE, ["QCho", "QFun", "QPla"], LEAF_SUBDB)
        self.assertIsNone(clade)
        self.assertAlmostEqual(frac, 1 / 3)

    def test_no_placements(self):
        self.assertEqual(subdb_vote(TREE, [], LEAF_SUBDB), (None, 0))


class test_determine_subdb(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        # minimal database: db_<sub>/db_info.json and db_<sub>/taxinfo.csv per sub-database
        # hmm_placement_place must exist, but is not used, so it points to the sub-database folder
        cls.db = tempfile.TemporaryDirectory()
        for leaf, subdb in LEAF_SUBDB.items():
            subdb_dir = os.path.join(cls.db.name, "db_{}".format(subdb))
            os.makedirs(subdb_dir, exist_ok=True)
            with open(os.path.join(subdb_dir, "db_info.json"), "w") as fout:
                fout.write('{"files": {"taxinfo": "taxinfo.csv", "hmm_placement_place": "."}}')
            with open(os.path.join(subdb_dir, "taxinfo.csv"), "a") as fout:
                fout.write("0,{}\n".format(leaf))

    @classmethod
    def tearDownClass(cls):
        cls.db.cleanup()

    def determine_subdb(self, places):
        eukcc_instance = eukcc.__new__(eukcc)
        eukcc_instance.state = {
            "db": self.db.name,
            "placement": {
                "tog": TREE,
                "placements": [{"n": [place]} for place in places],
            },
        }
        return eukcc_instance.determine_subdb()

    def test_single_placement(self):
        self.assertEqual(self.determine_subdb(["QCho"]), "protozoa")

    def test_majority(self):
        # QFun represents 2/3 of the placements, so is 67% of the votes, which is above the 60% threshold
        self.assertEqual(self.determine_subdb(["QFun", "QFun", "QPla"]), "fungi")

    def test_no_majority(self):
        # QCho, QFun and QPla each represent 1/3 of the placements, so no majority is reached
        self.assertEqual(self.determine_subdb(["QCho", "QFun", "QPla"]), "base")
