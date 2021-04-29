import unittest
import os
from eukcc.file import file

TESTDATA_SCMG = os.path.join(os.path.dirname(__file__), "testfiles", "SCMG_TEST.csv")
TESTDATA_dir = os.path.join(os.path.dirname(__file__), "testfiles")


class test_file_class(unittest.TestCase):
    def test_file_isdir(self):
        self.assertTrue(file.isdir(TESTDATA_dir))

    def test_file_isfile(self):
        self.assertTrue(file.isfile(TESTDATA_SCMG))
