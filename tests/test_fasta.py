import unittest
import os
from eukcc.fasta import Fasta, reduce_fasta, validate_fasta, N50, clean_metaeuk_fasta

TESTDATA_FASTA = os.path.join(os.path.dirname(__file__), "testfiles/test.fa")
TRICKY_FASTA = os.path.join(os.path.dirname(__file__), "testfiles/valid_but_tricky_fasta.fa")
NOT_A_FASTA = os.path.join(os.path.dirname(__file__), "testfiles/not_a_fasta.fasta")
CORRUPT_FILE = os.path.join(os.path.dirname(__file__), "testfiles/corrupt_file")
GMES_FILE = os.path.join(os.path.dirname(__file__), "testfiles/gmes_fail.fasta")
MULTILINE_FILE = os.path.join(os.path.dirname(__file__), "testfiles/multiline.fasta")

METAEUK_dirty = os.path.join(os.path.dirname(__file__), "testfiles/test_metaeuk_file.fa")
METAEUK_clean = os.path.join(os.path.dirname(__file__), "testfiles/test_metaeuk_file_clean.fa")


class test_base(unittest.TestCase):
    def test_Fasta(self):
        fa = Fasta(TESTDATA_FASTA)
        for seq in fa:
            if seq.name == "TEST":
                self.assertEqual(seq.seq, "HELLO")

    def test_Fasta_extract(self):
        fa = Fasta(TESTDATA_FASTA)
        for seq in fa:
            if seq.name == "TEST":
                break
        outfile = ".test_fasta"
        reduce_fasta(TESTDATA_FASTA, outfile, ["TEST"])
        fa2 = Fasta(outfile)
        for seq2 in fa2:
            if seq2.name == "TEST":
                break
        self.assertEqual(seq.seq, seq2.seq)
        self.assertEqual(seq.name, seq2.name)
        self.assertEqual(seq.long_name, seq2.long_name)
        os.remove(outfile)

    def test_validate_fasta(self):
        self.assertFalse(validate_fasta(NOT_A_FASTA))
        self.assertFalse(validate_fasta(CORRUPT_FILE))
        self.assertFalse(validate_fasta(GMES_FILE))
        self.assertTrue(validate_fasta(TESTDATA_FASTA))
        self.assertTrue(validate_fasta(TRICKY_FASTA))

    def test_open_invalid(self):
        with self.assertRaises(ValueError):
            x = Fasta(NOT_A_FASTA)
            for s in x:
                print(s)
        with self.assertRaises(ValueError):
            x = Fasta(CORRUPT_FILE)
            for s in x:
                print(s)

    def test_n50(self):
        # wikipedia example
        fasta_1 = os.path.join(os.path.dirname(__file__), "testfiles/N50_1.fa")
        self.assertEqual(N50(fasta_1), 8)

    def test_multiline(self):
        expected = ["AAGGCCTT", "CCTTAAGG", "AAGGCCTT"]
        for e, seq in zip(expected, Fasta(MULTILINE_FILE)):
            self.assertEqual(e, seq.seq)

    def test_metaeuk_cleaner(self):
        cleaned = Fasta(clean_metaeuk_fasta(METAEUK_dirty, ".test_tmp_metaeuk_file"))
        clean = Fasta(METAEUK_clean)

        for seq_1, seq_2 in zip(clean, cleaned):
            self.assertEqual(seq_1.name, seq_2.name)
            self.assertEqual(seq_1.seq, seq_2.seq)
        os.remove(".test_tmp_metaeuk_file")
