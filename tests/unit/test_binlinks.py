import os
import shutil
import tempfile
import unittest
from scripts.binlinks import bin_map, contig_map

FIXTURES = os.path.join(os.path.dirname(__file__), "../fixtures")


class test_binlinks(unittest.TestCase):
    def test_bins_like_eukcc_folder(self):
        with tempfile.TemporaryDirectory() as bindir:
            # test.fa.gz has the same contigs as test.fa, so each bin gets its own copy with new names
            for name, fixture in [("bin_1.fa", "test.fa"), ("bin_2.FASTA", "test.fa"), ("bin_3.fna.gz", "test.fa.gz")]:
                shutil.copy(os.path.join(FIXTURES, fixture), os.path.join(bindir, name))
            shutil.copy(os.path.join(FIXTURES, "test.fa"), os.path.join(bindir, "notes.txt"))

            contigs, contigs_per_bin = bin_map(bindir)
            self.assertEqual(set(contigs_per_bin), {"bin_1.fa", "bin_2.FASTA", "bin_3.fna.gz"})
            self.assertEqual(set(contigs_per_bin.values()), {3})
            self.assertEqual(contig_map(bindir)["sequence"], 9)
