import os
import unittest
from eukcc.refine import bins_without_proteins, write_table, write_no_marker_genes_output


class test_refine(unittest.TestCase):
    def test_bins_without_proteins(self):
        input_bins = ["bins/bin_1.fa", "bins/bin_2.fa.gz", "bins/bin_3.fna"]
        faas = ["workdir/faa/metaeuk_bin_1.fa", "workdir/faa/metaeuk_bin_3.fna"]
        self.assertEqual(bins_without_proteins(input_bins, faas), ["bins/bin_2.fa.gz"])

    def test_unestimated_bin_rows(self):
        table = ".test_tmp_refine.tsv"
        missing = ".test_tmp_missing.txt"
        write_table(None, table, header=True)
        write_no_marker_genes_output(None, missing, header=True)
        # a bin with zero proteins only has its fasta path
        write_table({"fasta": "bins/bin_2.fa.gz"}, table)
        write_no_marker_genes_output({"fasta": "bins/bin_2.fa.gz"}, missing)
        # a bin with proteins but no marker gene set
        write_table({"fasta": "workdir/faa/metaeuk_bin_4.fa", "quality": None}, table)
        with open(table) as fin:
            rows = [line.rstrip("\n").split("\t") for line in fin]
        with open(missing) as fin:
            names = [line.rstrip("\n") for line in fin]
        os.remove(table)
        os.remove(missing)

        self.assertEqual(rows[1], ["bin_2.fa.gz", "NA", "NA", "NA", "NA"])
        self.assertEqual(rows[2], ["bin_4.fa", "NA", "NA", "NA", "NA"])
        self.assertEqual(names, ["bin", "bin_2.fa.gz"])
