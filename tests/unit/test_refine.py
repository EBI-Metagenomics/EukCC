import os
import tempfile
import unittest
from unittest import mock
import eukcc.refine as refine_module
from eukcc.eukcc import eukcc_state
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


class test_no_merge(unittest.TestCase):
    def run_refine(self, no_merge):
        with tempfile.TemporaryDirectory() as out:
            state = eukcc_state(
                workdir=out,
                options={"out": out, "input_bins": [], "links": None, "no_merge": no_merge, "n_combine": 1},
            )
            # no metaeuk, no bins: only the steps around the merging are run
            with mock.patch.object(refine_module, "eukcc"), \
                    mock.patch.object(refine_module, "merge_fasta"), \
                    mock.patch.object(refine_module, "split_contig_faa", return_value=[]), \
                    mock.patch.object(refine_module, "evaluate_multiples"):
                refine_module.refine(state)
            return sorted(os.listdir(out))

    def test_no_merge(self):
        files = self.run_refine(no_merge=True)
        self.assertIn("eukcc.tsv", files)
        self.assertNotIn("merged_bins.tsv", files)
        self.assertNotIn("merged_bins", files)

    def test_merge(self):
        files = self.run_refine(no_merge=False)
        self.assertIn("merged_bins.tsv", files)
        self.assertIn("merged_bins", files)


class test_reuse(unittest.TestCase):
    def test_split_contig_faa_rebuilds_folder(self):
        with tempfile.TemporaryDirectory() as workdir:
            # a file left by an older run with another bin
            os.makedirs(os.path.join(workdir, "faa"))
            open(os.path.join(workdir, "faa", "metaeuk_old.fa"), "w").close()
            faa = os.path.join(workdir, "proteins.faa")
            with open(faa, "w") as fout:
                fout.write(">metaeuk_bin_1.fa_binsplit_contig_1\nMKV\n")
            faas = refine_module.split_contig_faa(faa, workdir, delim="_binsplit_")
            self.assertEqual([os.path.basename(path) for path in faas], ["metaeuk_bin_1.fa"])

    def test_bin_cache_needs_same_proteins(self):
        with tempfile.TemporaryDirectory() as workdir:
            faa = os.path.join(workdir, "metaeuk_bin_1.fa")
            with open(faa, "w") as fout:
                fout.write(">contig_1\nMKV\n")
            state = eukcc_state(options={"out": workdir})
            # the bin folder is made by run_eukcc, which is mocked here
            wd = os.path.join(workdir, "bin")
            os.makedirs(wd)
            with mock.patch.object(refine_module.bin, "run_eukcc") as run:
                refine_module.bin(state, wd, faa, protein=True)
                # same proteins: the cached estimate is used
                refine_module.bin(state, wd, faa, protein=True)
                self.assertEqual(run.call_count, 1)
                # other proteins with the same file name: estimated again
                with open(faa, "w") as fout:
                    fout.write(">contig_1\nMKVL\n")
                refine_module.bin(state, wd, faa, protein=True)
                self.assertEqual(run.call_count, 2)
