import unittest
import os
from eukcc.eukcc import eukcc, eukcc_state, check_set_size


class test_write_result(unittest.TestCase):
    def test_lng_names_column(self):
        result = [
            {
                "fasta": "bin.fa",
                "quality": {"completeness": 90, "contamination": 1},
                "lng": ["1", "131567", "2759"],
                "lng_names": "root;cellular organisms;Eukaryota",
            },
            {"fasta": "empty.fa"},
        ]
        eukcc_instance = eukcc.__new__(eukcc)
        eukcc_instance.state = eukcc_state(options={"out": "."})
        eukcc_instance.write_result(result, ".test_tmp_eukcc.tsv")
        with open(".test_tmp_eukcc.tsv") as fin:
            lines = [line.rstrip("\n").split("\t") for line in fin]
        os.remove(".test_tmp_eukcc.tsv")
        os.remove("savestate.json.gz")

        self.assertEqual(lines[0][-1], "ncbi_lng_names")
        self.assertEqual(lines[1], ["bin.fa", "90", "1", "1-131567-2759", "root;cellular organisms;Eukaryota"])
        self.assertEqual(lines[2], ["empty.fa", "NA", "NA", "NA", "NA"])


class test_check_set_size(unittest.TestCase):
    def test_valid_sizes(self):
        # 0 means no limit, so it is valid with any minimal set size
        for max_set_size in [0, 20, 500]:
            check_set_size({"set_size": 20, "max_set_size": max_set_size})

    def test_max_smaller_than_min(self):
        with self.assertLogs(level="ERROR") as logs, self.assertRaises(SystemExit) as stop:
            check_set_size({"set_size": 20, "max_set_size": 10})
        self.assertEqual(stop.exception.code, 1)
        self.assertIn("--max_set_size (10) cannot be smaller than --set_size (20)", logs.output[0])
