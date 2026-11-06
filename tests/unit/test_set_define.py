import unittest
from types import SimpleNamespace
from unittest import mock
import eukcc.set_define as set_define


class test_define_set(unittest.TestCase):
    def test_no_marker_set(self):
        args = SimpleNamespace(
            genomes=["genomes/a.fa", "genomes/b.fa"],
            db="db",
            out="out",
            name="set",
            set_size=20,
            max_set_size=500,
            marker_prevalence=95,
            seqtype=None,
        )
        # no marker gene set for these genomes
        with mock.patch.object(set_define, "eukcc") as eukcc:
            eukcc.return_value.pick_marker_set.return_value = None
            eukcc.return_value.terminate.side_effect = SystemExit
            with self.assertLogs(level="ERROR") as logs, self.assertRaises(SystemExit):
                set_define.define_set_w_genomes(args)
        eukcc.return_value.terminate.assert_called_once_with(1)
        self.assertIn("for the 2 genomes (a.fa, b.fa)", logs.output[0])
        self.assertIn("--set_size (now 20)", logs.output[0])
