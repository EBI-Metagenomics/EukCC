import unittest
from unittest import mock
from click.testing import CliRunner
from eukcc import __version__
import eukcc.__main__ as cli


class test_cli(unittest.TestCase):
    def setUp(self):
        self.runner = CliRunner()

    def test_version(self):
        result = self.runner.invoke(cli.main, ["--version"])
        self.assertEqual(result.exit_code, 0)
        self.assertEqual(result.output.strip(), "EukCC version {}".format(__version__))

    def test_help(self):
        for args in [["-h"], ["single", "-h"], ["folder", "-h"], ["define_set", "-h"], ["ncbi_update", "-h"]]:
            result = self.runner.invoke(cli.main, args)
            self.assertEqual(result.exit_code, 0, args)

    def test_no_command(self):
        result = self.runner.invoke(cli.main, [])
        self.assertEqual(result.exit_code, 0)
        self.assertIn("Commands:", result.output)

    def test_missing_argument(self):
        result = self.runner.invoke(cli.main, ["single"])
        self.assertEqual(result.exit_code, 2)
        self.assertIn("Missing argument 'FASTA'", result.output)

    def test_single_options(self):
        with self.runner.isolated_filesystem(), mock.patch.object(cli, "run_eukcc") as run:
            result = self.runner.invoke(
                cli.main,
                ["-d", "single", "bin.fa", "-o", "out", "--AA", "--genomes", "a.fa", "--genomes", "b.fa",
                 "--select_species", "--no_dynamic_root", "--keep", "--marker_prevalence", "97.5"],
            )
        self.assertEqual(result.exit_code, 0, result.output)
        args = vars(run.call_args.args[0])
        self.assertEqual(args["fasta"], "bin.fa")
        self.assertEqual(args["seqtype"], "AA")
        self.assertEqual(args["comparison_genomes"], ["a.fa", "b.fa"])
        self.assertIsNone(args["taxids"])
        self.assertEqual(args["set_selection"], "species")
        self.assertFalse(args["dynamic_root"])
        self.assertFalse(args["clean"])
        self.assertEqual(args["marker_prevalence"], 97.5)
        self.assertTrue(args["debug"])

    def test_single_defaults(self):
        with self.runner.isolated_filesystem(), mock.patch.object(cli, "run_eukcc") as run:
            result = self.runner.invoke(cli.main, ["single", "bin.fa"])
        self.assertEqual(result.exit_code, 0, result.output)
        args = vars(run.call_args.args[0])
        self.assertIsNone(args["seqtype"])
        self.assertIsNone(args["comparison_genomes"])
        self.assertEqual(args["set_selection"], "best_guess")
        self.assertTrue(args["dynamic_root"])
        self.assertTrue(args["clean"])
        self.assertEqual(args["out"], ".")
