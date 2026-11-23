import sys
import csv
import io
import unittest
from contextlib import redirect_stdout
from pathlib import Path
from shutil import copyfile
from tempfile import TemporaryDirectory


PROJECT_DIR = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(PROJECT_DIR / "src"))

from fasta_io import gc_content, main, n_invalid, normalize, read_fasta, validate, write_fasta


class TestFastaIO(unittest.TestCase):
    def test_clean_fasta_expected_stats(self):
        expected_stats = [
            ("TEM-1", 861, 422, 49.0128),
            ("TEM-2", 861, 424, 49.2451),
            ("TEM-3", 861, 422, 49.0128),
            ("TEM-4", 861, 421, 48.8966),
        ]

        for name, expected_length, expected_gc_count, expected_gc_percent in expected_stats:
            with self.subTest(variant=name):
                clean_path = PROJECT_DIR / "data" / "clean" / f"{name}.fasta"
                raw_path = PROJECT_DIR / "data" / "raw" / f"{name}_cds.fasta"

                self.assertTrue(clean_path.is_file(), f"Missing clean FASTA: {clean_path}")
                header, seq = read_fasta(clean_path)
                _, raw_seq = read_fasta(raw_path)

                self.assertTrue(header, "Header must not be empty")
                self.assertEqual(seq, normalize(seq), "Sequence must be normalized")
                validate(seq)
                self.assertEqual(seq, normalize(raw_seq), "Clean sequence differs from raw CDS")
                self.assertEqual(len(seq), expected_length)
                self.assertEqual(seq.count("G") + seq.count("C"), expected_gc_count)
                # Expected percentages are rounded to four decimal places.
                self.assertAlmostEqual(gc_content(seq), expected_gc_percent, delta=0.00005)

    def test_read_fasta(self):
        with TemporaryDirectory() as tmp_dir:
            path = Path(tmp_dir) / "sample.fasta"
            path.write_text(">TEM test\nATG\nCCT\n", encoding="utf-8")

            self.assertEqual(read_fasta(path), ("TEM test", "ATGCCT"))

    def test_read_fasta_rejects_invalid_header(self):
        with TemporaryDirectory() as tmp_dir:
            path = Path(tmp_dir) / "invalid.fasta"
            path.write_text("ATGC\n", encoding="utf-8")

            with self.assertRaises(ValueError):
                read_fasta(path)

    def test_read_fasta_rejects_multiple_records(self):
        with TemporaryDirectory() as tmp_dir:
            path = Path(tmp_dir) / "multiple.fasta"
            path.write_text(">one\nATG\n>two\nCCT\n", encoding="utf-8")

            with self.assertRaises(ValueError):
                read_fasta(path)

    def test_normalize_and_validate(self):
        seq = normalize(" atgc ")
        self.assertEqual(seq, "ATGC")
        self.assertIsNone(validate(seq))

    def test_validate_reports_invalid_character(self):
        # Only B is invalid in ABC; the current API uses zero-based positions.
        with self.assertRaises(ValueError) as caught:
            validate("ABC")
        self.assertIn("{1: 'B'}", str(caught.exception))


    def test_normalize_removes_all_whitespace_but_preserves_invalid_characters(self):
        self.assertEqual(normalize(" aT\tg\nC "), "ATGC")
        self.assertEqual(normalize(" aN-1\t "), "AN-1")
        self.assertEqual(normalize(" aßt- "), "AßT-")
        self.assertEqual(normalize(normalize(" aT\tg\nC ")), "ATGC")

    def test_read_fasta_rejects_empty_records(self):
        with TemporaryDirectory() as tmp_dir:
            path = Path(tmp_dir) / "empty.fasta"
            for content in ("", ">\nATGC\n", ">   \nATGC\n", ">one\n \t\n"):
                with self.subTest(content=content):
                    path.write_text(content, encoding="utf-8")
                    with self.assertRaises(ValueError):
                        read_fasta(path)

    def test_gc_percent_and_invalid_count(self):
        for seq, expected in (("ATAT", 0.0), ("ACGT", 50.0), ("GGCC", 100.0)):
            with self.subTest(seq=seq):
                self.assertEqual(gc_content(seq), expected)
                self.assertEqual(n_invalid(seq), 0)
        self.assertEqual(n_invalid("ATGN-"), 2)
        for seq in ("", "ATGN-"):
            with self.subTest(invalid=seq):
                with self.assertRaises(ValueError):
                    gc_content(seq)

    def test_write_fasta_round_trip(self):
        with TemporaryDirectory() as tmp_dir:
            path = Path(tmp_dir) / "round_trip.fasta"
            write_fasta(path, "TEM-1", "ATGC")
            self.assertEqual(path.read_text(encoding="utf-8"), ">TEM-1\nATGC\n")
            self.assertEqual(read_fasta(path), ("TEM-1", "ATGC"))


class TestFastaPipeline(unittest.TestCase):
    def setUp(self):
        tmp_dir = TemporaryDirectory()
        self.addCleanup(tmp_dir.cleanup)
        self.root = Path(tmp_dir.name)
        self.raw_dir = self.root / "data" / "raw"
        self.raw_dir.mkdir(parents=True)

    def write_input(self, name, content):
        (self.raw_dir / f"{name}_cds.fasta").write_text(content, encoding="utf-8")

    def run_pipeline(self):
        output = io.StringIO()
        with redirect_stdout(output):
            main(self.root)
        return output.getvalue()

    def test_csv_and_fasta_outputs_are_normalized_and_repeatable(self):
        self.write_input("TEM-2", ">original header\nggCC\n")
        self.write_input("TEM-1", ">original header\n aT\tg\nC \n")
        (self.raw_dir / "TEM-1_full.fasta").write_text("not a CDS", encoding="utf-8")
        originals = {path: path.read_bytes() for path in self.raw_dir.iterdir()}

        self.run_pipeline()
        stats_path = self.root / "results" / "seq_stats.csv"
        with stats_path.open(newline="", encoding="utf-8") as file:
            self.assertEqual(list(csv.reader(file)), [
                ["name", "length", "gc_percent", "n_invalid"],
                ["TEM-1", "4", "50.0000", "0"],
                ["TEM-2", "4", "100.0000", "0"],
            ])
        clean_dir = self.root / "data" / "clean"
        self.assertEqual((clean_dir / "TEM-1.fasta").read_text(), ">TEM-1\nATGC\n")
        self.assertEqual((clean_dir / "TEM-2.fasta").read_text(), ">TEM-2\nGGCC\n")
        outputs = {path: path.read_bytes() for path in (stats_path, *clean_dir.iterdir())}
        self.run_pipeline()
        for path, content in {**originals, **outputs}.items():
            self.assertEqual(path.read_bytes(), content)

    def test_invalid_sequence_is_reported_and_stale_clean_file_is_removed(self):
        self.write_input("TEM-1", ">one\nATGN-\n")
        self.write_input("TEM-2", ">two\nACGT\n")
        clean_dir = self.root / "data" / "clean"
        clean_dir.mkdir()
        stale = clean_dir / "TEM-1.fasta"
        stale.write_text(">TEM-1\nAAAA\n")
        unrelated = clean_dir / "notes.txt"
        unrelated.write_text("keep me")

        output = self.run_pipeline()
        with (self.root / "results" / "seq_stats.csv").open(newline="", encoding="utf-8") as file:
            self.assertEqual(list(csv.DictReader(file)), [
                {"name": "TEM-1", "length": "5", "gc_percent": "", "n_invalid": "2"},
                {"name": "TEM-2", "length": "4", "gc_percent": "50.0000", "n_invalid": "0"},
            ])
        self.assertFalse(stale.exists())
        self.assertEqual(read_fasta(clean_dir / "TEM-2.fasta"), ("TEM-2", "ACGT"))
        self.assertEqual(unrelated.read_text(), "keep me")
        self.assertIn("TEM-1", output)
        self.assertIn("N", output)

    def test_malformed_input_leaves_existing_outputs_unchanged(self):
        self.write_input("TEM-1", ">one\nATGC\n")
        clean_dir = self.root / "data" / "clean"
        clean_dir.mkdir()
        clean = clean_dir / "TEM-1.fasta"
        clean.write_text(">old\nAAAA\n")
        results = self.root / "results"
        results.mkdir()
        stats = results / "seq_stats.csv"
        stats.write_text("previous report\n")
        for content in ("ATGC\n", ">two\n \t\n", ">two\nATG\n>three\nCCC\n"):
            with self.subTest(content=content):
                self.write_input("TEM-2", content)
                with self.assertRaises(ValueError):
                    self.run_pipeline()
                self.assertEqual(clean.read_text(), ">old\nAAAA\n")
                self.assertEqual(stats.read_text(), "previous report\n")

    def test_missing_inputs_raise_without_creating_outputs(self):
        with self.assertRaises(FileNotFoundError):
            self.run_pipeline()
        self.assertFalse((self.root / "results").exists())
        self.assertFalse((self.root / "data" / "clean").exists())

    def test_non_ascii_invalid_characters_keep_their_length_and_positions(self):
        self.write_input("TEM-1", ">one\naßt-\n")
        output = self.run_pipeline()
        with (self.root / "results" / "seq_stats.csv").open(newline="", encoding="utf-8") as file:
            self.assertEqual(list(csv.DictReader(file)), [
                {"name": "TEM-1", "length": "4", "gc_percent": "", "n_invalid": "2"},
            ])
        self.assertIn("{1: 'ß', 3: '-'}", output)
        self.assertFalse((self.root / "data" / "clean" / "TEM-1.fasta").exists())

    def test_real_cds_produce_expected_csv(self):
        for name in ("TEM-1", "TEM-2", "TEM-3", "TEM-4"):
            copyfile(PROJECT_DIR / "data" / "raw" / f"{name}_cds.fasta",
                     self.raw_dir / f"{name}_cds.fasta")
        self.run_pipeline()
        self.assertEqual((self.root / "results" / "seq_stats.csv").read_text(),
                         "name,length,gc_percent,n_invalid\n"
                         "TEM-1,861,49.0128,0\n"
                         "TEM-2,861,49.2451,0\n"
                         "TEM-3,861,49.0128,0\n"
                         "TEM-4,861,48.8966,0\n")


if __name__ == "__main__":
    unittest.main(verbosity=2)
