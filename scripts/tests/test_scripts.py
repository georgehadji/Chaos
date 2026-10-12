"""Offline tests for the ICM scripts. Run: python -m unittest discover -s scripts/tests"""
import sys
import tempfile
import unittest
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT / "scripts"))

import bibcheck  # noqa: E402
import icm_status  # noqa: E402
import icm_validate  # noqa: E402
import promote  # noqa: E402

RECORD = {"title": ["A universal instability of many-dimensional oscillator systems"],
          "author": [{"family": "Chirikov", "given": "B. V."}, {"name": "Some Institute"}],
          "volume": "52", "page": "263–379", "DOI": "10.1016/0370-1573(79)90023-1",
          "container-title": ["Physics Reports"], "issued": {"date-parts": [[1979]]}}
ENTRY = {"key": "Chirikov1979", "title": "A universal instability of many-dimensional oscillator systems",
         "author": "Chirikov, B. V.", "volume": "52(5)", "pages": "263-379"}


class BibcheckTests(unittest.TestCase):
    def test_matching_record_has_no_issues(self):
        self.assertEqual(bibcheck.compare(ENTRY, RECORD), [])

    def test_wrong_title_and_pages_are_reported(self):
        issues = bibcheck.compare(dict(ENTRY, title="Something else", pages="1-2"), RECORD)
        self.assertTrue(any(i.startswith("title") for i in issues))
        self.assertTrue(any(i.startswith("pages") for i in issues))

    def test_key_rule(self):
        self.assertEqual(bibcheck.key_base(["Hénon", "Heiles"], 1964), "HenonHeiles1964")
        self.assertEqual(bibcheck.key_base(["A", "B", "C", "D"], 2000), "A2000")
        self.assertEqual(bibcheck.key_base(["FARANDA", "Udwadia", None, "None"], 2012), "FarandaUdwadia2012")
        self.assertEqual(bibcheck.assign_keys(["W1993", "X2000", "W1993"]), ["W1993a", "X2000", "W1993b"])

    def test_year_falls_back_when_issued_is_empty(self):
        rec = {"issued": {"date-parts": [[None]]}, "created": {"date-parts": [[2002, 11, 23]]},
               "event": {"name": "1997 1st International Conference, Control of Oscillations and Chaos"}}
        self.assertEqual(bibcheck.rec_year(rec), 1997)
        self.assertEqual(bibcheck.rec_year({"issued": {"date-parts": [[None]]}, "created": {"date-parts": [[2002]]}}), 2002)

    def test_arxiv_ids_new_and_old_style(self):
        text = "- arXiv:2304.08340\n- arXiv:acc-phys/9411001\n- arXiv:chao-dyn/9406010\n"
        self.assertEqual(bibcheck.ARXIV_RE.findall(text), ["2304.08340", "acc-phys/9411001", "chao-dyn/9406010"])

    def test_bibtex_round_trip(self):
        entry = bibcheck.parse_bib(bibcheck.to_bibtex("Chirikov1979", RECORD))[0]
        self.assertEqual(entry["key"], "Chirikov1979")
        self.assertEqual(bibcheck.compare(entry, RECORD), [])

    def test_em_dash_written_as_bibtex_markup(self):
        rec = dict(RECORD, title=["Chaos" + chr(0x2014) + "order"])
        text = bibcheck.to_bibtex("X2000", rec)
        self.assertNotIn(chr(0x2014), text)
        self.assertIn("Chaos---order", text)


class StatusTests(unittest.TestCase):
    def test_slug_detection_handles_hyphenated_slugs(self):
        with tempfile.TemporaryDirectory() as d:
            out = Path(d)
            (out / "p1-le-benchmark-preregistration.md").write_text("x")
            (out / "p1-le-benchmark-fig-1.pdf").write_text("x")
            pats = icm_status.stage_patterns("## Outputs\n| a | `output/[study]-preregistration.md` | f |\n"
                                             "| b | `output/[study]-fig-*.pdf` | f |\n")
            self.assertEqual(icm_status.slugs_in(out, pats), {"p1-le-benchmark"})


class PromoteTests(unittest.TestCase):
    def make_repo(self, d, new_key):
        r = Path(d)
        for sub in ["03-verify/output", "04-synthesize/output", "05-promote/output"]:
            (r / promote.LR / sub).mkdir(parents=True)
        (r / promote.BIB).mkdir(parents=True)
        (r / "_config/catalogs").mkdir(parents=True)
        (r / promote.BIB / "references.bib").write_text("@article{Old2000,\n  doi = {10.1/old}\n}\n")
        (r / promote.LR / "03-verify/output/t-verified.bib").write_text(f"@article{{{new_key},\n  doi = {{10.1/new}}\n}}\n")
        (r / promote.LR / "03-verify/output/t-verify-report.md").write_text("## UNVERIFIED\n\n- Lost1999: no record\n")
        (r / "_config/catalogs/x.md").write_text("| L1 | m | New2001? | M2 | planned |\n| L2 | m | Old2000? | M5 | planned |\n")
        return r

    def test_plan_resolves_keys_and_collects_unverified(self):
        with tempfile.TemporaryDirectory() as d:
            p = promote.plan(self.make_repo(d, "New2001"), "t")
            self.assertEqual(p["problems"], [])
            self.assertEqual(p["flips"], {"_config/catalogs/x.md": ["New2001"]})
            self.assertEqual(len(p["unverified"]), 1)

    def test_duplicate_title_blocks(self):
        with tempfile.TemporaryDirectory() as d:
            r = self.make_repo(d, "New2001")
            (r / promote.BIB / "references.bib").write_text("@article{Old2000,\n  title = {Same Paper},\n  doi = {10.1/old}\n}\n")
            (r / promote.LR / "03-verify/output/t-verified.bib").write_text(
                "@article{New2001,\n  title = {Same paper.},\n  doi = {10.1/other-form}\n}\n")
            self.assertIn("duplicate title New2001 = Old2000", promote.plan(r, "t")["problems"])

    def test_duplicate_key_blocks(self):
        with tempfile.TemporaryDirectory() as d:
            p = promote.plan(self.make_repo(d, "Old2000"), "t")
            self.assertIn("duplicate key Old2000", p["problems"])


class ValidateTests(unittest.TestCase):
    def test_repo_passes_validation(self):
        fails, _ = icm_validate.validate(ROOT)
        self.assertEqual(fails, [])


if __name__ == "__main__":
    unittest.main()
