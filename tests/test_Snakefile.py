import os
import unittest
from pathlib import Path

from . import utils

DIR_1 = Path(os.environ.get("DIR_1", "tests/Snakefile_results"))
DIR_2 = Path(os.environ.get("DIR_2", "results/tests/Snakefile/default"))

def sort_info_part(l: str) -> str:
    """Dirty way to sort INFO"""
    if l.startswith("#C"):
        return l

    parts = l.split("\t")
    info_part = parts[7]
    sorted_info_part = ";".join(sorted(x for x in info_part.split(";")))
    parts[7] = sorted_info_part
    return "\t".join(parts)

SORT_INFO = True

class BaseTestCase(unittest.TestCase):

    def assert_vcf_equal(self, file1: Path, file2: Path):
        lines1, lines2 = utils.read_lines_vcf(file1), utils.read_lines_vcf(file2)
        for l1, l2 in zip(lines1, lines2):
            if SORT_INFO:
                l1, l2 = sort_info_part(l1), sort_info_part(l2)
            if l1 != l2:
                self.fail(utils.inline_diff(l1, l2))

    def assert_text_equal(self, file1: Path, file2: Path):
        lines1, lines2 = utils.read_lines(file1), utils.read_lines(file2)
        for l1, l2 in zip(lines1, lines2):
            if l1 != l2:
                self.fail(utils.inline_diff(l1, l2))


class TestSnakefileOutput(BaseTestCase):

    def test_pgx_pilot_full_sample_data(self):
        filename = "pgx_pilot.full_sample_data.vcf.gz"
        self.assert_vcf_equal(DIR_1 / filename, DIR_2 / filename)

    def test_pgx_pilot_sites(self):
        filename = "pgx_pilot.sites.all.vcf.gz"
        self.assert_vcf_equal(DIR_1 / filename, DIR_2 / filename)

        filename = "pgx_pilot.sites.pass.vcf.gz"
        self.assert_vcf_equal(DIR_1 / filename, DIR_2 / filename)


class TestSnakefileTempOutput(BaseTestCase):

    def test_groups_final(self):
        filename = "temp/groups_final.txt"
        self.assert_text_equal(DIR_1 / filename, DIR_2 / filename)

    def test_groups_raw(self):
        filename = "temp/groups_raw.txt"
        self.assert_text_equal(DIR_1 / filename, DIR_2 / filename)

    def test_sex_map(self):
        filename = "temp/sex_map.txt"
        self.assert_text_equal(DIR_1 / filename, DIR_2 / filename)

    @unittest.skipIf(DIR_2.name == "stream", "Stream mode skips intermediate temp VCFs")
    def test_intermediate(self):
        intermediate_files = [
            "temp/01_selected.vcf",
            "temp/02_normalized.vcf",
            "temp/03_raw_stats.vcf",
            "temp/04_masked.vcf",
            "temp/05_ploidy_fixed.vcf",
            "temp/06_final_stats.vcf",
        ]
        for filename in intermediate_files:
            with self.subTest(filename=filename):
                self.assert_vcf_equal(DIR_1 / filename, DIR_2 / filename)

