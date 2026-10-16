"""Unit tests for pipeline helper scripts."""

import os
import sys
import tempfile
import unittest

import pandas as pd

# Ensure repository root is on sys.path
REPO_ROOT = os.path.abspath(os.path.join(os.path.dirname(__file__), ".."))
if REPO_ROOT not in sys.path:
    sys.path.insert(0, REPO_ROOT)

from scripts.generate_groups import generate_genomics_metadata  # noqa: E402
from scripts.generate_sex_map import generate_sex_map  # noqa: E402
from scripts.preprocess_vcf_for_pypgx import (  # noqa: E402
    preprocess_vcf_for_pypgx,
)
from scripts.tag_variant_qc import _extract_scalar  # noqa: E402


class TestMetadataParsing(unittest.TestCase):
    """Test header detection and parsing in metadata scripts."""

    def test_headerless_file_with_id_substring_in_sample(self):
        """Ensure sample IDs with 'id' or 'sample' are not taken as headers."""
        data = "patient_id_01\tM\tPT\nsample_id_02\tF\tES\nnormal_03\t1\tFR\n"
        with (
            tempfile.NamedTemporaryFile("w", delete=False) as f_in,
            tempfile.NamedTemporaryFile("w", delete=False) as f_groups,
            tempfile.NamedTemporaryFile("w", delete=False) as f_sex,
        ):
            f_in.write(data)
            in_path = f_in.name
            groups_path = f_groups.name
            sex_path = f_sex.name

        try:
            generate_genomics_metadata(in_path, groups_path)
            generate_sex_map(in_path, sex_path)

            groups_df = pd.read_csv(groups_path, sep="\t", header=None)
            sex_df = pd.read_csv(sex_path, sep="\t", header=None)

            # All 3 samples must be present, especially patient_id_01
            self.assertEqual(len(groups_df), 3)
            self.assertEqual(len(sex_df), 3)
            self.assertEqual(groups_df.iloc[0, 0], "patient_id_01")
            self.assertEqual(sex_df.iloc[0, 0], "patient_id_01")
            self.assertEqual(sex_df.iloc[0, 1], "M")
            self.assertEqual(sex_df.iloc[1, 1], "F")
            self.assertEqual(sex_df.iloc[2, 1], "M")
        finally:
            os.remove(in_path)
            os.remove(groups_path)
            os.remove(sex_path)

    def test_standard_header_file(self):
        """Ensure explicit headers are recognized and excluded from data rows."""
        data = "Sample_ID\tSex\tCountry_Code\nHG001\tFEMALE\tDE\nHG002\tMALE\tIT\n"
        with (
            tempfile.NamedTemporaryFile("w", delete=False) as f_in,
            tempfile.NamedTemporaryFile("w", delete=False) as f_groups,
            tempfile.NamedTemporaryFile("w", delete=False) as f_sex,
        ):
            f_in.write(data)
            in_path = f_in.name
            groups_path = f_groups.name
            sex_path = f_sex.name

        try:
            generate_genomics_metadata(in_path, groups_path)
            generate_sex_map(in_path, sex_path)

            groups_df = pd.read_csv(groups_path, sep="\t", header=None)
            sex_df = pd.read_csv(sex_path, sep="\t", header=None)

            self.assertEqual(len(groups_df), 2)
            self.assertEqual(len(sex_df), 2)
            self.assertEqual(sex_df.iloc[0, 0], "HG001")
            self.assertEqual(sex_df.iloc[0, 1], "F")
            self.assertEqual(sex_df.iloc[1, 0], "HG002")
            self.assertEqual(sex_df.iloc[1, 1], "M")
        finally:
            os.remove(in_path)
            os.remove(groups_path)
            os.remove(sex_path)


class TestTagVariantQcHelpers(unittest.TestCase):
    """Test helper functions in tag_variant_qc."""

    def test_extract_scalar(self):
        """Test extraction of scalar values from sequences or raw values."""
        self.assertEqual(_extract_scalar([1.5, 2.5]), 1.5)
        self.assertEqual(_extract_scalar((3.0, 4.0)), 3.0)
        self.assertEqual(_extract_scalar(42.0), 42.0)
        self.assertIsNone(_extract_scalar([]))
        self.assertIsNone(_extract_scalar(None))

    def test_streaming_cli(self):
        """Verify tag_variant_qc CLI accepts stdin and streams BCF to stdout."""
        import subprocess

        from cyvcf2 import VCF

        vcf_content = (
            "##fileformat=VCFv4.2\n"
            "##contig=<ID=chr1,length=10000>\n"
            '##INFO=<ID=DP,Number=1,Type=Integer,Description="Depth">\n'
            "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n"
            "chr1\t100\t.\tA\tG\t50\tPASS\tDP=30\n"
        )
        script = os.path.join(REPO_ROOT, "scripts", "tag_variant_qc.py")
        proc = subprocess.run(
            [sys.executable, script, "-", "-Ou"],
            input=vcf_content.encode("utf-8"),
            capture_output=True,
            check=True,
        )
        with tempfile.NamedTemporaryFile("wb", suffix=".bcf", delete=False) as f_bcf:
            f_bcf.write(proc.stdout)
            bcf_path = f_bcf.name

        try:
            vcf_obj = VCF(bcf_path)
            self.assertIn("QC_STATUS", vcf_obj)
            rec = next(vcf_obj)
            self.assertEqual(rec.INFO.get("QC_STATUS"), "PASS")
            vcf_obj.close()
        finally:
            if os.path.exists(bcf_path):
                os.remove(bcf_path)


class TestPreprocessVcf(unittest.TestCase):
    """Test VCF preprocessing for PyPGX."""

    def test_vcf_without_ad_header(self):
        """Ensure VCF files without FORMAT/AD in header are handled gracefully."""
        vcf_content = (
            "##fileformat=VCFv4.2\n"
            "##contig=<ID=chr1,length=10000>\n"
            '##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">\n'
            "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tSAMPLE1\n"
            "chr1\t100\t.\tA\tG\t.\tPASS\t.\tGT\t0/1\n"
        )
        with (
            tempfile.NamedTemporaryFile("w", suffix=".vcf", delete=False) as f_in,
            tempfile.NamedTemporaryFile("w", suffix=".vcf.gz", delete=False) as f_out,
        ):
            f_in.write(vcf_content)
            f_in_path = f_in.name
            f_out_path = f_out.name

        try:
            preprocess_vcf_for_pypgx(f_in_path, f_out_path)
            self.assertTrue(os.path.exists(f_out_path))
            self.assertGreater(os.path.getsize(f_out_path), 0)
        finally:
            if os.path.exists(f_in_path):
                os.remove(f_in_path)
            if os.path.exists(f_out_path):
                os.remove(f_out_path)


if __name__ == "__main__":
    unittest.main()
