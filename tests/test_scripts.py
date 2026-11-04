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
from scripts.tag_variant_qc import (  # noqa: E402
    _extract_scalar,
    tag_variant_qc,
)


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

    def test_chry_male_missingness_and_no_site_dp_or_ac_fail(self):
        """Verify chrY uses F_MISSING_M, INFO/DP is ignored, and AC=0 remains PASS."""
        from cyvcf2 import VCF

        vcf_content = (
            "##fileformat=VCFv4.2\n"
            "##contig=<ID=chr1,length=10000>\n"
            "##contig=<ID=chrY,length=10000>\n"
            '##INFO=<ID=AC,Number=A,Type=Integer,Description="Allele count">\n'
            '##INFO=<ID=DP,Number=1,Type=Integer,Description="Total Depth">\n'
            '##INFO=<ID=F_MISSING,Number=1,Type=Float,Description="Missing">\n'
            '##INFO=<ID=F_MISSING_M,Number=1,Type=Float,Description="Male missing">\n'
            "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n"
            "chrY\t100\t.\tA\tG\t50\t.\tAC=1;F_MISSING=0.65;F_MISSING_M=0.0\n"
            "chrY\t200\t.\tA\tG\t50\t.\tAC=1;F_MISSING=0.65;F_MISSING_M=0.25\n"
            "chr1\t300\t.\tA\tG\t50\t.\tAC=0;DP=2;F_MISSING=0.05\n"
        )
        thresholds = {
            "qual": 30.0,
            "qd": 2.0,
            "mq": 40.0,
            "fs": 60.0,
            "readpos": -8.0,
            "hwe": 1e-6,
            "max_missing": 0.1,
        }
        with (
            tempfile.NamedTemporaryFile("w", suffix=".vcf", delete=False) as f_in,
            tempfile.NamedTemporaryFile("w", suffix=".vcf", delete=False) as f_out,
        ):
            f_in.write(vcf_content)
            f_in_path = f_in.name
            f_out_path = f_out.name

        try:
            tag_variant_qc(f_in_path, f_out_path, thresholds)
            vcf_obj = VCF(f_out_path)
            statuses = {rec.POS: rec.INFO.get("QC_STATUS") for rec in vcf_obj}
            vcf_obj.close()
            self.assertEqual(
                statuses,
                {
                    100: "PASS",
                    200: "FAIL_MISSING",
                    300: "PASS",
                },
            )
        finally:
            os.remove(f_in_path)
            os.remove(f_out_path)

    def test_nullable_thresholds(self):
        """Verify setting thresholds to None disables those QC checks."""
        from cyvcf2 import VCF

        vcf_content = (
            "##fileformat=VCFv4.2\n"
            "##contig=<ID=chr1,length=10000>\n"
            '##INFO=<ID=QD,Number=1,Type=Float,Description="QD">\n'
            '##INFO=<ID=MQ,Number=1,Type=Float,Description="MQ">\n'
            '##INFO=<ID=FS,Number=1,Type=Float,Description="FS">\n'
            '##INFO=<ID=ReadPosRankSum,Number=1,Type=Float,Description="ReadPos">\n'
            '##INFO=<ID=F_MISSING,Number=1,Type=Float,Description="Missing">\n'
            "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n"
            "chr1\t100\t.\tA\tG\t10\t.\tQD=0.5;MQ=20.0;FS=100.0;ReadPosRankSum=-12.0;F_MISSING=0.05\n"
        )
        thresholds = {
            "qual": None,
            "qd": None,
            "mq": None,
            "fs": None,
            "readpos": None,
            "hwe": 1e-6,
            "max_missing": 0.1,
        }
        with (
            tempfile.NamedTemporaryFile("w", suffix=".vcf", delete=False) as f_in,
            tempfile.NamedTemporaryFile("w", suffix=".vcf", delete=False) as f_out,
        ):
            f_in.write(vcf_content)
            f_in_path = f_in.name
            f_out_path = f_out.name

        try:
            tag_variant_qc(f_in_path, f_out_path, thresholds)
            vcf_obj = VCF(f_out_path)
            rec = next(vcf_obj)
            self.assertEqual(rec.INFO.get("QC_STATUS"), "PASS")
            vcf_obj.close()
        finally:
            os.remove(f_in_path)
            os.remove(f_out_path)


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
