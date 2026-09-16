import difflib
import gzip
import re
import subprocess
from pathlib import Path


def read_lines_vcf(path: Path) -> list[str]:
    """Read lines from .bcf, vcf.gz, or .vcf without header information."""

    if not path.exists():
        raise FileNotFoundError(path)
        # path = _find_alternative_vcf_path(path)

    if path.suffix == ".bcf":
        raw = subprocess.check_output(["bcftools", "view", str(path)], text=True)
        lines = raw.splitlines()
    else:
        lines = read_lines(path)

    # ignore header lines
    return [line for line in lines if not line.startswith("##")]

    # ignore bcftools and commandline header lines
    # return [line for line in lines if not re.match(r"##(bcftools|commandline)", line)]



def _find_alternative_vcf_path(path: Path) -> Path:
    # try and find other extensions
    s = str(path)
    s = s.removesuffix(".bcf").removesuffix(".gz").removesuffix(".vcf")
    for ext in [".vcf", ".vcf.gz", ".bcf"]:
        if (alt_p := Path(s + ext)).exists():
            print("found it", alt_p)
            return alt_p

    raise FileNotFoundError(f"No .vcf like file can be found for: {path}")


def read_lines(path: Path) -> list[str]:
    if path.suffix == ".gz":
        with gzip.open(path, "rt") as f:
            return f.read().splitlines()

    return path.read_text().splitlines()


def inline_diff(a, b):
    matcher = difflib.SequenceMatcher(None, a, b)
    result = []

    for tag, i1, i2, j1, j2 in matcher.get_opcodes():
        if tag == 'equal':
            result.append(a[i1:i2])
        elif tag == 'delete':
            # Red text for deletions
            result.append(f"\033[91m{a[i1:i2]}\033[0m")
        elif tag == 'insert':
            # Green text for insertions
            result.append(f"\033[92m{b[j1:j2]}\033[0m")
        elif tag == 'replace':
            # Red for old, Green for new
            result.append(f"\033[91m{a[i1:i2]}\033[0m\033[92m{b[j1:j2]}\033[0m")

    return "".join(result)
