import sys
import shutil
from pathlib import Path

import pytest
from Bio.Seq import Seq
from Bio.SeqFeature import FeatureLocation, SeqFeature
from Bio.SeqRecord import SeqRecord

import ganflu
from ganflu import ganflu as ganflu_cli
from ganflu.scripts.gff3togbk import add_translations


DATA_DIR = Path(__file__).parent / "data"


def test_cli_version_uses_package_version():
    assert ganflu_cli._version() == ganflu.__version__


def test_cli_requires_isolate(monkeypatch):
    monkeypatch.setattr(sys, "argv", ["ganflu", "-i", "input.fa", "-t", "IBV"])
    with pytest.raises(SystemExit) as exc_info:
        ganflu_cli._get_args()
    assert exc_info.value.code == 2


def test_internal_stop_marks_cds_as_misc_feature():
    record = SeqRecord(Seq("ATGTAGAAATAA"), id="internal_stop")
    feature = SeqFeature(
        FeatureLocation(0, len(record.seq), strand=1),
        type="CDS",
        qualifiers={"product": "test"},
    )
    record.features = [feature]

    add_translations(record)

    assert record.features[0].type == "misc_feature"
    assert "nonfunctional due to mutation" in record.features[0].qualifiers["note"]


@pytest.mark.parametrize(
    ("input_file", "target", "isolate", "stem", "expected_cds"),
    [
        ("PR8.fasta", "IAV", "A/Puerto_Rico/8/1934", "PR8", 12),
        ("B_Victoria_2_87.fa", "IBV", "B/Victoria/2/1987", "B_Victoria_2_87", 11),
    ],
)
def test_cli_smoke_with_packaged_reference_data(
    tmp_path, input_file, target, isolate, stem, expected_cds
):
    if not shutil.which("miniprot"):
        pytest.skip("miniprot is required for CLI smoke tests")

    output_stem = tmp_path / stem
    argv = [
        "ganflu",
        "-i",
        str(DATA_DIR / input_file),
        "-o",
        str(output_stem),
        "-t",
        target,
        "--isolate",
        isolate,
    ]

    original_argv = sys.argv
    try:
        sys.argv = argv
        assert ganflu_cli.main() == 0
    finally:
        sys.argv = original_argv

    for suffix in (".gbk", ".gff3", ".cds.fna", ".faa"):
        output_file = tmp_path / f"{stem}{suffix}"
        assert output_file.is_file()
        assert output_file.stat().st_size > 0

    assert isolate in (tmp_path / f"{stem}.gbk").read_text(encoding="utf-8")
    cds_headers = [
        line
        for line in (tmp_path / f"{stem}.cds.fna").read_text(encoding="utf-8").splitlines()
        if line.startswith(">")
    ]
    assert len(cds_headers) == expected_cds
