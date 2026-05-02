import csv
import json
import sys
import shutil
from pathlib import Path

import pytest
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqFeature import FeatureLocation, SeqFeature
from Bio.SeqRecord import SeqRecord

import ganflu
from ganflu import ganflu as ganflu_cli
from ganflu.scripts import auto_mode, gff3togbk
from ganflu.scripts.gff3togbk import add_translations


DATA_DIR = Path(__file__).parent / "data"


def test_cli_version_uses_package_version():
    assert ganflu_cli._version() == ganflu.__version__


def test_cli_isolate_is_optional_and_defaults_to_output_prefix(monkeypatch):
    monkeypatch.setattr(sys, "argv", ["ganflu", "-i", "input.fa", "-t", "IBV"])
    args = ganflu_cli._get_args()
    assert args.isolate is None
    assert ganflu_cli.resolve_isolate(args.isolate, "/tmp/my sample") == "my_sample"


def test_cli_accepts_icv_target(monkeypatch):
    monkeypatch.setattr(
        sys,
        "argv",
        ["ganflu", "-i", "input.fa", "-t", "ICV", "--isolate", "C/Ann_Arbor/1/1950"],
    )
    args = ganflu_cli._get_args()
    assert args.target == "ICV"


def test_cli_accepts_idv_target(monkeypatch):
    monkeypatch.setattr(
        sys,
        "argv",
        ["ganflu", "-i", "input.fa", "-t", "IDV", "--isolate", "D/swine/Oklahoma/1334/2011"],
    )
    args = ganflu_cli._get_args()
    assert args.target == "IDV"


def test_cli_accepts_auto_target(monkeypatch):
    monkeypatch.setattr(
        sys,
        "argv",
        ["ganflu", "-i", "input.fa", "-t", "auto", "--isolate", "unknown"],
    )
    args = ganflu_cli._get_args()
    assert args.target == "auto"
    assert args.auto_targets == "IAV,IBV,ICV,IDV"


def test_cli_accepts_gui_command(monkeypatch):
    monkeypatch.setattr(
        sys,
        "argv",
        ["ganflu", "gui", "--host", "127.0.0.1"],
    )
    args = ganflu_cli._get_args()
    assert args.command == "gui"
    assert args.host == "127.0.0.1"
    assert args.port == 0


def test_cli_help_mentions_gui_subcommand(monkeypatch, capsys):
    monkeypatch.setattr(sys, "argv", ["ganflu"])
    with pytest.raises(SystemExit) as exc_info:
        ganflu_cli._get_args()
    assert exc_info.value.code == 1
    help_text = capsys.readouterr().err
    assert "Subcommands:" in help_text
    assert "ganflu gui" in help_text
    assert "ganflu gui --help" in help_text


def test_gui_webapp_assets_and_url_format():
    web_dir = ganflu_cli.get_webapp_dir()
    assert web_dir.joinpath("index.html").is_file()

    class DummyServer:
        server_address = ("0.0.0.0", 8765)

    assert ganflu_cli.get_server_url(DummyServer()) == "http://127.0.0.1:8765/"


def make_candidate(segment, score):
    return auto_mode.CandidateHit(
        contig_id="contig1",
        target="IAV",
        parent_id=f"{segment}_hit",
        product=segment,
        segment=segment,
        identity=1.0,
        positive=1.0,
        ref_aa_start=1,
        ref_aa_end=100,
        ref_aa_length=100,
        aa_coverage=1.0,
        query_start=1,
        query_end=300,
        query_coverage=1.0,
        strand="+",
        raw_score=score,
        normalized_score=score,
        cds_count=1,
        stop_codon_present=True,
    )


def test_auto_segment_close_score_is_accepted_not_review():
    contig = SeqRecord(Seq("ATG" * 100), id="contig1")
    call = auto_mode.classify_contig(
        contig,
        [make_candidate("PB2", 0.90), make_candidate("PB1", 0.85)],
        auto_mode.AutoThresholds(min_margin=0.10),
    )

    assert call.call == "accept"
    assert call.segment == "PB2"
    assert call.status == "complete"
    assert "segment_close_score" in call.flags


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


def test_gff3togbk_uses_gene_metadata_for_cds_qualifiers(tmp_path):
    fasta = tmp_path / "input.fa"
    gff3 = tmp_path / "input.gff3"
    toml = tmp_path / "reference.toml"
    output = tmp_path / "sample.gbk"

    fasta.write_text(">seq1\nATGAAATAG\n", encoding="utf-8")
    gff3.write_text(
        "\n".join(
            [
                "##gff-version 3",
                "seq1\tminiprot\tCDS\t1\t9\t.\t+\t0\tParent=MP000001;Identity=1.0000;Target=PB2 1 3",
            ]
        )
        + "\n",
        encoding="utf-8",
    )
    toml.write_text(
        "\n".join(
            [
                "[segments]",
                'PB2 = {id = 1, file = "PB2.fa", description = "{organism} segment 1 polymerase PB2 (PB2) gene, complete cds"}',
                "",
                "[serotype]",
                "",
                "[genes]",
                'PB2 = {product = "polymerase PB2", note = "polymerase basic protein 2", has_intron = false, ribosomal_slippage = false}',
                "",
                "[annotations]",
                'molecule_type = "cRNA"',
                'organism = "Influenza A virus ({isolate})"',
                'taxonomy = ["Viruses"]',
                'data_file_division = "VRL"',
                'topology = "linear"',
            ]
        )
        + "\n",
        encoding="utf-8",
    )

    gff3togbk.main(
        [
            "-i",
            str(fasta),
            "-g",
            str(gff3),
            "--toml",
            str(toml),
            "-o",
            str(output),
            "--isolate",
            "A/Test/1/2026",
        ]
    )

    record = SeqIO.read(output, "genbank")
    cds = next(feature for feature in record.features if feature.type == "CDS")
    assert cds.qualifiers["gene"] == ["PB2"]
    assert cds.qualifiers["product"] == ["polymerase PB2"]
    assert "polymerase basic protein 2" in cds.qualifiers["note"]

    cds_header = next(
        line
        for line in (tmp_path / "sample.cds.fna").read_text(encoding="utf-8").splitlines()
        if line.startswith(">")
    )
    assert cds_header.startswith(">sample_PB2 polymerase PB2 ")


@pytest.mark.parametrize(
    ("input_file", "target", "isolate", "stem", "expected_cds"),
    [
        ("PR8.fasta", "IAV", "A/Puerto_Rico/8/1934", "PR8", 12),
        ("B_Victoria_2_87.fa", "IBV", "B/Victoria/2/1987", "B_Victoria_2_87", 11),
        ("Ann_Arbor.fna", "ICV", "C/Ann_Arbor/1/1950", "Ann_Arbor", 9),
        ("swine_Oklahoma_1334.fna", "IDV", "D/swine/Oklahoma/1334/2011", "swine_Oklahoma_1334", 8),
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


def test_cli_auto_smoke_identifies_iav_and_writes_reports(tmp_path):
    if not shutil.which("miniprot"):
        pytest.skip("miniprot is required for CLI smoke tests")

    output_stem = tmp_path / "PR8_auto"
    argv = [
        "ganflu",
        "-i",
        str(DATA_DIR / "PR8.fasta"),
        "-o",
        str(output_stem),
        "-t",
        "auto",
        "--isolate",
        "A/Puerto_Rico/8/1934",
    ]

    original_argv = sys.argv
    try:
        sys.argv = argv
        assert ganflu_cli.main() == 0
    finally:
        sys.argv = original_argv

    for suffix in (
        ".IAV.gbk",
        ".IAV.gff3",
        ".IAV.cds.fna",
        ".IAV.faa",
        ".auto.tsv",
        ".auto.summary.json",
    ):
        output_file = tmp_path / f"PR8_auto{suffix}"
        assert output_file.is_file()
        assert output_file.stat().st_size > 0

    with open(tmp_path / "PR8_auto.auto.tsv", newline="", encoding="utf-8") as handle:
        rows = list(csv.DictReader(handle, delimiter="\t"))
    assert len(rows) == 8
    assert {row["call"] for row in rows} == {"accept"}
    assert {row["target"] for row in rows} == {"IAV"}
    assert {row["status"] for row in rows} == {"complete"}

    summary = json.loads((tmp_path / "PR8_auto.auto.summary.json").read_text(encoding="utf-8"))
    assert summary["counts"]["accepted"] == 8
    assert summary["by_target"]["IAV"]["accepted"] == 8
    assert "auto.summary_json" in summary["outputs"]
