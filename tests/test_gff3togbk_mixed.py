from collections import defaultdict

import pytest
from Bio import SeqIO

from ganflu.scripts import gff3togbk


def make_cds_feature(seqid, parent_id, product, identity):
    return gff3togbk.GFF3Feature(
        seqid,
        "miniprot",
        "CDS",
        "1",
        "9",
        ".",
        "+",
        "0",
        f"Parent={parent_id};Rank=1;Identity={identity:.4f};Target={product} 1 3",
    )


def note_values(feature):
    notes = feature.qualifiers.get("note", [])
    if isinstance(notes, str):
        return [notes]
    return notes


def test_to_seqfeatures_keeps_antigen_subtype_state_contig_local():
    features = [
        make_cds_feature("contig_h1", "MP1", "HA_H1", 0.99),
        make_cds_feature("contig_h5", "MP2", "HA_H5", 0.98),
        make_cds_feature("contig_h5", "MP3", "HA_H7", 0.50),
        make_cds_feature("contig_pb2", "MP4", "PB2", 0.99),
    ]

    seq_features = gff3togbk.to_seqfeatures(
        features,
        defaultdict(dict),
        ["HA", "NA"],
        [],
        {
            "HA": {"product": "hemagglutinin"},
            "PB2": {"product": "polymerase PB2"},
        },
    )

    h1 = seq_features["contig_h1"][0]
    h5 = seq_features["contig_h5"][0]
    pb2 = seq_features["contig_pb2"][0]

    assert "subtype: H1" in note_values(h1)
    assert "subtype: H5" in note_values(h5)
    assert not any("subtype:" in note for note in note_values(pb2))


def write_toml(path):
    path.write_text(
        "\n".join(
            [
                "[segments]",
                'HA = {id = 4, file = "HA.fa", description = "{organism} segment 4 hemagglutinin (HA) gene, complete cds"}',
                'NA = {id = 6, file = "NA.fa", description = "{organism} segment 6 neuraminidase (NA) gene, complete cds"}',
                "",
                "[serotype]",
                'HA = {segment = "HA", prefix = "H"}',
                'NA = {segment = "NA", prefix = "N"}',
                "",
                "[genes]",
                'HA = {product = "hemagglutinin", has_intron = false, ribosomal_slippage = false}',
                'NA = {product = "neuraminidase", has_intron = false, ribosomal_slippage = false}',
                "",
                "[annotations]",
                'molecule_type = "cRNA"',
                'organism = "Influenza A virus ({isolate}({subtype}))"',
                'taxonomy = ["Viruses"]',
                'data_file_division = "VRL"',
                'topology = "linear"',
            ]
        )
        + "\n",
        encoding="utf-8",
    )


@pytest.mark.parametrize(
    ("ha_products", "na_products", "expected_serotype"),
    [
        (["HA_H1", "HA_H5"], ["NA_N1", "NA_N1"], "HXN1"),
        (["HA_H1", "HA_H7"], ["NA_N1", "NA_N9"], "HXNX"),
    ],
)
def test_gff3togbk_aggregates_mixed_iav_serotype_from_ha_na_notes(
    tmp_path,
    ha_products,
    na_products,
    expected_serotype,
):
    fasta = tmp_path / "mixed.fa"
    gff3 = tmp_path / "mixed.gff3"
    toml = tmp_path / "IAV.toml"
    output = tmp_path / "mixed.gbk"
    seq = "ATGAAATAA"

    fasta.write_text(
        "\n".join(
            [
                f">ha1\n{seq}",
                f">ha2\n{seq}",
                f">na1\n{seq}",
                f">na2\n{seq}",
            ]
        )
        + "\n",
        encoding="utf-8",
    )
    gff3.write_text(
        "\n".join(
            [
                "##gff-version 3",
                f"na1\tminiprot\tCDS\t1\t9\t.\t+\t0\tParent=MP3;Identity=0.9700;Target={na_products[0]} 1 3",
                f"ha1\tminiprot\tCDS\t1\t9\t.\t+\t0\tParent=MP1;Identity=0.9900;Target={ha_products[0]} 1 3",
                f"na2\tminiprot\tCDS\t1\t9\t.\t+\t0\tParent=MP4;Identity=0.9600;Target={na_products[1]} 1 3",
                f"ha2\tminiprot\tCDS\t1\t9\t.\t+\t0\tParent=MP2;Identity=0.9800;Target={ha_products[1]} 1 3",
            ]
        )
        + "\n",
        encoding="utf-8",
    )
    write_toml(toml)

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
            "A/mixed/1/2026",
        ]
    )

    gbk_text = output.read_text(encoding="utf-8")
    assert expected_serotype in gbk_text
    records = list(SeqIO.parse(output, "genbank"))
    assert len(records) == 4
    for record in records:
        cds = next(feature for feature in record.features if feature.type == "CDS")
        gene = cds.qualifiers["gene"][0]
        notes = note_values(cds)
        assert any(note.startswith(f"subtype: {gene[0]}") for note in notes)
