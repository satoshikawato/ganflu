from ganflu.scripts import gff3_prune, gff3togbk


def mrna_line(parent_id, seqid, product, identity, *, start=1, end=90, score=100, rank=1, target_end=30):
    return (
        f"{seqid}\tminiprot\tmRNA\t{start}\t{end}\t{score}\t+\t.\t"
        f"ID={parent_id};Rank={rank};Identity={identity:.4f};Positive={identity:.4f};"
        f"Target={product} 1 {target_end}"
    )


def cds_line(parent_id, seqid, product, identity, *, start=1, end=90, phase=0, target_start=1, target_end=30):
    return (
        f"{seqid}\tminiprot\tCDS\t{start}\t{end}\t.\t+\t{phase}\t"
        f"Parent={parent_id};Rank=1;Identity={identity:.4f};Target={product} {target_start} {target_end}"
    )


def write_gff3(tmp_path, lines):
    path = tmp_path / "raw.gff3"
    path.write_text("##gff-version 3\n" + "\n".join(lines) + "\n", encoding="utf-8")
    return path


def feature_rows(path, feature_type):
    rows = []
    for line in path.read_text(encoding="utf-8").splitlines():
        if line.startswith("#"):
            continue
        columns = line.split("\t")
        if len(columns) == 9 and columns[2] == feature_type:
            rows.append(columns)
    return rows


def target_products(path, feature_type="mRNA"):
    products = []
    for columns in feature_rows(path, feature_type):
        attrs = gff3_prune.parse_gff3_attributes(columns[8])
        products.append(attrs["Target"].split()[0])
    return products


def test_prune_keeps_best_ha_subtype_candidate(tmp_path):
    raw_gff3 = write_gff3(
        tmp_path,
        [
            mrna_line("MP1", "contig_ha", "HA_H1", 0.80),
            cds_line("MP1", "contig_ha", "HA_H1", 0.80),
            mrna_line("MP2", "contig_ha", "HA_H5", 0.95),
            cds_line("MP2", "contig_ha", "HA_H5", 0.95),
        ],
    )
    output = tmp_path / "pruned.gff3"

    result = gff3_prune.prune_gff3(
        raw_gff3,
        output,
        protein_lengths={"HA_H1": 30, "HA_H5": 30},
        antigen_names={"HA"},
    )

    assert result.selected_parent_count == 1
    assert target_products(output) == ["HA_H5"]


def test_prune_keeps_best_na_subtype_candidate(tmp_path):
    raw_gff3 = write_gff3(
        tmp_path,
        [
            mrna_line("MP1", "contig_na", "NA_N2", 0.70),
            cds_line("MP1", "contig_na", "NA_N2", 0.70),
            mrna_line("MP2", "contig_na", "NA_N1", 0.96),
            cds_line("MP2", "contig_na", "NA_N1", 0.96),
        ],
    )
    output = tmp_path / "pruned.gff3"

    result = gff3_prune.prune_gff3(
        raw_gff3,
        output,
        protein_lengths={"NA_N1": 30, "NA_N2": 30},
        antigen_names={"NA"},
    )

    assert result.selected_parent_count == 1
    assert target_products(output) == ["NA_N1"]


def test_prune_deduplicates_same_product_same_coordinate_parent(tmp_path):
    raw_gff3 = write_gff3(
        tmp_path,
        [
            mrna_line("MP1", "contig_pb2", "PB2", 0.99, start=5, end=95),
            cds_line("MP1", "contig_pb2", "PB2", 0.99, start=5, end=95),
            mrna_line("MP2", "contig_pb2", "PB2", 0.99, start=5, end=95),
            cds_line("MP2", "contig_pb2", "PB2", 0.99, start=5, end=95),
        ],
    )
    output = tmp_path / "pruned.gff3"

    result = gff3_prune.prune_gff3(raw_gff3, output, protein_lengths={"PB2": 30})

    assert result.selected_parent_count == 1
    assert target_products(output) == ["PB2"]
    assert len(feature_rows(output, "CDS")) == 1


def test_prune_retains_distinct_products_spliced_children_and_fragments(tmp_path):
    raw_gff3 = write_gff3(
        tmp_path,
        [
            mrna_line("MP1", "contig_m", "M1", 0.99),
            cds_line("MP1", "contig_m", "M1", 0.99),
            mrna_line("MP2", "contig_m", "M2", 0.98),
            cds_line("MP2", "contig_m", "M2", 0.98),
            mrna_line("MP3", "contig_ns", "NS1", 0.99),
            cds_line("MP3", "contig_ns", "NS1", 0.99),
            mrna_line("MP4", "contig_ns", "NS2", 0.98, end=150, target_end=50),
            cds_line("MP4", "contig_ns", "NS2", 0.98, start=1, end=30, target_start=1, target_end=10),
            cds_line("MP4", "contig_ns", "NS2", 0.98, start=100, end=150, target_start=11, target_end=50),
            mrna_line("MP5", "contig_pa", "PA-X_fragment01", 0.99, target_end=20),
            cds_line("MP5", "contig_pa", "PA-X_fragment01", 0.99, target_end=20),
            mrna_line("MP6", "contig_pa", "PA-X_fragment02", 0.99, target_end=20),
            cds_line("MP6", "contig_pa", "PA-X_fragment02", 0.99, target_end=20),
        ],
    )
    output = tmp_path / "pruned.gff3"

    result = gff3_prune.prune_gff3(
        raw_gff3,
        output,
        protein_lengths={
            "M1": 30,
            "M2": 30,
            "NS1": 30,
            "NS2": 50,
            "PA-X_fragment01": 20,
            "PA-X_fragment02": 20,
        },
    )

    assert result.selected_parent_count == 6
    assert set(target_products(output)) == {"M1", "M2", "NS1", "NS2", "PA-X_fragment01", "PA-X_fragment02"}
    ns2_cds_rows = [
        columns
        for columns in feature_rows(output, "CDS")
        if gff3_prune.parse_gff3_attributes(columns[8])["Target"].startswith("NS2 ")
    ]
    assert len(ns2_cds_rows) == 2


def test_prune_applies_accepted_segment_filter_before_selection(tmp_path):
    raw_gff3 = write_gff3(
        tmp_path,
        [
            mrna_line("MP1", "contig1", "PB2", 0.99),
            cds_line("MP1", "contig1", "PB2", 0.99),
            mrna_line("MP2", "contig1", "HA_H1", 0.98),
            cds_line("MP2", "contig1", "HA_H1", 0.98),
        ],
    )
    output = tmp_path / "pruned.gff3"

    result = gff3_prune.prune_gff3(
        raw_gff3,
        output,
        protein_lengths={"PB2": 30, "HA_H1": 30},
        antigen_names={"HA"},
        accepted_segments={"contig1": "HA"},
        segment_keys=["PB2", "HA"],
        product_to_segment=gff3togbk.product_to_segment,
    )

    assert result.selected_parent_count == 1
    assert target_products(output) == ["HA_H1"]
