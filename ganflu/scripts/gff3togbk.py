#!/usr/bin/env python
# coding: utf-8

import sys
import os
import logging
import argparse
try:
    import tomllib
except ModuleNotFoundError:  # pragma: no cover - Python 3.10 fallback
    tomllib = None
    import toml
from collections import defaultdict
from datetime import datetime
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio.Data import CodonTable
from Bio.SeqFeature import SeqFeature, FeatureLocation, CompoundLocation

logger = logging.getLogger()
handler = logging.StreamHandler(sys.stdout)
logging.basicConfig(level=logging.INFO, format='%(asctime)s - %(levelname)s - %(message)s')


def parse_arguments(raw_args=None):
    parser = argparse.ArgumentParser(description="")
    parser.add_argument("-i", "--input", required=True, help="Input FASTA file")
    parser.add_argument("-g", "--gff", required=True, help="Input GFF3 file")
    parser.add_argument("-t", "--toml", required=True, help="Input TOML file")
    parser.add_argument("-o", "--output", required=True, help="Output GenBank file")
    parser.add_argument("-s", "--isolate", required=True, help="Isolate name")
    parser.add_argument("--preserve_original_id", "--preserve-original-id", dest="preserve_original_id", action="store_true", help="Preserve original FASTA record IDs in GenBank output")
    parser.add_argument("--cds-fna", dest="cds_fna", default=None, help="Output CDS nucleotide FASTA file (default: <output stem>.cds.fna)")
    parser.add_argument("--faa", dest="faa", default=None, help="Output amino acid FASTA file (default: <output stem>.faa)")
    if raw_args is None and len(sys.argv) == 1:
        parser.print_help(sys.stderr)
        sys.exit(1)
    return parser.parse_args(raw_args)


def load_toml_file(path):
    if tomllib is not None:
        with open(path, "rb") as handle:
            return tomllib.load(handle)
    with open(path, "r", encoding="utf-8") as handle:
        return toml.load(handle)

class GFF3Feature:
    def __init__(self, seqid, source, type_, start, end, score, strand, phase, attributes):
        self.seqid = seqid
        self.source = source
        self.type = type_
        self.start = int(start)
        self.end = int(end)
        self.score = None if score == '.' else float(score)
        self.strand = strand
        self.phase = None if phase == '.' else int(phase)
        self.attributes = attributes
        self.children = []  
        
    def _parse_attributes(self, attributes):
        attr_dict = {}
        if attributes != '.':
            attributes = attributes.split(";")
            for attribute in attributes:
                
                key, value = attribute.split('=')
                if key == "Target":
                    value = value.split(' ')[0]
                attr_dict[key] = value

        return attr_dict

    def to_seqfeature(self):
        if len(self.children) >1:
            locations = [FeatureLocation(child.start - 1, child.end, strand=1 if child.strand == '+' else -1) for child in self.children]
            location = CompoundLocation(locations)
        else:
            # Single-location feature
            location = FeatureLocation(self.start - 1, self.end, strand=1 if self.strand == '+' else -1)
        qualifiers = '' #self.attributes.copy()
        return SeqFeature(location=location, type=self.type, qualifiers=self.attributes)
    

class GFF3Parser:
    def __init__(self, file_path):
        self.file_path = file_path
        self.features = []

    def parse(self):
        feature_dict = {}
        with open(self.file_path, 'r') as file:
            for line in file:
                if line.startswith('#'):
                    continue

                columns = line.strip().split('\t')
                if len(columns) != 9:
                    raise ValueError('GFF3 line does not have 9 columns')

                feature = GFF3Feature(*columns)
                attrs = feature.attributes.split(";")
                attr_dict = {}
                if feature.type == "CDS":
                    if len(attrs) > 1:
                        for attr in attrs:
                            key, value = attr.split("=")
                            if key == "Target":
                                product = value.split(" ")[0]
                                attr_dict["product"] = product
                            elif key == "Identity":
                                attr_dict["identity"] = value
                        feature_dict[product] = GFF3Feature(feature.seqid, feature.source, feature.type, feature.start, feature.end, feature.score, feature.strand, feature.phase, attr_dict)
            
            for key in feature_dict.keys():
                self.features.append(feature_dict[key])

    def get_features(self):
        
        return self.features

    def get_feature_dict(self):
        return self.feature_dict

    def to_seqfeatures(self):
        seqfeatures = []
        for feature in self.features:
            seqfeature = feature.to_seqfeature()
            seqfeatures.append(seqfeature)
        return seqfeatures

def get_gff_features(file_path):
    features = []
    with open(file_path, 'r') as file:    
        for line in file:
            if line.startswith('#'):
                continue
            columns = line.strip().split('\t')
            if len(columns) != 9:
                raise ValueError('GFF3 line does not have 9 columns')
            feature = GFF3Feature(*columns)
            features.append(feature)
    return features


def parse_gff3_attributes(attributes):
    if attributes == '.':
        return {}
    return dict(attr.split("=", 1) for attr in attributes.split(";") if attr)


def qualifier_values(value):
    if not value:
        return []
    if isinstance(value, list):
        return value
    return [value]


def ensure_note_list(feature):
    notes = feature.qualifiers.get("note")
    if notes is None:
        feature.qualifiers["note"] = []
    elif not isinstance(notes, list):
        feature.qualifiers["note"] = [notes]
    return feature.qualifiers["note"]


def build_cds_qualifiers(gene_name, gene_configs=None, slip=False):
    gene_configs = gene_configs or {}
    gene_config = gene_configs.get(gene_name, {})
    qualifiers = {
        "gene": gene_name,
        "product": gene_config.get("product", gene_name),
    }

    note = gene_config.get("note")
    if note:
        qualifiers["note"] = qualifier_values(note)

    if slip or gene_config.get("ribosomal_slippage"):
        qualifiers["ribosomal_slippage"] = []

    return qualifiers

def update_existing_feature(product_name, gff3_feature, existing_feature, slip=False):
    added_location = FeatureLocation(gff3_feature.start - 1, gff3_feature.end, strand=1 if gff3_feature.strand == '+' else -1)
    if isinstance(existing_feature.location, CompoundLocation):
        new_location = CompoundLocation(list(existing_feature.location.parts) + [added_location])
    else:
        new_location = CompoundLocation([existing_feature.location, added_location])
    existing_feature.location = new_location
    if slip:
        existing_feature.qualifiers["ribosomal_slippage"] = None
    return existing_feature

def create_new_feature(product_name, gff3_feature, gene_configs=None, slip=False):
    attr_dict = build_cds_qualifiers(product_name, gene_configs, slip=slip)
    location = FeatureLocation(gff3_feature.start - 1, gff3_feature.end, strand=1 if gff3_feature.strand == '+' else -1)
    return SeqFeature(location=location, type=gff3_feature.type, qualifiers=attr_dict)

def format_organism(template, isolate, subtype):
    if subtype:
        return template.format(isolate=isolate, subtype=subtype)
    return template.replace("({subtype})", "").format(isolate=isolate, subtype=subtype)

def product_to_segment(product_name, segment_keys):
    segment_aliases = {
        "PB1-F2": "PB1",
        "PA-X": "PA",
        "NB": "NA",
        "M1": "M",
        "M2": "M",
        "BM2": "M",
        "CM2": "M",
        "P42": "M",
        "NS1": "NS",
        "NS2": "NS",
    }
    if product_name in segment_keys:
        return product_name
    if segment_aliases.get(product_name) in segment_keys:
        return segment_aliases[product_name]
    for segment_key in sorted(segment_keys, key=len, reverse=True):
        suffix = product_name.removeprefix(segment_key)
        if suffix != product_name and (suffix.startswith(("-", "_")) or suffix.isdigit()):
            return segment_key
    return None

def get_segment_key(contig_id, features_in_contig, segment_keys):
    segment_keys = list(segment_keys)
    candidates = []
    predicted_names = []
    for feature in features_in_contig:
        feature_names = qualifier_values(feature.qualifiers.get("gene"))
        feature_names.extend(qualifier_values(feature.qualifiers.get("product")))
        for feature_name in feature_names:
            predicted_names.append(feature_name)
            segment_key = product_to_segment(feature_name, segment_keys)
            if segment_key and segment_key not in candidates:
                candidates.append(segment_key)

    if len(candidates) == 1:
        return candidates[0]
    if len(candidates) > 1:
        raise KeyError(f"Ambiguous segment predictions for FASTA record ID '{contig_id}': {', '.join(candidates)}")
    expected_segments = ", ".join(segment_keys)
    predicted_products = ", ".join(predicted_names) if predicted_names else "none"
    raise KeyError(f"Could not infer segment from predicted CDS products for FASTA record ID '{contig_id}'. Products: {predicted_products}. Expected one of: {expected_segments}")

def get_output_id_prefix(output_path):
    prefix = os.path.splitext(os.path.basename(output_path))[0]
    return "".join(char if char.isalnum() or char in "_-" else "_" for char in prefix)

def get_fasta_output_paths(output_path, cds_fna=None, faa=None):
    stem, _ = os.path.splitext(output_path)
    return cds_fna or f"{stem}.cds.fna", faa or f"{stem}.faa"

def format_record_id(prefix, segment_key, segment_counts, segment_seen):
    if segment_counts[segment_key] == 1:
        return f"{prefix}_{segment_key}"
    segment_seen[segment_key] += 1
    return f"{prefix}_{segment_key}_{segment_seen[segment_key]}"

def process_cds_feature(gff3_feature, features_in_seq, antigen_dict, antigen_list, slip_list, gene_configs=None):
    attrs = parse_gff3_attributes(gff3_feature.attributes)
    
    product = attrs["Target"].split(" ")[0]
    product_name = product.split("_")[0]


    if any(antigen in product_name for antigen in antigen_list):
        
        if len(product.split("_")) > 1:
            antigen = product_name

            subtype = product.split("_")[1]
            identity = float(attrs["Identity"])
            
            location = FeatureLocation(gff3_feature.start - 1, gff3_feature.end, strand=1 if gff3_feature.strand == '+' else -1)
            antigen_dict[antigen][subtype] = {"subtype": subtype, "identity": identity, "location": location, "gff3_feature": gff3_feature}
            if product_name in features_in_seq:
                final_subtype = max(antigen_dict[antigen], key=lambda x: antigen_dict[antigen][x]["identity"])
                features_in_seq[product_name] = create_new_feature(product_name, antigen_dict[antigen][final_subtype]["gff3_feature"], gene_configs)
                ensure_note_list(features_in_seq[product_name]).append(f"subtype: {final_subtype}")
            else:
                features_in_seq[product_name] = create_new_feature(product_name, gff3_feature, gene_configs)
                ensure_note_list(features_in_seq[product_name]).append(f"subtype: {subtype}")
        else:
            if product_name in features_in_seq:
                features_in_seq[product_name] = update_existing_feature(product_name, gff3_feature, features_in_seq[product_name])
            else:
                features_in_seq[product_name] = create_new_feature(product_name, gff3_feature, gene_configs)
                ensure_note_list(features_in_seq[product_name])
    elif any(slip_gene in product_name for slip_gene in slip_list):

        if product_name in features_in_seq:
            features_in_seq[product_name] = update_existing_feature(product_name, gff3_feature, features_in_seq[product_name], slip=True)
        else:
            features_in_seq[product_name] = create_new_feature(product_name, gff3_feature, gene_configs, slip=True)
            ensure_note_list(features_in_seq[product_name])
    else:
        if product_name in features_in_seq:
            features_in_seq[product_name] = update_existing_feature(product_name, gff3_feature, features_in_seq[product_name])
        else:
            features_in_seq[product_name] = create_new_feature(product_name, gff3_feature, gene_configs)
            ensure_note_list(features_in_seq[product_name])
            
def to_seqfeatures(gff3_features, antigen_dict, antigen_list, slip_list, gene_configs=None):
    seqfeatures = defaultdict(dict)
    for gff3_feature in gff3_features:
        if gff3_feature.type == "CDS":
            process_cds_feature(gff3_feature, seqfeatures[gff3_feature.seqid], antigen_dict, antigen_list, slip_list, gene_configs)
    return {key: list(seqfeatures[key].values()) for key in seqfeatures}

def add_translations(seq_record):
    for feature in seq_record.features:
        ensure_note_list(feature)
        if feature.type == "CDS":
            try:
                feature.qualifiers["translation"] = feature.translate(seq_record.seq)
            except CodonTable.TranslationError as e:
                if "Extra in frame stop codon" in str(e):
                    feature.type = "misc_feature"
                    
                    feature.qualifiers["note"].append("nonfunctional due to mutation")
                        
                elif "is not a start codon" in str(e):
                    feature.qualifiers["translation"] = feature.translate(seq_record.seq, cds=False)[:-1]
                    feature.qualifiers["note"].append("start codon not found; possibly truncated")
                elif "Final codon" in str(e) and "is not a stop codon" in str(e):
                    feature.qualifiers["translation"] = feature.translate(seq_record.seq, cds=False)
                    feature.qualifiers["note"].append("stop codon not found; possibly truncated")
    return seq_record

def get_first_qualifier(feature, key, default=""):
    value = feature.qualifiers.get(key, default)
    if isinstance(value, list):
        return value[0] if value else default
    return value

def sanitize_fasta_id(value):
    value = str(value) if value else "unknown"
    return "".join(char if char.isalnum() or char in "._-" else "_" for char in value)

def format_cds_record_id(seq_record_id, feature, seen_ids):
    product = get_first_qualifier(feature, "gene") or get_first_qualifier(feature, "product", "CDS")
    seq_record_id = sanitize_fasta_id(seq_record_id)
    product_id = sanitize_fasta_id(product)
    if seq_record_id == product_id or seq_record_id.endswith(f"_{product_id}"):
        base_id = seq_record_id
    else:
        base_id = f"{seq_record_id}_{product_id}"
    seen_ids[base_id] += 1
    if seen_ids[base_id] == 1:
        return base_id
    return f"{base_id}_{seen_ids[base_id]}"

def get_feature_translation(feature, seq_record):
    translation = get_first_qualifier(feature, "translation")
    if translation:
        return Seq(str(translation))
    try:
        translation = feature.translate(seq_record.seq, cds=False)
    except CodonTable.TranslationError:
        return None
    if str(translation).endswith("*"):
        translation = translation[:-1]
    return Seq(str(translation))

def build_cds_fasta_records(seq_records):
    cds_records = []
    aa_records = []
    seen_ids = defaultdict(int)
    for seq_record in seq_records:
        for feature in seq_record.features:
            if feature.type != "CDS":
                continue
            cds_record_id = format_cds_record_id(seq_record.id, feature, seen_ids)
            product = get_first_qualifier(feature, "product", "CDS")
            location = str(feature.location)
            description = f"{product} {seq_record.id}:{location}"

            cds_seq = feature.extract(seq_record.seq)
            cds_records.append(SeqRecord(cds_seq, id=cds_record_id, name=cds_record_id, description=description))

            translation = get_feature_translation(feature, seq_record)
            if translation:
                aa_records.append(SeqRecord(translation, id=cds_record_id, name=cds_record_id, description=description))
    return cds_records, aa_records

def write_cds_fasta_files(seq_records, cds_fna_path, faa_path):
    cds_records, aa_records = build_cds_fasta_records(seq_records)
    with open(cds_fna_path, "w") as handle:
        SeqIO.write(cds_records, handle, "fasta")
    with open(faa_path, "w") as handle:
        SeqIO.write(aa_records, handle, "fasta")
    return len(cds_records), len(aa_records)



def main(raw_args=None):
    args = parse_arguments(raw_args)
    isolate = args.isolate
    
    antigen_dict = defaultdict(dict)
    gff_features = get_gff_features(args.gff)
    config = load_toml_file(args.toml)
    
    seq_records = [record for record in SeqIO.parse(args.input, "fasta")]
    antigen_list = list(config.get("serotype", {}).keys())
    # To store the keys of a nested dictionary whose value for "ribosomal_slippage" is true, i.e. ribosomal slippage is present in the gene:
    gene_configs = config.get("genes", {})
    slip_list = [key for key, value in gene_configs.items() if value.get("ribosomal_slippage")]
    seq_features = to_seqfeatures(gff_features, antigen_dict, antigen_list, slip_list, gene_configs)
    antigen_dict = defaultdict(list)
    for key in seq_features.keys():
        for feature in seq_features[key]:
            gene_name = get_first_qualifier(feature, "gene")
            if gene_name in antigen_list:
                subtype = ""
                if "note" in feature.qualifiers:
                    # if the list of notes contains an element with the word "subtype", then get the subtype information
                    if any("subtype" in s for s in feature.qualifiers["note"]):

                        # Get the element of the list that contains the subtype information
                        index_for_subtype = [i for i, s in enumerate(feature.qualifiers["note"]) if "subtype" in s][0]
                        subtype = feature.qualifiers["note"][index_for_subtype].split(": ")[1]
                antigen_dict[gene_name].append(subtype)
    for key, value in antigen_dict.items():
        subtypes = [subtype for subtype in value if subtype]
        unique_subtypes = list(dict.fromkeys(subtypes))
        if len(unique_subtypes) == 0:
            antigen_dict[key] = ""
        elif len(unique_subtypes) == 1:
            antigen_dict[key] = unique_subtypes[0]
        else:
            antigen_dict[key] = f"{key[0]}X"
    serotype = "".join([value for _, value in antigen_dict.items() if value])
    annotations = {}
    annotations["molecule_type"] = config["annotations"]["molecule_type"]
    annotations["isolate"] = args.isolate
    annotations["topology"] = config["annotations"]["topology"]
    annotations["taxonomy"] = config["annotations"]["taxonomy"]
    annotations["data_file_division"] = config["annotations"]["data_file_division"]
    annotations["serotype"] = serotype
    annotations["source"] = format_organism(config["annotations"]["organism"], isolate, serotype)
    annotations["organism"] = annotations["source"]
    annotations["date"] = datetime.now().strftime("%d-%b-%Y").upper()
    try:
        out_records = []
        record_segments = []
        for record in seq_records:
            contig_id = record.id
            features_in_contig = seq_features.get(contig_id, [])
            segment_key = get_segment_key(contig_id, features_in_contig, config["segments"].keys())
            record_segments.append((record, features_in_contig, segment_key))

        segment_counts = defaultdict(int)
        for _, _, segment_key in record_segments:
            segment_counts[segment_key] += 1
        segment_seen = defaultdict(int)
        id_prefix = get_output_id_prefix(args.output)
        for record, features_in_contig, segment_key in record_segments:
            contig_id = record.id
            contig_seq = record.seq
            record_id = contig_id if args.preserve_original_id else format_record_id(id_prefix, segment_key, segment_counts, segment_seen)
            
            description = config["segments"][segment_key]["description"].format(organism=annotations["organism"], subtype=serotype)
            out_record = SeqRecord(contig_seq, id=record_id, description = description, name=record_id, annotations=annotations, features=features_in_contig)
            out_record = add_translations(out_record)
            out_records.append(out_record)

        with open(args.output, 'w') as handle:  
            SeqIO.write(out_records, handle, 'genbank')
        cds_fna_path, faa_path = get_fasta_output_paths(args.output, args.cds_fna, args.faa)
        cds_count, aa_count = write_cds_fasta_files(out_records, cds_fna_path, faa_path)
        logger.info(f"CDS nucleotide FASTA output: {cds_fna_path} ({cds_count} records)")
        logger.info(f"Amino acid FASTA output: {faa_path} ({aa_count} records)")

    except FileNotFoundError as e:
        logger.error(f"File not found: {e}")
        raise
    except Exception as e:
        logger.error(f"An error occurred: {e}")
        raise

if __name__ == "__main__":
    main()
