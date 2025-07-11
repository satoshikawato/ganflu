#!/usr/bin/env python
# coding: utf-8

import sys
import os
import logging
import argparse
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
    if len(sys.argv) == 1:
        parser.print_help(sys.stderr)
        sys.exit(1)
    return parser.parse_args(raw_args)

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

def create_new_feature(product_name, gff3_feature, slip=False):
    attrs = dict(attr.split("=") for attr in gff3_feature.attributes.split(";"))
    attr_dict = {"product": product_name}

    if slip:
        attr_dict["ribosomal_slippage"] = []
    location = FeatureLocation(gff3_feature.start - 1, gff3_feature.end, strand=1 if gff3_feature.strand == '+' else -1)
    return SeqFeature(location=location, type=gff3_feature.type, qualifiers=attr_dict)

def process_cds_feature(gff3_feature, features_in_seq, antigen_dict, antigen_list, slip_list):
    attrs = dict(attr.split("=") for attr in gff3_feature.attributes.split(";"))
    
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
                features_in_seq[product_name] = create_new_feature(product_name, antigen_dict[antigen][final_subtype]["gff3_feature"])
                if features_in_seq[product_name].qualifiers.get("note"):
                    pass
                else:
                    features_in_seq[product_name].qualifiers["note"] = []                
                features_in_seq[product_name].qualifiers["note"].append(f"subtype: {final_subtype}")
            else:
                features_in_seq[product_name] = create_new_feature(product_name, gff3_feature)
                if features_in_seq[product_name].qualifiers.get("note"):
                    pass
                else:
                    features_in_seq[product_name].qualifiers["note"] = []
                features_in_seq[product_name].qualifiers["note"].append(f"subtype: {subtype}")
    elif any(slip_gene in product_name for slip_gene in slip_list):

        if product_name in features_in_seq:
            features_in_seq[product_name] = update_existing_feature(product_name, gff3_feature, features_in_seq[product_name], slip=True)
        else:
            features_in_seq[product_name] = create_new_feature(product_name, gff3_feature, slip=True)
            if features_in_seq[product_name].qualifiers.get("note"):
                pass
            else:
                features_in_seq[product_name].qualifiers["note"] = []
    else:
        if product_name in features_in_seq:
            features_in_seq[product_name] = update_existing_feature(product_name, gff3_feature, features_in_seq[product_name])
        else:
            features_in_seq[product_name] = create_new_feature(product_name, gff3_feature)
            if features_in_seq[product_name].qualifiers.get("note"):
                pass
            else:
                features_in_seq[product_name].qualifiers["note"] = []
            
def to_seqfeatures(gff3_features, antigen_dict, antigen_list, slip_list):
    seqfeatures = defaultdict(dict)
    for gff3_feature in gff3_features:
        if gff3_feature.type == "CDS":
            process_cds_feature(gff3_feature, seqfeatures[gff3_feature.seqid], antigen_dict, antigen_list, slip_list)
    return {key: list(seqfeatures[key].values()) for key in seqfeatures}

def add_translations(seq_record):
    for feature in seq_record.features:
        if feature.qualifiers.get("note"):
            pass
        else:
            feature.qualifiers["note"] = []
        if feature.type == "CDS":
            try:
                feature.qualifiers["translation"] = feature.translate(seq_record.seq)
            except CodonTable.TranslationError as e:
                if "Extra in frame stop codon" in str(e):
                    feature.type == "misc_feature"
                    
                    feature.qualifiers["note"].append("nonfunctional due to mutation")
                        
                elif "is not a start codon" in str(e):
                    feature.qualifiers["translation"] = feature.translate(seq_record.seq, cds=False)[:-1]
                    feature.qualifiers["note"].append("start codon not found; possibly truncated")
    return seq_record



def main(raw_args=None):
    args = parse_arguments(raw_args)
    isolate = args.isolate
    
    antigen_dict = defaultdict(dict)
    gff_features = get_gff_features(args.gff)
    config = toml.load(args.toml)
    
    seq_records = [record for record in SeqIO.parse(args.input, "fasta")]
    antigen_list = list(config["serotype"].keys())
    # To store the keys of a nested dictionary whose value for "ribosomal_slippage" is true, i.e. ribosomal slippage is present in the gene:
    slip_list = [key for key, value in antigen_dict.items() if value.get("ribosomal_slippage")]
    seq_features = to_seqfeatures(gff_features, antigen_dict, antigen_list, slip_list)
    antigen_dict = defaultdict(list)
    for key in seq_features.keys():
        for feature in seq_features[key]:
            if feature.qualifiers["product"] in antigen_list:
                if "note" in feature.qualifiers:
                    # if the list of notes contains an element with the word "subtype", then get the subtype information
                    if any("subtype" in s for s in feature.qualifiers["note"]):

                        # Get the element of the list that contains the subtype information
                        index_for_subtype = [i for i, s in enumerate(feature.qualifiers["note"]) if "subtype" in s][0]
                        subtype = feature.qualifiers["note"][index_for_subtype].split(": ")[1]
                    antigen_dict[feature.qualifiers["product"]].append(subtype)
                else:
                    antigen_dict[feature.qualifiers["product"]].append("Unknown")
    for key, value in antigen_dict.items():
        if len(value) == 1:
            antigen_dict[key] = value[0]
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
    annotations["source"] = config["annotations"]["organism"].format(isolate=isolate, subtype=serotype)
    annotations["organism"] = config["annotations"]["organism"].format(isolate=isolate, subtype=serotype)
    annotations["date"] = datetime.now().strftime("%d-%b-%Y").upper()
    try:
        out_records = []
        for record in seq_records:
            contig_id = record.id
            contig_seq = record.seq
            features_in_contig = seq_features[contig_id]
            
            description = config["segments"][contig_id.split("_")[0]]["description"].format(organism=annotations["organism"], subtype=serotype)
            out_record = SeqRecord(contig_seq, id=contig_id, description = description, name=contig_id, annotations=annotations, features=features_in_contig)
            out_record = add_translations(out_record)
            out_records.append(out_record)

        with open(args.output, 'w') as handle:  
            SeqIO.write(out_records, handle, 'genbank')

    except FileNotFoundError as e:
        logger.error(f"File not found: {e}")
    except Exception as e:
        logger.error(f"An error occurred: {e}")

if __name__ == "__main__":
    main()
