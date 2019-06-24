#!/usr/bin/env python3
# coding: utf8

"""
Description:
Using two fasta files (genome and proteome) and the functional GFF from an IMG run create a Genbank file.

Usage:
gbk_creator_from_gff.py -fg <Genome fasta file> -fp <Protein Fasta file> -g <GFF file> -s <Species name> -o <GBK Output file name>
"""

import argparse
import datetime
import gffutils
import numpy as np
import os
import pandas as pa
import re
import requests
import sys

from Bio import SeqFeature as sf
from Bio import SeqIO
from Bio.Alphabet import IUPAC
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from collections import OrderedDict


def create_taxonomic_data(species_name):
    """
    Query the EBI with the species name to create a dictionary containing taxon id,
    taxonomy and some other informations.
    """
    species_informations = {}
    species_name_url = species_name.replace(' ', '%20')

    url = 'https://www.ebi.ac.uk/ena/data/taxonomy/v1/taxon/scientific-name/' + species_name_url
    response = requests.get(url)
    temp_species_informations = response.json()[0]

    for temp_species_information in temp_species_informations:
        if temp_species_information == 'lineage':
            species_informations['taxonomy'] = temp_species_informations[temp_species_information].split('; ')[:-1]
        elif temp_species_information == 'division':
            species_informations['data_file_division'] = temp_species_informations[temp_species_information]
        elif temp_species_information == 'taxId':
            species_informations['db_xref'] = 'taxon:' + str(temp_species_informations[temp_species_information])
        else:
            species_informations[temp_species_information] = temp_species_informations[temp_species_information]

    compatible_species_name = species_name.replace('/', '_')
    species_informations['description'] = compatible_species_name + ' genome'
    species_informations['organism'] = compatible_species_name
    species_informations['keywords'] = [compatible_species_name]

    return species_informations


def contig_info(contig_id, contig_seq, species_informations):
    """
    Create contig information from species_informations dictionary and contig id and contig seq.
    """
    record = SeqRecord(contig_seq, id=contig_id, name=contig_id,
                    description=species_informations['description'])

    record.seq.alphabet = IUPAC.ambiguous_dna
    if 'data_file_division' in species_informations:
        record.annotations['data_file_division'] = species_informations['data_file_division']
    record.annotations['date'] = datetime.date.today().strftime('%d-%b-%Y').upper()
    if 'topology' in species_informations:
        record.annotations['topology'] = species_informations['topology']
    record.annotations['accessions'] = contig_id
    if 'organism' in species_informations:
        record.annotations['organism'] = species_informations['organism']
    # Use of literal_eval for taxonomy and keywords to retrieve list.
    if 'taxonomy' in species_informations:
        record.annotations['taxonomy'] = species_informations['taxonomy']
    if 'keywords' in species_informations:
        record.annotations['keywords'] = species_informations['keywords']
    if 'source' in species_informations:
        record.annotations['source'] = species_informations['source']

    new_feature_source = sf.SeqFeature(sf.FeatureLocation(1-1,
                                                        len(contig_seq)),
                                                        type="source")
    new_feature_source.qualifiers['scaffold'] = contig_id
    if 'isolate' in species_informations:
        new_feature_source.qualifiers['isolate'] = species_informations['isolate']
    # db_xref corresponds to the taxon NCBI ID.
    # Important if you want to use Pathway Tools after.
    if 'db_xref' in species_informations:
        new_feature_source.qualifiers['db_xref'] = species_informations['db_xref']
    if 'cell_type' in species_informations:
        new_feature_source.qualifiers['cell_type'] = species_informations['cell_type']
    if 'dev_stage' in species_informations:
        new_feature_source.qualifiers['dev_stage'] = species_informations['dev_stage']
    if 'mol_type' in species_informations:
        new_feature_source.qualifiers['mol_type'] = species_informations['mol_type']

    record.features.append(new_feature_source)

    return record


def strand_change(input_strand):
    """
    The input is strand in str ('-', '+') modify it to be a strand in int (-1, +1) to 
    be compatible with SeqIO strand reading.
    """
    if isinstance(input_strand, str):
        if input_strand == '-':
            new_strand = -1
        elif input_strand == '+':
            new_strand = +1
        if input_strand == '.':
            new_strand = None
        elif input_strand == '?':
            new_strand = 0
    elif isinstance(input_strand, int):
        if input_strand == -1:
            new_strand = input_strand
        elif input_strand == +1:
            new_strand = input_strand

    return new_strand


def gff_to_gbk(genome_fasta, prot_fasta, gff_file, species_name, gbk_out):
    """
    From a genome fasta (containing each contigs of the genome),
    a protein fasta (containing each protein sequence),
    an annotation table (containing gene name associated with GO terms, InterPro and EC),
    a gff file (containing gene, exon, mRNA, ncRNA, tRNA),
    a contig information table (containing species name, taxon ID, ..)
    create a genbank file.
    """

    print('Creating GFF database (gffutils)')
    # Create the gff database file.
    # gffutils use sqlite3 file-based database to access data inside GFF.
    # ':memory:' ask gffutils to keep database in memory instead of writting in a file.
    gff_database = gffutils.create_db(gff_file, ':memory:', force=True, keep_order=True, merge_strategy='merge', sort_attribute_values=True)

    # Length of your gene ID.
    # Catch it in the GFF database.
    # It's pretty dumb as we go into a loop for one information.
    # But I don't find another way to catch the length of gene_id.
    length_gene_id = 0

    for gene in gff_database.features_of_type('gene'):
        length_gene_id = len(gene.id.replace('gene:', ''))
        break

    # Get the longest contig ID to check if all contig IDs have the
    # same length, if not add 0 (at the supposed position of the number).
    longest_contig_id = ""

    for contig_for_length_id in gff_database.features_of_type('sequence_assembly'):
        if len(longest_contig_id) < len(contig_for_length_id.id):
            longest_contig_id = contig_for_length_id.id

    print('Formatting fasta and annotation file')
    # Dictionary with scaffold/chromosome id as key and sequence as value.
    contig_seqs = OrderedDict()

    for record in SeqIO.parse(genome_fasta, "fasta"):
        id_contig = record.id
        contig_seqs[id_contig] = record.seq


    # Dictionary with gene id as key and protein sequence as value.
    gene_protein_seq = {}

    for record in SeqIO.parse(prot_fasta, "fasta"):
        gene_protein_seq[record.id] = record.seq

    # Create a taxonomy dictionary querying the EBI.
    species_informations = create_taxonomic_data(species_name)

    # All SeqRecord objects will be stored in a list and then give to the SeqIO writer to create the genbank.
    seq_objects = []

    print('Assembling Genbank informations')

    # Iterate through each contig.
    # Then iterate through gene and throug RNA linked with the gene.
    # Then look if protein informations are available.

    contig_number = 1
    for contig_id in sorted(contig_seqs):
        # Data for each contig.
        record = contig_info(contig_id, contig_seqs[contig_id], species_informations)
        feature_number = 1
        for feature in gff_database.featuretypes():
            if feature == 'CDS':
                for gene in gff_database.features_of_type('CDS'):
                    feature_label = species_name[:4] + '_' + str(contig_number).zfill(len(str(len(contig_seqs)))) + '_' +  str(feature_number).zfill(len(str(gff_database.count_features_of_type())))
                    gene_contig = gene.chrom
                    if gene_contig == contig_id:
                        id_gene = gene.id
                        start_position = gene.start -1
                        end_position = gene.end
                        strand = strand_change(gene.strand)
                        new_feature_gene = sf.SeqFeature(sf.FeatureLocation(start_position,
                                                                            end_position,
                                                                            strand),
                                                                            type="gene")
                        new_feature_gene.qualifiers['locus_tag'] = id_gene
                        new_feature_gene.qualifiers['gene'] = feature_label
                        # Add gene information to contig record.
                        record.features.append(new_feature_gene)

                        # Create mRNA data.
                        new_feature_mRNA = sf.SeqFeature(sf.FeatureLocation(start_position,
                                                                            end_position,
                                                                            strand),
                                                                            type="mRNA")
                        new_feature_mRNA.qualifiers['locus_tag'] = feature_label
                        new_feature_mRNA.qualifiers['gene'] = feature_label
                        # Add gene information to contig record.
                        record.features.append(new_feature_mRNA)

                        # Create CDS and get annotations.
                        new_feature_cds = sf.SeqFeature(sf.FeatureLocation(start_position,
                                                                            end_position,
                                                                            strand),
                                                                        type="CDS")

                        new_feature_cds.qualifiers['translation'] = gene_protein_seq[id_gene]
                        new_feature_cds.qualifiers['locus_tag'] = feature_label
                        new_feature_cds.qualifiers['gene'] = feature_label

                        if 'product' in gene.attributes:
                            new_feature_cds.qualifiers['note'] = gene.attributes['product']
                        if 'ec_number' in gene.attributes:
                            new_feature_cds.qualifiers['EC_number'] = [ec.replace('EC:', '') for ec in gene.attributes['ec_number']]

                        # Add CDS information to contig record
                        record.features.append(new_feature_cds)
                        feature_number += 1
            else:
                for data_feature in gff_database.features_of_type(feature):
                    feature_label = species_name[:4] + '_' + str(contig_number).zfill(len(str(len(contig_seqs)))) + '_' +  str(feature_number).zfill(len(str(gff_database.count_features_of_type())))
                    feature_contig = data_feature.chrom
                    if feature_contig == contig_id:
                        id_feature = data_feature.id
                        start_position = data_feature.start -1
                        end_position = data_feature.end
                        strand = strand_change(data_feature.strand)
                    new_feature = sf.SeqFeature(sf.FeatureLocation(start_position,
                                                                        end_position,
                                                                        strand),
                                                                        type=feature)
                    new_feature.qualifiers['locus_tag'] = feature_label
                    if 'product' in data_feature.attributes:
                        new_feature.qualifiers['note'] = data_feature.attributes['product']
                    record.features.append(new_feature)
                    feature_number += 1

        seq_objects.append(record)
        contig_number += 1

    # Create Genbank with the list of SeqRecord.
    SeqIO.write(seq_objects, gbk_out, 'genbank')


def main(genome_fasta, prot_fasta, gff_file_folder, species_name, gbk_out):

    # Check if gff is a file or is multiple files in a folder.
    # If it's multiple files, it wil merge them in one.
    if os.path.isfile(gff_file_folder):
        gff_file = gff_file_folder
    else:
        sys.exit('No GFF file.')

    gff_to_gbk(genome_fasta, prot_fasta, gff_file, species_name, gbk_out)

def run():
    parser = argparse.ArgumentParser(prog = "gbk_creator_from_gff.py")
    parser.add_argument("-fg", "--fgen", dest = "genome_fasta", metavar = "FILE", help = "contig fasta file", required = True)
    parser.add_argument("-fp", "--fprot", dest = "prot_fasta", metavar = "FILE", help = "protein fasta file", required = True)
    parser.add_argument("-g", "--gff", dest = "gff_file_folder", metavar = "FILE or FOLDER", help = "gff file or folder containing multiple gff", required = True)
    parser.add_argument("-s", "--speciesname", dest = "species_name", metavar = "STRING", help = "species scientific name", required = True)
    parser.add_argument("-o", "--output", dest = "gbk_out", metavar = "FILE", help = "output file", default = "mygbk.gbk")
    args = parser.parse_args()

    main(genome_fasta=args.genome_fasta, prot_fasta=args.prot_fasta,
                gff_file_folder=args.gff_file_folder, species_name=args.species_name, gbk_out=args.gbk_out)

if __name__ == '__main__':
	run()
