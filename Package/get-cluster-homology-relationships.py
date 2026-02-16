#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# pylint: disable=broad-except
# pylint: disable=invalid-name
# pylint: disable=line-too-long
# pylint: disable=multiple-statements
# pylint: disable=too-many-lines
# pylint: disable=unnecessary-pass

#-------------------------------------------------------------------------------

'''
This program gets homology relationships corresponding to set of sequence clusters
yielded by process of searching for sequence homology using the comparative genomics
data of quercusTOA (Quercus Taxonomy-oriented Annotation) database.

This software has been developed by:

    GI en Desarrollo de Especies y Comunidades Leñosas (WooSp)
    Dpto. Sistemas y Recursos Naturales
    ETSI Montes, Forestal y del Medio Natural
    Universidad Politecnica de Madrid
    https://github.com/ggfhf/

Licence: GNU General Public Licence Version 3.
'''

#-------------------------------------------------------------------------------

import argparse
import gzip
import os
import sys

from collections import defaultdict

import genlib
import sqllib

#-------------------------------------------------------------------------------

def main():
    '''
    Main line of the program.
    '''

    # check the operating system
    genlib.check_os()

    # get and check the arguments
    parser = build_parser()
    args = parser.parse_args()
    check_args(args)

    # connect to the quercusTOA comparative genomics database
    conn = sqllib.connect_database(args.comparative_genomics_database)

    # attach to the quercusTOA functional annotations database
    sqllib.attach_database(conn, 'functional_annotations_database', args.functional_annotations_database)

    # get the homology relationships corresponding to set of sequence clusters yielded by the alignment process
    get_cluster_homology_relationships(conn, args.blastp_alignment_file, args.homology_relationships_file)

    # close connection to quercusTOA comparative genomics database
    conn.close()

#-------------------------------------------------------------------------------

def build_parser():
    '''
    Build the parser with the available arguments.
    '''

    # create the parser and add arguments
    description = 'Description: This program gets homology relationships corresponding to set of sequence clusters\n' \
       'yielded by process of searching for sequence homology using the comparative genomics data of\n' \
       f'{genlib.get_app_short_name()} database.'
    text = f'{genlib.get_app_short_name()} v{genlib.get_app_version()} - {os.path.basename(__file__)}\n\n{description}\n'
    usage = f'\r{text.ljust(len("usage:"))}\nUsage: {os.path.basename(__file__)} arguments'
    parser = argparse.ArgumentParser(usage=usage)
    parser._optionals.title = 'Arguments'    #pylint: disable=protected-access
    parser.add_argument('--comparative-db', dest='comparative_genomics_database', help=f'Path of the {genlib.get_app_short_name} comparative genomics database (mandatory).')
    parser.add_argument('--annotations-db', dest='functional_annotations_database', help=f'Path of the {genlib.get_app_short_name} functional annotations database (mandatory).')
    parser.add_argument('--blastp-alignments', dest='blastp_alignment_file', help='Path of the alignment file yielded by blastp (mandatory).')
    parser.add_argument('--homology', dest='homology_relationships_file', help='Path of the homology relationships file (mandatory).')
    parser.add_argument('--verbose', dest='verbose', help=f'Additional job status info during the run: {genlib.get_verbose_code_list_text()}; default: {genlib.Const.DEFAULT_VERBOSE}.')
    parser.add_argument('--trace', dest='trace', help=f'Additional info useful to the developer team: {genlib.get_trace_code_list_text()}; default: {genlib.Const.DEFAULT_TRACE}.')

    # return the paser
    return parser

#-------------------------------------------------------------------------------

def check_args(args):
    '''
    Check the input arguments.
    '''

    # initialize the control variable
    OK = True

    # check "comparative_genomics_database"
    if args.comparative_genomics_database is None:
        genlib.Message.print('error', f'*** The {genlib.get_app_short_name} comparative genomics database is not indicated in the input arguments.')
        OK = False
    elif not os.path.isfile(args.comparative_genomics_database):
        genlib.Message.print('error', f'*** The file {args.comparative_genomics_database} does not exist.')
        OK = False

    # check "functional_annotations_database"
    if args.functional_annotations_database is None:
        genlib.Message.print('error', f'*** The {genlib.get_app_short_name} functional annotations database is not indicated in the input arguments.')
        OK = False
    elif not os.path.isfile(args.functional_annotations_database):
        genlib.Message.print('error', f'*** The file {args.functional_annotations_database} does not exist.')
        OK = False

    # check "blastp_alignment_file"
    if args.blastp_alignment_file is None:
        genlib.Message.print('error', '*** The input blastp alignment file is not indicated in the input arguments.')
        OK = False
    elif not os.path.isfile(args.blastp_alignment_file):
        genlib.Message.print('error', f'*** The file {args.blastp_alignment_file} does not exist.')
        OK = False

    # check "homology_relationships_file"
    if args.homology_relationships_file is None:
        genlib.Message.print('error', '*** The homology relationships file is not indicated in the input arguments.')
        OK = False

    # check "verbose"
    if args.verbose is None:
        args.verbose = genlib.Const.DEFAULT_VERBOSE
    elif not genlib.check_code(args.verbose, genlib.get_verbose_code_list(), case_sensitive=False):
        genlib.Message.print('error', f'*** verbose has to be {genlib.get_verbose_code_list_text()}.')
        OK = False
    if args.verbose.upper() == 'Y':
        genlib.Message.set_verbose_status(True)

    # check "trace"
    if args.trace is None:
        args.trace = genlib.Const.DEFAULT_TRACE
    elif not genlib.check_code(args.trace, genlib.get_trace_code_list(), case_sensitive=False):
        genlib.Message.print('error', f'*** trace has to be {genlib.get_trace_code_list_text()}.')
        OK = False
    if args.trace.upper() == 'Y':
        genlib.Message.set_trace_status(True)

    # if there are errors, exit with exception
    if not OK:
        raise genlib.ProgramException('', 'P001')

#-------------------------------------------------------------------------------

def get_cluster_homology_relationships(conn, blastp_alignment_file, homology_relationships_file):
    '''
    Get the homology relationships corresponding to set of sequence clusters yielded by the alignment process.
    WARNING: The alignment file should contain only one record per query sequence.
    '''

    # open the alignment file yielded by blastp
    if blastp_alignment_file.endswith('.gz'):
        try:
            blastp_alignment_file_id = gzip.open(blastp_alignment_file, mode='rt', encoding='iso-8859-1')
        except Exception as e:
            raise genlib.ProgramException(e, 'F002', blastp_alignment_file)
    else:
        try:
            blastp_alignment_file_id = open(blastp_alignment_file, mode='r', encoding='iso-8859-1')
        except Exception as e:
            raise genlib.ProgramException(e, 'F001', blastp_alignment_file)

    # initialize the counter of records corresponding to the alignment file yielded by blastp
    blastp_alignment_record_counter = 0

    # open the homology relationships file
    if homology_relationships_file.endswith('.gz'):
        try:
            homology_relationships_file_id = gzip.open(homology_relationships_file, mode='wt', encoding='iso-8859-1', newline='\n')
        except Exception as e:
            raise genlib.ProgramException(e, 'F004', homology_relationships_file)
    else:
        try:
            homology_relationships_file_id = open(homology_relationships_file, mode='w', encoding='iso-8859-1', newline='\n')
        except Exception as e:
            raise genlib.ProgramException(e, 'F003', homology_relationships_file)

    # write head in the homology relationships file
    homology_relationships_record = 'Sequence id;Species id;Homologous gene id;Homologous protein isoforms'
    homology_relationships_file_id.write(f'{homology_relationships_record}\n')

    # read the first record of alignment file yielded by blastp
    (blastp_alignment_record, _, blastp_alignment_data_dict) = genlib.read_alignment_outfmt6_record(blastp_alignment_file, blastp_alignment_file_id, blastp_alignment_record_counter)

    # while there are records in the alignment file yielded by blastp
    while blastp_alignment_record != '':

        # get alignment data
        qseqid = blastp_alignment_data_dict['qseqid']
        sseqid = blastp_alignment_data_dict['sseqid']

        # get the reference protein identification
        relationships_dict = sqllib.get_mmseqs2_protein_clusters_dict(conn, sseqid)
        for key, data in relationships_dict.items():
            reference_protein_id = data['seq_id']

        # get homology relationships dictionary
        (homology_relationships_dict) = get_homology_relationships(conn, reference_protein_id)

        # write the homology relationships in the homology relationships file
        homology_relationships_record = ''
        if homology_relationships_dict:
            for key in sorted(homology_relationships_dict):
                data = homology_relationships_dict[key]
                species_id = data['species_id']
                gene_id = data['gene_id']
                protein_isoform_ids = data['protein_isoform_ids']
                homology_relationships_record = f'{qseqid};{species_id};{gene_id};{protein_isoform_ids}'
                homology_relationships_file_id.write(f'{homology_relationships_record}\n')
        else:
            homology_relationships_record = f'{qseqid};-;-;-'
            homology_relationships_file_id.write(f'{homology_relationships_record}\n')

        # print counters
        genlib.Message.print('verbose', f'\rblastp clade alignment file: {blastp_alignment_record_counter} processed records')

        # read the next record of clade alignment file yielded by blastp
        (blastp_alignment_record, _, blastp_alignment_data_dict) = genlib.read_alignment_outfmt6_record(blastp_alignment_file, blastp_alignment_file_id, blastp_alignment_record_counter)

    genlib.Message.print('verbose', '\n')

    # close files
    blastp_alignment_file_id.close()
    homology_relationships_file_id.close()

    genlib.Message.print('verbose', '\n')
    genlib.Message.print('info', f'The file {homology_relationships_file} is created.')

#-------------------------------------------------------------------------------

def get_homology_relationships(conn, reference_protein_id):
    '''
    Get the homology relationships of a protein identification.
    '''

    # initialize the homology relationhips dictionary
    homology_relationships_dict = genlib.NestedDefaultDict()

    # get the list of protein isoforms corresponding to the reference protein identification
    mmseqs2_protein_isoforms_list = sqllib.get_mmseqs2_protein_isoforms_list(conn, '', [reference_protein_id])

    # build the list of protein isoforms identification corresponding to the reference protein identification
    reference_protein_isoform_ids_list = []
    for _, protein_isoform_data in enumerate(mmseqs2_protein_isoforms_list):
        reference_protein_isoform_ids_list.append(protein_isoform_data['protein_id'])

    # check if there are items in the list of protein isoforms corresponding to the reference protein identification
    if reference_protein_isoform_ids_list:

        # add data of reference protein to the homology relationships dictionary
        species_id = mmseqs2_protein_isoforms_list[0]['species_id']
        species_name = genlib.get_quercus_species_name(species_id)
        gene_id = mmseqs2_protein_isoforms_list[0]['gene_id']
        protein_isoform_ids = '|'.join(sorted(reference_protein_isoform_ids_list))
        homology_relationships_dict[f'{species_id}-{gene_id}'] = {'species_id': species_id, 'species_name': species_name, 'gene_id':gene_id, 'protein_isoform_ids': protein_isoform_ids}

        # for each identification in the list of protein isoforms corresponding to the reference protein identification
        for reference_protein_isoform_id in reference_protein_isoform_ids_list:

            # get data of the orthologous proteins
            orthologous_protein_data_list = sqllib.get_orthologous_protein_data_list(conn, reference_protein_isoform_id)

            # check if there are data of orthologous proteins:
            if orthologous_protein_data_list:

                # cluster the list of target orthologous proteins by species identification and gene identification
                clustered_homologous_proteins_dict = defaultdict(list)
                for item in orthologous_protein_data_list:
                    key = (item['target_species_id'], item['gene_id'])
                    clustered_homologous_proteins_dict[key].append(item['target_protein_id'])

                # build the orthologous data list
                orthologous_data_list = [
                    {
                        'species_id': species_id,
                        'gene_id': gene_id,
                        'protein_isoform_ids': '|'.join(sorted(protein_ids_list))
                    }
                    for (species_id, gene_id), protein_ids_list in clustered_homologous_proteins_dict.items()
                ]

                # add orthologous data to the homology relationships dictionary
                for _, data in enumerate(orthologous_data_list):
                    species_id = data['species_id']
                    species_name = genlib.get_quercus_species_name(species_id)
                    gene_id = data['gene_id']
                    protein_isoform_ids = data['protein_isoform_ids']
                    homology_relationships_dict[f'{species_id}-{gene_id}'] = {'species_id': species_id, 'species_name': species_name, 'gene_id':gene_id, 'protein_isoform_ids': protein_isoform_ids}

    # return the homology relationships dictionary
    return homology_relationships_dict

#-------------------------------------------------------------------------------

if __name__ == '__main__':

    main()
    sys.exit(0)

#-------------------------------------------------------------------------------
