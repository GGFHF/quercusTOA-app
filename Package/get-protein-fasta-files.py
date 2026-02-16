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
This program gets the protein FASTA files corresponding to homology relationships
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

    # connect to the quercusTOA sequences database
    conn = sqllib.connect_database(args.sequences_database)

    # get the protein FASTA files corresponding to homology relationships yielded by process of searching for sequence homology
    get_protein_fasta_files(conn, args.homology_relationships_file, args.blastp_alignment_file, args.analysis_fasta_file, args.consensus_seqs_file, args.output_dir)

    # close connection to quercusTOA comparative genomics database
    conn.close()

#-------------------------------------------------------------------------------

def build_parser():
    '''
    Build the parser with the available arguments.
    '''

    # create the parser and add arguments
    description = 'Description: This program gets the protein FASTA files corresponding to homology relationships\n' \
       'yielded by process of searching for sequence homology.'
    text = f'{genlib.get_app_short_name()} v{genlib.get_app_version()} - {os.path.basename(__file__)}\n\n{description}\n'
    usage = f'\r{text.ljust(len("usage:"))}\nUsage: {os.path.basename(__file__)} arguments'
    parser = argparse.ArgumentParser(usage=usage)
    parser._optionals.title = 'Arguments'    #pylint: disable=protected-access
    parser.add_argument('--sequences-db', dest='sequences_database', help=f'Path of the {genlib.get_app_short_name} sequences database (mandatory).')
    parser.add_argument('--homology', dest='homology_relationships_file', help='Path of the homology relationships file (mandatory).')
    parser.add_argument('--blastp-alignments', dest='blastp_alignment_file', help='Path of the alignment file yielded by blastp (mandatory).')
    parser.add_argument('--analysis', dest='analysis_fasta_file', help='Path of the FASTA file with analysis protein sequences or NONE; default: NONE.')
    parser.add_argument('--consensus', dest='consensus_seqs_file', help='Path of the FASTA file with consensus sequences (mandatory).')
    parser.add_argument('--outdir', dest='output_dir', help='Path of output directoty where files with selected variant data are saved (mandatory).')
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

    # check "sequences_database"
    if args.sequences_database is None:
        genlib.Message.print('error', f'*** The {genlib.get_app_short_name} sequence database is not indicated in the input arguments.')
        OK = False
    elif not os.path.isfile(args.sequences_database):
        genlib.Message.print('error', f'*** The file {args.sequences_database} does not exist.')
        OK = False

    # check "homology_relationships_file"
    if args.homology_relationships_file is None:
        genlib.Message.print('error', '*** The homology relationships file is not indicated in the input arguments.')
        OK = False
    elif not os.path.isfile(args.homology_relationships_file):
        genlib.Message.print('error', f'*** The file {args.homology_relationships_file} does not exist.')
        OK = False

    # check "blastp_alignment_file"
    if args.blastp_alignment_file is None:
        genlib.Message.print('error', '*** The blastp alignment file is not indicated in the input arguments.')
        OK = False
    elif not os.path.isfile(args.blastp_alignment_file):
        genlib.Message.print('error', f'*** The file {args.blastp_alignment_file} does not exist.')
        OK = False

    # check "analysis_fasta_file"
    if args.analysis_fasta_file is None or args.analysis_fasta_file.upper() == 'NONE':
        args.analysis_fasta_file = 'NONE'
    elif not os.path.isfile(args.analysis_fasta_file):
        genlib.Message.print('error', f'*** The file {args.analysis_fasta_file} does not exist.')
        OK = False

    # check "consensus_seqs_file"
    if args.consensus_seqs_file is None:
        genlib.Message.print('error', '*** The consensus sequence file is not indicated in the input arguments.')
        OK = False
    elif not os.path.isfile(args.consensus_seqs_file):
        genlib.Message.print('error', f'*** The file {args.consensus_seqs_file} does not exist.')
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

def get_protein_fasta_files(conn, homology_relationships_file, blastp_alignments_file, analysis_fasta_file, consensus_seqs_file, output_dir):
    '''
    Get the protein FASTA files corresponding to homology relationships yielded by process
    of searching for sequence homology.
    '''

    # get the blastp alingments dictionary
    blastp_alignments_dict = get_blaspt_alignment_dict(blastp_alignments_file)

    # get the analysis FASTA dictionary
    analysis_fasta_dict = {}
    if analysis_fasta_file != 'NONE':
        analysis_fasta_dict = genlib.get_fasta_seq_dict(analysis_fasta_file, cutting_char=' ')

    # get the consensus sequences dictionary
    consensus_seqs_dict = genlib.get_fasta_seq_dict(consensus_seqs_file, cutting_char=' ')

    # open the homology relationships file
    if homology_relationships_file.endswith('.gz'):
        try:
            homology_relationships_file_id = gzip.open(homology_relationships_file, mode='rt', encoding='iso-8859-1')
        except Exception as e:
            raise genlib.ProgramException(e, 'F002', homology_relationships_file)
    else:
        try:
            homology_relationships_file_id = open(homology_relationships_file, mode='r', encoding='iso-8859-1')
        except Exception as e:
            raise genlib.ProgramException(e, 'F001', homology_relationships_file)

    # set the header indicator
    is_header = True

    # initialize the counter of records corresponding to the homology relationships file
    homology_relationships_record_counter = 0

    # read the first record of the homology relationships file
    (homology_relationships_record, _, homology_relationships_data_dict) = genlib.read_homology_relationships_record(homology_relationships_file, homology_relationships_file_id, homology_relationships_record_counter)

    # while there are records in the homology relationships file
    while homology_relationships_record != '':

        # when is the head record
        if is_header:

            # add 1 to the counter of records corresponding to the homology relationships file
            homology_relationships_record_counter += 1

            # set the header indicator
            is_header = False

            # read the next record of the homology relationships file
            (homology_relationships_record, _, homology_relationships_data_dict) = genlib.read_homology_relationships_record(homology_relationships_file, homology_relationships_file_id, homology_relationships_record_counter)

        # when is a data record
        else:

            # initialize the old sequence identification
            old_seq_id = homology_relationships_data_dict['seq_id']

            # set the file path of homologous protein FASTA sequences corresponding to the current sequence identification
            protein_sequence_fasta_file = f'{output_dir}{os.sep}{homology_relationships_data_dict['seq_id']}-homologous-proteins.fasta'

            # open the protein FASTA sequence file
            if protein_sequence_fasta_file.endswith('.gz'):
                try:
                    protein_sequence_fasta_file_id = gzip.open(protein_sequence_fasta_file, mode='wt', encoding='iso-8859-1', newline='\n')
                except Exception as e:
                    raise genlib.ProgramException(e, 'F004', protein_sequence_fasta_file)
            else:
                try:
                    protein_sequence_fasta_file_id = open(protein_sequence_fasta_file, mode='w', encoding='iso-8859-1', newline='\n')
                except Exception as e:
                    raise genlib.ProgramException(e, 'F003', protein_sequence_fasta_file)

            # set the file path of homologous gene FASTA sequences corresponding to the current sequence identification
            gene_sequence_fasta_file = f'{output_dir}{os.sep}{homology_relationships_data_dict['seq_id']}-homologous-genes.fasta'

            # open the gene sequence FASTA file
            if gene_sequence_fasta_file.endswith('.gz'):
                try:
                    gene_sequence_fasta_file_id = gzip.open(gene_sequence_fasta_file, mode='wt', encoding='iso-8859-1', newline='\n')
                except Exception as e:
                    raise genlib.ProgramException(e, 'F004', gene_sequence_fasta_file)
            else:
                try:
                    gene_sequence_fasta_file_id = open(gene_sequence_fasta_file, mode='w', encoding='iso-8859-1', newline='\n')
                except Exception as e:
                    raise genlib.ProgramException(e, 'F003', gene_sequence_fasta_file)

            # get the analysis sequence
            analysis_seq = analysis_fasta_dict.get(homology_relationships_data_dict['seq_id'], '')

            # write the analysis sequence in the  the protein FASTA sequence file
            if analysis_seq != '':
                protein_sequence_fasta_file_id.write(f'>{homology_relationships_data_dict['seq_id']}\n')
                protein_sequence_fasta_file_id.write(f'{analysis_seq}\n')

            # the the cluster identification
            cluster_id = blastp_alignments_dict[homology_relationships_data_dict['seq_id']]

            # get the consensus sequence
            consensus_seq = consensus_seqs_dict[cluster_id]

            # write the consensus sequence in the  the protein FASTA sequence file
            protein_sequence_fasta_file_id.write(f'>quercusTOA-{cluster_id}\n')
            protein_sequence_fasta_file_id.write(f'{consensus_seq}\n')

            # while there are records in the homology relationships file and same sequence identification
            while homology_relationships_record != '' and homology_relationships_data_dict['seq_id'] ==  old_seq_id:

                # add 1 to the counter of records corresponding to the homology relationships file
                homology_relationships_record_counter += 1

                # set the protein isoform identification list
                protein_isoform_id_list = homology_relationships_data_dict['protein_isoform_ids'].split('|')

                # for each protein isoform identification
                for protein_isoform_id in protein_isoform_id_list:

                    # get the protein isoform sequence data
                    protein_isoform_seq_dict = sqllib.get_protein_seq_dict(conn, protein_isoform_id)
                    protein_isoform_seq = protein_isoform_seq_dict['seq']

                    # write the protein isoform sequence in the  the protein FASTA sequence file
                    protein_sequence_fasta_file_id.write(f'>{protein_isoform_id}[{homology_relationships_data_dict['species_id']}]\n')
                    protein_sequence_fasta_file_id.write(f'{protein_isoform_seq}\n')

                    # get the gene sequence data of protein isoform
                    gene_seq_dict = sqllib.get_gene_seq_dict(conn, homology_relationships_data_dict['gene_id'])
                    gene_isoform_seq = gene_seq_dict['seq']

                    # write the gene sequence of the protein isoform in the the gene FASTA sequence file
                    gene_sequence_fasta_file_id.write(f'>{homology_relationships_data_dict['gene_id']}[{homology_relationships_data_dict['species_id']}]\n')
                    gene_sequence_fasta_file_id.write(f'{gene_isoform_seq}\n')

                # print the record counters
                genlib.Message.print('verbose', f'\rHolology relationships records: {homology_relationships_record_counter:5d}')

                # read the next record of the homology relationships file
                (homology_relationships_record, _, homology_relationships_data_dict) = genlib.read_homology_relationships_record(homology_relationships_file, homology_relationships_file_id, homology_relationships_record_counter)

            # close protein FASTA sequence file of the current sequence identification
            protein_sequence_fasta_file_id.close()

            # close gene FASTA sequence file of the current sequence identification
            gene_sequence_fasta_file_id.close()

    genlib.Message.print('verbose', '\n')

    # close the homology relationships file
    homology_relationships_file_id.close()

    genlib.Message.print('verbose', '\n')
    genlib.Message.print('info', 'The protein FASTA files corresponding to homology relationships are created.')

#-------------------------------------------------------------------------------


def get_blaspt_alignment_dict(file_path):
    '''
    Get the alignment dictionary from the corresponding blastp alignment file yielded by
    the search homology of sequences.
    '''

    # initialize the alignment dictionary
    alignments_dict = {}

    # open the alignment file
    if file_path.endswith('.gz'):
        try:
            file_id = gzip.open(file_path, mode='rt', encoding='iso-8859-1')
        except Exception as e:
            raise genlib.ProgramException(e, 'F002', file_path) from  e
    else:
        try:
            file_id = open(file_path, mode='r', encoding='iso-8859-1')
        except Exception as e:
            raise genlib.ProgramException(e, 'F001', file_path) from e

    # initialize the counter of records corresponding to the alignment file
    alignment_record_counter = 0

    # read the first record of alignment file
    (blastp_alignment_record, _, blastp_alignment_data_dict) = genlib.read_alignment_outfmt6_record(file_path, file_id, alignment_record_counter)

    # while there are records in the alignment file
    while blastp_alignment_record != '':

        # get alignment data
        qseqid = blastp_alignment_data_dict['qseqid']
        sseqid = blastp_alignment_data_dict['sseqid']

        # insert data in the alignment dictionary
        alignments_dict[qseqid] = sseqid

        # print counters
        genlib.Message.print('verbose', f'\rblastp clade alignment file: {alignment_record_counter} processed records')

        # read the next record of alignment file
        (blastp_alignment_record, _, blastp_alignment_data_dict) = genlib.read_alignment_outfmt6_record(file_path, file_id, alignment_record_counter)

    genlib.Message.print('verbose', '\n')

    # close files
    file_id.close()

    # return the BLASP alignment dictionary
    return alignments_dict

#-------------------------------------------------------------------------------

if __name__ == '__main__':

    main()
    sys.exit(0)

#-------------------------------------------------------------------------------
