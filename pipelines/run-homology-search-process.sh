#!/bin/bash

#-------------------------------------------------------------------------------

# This script processes a sequence homology search using quercusTOA-app.
#
# This software has been developed by:
#
#    GI en Desarrollo de Especies y Comunidades Leñosas (WooSp)
#    Dpto. Sistemas y Recursos Naturales
#    ETSI Montes, Forestal y del Medio Natural
#    Universidad Politecnica de Madrid
#    https://github.com/ggfhf/
#
# Licence: GNU General Public Licence Version 3.

#-------------------------------------------------------------------------------

# Infrastructure
# ==============
#
# Before executing this script it is necessary that the software used is installed
# (Miniforge3, CodAn, BLAST+ and DIAMOND) and the quercusTOA-db is downloaded. Check
# how to perform these actions in the capter "Functional annotation and enrichment
# analysis on Linux servers" of quercusTOA-app manual.

#-------------------------------------------------------------------------------

# Control parameters

PARAM_NUM=7
if [ "$#" -ne "$PARAM_NUM" ]; then
    echo "*** ERROR: The following $PARAM_NUM parameters are required:"
    echo
    echo '    quercustoa_app_dir <- path of the quercusTOA-app directory.'
    echo '    quercustoa_db_dir <- path of the quercusTOA-db directory.'
    echo '    analysis_type <- FASTA file type: "TRANSCRIPTS" or "PROTEINS".'
    echo '    analysis_file <- FASTA file with analysis protein sequences.'
    echo '    model <- CodAn model: "PLANTS_full" or "PLANTS_partial".'
    echo '    threads <- threads number.'
    echo '    homology_dir <- path of the homology search output directory.'
    echo
    echo "Use: ${0##*/} quercustoa_app_dir quercustoa_db_dir analysis_type analysis_file model threads homology_dir"
    exit 1
fi

QUERCUSTOA_APP_DIR=${1}
QUERCUSTOA_DB_DIR=${2}
ANALYSIS_TYPE=${3}
ANALYSIS_FILE=${4}
MODEL=${5}
THREADS=${6}
HOMOLOGY_DIR=${7}

#-------------------------------------------------------------------------------

# set other variables

TEMP=$HOMOLOGY_DIR/temp
SEQS_ALIGNMENTS=$HOMOLOGY_DIR/seqs-alignments
SEP="#########################################"

#-------------------------------------------------------------------------------

# create directories

if [ -d "$HOMOLOGY_DIR" ]; then rm -rf $HOMOLOGY_DIR; fi; mkdir --parents $HOMOLOGY_DIR 
if [ -d "$TEMP" ]; then rm -rf $TEMP; fi; mkdir --parents $TEMP 
if [ -d "$SEQS_ALIGNMENTS" ]; then rm -rf $SEQS_ALIGNMENTS; fi; mkdir --parents $SEQS_ALIGNMENTS 

#-------------------------------------------------------------------------------

function init
{
    INIT_DATETIME=`date +%s`
    FORMATTED_INIT_DATETIME=`date "+%Y-%m-%d %H:%M:%S"`
    echo "$SEP"
    echo "Enrichment analysis process started at $FORMATTED_INIT_DATETIME."
}

#-------------------------------------------------------------------------------

function predict_orfs
{
    echo "$SEP"
    echo "Predicting ORFs and getting peptide sequences ..."
    if [[ "$ANALYSIS_TYPE" == "TRANSCRIPTS" ]]; then
        cd $HOMOLOGY_DIR
        source activate quercustoa-codan
        CODAN_DIR=`echo $CONDA_PREFIX`
        /usr/bin/time \
            codan.py \
                --cpu=$THREADS \
                --model=$CODAN_DIR/models/$MODEL \
                --transcripts=$ANALYSIS_FILE \
                --output=$HOMOLOGY_DIR/codan_output
        RC=$?
        if [ $RC -ne 0 ]; then manage_error codan.py $RC; fi
        conda deactivate
        echo "ORFs are predicted and peptide sequences are gotten."
    elif [[ "$FASTA_TYPE" == "PROTEINS" ]]; then
        echo "This step is not run with a proteins file."
    fi
}

#-------------------------------------------------------------------------------

function align_peptides_2_alignment_tool_quercus_db
{
    echo "$SEP"
    echo "Aligning peptides to the aligner Quercus database ..."
    cd $HOMOLOGY_DIR
    if [[ "$ANALYSIS_TYPE" == "TRANSCRIPTS" ]]; then
        PEPTIDE_FILE=$HOMOLOGY_DIR/codan_output/PEP_sequences.fa
    elif [[ "$ANALYSIS_TYPE" == "PROTEINS" ]]; then
        PEPTIDE_FILE=$ANALYSIS_FILE
    fi
    source activate quercustoa-diamond
    /usr/bin/time \
        diamond blastp \
            --threads $THREADS \
            --db $QUERCUSTOA_DB_DIR/Quercus-consensus-diamond-db/Quercus-consensus-diamond-db \
            --query $PEPTIDE_FILE \
            --evalue 1e-06 \
            --max-target-seqs 1 \
            --max-hsps 1 \
            --outfmt 6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore \
            --out $TEMP/blastp-clade-alignments.csv
    RC=$?
    if [ $RC -ne 0 ]; then manage_error diamond-blastp $RC; fi
    conda deactivate
    echo "Alignment is done."
}
#-------------------------------------------------------------------------------

function get_cluster_homology_relationships
{
    echo "$SEP"
    echo "Getting cluster homology relationships ..."
    cd $HOMOLOGY_DIR
    source activate quercustoa
    /usr/bin/time \
        $QUERCUSTOA_APP_DIR/get-cluster-homology-relationships.py \
            --comparative-db=$QUERCUSTOA_DB_DIR/comparative-genomics.db \
            --annotations-db=$QUERCUSTOA_DB_DIR/functional-annotations.db \
            --blastp-alignments=$TEMP/blastp-clade-alignments.csv \
            --homology=$HOMOLOGY_DIR/homology-relationships-file.csv \
            --verbose=N \
            --trace=N
    RC=$?
    if [ $RC -ne 0 ]; then manage_error get-cluster-homology-relationships.py $RC; fi
    conda deactivate
    echo "Cluster homology relationships are gotten."
}

#-------------------------------------------------------------------------------

function get_protein_fasta_files
{
    echo "$SEP"
    echo "Getting protein FASTA files ..."
    cd $HOMOLOGY_DIR
    source activate quercustoa
    /usr/bin/time \
        $QUERCUSTOA_APP_DIR/get-protein-fasta-files.py \
            --sequences-db=$QUERCUSTOA_DB_DIR/sequences.db \
            --homology=$HOMOLOGY_DIR/homology-relationships-file.csv \
            --blastp-alignments=$TEMP/blastp-clade-alignments.csv \
            --analysis=$PEPTIDE_FILE \
            --consensus=$QUERCUSTOA_DB_DIR/Quercus-consensus-seqs.fasta \
            --outdir=$SEQS_ALIGNMENTS \
            --verbose=N \
            --trace=N
    RC=$?
    if [ $RC -ne 0 ]; then manage_error get-cluster-homology-relationships.py $RC; fi
    conda deactivate
    echo "Protein FASTA files are gotten."
}

#-------------------------------------------------------------------------------

function get_alignment_files
{
    echo "$SEP"
    echo "Getting alignment files and their plots ..."
    cd $HOMOLOGY_DIR
    source activate quercustoa-mafft
    FASTA_FILE_LIST=$SEQS_ALIGNMENTS/fasta_file_list.txt
    ls $SEQS_ALIGNMENTS/*.fasta > $FASTA_FILE_LIST
    while read FASTA_FILE; do
        ALIGNMENT_FILE=`echo $FASTA_FILE | sed "s/.fasta/.aln/g"`
        if [[ "$FASTA_FILE" == *homologous-proteins.fasta ]]; then
            TREE=Y
        else
            TREE=N
        fi
        PLOT_FILE=`echo $FASTA_FILE | sed "s/.fasta/.png/g"`
        /usr/bin/time \
            $QUERCUSTOA_APP_DIR/align-fasta-seqs.py \
                --seqs=$FASTA_FILE \
                --tree=$TREE \
                --verbose=N \
                --trace=N
        RC=$?
        if [ $RC -ne 0 ]; then manage_error align-fasta-seqs.py $RC; fi
    done < $FASTA_FILE_LIST
    conda deactivate
    echo "FASTA files are aligned."
}

#-------------------------------------------------------------------------------

function end
{
    END_DATETIME=`date +%s`
    FORMATTED_END_DATETIME=`date "+%Y-%m-%d %H:%M:%S"`
    calculate_duration
    echo "$SEP"
    echo "Homology search process ended OK at $FORMATTED_END_DATETIME with a run duration of $DURATION s ($FORMATTED_DURATION)."
    echo "$SEP"
    exit 0
}

#-------------------------------------------------------------------------------

function manage_error
{
    END_DATETIME=`date +%s`
    FORMATTED_END_DATETIME=`date "+%Y-%m-%d %H:%M:%S"`
    calculate_duration
    echo "$SEP"
    echo "ERROR: $1 returned error $2"
    echo "Homology search process ended WRONG at $FORMATTED_END_DATETIME with a run duration of $DURATION s ($FORMATTED_DURATION)."
    echo "$SEP"
    exit 3
}

#-------------------------------------------------------------------------------

function calculate_duration
{
    DURATION=`expr $END_DATETIME - $INIT_DATETIME`
    HH=`expr $DURATION / 3600`
    MM=`expr $DURATION % 3600 / 60`
    SS=`expr $DURATION % 60`
    FORMATTED_DURATION=`printf "%03d:%02d:%02d\n" $HH $MM $SS`
}

#-------------------------------------------------------------------------------

init
predict_orfs
align_peptides_2_alignment_tool_quercus_db
get_cluster_homology_relationships
get_protein_fasta_files
get_alignment_files
end
