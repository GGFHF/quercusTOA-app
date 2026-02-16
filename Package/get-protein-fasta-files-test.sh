#!/bin/bash

#-------------------------------------------------------------------------------

# This script executes a test of the program get-cluster-homology-relationships.py
# in a Linux environment.
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

# Control parameters

if [ -n "$*" ]; then echo 'This script does not have parameters'; exit 1; fi

#-------------------------------------------------------------------------------

# Set run environment

APP_DIR=$QUERCUSTOA
DATA_DIR=$APP_DIR/data
OUTPUT_DIR=$APP_DIR/output

if [ ! -d "$OUTPUT_DIR" ]; then mkdir --parents $OUTPUT_DIR; fi

INITIAL_DIR=$(pwd)
cd $APP_DIR

#-------------------------------------------------------------------------------

# Run the program get-cluster-homology-relationships.py

/usr/bin/time \
    ./get-cluster-homology-relationships.py \
        --sequences-db=$DATA_DIR/sequences.db \
        --homology=$DATA_DIR/cluster-homology-relationships.csv \
        --blastp-alignments=$DATA_DIR/blastp-alignments.csv \
        --analysis=$DATA_DIR/analysis-sequences.fasta \
        --consensus=$DATA_DIR/Quercus-consensus-seqs.fasta \
        --outdir=$OUTPUT_DIR \
        --verbose=Y  \
        --trace=N
if [ $? -ne 0 ]; then echo 'Script ended with errors.'; exit 1; fi

#-------------------------------------------------------------------------------

# End

cd $INITIAL_DIR

exit 0

#-------------------------------------------------------------------------------
