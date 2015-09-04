#!/usr/bin/env bash
# Splits all_SRA_introns.tsv.gz into 100 partitions, compiles asymptote.cpp
# $1: path to all_SRA_introns.tsv.gz
# $2: where to dump split files and results
set -e

ALLSRAINTRONS=$1
DUMPDIR=$2
DIR=$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )
cd $DIR
g++ asymptote.cpp -o $DUMPDIR/asymptote
cd $DUMPDIR
gzip -cd $ALLSRAINTRONS | split -b 428821 - srasplit
# Now execute $DUMPDIR/asymptote on each srasplit file