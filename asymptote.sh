#!/usr/bin/env bash
# Executes asymptote.cpp over all_SRA_introns.tsv.gz for random seeds on [0, 50]; merges result
# $1: path to all_SRA_introns.tsv.gz
# $2: where to dump split files and results
set -e

ALLSRAINTRONS=$1
DUMPDIR=$2
DIR=$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )
cd $DIR
mkdir -p $DUMPDIR
g++ asymptote.cpp -o $DUMPDIR/asymptote
cd $DUMPDIR
gzip -cd $ALLSRAINTRONS | $DUMPDIR/asymptote 0 >0.res &
gzip -cd $ALLSRAINTRONS | $DUMPDIR/asymptote 1 >1.res &
gzip -cd $ALLSRAINTRONS | $DUMPDIR/asymptote 2 >2.res &
gzip -cd $ALLSRAINTRONS | $DUMPDIR/asymptote 3 >3.res &
wait
gzip -cd $ALLSRAINTRONS | $DUMPDIR/asymptote 4 >4.res &
gzip -cd $ALLSRAINTRONS | $DUMPDIR/asymptote 5 >5.res &
gzip -cd $ALLSRAINTRONS | $DUMPDIR/asymptote 6 >6.res &
gzip -cd $ALLSRAINTRONS | $DUMPDIR/asymptote 7 >7.res &
wait
gzip -cd $ALLSRAINTRONS | $DUMPDIR/asymptote 8 >8.res &
gzip -cd $ALLSRAINTRONS | $DUMPDIR/asymptote 9 >9.res &
gzip -cd $ALLSRAINTRONS | $DUMPDIR/asymptote 10 >10.res &
gzip -cd $ALLSRAINTRONS | $DUMPDIR/asymptote 11 >11.res &
wait
cat 0.res 1.res 2.res 3.res 4.res 5.res 6.res 7.res 8.res 9.res 10.res 11.res 12.res >all_asymptote_results.tsv
rm -rf *.res