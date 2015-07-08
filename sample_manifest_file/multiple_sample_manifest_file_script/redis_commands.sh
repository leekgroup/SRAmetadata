#!/usr/bin/env bash
# Commands for loading the redis database with introns
DIR=$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )
cd $DIR
REDIS=/scratch2/langmead-fs1/redis/redis-3.0.2/src
SRAINTRONS=/scratch0/langmead-fs1/sraintrons
MAWK=/scratch0/langmead-fs1/sraintrons/mawk-1.3.4-20150503/mawk
PYTHON=pypy

$MAWK -F ',' '{print "zadd numsamples " (NF+1)/2 " " NR}' $SRAINTRONS/all_SRA_introns.tsv | $PYTHON proto.py | $REDIS/redis-cli --pipe --pipe-timeout 0
cut -f8 $SRAINTRONS/all_SRA_introns.tsv | $MAWK -F ',' '{maxval=$1;for(i=2;i<=NF;i++){if($i>maxval){maxval=$i;}} print "zadd maxcoverage " maxval " " NR}' | $PYTHON proto.py | $REDIS/redis-cli --pipe --pipe-timeout 0
cut -f8 $SRAINTRONS/all_SRA_introns.tsv | $MAWK -F ',' '{minval=$1;for(i=2;i<=NF;i++){if($i<minval){minval=$i;}} print "zadd mincoverage " minval " " NR}' | $PYTHON proto.py | $REDIS/redis-cli --pipe --pipe-timeout 0
$MAWK '{print "hmset " NR " chrom " $1 " start " $2 " end " $3 " strand " $4 " motif " $5$6}' $SRAINTRONS/all_SRA_introns.tsv | $PYTHON proto.py | $REDIS/redis-cli --pipe --pipe-timeout 0