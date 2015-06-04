#!/usr/bin/env bash
# Used PyPy 2.4.0 with GCC 4.8.2
gzip -cd collected_introns.tsv.gz | pypy asymptote.py --sample-fraction 0.05 -p 31 >intron_accuracy_sample_fraction_0.05.tsv
gzip -cd collected_introns.tsv.gz | pypy asymptote.py --sample-fraction 0.03 -p 31 >intron_accuracy_sample_fraction_0.03.tsv
gzip -cd collected_introns.tsv.gz | pypy asymptote.py --sample-fraction 0.10 -p 31 >intron_accuracy_sample_fraction_0.10.tsv
gzip -cd collected_introns.tsv.gz | pypy asymptote.py --sample-fraction 0.05 --read-min 5 -p 31 >intron_accuracy_sample_fraction_0.05_read_min_5.tsv
gzip -cd collected_introns.tsv.gz | pypy asymptote.py --sample-fraction 0.05 --read-min 10 -p 31 >intron_accuracy_sample_fraction_0.05_read_min_10.tsv
gzip -cd collected_introns.tsv.gz | pypy asymptote.py --sample-fraction 0.05 --read-min 15 -p 31 >intron_accuracy_sample_fraction_0.05_read_min_15.tsv