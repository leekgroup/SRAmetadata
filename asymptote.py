"""
asymptote.py
Part of SRA project

Consider a filter F that selects introns either a) occurring in some proportion
P of S RNA-seq samples or b) covered by at least R reads in any one sample.
This script operates on a set of introns initially detected across the S
samples, applying F to both S and random samples of S. For each random sample
s of S, the following tab-delimited quantities are output.
1) Size of s = |s|
2) Number of introns in s that passed filter
3) Number of introns N in common between s and S.
4) Precision = N / |s|
5) Recall = N / |S|

Input (tab-delimited fields read from stdin):
1) Strand (e.g., chr1+)
2) Intron start position
3) Intron end position
4) Comma-separated list of sample indexes
5) Comma-separated list of numbers of reads in which intron was initially
detected in samples from 4)

The input is in exactly the same format as Rail-RNA's "itn" deliverable.
"""
import sys
import random
import multiprocessing # faster faster
import time
import signal

def init_worker():
    """ Prevents KeyboardInterrupt from reaching a pool's workers.

        Exiting gracefully after KeyboardInterrupt or SystemExit is a
        challenge. The solution implemented here is by John Reese and is from
        http://noswap.com/blog/python-multiprocessing-keyboardinterrupt .

        No return value.
    """
    signal.signal(signal.SIGINT, signal.SIG_IGN)

def counts(populations, lines, min_reads=None, sample_fraction=None,
            indexes=range(3000)):
    """ Finds counts for computing intron precisions and recalls

        populations: list of sets of sample indexes randomly sampled from
            the set of possible sample indexes
        lines: iterable of some subset of lines input from stdin
        min_reads: minimum number of reads in which intron should be detected
            in any one sample to guarantee passing filter OR None if no
            such criterion should be used
        sample_fraction: minimum proportion of samples in which intron should
            be detected to guarantee passing filter OR None if no such
            criterion should be used
        indexes: list of all sample indexes

        Return value: tuple (total number of introns evaluated,
                number of introns passing filter for full range of sample
                indexes,
                list of numbers of introns passing filter for each
                corresponding sample index subset in the list populations,
                list of numbers of introns passing both filters (for
                    full range of sample indexes and sample index
                    subset corresponding to list item from populations)
            )
    """
    population_size = len(populations)
    passes, intersections = [0]*population_size, [0]*population_size
    unfiltered, sample_index_count = 0, len(indexes)
    if sample_fraction is not None:
        sample_fraction = float(sample_fraction)
        sample_threshold = round(sample_index_count * sample_fraction)
    else:
        sample_threshold = 0
    for line in lines:
        if not line: continue
        tokens = line.strip().split('\t')
        found_indexes = [int(index) for index in tokens[3].split(',')]
        found_index_set = set(found_indexes)
        coverages = [int(coverage) for coverage in tokens[4].split(',')]
        read_criterion = (min_reads is not None
                            and max(coverages) >= min_reads)
        sample_criterion = (sample_fraction is not None
                            and len(found_indexes) >= sample_threshold)
        no_filter = (min_reads is None and sample_fraction is None)
        intersect_it = False
        if no_filter or read_criterion or sample_criterion:
            unfiltered += 1
            intersect_it = True
        for i in xrange(population_size):
            if ((no_filter and populations[i].intersection(found_index_set))
                or (sample_fraction is not None
                        and len(populations[i].intersection(found_index_set))
                            >= round(len(populations[i]) * sample_fraction))
                or (min_reads is not None
                        and max([coverage for j, coverage
                                    in enumerate(coverages)
                                    if found_indexes[j] in populations[i]])
                            >= min_reads)):
                passes[i] += 1
                if intersect_it:
                    intersections[i] += 1
    return (len(lines), unfiltered, passes, intersections)

if __name__ == '__main__':
    start_time = time.time()
    import argparse
    # Print file's docstring if -h is invoked
    parser = argparse.ArgumentParser(description=__doc__, 
                formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument('--min-reads', type=int, required=False,
        default=None,
        help=('minimum number of reads in which a given intron should '
              'appear in any sample to pass filter; leave unspecified '
              'for no minimum'))
    parser.add_argument('--sample-fraction', type=float, required=False,
        default=None,
        help=('minimum proportion of samples in which a given intron '
              'should occur to pass filter; leave unspecified for '
              'no minimum; sample threshold is rounded to nearest integer'))
    parser.add_argument('--sample-index-count', type=int, required=False,
        default=3000,
        help=('number of samples spanned by sample indexes; this is 3000 for '
              'the SRA project'))
    parser.add_argument('--interval', type=int, required=False,
        default=50,
        help='interval between successive sample sizes')
    parser.add_argument('--sample-count', type=int, required=False,
        default=100,
        help=('number of samples of all sample indexes to take at each'
              'sample size'))
    parser.add_argument('--chunk-size', type=int, required=False,
        default=50000,
        help='number of input lines a thread should analyze at a time')
    parser.add_argument('--num-processes', '-p', type=int, required=False,
        default=None,
        help=('number of processes run simultaneously to analyze chunks; '
              'defaults to number of available processing cores'))
    args = parser.parse_args()
    sys.stderr.write('\x1b[KPreparing chunks for processing...\r')
    indexes = range(args.sample_index_count)
    random.seed(args.sample_fraction + args.min_reads)
    populations = [
                set(random.sample(indexes, sample_size))
                for sample_size in xrange(
                        args.interval, args.sample_index_count, args.interval
                    )
                for _ in xrange(args.sample_count)
            ]
    if args.num_processes is None:
        num_processes = multiprocessing.cpu_count()
    else:
        num_processes = args.num_processes
    pool = multiprocessing.Pool(num_processes, init_worker, maxtasksperchild=5)
    results = []
    to_dispatch = filter(lambda x: x is not None,
                         [next(sys.stdin, None)
                            for _ in xrange(args.chunk_size)])
    dispatch_count, active = 0, 0
    while to_dispatch:
        len_results_before = len(results)
        pool.apply_async(counts, [populations, to_dispatch, args.min_reads,
                                    args.sample_fraction, indexes],
                            callback=results.append)
        dispatch_count += 1
        active += 1
        while active >= num_processes:
            sys.stderr.write('\x1b[KChunks processed: %d\r' % len(results))
            active -= len(results) - len_results_before
            time.sleep(0.4)
        to_dispatch = filter(lambda x: x is not None,
                             [next(sys.stdin, None)
                                for _ in xrange(args.chunk_size)]
                        )
    while len(results) < dispatch_count:
        sys.stderr.write('\x1b[KChunks processed: %d\r' % len(results))
        time.sleep(0.4)
    process_time = time.time()
    print >>sys.stderr, '\x1b[KChunks processed in %02f s.' % (
                process_time - start_time
            )
    sys.stderr.write('\x1b[KComputing and writing output...\r')
    results = zip(*results)
    total_introns = sum(results[0])
    total_filtered_introns = sum(results[1])
    sample_filtered_introns = [sum(chunks) for chunks in zip(*results[2])]
    sample_intersects = [sum(chunks) for chunks in zip(*results[3])]
    line_format = '%d\t%d\t%d\t%08f\t%08f'
    print >>sys.stderr, ('\x1b[K%d intron(s) from %d samples were processed, '
                         'and %d intron(s) made it through the filter.') % (
                                total_introns,
                                args.sample_index_count,
                                total_filtered_introns
                            )
    for i in xrange(len(populations)):
        print line_format % (
                    len(populations[i]),
                    sample_filtered_introns[i],  
                    sample_intersects[i],
                    (float(sample_intersects[i]) / sample_filtered_introns[i])
                    if sample_filtered_introns[i] else 0.0,
                    (float(sample_intersects[i]) / total_filtered_introns)
                    if total_filtered_introns else 0.0
                )
    end_time = time.time()
    print >>sys.stderr, ('\x1b[KOutput computed and written in %02f s. Entire '
                         'job finished in %02f s.') % (end_time - process_time,
                                                       end_time - start_time)