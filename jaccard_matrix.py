"""
jaccard_matrix.py
Part of SRA project

Constructs Jaccard similarity matrix for intron sets across sample indexes.
Uses twice as much memory as it needs to constructing a symmetric matrix, but
who cares?

Input (tab-delimited fields read from stdin):
1) Strand (e.g., chr1+)
2) Intron start position
3) Intron end position
4) Comma-separated list of sample indexes
5) Comma-separated list of numbers of reads in which intron was initially
detected in samples from 4)

Tab-delimited output:
1) First sample index
2) Second sample index
3) Jaccard index

The input is in exactly the same format as Rail-RNA's "itn" deliverable.
"""
import sys
import multiprocessing # faster faster
import time
import signal
from collections import defaultdict
from itertools import combinations_with_replacement

def init_worker():
    """ Prevents KeyboardInterrupt from reaching a pool's workers.

        Exiting gracefully after KeyboardInterrupt or SystemExit is a
        challenge. The solution implemented here is by John Reese and is from
        http://noswap.com/blog/python-multiprocessing-keyboardinterrupt .

        No return value.
    """
    signal.signal(signal.SIGINT, signal.SIG_IGN)

def intersections_and_unions(lines, sample_count=3000):
    """ Finds intersections and unions across samples for set of intron lines.

        lines: iterable of some subset of lines input from stdin
        sample_count: number of samples

        Return value: 2-tuple whose elements are lists of lists i and u.
        i[k][l] is the number of introns in common between samples
        k and l, and u[k][l] is the number of introns found in either sample k
        or sample l.
    """
    intersections = [[0 for _ in xrange(sample_count)]
                        for __ in xrange(sample_count)]
    unions = [[0 for _ in xrange(sample_count)]
                 for __ in xrange(sample_count)]
    for line in lines:
        if not line: continue
        tokens = line.strip().split('\t')
        found_indexes = [int(index) for index in tokens[3].split(',')]
        for i, j in combinations_with_replacement(found_indexes, 2):
            intersections[i][j] += 1
            if i != j:
                intersections[j][i] += 1
        for i in xrange(sample_count):
            for j in found_indexes:
                unions[i][j] += 1
                if i != j:
                    unions[j][i] += 1
    return (intersections, unions)

if __name__ == '__main__':
    start_time = time.time()
    import argparse
    # Print file's docstring if -h is invoked
    parser = argparse.ArgumentParser(description=__doc__, 
                formatter_class=argparse.RawDescriptionHelpFormatter)
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
        pool.apply_async(intersections_and_unions,
                            [to_dispatch, args.sample_count],
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
    intersections = [[0 for _ in xrange(args.sample_count)]
                        for __ in xrange(args.sample_count)]
    unions = [[0 for _ in xrange(args.sample_count)]
                 for __ in xrange(args.sample_count)]
    result_count = len(results)
    for i in xrange(args.sample_count):
        for j in xrange(args.sample_count):
            intersections[i][j] = sum([results[k][0][i][j]
                                        for k in xrange(result_count)])
            unions[i][j] = sum([results[k][1][i][j]
                                    for k in xrange(result_count)])
    for i in xrange(args.sample_count):
        for j in xrange(i, args.sample_count):
            assert intersections[i][j] == intersections[j][i]
            assert unions[i][j] == unions[j][i]
            assert i != j or intersections[i][j] == unions[i][j]
            try:
                print '%d\t%d\t%.15f' % (
                        i, j, float(intersections[i][j]) / unions[i][j]
                    )
            except ZeroDivisionError:
                print '%d\t%d\tNA' % (i, j)
    end_time = time.time()
    print >>sys.stderr, ('\x1b[KOutput computed and written in %02f s. Entire '
                         'job finished in %02f s.') % (end_time - process_time,
                                                       end_time - start_time)