"""
jaccard_matrix.py
Part of SRA project

Constructs Jaccard similarity matrix for intron sets across sample indexes.
Uses twice as much memory as it needs to constructing a symmetric matrix, but
who cares?

Input (tab-delimited fields; read from file):
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

def counts(lines, intron_counts, mutex):
    """ Finds intron counts across samples

        mutex: for acquiring lock
        intron_counts: defaultdict that maps sample indexes to number of 
            introns
        lines: iterable of some subset of lines input from stdin

        Return value: 0 if successful
    """
    sample_counts = defaultdict(int)
    for line in lines:
        if not line: continue
        tokens = line.strip().split('\t')
        for index in tokens[3].split(','):
            mutex.acquire()
            sample_counts[int(index)] += 1
            mutex.release()
    return sample_counts

def intersections_and_unions(lines, intersections, unions, mutex,
                                forbidden_samples=set(),
                                sample_count=3000):
    """ Updates intersections/unions across samples for set of intron lines.

        lines: iterable of some subset of lines input from stdin
        intersections: matrix i for which i[k][l] is the number of introns in
            common between samples k and l
        unions: matrix u for which u[k][l] is the number of introns found in
            either sample k or sample l
        mutex: for acquiring lock
        forbidden_samples: set of samples to exclude from consideration
        sample_count: number of samples
        
        Return value: 2-tuple whose elements are lists of lists i and u.
        i[k][l] is the number of introns in common between samples
        k and l, and u[k][l] is the number of introns found in either sample k
        or sample l.
    """
    for line in lines:
        if not line: continue
        tokens = line.strip().split('\t')
        found_indexes = [int(index) for index in tokens[3].split(',')
                            if index not in forbidden_samples]
        for i, j in combinations_with_replacement(found_indexes, 2):
            intersections[i][j] += 1
            if i != j:
                mutex.acquire()
                intersections[j][i] += 1
                mutex.release()
        for i in xrange(sample_count):
            if i in forbidden_samples: continue
            for j in found_indexes:
                mutex.acquire()
                unions[i][j] += 1
                if i != j:
                    unions[j][i] += 1
                mutex.release()
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
    parser.add_argument('--input', '-i' type=str, required=True,
        help='path to input file')
    parser.add_argument('--filter', type=int, required=False,
        default=10000,
        help='filter out samples with fewer than this many detected introns')
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
    mutex = multiprocessing.Lock()
    intron_counts = defaultdict(int)
    results = []
    print >>sys.stderr, '\x1b[Counting introns...'
    with open(args.input) as input_stream:
        to_dispatch = filter(lambda x: x is not None,
                             [next(input_stream, None)
                                for _ in xrange(args.chunk_size)])
        dispatch_count, active = 0, 0
        while to_dispatch:
            len_results_before = len(results)
            pool.apply_async(counts, [to_dispatch, intersections, unions,
                                      mutex, forbidden_samples, ],
                                callback=results.append)
            dispatch_count += 1
            active += 1
            while active >= num_processes:
                sys.stderr.write('\x1b[KChunks processed: %d\r' % len(results))
                active -= len(results) - len_results_before
                time.sleep(0.4)
            to_dispatch = filter(lambda x: x is not None,
                                 [next(input_stream, None)
                                    for _ in xrange(args.chunk_size)]
                            )
    while len(results) < dispatch_count:
        sys.stderr.write('\x1b[KChunks processed: %d\r' % len(results))
        time.sleep(0.4)
    process_time = time.time()
    print >>sys.stderr, '\x1b[KChunks processed in %02f s.' % (
                process_time - start_time
            )
    forbidden_samples = set(
                [int(sample_index) for sample_index in intron_counts
                    if intron_counts[sample_index] < args.filter]
            )
    print >>sys.stderr, '\x1b[%d samples filtered out.' % len(
                                                            forbidden_samples
                                                        )
    print >>sys.stderr, '\x1b[Computing Jaccard matrix...'
    print ';'.join(map(str, sorted(forbidden_samples)))
    intersections = [[0 for _ in xrange(args.sample_count)]
                        for __ in xrange(args.sample_count)]
    unions = [[0 for _ in xrange(args.sample_count)]
                 for __ in xrange(args.sample_count)]
    results = []
    with open(args.input) as input_stream:
        to_dispatch = filter(lambda x: x is not None,
                             [next(input_stream, None)
                                for _ in xrange(args.chunk_size)])
        dispatch_count, active = 0, 0
        while to_dispatch:
            len_results_before = len(results)
            pool.apply_async(intersections_and_unions,
                                [to_dispatch, intron_counts, mutex],
                                callback=results.append)
            dispatch_count += 1
            active += 1
            while active >= num_processes:
                sys.stderr.write('\x1b[KChunks processed: %d\r' % len(results))
                active -= len(results) - len_results_before
                time.sleep(0.4)
            to_dispatch = filter(lambda x: x is not None,
                                 [next(input_stream, None)
                                    for _ in xrange(args.chunk_size)]
                            )
    while len(results) < dispatch_count:
        sys.stderr.write('\x1b[KChunks processed: %d\r' % len(results))
        time.sleep(0.4)
    second_process_time = time.time()
    print >>sys.stderr, '\x1b[KChunks processed in %02f s.' % (
                second_process_time - process_time
            )
    sys.stderr.write('\x1b[KComputing and writing output...\r')
    for i in xrange(args.sample_count):
        if i in forbidden_samples: continue
        for j in xrange(i, args.sample_count):
            if j in forbidden_samples: continue
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