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
from threading import Thread
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

def count_listener(q, intron_counts):
    """ Accumulates intron counts.

        q: queue with sample indexes
        intron_counts: defaultdict mapping sample indexes to intron counts

        No return value.
    """
    last_index = q.get()
    while last_index != -1:
        intron_counts[last_index] += 1
        last_index = q.get()

def matrix_listener(q, matrix):
    """ Accumulates intersections across samples.

        q: queue with sample indexes
        matrix: matrix i for which i[k][l] is the number of introns in
            common between samples k and l, or matrix u for which u[k][l] is
            the number of introns found in either sample k or sample l

        No return value.
    """
    last_record = q.get()
    while last_record != -1:
        matrix[last_record[1]][last_record[2]] += 1
        last_record = q.get()

def counts(q, lines):
    """ Finds intron counts across samples

        q: queue for accumulating results
        lines: iterable of some subset of lines input from stdin

        Return value: 0 if successful
    """
    for line in lines:
        if not line: continue
        tokens = line.strip().split('\t')
        for index in tokens[3].split(','):
            q.put(int(index))
    return 0

def intersections_and_unions(q_i, q_u, lines, forbidden_samples=set(),
                                sample_count=3000):
    """ Updates intersections/unions across samples for set of intron lines.

        q_i: queue for accumulating intersections
        q_u: queue for accumulating unions
        lines: iterable of some subset of lines input from stdin
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
            q_i.put((True, i, j))
            if i != j:
                q_i.put((True, j, i))
        for i in xrange(sample_count):
            if i in forbidden_samples: continue
            for j in found_indexes:
                q_u.put((False, i, j))
                if i != j:
                    q_u.put((False, j, i))
    return 0

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
    parser.add_argument('--input', '-i', type=str, required=True,
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
    manager = multiprocessing.Manager()
    q = manager.Queue()
    intron_counts = defaultdict(int)
    results = []
    print >>sys.stderr, '\x1b[KCounting introns...'
    listener = Thread(target=count_listener, args=(q, intron_counts))
    listener.daemon = True
    listener.start()
    with open(args.input) as input_stream:
        to_dispatch = filter(lambda x: x is not None,
                             [next(input_stream, None)
                                for _ in xrange(args.chunk_size)])
        dispatch_count, active = 0, 0
        while to_dispatch:
            len_results_before = len(results)
            pool.apply_async(counts, [q, to_dispatch], callback=results.append)
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
    q.put(-1)
    listener.join()
    process_time = time.time()
    print >>sys.stderr, '\x1b[KChunks processed in %02f s.' % (
                process_time - start_time
            )
    forbidden_samples = set(
                [int(sample_index) for sample_index in intron_counts
                    if intron_counts[sample_index] < args.filter]
            )
    print >>sys.stderr, '\x1b[K%d samples filtered out.' % len(
                                                            forbidden_samples
                                                        )
    print >>sys.stderr, '\x1b[KComputing Jaccard matrix...'
    print ';'.join(map(str, sorted(forbidden_samples)))
    intersections = [[0 for _ in xrange(args.sample_count)]
                        for __ in xrange(args.sample_count)]
    unions = [[0 for _ in xrange(args.sample_count)]
                 for __ in xrange(args.sample_count)]
    results = []
    q_i, q_u = manager.Queue(), manager.Queue()
    listener_i = Thread(target=matrix_listener, args=(q, intersections))
    listener_i.daemon = True
    listener_i.start()
    listener_u = Thread(target=matrix_listener, args=(q, unions))
    listener_u.daemon = True
    listener_u.start()
    with open(args.input) as input_stream:
        to_dispatch = filter(lambda x: x is not None,
                             [next(input_stream, None)
                                for _ in xrange(args.chunk_size)])
        dispatch_count, active = 0, 0
        while to_dispatch:
            len_results_before = len(results)
            pool.apply_async(intersections_and_unions,
                [q_i, q_u, to_dispatch, forbidden_samples, args.sample_count],
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
    q_i.put(-1)
    q_u.put(-1)
    listener_i.join()
    listener_u.join()
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