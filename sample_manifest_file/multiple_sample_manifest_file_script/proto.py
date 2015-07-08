"""
proto.py

Generates Redis protocol.
Reads from stdin, writes to stdout.
"""
#!/usr/bin/env python
import sys

for line in sys.stdin:
    args = line.rstrip().split('\t')
    sys.stdout.write('\r\n'.join(['*{}'.format(len(args))] + [
                                '${}\r\n{}'.format(len(arg), arg)
                                for arg in args
                            ] + ['\r\n']))
