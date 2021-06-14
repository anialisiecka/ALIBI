# -*- coding: utf-8 -*-

from collections import defaultdict
from block import Block

def strand(x):
    if x == '+': return 1
    else: return -1

def connect_nodes(v1, v2):
    e = (v1[0], v1[1], v2[0], -v2[1])
    if e[0] > e[2]: e=e[2:]+e[:2]
    elif e[0] is e[2]: e=(e[0], max(e[1], e[3]), e[0], min(e[1], e[3]))
    return e

def function1(line, s):
    for i in range(len(line) - 1):
        v1 = (int(line[i][:-1]) - 1, strand(line[i][-1]))
        v2 = (int(line[i+1][:-1]) - 1, strand(line[i+1][-1]))
        e = connect_nodes(v1, v2)
        s[e] += 1

def read(gfa_file):
    blocks = []
    s = defaultdict(int)
    with open(gfa_file) as f:
        i = 0
        for line in f:
            if line.strip().startswith('S'):
                blocks.append(Block(i,line.split()[2]))
                i += 1
            elif line.strip().startswith('P'):
                line = line.split()[2].split(',')
                function1(line, s)
    edges=sorted([(x, s[x]) for x in s], key=lambda e: e[1], reverse=True)
    return blocks, edges
