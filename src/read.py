from Bio import AlignIO
from collections import defaultdict
from block import Block

class SequenceInfo:
    def __init__(self, block_id, seq_id, start_pos, strand):
        self.block = block_id
        self.id = seq_id
        self.start = start_pos
        self.strand = strand

def start_position(sequence):
    # Return start position relative to the plus strand
    if sequence.annotations["strand"] == -1:
        return sequence.annotations["srcSize"] - sequence.annotations["start"] - sequence.annotations["size"]
    else:
        return sequence.annotations["start"]

def connect_maf_blocks(seq1, seq2, d):
    e = (seq1.block, seq1.strand, seq2.block, -seq2.strand)
    if seq1.block > seq2.block:
        e = e[2:] + e[:2]
    elif seq1.block == seq2.block:
        e = (e[0], max(e[1], e[3]), e[2], min(e[1], e[3]))
    d[e] += 1

def weight_maf(seq):
    d = defaultdict(int)
    for i in range(len(seq) - 1):
        if seq[i].id == seq[i+1].id:
            connect_maf_blocks(seq[i], seq[i+1], d)
    edges = sorted([(x, d[x]) for x in d], key=lambda e: e[1], reverse=True)
    return edges

def read_maf(maf_file):
    blocks, seq = [], []
    for i, mafblock in enumerate(AlignIO.parse(maf_file, "maf")):
        blocks.append(Block(i, mafblock))
        for sequence in mafblock:
            seq.append(SequenceInfo(i, sequence.id, start_position(sequence), sequence.annotations["strand"]))
    seq.sort(key = lambda s: (s.id, s.start))
    edges = weight_maf(seq)
    return blocks, edges

##############

def strand(x):
    if x == '+': return 1
    else: return -1

def connect_gfa_nodes(v1, v2):
    e = (v1[0], v1[1], v2[0], -v2[1])
    if e[0] > e[2]:
        e = e[2:] + e[:2]
    elif e[0] is e[2]:
        e = (e[0], max(e[1], e[3]), e[0], min(e[1], e[3]))
    return e

def weight_gfa(line, s):
    for i in range(len(line) - 1):
        v1 = (int(line[i][:-1]) - 1, strand(line[i][-1]))
        v2 = (int(line[i+1][:-1]) - 1, strand(line[i+1][-1]))
        e = connect_gfa_nodes(v1, v2)
        s[e] += 1

def read_gfa(gfa_file):
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
                weight_gfa(line, s)
    edges=sorted([(x, s[x]) for x in s], key=lambda e: e[1], reverse=True)
    return blocks, edges
