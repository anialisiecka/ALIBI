from Bio import AlignIO
from Bio.Seq import Seq
from Bio.AlignIO.MafIO import MafWriter

def write_maf(outfile, blocks):
    min_order = blocks[0].minimum() # min order
    records = [None]*len(blocks)
    for i, block in enumerate(blocks):
        block.orient_maf_block()
        records[blocks[i].order()-min_order] = block.alignment
    handle = open(outfile, 'w')
    mw = MafWriter(handle)
    mw.write_header()
    for multiple_alignment in records:
        mw.write_alignment(multiple_alignment)
    handle.close()

#######

def strand(x):
    if x == '+': return 1
    else: return -1

def reverseCigar(cigar):
    new_cigar = []
    x = ''
    for i in range(len(cigar)):
        if cigar[i].isdigit():
            x += cigar[i]
        else:
            x += cigar[i]
            new_cigar.append(x)
            x = ''
    new_cigar.reverse()
    return ''.join(new_cigar)

def write_gfa(infile, outfile, blocks):
    links = []
    paths = []
    with open(infile) as f:
        for line in f:
            if line.strip().startswith('L'):
                links.append(line.split()[1:])
            elif line.strip().startswith('P'):
                paths.append(line.split()[1:])
    sortedBlocks = sorted(blocks, key=lambda b: b.order()) # sorted list of nodes
    with open(outfile, 'w') as f:
        f.write('H'+'\t'+'VN:Z:1.0'+'\n')
        for b in sortedBlocks:
            if b.orientation() == -1:
                b.alignment = str(Seq(b.alignment).reverse_complement())
            f.write('S' + '\t' + str(b.id+1) + '\t' + b.alignment + '\n')
        for p in paths:
            f.write('P'+'\t'+p[0]+'\t')
            x = p[1].split(',')
            cigars = p[2].split(',')
            for i in range(len(x)):
                v = x[i][:-1]
                if blocks[int(v)-1].orientation() == strand(x[i][-1]):
                    x[i] = v + '+'
                else:
                    x[i] = v + '-'
                if cigars[0] != '*' and blocks[int(v)-1].orientation() == -1:
                    cigars[i] = reverseCigar(cigars[i])
            f.write(','.join(x) + '\t' + ','.join(cigars) + '\n')
        d = {1:'+', -1: '-'}
        for l in links:
            if blocks[int(l[0])-1].order() <= blocks[int(l[2])-1].order():
                l[1] = d[blocks[int(l[0])-1].orientation()*strand(l[1])]
                l[3] = d[blocks[int(l[2])-1].orientation()*strand(l[3])]
            else:
                s1 = d[-1*blocks[int(l[2])-1].orientation()*strand(l[3])]
                s2 = d[-1*blocks[int(l[0])-1].orientation()*strand(l[1])]
                l[0], l[2] = l[2], l[0]
                l[1], l[3] = s1, s2
        links.sort(key=lambda x: (blocks[int(x[0])-1].order(), blocks[int(x[2])-1].order()))
        for l in links:
            f.write('L' + '\t' + '\t'.join(l) + '\n')
