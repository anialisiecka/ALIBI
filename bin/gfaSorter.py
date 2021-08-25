import os
import sys
from gfaReader import read
from graph import Graph

def strand(x):
    if x == '+': return 1
    else: return -1

def reverseComplement(seq):
    d = {'A':'T', 'T':'A', 'C':'G', 'G':'C', 'a':'t', 't':'a', 'c':'g', 'g':'c'}
    rc=''
    for i in range(len(seq)-1, -1, -1):
        rc+=d[seq[i]]
    return rc

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

def write(infile, outfile, blocks):
    links=[]
    paths=[]
    with open(infile) as f:
        for line in f:
            if line.strip().startswith('L'):
                links.append(line.split()[1:])
            elif line.strip().startswith('P'):
                paths.append(line.split()[1:])
    sortedBlocks=sorted(blocks, key=lambda b: b.order()) # sorted list of nodes
    with open(outfile, 'w') as f:
        f.write('H'+'\t'+'VN:Z:1.0'+'\n')
        for b in sortedBlocks:
            if b.orientation() == -1: b.seq=reverseComplement(b.seq)
            f.write('S'+'\t'+str(b.id+1)+'\t'+b.seq+'\n')
        for p in paths:
            f.write('P'+'\t'+p[0]+'\t')
            x=p[1].split(',')
            cigars = p[2].split(',')
            for i in range(len(x)):
                v = x[i][:-1]
                if blocks[int(v)-1].orientation() == strand(x[i][-1]):
                    x[i] = v + '+'
                else:
                    x[i] = v + '-'
                if blocks[int(v)-1].orientation() == -1:
                    cigars[i] = reverseCigar(cigars[i])
            f.write(','.join(x)+'\t'+','.join(cigars)+'\n')
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
            f.write('L'+'\t'+'\t'.join(l)+'\n')
            
def reorder(R_f, R_b, blocks):
    R_f.sort(key = lambda x: blocks[x].order())
    R_b.sort(key = lambda x: blocks[x].order())
    L = R_b + R_f
    O = sorted(blocks[x].order() for x in L)
    for i in range(len(L)):
        blocks[L[i]].reorder(O[i])

def addEdgeWithinComponent(e, w, G, blocks):
    x, y = e[0], e[2]
    lb = blocks[y].order()
    ub = blocks[x].order()
    if lb is ub: G.fa += w
    elif lb < ub:
        R_f = G.dfsF(y, blocks, ub)
        if R_f:
            R_b = G.dfsB(x, blocks, lb)
            reorder(R_f, R_b, blocks)
            G.addEdge(x,y)
        else:  G.fa.add(e)
    else: G.addEdge(x,y)
    
def addEdgeBetweenComponents(e, blocks):
    # Add edge between blocks of different components
    reverse, flank = 1, 1
    if blocks[e[0]].size() < blocks[e[2]].size():
        e = e[2:] + e[:2]
    if blocks[e[0]].orientation()*blocks[e[2]].orientation() is e[1]*e[3]:
        reverse = -1
    if blocks[e[0]].orientation()*e[1] < 0:
        flank = -1
    blocks[e[2]].unionto(blocks[e[0]], reverse, flank)

def connect_components(blocks):
    if len(blocks) != blocks[0].size(): 
        d = {blocks[0].find(): 0}
        n = blocks[0].maximum()
        for block in blocks:
            if block.find() not in d:
                i = n - block.minimum() + 1
                d[block.find()] = i 
                block.reorder(i + block.order())
                n += block.size()
            else:
                block.reorder(d[block.find()] + block.order())

def linSort(infile): 
    # blocks - list of Block instances
    # edges - list of edges sorted by the weight
    blocks, edges = read(infile)
    G = Graph(len(blocks))
    outfile = os.path.splitext(infile)[0]+"_sorted.gfa"
    for e, w in edges:
        if blocks[e[0]].find() is blocks[e[2]].find():
            if blocks[e[0]].orientation()*blocks[e[2]].orientation() is e[1]*e[3]:
                G.rj.add(e)
                continue
            elif blocks[e[0]].orientation()*e[1] > 0:
                addEdgeWithinComponent(e, w, G, blocks)
            else:
                addEdgeWithinComponent(e[2:]+e[:2], w, G, blocks)
        else:
            addEdgeBetweenComponents(e, blocks)
            G.addEdge(e[0], e[2])
    connect_components(blocks)
    write(infile, outfile, blocks)

def main():
    infile = sys.argv[1]
    linSort(infile)

if __name__ == "__main__":
    main()
