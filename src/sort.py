import os
import sys
from read import read_maf, read_gfa
from write import write_maf, write_gfa
from graph import Graph


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
    if lb is ub: G.fa.add(e)
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
    
    if flank == 1 and blocks[e[0]].order() == blocks[e[0]].maximum():
        blocks[e[2]].unionto(blocks[e[0]], reverse, flank)
    elif flank == -1 and blocks[e[0]].order() == blocks[e[0]].minimum():
        blocks[e[2]].unionto(blocks[e[0]], reverse, flank)
    else:
        c_end = 0
        c_midst = 0
        if flank==1:
            c_end += blocks[e[0]].maximum() - blocks[e[0]].order()
        else:
            c_end += blocks[e[0]].order() - blocks[e[0]].minimum()

        n = blocks[e[0]].order()
        k = blocks[e[2]].size()
        for u in G.graph[e[0]]:
            c_midst += k

        if c_end <= c_midst:
            blocks[e[2]].unionto(blocks[e[0]], reverse, flank)
        else:
            blocks[e[2]].uniontoMidst(blocks[e[0]], reverse, flank, blocks)

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

def format_file(infile):
    with open(infile) as f:
        first_line = f.readline().strip()
        if first_line.startswith('##maf'):
            return 'maf'
        elif first_line.startswith('H'):
            return 'gfa'
        else:
            return 'bad'

def linSort(infile):  
    # blocks - list of Block instances
    # edges - list of edges sorted by the weight
    blocks, edges = [], []
    in_format = format_file(infile)
    if in_format == 'maf':
        blocks, edges = read_maf(infile)
    elif in_format == 'gfa':
        blocks, edges = read_gfa(infile)
    else:
        print("Only GFA and MAF file are accepted.")    
        return
    
    G = Graph(len(blocks))
    outfile = os.path.splitext(infile)[0] + '_sorted.' + in_format
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
    if in_format == 'maf':
        write_maf(outfile, blocks)
    elif in_format == 'gfa':
        write_gfa(infile, outfile, blocks)

def main():
    infile = sys.argv[1]
    linSort(infile)

if __name__ == "__main__":
    main()
