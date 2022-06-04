import sys

def metrics(infile):

    n = 0 # number of nodes
    order = {}
    wfa = 0.0 # weighted feedback arcs
    wrj = 0.0 # weighted reversing joins

    with open(infile) as f:
        for line in f:
            line = line.strip().split()
            if line[0] == "S":
                order[line[1]] = n
                n += 1
            elif line[0] == "P":
                path=line[2].split(',')
                for i in range(len(path)-1):
                    if order[path[i][:-1]] >= order[path[i+1][:-1]] and path[i][-1]==path[i+1][-1]=='+':
                        wfa+=1
                    elif order[path[i][:-1]] <= order[path[i+1][:-1]] and path[i][-1]==path[i+1][-1]=='-':
                        wfa+=1
                    elif path[i][-1]!=path[i+1][-1]:
                        wrj+=1

    print ("wfa: ", wfa)
    print ("wrj: ", wrj)

def main():
    infile = sys.argv[1]
    metrics(infile)

if __name__ == "__main__":
    main()
