from collections import defaultdict

class Graph:
    def __init__(self, n):
        self.size=n     # number of vertices
        self.graph=defaultdict(list)    # default dictionary to store graph
        self.fa=set()   # set of feedback arcs
        self.rj=set()    # set of reversing joins

    def getSize(self):
        return self.size

    def getAdjacent(self, v):
        return self.graph[v]

    # add edge to the graph
    def addEdge(self, u, v):
        self.graph[u].append(v)
        self.graph[v].append(u)

    # Return a list of nodes connected to node v by edge leaving the vertex v
    def outEdges(self, v, blocks, ub):
        return [u for u in self.graph[v] if ub >= blocks[u].order() > blocks[v].order()]
    # Return a list of nodes connected to node v by edge coming in to the vertex v
    def inEdges(self, v, blocks, lb):
        return [u for u in self.graph[v] if lb < blocks[u].order() < blocks[v].order()]

    # use DFS to find all vertices connected with root and whose order is between
    # the order of root and ub
    def dfsF(self, root, blocks, ub):
        visited = set()
        stack = [root,]
        while stack:
            node = stack.pop()
            if blocks[node].order()==ub: return []
            if node not in visited:
                visited.add(node)
                stack.extend([x for x in self.outEdges(node, blocks, ub) if x not in visited])
        return list(visited)
    
    # use DFS to find all vertices connected with root and whose order is between
    # lb and the order of root
    def dfsB(self, root, blocks, lb):
        visited = set()
        stack = [root,]
        while stack:
            node = stack.pop()
            if node not in visited:
                visited.add(node)
                stack.extend([x for x in self.inEdges(node, blocks, lb) if x not in visited])
        return list(visited)
