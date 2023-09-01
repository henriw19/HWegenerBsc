import networkx as nx
import matplotlib.pyplot as plt
import numpy as np
from hwbsc.stab_codes import *

def edge_matrix(n):
    rep = np.zeros((n-1,n),dtype=int)
    for i in range(n-1):
        for j in range(n):
            rep[i,i] = 1
            rep[i,i+1] = 1
    hyprep = hypergraph_prod(rep).H
    a = hyprep.shape[0]
    b = hyprep.shape[1]
    mat = np.block([[hyprep[a//2:,:b//2]],[hyprep[:a//2,b//2:]]])
    edge_matrix = np.block([[np.zeros((mat.shape[0],mat.shape[0]),dtype=int),mat],[mat.T,np.zeros((mat.shape[1],mat.shape[1]),dtype=int)]])
    return edge_matrix, hyprep

class Node:
    def __init__(self, index):
        self.index = index
        self.neighbours=set()

    def __hash__(self):
        return hash(self.index)
    
    def __eq__(self, other):
        if isinstance(other, Node):
            return self.index == other.index
        return False

    def add_neighbour(self, node):
        self.neighbours.add(node)

class Cluster:
    def __init__(self, node):

        self.nodes= {node}
        print("zallo")

    def grow(self):
        
        erweiterungsliste = []
        for node in self.nodes:
            for neighbour in node.neighbours:
                erweiterungsliste.append(neighbour)
        
        for nachbar in erweiterungsliste:
            self.nodes.add(nachbar)

def main():
    H = np.array([[1,1,0],[0,1,1]])
    hz = np.vstack((np.kron(np.eye(H.shape[0],dtype = int),H),np.kron(H.T,np.eye(H.shape[1],dtype = int))))
    matrix = np.block([[np.zeros((hz.shape[0],hz.shape[0]),dtype=int),hz],[hz.T,np.zeros((hz.shape[1],hz.shape[1]),dtype=int)]])
    G=nx.Graph(matrix)
    nx.draw(G, with_labels = True)
    plt.show()

    liste = []
    for index in range(len(matrix)):
        liste.append(Node(index))
    for node in liste:
        i = node.index
        nachbarliste = list(np.where(matrix[i]==1)[0])
        for element in nachbarliste:
            node.add_neighbour(liste[element])

    cluster = Cluster(liste[0])
    for node in cluster.nodes:
        print(node.index)
    cluster.grow()
    cluster.grow()

    for node in cluster.nodes:
        print(node.index)

if __name__ == '__main__':
    main()
