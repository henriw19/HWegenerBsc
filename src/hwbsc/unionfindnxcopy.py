import networkx as nx
import numpy as np
from pymatching import Matching
import matplotlib.pyplot as plt
from hwbsc.peelingdecoder import *
from hwbsc.stabilizers_utils import find_log_op



n = 3
rep = np.zeros((n-1,n),dtype=int)
for i in range(n-1):
    for j in range(n):
        rep[i,i] = 1
        rep[i,i+1] = 1
H = rep
hz = np.vstack((np.kron(np.eye(H.shape[0],dtype = int),H),np.kron(H.T,np.eye(H.shape[1],dtype = int))))
hx = np.vstack((np.kron(H,np.eye(H.shape[0],dtype = int)),np.kron(np.eye(H.shape[1],dtype = int),H.T)))
hi = hx
matrix = np.block([[np.zeros((hi.shape[0],hi.shape[0]),dtype=int),hi],[hi.T,np.zeros((hi.shape[1],hi.shape[1]),dtype=int)]])

xGr = nx.Graph(matrix)

def union_decoder(matrix,syndrome):
    
    error = []
    clusterlist = []
    syndrome_nodes = []

    clustergraph = nx.Graph()

    xGr = nx.Graph(matrix)

    for i in range(len(syndrome)):
        if syndrome[i] == 1:
            xGr.nodes[i]['syndrome'] = True
        else:
            xGr.nodes[i]['syndrome'] = False

    for i in range(len(syndrome)):
        if syndrome[i] == 1:
            syndrome_nodes.append(i)
            clusterlist.append({list(xGr.nodes())[i]})

    #len = len(clusterlist)
    counter = 0
    predicted_error = []
    
    
    while clusterlist:
        mergelist = []
        length = len(clusterlist)
        for i in range(length):
            for j in range(length):
                if i != j:
                    if clusterlist[i].intersection(clusterlist[j]) != set():   
                        # cluster schneiden sich: mergen
                        if [i,j] not in mergelist and [j,i] not in mergelist:
                            mergelist.append([i,j])
                        
      
        if mergelist != []:
            # merge cluster aus mergelist                       
            for l,clusters_to_merge in enumerate(mergelist):
                clusterlist[clusters_to_merge[0]] = clusterlist[clusters_to_merge[0]].union(clusterlist[clusters_to_merge[1]])

                clusterlist.remove(clusterlist[clusters_to_merge[1]])

                # check that the number of ancillas in the merged cluster is even
                if len(np.intersect1d(list(clusterlist[l]), syndrome_nodes)) == 2:
                    # einfach verbinden      
                    predicted_error.append(decode_cluster(np.intersect1d(list(clusterlist[l]), syndrome_nodes),xGr))
                    
                    # re-create a graph from the cluster   
                    #subgraph = xGr.subgraph(list(clusterlist[l]))
                    # call peeling decoder:
                    #peeling_decoder(subgraph)
                    
                            



                    clusterlist.remove(clusterlist[l])

                elif len(np.intersect1d(list(clusterlist[l]), syndrome_nodes)) % 2 == 0:
                    #pymatching
                    sub = xGr.subgraph(clusterlist[i])
                    matching = Matching(xGr)
                    prediction = matching.decode_to_edges_array(syndrome)
                    predicted_error.append(prediction)
                    
        # grow all:
        for i in range(len(clusterlist)):
            clusterlist[i] = grow(clusterlist[i], xGr)

    bin_pred_err = np.zeros((len(predicted_error),matrix.shape[0]),dtype=int)

    for k,pred_error in enumerate(predicted_error):
        for t in range(H.shape[0]):
            if t in pred_error:
                bin_pred_err[k][t] = 1
   
    return predicted_error

def grow(cluster, graph):
    neighbor_indices = set()
    for index in cluster:
        neighbor_indices.add(index)

    for clusterindex in cluster:
        for node_index in xGr.neighbors(clusterindex):
            node = list(xGr.nodes())[node_index]
            neighbor_indices.add(node)
    return neighbor_indices

def decode_cluster(syndrome_nodes, graph):
    connection_indices = nx.shortest_path(graph, syndrome_nodes[0], syndrome_nodes[1])
    return connection_indices

# syndrome bei 0 und 6, 12, 15
syndrome = np.array([1,0,0,0,0,0,1,0,0,0,0,0,1,0,0,1,0,0,0])
print(union_decoder(matrix,syndrome))



nx.draw_networkx(xGr)
#plt.savefig('delplot')
