import networkx as nx
import numpy as np
from pymatching import Matching
import matplotlib.pyplot as plt
from hwbsc.color_code import color_graph, hexa_color_pcm_generator

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

# xGr = nx.Graph(matrix)

def color_union_find(facegraph,syndrome):
    error = []
    clusterlist = []
    syndrome_nodes = []


    for node in facegraph.nodes():
        facegraph.nodes[node]['syndrome'] = 'False'
    # initialize syndrome on facegraph
    for syn_node,syn in enumerate(syndrome):
        if syn == 1:
            facegraph.nodes[syn_node]['syndrome'] = 'True'
            syndrome_nodes.append(syn_node)
            clusterlist.append({syn_node})
    

    count = 1
    predicted_error = []
    numclusters = len(clusterlist)

    # clusterlist[0].update({3})
    # clusterlist.setdefault(0, set()).add(9)
    flag = False

    while count < 2:
        
        mergeset = set()
        
        # for i,(cluster1,nodesincluster1) in enumerate(clusterlist.items()):
        for i in range(numclusters):
            
            # for j,(cluster2,nodesincluster2) in enumerate(clusterlist.items()):
            for j in range(numclusters):
                if i != j and j <= len(clusterlist) and i <= len(clusterlist):
                    print("i,j:",i,j)
                    if clusterlist[i].intersection(clusterlist[j]) == set():
                        # cluster haben leere Schnittmenge: vergrößern 
                        

                        print("grow j =",j)
                        print(clusterlist)
                        clusterlist[j].update(grow(clusterlist[j], facegraph))

                        # clusterlist[cluster2].union({grow(clusterlist[cluster2], facegraph)})
                        print("j gegrow",clusterlist)

                       #### WHY
                        if j == numclusters-1:
                            print("grow i =",i)
                            print(clusterlist)
                            clusterlist[i].update(grow(clusterlist[i], facegraph))
                            print("i gegrow",clusterlist)
                        
                        
                    for c1 in range(len(clusterlist)):
                        for c2 in range(len(clusterlist)):
                            if c1 != c2:
                                if clusterlist[c1].intersection(clusterlist[c2]) != set():
                                    mergeset = [c1,c2]
                    print("mergelist",mergeset)

                    
                    # cluster schneiden sich: mergen
                    # clusterlist[i] = clusterlist[i].union(clusterlist[j])
                    # merge(clusterlist, mergelist)
                    # for merge in mergeset:
                    print("merge",mergeset[0],mergeset[1])
                    clusterlist[mergeset[0]].update(clusterlist[mergeset[1]])
                    del clusterlist[mergeset[1]]


                    count +=1
                    print("counter",count)
                    
                    print(clusterlist)

                    # check that parity of each color in the cluster is even:
                    r = 0
                    g = 0
                    b = 0
                    
                    b = i
                    if len(clusterlist) == 1:
                        b = 0
                    for nodeincluster in clusterlist[b]:
                        if facegraph.nodes[nodeincluster]['syndrome'] == 'True':
                            if facegraph.nodes[nodeincluster]['color'] == 'r':
                                r += 1
                            if facegraph.nodes[nodeincluster]['color'] == 'g':
                                g += 1
                            if facegraph.nodes[nodeincluster]['color'] == 'b':
                                b += 1 
                        # gerade parität:
                    if r%2 == 0 and g%2 == 0 and b%2 == 0:
                        print("decode even")
                    if r%2 == 1 and g%2 == 1 and b%2 == 1:
                        print("decode UNeven")
                    else:
                        pass
                        # weder noch: noch nicht decoden
             
                    

                # for i,j in mergelist:
                #     print("decode")
                #     clusterlist.remove(clusterlist[i])
                        

                        
                        

    return predicted_error

def grow(cluster, graph):
    neighbor_ancillas = set()
    for index in cluster:
        neighbor_ancillas.add(index)

    for clusterindex in cluster:
        for neighbordata in list(graph.neighbors(clusterindex)):
            ancilla = list(graph.nodes())[neighbordata]
            neighbor_ancillas.add(ancilla)
    return neighbor_ancillas

def merge(clusterlist, mergelist):
    clusterlistnew = clusterlist
    todelete = []
    for tomerge in mergelist:
        clusterlistnew[tomerge[0]].update(tomerge[1])
        todelete.append(tomerge[1])
        
    return clusterlistnew, todelete


def decode_cluster(syndrome_nodes, graph):
    connection_indices = nx.shortest_path(graph, syndrome_nodes[0], syndrome_nodes[1])
    return connection_indices

colG, facegraph, node_colors, edge_colors = color_graph(hexa_color_pcm_generator(6,4)[0])

# nx.draw(facegraph, node_color=node_colors, edge_color=edge_colors, with_labels = True)
# plt.show()


syndrome = np.zeros(24).astype(int)
syndrome[0] = 1
syndrome[1] = 1
syndrome[3] = 1
# syndrome[13] = 1
# syndrome[19] = 1
# syndrome[21] = 1
print(color_union_find(facegraph,syndrome))