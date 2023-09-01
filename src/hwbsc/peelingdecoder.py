import numpy as np
import networkx as nx
import matplotlib.pyplot as plt

def peeling_decoder(cluster, xGr, syndrome):
    """
    Perform peeling decoding on a cluster of a networkx graph representing a quantum code.

    Parameters
    ----------
    cluster : networkx.Graph
        The graph representing the cluster of the quantum correction code.
        Lit-up syndrome nodes are flagged with the attribute [syndrome] = True

    Returns
    -------
    A : list
        A list of data nodes representing the proposed error correction based on the peeling decoding.

    """
    # print("cluster",cluster)
    # print("syndrome",syndrome)
    num_data = len(list(xGr.nodes())) - len(syndrome)

    # light up the syndrome nodes: NO
    # for i in range(num_data+len(syndrome)): 
    #     xGr.nodes[i]['syndrome'] = False
    #     if i >= num_data:
    #         if syndrome[i-num_data] == 1:
    #             xGr.nodes[i]['syndrome'] = True

    # initialize the proposed error correction
    
    A = set()
    
    cluster_to_solve = xGr.subgraph(cluster)
    spanning_tree = nx.minimum_spanning_tree(cluster_to_solve)
    # iterate over leaf nodes
    while len(list(spanning_tree.nodes)) > 1:
        us_to_remove = []
        
        for u in spanning_tree:
            v = list(spanning_tree.neighbors(u))
            if len(v) == 1:
                v = v[0]              
                # if u is a lit-up syndrome, append its edge to A, and light up neighbor node
                if spanning_tree.nodes[u]['syndrome'] == True:
                    # check if node is not a check
                    if v not in list(xGr.nodes())[-len(syndrome):]:
                        A.add(v)
                        
                    spanning_tree.nodes[u]['syndrome'] = False
                    # xGr.nodes[u]['syndrome'] = False
                    spanning_tree.nodes[v]['syndrome'] += True 
                    # xGr.nodes[v]['syndrome'] += True
                us_to_remove.append(u)
        # remove all leaf nodes after each iteration
        spanning_tree.remove_nodes_from(us_to_remove)
    #print("A",A)
    return A#,xGr






n = 5
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


# nx.draw(xGr, with_labels = True)
#plt.show()




# syndrome bei 13, 15, 16, 17
# syndrome = np.zeros(61)
# syndrome[58] = 1
# syndrome[49] = 1
# print(peeling_decoder(xGr,xGr, syndrome, 20))   