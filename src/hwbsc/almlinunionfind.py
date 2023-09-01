import networkx as nx
import numpy as np
import matplotlib.pyplot as plt
from hwbsc.peelingdecoder import peeling_decoder


n = 3
rep = np.zeros((n,n),dtype=int)
for i in range(n):
    for j in range(n):
        rep[i,i] = 1
        rep[i,(i+1)%n] = 1
H = rep
hz = np.vstack((np.kron(np.eye(H.shape[0],dtype = int),H),np.kron(H.T,np.eye(H.shape[1],dtype = int))))
hx = np.vstack((np.kron(H,np.eye(H.shape[0],dtype = int)),np.kron(np.eye(H.shape[1],dtype = int),H.T)))
#print(hz.T)
hi = hx
matrix = np.block([[np.zeros((hi.shape[0],hi.shape[0]),dtype=int),hi],[hi.T,np.zeros((hi.shape[1],hi.shape[1]),dtype=int)]])
# xGr = nx.Graph(matrix)



class Alm_lin_union_find():
       

    def __init__(self, hi, syndrome):
        # key =  root; value = [ parity, size, set(cluster nodes), boundary list  ]
        # initialize the cluster tree
        
        matrix = np.block([[np.zeros((hi.shape[0],hi.shape[0]),dtype=int),hi],[hi.T,np.zeros((hi.shape[1],hi.shape[1]),dtype=int)]])
        self.xGr = nx.Graph(matrix)
               
        self.find_LUT = {}
        self.cluster_roots = {}
        self.support = {}
        self.counter = 0
        self.odd_cluster_list = []
        self.decode_prediction = set()
        num_data = hi.shape[0]
        num_check = len(syndrome)

        if sum(syndrome) == 0:
            self.decodebin = np.zeros(num_data,dtype=int)


        # initialize find_LookUpTable
        for node in self.xGr.nodes():
            self.xGr.nodes[node]['syndrome'] = False
            self.find_LUT[node] = -1

        # flag syndromes in nx graph
        for i in range(len(syndrome)):
            if syndrome[i] == 1:
                self.xGr.nodes[i+num_data]['syndrome'] = True
            else:
                self.xGr.nodes[i+num_data]['syndrome'] = False

        for i,node in enumerate(list(self.xGr.nodes())[-len(syndrome):]):
            if syndrome[i] == 1:
                self.cluster_roots[i+ num_data] = [1, 1, {node}, {node}]
                self.find_LUT[i+num_data] = i + num_data    
                self.odd_cluster_list.append(i+num_data)    
        # nx.draw(self.xGr,with_labels=True)
        plt.show()
        for edge in self.xGr.edges():
            self.support[edge] = 0
            self.support[(edge[1],edge[0])] = 0
        

        print(self.xGr.nodes(data=True))

        # break condition: no odd-parity clusters left
        while self.odd_cluster_list:
            self.counter +=1
            if self.counter > 30:
                print("fail")

            for key, value in self.cluster_roots.items():
                if key in self.odd_cluster_list:
                    print(self.odd_cluster_list)
                    fusion_list = []

                    # fusion list contains all new-grown neighbor edges
                    fusion_list.append(self.grow(key))

                    for edge in fusion_list[0]:                   
                     # verify if new grown node belongs to an odd cluster:
                        
                        c1 = self.find(edge[0])
                        c2 = self.find(edge[1])
                        
                        if c1 != c2 and c1 in self.odd_cluster_list and c2 in self.odd_cluster_list:
                
                            # compare cluster sizes:
                            if self.cluster_roots[c1][1] >= self.cluster_roots[c2][1]:
                                # add c2 to c1
                                self.union(c1, c2)
                                # if parity even: decode!
                                if self.cluster_roots[c1][0] % 2 == 0 and len(self.odd_cluster_list) > 0:
                                    self.odd_cluster_list.remove(c1)
                                    
                                    self.decode_prediction.update(peeling_decoder(self.cluster_roots[c1][2], self.xGr, syndrome))
                                    # self.xGr = peeling_decoder(self.cluster_roots[c1][2], self.xGr, syndrome)[1]
                                    self.xGr.nodes[c1]['syndrome'] = False
                                    self.xGr.nodes[c2]['syndrome'] = False

                            elif self.cluster_roots[c1][1] < self.cluster_roots[c2][1]:
                                # add c1 to c2
                                self.union(c2, c1)
                                # if parity even: decode!
                                if self.cluster_roots[c2][0] % 2 == 0 and len(self.odd_cluster_list) > 0:
                                    self.odd_cluster_list.remove(c2)
                                    self.decode_prediction.update(peeling_decoder(self.cluster_roots[c2][2], self.xGr, syndrome))
                                    self.xGr.nodes[c1]['syndrome'] = False
                                    self.xGr.nodes[c2]['syndrome'] = False
        # create output:
        self.decode_prediction = list(self.decode_prediction)   
        self.decode = list(self.decode_prediction)
        self.decodebin = np.zeros(num_data,dtype=int)
        for index in self.decode_prediction:
            self.decodebin[index] = 1

    def grow(self, key):
        """
        Increase the radius of the clusters by one node, 
        fuse the clusters that meet and update the clusters 
        and their boundary lists.

        Parameters:
        -----------
        key : np.ndarray
            Root of the clusters to be grown

        Returns:
        --------
        merged_cluster : ???
            The merged cluster


        """        
        #print("grow",key)

        fusion_edges = set()
        neighbor_nodes = []
        cluster = self.cluster_roots[key][3]
        # the roots from where the cluster grows
        roots = []

        # iterate over the cluster boundary nodes
        for clusterbound in cluster:
            roots.append(clusterbound)
            # iterate over all neighbors
            for new_neighbors in self.xGr.neighbors(clusterbound):
                neighbor = list(self.xGr.nodes())[new_neighbors]
                edge = (clusterbound,neighbor)
                # verify if the edge has not been grown already and if the new node belongs to an odd cluster
                # if both is True, then clusters should be fused
                if (self.support[edge] == 0) and self.find(edge[1]) in self.odd_cluster_list : #(self.find(edge[1]) != key and self.find(edge[1]) != -1):
                    fusion_edges.add(edge)

                # if False, add the new node to the cluster in the find_LUT
                else:
                    self.find_LUT[edge[1]] = key
                
                # either way, mark the new edges as grown in support
                self.support[edge] = 1
                neighbor_nodes.append(new_neighbors)
        
        self.cluster_roots[key][2].update(neighbor_nodes)
        self.cluster_roots[key][1] = len(self.cluster_roots[key][2])

        # update cluster boundary list
        for neighbor in neighbor_nodes:
            for root in roots:
                if root in self.cluster_roots[key][3]:
                    self.cluster_roots[key][3].remove(root)
            if self.find(neighbor) in self.odd_cluster_list or self.find(neighbor) == -1:
                self.cluster_roots[key][3].add(neighbor)
            # if any(self.find_LUT[num] == -1 for num in list(self.xGr.neighbors(neighbor))):
            #     self.cluster_roots[key][3].add(neighbor)
        
        return fusion_edges

    def union(self, u, v):
        """
        Perform a fusion operation on the clusters with indices u and v

        Parameters:
        -----------
        u : ???
            The root of the first cluster to be merged
        v : ???
            The root of the second cluster to be merged

        Returns:
        --------
        merged_cluster : ???
            The merged cluster


        """
        #weighted union
        self.cluster_roots[u][2] |= self.cluster_roots[v][2]
        self.cluster_roots[u][1] += self.cluster_roots[v][1]
        if v in self.odd_cluster_list:
            self.odd_cluster_list.remove(v)

        # update find_LUT:
        for node in self.cluster_roots[u][2]:
            self.find_LUT[node] = u
        for node in self.cluster_roots[v][2]:
            self.find_LUT[node] = u

        # update parity
        self.cluster_roots[u][0] += self.cluster_roots[v][0] % 2

        return self.cluster_roots




    def find(self, v):
        """
        Identify the cluster to which node v belongs.

        Parameters:
        -----------
        v : int
            The node of which we want to find the cluster
        
        Returns:
        --------
        cluster_index : int
            The index of the cluster to which v belongs

        """
        return self.find_LUT[v]
    

# lit up ancillas: 13, 14, 15, 16, 17, 18
#syndrome = np.array([1,0,1,1,0,0])
# logical = np.array([1,0,0,1,0,0,1,0,0,0,0,0,0,0,0,0,0,0])
# import itertools
# for i in range(512,2**18):
#     print(i/2**18)
#     error = np.array(list(itertools.product([0,1], repeat = 18))[i])
#     syndrome = hx.T @ error.T %2
#     if not np.array_equal(hx.T@Alm_lin_union_find(hx, syndrome).decodebin %2, syndrome):
#         print(hx.T@Alm_lin_union_find(hx, syndrome).decodebin, syndrome)
#         break
    
    #result = [0 0 0 0 0 0 0 0 1 0 0 0 0 0 0 0 0 1]

    # if not np.array_equal(syndrome,hx.T @ (Alm_lin_union_find(hx, syndrome).decodebin) %2):
    #     print()
        # print(i,"\n" ,syndrome,"\n", hx.T @ (Alm_lin_union_find(hx, syndrome).decodebin) %2)
        # print("error",error)
        # print(Alm_lin_union_find(hx, syndrome).decode)
        # residual =(error + Alm_lin_union_find(hx, syndrome).decodebin) %2
        # print("residual", residual)
        # print( hx.T@residual %2)
        

# error = np.array([1, 0, 0, 1, 0, 0, 0 ,1 ,1 ,0 ,0 ,0 ,0 ,0 ,0 ,0, 0, 1])
# print("err",error)
# syndrome = hx.T @ error.T %2

# # # syndrome = np.zeros(20)
# # # # syndrome[61-51] = 1
# # # # syndrome[61-55] = 1
# # # # # syndrome[57] = 1
# # # syndrome[59] = 1
# decode = Alm_lin_union_find(hx, syndrome).decodebin
# print("decode", decode)
# print("decode syn", hx.T@decode%2 )
# # residual =(error + Alm_lin_union_find(hx, syndrome).decodebin) %2
