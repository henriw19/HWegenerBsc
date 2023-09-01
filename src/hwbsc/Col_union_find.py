import networkx as nx
import numpy as np
import matplotlib.pyplot as plt
from hwbsc.inversion_decoder import gfinversion



# n = 3
# rep = np.zeros((n,n),dtype=int)
# for i in range(n):
#     for j in range(n):
#         rep[i,i] = 1
#         rep[i,(i+1)%n] = 1
# # H = rep
# hz = np.vstack((np.kron(np.eye(H.shape[0],dtype = int),H),np.kron(H.T,np.eye(H.shape[1],dtype = int))))
# hx = np.vstack((np.kron(H,np.eye(H.shape[0],dtype = int)),np.kron(np.eye(H.shape[1],dtype = int),H.T)))
# #######print(hz.T)
# hi = hx
# matrix = np.block([[np.zeros((hi.shape[0],hi.shape[0]),dtype=int),hi],[hi.T,np.zeros((hi.shape[1],hi.shape[1]),dtype=int)]])
# # facegraph = nx.Graph(matrix)



class Alm_lin_union_find():
       

    def __init__(self, colG, facegraph, syndrome, err = []):
        # key =  root; value = [ parity, size, set(cluster nodes), boundary list  ]
        # initialize the cluster tree
        colors = []
        
        for i in syndrome:
            if i == 0:
                colors.append('b')
            if i == 1:
                colors.append('r')
        for i in range(144):
            colors.append('b')
        
        # nx.draw(colG, with_labels = True, node_color = colors)
        # plt.show()
        facegraph = facegraph
        facegraphcopy = facegraph
               
        self.find_LUT = {}
        self.cluster_roots = {}
        self.support = {}
        self.counter = 0
        self.odd_cluster_list = []
        num_data = int(colG.number_of_nodes()-facegraph.number_of_nodes())
        
        self.decode_prediction = np.zeros(facegraph.number_of_nodes()*2).astype(int)

        num_check = len(syndrome)

        self.decode = np.zeros(num_data,dtype=int)

        # initialize find_LookUpTable
        for node in facegraph.nodes():
            facegraph.nodes[node]['syndrome'] = False
            self.find_LUT[node] = -1

        # flag syndromes in nx graph
        for i in range(len(syndrome)):
            if syndrome[i] == 1:
                facegraph.nodes[i]['syndrome'] = True
            else:
                facegraph.nodes[i]['syndrome'] = False

        for i,node in enumerate(list(facegraph.nodes())):
            if syndrome[i] == 1:
                self.cluster_roots[i] = [False, 1, {node}, {node}]
                self.find_LUT[i] = i   
                self.odd_cluster_list.append(i)    
        # nx.draw(facegraph,with_labels=True)
        # plt.show()
        for edge in facegraph.edges():
            self.support[edge] = 0
            self.support[(edge[1],edge[0])] = 0
        
        ######print(facegraph.nodes(data=True))

        # break condition: no odd-parity clusters left
        while self.odd_cluster_list:
            self.counter +=1
            if self.counter > 100:
                # print("fail")
                break

            for key, value in self.cluster_roots.items():
                if key in self.odd_cluster_list:
                    ######print(self.odd_cluster_list)
                    
                    fusion_list = []

                    # fusion list contains all new-grown neighbor edges
                    fusion_list.append(self.grow(key,facegraph))
                    ####print(self.cluster_roots)
                    for edge in fusion_list[0]:                   
                     # verify if new grown node belongs to an odd cluster:


                     
                        c1 = self.find(edge[0])
                        c2 = self.find(edge[1])
                        
                        if c1 != c2 and c1 in self.odd_cluster_list and c2 in self.odd_cluster_list:
                
                            # compare cluster sizes:
                            if self.cluster_roots[c1][1] >= self.cluster_roots[c2][1]:
                                # add c2 to c1
                                self.union(c1, c2,facegraph)
                                # if parity even: decode!
                                # self.grow(c1,facegraph)
                                # self.parity_check(facegraph,c1)
                                if self.cluster_roots[c1][0]  == True and len(self.odd_cluster_list) > 0:
                                    self.odd_cluster_list.remove(c1)
                                    

                                    ######print("decode1:",self.cluster_roots[c1][2],self.cluster_roots[c2][2])
                                    subcluster = self.cluster_roots[c1][2].union(self.cluster_roots[c2][2])
                                    ####print(self.cluster_roots)
                                    
                                    # disable syndrome colors:
                                    for ancilla in self.cluster_roots[c1][2]:
                                        facegraph.nodes[ancilla]['color'] = '0'

                                    ##print("decode", self.cluster_roots[c1][2])
                                    predo = gfinversion(colG, subcluster, facegraph, syndrome, err)
                                    self.decode_prediction ^= predo %2

                                    for ancilla in subcluster:
                                        syndrome[ancilla] = 0

                                    facegraph.nodes[c1]['syndrome'] = False
                                    facegraph.nodes[c2]['syndrome'] = False

                            elif self.cluster_roots[c1][1] < self.cluster_roots[c2][1]:
                                # add c1 to c2
                                self.union(c2, c1,facegraph)
                                # if parity even: decode!
                                # self.grow(c2,facegraph)
                                # self.parity_check(facegraph,c2)
                                if self.cluster_roots[c2][0] == True and len(self.odd_cluster_list) > 0:
                                    self.odd_cluster_list.remove(c2)
                                    ######print("decode2")
                                    subcluster = self.cluster_roots[c1][2].union(self.cluster_roots[c2][2])
                                    ####print(self.cluster_roots)
                                    ##print("decode", c2)

                                    # disable syndrome colors:
                                    for ancilla in self.cluster_roots[c2][2]:
                                        facegraph.nodes[ancilla]['color'] = '0'



                                    # self.decode_prediction.update(gfinversion(colG, subcluster, facegraph, syndrome))
                                    ##print("decode", self.cluster_roots[c2][2])
                                    predo =gfinversion(colG, subcluster, facegraph, syndrome, err)
                                    self.decode_prediction ^= predo %2
                                    
                                    for ancilla in subcluster:
                                        syndrome[ancilla] = 0
                                    

                                    facegraph.nodes[c1]['syndrome'] = False
                                    facegraph.nodes[c2]['syndrome'] = False
        
        
        
        # create output:

        
        self.decode_prediction = list(self.decode_prediction)   
        self.decode = list(self.decode_prediction)
        # self.decodebin = np.zeros(num_data,dtype=int)
        # for index in self.decode_prediction:
        #     self.decodebin[index] = 1

    def grow(self, key,facegraph):
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

        fusion_edges = set()
        neighbor_nodes = []
        cluster = self.cluster_roots[key][3]
        # the roots from where the cluster grows
        roots = []

        # iterate over the cluster boundary nodes
        for clusterbound in cluster:
            roots.append(clusterbound)
            # iterate over all neighbors
            for new_neighbors in facegraph.neighbors(clusterbound):
                neighbor = list(facegraph.nodes())[new_neighbors]
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
            # if any(self.find_LUT[num] == -1 for num in list(facegraph.neighbors(neighbor))):
            #     self.cluster_roots[key][3].add(neighbor)
        
        return fusion_edges

    def union(self, u, v, facegraph):
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
        # print("union",u,v)
        for clusterroot, nodes in self.cluster_roots.items():
            
                if clusterroot != u and clusterroot != v and clusterroot in self.odd_cluster_list:
                    if nodes[2].intersection(self.cluster_roots[u][2]) != set():
                        ####print("submerge", u, clusterroot)
                        self.cluster_roots[u][2] |= self.cluster_roots[clusterroot][2]
                        self.cluster_roots[u][1] += self.cluster_roots[clusterroot][1]
                        if clusterroot in self.odd_cluster_list:
                            self.odd_cluster_list.remove(clusterroot)

                        # update find_LUT:
                        for node in self.cluster_roots[u][2]:
                            self.find_LUT[node] = u
                        for node in self.cluster_roots[clusterroot][2]:
                            self.find_LUT[node] = u

                    elif nodes[2].intersection(self.cluster_roots[v][2]) != set():
                        ####print("submerge", v, clusterroot)
                        self.cluster_roots[v][2] |= self.cluster_roots[clusterroot][2]
                        self.cluster_roots[v][1] += self.cluster_roots[clusterroot][1]
                        if clusterroot in self.odd_cluster_list:
                            self.odd_cluster_list.remove(clusterroot)

                        # update find_LUT:
                        for node in self.cluster_roots[v][2]:
                            self.find_LUT[node] = v
                        for node in self.cluster_roots[clusterroot][2]:
                            self.find_LUT[node] = v
                
        ######print(self.cluster_roots)

        
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


        


        # # update parity
        # # check that parity of each color in the cluster is even:
        # r = 0
        # g = 0
        # b = 0
        # ######print("2",self.cluster_roots)
        
        # for nodeincluster in self.cluster_roots[u][2]:
        #     if facegraph.nodes[nodeincluster]['syndrome'] == True:
        #         if facegraph.nodes[nodeincluster]['color'] == 'r':
        #             r += 1
        #         if facegraph.nodes[nodeincluster]['color'] == 'g':
        #             g += 1
        #         if facegraph.nodes[nodeincluster]['color'] == 'b':
        #             b += 1 
        #     # gerade parität:
        # ###print("cluster",u,"r,g,b",r,g,b)
        # if r%2 == 0 and g%2 == 0 and b%2 == 0 and r+g+b >= 3:
        #     self.cluster_roots[u][0] = True
        #     ######print("gerade")
        # # elif r%2 == 0 and g%2 == 0 and b%2 == 0 and (r+g+b)%2 ==2 and (r == 0 or g == 0 or b == 0):
        # #     self.cluster_roots[u][0] = True
        # #     ######print("gerade")
        # elif r%2 == 1 and g%2 == 1 and b%2 == 1 and r+g+b >= 3:
        #     self.cluster_roots[u][0] = True
        #     ######print("ungerade")

        self.parity_check(facegraph,u)



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
    
    def parity_check(self, facegraph, u):
        # check that parity of each color in the cluster is even:
        r = 0
        g = 0
        b = 0
        ######print("2",self.cluster_roots)
        for nodeincluster in self.cluster_roots[u][2]:
            if facegraph.nodes[nodeincluster]['syndrome'] == True:
                if facegraph.nodes[nodeincluster]['color'] == 'r':
                    r += 1
                if facegraph.nodes[nodeincluster]['color'] == 'g':
                    g += 1
                if facegraph.nodes[nodeincluster]['color'] == 'b':
                    b += 1 
            # gerade parität:
        ###print("cluster",u,"r,g,b",r,g,b)
        if r%2 == 0 and g%2 == 0 and b%2 == 0 and r+g+b >= 3:
            self.cluster_roots[u][0] = True
            ######print("gerade")
        # elif r%2 == 0 and g%2 == 0 and b%2 == 0 and (r+g+b)%2 ==2 and (r == 0 or g == 0 or b == 0):
        #     self.cluster_roots[u][0] = True
        #     ######print("gerade")
        elif r%2 == 1 and g%2 == 1 and b%2 == 1 and r+g+b >= 3:

            self.cluster_roots[u][0] = True
            ######print("ungerade")


# # output = hexa_color_pcm_generator(12,6)[0][60]
# pcm = hexa_color_pcm_generator(12,6)[0]
# err = np.zeros(144).astype(int)


# # ### 400.000 combinations:
# err[0] = 1  # 72
# err[7] = 1  #79
# err[122] = 1    #194
# err[123] = 1 #195
# err[140] = 1 #212
# err[2] = 1 #74

# err[55] = 1
# err[72] = 1
# # err[89] = 1
# # err[91] = 1
# # err[12] = 1
# # err[23] = 1
# err[45] = 1
# err[34] = 1
# #########

# [0,2,7,32,33,50]
# #small: [0,2,5,22]








# #1: 90,114
# #2: 110, 113
# #3: 161, 162

# # 0 0 1 0 0 0 0 0 0 0 1 0 0 0




# err = [0, 1, 1, 0, 1 ,1 ,0 ,1 ,0 ,0 ,0 ,0 ,0 ,0 ,1 ,0 ,0 ,1 ,0 ,0 ,1 ,0 ,1 ,0 ,1 ,0 ,1 ,0 ,0 ,1 ,1, 0, 1, 0, 0 ,1 ,1,
#  1, 0, 1, 0, 1, 0, 1, 1, 1, 1, 1, 1, 0, 1, 0, 0, 1, 1, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 ,1 ,0 ,1
# , 0, 0, 1, 0, 1, 0, 0, 0 ,0 ,0 ,0, 0, 1, 0, 0, 1, 0 ,1 ,0 ,1, 0, 0 ,0, 0, 1 ,1 ,1 ,0 ,1, 0, 0, 1, 0, 1, 1, 0, 1,
#  1 ,0, 0, 0, 0, 0, 0, 0, 1, 1 ,0, 0, 0, 0, 0, 1 ,0, 1, 0, 0, 1, 1 ,0, 0, 1 ,1 ,0 ,1 ,0, 1, 0, 0, 0]
# ogsyndrome = pcm @ err %2
# # ######print("syndrome",syndrome.shape,"pcm",pcm.shape,"err",err.shape)
# colG, facegraph, node_colors, edge_colors = color_graph(pcm)
# ####print(  "COLORS"   ,node_colors)

# #print("syndrome",ogsyndrome) #11,12.21.22.23.24,45,71
# bbbb = []
# for i,syn in enumerate(ogsyndrome):
#     if syn == 1:
#         bbbb.append(i)
# #print("syndromes",bbbb)
# prediction = Alm_lin_union_find(colG,facegraph,ogsyndrome,err).decode
# # #print(pcm @ prediction%2,"\n", pcm @ err %2)
# if np.array_equal(pcm @ prediction%2, pcm @ err %2):
    
#     print("Success")   
#     pass 
# else:
#     #print(colG.nodes())
#     print("Failure")
#     pass

















# # # # lit up ancillas: 13, 14, 15, 16, 17, 18
# # # #syndrome = np.array([1,0,1,1,0,0])
# # # # logical = np.array([1,0,0,1,0,0,1,0,0,0,0,0,0,0,0,0,0,0])
# # # # import itertools
# # # # for i in range(512,2**18):
# # # #     ######print(i/2**18)
# # # #     error = np.array(list(itertools.product([0,1], repeat = 18))[i])
# # # #     syndrome = hx.T @ error.T %2
# # # #     if not np.array_equal(hx.T@Alm_lin_union_find(hx, syndrome).decodebin %2, syndrome):
# # # #         ######print(hx.T@Alm_lin_union_find(hx, syndrome).decodebin, syndrome)
# # # #         break
    
# # #     #result = [0 0 0 0 0 0 0 0 1 0 0 0 0 0 0 0 0 1]

# # #     # if not np.array_equal(syndrome,hx.T @ (Alm_lin_union_find(hx, syndrome).decodebin) %2):
# # #     #     ######print()
# # #         # ######print(i,"\n" ,syndrome,"\n", hx.T @ (Alm_lin_union_find(hx, syndrome).decodebin) %2)
# # #         # ######print("error",error)
# # #         # ######print(Alm_lin_union_find(hx, syndrome).decode)
# # #         # residual =(error + Alm_lin_union_find(hx, syndrome).decodebin) %2
# # #         # ######print("residual", residual)
# # #         # ######print( hx.T@residual %2)
        

# # # # error = np.array([1, 0, 0, 1, 0, 0, 0 ,1 ,1 ,0 ,0 ,0 ,0 ,0 ,0 ,0, 0, 1])
# # # # ######print("err",error)
# # # # syndrome = hx.T @ error.T %2

# # # # # # syndrome = np.zeros(20)
# # # # # # # syndrome[61-51] = 1
# # # # # # # syndrome[61-55] = 1
# # # # # # # # syndrome[57] = 1
# # # # # # syndrome[59] = 1
# # # # decode = Alm_lin_union_find(hx, syndrome).decodebin
# # # # ######print("decode", decode)
# # # # ######print("decode syn", hx.T@decode%2 )
# # # # # residual =(error + Alm_lin_union_find(hx, syndrome).decodebin) %2
