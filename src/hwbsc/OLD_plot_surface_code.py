import networkx as nx
import numpy as np
import math
import matplotlib.pyplot as plt
from hwbsc.stab_codes import *
from hwbsc. surface_code import *
# from hwbsc.decoders import *

#### create linear code
n = 5

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
    
####



#### test error
error = np.zeros(hyprep.shape[1])
error[12] = 1
error[13] = 1
error[17] = 1
error[24] = 1
error_indices = np.where(error == 1)[0]
syndrome = hyprep@error%2
print(syndrome)
syndrome_indices = np.where(syndrome == 1)[0]
print(syndrome_indices)
####



#### logicals
logicals = StabCode(hyprep).logicals
logicalx = logicals[0][:b//2]
logicaly = logicals[1][b//2:]
logicalx_indices = np.where(logicalx == 1)[0]
logicaly_indices = np.where(logicaly == 1)[0]
####


#### 
show_syndrome = True
show_logicals = True
show_errors = True
####






G = nx.Graph()
length = edge_matrix.shape[0]
pos = {}
node_type = [0]*length
node_color = ['white']*length
sqrt = int(edge_matrix.shape[0]**0.5)


# ancillas
node = 0
while node < length//2:
    for zeile1 in range(math.ceil(sqrt/2),math.ceil(sqrt/2)+sqrt):
        if zeile1 >=sqrt:
            k = int(math.ceil(sqrt/2))
            l = sqrt
        else:
            k =0
            l = int(np.ceil(sqrt/2))
        for spalte1 in range(k,l):
            pos[node] = (zeile1*2%sqrt,spalte1*2%sqrt)

            #syndrome
            if (node-length//4 in syndrome_indices) and show_syndrome == True:
              
                node_color[node] = 'blue'

            node_type[node] = 's'
            G.add_node(node)
            node +=1


# data qubits
node = length//2
while node < length:
    for zeile2 in range(0,sqrt):
        if zeile2 >= int(math.ceil(sqrt/2)):
            k = math.ceil(sqrt/2)
            l = sqrt
        else:
            k =0
            l = math.ceil(sqrt/2)
        for spalte2 in range(k,l):
            pos[node] = (zeile2*2%sqrt,spalte2*2%sqrt)

            #faulty qubit
            if (node-length//2 in error_indices) and show_errors == True:
                node_color[node] = 'red'

            # logicals
            elif (node-length//2 in logicalx_indices or node-length//2 in logicaly_indices) and show_logicals == True:
                node_color[node] = 'green'

            node_type[node] = 'o'
            G.add_node(node)
            node +=1
               

#plotting options
options = {
    #'node_color': 'white',
    'edgecolors': 'black',
    'linewidths': 2,
    'node_size': 200,
}


# generate graph
for i in range(len(edge_matrix)):
    for j in range(i + 1, len(edge_matrix)):
        if edge_matrix[i][j] == 1:
            G.add_edge(i, j)

for i, node in enumerate(G.nodes):
    shape = node_type[i]
    color = node_color[i]
    nx.draw_networkx_nodes(G,pos,nodelist=[node], node_shape=shape,node_color =color, **options)
nx.draw_networkx_edges(G, pos, edge_color= 'black')

#plt.savefig('testtest')



print(SurfCode(hyprep,[9,11,13,14,19],error_indices).draw(show_logicals=True,show_errors=True, show_syndrome=True))



    





        
