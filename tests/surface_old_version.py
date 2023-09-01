import networkx as nx
import numpy as np
import matplotlib.pyplot as plt
from hwbsc.stab_codes import *
from hwbsc.decoders import *

#### code
n = 3

rep = np.zeros((n-1,n),dtype=int)
for i in range(n-1):
    for j in range(n):
        rep[i,i] = 1
        rep[i,i+1] = 1
rep = np.array([[1,1,0],[0,1,1]])

hyprep = hypergraph_prod(rep).H
a = hyprep.shape[0]
b = hyprep.shape[1]
mat = np.block([[hyprep[a//2:,:b//2]],[hyprep[:a//2,b//2:]]])
edge_matrix = np.block([[np.zeros((mat.shape[0],mat.shape[0]),dtype=int),mat],[mat.T,np.zeros((mat.shape[1],mat.shape[1]),dtype=int)]])

####


#### error
error = np.zeros(12)
error[3] = 1
error[9] = 1
error[10] = 1

error_indices = np.where(error == 1)[0]

syndrome_indices = np.array([9,1,3])
####

#### logicals
logicals = StabCode(hyprep).logicals
logicalx = logicals[0][:b//2]
logicaly = logicals[1][b//2:]
logicalx_indices = np.where(logicalx == 1)[0]
logicaly_indices = np.where(logicaly == 1)[0]
####

####
show_syndrome = False
show_logicals = False
show_errors = False
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
    for zeile1 in range(1,sqrt,2):
        for spalte1 in range(0,sqrt,2):
            pos[node] = (zeile1,spalte1)
            print(zeile1,spalte1)

            #syndrome
            if (node-length//4 in syndrome_indices) and show_syndrome == True:
                node_color[node] = 'blue'

            node_type[node] = 's'
            G.add_node(node)
            node +=1

    for zeile2 in range(0,sqrt,2):
        for spalte2 in range(1,sqrt,2):
            pos[node] = (zeile2,spalte2)
            print(zeile2,spalte2)

            #syndrome
            if (node-length//4 in syndrome_indices) and show_syndrome == True:
                node_color[node] = 'blue'

            node_type[node] = 's'
            G.add_node(node)
            node +=1









# data qubits
node = length//2
while node < length:
    for zeile3 in range(0,sqrt,2):
        for spalte3 in range(0,sqrt,2):
            pos[node] = (zeile3,spalte3)

            #faulty qubit
            if (node-length//2 in error_indices) and show_errors == True:
                node_color[node] = 'red'

            elif (node-length//2 in logicalx_indices or node-length//2 in logicaly_indices) and show_logicals == True:
                print(node-length//2)
                node_color[node] = 'green'

            node_type[node] = 'o'
            G.add_node(node)
            node +=1
                    
    for zeile4 in range(1,sqrt,2):
        for spalte4 in range(1,sqrt,2):
            pos[node] = (zeile4,spalte4)

            #faulty qubit
            if (node-length//2 in error_indices) and show_errors == True:
                node_color[node] = 'red'
            elif (node-length//2 in logicalx_indices or node-length//2 in logicaly_indices) and show_logicals == True:
                

                node_color[node] = 'green'

            node_type[node] = 'o'
            G.add_node(node)
            node +=1





options = {
    #'node_color': 'white',
    'edgecolors': 'black',
    'linewidths': 2,
    'node_size': 200,
}


for i in range(len(edge_matrix)):
    for j in range(i + 1, len(edge_matrix)):
        if edge_matrix[i][j] == 1:
            G.add_edge(i, j)


for i, node in enumerate(G.nodes):
    shape = node_type[i]
    color = node_color[i]
    nx.draw_networkx_nodes(G,pos,nodelist=[node], node_shape=shape,node_color =color, **options)

nx.draw_networkx_edges(G, pos, edge_color= 'black')


plt.savefig('test2')

print(G)



for i1 in range(3,8):
    print("hier")
    if i1 >=5:
        k = 3
        l = 5
    else:
        k =0
        l = 3
    for i2 in range(k,l):
        
        print((i1 *2 % 5, i2*2%5))
    




