import numpy as np
import matplotlib.pyplot as plt
from hwbsc.stab_codes import *

# n = 3
# rep = np.zeros((n-1,n),dtype=int)
# for i in range(n-1):
#     for j in range(n):
#         rep[i,i] = 1
#         rep[i,i+1] = 1

# H = hypergraph_prod(rep).H

# hz = np.vstack((np.kron(np.eye(H.shape[0],dtype = int),H),np.kron(H.T,np.eye(H.shape[1],dtype = int))))
# hx = np.vstack((np.kron(H,np.eye(H.shape[0],dtype = int)),np.kron(np.eye(H.shape[1],dtype = int),H.T)))
# hi = hx
# matrix = np.block([[np.zeros((hi.shape[0],hi.shape[0]),dtype=int),hi],[hi.T,np.zeros((hi.shape[1],hi.shape[1]),dtype=int)]])

# print(matrix)
#SurfCode().graph(hyprep)


# G = nx.grid_2d_graph(4, 4)
# nx.draw(G)
# plt.show()



# mat = np.block([[hyprep[6:,:13]],[hyprep[:6,13:]]])
# factor_graph_matrix = np.block([[np.zeros((mat.shape[0],mat.shape[0]),dtype=int),mat],[mat.T,np.zeros((mat.shape[1],mat.shape[1]),dtype=int)]])


# G = nx.Graph(factor_graph_matrix)

# # Define node shapes and colors
# node_shapes = {0: 'o', 1: 's'}
# node_colors = {0: 'red', 1: 'blue'}

# # Set node attributes for shape and color
# #node_attributes = {'shape': [node_shapes[node] for node in G.nodes()],
# #                   'color': [node_colors[node] for node in G.nodes()]}

# # Plot the graph matrix
# pos = nx.spring_layout(G)
# nx.draw_networkx(G, pos, node_size=500, **node_attributes)

# # Customize the plot
# plt.axis('off')
# plt.title("Graph Matrix")

# # Show the plot
# #plt.show()


# get data:

colors = ['cornflowerblue','khaki','plum','mediumseagreen','lightsteelblue','moccasin', 'yellowgreen','salmon']
skiplist = [4,7,10,13,5,1,1,1]

fig, ax = plt.subplots()
for distance in range(3,10,2):
    
    data = np.loadtxt(f"uniondistance{distance}data.txt", delimiter=',', skiprows=1)
    p = data[:, 0]
    log = data[:, 1]
    ax.plot(p, log, color = colors[int((distance-3)/2)])

# ax.set_xscale('log')
# ax.set_yscale('log')
ax.set_xlabel('x')
ax.set_ylabel('y')
ax.legend()
plt.axvline(x=0.16, color='0', linestyle='--')

plt.xscale('log') 
plt.yscale('log')
# plt.xticks(np.logspace(-1,-0.8,6))
# plt.yticks(np.logspace(-1,-0,4,20))
plt.xlabel("phys err rate")
plt.ylabel("logical err rate")
plt.legend(loc = 0)
# plt.show()

# syndrome = [1,0,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0]
# pred_error = union_decoder(H, syndrome)
# print(pred_error)
# SurfCode(H,syndrome, predicted_error = pred_error).draw(H, show_syndrome=True, show_errors = True)

