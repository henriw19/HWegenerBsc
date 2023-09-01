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

fig, ax = plt.subplots()

data1 = np.loadtxt(f"thresholddata/coloruniondistance12data.txt", delimiter=',', skiprows=1)
# data2 = np.loadtxt(f"coloruniondistance18_6bestdata.txt", delimiter=',', skiprows=1)
# data3 = np.loadtxt(f"coloruniondistance18_12bestdata.txt", delimiter=',', skiprows=1)
# data4 = np.loadtxt(f"coloruniondistance24_6bestdata.txt", delimiter=',', skiprows=1)
# data5 = np.loadtxt(f"coloruniondistance24_18bestdata.txt", delimiter=',', skiprows=1)
# data6 = np.loadtxt(f"coloruniondistance30_6bestdata.txt", delimiter=',', skiprows=1)
# data7 = np.loadtxt(f"coloruniondistance30_24bestdata.txt", delimiter=',', skiprows=1)
# data8 = np.loadtxt(f"coloruniondistance24_12data.txt", delimiter=',', skiprows=1)
# data9 = np.loadtxt(f"coloruniondistance36data.txt", delimiter=',', skiprows=1)



p1 = data1[5:, 0]
# p2 = data2[7:, 0]
# p3 = data3[7:, 0]
# p4 = data4[5:, 0]
# p5 = data5[8:, 0]
# p6 = data6[5:, 0]
# p7 = data7[8:, 0]
# p8 = data8[:, 0]
# p9 = data9[:, 0]


log1 = data1[5:, 1]
# log2 = data2[7:, 1]
# log3 = data3[7:, 1]
# log4 = data4[5:, 1]
# log5 = data5[8:, 1]
# log6 = data6[5:, 1]
# log7 = data7[8:, 1]
# log8 = data8[:, 1]
# log9 = data9[:, 1]



ax.plot(p1,log1, label= '[[144,4,8]]', marker='o')
# ax.plot(p2,log2, label= '[[216,4,12]]',marker='o')
# ax.plot(p4,log4, label = '[[288,4,8]]',marker='o')
# ax.plot(p6,log6, label = '[[360,4,8]]',marker='o')
# ax.plot(p3,log3, label = '[[432,4,16]]',marker='o')
# ax.plot(p5,log5, label = '[[864,4,24]]',marker='o')
# ax.plot(p7,log7, label = '[[1440,4,32]]',marker='o')
# ax.plot(p8,log8, label = '[[xx]]',marker='o')
# ax.plot(p9,log9, label = '[[xxxx]]',marker='o')





# ax.plot(np.sort(np.hstack((p1,p2))), np.sort(np.hstack((log1,log2))), color = 'plum')
# ax.plot(np.sort(np.hstack((p3,p4))), np.sort(np.hstack((log3,log4))), color = 'mediumseagreen')


# ax.set_xscale('log')
# ax.set_yscale('log')
ax.set_xlabel('x')
ax.set_ylabel('y')
ax.legend()
plt.axvline(x=0.077, color='0', linestyle='--')


plt.xscale('log') 
plt.yscale('log')
# plt.xticks(np.logspace(-1,-0.8,6))
# plt.yticks(np.logspace(-1,-0,4,20))
plt.xlabel("phys error rate")
plt.ylabel("logical error rate")
plt.legend(loc = 0)
plt.title('Hexagonal Toric Color Code Threshold')
# plt.show()
plt.savefig('img/hexacolor.jpg')

# syndrome = [1,0,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0]
# pred_error = union_decoder(H, syndrome)
# print(pred_error)
# SurfCode(H,syndrome, predicted_error = pred_error).draw(H, show_syndrome=True, show_errors = True)









# rescale the logical error rate (word err rate?)
# more iterations on lower phys err rate
# same ratio on toric code