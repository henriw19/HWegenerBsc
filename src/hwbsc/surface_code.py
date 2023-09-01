import numpy as np
import networkx as nx
import math
import matplotlib.pyplot as plt
from hwbsc.stab_codes import *
from hwbsc.codes import *

class SurfCode():

    """
    
    """

    def __init__(self, stabiliser_check_matrix:np.ndarray, syndrome = [], predicted_error = []):
        # ....
        self.stabiliser_check_matrix = stabiliser_check_matrix
        self.a = stabiliser_check_matrix.shape[0]
        self.b = stabiliser_check_matrix.shape[1]
        logicals = StabCode(stabiliser_check_matrix).logicals
        logicalx = logicals[0][:self.b//2]
        logicaly = logicals[1][self.b//2:]
        self.logicalx = logicalx
        self.logicaly = logicaly
        self.logicals = np.block([[logicalx],[logicaly]])
        self.syndrome = syndrome
        self.predicted_error = predicted_error
        print(syndrome)
        print(predicted_error)
        
    

    def draw(self, show_logicals = False, show_syndrome = False, show_errors = False):
        # a,b = self.stabiliser_check_matrix.shape
        a = self.a
        b= self.b
        mat = np.block([[self.stabiliser_check_matrix[a//2:,:b//2]],[self.stabiliser_check_matrix[:a//2,b//2:]]])
        edge_matrix = np.block([[np.zeros((mat.shape[0],mat.shape[0]),dtype=int),mat],[mat.T,np.zeros((mat.shape[1],mat.shape[1]),dtype=int)]])
        logicalx_indices = np.where(self.logicalx == 1)[0]
        logicaly_indices = np.where(self.logicaly == 1)[0]
        

        #### test error
        # error = np.zeros(12)
        # error[3] = 1
        # error[9] = 1
        # error[10] = 1
        error_indices = self.predicted_error
        syndrome_indices = self.syndrome
        print(error_indices)
        ####


        G = nx.Graph()
        
        length = edge_matrix.shape[0]
        pos = {}
        node_type = [0]*length
        node_colors = ['white']*length
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
                        node_colors[node] = 'blue'
                        

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
                        print("err")
                        node_colors[node] = 'red'

                    # logicals
                    elif (node-length//2 in logicalx_indices or node-length//2 in logicaly_indices) and show_logicals == True:
                        node_colors[node] = 'green'

                    node_type[node] = 'o'
                    G.add_node(node)
                    node +=1
                    
        #plotting options
        options = {
            'edgecolors': 'black',
            'linewidths': 2,
            'node_size': 200,
        }

        # generate graph
        for i in range(len(edge_matrix)):
            for j in range(i + 1, len(edge_matrix)):
                if edge_matrix[i][j] == 1:
                    G.add_edge(i, j)
        node_colors[7] = 'blue'
        node_colors[8] = 'blue'
        node_colors[13] = 'blue'
        node_colors[17] = 'blue'
        node_colors[19] = 'blue'



        for i, node in enumerate(G.nodes):
            shape = node_type[i]
            color = node_colors[i]
            
            nx.draw_networkx_nodes(G,pos,nodelist=[node], node_shape=shape,node_color = color, **options)
            
            # G.add_node(node, node_shape=shape,node_color =color, **options)
        # print(nx.get_node_attributes(G, 'color'))
        
        nx.draw_networkx_edges(G, pos, edge_color= 'black')

        plt.savefig('surfaceplot',dpi=1000)
        print("yay")
        plt.show()
        return 0
        

    def __str__(self) -> str:
        """
        The string that is returned when the object is printed
        """
        return '-'

