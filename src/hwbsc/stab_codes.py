import numpy as np
import networkx as nx
import matplotlib.pyplot as plt
from hwbsc.stabilizers_utils import compute_distance,check_H_commutativity, check_commutativity, find_log_op, check_css, hexa_color_pcm_generator
from hwbsc.codes import steane_code, hamming_code
from ldpc.mod2 import rank, nullspace

class StabCode():

    """
    Asserts if a given check matrix fulfills all condictions to be a stabilizer matrix
    and returns the code parameters

    Parameters:
    -----------
    H : np.ndarray
        check matrix
    name : str
        name of the code
    distance : bool
        function also returns code distance if True

    Returns:
    --------
    name : str
        name of the code
    n : int
        number of qubits of the code
    k : int
        number of logical qubits of the code
    d : int
        distance of the code
    logicals : np.ndarray
        the logical operators of the code

    """

    def __init__(self, stabiliser_check_matrix: np.ndarray, name: str = "Unnamed code", distance: bool = False):
        self.H = stabiliser_check_matrix

        #check that the matrix is a valid GF2 stabiliser check matrix.
        try:
            assert self.H.shape[1] % 2 == 0
        except AssertionError:
            raise AssertionError("Stabiliser check matrix is not a valid GF2 stabiliser check matrix.")

        self.name = name
        print(stabiliser_check_matrix.shape)
        self.n = int(self.H.shape[1]/2) #number of qubits
        self.k = self.n - rank(stabiliser_check_matrix[:18,37:]) #number of logical qubits
        self.d = np.NAN
        if distance == True:
            self.d = compute_distance(stabiliser_check_matrix) #distance of the code

        ###
        # check that stabilisers commutes (ie. that self.H is a valid set of stabilisers)
        try:
            assert check_H_commutativity(self.H) == True
        except AssertionError:
            raise AssertionError("stabilisers do not commute")

        ##compute the logicals
        self.logicals = find_log_op(self.H)
        

        #check the logical commute with all the stabilisers. Throw an error if they don't.
        try:
            for logical in self.logicals:
                for stabilizer in stabiliser_check_matrix:
                    assert check_commutativity(logical, stabilizer) == 0
        except:
            raise AssertionError("logicals do not commute with all the stabilisers")

        #check that the logicals anti-commute with at least one other logical. Throw an error if they don't.
        try:
            for logical in self.logicals:
                for other_logical in self.logicals:
                    if check_commutativity(logical, other_logical) == 1:
                        break
                else:
                    assert False
        except:
            raise AssertionError("logicals do not anti-commute with at least one other logical")


    def code_parameters(self):
        return self.n, self.k, self.d
    
    def __str__(self) -> str:
        """
        The string that is returned when the object is printed
        """
        return f"{self.name}:\n\tParameters: [[n={self.n}, k={self.k}, d={self.d}]]"


class CssCode(StabCode):

    """
    This class inherits from StabCode. Asserts if a given check matrix fulfills all condictions to be a CSS code
    and returns the code parameters

    Parameters:
    -----------
    hx : np.ndarray
        The hx matrix of a CSS code
    hz : np.ndarray
        The hz matrix of a CSS code
    name : str
        name of the CSS code

    Returns:
    --------
    code parameters : str
        The [[n,k,d]] code parameters

    """

    def __init__(self, hx: np.ndarray, hz: np.ndarray, name: str = "Unnamed code", distance: bool = False):

        try:
            assert hx.shape[1] == hz.shape[1]
        except AssertionError:
            raise AssertionError("The hx and hz matrices must have the same number of columns")

        #check that hx and hz are valid stabiliser check matrices. Eg hx@hz.T = 0. Throw an error if they aren't.

        zero = np.zeros((hx.shape[0],hx.shape[0]))
        try:
            assert np.array_equal((hx @ hz.T) % 2 , zero)
        except AssertionError:
            raise AssertionError("hx and hz are not valid stabiliser check matrices")

        stabiliser_check_matrix = np.block([[np.zeros_like(hz),hz],[hx, np.zeros_like(hx)]]).astype(int)

        super().__init__(stabiliser_check_matrix,name,distance)

class SteaneCode(CssCode):

    """
    This class inherits from CssCode and returns the [[7,1,3]] Steane code.
    
    Returns:
    --------
    Steane code : np.ndarray
        The check matrix of the Steane code
    """

    def __init__(self):
        hx = hz = hamming_code(3) #set this to the Hamming Code parity check matrix
        super().__init__(hx, hz, "Steane Code")
        self.d = 3 #we can set the distance manually here

class hypergraph_prod(CssCode):
    """
    Compute the hypergraph product for a classical parity check matrix

    Parameters:
    -----------
    H : np.ndarray
        parity check matrix
    
    Returns:
    --------
    hypergraph_product : np.ndarray
        the hypergraph product of the classical parity check matrix
    
    """
    
    def __init__(self, classical_code: np.ndarray, name: str = "Unnamed code"):

        m = classical_code.shape[0]
        n = classical_code.shape[1]
        
        hx = np.hstack((np.kron(classical_code,np.eye(n)), np.kron(np.eye(m), classical_code.T) ))
        hz = np.hstack((np.kron(np.eye(n),classical_code),np.kron(classical_code.T,np.eye(m))) )

        super().__init__(hx, hz, name)

class Color_code(CssCode):
    def __init__(self, m:int, n:int, name: str = "Color Code"):
        pcm = hexa_color_pcm_generator(m,n)[0]
        mat = pcm
     
        edge_matrix = np.block([[np.zeros((mat.shape[0],mat.shape[0]),dtype=int),mat],
                                [mat.T,np.zeros((mat.shape[1],mat.shape[1]),dtype=int)]])
        colG = nx.from_numpy_array(edge_matrix)
        faces = []

        for i, node in enumerate(edge_matrix):
            if len(list(colG.neighbors(i))) == 6:
                faces.append(i)
        
        
        
        
        

        distance = nx.shortest_path_length(colG, source=faces[4], target=faces[5])
        # initialize all face colors
        nx.set_node_attributes(colG, {n: '0' for n in faces}, name='color')
        colG.nodes[faces[0]]['color'] = 'r'
        colG.nodes[faces[1]]['color'] = 'g'

        colors = {'r','g','b'}


        for face in faces:
            closest_neighbor_faces = []
            if colG.nodes[face]['color'] != '0':
                for face2 in faces:
                    if face != face2 and nx.shortest_path_length(colG, source=face, target=face2) == 2 and colG.nodes[face2]['color'] != '0':
                        for face3 in faces:
                            if nx.shortest_path_length(colG, source=face, target=face3) == 2  and nx.shortest_path_length(colG, source=face2, target=face3) == 2:
                                if face != face3 and face2 != face3:
                                    newcolor = next(iter(colors.difference({colG.nodes[face]['color'], colG.nodes[face2]['color']})))
                                    colG.nodes[face3]['color'] = newcolor


        # print(faces)
        facegraph = colG.subgraph(faces).copy()

        # add edges to facegraph
        edges = []
        for face in faces:
            for neighbor in faces:
                if face != neighbor and nx.shortest_path_length(colG, source=face, target=neighbor) == 2:
                    facegraph.add_edge(face, neighbor)
                    edgecolor = next(iter(colors.difference({facegraph.nodes[face]['color'], facegraph.nodes[neighbor]['color']})))
                    facegraph.edges[(face,neighbor)]['color'] = edgecolor

    
        node_colors = [data['color'] for _, data in facegraph.nodes(data=True)] 
        edge_colors = [data['color'] for _,_, data in facegraph.edges(data=True)] 

        self.facegraph = facegraph, node_colors, edge_colors
        # nx.draw(facegraph, node_color=node_colors, edge_color=edge_colors, with_labels = True)
        
        # plt.show()
        # print(facegraph.nodes(data=True))


        for face in colG:
            if face in faces:           
                rc = 0
                gc = 0
                bc = 0
                c=0
                for neighbor in colG:
                    if nx.shortest_path_length(colG, source=face, target=neighbor) < 3:
                        if neighbor in faces and face != neighbor:
                            c += 1
                            # print(face,neighbor)
                            if colG.nodes[neighbor]['color'] == 'r':
                                rc += 1
                            if colG.nodes[neighbor]['color'] == 'g':
                                gc += 1
                            if colG.nodes[neighbor]['color'] == 'b':
                                bc += 1
                if c!= 6:
                    print("eeerror")

        
        # check if legit
        for node in facegraph.nodes:
            colorsdic = {'r':0,'g':0,'b':0}
            for neighbor in facegraph.neighbors(node):
                #print(neighbor,"\n")
                othercolors = []
                if facegraph.nodes[neighbor]['color'] == 'g':
                    colorsdic['g'] += 1
                elif facegraph.nodes[neighbor]['color'] == 'b':
                    colorsdic['b'] += 1
                elif facegraph.nodes[neighbor]['color'] == 'r':
                    colorsdic['r'] += 1
            del colorsdic[facegraph.nodes[node]['color']]
            keys = list(colorsdic.keys())
            if not colorsdic[keys[0]] == 3 and colorsdic[keys[1]] == 3:
                raise ValueError("Failed to color graph")

        super().__init__(pcm, pcm, name)



