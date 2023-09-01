import networkx as nx
from pymatching import Matching
from numpy import zeros, uint8, linspace, array
from random import random
from random import sample
from typing import List, FrozenSet, Set
from os import makedirs
from os.path import exists
from time import time
from matplotlib.pyplot import figure, savefig, title, show, xlabel, ylabel, legend, errorbar, close, plot
"""
The main function of this file will generate a set of images
of the pertaining color code graph, its dual, and its respective 
2- colored subgraphs and print an error prediction on the subgraphs.

A folder "img/hexcolor/" will be created to save image files if it does not exist.
"""

def colorize_graph_black(G: nx.Graph) -> bool:
    """
    Args:
        G(nx.Graph): some graph
    Returns:
        bool: wether it was successful in changing the graph object 
              to all black edges.
    """
    for u, v, attr in G.edges(data=True):
        G[u][v]['color'] = 'black'
    return True

def tor_hex48_color_encode(G: nx.Graph,m: int=6,n: int=4) -> bool:
    """
    Args:
        G(nx.Graph): graph we want to encode with three colored faces
        n,m: how many by how many hexagon, default to 6 and 4 like in delfosse
    Returns:
        bool: Success of graph object modification procedure
    """
    rgb_list = ['r', 'g', 'b']
    # initialize all edge colors to black
    for u, v, attr in G.edges(data=True):
        G[u][v]['color'] = 'black'
    
    # colorizing algorithm
    
    # horizontal edges
    for i in range(int(n/2)):
            for j in range(m):
                first_coordinate = (2*i,2*j)
                second_coordinate = (2*i+1,2*j)
                G[first_coordinate][second_coordinate]['color'] = rgb_list[j%3]
    for i in range(int(n/2)):
        for j in range(m):
            first_coordinate = (2*i+1,2*j+1)
            second_coordinate = ((2*i+2)%n,2*j+1)
            G[first_coordinate][second_coordinate]['color'] = rgb_list[(j-1)%3]
    # left ladder edges
    for i in range(int(n/2)):
        for j in range(2*m):
            first_coordinate = (2*i,j)
            second_coordinate = (2*i,(j+1)%(2*m))
            G[first_coordinate][second_coordinate]['color'] = rgb_list[(1-j)%3]
    # right ladder edges
    for i in range(int(n/2)):
        for j in range(2*m):
            first_coordinate = (2*i+1,j)
            second_coordinate = (2*i+1,(j+1)%(2*m))
            G[first_coordinate][second_coordinate]['color'] = rgb_list[(1-j)%3]
    return True

def make_a_base_graph(m: int=6,n: int=4) -> nx.Graph:
    """
    Args:
        m(int), n(int): desired dimension of faces on graph (m by n)
    Returns:
        G(nx.Graph): A basic color code graph of dimensions m by n
                     with colored edges encircling opposite colored faces.
    """
    G = nx.hexagonal_lattice_graph(m, n, periodic=True)
    colorize_graph_black(G)
    if not tor_hex48_color_encode(G,m,n):
        raise ValueError("Graph ist dumm")
    for node in G.nodes:
        G.nodes[node]['color'] = 'yellow'
        G.nodes[node]['fault_ids'] = 0
    return G

def draw_graph_with_colored_edges_and_nodes(G: nx.Graph, file: str=None, name: str=None) -> bool:
    """
    Draws a graph who's nodes and edges have colors.
    Options:
        filename (str): save file to specified name (will plt.show() otherwise)
        name (str): will create figure with specified name
    """
    pos = nx.get_node_attributes(G, 'pos')
    node_colors = [data['color'] for _, data in G.nodes(data=True)] 
    edge_colors = [G[u][v]['color'] for u, v in G.edges()]

    figure()
    if name:
        title(name)
    if pos: 
        nx.draw(G, pos, with_labels=True, node_color=node_colors, edge_color=edge_colors)
    elif not pos:
        nx.draw(G, with_labels=True, node_color=node_colors, edge_color=edge_colors)
    if nx.get_edge_attributes(G,"fault_ids"):
        nx.draw_networkx_edge_labels(G, pos, edge_labels=nx.get_edge_attributes(G, "fault_ids"))
    if file:
        savefig(file)
    else:
        show()
    close()
    return True

def flag_color_graph(graph: nx.Graph, per: float=0.1) -> Set[any]:
    """
    Args:
        graph(nx.Graph): graph to be altered with errors on nodes
        per(float): probability on error occuring on each node
    Returns:
        set(node): Actually occured errors
    """
    error_nodes = set()
    for node in graph.nodes:
        if random()< per:
            graph.nodes[node]['fault_ids'] = 1
            graph.nodes[node]['color'] = 'y'
            error_nodes.add(node)
    return error_nodes

def find_6_loops(graph: nx.Graph) -> List[FrozenSet[any]]:
    """
    Args:
        nx.Graph: input graph
    Returns:
        set[frozenset]: Topology of nodes comprising faces on input graph """
    cycles = set()
    for node in graph.nodes:
        for node1 in graph.neighbors(node):
            for node2 in graph.neighbors(node1):
                for node3 in graph.neighbors(node2):
                    for node4 in graph.neighbors(node3):
                        for node5 in graph.neighbors(node4):
                            for node6 in graph.neighbors(node5):
                                if node6 == node:
                                    cycles.add(frozenset([node, node1, node2, node3, node4, node5]))
    faces = [cycle for cycle in cycles if len(cycle) == 6]
    return faces

def find_face_color(graph: nx.Graph, face: FrozenSet) -> str:
    """
    Args:
        graph: graph on which face lies
        face(set[nodes]): face to analyze
    Returns:
        color(str): color of face
    """
    rgb = set(['r','g','b'])
    boundary_colors = set()
    for node in face:
        for node2 in face:
            if node2 in graph.neighbors(node):
                boundary_colors.add(graph[node][node2]['color'])
    face_color = rgb - boundary_colors
    face_color = face_color.pop()
    return face_color
        
def dual_of_three_colored_graph(graph: nx.Graph):# -> nx.Graph:
    """
    Args:
        graph(nx.Graph): graph of which we want the dual 
    Returns:
        dual_graph(nx.Graph): the dual of that graph
    """
    dual_graph = nx.Graph()
    faces = find_6_loops(graph)
    # init nodes
    for i, face in enumerate(faces):
        dual_graph.add_node(i, color = 'black')
        color_of_face = find_face_color(graph, face)
        dual_graph.nodes[i]['color'] = color_of_face
        # This is the part for error -> syndrome inheritance to the dual graph
        dual_graph.nodes[i]['fault_ids'] = 0
        for node in face:
            if graph.nodes[node]['fault_ids'] == 1:
                dual_graph.nodes[i]['fault_ids'] = (dual_graph.nodes[i]['fault_ids']+1)%2
    # connect nodes
    for i, face in enumerate(faces):
        otherfaces = faces[:i]+faces[((i+1)%(len(faces)+1)):]
        for j, face2 in enumerate(otherfaces):
            lap_nodes = set(face & face2)
            if lap_nodes:
                # A three-colorable graph will only ever have two nodes between two faces
                node1 = lap_nodes.pop()
                node2 = lap_nodes.pop()
                connecting_color = graph[node1][node2]['color']
                # we are iterating not over the list of faces, but over the list of otherfaces
                # what is j?
                second_face_pos = [k for k in range(len(faces)) if faces[k] == otherfaces[j]].pop()
                #
                dual_graph.add_edge(i,second_face_pos, color = connecting_color)
                
    
    return dual_graph, faces

def subtile(Graph: nx.Graph, color: str) -> nx.Graph:
    """
    Args:
        Graph(nx.Graph): graph we want to subtile
        color(str): color in format "r","g", "b" of which all edges
                    in the subtiling will be comprised
    Returns: 
        G(nx.Graph): subtiled graph (does not edit original object)
    """
    G = Graph.copy()
    for edge in G.edges:
        u, v = edge[0], edge[1]
        if G.edges[u,v]['color'] != color:
            G.remove_edge(u,v)
    G.remove_nodes_from(list(nx.isolates(G)))
    return G

def decode_subtile(graph: nx.Graph) -> List[any]:
    """
    Args:
        graph(nx.Graph): graph with "fault_ids" property on some nodes
    Returns:
        prediction(List[edges]): predicted error edges on graph
    """
    for i, node in enumerate(graph.nodes):
        graph.nodes[node]['og_name'] = node
        graph = nx.relabel_nodes(graph,{node: i})
    matching = Matching(graph)
    # generate syndrome on renamed_copy
    syndrome = zeros(len(graph.nodes), dtype=uint8)
    for node in graph.nodes:
        if graph.nodes[node]['fault_ids'] == 1:
            syndrome[node] = 1
    # predict edges on the renamed_copy
    prediction = matching.decode_to_edges_array(syndrome)
    # rename nodes to be actually usable
    for edge in prediction:
        for i in range(len(edge)):
            edge[i] = graph.nodes[edge[i]]['og_name']
    # revert the graph back to normal
    for node in graph.nodes:
        graph = nx.relabel_nodes(graph, {node: graph.nodes[node]['og_name']})
    
    return prediction

def make_a_shower(graph: nx.Graph) -> nx.Graph:
    """
    Args:
        graph(nx.Graph): graph we want a yellow syndrome flagged copy of
    Returns:
        shower(nx.Graph): graph with yellow marked syndrome nodes
    """
    shower = graph.copy()
    for node in shower.nodes:
        if shower.nodes[node]['fault_ids'] == 1:
            shower.nodes[node]['color'] = 'y'
    return shower

def find_hyper_edges(dual_graph: nx.Graph, edges_array_r: List[any],
                     edges_array_g: List[any], edges_array_b: List[any]) -> List[any]:
    """ 
    Takes: a dualgraph and its subgraph matching edges arrays
    Returns: list of cycles on the dual graph
    """
    # generate the surrounding edges of cycles
    set_of_all_edges_bounding_hyperedge = set()
    for color in [edges_array_r, edges_array_g, edges_array_b]:
        for edge in color:
            addable_edge = tuple(sorted(edge))
            set_of_all_edges_bounding_hyperedge.add(addable_edge)
    # make a graph of only error cycles
    error_bound_graph = dual_graph.copy()
    bad_edges = []
    # this is necessary because we can't modify edges during iteration
    for edge in error_bound_graph.edges:
        if edge not in set_of_all_edges_bounding_hyperedge:
            bad_edges.append(edge)
    error_bound_graph.remove_edges_from(bad_edges)
    isolates = list(nx.isolates(error_bound_graph))
    error_bound_graph.remove_nodes_from(isolates)
    cycles = nx.cycle_basis(error_bound_graph)
    return cycles     

def flag_c_graph_specific(graph: nx.Graph, nodes: List[any]) -> bool:
    """ 
    Flags down specific nodes on a graph from a list of nodes
    """
    for node in nodes:
        graph.nodes[node]['fault_ids'] = 1
        graph.nodes[node]['color'] = 'y'
    return True

def lift(dual_edge_cycles: List[any], faces: List[FrozenSet]) -> Set[any]:
    """
    Takes: List of dual graph cycles, facenodes to face map
    Returns: List of enclosed og nodes
    """
    ## Find corresponding nodes in original graph
    enc_nodes = set()
    for dual_edge_cycle in dual_edge_cycles:
        face_on_dec = dual_edge_cycle.pop()
        bounded_nodes = faces[face_on_dec]
        for face in dual_edge_cycle:
            bounded_nodes = bounded_nodes & faces[face]
        bounded_nodes = frozenset(bounded_nodes)
        enc_nodes.add(bounded_nodes)
    # clear up empty frozenset and pop the items to a set
    if frozenset() in enc_nodes:
        enc_nodes.remove(frozenset())
    res = set()
    for enc_node in enc_nodes:
        res.add(next(iter(set(enc_node))))
    return res

def total_decoder(graph: nx.Graph, per: float) -> bool:
    """ 
    Takes: a color code graph and physical error rate
    Returns: Success of correction operation
    """
    actual_errors = flag_color_graph(graph, per)
    #### dualizing and subtiling
    dual, faces = dual_of_three_colored_graph(graph)
    subr, subg, subb = subtile(dual, 'r'), subtile(dual, 'g'), subtile(dual, 'b')
    #### decoding part
    pred_r, pred_g, pred_b = decode_subtile(subr), decode_subtile(subg), decode_subtile(subb)
    hyper_edge_cycles = find_hyper_edges(dual, pred_r, pred_g, pred_b)
    ## get back to og nodes from dual nodes/ faces
    og_enc_nodes_by_dual_cycles = lift(hyper_edge_cycles, faces)
    return og_enc_nodes_by_dual_cycles == actual_errors

def cc_ler_calc(graph: nx.Graph, per: float, nr: int) -> float:
    numErrors = 0
    for _ in range(nr):
        if not total_decoder(graph, per):
            numErrors += 1
    return numErrors/nr

def cc_threshold_plotter(dists: List[any], pers: List[float], nr:int, file=None) -> bool:
    log_errors_all_dist = []
    for d in dists:
        print("Simulating d = {}".format(d))
        origG = make_a_base_graph(d[0],d[1])
        lers = []
        for per in pers:
            print(f"per={per}")
            graph = origG.copy()
            lers.append(cc_ler_calc(graph, per, nr))
        log_errors_all_dist.append(array(lers))
    figure()
    for dist, logical_errors in zip(dists, log_errors_all_dist):
        std_err = (logical_errors*(1-logical_errors)/nr)**0.5
        errorbar(pers, logical_errors, yerr=std_err, label="L={}".format(dist))
    plot(pers, pers, label = 'Threshold')
    xlabel("Physical error rate")
    ylabel("Logical error rate")
    legend(loc=0)
    if file:
        if not exists("img/hexcolor"):
            makedirs("img/hexcolor")
        savefig("img/hexcolor/"+file)
    else:
        show()  
    close()   
    return True

def main() -> bool:

   
    #### just making sure image filesaves work
    if not exists("img/hexcolor"):
        makedirs("img/hexcolor")
    #### initialize color code graph with errors
    origG = make_a_base_graph()
    actual_errors = flag_color_graph(origG, 0.05)
    ## This is for manually setting faults
    # actual_errors = [(0,0),(1,2)]
    # flag_c_graph_specific(origG, actual_errors)
    #### dualizing and subtiling
    dual, faces = dual_of_three_colored_graph(origG)
    subr, subg, subb = subtile(dual, 'r'), subtile(dual, 'g'), subtile(dual, 'b')
    #### flag syndromes yellow for better visualizing
    dual_syn, subr_syn, subg_syn, subb_syn = make_a_shower(dual), make_a_shower(subr), make_a_shower(subg), make_a_shower(subb)
    #### decoding part
    start = time()
    pred_r, pred_g, pred_b = decode_subtile(subr), decode_subtile(subg), decode_subtile(subb)
    hyper_edge_cycles = find_hyper_edges(dual, pred_r, pred_g, pred_b)
    ## get back to og nodes from dual nodes/ faces
    print("hyp_edge_cycles are: ", hyper_edge_cycles)
    og_enc_nodes_by_dual_cycles = lift(hyper_edge_cycles, faces)
    end = time()
    #### visualizing part
    print(f"This decoding and lifting took {end-start} seconds.")
    draw_graph_with_colored_edges_and_nodes(origG, "img/hexcolor/original.png")
    draw_graph_with_colored_edges_and_nodes(dual_syn, "img/hexcolor/dual.png")
    for i, graph in enumerate([subr_syn, subg_syn, subb_syn]):      
        draw_graph_with_colored_edges_and_nodes(graph, f"img/hexcolor/{i}.png")
    if og_enc_nodes_by_dual_cycles:
        for i in range(len(og_enc_nodes_by_dual_cycles)):
            print(f"The {i}th error node on the graph is {og_enc_nodes_by_dual_cycles.pop()}")
    print("The actual errors were: ", actual_errors)
    return True

if __name__ == "__main__":
    main()