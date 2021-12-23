from collections import defaultdict
import networkx as nx
from networkx import isomorphism
from networkx import readwrite
from networkx.algorithms.connectivity import build_auxiliary_node_connectivity
from networkx.algorithms.flow import build_residual_network
import itertools
import matplotlib.pyplot as plt
from networkx.algorithms.operators.product import power
from networkx.linalg.algebraicconnectivity import _rcm_estimate
import numpy as np
# from primefac import *
from networkx.algorithms.flow import shortest_augmenting_path
import math
import random
from ortools.sat.python import cp_model

def one(n, idx, v):
    arr = [0 for _ in range(n)]
    arr[idx] = v
    return tuple(arr)

def clamp(x, m, M):
    return max(min(x, M), m)

def inrange(X, lims):
    return all([0 <= X[i] < lims[i] for i in range(len(X))])

def rect_graph(x, y, z):
    lims = [x, y, z]
    node_mapping = dict(zip(range(x*y*z), itertools.product(range(x), range(y), range(z))))
    G = nx.Graph()
    G.add_nodes_from(itertools.product(range(x), range(y), range(z)))
    directions = [one(3, idx, s) for idx in range(3) for s in [1, -1]]
    for i in range(x*y*z):
        neighbours = [tuple(map(lambda m, n: m+n, node_mapping[i], direction)) for direction in directions]
        neighbours = list(map(lambda x: (node_mapping[i], x), filter(lambda L: inrange(L, lims), neighbours)))
        G.add_edges_from(neighbours)
    return G, node_mapping

def manhattan(p1, p2):
    return sum(map(lambda x, y: abs(x-y), p1, p2))

def compactible(H1 : nx.DiGraph, dim, upper):
    # is H a topological minor of a 3D grid graph with dimensions dim, with extra constraints
    # manhattan distance between neighbouring nodes of G should be leq than upper (rough bound for path length)
    # input/output nodes must be on the surface of the grid 
    # SAT solving method partially based on "topological_minor" function from SAGE
    model = cp_model.CpModel()
    G, _ = rect_graph(*dim)
    H = H1.copy()
    H.add_node('power') #add power providing node
    M = dict([(hvert, dict([(gvert, model.NewBoolVar(f'{gvert}{hvert}')) for gvert in G])) for hvert in H]) #simple degree pruning
    # print(M['n0'])
    boundary_nodes = [node for node, data in H.nodes(data=True) if 'description' in data and ('input' in data['description'] or 'output' in data['description'])] + ['power']
    for hvert in H:
        for gvert in G:
            if G.degree[gvert] < H.degree[hvert]:
                model.Add(M[hvert][gvert] == 0)

    for boundary_node in boundary_nodes:
        for gvert in G:
            if G.degree[gvert] == 6: #boundary nodes must have open face (not degree 6)
                model.Add(M[boundary_node][gvert] == 0)
    for hvert in H:
        model.Add(sum(M[hvert].values()) == 1) #each vert in H must be mapped to a vert in G
    for gvert in G:
        model.Add(sum([M[hvert][gvert] for hvert in H]) <= 1) #a vert in G must only be mapped to at most one H
    # print(status == cp_model.OPTIMAL)
    for gvert in G:
        for gvert2 in G:
            if manhattan(gvert, gvert2) > upper: 
                for edge in H.edges():
                    model.Add(M[edge[0]][gvert] + M[edge[1]][gvert2] < 2) #no neighbours with manhattan distance > upper
    is_mapped = dict([(gvert, model.NewBoolVar(f'{gvert}mapped')) for gvert in G])
    for gvert in G:
        model.Add(is_mapped[gvert] == sum([M[hvert][gvert] for hvert in H]))
    #treat each edge of H as a flow problem
    #an edge (h1, h2) is defined as a "commodity" in the problem
    #g(h1) is the mapped version of h1 onto g. it will be a source node
    #g(h2) is the mapped version of h2 onto g. it will be a sink node
    #source nodes have flow_out - flow_in = 1
    #sink nodes have flow_out - flow_in = -1
    #other nodes have flow_out - flow_in = 0

    flow = dict([(hedge, dict([(gedge[::order], model.NewBoolVar(f'{gedge}{hedge}')) for gedge in G.edges() for order in [1, -1]])) for hedge in H.edges()])
    flow_in = lambda hedge, gvert: sum([flow[hedge][(neigh, gvert)] for neigh in G.neighbors(gvert)]) #assuming G is undirected
    flow_out = lambda hedge, gvert: sum([flow[hedge][(gvert, neigh)] for neigh in G.neighbors(gvert)])
    
    #add power flow (source is some node named "power" with <<number of other nodes in H>> flow, each other node in H is sink with flow 1)
    #special commodity called "power", flow network where nodes have more capacity than just 1
    #still cant intersect with other flows
    power_flow = dict([(gedge[::order], model.NewIntVar(0, H.number_of_nodes()-1, f'{gedge}power')) for gedge in G.edges() for order in [1, -1]])
    power_flow_in = dict([(gvert, model.NewIntVar(0, H.number_of_nodes()-1, f'{gvert}power_in')) for gvert in G])
    power_flow_out = dict([(gvert, model.NewIntVar(0, H.number_of_nodes()-1, f'{gvert}power_out')) for gvert in G])
    is_mapped_not_power = dict([(gvert, model.NewBoolVar(f'{gvert}notpower')) for gvert in G])
    
    for gvert in G:
        model.Add(power_flow_in[gvert] == sum([power_flow[(neigh, gvert)] for neigh in G.neighbors(gvert)]))
        model.Add(power_flow_out[gvert] == sum([power_flow[(gvert, neigh)] for neigh in G.neighbors(gvert)]))
   
        model.Add(power_flow_out[gvert] - power_flow_in[gvert] == (H.number_of_nodes())*M['power'][gvert] - is_mapped[gvert]) #source with flow number of nodes if is mapped to power node, else if mapped to something else then sink with flow 1, else 0

        # temp_bool = model.NewBoolVar(f'TEMP{gvert}')
        model.Add(is_mapped_not_power[gvert] == is_mapped[gvert] - M['power'][gvert])
        model.Add(power_flow_out[gvert] == 0).OnlyEnforceIf(is_mapped_not_power[gvert]) #forces power flow out from sinks to be 0
        model.Add(power_flow_in[gvert] == 0).OnlyEnforceIf(M['power'][gvert]) #forces power flow into sources to be 0
    
    power_internal = dict([(gvert, model.NewBoolVar(f'{gvert}powint')) for gvert in G])
    for gvert in G:
        temp_mul = model.NewIntVar(0, (H.number_of_nodes()-1)**2, f'TEMP2{gvert}')
        model.AddMultiplicationEquality(temp_mul, [power_flow_in[gvert], power_flow_out[gvert]])
        model.Add(temp_mul - power_internal[gvert]*(H.number_of_nodes()-1)**2 <= 0)


    # for gedge in G.edges():
    #     model.AddDivisionEquality(scaled_power_flow[gedge], power_flow[gedge] + H.number_of_nodes() - 1, H.number_of_nodes())  

    for hedge in H.edges():
        for gvert in G:
            model.Add(flow_out(hedge, gvert) - flow_in(hedge, gvert) == M[hedge[0]][gvert] - M[hedge[1]][gvert])
    
    # is_internal = lambda hedge, gvert: flow_in(hedge, gvert) + flow_out(hedge, gvert) == 2
    is_internal = dict([(hedge, dict([(gvert, model.NewBoolVar(f'{hedge}{gvert}')) for gvert in G])) for hedge in H.edges()])
    for hedge in H.edges():
        for gvert in G:
            model.Add((flow_in(hedge, gvert) + flow_out(hedge, gvert) - is_internal[hedge][gvert]) <= 1) 

    for gvert in G: #each node maps to at most one edge path
        model.Add(sum([is_internal[hedge][gvert] for hedge in H.edges()]) + power_internal[gvert] + is_mapped[gvert] <= 1)

    for gedge in G.edges():
        model.Add(sum([flow[hedge][gedge] for hedge in H.edges()]) + sum([flow[hedge][gedge[::-1]] for hedge in H.edges()]) <= 1)

    solver = cp_model.CpSolver()
    status = solver.Solve(model)
    # print(status == cp_model.INFEASIBLE)
    if status == cp_model.OPTIMAL or status == cp_model.FEASIBLE:
        edge_mapping = defaultdict(list)
        power_edges = []
        for hedge in H.edges():
            for gedge in G.edges():
                for i in [1, -1]:
                    if solver.Value(flow[hedge][gedge[::i]]) == 1:
                        edge_mapping[hedge].append(gedge)
                    if solver.Value(power_flow[gedge[::i]]) != 0:
                        power_edges.append(gedge)

        node_mapping = dict([(hvert, gvert) for hvert in H for gvert in G if solver.Value(M[hvert][gvert]) == 1])
        return edge_mapping, power_edges, node_mapping
    return False
    
def draw_mapping(edge_mapping, power_edges, node_map, dim, H):
    G, _ = rect_graph(*dim)
    inv_node_map = dict([(v, H.nodes[k]['label'].strip() if k in H.nodes and len(H.nodes[k]['label'].strip()) > 0 else k) for k, v in node_map.items()])
    colored_edges = list(itertools.chain.from_iterable(edge_mapping.values()))
    # print(edge_mapping.keys())
    # print(colored_edges)
    nx.draw(G, labels=inv_node_map, edge_color=['red' if gedge in colored_edges else ('purple' if gedge in power_edges else 'blue') for gedge in G.edges()])
    plt.show()
    
H = readwrite.read_graphml('./test_graphs/power_test.graphml')
print([data for _, data in H.nodes(data=True)])
edge_map, power_edges, node_map = compactible(H, (2, 2, 2), 2)
print(edge_map)
print(node_map)
# edge_map = {('n0', 'n1'): [((2, 0, 1), (1, 0, 1))], ('n1', 'n4'): [((1, 0, 1), (0, 0, 1))], ('n2', 'n1'): [((1, 0, 0), (1, 0, 1))], ('n3', 'n4'): [((0, 1, 1), (0, 0, 1)), ((0, 1, 0), (0, 1, 1))], ('n4', 'n5'): [((0, 0, 1), (0, 0, 0))], ('n6', 'n0'): [((2, 1, 0), (2, 0, 0)), ((2, 0, 0), (2, 0, 1))], ('n7', 'n6'): [((1, 1, 0), (2, 1, 0))], ('n8', 'n6'): [((2, 2, 0), (2, 1, 0))], ('n8', 'n12'): [((2, 2, 0), (2, 2, 1))], ('n9', 'n8'): [((1, 2, 0), (2, 2, 0))], ('n10', 'n11'): [((1, 1, 1), (2, 1, 1))], ('n11', 'n6'): [((2, 1, 1), (2, 1, 0))], ('n12', 'n11'): [((2, 2, 1), (2, 1, 1))], ('n13', 'n12'): [((1, 2, 1), (2, 2, 1))], ('n14', 'n13'): [((0, 2, 0), (0, 2, 1)), ((0, 2, 1), (1, 2, 1))]}
draw_mapping(edge_map, power_edges, node_map, (2, 2, 2), H)