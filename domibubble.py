import networkx as nx
from collections import defaultdict
import subprocess
from functools import reduce

def break_cycles_tarjan(G):
    index_counter = [0]
    index = {}
    lowlink = {}
    stack = []
    new_G = G.copy()

    def strongconnect(node):
        index[node] = index_counter[0]
        lowlink[node] = index_counter[0]
        index_counter[0] += 1
        stack.append(node)

        for successor in G.successors(node):
            if successor not in index:
                strongconnect(successor)
                lowlink[node] = min(lowlink[node], lowlink[successor])
            elif successor in stack:
                lowlink[node] = min(lowlink[node], index[successor])

        if lowlink[node] == index[node]:
            connected_component = set()

            while True:
                successor = stack.pop()
                connected_component.add(successor)
                if successor == node: break

            if len(connected_component) > 1:
                max_index_edge = None
                max_index = -1

                for node in connected_component:
                    for successor in new_G.successors(node):
                        if successor in connected_component:
                            edge_index = index[node] + index[successor]
                            if edge_index > max_index:
                                max_index = edge_index
                                max_index_edge = (node, successor)

                if max_index_edge:
                    # reverse the edge with the highest index
                    new_G.remove_edge(*max_index_edge)
                    new_G.add_edge(max_index_edge[1], max_index_edge[0])

    for node in G:
        if node not in index:
            strongconnect(node)

    return new_G

def dominator_tree(graph):
    """
    Generate the dominator tree of a given DAG (Directed Acyclic Graph).

    Parameters:
    graph (nx.DiGraph): A directed acyclic graph represented using networkx.

    Returns:
    dtree (nx.DiGraph): The dominator tree of the input graph.
    """
    # Step 1: Find all the head nodes (nodes with in-degree 0)
    head_nodes = [node for node, in_degree in graph.in_degree() if in_degree == 0]
    # add a new node to the graph that dominates all the head nodes
    graph.add_node(-1)
    for head_node in head_nodes:
        graph.add_edge(-1, head_node)

    # Step 2: Initialize the dominator tree with the same nodes as the input graph
    dtree = nx.DiGraph()
    dtree.add_nodes_from(graph.nodes)

    # Step 3: Create a mapping from each node to its immediate dominator
    idom = defaultdict(lambda: None)

    # Step 4: Define a function to compute the intersection of dominators
    def intersect(b1, b2):
        finger1, finger2 = b1, b2
        while finger1 != finger2:
            while (finger1 is not None) and (finger2 is None or finger1 > finger2):
                finger1 = idom[finger1]
            while (finger2 is not None) and (finger1 is None or finger2 > finger1):
                finger2 = idom[finger2]
        return finger1

    # Step 5: Compute the immediate dominators
    for node in nx.topological_sort(graph):
        if node != -1:
            idom[node] = reduce(intersect, [pred for pred in graph.predecessors(node)])

    # Step 6: Add edges between nodes and their immediate dominators in the dominator tree
    for node, dom_node in idom.items():
        if dom_node is not None:
            dtree.add_edge(dom_node, node)

    return dtree
 
def find_bubbles(G):
    G = break_cycles_tarjan(G)
    # write the graph to a file as dot
    nx.drawing.nx_pydot.write_dot(G, 'broken.dot')
    # run dot to generate a png
    subprocess.call(['dot', '-Tpng', 'broken.dot', '-o', 'broken.png'])
    topological_order = list(nx.topological_sort(G))
    dom_tree = dominator_tree(G)
    bubbles = []
    # write the dominator tree to a dot file
    nx.drawing.nx_pydot.write_dot(dom_tree, 'dom_tree.dot')
    # use dot to convert the dot file to a png file
    subprocess.call(['dot', '-Tpng', 'dom_tree.dot', '-o', 'dom_tree.png'])
    for parent in dom_tree:
        if G.in_degree(parent) > 1:
            children = list(dom_tree.successors(parent))
            max_child_order = max(topological_order.index(c) for c in children)
            bubble_end = topological_order[max_child_order]
            bubbles.append((parent, bubble_end))

    return bubbles

# Example usage:

# Create a directed graph with cycles
G = nx.DiGraph()
G.add_edges_from([(0, 1), (1, 2), (2, 3), (3, 4), (4, 5), (2, 5), (5, 6), (6, 7), (7, 8), (8, 9), (9, 6), (9, 10)])
# write the graph to a dot file
nx.drawing.nx_pydot.write_dot(G, 'cyclic.dot')
# run dot to create a png
subprocess.call(['dot', '-Tpng', 'cyclic.dot', '-o', 'cyclic.png'])

# Find bubbles in the graph
bubbles = find_bubbles(G)
print("Bubbles found:", bubbles)

