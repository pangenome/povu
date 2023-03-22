import networkx as nx
from collections import defaultdict
import subprocess

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

def dominator_tree(G, entry):
    # Initialize data structures for the algorithm
    semi = defaultdict(set)
    idom = {}
    s = {}
    bucket = defaultdict(list)

    # Depth-first search to populate 'semi' and 's' dictionaries
    def dfs(n):
        if n not in semi:
            semi[n] = len(semi)
            s[semi[n]] = n
            for m in G.successors(n):
                dfs(m)

    # Function to find the ancestor of 'v' with the lowest semi value
    def ancestor_with_lowest_semi(v):
        if v not in parent:
            return None
        u = v
        while v in parent:
            if semi[v] < semi[u]:
                u = v
            v = parent.get(v, None)
        return u

    # Function to set the parent of 'w' to 'v'
    def link(v, w):
        parent[w] = v

    # Perform the depth-first search starting from the entry node
    dfs(entry)
    n = len(semi) - 1
    parent = {}

    # Iterate through nodes in reverse order of the depth-first search
    while n > 0:
        w = s[n]
        
        # Calculate semi values for 'w' by iterating through predecessors of 'w'
        for v in G.predecessors(w):
            u = ancestor_with_lowest_semi(v)
            if u is not None and semi[u] < semi[w]:
                semi[w] = semi[u]

        # Add 'w' to the bucket of nodes with the same semi value
        bucket[s[semi[w]]].append(w)

        # Process nodes in the bucket of the parent of 'w'
        while bucket[s[n - 1]]:
            v = bucket[s[n - 1]].pop()
            u = ancestor_with_lowest_semi(v)
            if u is not None:
                # Set the immediate dominator of 'v'
                idom[v] = u if semi[u] < semi[v] else s[n - 1]

        n -= 1

    # Update the immediate dominator values for all nodes
    for n in range(1, len(s)):
        print("n: ", n)
        print("s[n]: ", s[n])
        w = s[n]
        if w in idom and idom[w] != s[semi[w]]:
            idom[w] = idom[idom[w]]

    # Build the dominator tree as a NetworkX DiGraph using the calculated idom values
    dom_tree = nx.DiGraph()
    for node, dominator in idom.items():
        if dominator is not None:
            dom_tree.add_edge(dominator, node)

    return dom_tree
 
def find_bubbles(G, entry):
    G = break_cycles_tarjan(G)
    topological_order = list(nx.topological_sort(G))
    dom_tree = dominator_tree(G, entry)
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

# Entry node (source node) of the graph
entry = 0

# Find bubbles in the graph
bubbles = find_bubbles(G, entry)
print("Bubbles found:", bubbles)

