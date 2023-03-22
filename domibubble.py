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
    semi = defaultdict(set)
    idom = {}
    s = {}
    bucket = defaultdict(list)

    def dfs(n):
        print("dfs", n)
        if n not in semi:
            semi[n] = len(semi)
            s[semi[n]] = n
            for m in G.successors(n):
                dfs(m)

    def ancestor_with_lowest_semi(v):
        print("ancestor_with_lowest_semi", v)
        u = v
        while v in s:
            if semi[v] < semi[u]:
                u = v
            v = parent[v]
        return u

    def link(v, w):
        parent[w] = v

    dfs(entry)
    n = len(semi) - 1
    parent = {}

    while n > 0:
        w = s[n]
        for v in G.predecessors(w):
            u = ancestor_with_lowest_semi(v)
            if semi[u] < semi[w]:
                semi[w] = semi[u]

        bucket[s[semi[w]]].append(w)

        while bucket[s[n - 1]]:
            v = bucket[s[n - 1]].pop()
            u = ancestor_with_lowest_semi(v)
            idom[v] = u if semi[u] < semi[v] else s[n - 1]

        n -= 1

    for n in range(1, len(s)):
        w = s[n]
        if idom[w] != s[semi[w]]:
            idom[w] = idom[idom[w]]

    return idom

def find_bubbles(G, entry):
    G = break_cycles_tarjan(G)
    nx.drawing.nx_pydot.write_dot(G, 'broken.dot')
    # run dot to create a png
    subprocess.call(['dot', '-Tpng', 'broken.dot', '-o', 'broken.png'])
    topological_order = list(nx.topological_sort(G))
    dom_tree = dominator_tree(G, entry)
    bubbles = []

    for parent, child in dom_tree.items():
        if G.in_degree(parent) > 1:
            max_child_order = max(topological_order.index(c) for c in G.successors(parent))
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

