import random
import time
from typing import List
from typing import Tuple
from collections import defaultdict, deque

def main():
    seq_truth = "aaaaaaaaaaa"
    kmers = get_kmers(seq_truth, 3, True)
    g = debrujin_graph_from_kmers(kmers)
    print('debruijn graph finished')
    if has_Eulerian_path(g):
        print('Has Eularian path')
        path = find_Eulerian_path(g)
    else:
        print('does not have Eularian path')


def random_DNA_sequence(min_length=10, max_length=10000):
    DNA = ""
    length = random.randint(min_length, max_length)
    nucleotides = "atgc"
    for n in range(length):
        DNA += nucleotides[random.randint(0, 3)]
    return DNA

def get_kmers(seq, k, randomized=True):
    kmers = [seq[i:i+k] for i in range(len(seq) - k + 1)]
    if randomized:
        nkmers = len(kmers)
        for i in range(nkmers-1):
            j = random.randint(i, len(kmers)-1)
            kmers[i], kmers[j] = kmers[j], kmers[i]
    return kmers

def compare_composition(s1, s2, k):
    same_composition = True
    if len(s1) != len(s2):
        same_composition = False
    elif s1 == s2:
        same_composition = True
    else:
        composition1 = get_kmers(s1, k)
        composition1.sort()
        composition2 = get_kmers(s2, k)
        composition2.sort()
        same_composition = composition1 == composition2
    return same_composition

def debrujin_graph_from_kmers(patterns):
    kmers = set()
    for pattern in patterns:
        kmers = set(suffix_composition(len(pattern), pattern))
    graph = defaultdict(deque)
    for kmer in patterns:
        graph[prefix(kmer)].append(suffix(kmer))
    return graph

def create_deBruijn_graph_by_string_comp(kmers):
    graph = {}
    for kmer in kmers:
        # make the prefix and suffix
        prefix = kmer[:-1]
        suffix = kmer[1:]
        # if prefix isn't in the graph add it and add the suffix, else add the suffix to the right prefix
        if prefix not in graph:
            graph[prefix] = []
        graph[prefix].append(suffix)
    return graph

def suffix_composition(k, text):
    kmers = []
    for i in range(len(text) + 1 - k):
        kmers.append(text[i:i+k-1])
    return sorted(kmers)
    
def suffix(string):
    return string[1:]

def prefix(string):
    return string[:-1]

def find_Eulerian_path(graph):
    stack = deque()
    balanced_count = balanceCount(graph)
    stack.appendleft([k for k, v in balanced_count.items() if v == -1][0])
    path = deque()
    while stack:
        u_v = stack[0]
        try:
            w = graph[u_v][0]
            stack.appendleft(w)
            graph[u_v].popleft()
        except:
            path.appendleft(stack.popleft())
    return path

def balanceCount(adjacentList):
    balanced_count = defaultdict(int)
    # Look for nodes balancing
    for node in adjacentList.keys():
        for out in adjacentList[node]:
            balanced_count[node] -= 1
            balanced_count[out] += 1
    return balanced_count

def build_sequence(path):
    return ''.join(path)

def has_Eulerian_path(graph):
    # Check if the graph is strongly connected
    visited = set()
    start = list(graph.keys())[0]
    dfs(start, graph, visited)
    if len(visited) != len(graph):
        return False

    # Count the in-degree and out-degree of each node
    indegree = {node: 0 for node in graph}
    outdegree = {node: 0 for node in graph}
    for node, neighbors in graph.items():
        outdegree[node] = len(neighbors)
        for neighbor in neighbors:
            indegree[neighbor] += 1

    # Check if the graph has at most two nodes with odd degree
    num_odd = sum(deg % 2 == 1 for node, deg in outdegree.items())
    num_odd += sum(deg % 2 == 1 for node, deg in indegree.items())
    return num_odd <= 2

def dfs(node, graph, visited):
    if node not in visited:
        visited.add(node)
        for neighbor in graph[node]:
            dfs(neighbor, graph, visited)


if __name__ == "__main__":
    main()