import random
import time
from typing import List
from typing import Tuple
from collections import defaultdict, deque, Counter

def main():
    seq_truths = ["aaaaaaaaaaa", "agcagctcagc", "agcagctcag", random_DNA_sequence(11, 15)]
    for seq_truth in seq_truths:
        kmers = get_kmers(seq_truth, 3, True)
        print(kmers)
        g = debrujin_graph_from_kmers(kmers)
        print('      Debruijn graph finished')
        print(g)
        if has_Eulerian_path(g):
            print('      Has Eularian path')
            path = find_Eulerian_path(g)
            print(path)
            print('      Eularian path:', build_sequence(path))
            print('      Composition:', compare_composition(seq_truth, build_sequence(path), 3))
            print('\n')
        else:
            print('      does not have Eularian path')    
            print('\n')

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

# by doing join(path) the path will have double of every nucleotide except for the last one.
def build_sequence(path):
    seq1 = ''.join(path)
    seq2 = ''
    for i in range(0, len(seq1)):
        if (i % 2) == 0:
            seq2 += seq1[i]
    seq2 += seq1[len(seq1)-1]
    return seq2

def has_Eulerian_path(graph):
    # Calculate in-degree and out-degree for each vertex
    in_degrees = defaultdict(int)
    out_degrees = defaultdict(int)
    for u, neighbors in graph.items():
        for v in neighbors:
            out_degrees[u] += 1
            in_degrees[v] += 1
    
    # Check if the graph is connected
    visited = set()
    stack = deque(['ag'])
    while stack:
        vertex = stack.pop()
        visited.add(vertex)
        for neighbor in graph[vertex]:
            if neighbor not in visited:
                stack.append(neighbor)
    if len(visited) != len(graph):
        return False
    
    # Count the number of vertices with odd degree
    odd_degree_count = 0
    for vertex in graph:
        if out_degrees[vertex] - in_degrees[vertex] == 1:
            odd_degree_count += 1
        elif in_degrees[vertex] - out_degrees[vertex] == 1:
            odd_degree_count += 1
            start_vertex = vertex
        elif abs(out_degrees[vertex] - in_degrees[vertex]) > 1:
            return False
    
    # If the count is 0 or 2, then the graph has an Eulerian path
    return odd_degree_count == 0 or odd_degree_count == 2

if __name__ == "__main__":
    main()