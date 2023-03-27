import random
import time
from typing import List
from typing import Tuple
from collections import defaultdict, deque

def main():
    test_seq_assembly()

def test_and_print_message(seq, seq_truth, k, message):
    if seq == seq_truth:
        print(f"Passed {message} (assembled original sequence). Congratulations!")
    elif compare_composition(seq, seq_truth, k):
        print(f"Passed {message} (assembled a sequence of the same composition with the original sequence). Congratulations!")
    else:
        print("FAILED test 1!")

def test_1(method):
    print(f"Testing k-assembler by {method}")
    seqs_truth = ["aaaaaaaaaaa","agcagctcagc","agcagctcagg",random_DNA_sequence(10000, 20000)]
    ks = [5, 3, 3, 20]
    for i in range(len(seqs_truth)):     
        print(f"\nExample {i}:")
        seq_truth = seqs_truth[i]
        k = ks[i]
        kmers = get_kmers(seq_truth, k, True)
        begin = time.time()
        if method == "k-mer pairwise comparison":
            g = create_deBruijn_graph_by_string_comp(kmers)
        elif method == "k-mer hashing":
            g = debrujin_graph_from_kmers(kmers)
        else:
            print("ERROR: unknown methods!") 
        end = time.time()
        elapsed_secs = end - begin
        print(f"Elapsed time for building de Bruijn graph: {elapsed_secs}")
        if has_Eulerian_path(g):
            print("Passed test for existence of Eulerian path. Congratulations!")
        else:
            print("Failed test for existence of Eulerian path!")
            return 0
        try:
            path = eulPath(g)
            seq = build_sequence(path)
            message = f"Test 1 Example {i}"
            test_and_print_message(seq, seq_truth, k, message)
        except Exception as e:
            print(f"ERROR: {e}")

def test_2(method):
    seq_truth = random_DNA_sequence()
    k = 10
    kmers = get_kmers(seq_truth, k)
    try:
        seq = assemble_kmers(kmers, method)
        test_and_print_message(seq, seq_truth, k, "Test 2")
    except Exception as e:
        print(e)

def test_3(method):
    seq_truth = random_DNA_sequence(15, 15)
    k = 4
    print(f"Sequence: {seq_truth}")
    kmers = get_kmers(seq_truth, k)
    print("kmers:")
    for kmer in kmers:
        print(kmer)
    try:
        seq = assemble_kmers(kmers, method)
        test_and_print_message(seq, seq_truth, k, "Test 3")
    except Exception as e:
        print(e)

def test_seq_assembly():
    methods = [
        "k-mer pairwise comparison",
        "k-mer hashing"
    ]
    for method in methods:   
        print("-----------")
        test_1(method)
        print()
        test_2(method)
        print()
    print("-----------")
    test_3("k-mer hashing")

def assemble_kmers(kmers, method):
    seq = ""
    if method == "k-mer pairwise comparison":
        g = create_deBruijn_graph_by_string_comp(kmers)
    elif method == "k-mer hashing":
        g = debrujin_graph_from_kmers(kmers)
    else:
        raise Exception("ERROR: unknown methods!")
    if not has_Eulerian_path(g):
        raise Exception("ERROR: Eulerian path does not exist!")
    else:
        path = eulPath(g)
        seq = build_sequence(path)
    return seq

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

def genomePath(kmers, apppend_last=True):
    genome = ''
    for kmer in kmers:
        genome += kmer[0]
    if apppend_last:
        genome += kmer[1:]
    return genome

def debrujin_graph_from_kmers(patterns):
    kmers = set()
    for pattern in patterns:
        kmers = set(suffix_composition(len(pattern), pattern))
    graph = defaultdict(deque)
    for kmer in patterns:
        graph[prefix(kmer)].append(suffix(kmer))
    return graph

def create_deBruijn_graph_by_string_comp(kmers):
    graph = defaultdict(deque)
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

def eulPath(graph):
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

def has_Eulerian_path(graph):
    # Check if the graph is strongly connected
    visited = set()
    start = list(graph.keys())[0]
    visited = dfs(start, graph, visited)
    #print(visited)
    #print(graph)
    if len(visited) != len(graph.keys()):
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
            visited = dfs(neighbor, graph, visited)
        return visited
    else:
        return visited


if __name__ == "__main__":
    main()