# TODO: ADD IN CODE TO SKELETON
# (TASK 1) DNAHasher: lineno 316
# (TASK 2) debrujin_graph_from_kmers: lineno 355
# (TASK 3) eulPath: lineno 208

# all imports
import time
import random
from typing import List
from typing import Tuple
from collections import defaultdict, deque

def main():
    test_seq_assembly();

#########################################################################################################################################
# test.cpp

class Node:
    def __init__(self, label):
        self.label = label
        self.outgoing = []
        self.num_of_incoming = 0

class DiGraph:
    def __init__(self):
        self.nodes = []

def test_and_print_message(seq, seq_truth, k, message):
    if seq == seq_truth:
        print(f"Passed {message} (assembled original sequence). Congratulations!")
    elif compare_composition(seq, seq_truth, k):
        print(f"Passed {message} (assembled a sequence of the same composition with the original sequence). Congratulations!")
    else:
        print("FAILED test 1!")

def test_1(method):
    print(f"Testing k-assembler by {method}")
    seqs_truth = [
        "aaaaaaaaaaa",
        "agcagctcagc",
        "agcagctcagg",
        random_DNA_sequence(10000, 20000)
    ]
    ks = [5, 3, 3, 20]
    for i in range(len(seqs_truth)):     
        print(f"\nExample {i}:")
        seq_truth = seqs_truth[i]
        k = ks[i]
        kmers = get_kmers(seq_truth, k, True)
        g = DiGraph()
        begin = time.time()
        if method == "k-mer pairwise comparison":
            create_deBruijn_graph_by_string_comp(kmers, g)
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
        try:
            path = find_Eulerian_path(g)
            seq = build_sequence(path, g)
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
        seq = assemble_kmers(kmers, method, "deBruijn.dot")
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

########################################################################################################################################

###############################################################################
# sequence.cpp

def random_DNA_sequence(min_length=10, max_length=10000) -> str:
    DNA = ""
    length = random.randint(min_length, max_length)
    nucleotides = "atgc"
    for n in range(length):
        DNA += nucleotides[random.randint(0, 3)]
    return DNA

def get_kmers(seq: str, k: int, randomized=True) -> List[str]:
    kmers = [seq[i:i+k] for i in range(len(seq) - k + 1)]
    if randomized:
        nkmers = len(kmers)
        for i in range(nkmers-1):
            j = random.randint(i, len(kmers)-1)
            kmers[i], kmers[j] = kmers[j], kmers[i]
    return kmers

def compare_composition(s1: str, s2: str, k: int) -> bool:
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
#################################################################

#################################################################
# k-assembler.cpp

def assemble_kmers(kmers, method, dotfile=""):
    seq = ""
    g = DiGraph()
    if method == "k-mer pairwise comparison":
        create_deBruijn_graph_by_string_comp(kmers, g)
    elif method == "k-mer hashing":
        g = debrujin_graph_from_kmers(kmers)
    else:
        raise Exception("ERROR: unknown methods!")
    if dotfile != "":
        printDOTFile(g, dotfile)
    if not has_Eulerian_path(g):
        raise Exception("ERROR: Eulerian path does not exist!")
    else:
        path = find_Eulerian_path(g)
        seq = build_sequence(path, g)
    return seq

def build_sequence(path, g):
    nodes = g.nodes
    k = len(nodes[path[0]].m_label) + 1
    seq = [""] * (k - 1 + len(path) - 1)
    seq[:k-1] = list(nodes[path[0]].m_label)
    i = k - 1
    for pos in path[1:]:
        seq[i] = nodes[pos].m_label[-1]
        i += 1
    return "".join(seq)

def printDOTFile(g, file):
    with open(file, 'w') as ofs:
        ofs.write("digraph {\n")
        ofs.write("label=\"de Bruijn graph\"\n")
        nodes = g.nodes
        for node in nodes:
            for to in node.outgoing:
                prefix = node.m_label
                suffix = nodes[to].m_label
                ofs.write(f"{prefix}->{suffix} [label={prefix}{suffix[-1]}];\n")
        ofs.write("}\n")

###########################################################################################################

##########################################################################################
# DeBruijnByStringComp.cpp

def create_deBruijn_graph_by_string_comp(kmers: List[str], g) -> None:
    nodes = []
    for kmer in kmers:
        k = len(kmer)
        prefix = kmer[:k-1]
        i = 0
        for node in nodes:
            if node.label == prefix:
                break
            i += 1
        else:
            from_node = Node(prefix)
            #from_node.m_label = prefix
            from_node.m_num_of_incoming = 0
            nodes.append(from_node)
            i = len(nodes) - 1
        from_node = nodes[i]
        suffix = kmer[1:k]
        j = 0
        for node in nodes:
            if node.label == suffix:
                break
            j += 1
        else:
            to_node = Node(suffix)
            #to_node.label = suffix
            to_node.m_num_of_incoming = 0
            nodes.append(to_node)
            j = len(nodes) - 1
        to_node = nodes[j]
        from_node.outgoing.append(j)
        to_node.m_num_of_incoming += 1
    g.nodes = nodes

##########################################################################################

##########################################################################################
# deBruijnByHash.cpp

def genomePath(kmers, apppend_last=True):
    genome = ''
    for kmer in kmers:
        genome += kmer[0]
    if apppend_last:
        genome += kmer[1:]
    return genome

#class DNAHasher:
#    def __call__(self, seq: str) -> int:
#        val = 0
        # Write a DNA sequence hash function here
        # BEGIN your code here:

        # END your code above
#        return val

class AlphabetHasher:
    def __call__(self, seq: str) -> int:
        val = 0
        max_width = 20
        for i in range(min(len(seq), max_width)):
            val = val << 5
            val += ord(seq[i].lower()) - ord('a')
        return val

class CSeqHash(defaultdict):
    def __init__(self):
        super().__init__(list)

    def insert(self, key: str, node_id: int):
        self[key].append(node_id)

def create_hash_table(kmers: List[str]) -> CSeqHash:
    ht = CSeqHash()
    node_id = 0
    for kmer in kmers:
        for j in range(2):
            key = kmer[j: len(kmer)-1]
            if key not in ht:
                ht.insert(key, node_id)
                node_id += 1
    return ht
  
def debrujin_graph_from_kmers(patterns):
    kmers = set()
    for pattern in patterns:
        kmers = set(suffix_composition(len(pattern), pattern))
    graph = defaultdict(deque)
    for kmer in patterns:
        graph[prefix(kmer)].append(suffix(kmer))
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

##########################################################################################

###########################################################################################################
# EulerPath.cpp

def source(g) -> int:
    for i in range(len(g.nodes)):
        if len(g.nodes[i].outgoing) == g.nodes[i].m_num_of_incoming + 1:
            return i
    return len(g.nodes)

def sink(g) -> int:
    for i in range(len(g.nodes)):
        if len(g.nodes[i].outgoing) + 1 == g.nodes[i].m_num_of_incoming:
            return i
    return len(g.nodes)


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
    
def balanceCount(adjacentList):
    balanced_count = defaultdict(int)
    # Look for nodes balancing
    for node in adjacentList.keys():
        for out in adjacentList.nodes[node]:
            balanced_count[node] -= 1
            balanced_count[out] += 1
    return balanced_count

def find_Eulerian_path(g):
    path = []
    cycle = []
    src = source(g)
    dest = sink(g)
    src = 0 if src >= len(g.nodes) else src
    dest = 0 if dest >= len(g.nodes) else dest
    nodes = g.nodes
    nodes[dest].outgoing.append(src)
    nodes[src].m_num_of_incoming += 1
    cycle = eulPath(g)
    pos_src, pos_dest = None, None
    for pos_dest in range(len(cycle)):
        if pos_dest != len(cycle) - 1:
            pos_src = pos_dest + 1
            if cycle[pos_src] == src and cycle[pos_dest] == dest:
                break
        else:
            break
    if pos_src is not None and pos_dest is not None:
        path += cycle[pos_src:pos_dest+1]
        path += cycle[:pos_src+1]
    else:
        raise Exception("Searching for Eulerian path has failed!")
    return path

def has_Eulerian_path(g):
    exist = True
    numSources = 0
    numSinks = 0
    for node in g.nodes:
        out = len(node.outgoing)
        _in = node.m_num_of_incoming
        if out == _in:
            continue
        elif out == _in + 1:
            numSources += 1
            if numSources > 1:
                exist = False
                break
        elif out + 1 == _in:
            numSinks += 1
            if numSinks > 1:
                exist = False
                break
        else:
            exist = False
            break
    return exist

##########################################################################################

if __name__ == "__main__":
    main()