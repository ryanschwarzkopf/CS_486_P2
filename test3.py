import random
import time
from typing import List
from typing import Tuple
from collections import defaultdict, deque

def main():
    string_comp_test()

def string_comp_test():
    seqs_truth = ["aaaaaaaaaaa","agcagctcagc","agcagctcagg",random_DNA_sequence(10000, 20000)]
    ks = [5, 3, 3, 20]
    for i in range(len(seqs_truth)):
        curr = seqs_truth[i]
        k = ks[i]
        kmers = get_kmers(curr, k, True)
        begin = time.time()
        create_deBruijn_graph_by_string_comp(kmers)
        end = time.time()
        print(end - begin)

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

class Node:
    def __init__(self, label: str):
        self.m_label = label
        self.m_num_of_incoming = 0
        self.m_outgoing = []

def create_deBruijn_graph_by_string_comp(kmers):
    nodes = []
    nodesFound = []
    i = 0
    for kmerNode in kmers:
        nodesFound.append(kmerNode)
        nodes.append(Node(kmerNode))
        for kmerCompare in kmers:
            if len(kmerNode) > 1 and len(kmerCompare) > 1:
                if kmerNode[:-1] == kmerCompare[1:]:
                    nodes[i].m_num_of_incoming += 1
                if kmerNode[1:] == kmerCompare[:-1]:
                    nodes[i].m_outgoing.append(kmerCompare)
        i += 1
    return nodes

if __name__ == "__main__":
    main()