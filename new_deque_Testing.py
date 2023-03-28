#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Mar 20 13:01:56 2023

@author: cynthiavannesatovar
"""

# This Python code defines several functions to construct a genome sequence from 
# a list of k-mers using a De Bruijn graph approach


# importing dependencies ******************************************************
import random
from collections import defaultdict, deque

'''
Get to know deque:
Appending and popping items from the right end of a Python list is usually 
efficient (O(1)), but can become slower (O(n)) when the list needs to 
reallocate memory. Appending and popping items from the left end of a list is 
generally inefficient (O(n)). While Python lists can be used as stacks and 
queues, their performance issues can impact overall application performance. 
To address these issues, Python's deque data type was introduced, which provides 
memory-efficient and fast append and pop operations on both ends of the data
structure.
'''

'''
The main function is construct_sequence(patterns), which takes a list of k-mers
as input and returns the reconstructed genome sequence. It does this by first
constructing the De Bruijn graph from the k-mers using the function 
debrujin_graph_from_kmers(patterns), and then finding the Eulerian path through
the graph using the function eulPath(dict). Finally, the reconstructed genome
sequence is obtained by concatenating the first k-1 characters of each k-mer in 
the Eulerian path using the function genomePath(kmers).
'''

def construct_sequence(patterns):
    graph = debrujin_graph_from_kmers(patterns)
    balanced_count = balanceCount(graph)
    path = deque()
    if has_Eulerian_path(balanced_count):
        path = eulPath(graph,balanced_count)
        print("Eulerian Path Exists")
    else:
        return("No Eulerian Path Exists")        
    genome = genomePath(path)
    return genome

'''
The function debrujin_graph_from_kmers(patterns) constructs the De Bruijn graph
by first generating all possible k-mers from the input list and then creating a 
dictionary where each key is a k-mer and the corresponding value is a list of
k-mers that overlap with the key by k-1 characters.
'''

'''
## To optimize the debrujin_graph_from_kmers(patterns) function, we can:

Use a defaultdict object from the collections module to replace the dict object. 
This will automatically initialize keys with an empty deque object, enabling us 
to use a deque to store the list of outgoing edges for each node in the graph.
Use the append() method to add each outgoing edge to the deque instead of 
dict[prefix(kmer)].append(suffix(kmer)). This will append a new item to a deque 
object instead of a Python list object, which can improve performance.
Change the loop to add all suffixes to the set of kmers using 
kmers.add(suffix(kmer)) instead of creating a new list kmers and concatenating 
it with the suffix composition of each pattern. This approach will avoid 
creating a list of unnecessary strings.
'''


def debrujin_graph_from_kmers(patterns):
    kmers = set()
    for pattern in patterns:
        kmers |= set(suffix_composition(len(pattern), pattern))
    graph = defaultdict(deque)
    for kmer in patterns:
        graph[prefix(kmer)].append(suffix(kmer))
    return graph

def genomePath(kmers, apppend_last=True):
    genome = ''
    for kmer in kmers:
        genome += kmer[0]
    if apppend_last:
        genome += kmer[1:]
    return genome

'''
The function eulPath(dict) takes a dictionary representation of the 
De Bruijn graph as input and returns the Eulerian path through the graph as a 
list of nodes. This is done by iteratively traversing the graph and adding nodes
to a stack until a dead end is reached, at which point the last node is added 
to the path and removed from the stack.
'''

'''
## To optimize the eulPath(dict) function, we can take the following steps:

Replace the stack object with a deque object from the collections module, 
which will allow us to use a deque as a stack.
Change stack.append() to stack.appendleft() to add elements to the left end of 
the deque instead of the right end, allowing the algorithm to use a deque as a
queue. Replace the path object with a deque object, which will allow us to use 
a deque as a queue. Replace path.append() with path.appendleft() to add
 elements to the left end of the deque instead of the right end, allowing the 
 algorithm to use a deque as a stack. Remove the [::-1] slice at the end of
 return path[::-1], as we will be using a deque to store the final path.
'''

def eulPath(graph, balanced_count):
    dictionary = deque()
    #print("BALANCED COUNT ITEMS")
    #print(balanced_count.items())
    dictionary.appendleft([k for k, v in balanced_count.items() if v == -1][0])
    path = deque()
    while dictionary:
        u_v = dictionary[0]
        try:
            w = graph[u_v][0]
            dictionary.appendleft(w)
            graph[u_v].popleft()
        except:
            path.appendleft(dictionary.popleft())
    return path

'''
suffix_composition(k, text, uniq=False), which generates all k-length suffixes 
of a given string
'''

def suffix_composition(k, text):
    kmers = []
    for i in range(len(text) + 1 - k):
        kmers.append(text[i:i+k-1])
    return sorted(kmers)
'''
 balanceCount(adjacentList), which computes the in-degree and out-degree of
 each node in the graph.
'''

def balanceCount(adjacentList):
    # create a set of all nodes in the graph
    all_nodes = set(adjacentList.keys())
    #print(all_nodes)
    #print(adjacentList)
    
    for nodes in adjacentList.values():
        all_nodes.update(nodes)

    # create a dictionary to hold the balanced counts
    balanced_count = dict.fromkeys(all_nodes, 0)
    
    if len(all_nodes) == 1:
        #print("Single Node Exception")
        for node in all_nodes:
            balanced_count[node] = -1
        return balanced_count

    # iterate over each node in the adjacency list
    for node in all_nodes:
        # subtract the out-degree of the node from its balanced count
        balanced_count[node] -= len(adjacentList[node])
        # iterate over each outgoing edge from the node
        for out in adjacentList[node]:
            # add 1 to the balanced count of the node at the other end of the edge
            try:
                balanced_count[out] += 1
            except KeyError:
                balanced_count[out] = 1
    return balanced_count
    

def suffix(string):
    return string[1:]

def prefix(string):
    return string[:-1]
    
def get_kmers(seq, k, randomized=True):
    kmers = [seq[i:i+k] for i in range(len(seq) - k + 1)]
    #print(kmers)
    if randomized:
        nkmers = len(kmers)
        for i in range(nkmers-1):
            j = random.randint(i, len(kmers)-1)
            kmers[i], kmers[j] = kmers[j], kmers[i]
    return kmers
    
def random_DNA_sequence(min_length=10, max_length=10000):
    DNA = ""
    length = random.randint(min_length, max_length)
    nucleotides = "atgc"
    for n in range(length):
        DNA += nucleotides[random.randint(0, 3)]
    return DNA
    
def has_Eulerian_path(balanced_count):
    for k, v in balanced_count.items():
        if v == -1:
            return True
    return False    
    
def main():
    seq_truths = ["aaaaaaaaaaa", "agcagctcagc", "agcagctcag", random_DNA_sequence(11, 15)]
    
    for i in seq_truths:
        print(i)
        kmers = get_kmers(i,3,True)
        print(kmers)
        print(construct_sequence(kmers))

if __name__ == "__main__":
    main()
    #data = "".join(open('./text1.txt')).split()
    #print(data)
    #print(data[1:])
    #print(construct_sequence(data[1:]))