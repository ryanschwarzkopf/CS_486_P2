#!/usr/bin/env python3
# -*- coding: utf-8 -*-
'''
title: "Bioinformatics P2 "
author: "Cynthia Tovar"
date: "2023-03-05"

'''

# This Python code defines several functions to construct a genome sequence from 
# a list of k-mers using a De Bruijn graph approach



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
    return genomePath(eulPath(debrujin_graph_from_kmers(patterns)))

'''
The function debrujin_graph_from_kmers(patterns) constructs the De Bruijn graph
by first generating all possible k-mers from the input list and then creating a 
dictionary where each key is a k-mer and the corresponding value is a list of
k-mers that overlap with the key by k-1 characters.
'''
def debrujin_graph_from_kmers(patterns):
    kmers = []
    for pattern in patterns:
        kmers = kmers + suffix_composition(len(pattern), pattern, uniq=True)
    kmers = set(kmers)
    dict = {}
    for kmer1 in kmers:
        dict[kmer1] = []
    for kmer in patterns:
        dict[prefix(kmer)].append(suffix(kmer))
    return dict


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
def eulPath(dict):
    stack=[]
    balanced_count = balanceCount(dict)
    stack.append([k for k, v in balanced_count.items() if v==-1][0])
    path = []
    while stack != []:
        u_v = stack[-1]
        try:
            w = dict[u_v][0]
            stack.append(w)
            dict[u_v].remove(w)
        except:
            path.append(stack.pop())
    return path[::-1]
'''
suffix_composition(k, text, uniq=False), which generates all k-length suffixes 
of a given string
'''

def suffix_composition(k, text, uniq=False):
    kmers = []
    for i in range(len(text)+1-k):
        kmers.append(text[i:i+k-1])
    if uniq:
        return sorted(list(kmers))
    else:
        return sorted(kmers)
'''
 balanceCount(adjacentList), which computes the in-degree and out-degree of
 each node in the graph.
'''

def balanceCount(adjacentList): # is the graph
    balanced_count = dict.fromkeys(adjacentList.keys(), 0)
    # Look for nodes balancing
    for node in adjacentList.keys():
        for out in adjacentList[node]:
            balanced_count[node] -= 1
            try:
                balanced_count[out] += 1
            except:
                balanced_count[out] = 1
    return balanced_count


def suffix(string):
    return string[1:]


def prefix(string):
    return string[0:-1]


if __name__ == "__main__":
    data = "".join(open('/Users/cynthiavannesatovar/Desktop/NMSU/spring2023/bioinformaticCourse/Project2/text1.txt')).split()
    print(construct_sequence(data[1:]))

