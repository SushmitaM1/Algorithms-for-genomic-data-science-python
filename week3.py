# -*- coding: utf-8 -*-
"""
Created on Sat Dec 26 01:08:44 2020

@author: myPC
"""
import kmer_index

def readFastq(filename):
    sequences = []
    qualities = []
    with open(filename) as fh:
        while True:
            fh.readline()  # skip name line
            seq = fh.readline().rstrip()  # read base sequence
            fh.readline()  # skip placeholder line
            qual = fh.readline().rstrip() # base quality line
            if len(seq) == 0:
                break
            sequences.append(seq)
            qualities.append(qual)
    return sequences, qualities

#Programmin homework week 3, Q3, this function is too slow for a long data with many reads, write a faster function
#Say we are concerned only with overlaps that 
#(a) are exact matches (no differences allowed), and (b) are at least k bases long. To make an overlap graph, we could call overlap(a, b, min_length=k) this function on every possible pair of reads from the dataset.  Unfortunately, that will be very slow!
def overlap(a, b, min_length):
    """ Return length of longest suffix of 'a' matching
        a prefix of 'b' that is at least 'min_length'
        characters long.  If no such overlap exists,
        return 0. """
    start = 0  # start all the way at the left
    while True:
        start = a.find(b[:min_length], start)  # look for b's prefix in a
        if start == -1:  # no more occurrences to right
            return 0
        # found occurrence; check for full suffix/prefix match
        if b.startswith(a[start:]):
            return len(a)-start
        start += 1  # move just past previous match
        
        
#Too slow, use below this code for faster one
def overlap_modified(reads, min_length):
    
    kmer_reads = {}
    
    #r_suffix = list()
    final_reads = set()
    read_suffix = set()    
    for r in reads:
        index = kmer_index.Index(r, min_length)
        for i in range(len(index.index)):
            if index.index[i][0] not in kmer_reads:
                kmer_reads[index.index[i][0]] = set([r])
            else:
                kmer_reads[index.index[i][0]].update([r])
    
    #return kmer_reads.keys()
    #(2) Now, for each read a, we find all overlaps involving a suffix of a.
    #we take a's length-k (min_length) suffix, find all reads containing that k-mer
    
    for r in reads:
                
        read_list = list(kmer_reads.get(r[-min_length : ]))
              
        if len(read_list) > 1:
        
            for i in range(len(read_list)):
                for j in range(len(read_list)):
                    if j != i:                                           
                        overlap_result = overlap(read_list[i], read_list[j], min_length)
                                                
                        if overlap_result > 0:
                            #(read_list[i], read_list[j])   #tuple with pair of reads
                                #print(reads_pair)
                            final_reads.update([(read_list[i], read_list[j])])
                            read_suffix.update([read_list[i]])     
           
    return read_suffix


#copied incomplete code from forum, modified it to give answers required, so much faster!!:
from collections import *
def Correct_overlap_all_pairs(reads, min_length):
    
    #kmer_dict = createKmersFromReadS(reads, min_length) #create dict with kmer key/and set of reads with that kmer value. User has to create this function
    kmer_dict = {}
    #r_suffix = list()
    #final_reads = set()
    #read_suffix = set()    
    for r in reads:
        index = kmer_index.Index(r, min_length)
        for i in range(len(index.index)):
            if index.index[i][0] not in kmer_dict:
                kmer_dict[index.index[i][0]] = set([r])
            else:
                kmer_dict[index.index[i][0]].update([r])
    
    overlap_graph = defaultdict(set) #read is key, set of overlapping reads are value
    val_list_total = []
    
    for read in reads:
        #create suffix for this read
        read_suffix = read[-min_length: ]
        
        
        #extract set of all reads containing this kmer/suffix
        read_set = kmer_dict[read_suffix].copy()
        
        assert(len(read_set) > 0) # check that the set isnt empty
        
        read_set.remove(read) #remove the read so we dont compare it with itself
        

        #THIS WORKS
        for compar_read in read_set:
            if overlap(read, compar_read, min_length) > 0:
                overlap_graph[read].add(compar_read)
         
        val_list_total.append(len(overlap_graph[read]))
    
    total_edges = 0
    for i in val_list_total:
        total_edges += i
    
    total_nodes_withasuffix = 0
    for i in val_list_total:
        if i != 0:
            total_nodes_withasuffix += 1
        
            
        
    
    return total_edges, total_nodes_withasuffix       
        