# -*- coding: utf-8 -*-
"""
Created on Mon Dec 14 05:14:59 2020

@author: myPC
"""
#!/usr/bin/python

import Bio

from Bio import SeqIO
from Bio.Seq import Seq
from collections import Counter


def no_of_records(fasta_file):
    "This function computes the number of records or sequences in the given fasta file"
    records = list(SeqIO.parse(fasta_file, "fasta"))
    return len(records)

def Seq_lengths(fasta_file):
    "This function computes the length of all sequences in a fasta file and sorts them in descending order"
    
    record_dict = {}
    
    #for loop to enter the key:value pair in above created dictionary, record_dict
    for seq_record in SeqIO.parse(fasta_file, "fasta"): 
        record_dict[seq_record.id] = len(seq_record)
    
    "to sort the values (sequence lengths) in ascending order"
    
    # Create a list of tuples sorted by index 1 i.e. value field     
    listofTuples = sorted(record_dict.items() ,  key=lambda x: x[1])
    
    # Iterate over the sorted sequence
    for elem in listofTuples :
        print(elem[0] , " ::" , elem[1] )
        
        
def orf(fasta_file):
    "This function is to answer ORF related questions about a fasta file"
    record_dict = {}
    #for loop to enter the key:value pair in above created dictionary, record_dict
    
    
    def orf_1seq(sequence):
         #"function to identify all orf and their lengths in a single sequence (for all three reading frames in a forward strand)
        all_orf=[]
        orf_length=[]
        orf_start=[]
        for rf in range(0,3):
            for i in range(rf, len(sequence), 3):
                if sequence[i:i+3] == "ATG":
                    for j in range(i+3, len(sequence), 3):
                        if sequence[j:j+3] in ["TAA","TAG","TGA"]:
                            all_orf.append(sequence[i:j+3])
                            orf_length.append(j+3-i)
                            orf_start.append(i+1)
                            break     
                       
        return orf_length
    
    
    #adding all ORFs and orf lengths (list returned from above function orf_1sequence), to each sequence in file
   
    for seq_record in SeqIO.parse(fasta_file, "fasta"): 
        
        allorf = orf_1seq(str(seq_record.seq)) 
        longest_orf = max(allorf, default=0)
        #start_orf_index =                 
        record_dict[seq_record.id] = longest_orf
        
    
    #print(max(record_dict.values(), default=0))
    return record_dict 

         
def repeats(fasta_file):
    "This function is to answer repeats related questions (4th set of problem)"
    "This function answers questions for 1 sequence, will call this function within whole fasta later"
    
    seq_file=[]
    
    def repeat_1seq(sequence):
        
        all_repeats=[]
        dict_freq = {}
        n=7
        
        for i in range(len(sequence)-n+1):
            all_repeats.append(sequence[i:i+n])
            
        for repeat in all_repeats:
            freq = all_repeats.count(repeat)
            if freq > 0:
                dict_freq[repeat] = freq
                 
        return dict_freq   
    
    for seq_record in SeqIO.parse(fasta_file, "fasta"):
        
        dict_freq_1seq = repeat_1seq(str(seq_record.seq)) 
        #max_repeat = max(dict_freq_1seq, default=0)                 
        #record_dict[seq_record.id] = 
        
        seq_file.append(dict_freq_1seq)
    
        #freq_file = repeat_1seq(str(seq_record.seq))
        
    mergeD = Counter(seq_file[0])
    
    for l in range(1, len(seq_file)):
        mergeD = Counter(mergeD) + Counter(seq_file[l])       
                                   
    
    return mergeD       
                        
       
    
    
   
      
        
                
                        
                        
            
            

    
               