#!/usr/bin/env python3
# -*- coding: utf-8 -*-
import sys
from Bio import SeqIO


# get parser
#parser   = SeqIO.parse( open("/run/user/4766/gvfs/smb-share:server=store.intra.i2bc.paris-saclay.fr,share=equipes/BIM/MEMBERS/paul.roginski/DNDS/TEST/rna-NM_001180142_one_line.fa"),'fasta' )
parser   = SeqIO.parse( open(sys.argv[1]),'fasta' )
# get dictionary of sequences
seqDict = SeqIO.to_dict( parser )


# Clean codons
clean_codons = []
#assume sequence are the same length
for i in range(0,len(seqDict[list(seqDict.keys())[0]].seq),3):
    
    nucleotides = [ str(seqDict[seq].seq)[i:i+3] for seq in seqDict ]
    
    if nucleotides.count("---") / len(nucleotides) < 0.75 :
        clean_codons += [i,i+1,i+2]
          
        
with open(sys.argv[2], 'w') as out:
    
    for entry in seqDict :
        out.write(">"+entry+"\n")
        out.write((''.join( [seqDict[entry].seq[i] for i in clean_codons] )+"\n").upper())
    
    
        
    
        
