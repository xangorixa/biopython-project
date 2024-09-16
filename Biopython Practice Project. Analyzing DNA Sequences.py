# -*- coding: utf-8 -*-
"""
Created on Mon Sep 16 15:38:32 2024

@author: xvsje
"""
from Bio.Seq import Seq
print(Seq("AGTACACTGGT"))
from Bio import Entrez
Entrez.email= "xvsj.ehs@outlook.com"
def fetch_sequence(record_id):
    handle = Entrez.efetch(db="nucleotide", id=record_id, rettype="gb", retmode="text")
    record = handle.read()
    handle.close()
    return record
    
sequence = fetch_sequence("NM_003121.5")
from Bio.Seq import Seq
from Bio.SeqUtils import gc_fraction
dna_seq = Seq(sequence)
gc_content= gc_fraction(dna_seq) * 100
motifs = ["ATG","TAG"]
found_motifs = {motif: dna_seq.find(motif) for motif in motifs}
print(f"GC content: {gc_content:.2f}%")
print("Motif positions:")
for motif, position in found_motifs.items():
    print(f"Motif {motif} found at position: {position}")
    
from Bio.Seq import Seq

sequence = """ggcaaacagcccgcccggcaccaccatgctcgccctggaggctgcacagctcgacgggccacacttcagctgtctgtacccagatggcgtcttctatgacctggacagctgcaagcattccagctaccctgattcagagggggctcctgactccctgtgggactggactgtggccccaccgtcccagccaccccctatgaagccttcgacccggcagcagccgcttttagccacccccaggctgcccagctctgctacgaaccccccacctacagccctgcagggaacctcgaactggcccccagcctggaggccccggggcctggcctccctgcataccccacggagaacttcgctagccagaccctggttcccccggcatatgccccgtaccccagccctgtgctatcagaggaggaagacttaccgttggacagccctgccctggaggtctcggacagcgagtcggatgaggccctcgtggctggccccgaggggaagggatccgaggcagggactcgcaagaagctgcgcctgtaccagttcctgctggggctactgacgcgcggggacatgcgtgagtgcgtgtggtgggtggagccaggcgccggcgtcttccagttctcctccaagcacaaggaactcctggcgcgccgctggggccagcagaaggggaaccgcaagcgcatgacctaccagaagctggcgcgcgccctccgaaactacgccaagaccggcgagatccgcaaggtcaagcgcaagctcacctaccagttcgacagcgcgctgctgcctgcagtccgccgggcctgagcacacccgaggctcccacctgcggagccgctgggggacctcacgtcccagccaggatccccctggaagaaaaagggcgtccccacactctaggtgataggacttacgcatccccaccttttggggtaaggggagtgctgccctgccataatccccaagcccagcccgggcctgtctgggattccccacttgtgcctggggtcctctgggatttctttgtcatgtacagactccctgggatcctcatgttttgggtgacaggacctatggaccactatactcggggaggcagggtagcagttcttccagaatcccaagagcttctctgggattttcttgtgatatctgattccccagtgaggcctgggacgtttttaagatcgctgtgtgtctgtaaaccctgaatctcatctggggtgggggccctgctggcaaccctgagccctgtccaaggttccctcttgtcagatctgagatttcctagttatgtctggggccctctgggagctgttatcatctcagatctcttcgcccatctatggctgtgttgtcacatctgtcccctcatttttgagatcccccaattctctggaactattctgctgcccctttttatgtgtctggagttccccaatcacatctagggctcctccaagatccttttgtcatgtctgaaatcactcttgagaggtctggggtggaggatggggagtcagtgaaatgtgtcatgtctgggccctgtcagggacacccttgttatatctgggatcctccaatcacatctgagacctcctaggctctccatctgatatgccctttcagggaccccacaaagactgagttctcatggggatcctacccttcctagtgccactccctatggccatgctgaagaccactctggccacgcgactgattttgggtgatcatggcagctccccacccatgtcatttctaaccagaagtctcaaggtcgtcacccccctgccccccaaccgaggccccggtcgctggtggtggtctctttagtgcactgtagcacttggtggtggaggtgtgagggatccacattaacagcaggccatcagctgggcaatggctcacacctgtaatcccagcactttgggaggcgaggcagggggaatggcttgaacccaggcattcaagaccagcctgggcaacataatgagacctcgtctctacaaaacataacaaaaacaattagccgagcgtgggggtgaacacctgtggtcccagctgctcaggaggctgaggtgggaggatctcttgagcccaggaagtaggaggctgtagtgagctgtaatcgtgccactgcactccagcctgggcgacagagtgagacaccgtcttaaaaacaaaaacaaggccgggcacggtggctcatgcctgttgtcccagcactttgggaggccgaggcaggcggatcacgaggtcgagagatcgagaccatcctggccaacatggtgaaaccctgtctctactaaaaatacagaaattagctgggcgtggtggcacgtgcctgtagtcccagctactcgggaggctgaggcaagagaatcgcttgaacgtgggaggcagaggttgcagtgagcctagattgtgccactgcactccagcctgggggacagagcgagactccgtctgaaaataaaaacaacaaaaacagcagaccattcaaaatagggagactttgcataatccagatttctgccttcacttaaaactttggacggtctggagagagtcggccagttttcggtggggggtggggagctggaacaggacagtagcctttcctaatgaggcatttgttctccaatctgccccagtcgctgccatccctggctatctcaccctagcagcttctcaagcctgttggctttagaccactgtataaacccagctggaactgaagcctgggtggactatggagccctggttgggacccccagggagtcaaaggctgcgggccaagaggccagaggtccttgagcctgggtgggcaggtggatctagggtgcatgacttgctgcttcccaaccttagtttgtcccttctgtgaaaaagggagagaaggaggaggaagatctcaaaaagactttccagcccagtgcggtggctcacgcctgtaatcccagcactttgggaggccgatgcaggtggatcacctgaggtaggagttcaagaccagcctgaccaacatagtgaagccccttctctactaaaaatacaaaattagctgggcgtggtggcatgtgcctgtactcccagctacttgggaggctgaggcaggagaatcgcttgaacctgggaggcggaggttgtagtgagctgagatcacaccactgcacaccagcctgggcgacaagagcgaaactccgtctcaaaaaaaaaaaactgttgcagccccgttgagcctttgacaccgcctgaaatccaccccactcccaggaggaggaggaggaaggaatgccaatgacctagagacacgagaagtccatgtggaggcacacagcagctgatggcagagcccaggctgggacctgcccttaagagaatgagtgggaagggggagggaggaagggcaggtaaaacgtcctccccagggccccctgcaacggggaaggtactttttacaaaagctatcattgtcaccctaaatgtggaataaaataagatgcatcgacgtagacaaacctcctgggaccttttgtcagggactgcaatcctgcccctccactgaggccgctggctctcagagacaccgtgacatcacgggtgatgatgagaggagttcaaagagagaattatatgctggcgcggtggctctgtaatcccaacactttggggggccaaggcaggaggatcgcttgagtacaggagtttgaaaccagcctgggcaagatagtgagatccccttcccacccgtctacaaaaaaaataaaaaattagcgggg"""
# Create a Seq object from the sequence
dna_seq = Seq(sequence)

# Pad the sequence if necessary
if len(dna_seq) % 3 != 0:
    dna_seq = dna_seq + Seq('N') * ((3 - (len(dna_seq) % 3)) % 3)

# Transcribe DNA to RNA
rna_seq = dna_seq.transcribe()

# Translate RNA to protein sequence
protein_seq = rna_seq.translate()

# Print results
print(f"RNA sequence: {rna_seq}")
print(f"Protein sequence: {protein_seq}")
import matplotlib.pyplot as plt

# Plot GC content
plt.figure(figsize=(10, 6))
plt.bar(['GC Content'], [gc_content], color='blue')
plt.ylabel('GC Content (%)')
plt.title('GC Content of DNA Sequence')
plt.show()

# Plot Motif Locations
plt.figure(figsize=(10, 6))
plt.bar(found_motifs.keys(), found_motifs.values(), color='green')
plt.xlabel('Motif')
plt.ylabel('Position in Sequence')
plt.title('Motif Positions in DNA Sequence')
plt.show()
