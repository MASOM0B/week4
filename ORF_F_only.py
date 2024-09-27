#This file will create ORFs only considering the forward strand i.e 3 reading frames

import argparse #ChatGPT3.5 was used to ask which module is used in Bio-Python to search for a fasta file type as an input to be read by the code itself
from Bio import SeqIO

import argparse
from Bio import SeqIO

def find_orfs(sequence):
    stop_codons = ['TAA', 'TAG', 'TGA']
    orfs = []

    for frame in range(3):
        for i in range(frame, len(sequence), 3):
            codon = sequence[i:i+3]
            if codon == 'ATG':
                for j in range(i + 3, len(sequence), 3):
                    stop_codon = sequence[j:j+3]
                    if stop_codon in stop_codons:
                        orfs.append(sequence[i:j+3])
                        break
    return orfs

def ORFs_F_only():
    #ChatGPT 3.5 was used to ask the correct way to call arg.parse as I was facing a syntax error
    parser = argparse.ArgumentParser()
    parser.add_argument('file', type=str)
    args = parser.parse_args()

    for record in SeqIO.parse(args.file, "fasta"):
        genome_sequence = str(record.seq)

        orfs = find_orfs(genome_sequence)
        for orf in orfs:
            print(orf)

if __name__ == "__main__":
    ORFs_F_only()
