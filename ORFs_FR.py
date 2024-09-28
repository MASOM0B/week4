#This code is the same as the last code but it also considers the Reverse complements so Forward and Reverse Strand both. So, we get 3+3=6 frames in total

import argparse
from Bio import SeqIO
from Bio.Seq import Seq

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

def ORFs_F_and_R():
    parser = argparse.ArgumentParser()
    parser.add_argument('file', type=str)
    args = parser.parse_args()

    for record in SeqIO.parse(args.file, "fasta"):
        genome_sequence = str(record.seq)

        forward_orfs = find_orfs(genome_sequence)

        reverse_comp_seq = str(Seq(genome_sequence).reverse_complement())
        reverse_orfs = find_orfs(reverse_comp_seq)

        all_orfs = forward_orfs + reverse_orfs

        for orf in all_orfs:
            print(orf)

if __name__ == "__main__":
    ORFs_F_and_R()

