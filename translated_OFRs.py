elp="The path to the input FASTA file"#This code is the same as the last one. It just translates all the ORFs in 6 frames in such a way that there is no redundancy and there is no * to signify the stop codon

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

def translated_ORFs_FR():

    parser = argparse.ArgumentParser()
    parser.add_argument('file', type=str)

    args = parser.parse_args()
    file_name = args.file

    fasta_seq = next(SeqIO.parse(file_name, "fasta"))
    dna_seq = str(fasta_seq.seq)

    forward_orfs = find_orfs(dna_seq)

    reverse_comp_seq = str(Seq(dna_seq).reverse_complement())
    reverse_orfs = find_orfs(reverse_comp_seq)

    all_orfs = forward_orfs + reverse_orfs

    translated_orfs = set() 
    for orf in all_orfs:
        protein = str(Seq(orf).translate(to_stop=True))
        translated_orfs.add(protein)

    for protein in translated_orfs:
        print(protein)

if __name__ == "__main__":
    translated_ORFs_FR()

