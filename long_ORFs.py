import argparse
from Bio import SeqIO
from Bio.Seq import Seq

def find_orfs(sequence, min_codon_length):
    stop_codons = ['TAA', 'TAG', 'TGA']
    orfs = []
    min_length = min_codon_length * 3  

    for frame in range(3):
        for i in range(frame, len(sequence), 3):
            codon = sequence[i:i+3]
            if codon == 'ATG':  
                for j in range(i + 3, len(sequence), 3):
                    stop_codon = sequence[j:j+3]
                    if stop_codon in stop_codons:
                        orf = sequence[i:j+3]
                        if len(orf) >= min_length: 
                            orfs.append(orf)
                        break
    return orfs

def long_ORFs():
    parser = argparse.ArgumentParser()
    parser.add_argument('file', type=str)
    parser.add_argument('--min_length', type=int, default=100)
    args = parser.parse_args()

    for record in SeqIO.parse(args.file, "fasta"):
        genome_sequence = str(record.seq)

        forward_orfs = find_orfs(genome_sequence, args.min_length)
        reverse_comp_seq = str(Seq(genome_sequence).reverse_complement())
        reverse_orfs = find_orfs(reverse_comp_seq, args.min_length)

        all_orfs = forward_orfs + reverse_orfs

        for orf in all_orfs:
            print(orf)

if __name__ == "__main__":
    long_ORFs()
