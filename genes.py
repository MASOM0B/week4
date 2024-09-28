import argparse
from Bio import SeqIO
from Bio.Seq import Seq

def find_orfs(sequence, min_codon_length, upstream_range, rbs_sequence):
    stop_codons = ['TAA', 'TAG', 'TGA']
    orfs = []
    min_length = min_codon_length * 3 

    for frame in range(3):
        for i in range(frame, len(sequence), 3):
            codon = sequence[i:i+3]
            if codon == 'ATG': 
                upstream_region = sequence[max(0, i - upstream_range):i]
                if rbs_sequence in upstream_region:
                    for j in range(i + 3, len(sequence), 3):
                        stop_codon = sequence[j:j+3]
                        if stop_codon in stop_codons:
                            orf = sequence[i:j+3]
                            if len(orf) >= min_length:
                                orfs.append(orf)
                            break
    return orfs

def genes():
    parser = argparse.ArgumentParser()
    parser.add_argument('file', type=str)
    parser.add_argument('--min_length', type=int, default=100)
    parser.add_argument('--upstream_range', type=int, default=20)
    parser.add_argument('--rbs_sequence', type=str, default='AGGAGG')
    args = parser.parse_args()

    for record in SeqIO.parse(args.file, "fasta"):
        genome_sequence = str(record.seq)

        forward_orfs = find_orfs(genome_sequence, args.min_length, args.upstream_range, args.rbs_sequence)
        reverse_orfs = find_orfs(str(Seq(genome_sequence).reverse_complement()), args.min_length, args.upstream_range, args.rbs_sequence)

        all_orfs = forward_orfs + reverse_orfs

        for orf in all_orfs:
            print(orf)

if __name__ == "__main__":
    genes()

