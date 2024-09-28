# Week 4: Workflow to create a genefinder tool

## 1. ORFs Forward only

Command: python3 ORFs_F_only.py /home/masom0b/ncbi_dataset/week_4/ncbi_dataset/data/GCF_000006745.1/GCF_000006745.1_ASM674v1_genomic.fna > ORFs_F_only_output.txt 

Code: 

```
import argparse #ChatGPT3.5 was used to ask which module is used in Bio-Python to search for a fasta file type as an input to be read by the code itself
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

```

Output = ORFs_F_only_output.txt

## 2. ORFs Forward and Reverse Both

Command: python3 ORFs_FR.py /home/masom0b/ncbi_dataset/week_4/ncbi_dataset/data/GCF_000006745.1/GCF_000006745.1_ASM674v1_genomic.fna > ORFs_FR_output.txt

Code:

```
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

def ORFs_FR():
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
    ORFs_FR()

```

Output= ORF_FR_output.txt

## 3. Rosalind Problem 72

Command:

Code:
 
```

```

Output=

## 4. Running ORFs Forward and Reverse for all 14 Genomes 
Command: 
 
```


```

Code: 

Output=
 
## 5. Only Lengthy ORFs
Command:

Code:


```

```

Output= 

## 6. Lengthy ORFs with a RBS 
Command: 

```

```
Code: 

Output= 
