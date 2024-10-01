# Week 4: Workflow to create a genefinder tool

This tool finds Open Reading Frames (ORFs) in FASTA files.

## Requirements
- Python 3.6+
- BioPython

## Installation
1. Clone this repository
2. Install dependencies: `pip install biopython`

Creation of git files, git repo and pushing the how I generally pushed the files
```
mkdir week_4
cd week_4
git init
touch ORFs_F_only.py README.md
nano ORFs_F_only.py
git add ORFs_F_only.py README.md
git commit -m "ORFs_F_only.py"
Username and Password (Instead of Password, had to create a token following by going into my profile, settings, Developer Settings, Personal Access Tokens, Fine-grained tokens, create token after adding a description)
git remote add origin https://github.com/MASOM0B/week4.git
git branch -M main
git push -u origin main
```

## 1. ORFs Forward only

```
Command: python3 ORFs_F_only.py /home/masom0b/ncbi_dataset/week_4/ncbi_dataset/data/GCF_000006745.1/GCF_000006745.1_ASM674v1_genomic.fna > ORFs_F_only_output.txt 
```

Python Code: 

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
```
git add ORFs_F_only.py 
Commit -m ORFs_F_only.py 
git push 
``` 

Output = ORFs_F_only_output.txt

This output was pushed along with 2 other outputs as the 8th commit

## 2. ORFs Forward and Reverse Both

Command: 

```
python3 ORFs_FR.py /home/masom0b/ncbi_dataset/week_4/ncbi_dataset/data/GCF_000006745.1/GCF_000006745.1_ASM674v1_genomic.fna > ORFs_FR_output.txt
```

Python Code:

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

Commands:  
```
nano Rosalind72.txt
python3 translated_ORFs.py Rosalind72.txt > Rosalind72_output.txt
```

Python Code:
 
```
#This code is the same as the last one. It just translates all the ORFs in 6 frames in such a way that there is no redundancy and there is no * to signify the stop codon

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

```

Output = 
```
A file Rosalind72_output.txt containing the following translations of FR_ORFs
MGWAVSPCTKLNAPNLAVGAVWCFG
MVARNTIFQEYGLRMRFLLRLYAVSFGHVQAYRCKLVSFRYS
MISLRSQTNDSTTKNRLY
MAYV
MSKRNRVKPEQKSHP
MC
MNCAVRIHSTPNIIK
MRPIWLLALYGVSADRPKAGTP
MRFLLRLYAVSFGHVQAYRCKLVSFRYS
MI
MCRHIGAS
MYSYCTVHQRQAFDTLSGVGLTSPLRPRTASPSRPHIMGWAVSPCTKLNAPNLAVGAVWCFG
MAHILGCTRRGGIPLAVLVASCTGKRLTCTDMPAHVQTKPRKAGAEISSVSRTPERWYF
MLGVLCILTAQFINGRPSTPSLASALRHP
MCRVRITPLDSRTASRSKVRITASSWSSNLASH
MPLNCVWGKLRRTRFGAISRNTIQRQQPDWAHLILCTGLLPTP
MTAPLRIVSTS
MVFRLIAPKRVRRSLPQTQFNGIRLMCRVRITPLDSRTASRSKVRITASSWSSNLASH
MYASWRHTAGCSCC
MCHGEPPDLPCERALHTRGMSLEMQG
MEG
MVFLATIRCGRGWLTRWR
MQG
MPAHVQTKPRKAGAEISSVSRTPERWYF
MGSRLIYLANGRYTPGE
MPPRRVHPRICAMGSRLIYLANGRYTPGE
MRPRGRGSSRP
MSLEMQG 
```
## 4. Running ORFs Forward and Reverse for all 14 Genomes 

Commands:

``` 
nano forloop_ORFs_FR.sh
chmod +x forloop_ORFs_FR.sh (Activating for execution)
./forloop_ORFs_FR.sh
```

For loop Code Bash Script:

```
base_dir="/home/masom0b/ncbi_dataset/week_4/ncbi_dataset/data"

output_dir="/home/masom0b/ncbi_dataset/week_4/14_Genomes_ORF_FR_output"

mkdir -p "$output_dir"  

for dir in "$base_dir"/*/; do 

    fna_file=$(find "$dir" -name "*GCF*.fna")

    if [[ -f "$fna_file" ]]; then

        base_name=$(basename "$fna_file" .fna) #

        output_file="$output_dir/${base_name}_all_ORFs_FR.txt"

        python3 /home/masom0b/ncbi_dataset/week_4/ORFs_FR.py "$fna_file" 100 > "$output_file"
    fi
done

```

Output = ORFS_FRs made in 6 frames for all the 14 Genomes having GCF strains in an output directory named 14_Genomes_ORFs_FR_output with 14 .txt files each having the same strain name that was run via the forloop embedding the python code

```
git add forloop_ORFs_FR_.sh
git add README.md
git add 14_Genomes_ORF_FR_output/
git commit -m "Task 4"
git push 
```

## 5. Only Lengthy ORFs


Command: 
```
nano long_ORFs.py
chmod +x long_ORFs.py
./forloop_longORFs.sh
```

For Loop Code Bash Script: 

```
base_dir="/home/masom0b/ncbi_dataset/week_4/ncbi_dataset/data"

output_dir="/home/masom0b/ncbi_dataset/week_4/14_Genomes_LongORFs_output"

mkdir -p "$output_dir"

for dir in "$base_dir"/*/; do

    fna_file=$(find "$dir" -name "*GCF*.fna")

    if [[ -f "$fna_file" ]]; then

        base_name=$(basename "$fna_file" .fna)

        output_file="$output_dir/${base_name}_all_long_ORFs.txt"

        python3 /home/masom0b/ncbi_dataset/week_4/long_ORFs.py "$fna_file" > "$output_file"
    fi
done
```


Python Code:

```
# 
import argparse
from Bio import SeqIO
from Bio.Seq import Seq

def find_orfs(sequence, min_codon_length):
    stop_codons = ['TAA', 'TAG', 'TGA']
    orfs = []
    min_length = min_codon_length * 3  # Minimum nucleotide length for the ORF (codons * 3)

    for frame in range(3):
        for i in range(frame, len(sequence), 3):
            codon = sequence[i:i+3]
            if codon == 'ATG':  # Start codon
                for j in range(i + 3, len(sequence), 3):
                    stop_codon = sequence[j:j+3]
                    if stop_codon in stop_codons:
                        orf = sequence[i:j+3]
                        if len(orf) >= min_length:  # Only include ORFs that meet the length requirement
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
```
Output = 14_Genomes_LongORFs_output

## 6. Lengthy ORFs with a RBS 

Commands:

```
cd /home/masom0b/ncbi_dataset/week_4
nano forloop_genes.sh
chmod +x forloop_genes.sh
./forloop_genes.sh

```

Code Bash Script: 

```
base_dir="/home/masom0b/ncbi_dataset/week_4/ncbi_dataset/data"

output_dir="/home/masom0b/ncbi_dataset/week_4/14_Genomes_genes_output"

mkdir -p "$output_dir"

for dir in "$base_dir"/*/; do

    fna_file=$(find "$dir" -name "*GCF*.fna")

    if [[ -f "$fna_file" ]]; then

        base_name=$(basename "$fna_file" .fna)

        output_file="$output_dir/${base_name}_all_genes.txt"

        python3 /home/masom0b/ncbi_dataset/week_4/genes.py "$fna_file" 100 20 > "$output_file"
    fi
done


```

Python Code: 

```
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

```

Output = 14_Genomes_genes_output dir

Adding to github:

Commands:

```
git remote add origin https://github.com/MASOM0B/week4.git
git branch -M main
git push -u origin main
```

## Pushing 
```
git remote add origin https://github.com/MASOM0B/week4.git
git branch -M main
git push -u origin main
```
