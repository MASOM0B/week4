base_dir="/home/masom0b/ncbi_dataset/week_4/ncbi_dataset/data"

output_dir="/home/masom0b/ncbi_dataset/week_4/14_Genomes_genes_output"

mkdir -p "$output_dir"

for dir in "$base_dir"/*/; do

    fna_file=$(find "$dir" -name "*GCF*.fna")

    if [[ -f "$fna_file" ]]; then

        base_name=$(basename "$fna_file" .fna)

        output_file="$output_dir/${base_name}_all_genes.txt"

        python3 /home/masom0b/ncbi_dataset/week_4/genes.py "$fna_file" > "$output_file"
    fi
done
