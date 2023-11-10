#!/bin/bash


# Function to add symbols to the eand of the sequence
add_symbol_if_missing() {
    local input_string="$1"  # Получаем строку в качестве аргумента
    local symbol="$2"       # Получаем символ, который нужно добавить

    # Проверяем, что строка не пустая и символ не пустой
    if [ -n "$input_string" ] && [ -n "$symbol" ]; then
        # Проверяем, не совпадает ли последний символ с символом, который нужно добавить
        if [ "${input_string: -1}" != "$symbol" ]; then
            input_string="$input_string$symbol"
        fi
    fi

    echo "$input_string"
}

# Parse command-line parameters
while [[ $# -gt 0 ]]; do
    case "$1" in
        -r|--path-results) path_results="$2"; shift 2;;
        -g|--path-genomes) path_genomes="$2"; shift 2;;
        -t|--file-type) fasta_type="$2"; shift 2;;
        -p|--result-prefix) result_pref="$2"; shift 2;;
        -q|--query-file) query_file="$2"; shift 2;;
        *)
            echo "Unknown option: $1"
            exit 1;;
    esac
done

# Check if required parameters are missing and produce an error if they are
if [ -z "$path_results" ] || [ -z "$path_genomes" ] || [ -z "$fasta_type" ] || [ -z "$query_file" ]; then
    echo "Error: Missing required parameter(s)"
    exit 1
fi

# Add trailing slashes to path variables if missing
path_results=$(add_symbol_if_missing "$path_results")
path_genomes=$(add_symbol_if_missing "$path_genomes")

# Create the rulult folder
if [ ! -d "$path_results" ]; then
    mkdir -p "$path_results"
fi

# Check if result_pref is empty, and if so, assign the default value
if [ -z "$result_pref" ]; then
    result_pref="${path_results}bl_"
fi
# Ensure result_pref ends with an underscore if provided
result_pref=$(add_symbol_if_missing "$result_pref" "_")

# -------------------------------------------------------------
# Create an array to store blast databases
blast_databases=()

# Create blast databases for all files in the path_genomes directory with the specified file type
# and remember the names of all files
for fasta_file in "${path_genomes}"*."${fasta_type}"
do
    # Get the file name without extension to use as the database name
    db_name=$(basename "$fasta_file" .fasta)
    db_path="${path_genomes}${db_name}"

    # Check if the BLAST database already exists
    if [ ! -e "${db_path}.nhr" ] && [ ! -e "${db_path}.nin" ] && [ ! -e "${db_path}.nsq" ]; then
        # Create a BLAST database for the file only if it doesn't already exist
        makeblastdb -in "$fasta_file" -dbtype nucl -out "$db_path"
    fi
    blast_databases+=("$db_path")
done
echo ${blast_databases}

# Perform BLAST on all databases
query_file_base=$(basename "$query_file" .fasta)
for db_path in "${blast_databases[@]}"
do
    echo "${result_pref}${query_file_base}${db_name}.txt"
    blastn -query "$query_file" -db "$db_path" -out "${result_pref}${query_file_base}${db_name}.txt"
done


merged_result_file="${result_pref}${query_file_base}_merged.txt"

# Сливаем все файлы результатов в один
cat "${result_pref}${query_file_base}"*.txt > "$merged_result_file"
rm  "${result_pref}${query_file_base}"*.txt

echo "$merged_result_file"