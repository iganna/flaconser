#!/bin/bash


# ./one_search.sh -r ..//tir -g ../../pacbio/pb_updated/ -t fasta -q ../candidates/tir.fasta -m m.txt -n 30                                                                         

# -------------------------------------------------------------
#            ERROR HANDLING BLOCK
# -------------------------------------------------------------

# Exit immediately if any command returns a non-zero status
set -e

# Keep track of the last executed command
trap 'last_command=$current_command; current_command=$BASH_COMMAND' DEBUG

# Define a trap for the EXIT signal
trap 'catch $?' EXIT

# Function to handle the exit signal
catch() {
    # Check if the exit code is non-zero
    if [ $1 -ne 0 ]; then
        echo "\"${last_command}\" command failed with exit code $1."
    fi
}

# -------------------------------------------------------------
#             FUNCTIONS
# -------------------------------------------------------------

source func_search.sh

# -------------------------------------------------------------
#                PARAMETERS
# -------------------------------------------------------------

# Parse command-line parameters
while [[ $# -gt 0 ]]; do
    case "$1" in
        -r|--path-results) path_results="$2"; shift 2;;
        -g|--path-genomes) path_genomes="$2"; shift 2;;
        -t|--file-type) fasta_type="$2"; shift 2;;
        -p|--result-prefix) result_pref="$2"; shift 2;;
        -q|--query-file) query_file="$2"; shift 2;;
        -m|--merged-file) merged_file="$2"; shift 2;;
        -n|--n-cores) n_cores="$2"; shift 2;;
        *)
            echo "Unknown option: $1"
            exit 1;;
    esac
done


# Number of cores
if ! [[ "$n_cores" =~ ^[0-9]+$ ]]; then
    n_cores=1
fi

# ---- Check missimg parameters

check_missing_variable "path_results"
check_missing_variable "path_genomes"
check_missing_variable "fasta_type"
check_missing_variable "query_file"
check_missing_variable "merged_file"


# ---- Add trailing slashes to path variables if missing

path_results=$(add_symbol_if_missing "$path_results" "/")
path_genomes=$(add_symbol_if_missing "$path_genomes" "/")


# ---- Create the rulult folder
if [ ! -d "$path_results" ]; then
    mkdir -p "$path_results"
fi

# ---- Check if result_pref is empty, and if so, assign the default value
if [ -z "$result_pref" ]; then
    result_pref="${path_results}blast_"
fi
# Ensure result_pref ends with an underscore if provided
result_pref=$(add_symbol_if_missing "$result_pref" "_")



# -------------------------------------------------------------
#     Create an array to store blast databases
# -------------------------------------------------------------

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
        echo "Creating database: ${db_path}"
        makeblastdb -in "$fasta_file" -dbtype nucl -out "$db_path" > /dev/null
    fi
    blast_databases+=("$db_path")
done

# echo "${blast_databases[@]}"


# -------------------------------------------------------------
#     Perform BLAST on all databases
# -------------------------------------------------------------
query_file_base=$(basename "$query_file" .fasta)
export result_pref
export query_file_base
export query_file

# Create a function that you will call in parallel tasks
process_db() {
    db_path="$1"
    db_name=$(basename "$db_path")
    # echo "${result_pref}${query_file_base}_${db_name}.txt"

    blastn -query "$query_file" -db "$db_path" -out "${result_pref}${query_file_base}_${db_name}.txt" -outfmt "6 qseqid qstart qend sstart send pident length sseqid sseq" 2>/dev/null 

    # blastn -query "$query_file" -db "$db_path" -out "${result_pref}${query_file_base}_${db_name}.txt"
}

# Export the function to make it available inside parallel tasks
export -f process_db

# Run parallel execution of the loop with the specified number of cores
parallel -j "$n_cores" process_db ::: "${blast_databases[@]}"


# -------------------------------------------------------------
#     Merging the results
# -------------------------------------------------------------

file_array=("${result_pref}${query_file_base}"*.txt)
for i in "${!file_array[@]}"; do
   if [[ ${file_array[i]} == "$merged_file" ]]; then
      unset 'file_array[i]'
   fi
done

# echo "${file_array[@]}"
# echo "$merged_file"

cat "${file_array[@]}" > "$merged_file"
rm "${file_array[@]}"
