#!/bin/bash

# ----------------------------------------------------------------------------
#            ERROR HANDLING BLOCK
# ----------------------------------------------------------------------------

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


# ----------------------------------------------------------------------------
#             FUNCTIONS
# ----------------------------------------------------------------------------

source func_search.sh

show_help() {
    cat << EOF
Usage: $(basename $0) [OPTIONS]

DESCRIPTION:
    This script is designed for bioinformatics data analysis, particularly for the analysis of genomic sequences. It includes error handling, parameter parsing, and execution of sequence analysis tools.

OPTIONS:
    -r, --path-results      Specify the path where results will be saved.
    -p, --result-prefix     Prefix for the result files.
    -q, --file-query        Path to the query file.
    -g, --path-genomes      Define the path to the genome files.
    -t, --file-type         Set the type of fasta files to be used.
    -d, --depth             Number of search cycles to perform.
    -s, --sim-cover         Similarity coverage parameter for analysis.
    -n, --n-cores           Number of cores to use for the analysis.
    -a, --aln-mafft         Flag to perform alignment using MAFFT for the gound hits.
    --pal                   Flag to search for palindromes.
    --ltr                   Flag to search for long terminal repeats (LTR).

EXAMPLES:
    $(basename $0) -r /path/to/results -g /path/to/genomes -q /path/to/query.fasta -s 0.85 -n 4

EOF
}

# ----------------------------------------------------------------------------
#                MAIN
# ----------------------------------------------------------------------------

# ---- Parse command-line parameters
while [[ $# -gt 0 ]]; do
    case "$1" in
        -h|--help) show_help; exit 0;;
        -r|--path-results) path_results="$2"; shift 2;;
        -g|--path-genomes) path_genomes="$2"; shift 2;;
        -t|--file-type) fasta_type="$2"; shift 2;;
        -p|--result-prefix) pref_result="$2"; shift 2;;
        -q|--file-query) file_query="$2"; shift 2;;
		-d|--depth) n_depth="$2"; shift 2;;
        -s|--sim-cover) sim_cover="$2"; shift 2;;
        -n|--n-cores) n_cores="$2"; shift 2;;
        -a|--aln-mafft) aln_mafft=true; shift;;
        --pal) search_palindromes=true; shift;;
        --ltr) search_ltr=true; shift;;
        *)
            echo "Unknown option: $1"
            exit 1;;
    esac
done

# ---- Number of search cycles
if [ -z "$n_depth" ]; then
    n_depth=1
fi

# ---- Number of cores
if ! [[ "$n_cores" =~ ^[0-9]+$ ]]; then
    n_cores=1
fi


# ---- Check of missimg parameters

check_missing_variable "path_results"
check_missing_variable "path_genomes"
check_missing_variable "fasta_type"
check_missing_variable "file_query"
check_missing_variable "sim_cover"


# ---- Add trailing slashes to path variables if missing

path_results=$(add_symbol_if_missing "$path_results" "/")
path_genomes=$(add_symbol_if_missing "$path_genomes" "/")


# ---- Create the rulult folder
if [ ! -d "$path_results" ]; then
    mkdir -p "$path_results"
fi


# ---- Create files in the the rulult folder
file_query_new=${path_results}new_query.fasta  # new query file to work with
file_merged=${path_results}merged.txt # file with merged blast results
file_out=${path_results}out.rds  # file with table of best blast hits


# Calling the result directory for each file
remove_file_if_exists "${file_query_new}"
remove_file_if_exists "${file_out}"
remove_file_if_exists "${file_merged}"

# Copy query file to work with it locally
cp ${file_query} ${file_query_new}

# ---- Looping for a specified number of search rounds
for i in $(seq 1 $n_depth)
do
    echo "Round ${i}"
    command="./one_search.sh"

    # Adding the -r parameter if ${pref_result} is set
    if [ -n "${pref_result}" ]; then
        command+=" -p ${pref_result}"
    fi

    # Adding other parameters to the command
    command+=" -g ${path_genomes} \
               -r ${path_results} \
               -t ${fasta_type} \
               -q ${file_query_new} \
               -m ${file_merged} \
               -n 30"

    
    # Executing the constructed command
    # echo "$command"
    eval "$command"


    # Running analysis of results
    final_status=$(Rscript one_preparation.R -q ${file_query_new} \
                              -m ${file_merged} \
                              -o ${file_out} \
                              -s ${sim_cover} )

    # echo ${final_status}

    if echo "$final_status" | grep -q 'Final'; then
        echo "DONE: The search have converged"
        break
    fi

done


# ---- Remove intermediate files

remove_file_if_exists "${file_query_new}"
remove_file_if_exists "${file_merged}"

# ---- save files into the file, in correct orientation
file_seqs=${path_results}seqs.fasta
remove_file_if_exists "${file_seqs}"

Rscript get_sequences.R -i ${file_out} -o ${file_seqs}


# ----------------------------------------------------------------------------
#                ADDITIONAL ANALYSIS
# ----------------------------------------------------------------------------

# ---- Run MAFFT

if [ -n "${aln_mafft}" ]; then
    file_aln=${path_results}aln_mafft.fasta
    remove_file_if_exists "${file_aln}"

    mafft --quiet  ${file_seqs} > ${file_aln} 
fi


# ---- Run Palindrom search

min_len=10000

if [ -n "${search_palindromes}" ]; then
    file_pal_pair=${path_results}palind_pair.rds 
    remove_file_if_exists "${file_pal_pair}"

    file_pal_single=${path_results}palind_single.rds 
    remove_file_if_exists "${file_pal_single}"

    Rscript find_palindromes.R -i ${file_out} -o ${file_pal_pair} -r ${file_pal_single} --min_len ${min_len}
fi


# ---- Run LTR search

if [ -n "${search_ltr}" ]; then
    file_ltr_pair=${path_results}ltr_pair.rds 
    remove_file_if_exists "${file_ltr_pair}"

    file_ltr_single=${path_results}ltr_single.rds 
    remove_file_if_exists "${file_ltr_single}"

    Rscript find_ltr.R -i ${file_out} -o ${file_ltr_pair} -r ${file_ltr_single} --min_len ${min_len}
fi
