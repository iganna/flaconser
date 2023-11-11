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
        -m|--merged-file) merged_file="$2"; shift 2;;
		-d|--depth) n_depth="$2"; shift 2;;
        -n|--n-cores) n_cores="$2"; shift 2;;
        *)
            echo "Unknown option: $1"
            exit 1;;
    esac
done

if [ -z "$n_depth" ]; then
    n_depth=1
fi

# Number of cores
if ! [[ "$n_cores" =~ ^[0-9]+$ ]]; then
    n_cores=1
fi

# Check if required parameters are missing and produce an error if they are
if [ -z "$path_results" ] || [ -z "$path_genomes" ] || [ -z "$fasta_type" ] || [ -z "$query_file" ]; then
    echo "Error: Missing required parameter(s)"
    exit 1
fi

# Add trailing slashes to path variables if missing
path_results=$(add_symbol_if_missing "$path_results" "/")
path_genomes=$(add_symbol_if_missing "$path_genomes" "/")

echo ${path_results}
echo ${path_genomes}


# Create the rulult folder
if [ ! -d "$path_results" ]; then
    mkdir -p "$path_results"
fi

# Ensure result_pref ends with an underscore if provided
result_pref=$(add_symbol_if_missing "$result_pref" "_")


query_file_new=${path_results}new_query.fasta
cp ${query_file} ${query_file_new}


for i in $(seq 1 $n_depth)
do
	./one_search.sh -r ${result_pref} -g ${path_genomes} -t ${fasta_type} -q ${query_file_new} -n 30 -m ${merged_file}
	Rscript your_script.R --merged_file ${merged_file} --target_file ${query_file_new} --seq_len 140 --seq_cover 0.9

	cat ${query_file} ${query_file_new} > ${query_file_new}
done



rm ${query_file_new}
