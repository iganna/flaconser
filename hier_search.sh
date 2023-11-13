# ==============================================================================

# exit when any command fails
set -e

# keep track of the last executed command
trap 'last_command=$current_command; current_command=$BASH_COMMAND' DEBUG
# echo an error message before exiting
#trap 'echo "\"${last_command}\" command filed with exit code $?."' EXIT

trap 'catch $?' EXIT
catch() {
if [ $1 -ne 0 ]; then
        echo "\"${last_command}\" command filed with exit code $1."
fi
#  echo "\"${last_command}\" command filed with exit code $1."
}
# ==============================================================================


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
        -p|--result-prefix) pref_result="$2"; shift 2;;
        -q|--file-query) file_query="$2"; shift 2;;
        -m|--file-merged) file_merged="$2"; shift 2;;
        -o|--file-out) file_out="$2"; shift 2;;
        -f|--file-final) file_final="$2"; shift 2;;
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

declare -A params=(
    ["path_results"]="$path_results"
    ["path_genomes"]="$path_genomes"
    ["fasta_type"]="$fasta_type"
    ["file_query"]="$file_query"
    ["file_out"]="$file_out"
    ["file_final"]="$file_final"
)

# Цикл для проверки каждого параметра
missing_params=false
for param in "${!params[@]}"; do
    if [ -z "${params[$param]}" ]; then
        echo "Error: Missing required parameter - $param"
        missing_params=true
    fi
done

# Выход из скрипта, если есть незаданные параметры
if [ "$missing_params" = true ]; then
    exit 1
fi

# Check if required parameters are missing and produce an error if they are
if [ -z "$path_results" ] || [ -z "$path_genomes" ] || [ -z "$fasta_type" ] ||  \
    [ -z "$file_query" ] || [ -z "$file_out" ] || [ -z "$file_final" ]; then
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

# Ensure pref_result ends with an underscore if provided
pref_result=$(add_symbol_if_missing "$pref_result" "_")


file_query_new=${path_results}new_query.fasta
cp ${file_query} ${file_query_new}

# Crean the directory before the work
file_out=${path_results}out.rds
file_merged=${path_results}merged.fasta
file_final=${path_results}final.fasta

[ -f "${file_out}" ] && rm "${file_out}"
[ -f "${file_merged}" ] && rm "${file_merged}"
[ -f "${file_final}" ] && rm "${file_final}"

for i in $(seq 1 $n_depth)
do


    echo "./one_search.sh -r ${pref_result} \
                    -g ${path_genomes} \
                    -t ${fasta_type} \
                    -q ${file_query_new} \
                    -m ${file_merged} \
                    -n 30"


	./one_search.sh -r ${pref_result} \
                    -g ${path_genomes} \
                    -t ${fasta_type} \
                    -q ${file_query_new} \
                    -m ${file_merged} \
                    -n 30 

	Rscript one_preparation.R -q ${file_query_new} \
                              -m ${file_merged} \
                              -o ${file_out} \
                              -f ${file_final} \
                              -s 0.9


	if [ -f ${file_final} ]; then
        break  
    fi
done



rm ${file_query_new}
