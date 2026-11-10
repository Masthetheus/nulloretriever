#!/bin/bash

FILE_DIR="$1"
BIN="scripts_opt/bin/null_order"          # adjust relative path as needed

# Identificador do FILE_DIR: assume terminar em _<numero>
FILE_DIR="${FILE_DIR%/}"                  # remove barra final, se houver
FILE_DIR_ID="${FILE_DIR##*_}"             # ex.: results_shuffled_42 -> 42

OUT_CSV="total_order_${FILE_DIR_ID}.csv"

shopt -s extglob

count=0

echo "organism,k,order,origin,counter" >> "$OUT_CSV"

for file in "$FILE_DIR"/*/k+([0-9])/null_bit_format
do
    [[ -f "$file" ]] || continue
    DIR_PATH=$(dirname "$file")           # .../organism/k<N>
    k="${DIR_PATH##*/}"                   # k<N>
    k="${k#k}"                            # <N>
    organism=$(basename "$(dirname "$DIR_PATH")")

    echo "$organism"
    echo "Processing: $file"
    $BIN "$file" 4
    ((count++))

    [[ -d order_results ]] || mkdir order_results

    for order_file in *.txt
    do
        order="${order_file: (-5):1}"
        [[ -f "$order_file" ]] || continue
        echo -n "$organism","$k","$order",0, >> "$OUT_CSV"
        wc -l < "$order_file" >> "$OUT_CSV"
        #[[ -d order_results/"$k"/"$organism" ]] || mkdir -p order_results/"$k"/"$organism"
        #mv "$order_file" order_results/"$k"/"$organism"/"$order_file"
        rm "$order_file"
    done
done

echo $count
echo "All files processed."
