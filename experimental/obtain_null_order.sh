#!/bin/bash

FILE_DIR="workflow/results_shuffled/"
BIN="scripts_opt/bin/null_order"          # adjust relative path as needed

shopt -s extglob

count=0

echo "organism,k,order,origin,counter" >> total_order.csv

for file in "$FILE_DIR"k+([0-9])*/*/null_bit_format
do
    [[ -f "$file" ]] || continue
    DIR_PATH=$(dirname "$file")
    organism="${DIR_PATH##*/}"
    if [[ "$file" =~ k([0-9]+) ]]; then
        k="${BASH_REMATCH[1]}"
    else
        k="null"
    fi
    echo "$organism"
    echo "Processing: $file"
    $BIN "$file" 4
    ((count++))
    [[ -d order_results ]] || mkdir order_results
    for order_file in *.txt
    do
        order="${order_file: (-5):1}"
        [[ -f "$order_file" ]] || continue
        echo -n "$organism","$k","$order",0, >> total_order.csv
        wc -l < "$order_file" >> total_order.csv
        #[[ -d order_results/"$k"/"$organism" ]] || mkdir -p order_results/"$k"/"$organism"
        #mv "$order_file" order_results/"$k"/"$organism"/"$order_file"
        rm "$order_file"
    done
done

echo $count
echo "All files processed."
