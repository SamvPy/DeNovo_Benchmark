#!/bin/bash

# Check if correct number of arguments are provided
if [ "$#" -ne 2 ]; then
    echo "Usage: $0 <file_containing_paths> <output_directory>"
    exit 1
fi

# Assigning arguments to variables
file_list=$1
output_dir=$2

# Iterate over each path in the file
while IFS= read -r file_path; do
    # Extracting filename without extension
    filename=$(basename "$file_path")
    filename_no_ext="${filename%.*}"

    # Running casanovo sequence command
    python -m ContraNovo.ContraNovo --mode=denovo --peak_path="${file_path}" --output="${output_dir}/${filename}" --model=/home/sam/denovo_project/model_checkpoints/ContraNovo/ControNovo.ckpt
done < "$file_list"
