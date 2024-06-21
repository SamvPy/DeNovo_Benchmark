#!/bin/bash

# Check if required arguments are provided
if [ "$#" -ne 2 ]; then
    echo "Usage: $0 <target_folder> <input_file>"
    exit 1
fi

# Assign input arguments to variables
target_folder="$1"
input_file="$2"

# Check if target folder exists
if [ ! -d "$target_folder" ]; then
    echo "Error: Target folder '$target_folder' does not exist."
    exit 1
fi

# Check if input file exists
if [ ! -f "$input_file" ]; then
    echo "Error: Input file '$input_file' does not exist."
    exit 1
fi

# Loop through each line in the input file
while IFS= read -r line; do
    # Extract filename from path
    filename=$(basename "$line")
    # Remove extension from filename
    filename_no_ext="${filename%.*}"
    # Construct target file path
    target_file="$target_folder/$filename_no_ext.ipc"
    
    # Run the conversion command
    python -m instanovo.utils.convert_to_ipc "$line" "$target_file" \
        --source_type mzml \
        --verbose
    
    echo "Converted $line to $target_file"
done < "$input_file"