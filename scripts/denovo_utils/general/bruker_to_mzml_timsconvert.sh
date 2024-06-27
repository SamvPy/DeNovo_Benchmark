#!/bin/bash
# Check if filename is provided as the first argument
if [ $# -eq 0 ]; then
    echo "Usage: $0 <filename>"
    exit 1
fi

# Read the filename from the first argument
filename=$1

# Check if the file exists
if [ ! -f "$filename" ]; then
    echo "File '$filename' not found."
    exit 1
fi

# Process each line in the file
while IFS= read -r line; do
    # Trim leading and trailing whitespace from the line
    line=$(echo "$line" | sed -e 's/^[[:space:]]*//' -e 's/[[:space:]]*$//')
    
    # Check if the line is not empty
    if [ -n "$line" ]; then
        # Run timsconvert command for each mzml file
        timsconvert --input "$line" --outdir /public/compomics3/Sam/CMB1473/mzml/timsconvert --compression none --exclude_mobility --verbose
    fi
done < "$filename"