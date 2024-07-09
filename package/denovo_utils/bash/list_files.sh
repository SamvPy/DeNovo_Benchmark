#!/bin/bash

# Check if required argument is provided
if [ "$#" -ne 1 ]; then
    echo "Usage: $0 <directory_path>"
    exit 1
fi

# Assign input argument to variable
directory="$1"

# Check if directory exists
if [ ! -d "$directory" ]; then
    echo "Error: Directory '$directory' does not exist."
    exit 1
fi

# List all files in the directory recursively and save paths to a text file
find "$directory" -mindepth 1 -maxdepth 1 -type f > /home/sam/denovo_project/util_files/file_paths.txt

echo "File paths listed in file_paths.txt"
