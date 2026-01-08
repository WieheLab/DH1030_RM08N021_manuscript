#!/bin/bash

# Check if a file path is provided
if [ -z "$1" ]; then
    echo "Usage: $0 <path_to_csv_file>"
    exit 1
fi

# Get the first line (header) of the CSV file
header=$(head -n 1 "$1")

# Split the header into columns
IFS=',' read -r -a columns <<< "$header"

# Print column names
for column in "${columns[@]}"
do
    echo "$column"
done

