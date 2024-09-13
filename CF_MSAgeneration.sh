#!/bin/bash

# Load necessary modules
module load colabfold/1.5.2/gcc

# Arguments
fasta_file="$1"
output_dir="$2"

# Check if arguments are provided
if [ -z "$fasta_file" ] || [ -z "$output_dir" ]; then
  echo "Usage: $0 <fasta_file> <output_dir>"
  exit 1
fi

# Ensure the output directory exists
mkdir -p "$output_dir"

# Set environment variables if needed
export RESULT_FOLDER="$output_dir"

# Run colabfold_search
colabfold_search --threads=32 --db-load-mode=2 "$fasta_file" "${COLABFOLD_DB}" "$RESULT_FOLDER"