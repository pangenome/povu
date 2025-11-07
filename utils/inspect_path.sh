#!/bin/bash

# Check if required arguments are provided
if [ "$#" -ne 4 ]; then
  echo "Usage: $0 <input_file> <left> <middle> <right>"
  echo "Example: $0 lpa.refs 13 186 10"
  exit 1
fi

# Assign command-line arguments to variables
input_file="$1"
left="$2"
middle="$3"
right="$4"

# Construct grep pattern
pattern=".{0,$left},$middle[+-].{0,$right}"

# Process the file
awk '{print NR-1 "\t" $0}' "$input_file" | while IFS=$'\t' read -r zero_index _ second_col rest; do
  grep -oP "$pattern" <<< "$rest" | awk -v col="$second_col" -v idx="$zero_index" '{print idx "\t" col "\t" $0}'
done
