#!/bin/bash

# Specify the folder path
folder_path="$1"

# Check if the folder exists
if [ ! -d "$folder_path" ]; then
  echo "The folder $folder_path does not exist."
  exit 1
fi

# Loop through each file in the folder
for file in "$folder_path"/*; do
  # Check if it's a file (skip directories)
  if [ -f "$file" ]; then
    echo "Processing file: $file"
    # Execute the Python script for the current file
    python3 gsm_automated_plots.py --json "$file"
  fi
done
