#!/usr/bin/env python3
import pandas as pd
import os

# List of CSV files and their column mappings
csv_files = {
    "lookup_times.csv": {"MemCache": "Baseline"},
    "range_scan_times.csv": {"MemCache": "Baseline"},
    "mixed_workload_times.csv": {"MemCache": "Baseline"},
    "hit_ratios.csv": {"MemCacheHitRatio": "BaselineHitRatio"},
    "eviction_stats.csv": {"MemCacheEvictions": "BaselineEvictions"},
    "summary_stats.csv": {"MemCacheMedian": "BaselineMedian"}
}

def rename_csv_columns():
    for file_name, column_mapping in csv_files.items():
        if os.path.exists(file_name):
            print(f"Processing {file_name}...")
            try:
                # Read the CSV
                df = pd.read_csv(file_name)
                
                # Rename columns
                df = df.rename(columns=column_mapping)
                
                # Save the updated CSV
                df.to_csv(file_name, index=False)
                print(f"Updated {file_name} with new column names: {column_mapping}")
            except Exception as e:
                print(f"Error processing {file_name}: {e}")
        else:
            print(f"Warning: {file_name} not found, skipping...")

if __name__ == "__main__":
    print("Renaming columns in CSV files...")
    rename_csv_columns()
    print("Done!")