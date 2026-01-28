#!/usr/bin/env python3
"""
getTables.py - Extract and tabulate benchmark results from the Results folder.

Extracts:
- preproc time (party 2)
- preproc sent (from party 0)
- online time (party 2)
- online sent (party 1 + party 2)

Outputs tables for varying num-verts, num-edges, and num-clients.
Times are in seconds, communication in MB.
"""

import os
import re
from pathlib import Path
from collections import defaultdict
import csv


def parse_log_file(filepath):
    """Parse a party log file and extract metrics."""
    metrics = {}
    try:
        with open(filepath, 'r') as f:
            content = f.read()
            
            # Extract preproc time (in ms)
            match = re.search(r'preproc time:\s+([\d.]+)\s+ms', content)
            if match:
                metrics['preproc_time_ms'] = float(match.group(1))
            
            # Extract preproc sent (in bytes)
            match = re.search(r'preproc sent:\s+([\d.]+)\s+bytes', content)
            if match:
                metrics['preproc_sent_bytes'] = float(match.group(1))
            
            # Extract online time (in ms)
            match = re.search(r'online time:\s+([\d.]+)\s+ms', content)
            if match:
                metrics['online_time_ms'] = float(match.group(1))
            
            # Extract online sent (in bytes)
            match = re.search(r'online sent:\s+([\d.]+)\s+bytes', content)
            if match:
                metrics['online_sent_bytes'] = float(match.group(1))
                
    except FileNotFoundError:
        pass
    
    return metrics


def extract_benchmark_data(results_dir):
    """
    Extract data from all benchmarks in the Results folder.
    
    Returns a nested dictionary:
    {
        'benchmark_name': {
            (num_verts, num_edges, num_clients): {
                'preproc_time_sec': float,
                'preproc_sent_mb': float,
                'online_time_sec': float,
                'online_sent_mb': float
            }
        }
    }
    """
    results_path = Path(results_dir)
    data = defaultdict(dict)
    
    # Iterate through all benchmark folders
    for benchmark_dir in results_path.iterdir():
        if not benchmark_dir.is_dir():
            continue
            
        benchmark_name = benchmark_dir.name
        print(f"Processing benchmark: {benchmark_name}")
        
        # Iterate through num_verts folders
        for verts_dir in benchmark_dir.iterdir():
            if not verts_dir.is_dir():
                continue
            
            num_verts = int(verts_dir.name)
            
            # Iterate through num_edges folders
            for edges_dir in verts_dir.iterdir():
                if not edges_dir.is_dir():
                    continue
                
                num_edges = int(edges_dir.name)
                
                # Iterate through num_clients folders
                for clients_dir in edges_dir.iterdir():
                    if not clients_dir.is_dir():
                        continue
                    
                    num_clients = int(clients_dir.name)
                    
                    # Parse party log files
                    party_0_log = clients_dir / 'party_0.log'
                    party_1_log = clients_dir / 'party_1.log'
                    party_2_log = clients_dir / 'party_2.log'
                    
                    party_0_data = parse_log_file(party_0_log)
                    party_1_data = parse_log_file(party_1_log)
                    party_2_data = parse_log_file(party_2_log)
                    
                    # Extract required metrics
                    if party_2_data and party_0_data:
                        preproc_time_sec = party_2_data.get('preproc_time_ms', 0) / 1000.0
                        preproc_sent_mb = party_0_data.get('preproc_sent_bytes', 0) / (1024 * 1024)
                        online_time_sec = party_2_data.get('online_time_ms', 0) / 1000.0
                        online_sent_mb = (
                            party_1_data.get('online_sent_bytes', 0) + 
                            party_2_data.get('online_sent_bytes', 0)
                        ) / (1024 * 1024)
                        
                        key = (num_verts, num_edges, num_clients)
                        data[benchmark_name][key] = {
                            'preproc_time_sec': preproc_time_sec,
                            'preproc_sent_mb': preproc_sent_mb,
                            'online_time_sec': online_time_sec,
                            'online_sent_mb': online_sent_mb
                        }
    
    return data


def print_table(benchmark_name, benchmark_data):
    """Print a formatted table for a single benchmark."""
    print(f"\n{'='*100}")
    print(f"Benchmark: {benchmark_name}")
    print(f"{'='*100}")
    
    # Sort by num_verts, num_edges, num_clients
    sorted_keys = sorted(benchmark_data.keys())
    
    # Print header
    header = f"{'Num Verts':<12} {'Num Edges':<12} {'Num Clients':<12} {'Preproc Time (s)':<18} {'Preproc Sent (MB)':<18} {'Online Time (s)':<18} {'Online Sent (MB)':<18}"
    print(header)
    print('-' * len(header))
    
    # Print data rows
    for key in sorted_keys:
        num_verts, num_edges, num_clients = key
        metrics = benchmark_data[key]
        
        row = (f"{num_verts:<12} {num_edges:<12} {num_clients:<12} "
               f"{metrics['preproc_time_sec']:<18.6f} {metrics['preproc_sent_mb']:<18.6f} "
               f"{metrics['online_time_sec']:<18.6f} {metrics['online_sent_mb']:<18.6f}")
        print(row)


def save_to_csv(benchmark_name, benchmark_data, output_dir='./'):
    """Save benchmark data to a CSV file."""
    output_path = Path(output_dir)
    output_path.mkdir(exist_ok=True)
    
    csv_filename = output_path / f"{benchmark_name}_results.csv"
    
    with open(csv_filename, 'w', newline='') as csvfile:
        fieldnames = ['Num Verts', 'Num Edges', 'Num Clients', 
                     'Preproc Time (s)', 'Preproc Sent (MB)', 
                     'Online Time (s)', 'Online Sent (MB)']
        writer = csv.DictWriter(csvfile, fieldnames=fieldnames)
        
        writer.writeheader()
        
        sorted_keys = sorted(benchmark_data.keys())
        for key in sorted_keys:
            num_verts, num_edges, num_clients = key
            metrics = benchmark_data[key]
            
            writer.writerow({
                'Num Verts': num_verts,
                'Num Edges': num_edges,
                'Num Clients': num_clients,
                'Preproc Time (s)': f"{metrics['preproc_time_sec']:.6f}",
                'Preproc Sent (MB)': f"{metrics['preproc_sent_mb']:.6f}",
                'Online Time (s)': f"{metrics['online_time_sec']:.6f}",
                'Online Sent (MB)': f"{metrics['online_sent_mb']:.6f}"
            })
    
    print(f"  Saved to: {csv_filename}")


def main():
    # Get the script directory
    script_dir = Path(__file__).parent
    results_dir = script_dir / 'Results'
    
    if not results_dir.exists():
        print(f"Error: Results directory not found at {results_dir}")
        return
    
    print(f"Extracting data from: {results_dir}")
    print("="*100)
    
    # Extract all benchmark data
    all_data = extract_benchmark_data(results_dir)
    
    if not all_data:
        print("No data found in Results directory!")
        return
    
    # Print tables for each benchmark
    for benchmark_name in sorted(all_data.keys()):
        print_table(benchmark_name, all_data[benchmark_name])
    
    # Save to CSV files
    print(f"\n{'='*100}")
    print("Saving results to CSV files...")
    print(f"{'='*100}")
    
    csv_output_dir = script_dir / 'results_tables'
    for benchmark_name in sorted(all_data.keys()):
        save_to_csv(benchmark_name, all_data[benchmark_name], csv_output_dir)
    
    print(f"\nAll tables saved to: {csv_output_dir}")


if __name__ == '__main__':
    main()
