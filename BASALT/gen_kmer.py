#!/usr/bin/env python
# -*- coding: UTF-8 -*-

"""
gen_kmer.py - Generate k-mer frequency profiles for MetaBinner

This script is a wrapper around calc.kmerfreq.pl that generates k-mer
frequency profiles in the format expected by MetaBinner.

Usage:
    gen_kmer.py <assembly_file> <min_length> <kmer_length>

Arguments:
    assembly_file : Path to the assembly FASTA file
    min_length    : Minimum contig length (default: 1000)
    kmer_length   : K-mer length (default: 4)
"""

import sys
import os
import subprocess

def main():
    if len(sys.argv) < 4:
        print("Usage: gen_kmer.py <assembly_file> <min_length> <kmer_length>")
        sys.exit(1)

    assembly_file = sys.argv[1]
    min_length = sys.argv[2]
    kmer_length = sys.argv[3]

    # Get the directory where this script is located
    script_dir = os.path.dirname(os.path.abspath(__file__))

    # Path to calc.kmerfreq.pl
    calc_kmer_script = os.path.join(script_dir, 'calc.kmerfreq.pl')

    # Check if calc.kmerfreq.pl exists
    if not os.path.exists(calc_kmer_script):
        print(f"Error: calc.kmerfreq.pl not found at {calc_kmer_script}")
        sys.exit(1)

    # Extract assembly name without extension
    assembly_name_parts = os.path.basename(assembly_file).split('.')
    if len(assembly_name_parts) > 1:
        assembly_name_parts.pop()  # Remove extension
    assembly_name = '.'.join(assembly_name_parts)

    # Output file name
    output_file = f"{assembly_name}_kmer_{kmer_length}_f{min_length}.csv"

    # Build the command
    cmd = [
        'perl',
        calc_kmer_script,
        '-i', assembly_file,
        '-o', output_file,
        '-m', min_length,
        '-k', kmer_length
    ]

    # Run the command
    try:
        result = subprocess.run(cmd, check=True, capture_output=True, text=True)
        if result.stdout:
            print(result.stdout)
        if result.stderr:
            print(result.stderr, file=sys.stderr)
        print(f"K-mer profile generated: {output_file}")
    except subprocess.CalledProcessError as e:
        print(f"Error running calc.kmerfreq.pl: {e}")
        if e.stdout:
            print(e.stdout)
        if e.stderr:
            print(e.stderr, file=sys.stderr)
        sys.exit(1)
    except FileNotFoundError:
        print("Error: perl not found. Please ensure perl is installed and in your PATH.")
        sys.exit(1)

if __name__ == '__main__':
    main()
