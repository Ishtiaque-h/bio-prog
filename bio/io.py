#!/usr/bin/env python3
'''
------------------------------------------------
 I/O files (FASTA/ FASTQ readers & writers)
------------------------------------------------
'''
<<<<<<< bio/io.py

=======
import os
from pathlib import Path

def _open(filename, mode: str):
    """Open a sequence file."""
    if not os.path.exists(filename):
        print('Error - Invalid filename.')
        return None
    return open(filename, mode)

def detect_extension(filename):
    # Check file extension to know if it's FASTA or FASTQ
    if filename.endswith((".fasta", ".fa",".fsa", ".fna", ".seq",".pep"):
        return "FASTA"
    elif filename.endswith(".fastq",".fq"):
        return "FASTQ"
    else:
        return None


