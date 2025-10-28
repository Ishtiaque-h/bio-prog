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
=======
def read_fasta(file):
 f= _open(file, "r")   
    lines = []
    for line in f:
        line = line.strip()
        if line != "":
            lines.append(line)

    if len(lines) == 0:
        print("Error: The file is empty.")
        return None

    if not lines[0].startswith(">"):
        print("Error: FASTA file must start with '>'.")
        return None

    seqs = []
    header = ""
    seq = ""

    for line in lines:
        if line.startswith(">"):
            if header != "":
                seqs.append((header, seq))
            header = line[1:]
            seq = ""
        else:
            seq = seq + line

    if header != "" and seq != "":
        seqs.append((header, seq))

    return seqs

