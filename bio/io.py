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
=======
def read_fastq(file):
    f = _open(file, "r")
    lines = []
    for line in f:
        line = line.strip()
        if line != "":
            lines.append(line)

    if len(lines) == 0:
        print("Error: The file is empty.")
        return None

    if len(lines) % 4 != 0:
        print("Error: FASTQ files must have groups of 4 lines per record.")
        return None

    seqs = []
    i = 0
    while i < len(lines):
        header = lines[i]
        seq = lines[i + 1]
        plus = lines[i + 2]
        qual = lines[i + 3]

        if not header.startswith("@"):
            print("Error: FASTQ header must start with '@'.")
            return None
        if not plus.startswith("+"):
            print("Error: Third line must start with '+'.")
            return None
        if len(seq) != len(qual):
            print("Error: Sequence and quality must be same length.")
            return None

        seqs.append((header[1:], seq, qual))
        i = i + 4

    return seqs