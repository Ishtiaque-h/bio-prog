#!/usr/bin/env python3
'''
------------------------------------------------
 I/O files (FASTA/ FASTQ readers & writers)
------------------------------------------------
'''
def detect_format(filename):
    # Check file extension to know if it's FASTA or FASTQ
    if filename.endswith((".fasta", ".fa", ".fna", ".txt")):
        return "FASTA"
    elif filename.endswith(".fastq"):
        return "FASTQ"
    else:
        return None



