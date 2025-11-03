#!/usr/bin/env python3
'''
------------------------------------------------
 Final OPS (extract, filter & covert)
 This module implements high-level sequence operations:
   • Extract: randomly sample sequences from a FASTA/FASTQ file.
   • Filter: retain only sequences longer than a user-specified threshold.
   • Convert: transform FASTQ records into FASTA format.
------------------------------------------------
'''

#---------------imports------------------------------
import random
from pathlib import Path
from typing import List, Tuple
from bio import io

#---------------Extracting------------------------------
def extract_random(input_path: Path, fmt: str, out_path: Path, k: int = 25) -> None:
    if fmt == "fasta":
        seqs = list(io.read_fasta(input_path))  # load all records into memory
        chosen = random.sample(seqs, k=min(k, len(seqs)))   # randomly choose up to k
        io.write_fasta(chosen, out_path)
    else:
        seqs = list(io.read_fastq(input_path))
        chosen = random.sample(seqs, k=min(k, len(seqs)))
        io.write_fastq(chosen, out_path)

#---------------Filtering------------------------------
def filter_by_min_len(input_path: Path, fmt: str, out_path: Path, min_len: int) -> int:
    kept = 0
    if fmt == "fasta":
        def it():
            nonlocal kept
            for h, s in io.read_fasta(input_path):
                if len(s) >= min_len:   # skip empty or short sequences
                    kept += 1
                    yield h, s
        io.write_fasta(it(), out_path)
    else:
        def it():
            nonlocal kept
            for h, s, q in io.read_fastq(input_path):
                if len(s) >= min_len:
                    kept += 1
                    yield h, s, q
        io.write_fastq(it(), out_path)
    return kept # number of sequences retained
    
#---------------Converting------------------------------
'''
------------------------------------------------
Behavior:
        - Drops quality lines ('+' and quality string).
        - Replaces '@' headers with '>'.
        - Writes sequences in standard FASTA format.
        - Raises ValueError if input file is malformed.
------------------------------------------------ '''
def convert_fastq_to_fasta(input_path: Path, out_path: Path) -> int:
    count = 0
    def it():
        nonlocal count
        for h, s, _q in io.read_fastq(input_path):
            count += 1
            yield h, s
    io.write_fasta(it(), out_path)
    return count    # number of sequences successfully converted