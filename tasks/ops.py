#!/usr/bin/env python3
'''
------------------------------------------------
 Final OPS (extract, filter & covert)
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
        seqs = list(io.read_fasta(input_path))
        chosen = random.sample(seqs, k=min(k, len(seqs)))
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
                if len(s) >= min_len:
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
    return kept
    
#---------------Converting------------------------------
def convert_fastq_to_fasta(input_path: Path, out_path: Path) -> int:
    count = 0
    def it():
        nonlocal count
        for h, s, _q in io.read_fastq(input_path):
            count += 1
            yield h, s
    io.write_fasta(it(), out_path)
    return count