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

