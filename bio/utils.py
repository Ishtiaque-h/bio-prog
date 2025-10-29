#!/usr/bin/env python3
'''
------------------------------------------------
 Helper/ utility functions
------------------------------------------------
'''
from .io import _open
from typing import Tuple

#----------------extension checker---------------------
def detect_extension(filename):
    # Check file extension to know if it's FASTA or FASTQ
    if filename.endswith(".fasta", ".fa",".fsa", ".fna", ".seq",".pep"):
        return "FASTA"
    elif filename.endswith(".fastq",".fq"):
        return "FASTQ"
    else:
        return None

#----------------fasta/fastq format checker---------------------    
def sniff_format(path_or_file) -> str:
    """
    Return 'fasta' or 'fastq' by peeking first non-empty character.
    Raises ValueError if format cannot be determined.
    """
    fh = _open(path_or_file, "r")
    pos = fh.tell()
    for line in fh:
        s = line.strip()
        if not s:
            continue
        if s.startswith(">"):
            fh.seek(pos)  # reset if it's a real file handle
            return "fasta"
        if s.startswith("@"):
            fh.seek(pos)
            return "fastq"
        break
    fh.seek(pos)
    raise ValueError("Unable to determine file format (expected FASTA or FASTQ).")

#----------------gc and atgc estimator-------------------------
def gc_and_counts(seq: str) -> Tuple[int, int, int]:
    """
    Return (gc_count, atgc_count, n_count) for a sequence.
    Ns are counted separately; atgc_count excludes Ns.
    """
    s = seq.upper()
    g = s.count("G"); c = s.count("C")
    a = s.count("A"); t = s.count("T")
    n = s.count("N")
    gc = g + c
    atgc_total = a + t + g + c
    return gc, atgc_total, n
