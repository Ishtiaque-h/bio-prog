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

=====================

def _open(filename, mode: str):
    """Open a sequence file."""
    if not os.path.exists(filename):
        print('Error - Invalid filename.')
        return None
    return open(filename, mode)

======

def detect_extension(filename):
    # Check file extension to know if it's FASTA or FASTQ
    if filename.endswith((".fasta", ".fa",".fsa", ".fna", ".seq",".pep"):
        return "FASTA"
    elif filename.endswith(".fastq",".fq"):
        return "FASTQ"
    else:
        return None

=======

=======
def read_fasta(path_or_file) -> Iterator[Tuple[str, str]]:
    """
    Yield (header, sequence) from a FASTA file.
    Header returned without leading '>'.
    WHole file will not be loaded into memory.
    """
    fh = _open(path_or_file, "r")
    header: Optional[str] = None
    seq_chunks = []
    for line in fh:
        line = line.strip()
        if not line:
            continue
        if line.startswith(">"):
            if header is not None:
                yield header, "".join(seq_chunks)
            header = line[1:].strip()
            seq_chunks = []
        else:
            seq_chunks.append(line)
    if header is not None:
        yield header, "".join(seq_chunks)
=======
def read_fastq(path_or_file) -> Iterator[Tuple[str, str, str]]:
    """
    Yield (header, sequence, quality) from a FASTQ file.
    Header returned without leading '@'.
    """
    fh = _open(path_or_file, "r")
    line_num = 0
    while True:
        h = fh.readline(); line_num += 1
        if not h:
            break
        seq = fh.readline(); line_num += 1
        plus = fh.readline(); line_num += 1
        qual = fh.readline(); line_num += 1
        if not (seq and plus and qual):
            raise ValueError(f"Truncated FASTQ record near line {line_num}.")
        if not h.startswith("@") or not plus.startswith("+"):
            raise ValueError(f"Invalid FASTQ format near line {line_num-3}.")
        s = seq.strip(); q = qual.strip()
        if len(s) != len(q):
            raise ValueError(f"FASTQ sequence/quality length mismatch near line {line_num-3}.")
        yield h[1:].strip(), s, q


