#!/usr/bin/env python3
'''
------------------------------------------------
 I/O files (FASTA/ FASTQ readers & writers)
------------------------------------------------
'''
<<<<<<< bio/io.py

#---------------imports-----------------------
import os
from pathlib import Path
from typing import Iterator, Tuple, Iterable, TextIO, Optional

#-----------------file opening handler----------------
def _open(path_or_file, mode: str = "r") -> TextIO:
    if hasattr(path_or_file, "read"):
        return path_or_file  # already a file-like
    return open(Path(path_or_file), mode, encoding="utf-8")


#----------------extension checker---------------
def detect_extension(filename):
    # Check file extension to know if it's FASTA or FASTQ
    if filename.endswith((".fasta", ".fa",".fsa", ".fna", ".seq",".pep"):
        return "FASTA"
    elif filename.endswith(".fastq",".fq"):
        return "FASTQ"
    else:
        return None

#-------------------fasta parser---------------------
def read_fasta(path_or_file) -> Iterator[Tuple[str, str]]:
    """
    Yield (header, sequence) from a FASTA file.
    Header returned without leading '>'.
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

#-------------------write fasta sequences------------------
def write_fasta(records: Iterable[Tuple[str, str]], path_or_file, line_width: int = 80) -> None:
    fh = _open(path_or_file, "w")
    for header, seq in records:
        fh.write(f">{header}\n")
        for i in range(0, len(seq), line_width):
            fh.write(seq[i:i+line_width] + "\n")

#-------------------fastq parser--------------------------
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

#---------------------write fastq sequences----------------------
def write_fastq(records: Iterable[Tuple[str, str, str]], path_or_file) -> None:
    fh = _open(path_or_file, "w")
    for header, seq, qual in records:
        fh.write(f"@{header}\n{seq}\n+\n{qual}\n")