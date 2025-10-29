#!/usr/bin/env python3
'''
------------------------------------------------
 I/O files (FASTA/ FASTQ readers & writers)
------------------------------------------------
'''

#---------------imports------------------------------
import os
from pathlib import Path
from typing import Iterator, Tuple, Iterable, TextIO, Optional

#-----------------file opening handler----------------
def _open(path_or_file, mode: str = "r") -> TextIO:
    """
    Opens a file or returns a file-like object.

    - If path_or_file is a path, it opens it.
    - If it's already file-like, it returns it.
    """
    if hasattr(path_or_file, "read"):
        return path_or_file  # already a file-like
    return open(Path(path_or_file), mode, encoding="utf-8")

#----------------extension checker---------------------
def detect_extension(path_or_file: Path, checked_fmt: str) -> None:
    """
    Non-fatal: warn if extension doesn’t match typical sets.
    """
    ext = path_or_file.suffix.lower()
    if checked_fmt == "fasta" and ext and ext not in FASTA_EXTS:
        print(f"Wrning!!!: Detected FASTA content but extension '{ext}' is uncommon for FASTA.")
    if checked_fmt == "fastq" and ext and ext not in FASTQ_EXTS:
        print(f"WARNING!!!: Detected FASTQ content but extension '{ext}' is uncommon for FASTQ.")

#----------------fasta/fastq format checker---------------------    
def check_format(path_or_file) -> str:
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

#-------------------fasta parser---------------------
def read_fasta(path_or_file) -> Iterator[Tuple[str, str]]:
    """
    Yield (header, sequence) from a FASTA file.
    Header returned without leading '>'.
    """
    fh = _open(path_or_file, "r")
    header = None
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
def write_fasta(records: Iterable[Tuple[str, str]], path_or_file) -> None:
    fh = _open(path_or_file, "w")
    for header, seq in records:
        fh.write(f">{header}\n")
        fh.write(seq + "\n")

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