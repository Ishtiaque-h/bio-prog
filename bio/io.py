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

#---------------global variables------------------------------
FASTA_EXTS = {".fasta", ".fa", ".fsa", ".fna", ".pep", ".seq"}
FASTQ_EXTS = {".fastq"}

# Uppercase only; we’ll .upper() sequences before checking
DNA_CHARS     = set("ACGTNRYSWKMBDHV-")     # IUPAC DNA + gap
RNA_CHARS     = set("ACGUNRYSWKMBDHV-")     # IUPAC RNA + gap
PROTEIN_CHARS = set("ACDEFGHIKLMNPQRSTVWYBXZJUO*-")  # 20 AA + common ambiguous + stop/gap
# Union used for "DNA or RNA or peptide" as per your spec
FASTA_ALLOWED = DNA_CHARS | RNA_CHARS | PROTEIN_CHARS

# FASTQ typically nucleotides; allow IUPAC DNA/RNA (no peptides in FASTQ by convention)
FASTQ_ALLOWED = DNA_CHARS | RNA_CHARS

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

#----------------fasta validator---------------------  
def validate_fasta(path: Path) -> None:
    """
    Enforces:
      1) Header lines start with '>' and immediately a non-space name char.
      2) The immediate next line after a header must begin the sequence (non-empty).
      3) Sequence names must be unique (name = token until first whitespace).
      4) Each record must have at least 2 lines (header + >=1 sequence line).
      5) All sequence lines contain only DNA/RNA/protein characters (uppercase-insensitive).
    """
    seen_names = set()
    with open(path, "r", encoding="utf-8") as fh:
        lines = fh.readlines()

    if not lines:
        raise ValueError("No content found.")

    i = 0
    record_count = 0
    while i < len(lines):
        line = lines[i].rstrip("\n")
        i += 1

        if not line.strip():
            continue  # skip stray empty lines between records

        # Expect a header
        if not line.startswith(">"):
            raise ValueError(f"FASTA error at line {i}: expected '>' header, found: {line[:20]!r}")

        # Rule 1: '>' immediately followed by non-space
        if len(line) == 1 or line[1].isspace():
            raise ValueError(f"FASTA error at line {i}: header must have name immediately after '>' with no leading space.")

        # Extract name until first whitespace
        header = line[1:]
        name = header.split()[0]
        if name in seen_names:
            raise ValueError(f"FASTA error at line {i}: duplicate sequence name '{name}'.")
        seen_names.add(name)

        # Rule 2 & 4: immediate next non-empty line must be sequence, and there must be at least one
        if i >= len(lines):
            raise ValueError(f"FASTA error after line {i}: header '{name}' has no sequence line.")
        # Allow empty lines, but the first non-empty after header must be sequence content
        j = i
        while j < len(lines) and not lines[j].strip():
            j += 1
        if j >= len(lines):
            raise ValueError(f"FASTA error after line {i}: header '{name}' has no sequence line.")
        if lines[j].startswith(">"):
            raise ValueError(f"FASTA error at line {j+1}: expected sequence line, found another header.")

        # Now collect sequence lines until next header or EOF; validate chars
        seq_len_this = 0
        while j < len(lines):
            s = lines[j].strip()
            if not s:
                j += 1
                continue
            if s.startswith(">"):
                break
            # Rule 5: validate allowed chars
            su = s.upper()
            if any(ch not in FASTA_ALLOWED for ch in su):
                raise ValueError(f"FASTA error at line {j+1}: invalid character(s) in sequence.")
            seq_len_this += len(su)
            j += 1

        if seq_len_this == 0:
            raise ValueError(f"FASTA error: header '{name}' has no sequence content.")
        record_count += 1
        i = j  # continue from next header or EOF

    if record_count == 0:
        raise ValueError("No FASTA records found.")

#----------------fastq validator---------------------  
def validate_fastq(path: Path) -> None:
    """
    Enforces:
      1) Records are in 4-line groups.
      2) Line 1: '@' immediately followed by name (no space right after '@'); name is until first space.
      3) Line 2: sequence (non-empty), DNA/RNA characters only (IUPAC allowed).
      4) Line 3: starts with '+' (free text allowed after).
      5) Line 4: quality string, same length as sequence.
      6) Sequence names must be unique.
    """
    seen = set()
    with open(path, "r", encoding="utf-8") as fh:
        line_no = 0
        while True:
            h = fh.readline(); line_no += 1
            if not h:
                break  # EOF cleanly on record boundary
            s = fh.readline(); line_no += 1
            p = fh.readline(); line_no += 1
            q = fh.readline(); line_no += 1

            if not (s and p and q):
                raise ValueError(f"FASTQ error near line {line_no}: truncated record.")

            if not h.startswith("@"):
                raise ValueError(f"FASTQ error at line {line_no-3}: header must start with '@'.")
            if len(h.strip()) < 2 or h[1].isspace():
                raise ValueError(f"FASTQ error at line {line_no-3}: name must immediately follow '@' (no leading space).")
            header = h[1:].strip()
            name = header.split()[0]
            if name in seen:
                raise ValueError(f"FASTQ error at line {line_no-3}: duplicate sequence name '{name}'.")
            seen.add(name)

            seq = s.strip().upper()
            if not seq:
                raise ValueError(f"FASTQ error at line {line_no-2}: empty sequence line.")
            if any(ch not in FASTQ_ALLOWED for ch in seq):
                raise ValueError(f"FASTQ error at line {line_no-2}: invalid character(s) in sequence.")

            if not p.startswith("+"):
                raise ValueError(f"FASTQ error at line {line_no-1}: third line must start with '+'.")

            qual = q.strip()
            if len(qual) != len(seq):
                raise ValueError(f"FASTQ error at line {line_no}: quality length {len(qual)} != sequence length {len(seq)}.")

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