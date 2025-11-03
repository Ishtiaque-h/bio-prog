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
FASTA_EXTS = {".fasta", ".fa", ".fsa", ".fna", ".pep", ".seq"}  # valid FASTA extensions
FASTQ_EXTS = {".fastq", ".fq"}                                  # valid FASTQ extensions

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
    - This allows flexibility in accepting both file paths and file handles
    """
    if hasattr(path_or_file, "read"):
        return path_or_file  # already a file-like
    return open(Path(path_or_file), mode, encoding="utf-8")

#----------------extension checker---------------------
def detect_extension(path_or_file: Path, checked_fmt: str) -> None:
    """
    Warn if file extension doesn’t match detected format.
    Non-fatal — continues execution.
    Helps catch potential file naming issues early.
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
    Auto-detects file format by examining the first non-blank line.
    Uses '>' for FASTA and '@' for FASTQ as distinguishing markers.
    """
    fh = _open(path_or_file, "r")
    pos = fh.tell() # Save position to reset file pointer after peeking
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
        break   # Stop after first non-empty line
    fh.seek(pos)
    raise ValueError("Unable to determine file format (expected FASTA or FASTQ).")

#----------------fasta validator---------------------  
def validate_fasta(path: Path) -> None:
    """
    FASTA validation:
      1) Each record begins with '>' immediately followed by a non-space name.
      2) No blank lines anywhere in the file.
      3) Sequence lines must contain only residue characters (DNA/RNA/protein),
      with NO spaces or tabs anywhere in the sequence.
      4) Sequence names must be unique (first token after '>').
      5)  ≥1 sequence line per record.
      # Comprehensive validation ensures file integrity before processing.
      # Tracks sequence names to enforce uniqueness constraint.
    """
    seen_names = set()  # Track unique names to prevent duplicates*
    line_no = 0
    in_seq = False
    seq_len_this = 0
    current_name = None

    with open(path, "r", encoding="utf-8") as fh:
        for raw in fh:
            line_no += 1
            line = raw.rstrip("\n")

            # Rule 2: no empty lines allowed anywhere
            if not line.strip():
                raise ValueError(f"FASTA error at line {line_no}: blank lines are not permitted.")

            if line.startswith(">"):  # header line
                # Before starting new record, verify previous record had sequence data
                if in_seq and seq_len_this == 0:
                    raise ValueError(f"FASTA error at line {line_no}: header '{current_name}' had no sequence lines.")

                # Rule 1: '>' immediately followed by non-space name
                if len(line) == 1 or line[1].isspace():
                    raise ValueError(f"FASTA error at line {line_no}: header must have sequence name immediately after '>' (no leading space).")

                header = line[1:].strip()
                name = header.split()[0]    # Extract first token as sequence name

                # Rule 4: unique names
                if name in seen_names:
                    raise ValueError(f"FASTA error at line {line_no}: duplicate sequence name '{name}'.")
                seen_names.add(name)

                # Start new record
                current_name = name
                in_seq = True
                seq_len_this = 0    # Reset counter for new sequence
                continue

            # Sequence line (multi-line allowed)
            if not in_seq:
                raise ValueError(f"FASTA error at line {line_no}: expected header starting with '>' before sequence data.")

            # Rule 3: No whitespace within sequence lines
            if any(ch.isspace() for ch in line):
                raise ValueError(f"FASTA error at line {line_no}: whitespaces are not allowed within sequence lines.")

            # Allowed alphabet check (case-insensitive)
            su = line.upper()
            if any(ch not in FASTA_ALLOWED for ch in su):
                raise ValueError(f"FASTA error at line {line_no}: invalid character(s) found in sequence.")

            seq_len_this += len(su)

    # EOF check - if file ended while inside a record, ensure we saw at least one sequence line
    if in_seq and seq_len_this == 0:
        raise ValueError(f"FASTA error: header '{current_name}' has no sequence content.")

    if not seen_names:
        raise ValueError("FASTA error: no valid records found.")

#----------------fastq validator---------------------  
def validate_fastq(path: Path) -> None:
    """
    Enforces:
      1) Records are in 4-line groups.
      2) Line 1: '@' immediately followed by name (no space right after '@'); name is until first space.
      3) Line 2: sequence (non-empty, non-spaced, single-line), DNA/RNA characters only (IUPAC allowed).
      4) Line 3: starts with '+' (additional text allowed after).
      5) Line 4: quality string — same length as seq, printable ASCII 33–126 only, no whitespace
      6) Sequence names must be unique.
      # Strict 4-line format validation ensures FASTQ integrity
      # Quality scores validated for proper Phred+33 encoding
    """
    seen_names = set()  # Track unique names to prevent duplicates
    with open(path, "r", encoding="utf-8") as fh:
        line_no = 0
        while True:
            h = fh.readline(); line_no += 1

            if not h:
                break       # EOF cleanly on record boundary
            
            # No empty lines allowed anywhere.
            if not h.rstrip('\n'):
                raise ValueError(f"FASTQ error at line {line_no}: blank lines are not permitted.")

            s = fh.readline(); line_no += 1     # Sequence line
            p = fh.readline(); line_no += 1     # Plus line (+)
            q = fh.readline(); line_no += 1     # Quality line

            # Rule 1: Records are 4-line groups. Handle truncated records.
            if not (s and p and q):
                raise ValueError(f"FASTQ error near line {line_no}: incomplete records are not permitted.")

            # Rule 2: line 1 includes header
            if not h.startswith("@"):
                raise ValueError(f"FASTQ error at line {line_no-3}: header must start with '@'.")
            
            # seq name must immediately follow '@' (no space)
            if len(h.rstrip('\n')) == 1 or h[1].isspace():
                raise ValueError(f"FASTQ error at line {line_no-3}: header must have sequence name immediately after '@' (no leading space).")
            header = h[1:].strip()
            name = header.split()[0]    # Extract first token as sequence name

            # Rule 6: Sequence names must be unique
            if name in seen_names:
                raise ValueError(f"FASTQ error at line {line_no-3}: duplicate sequence name '{name}'.")
            seen_names.add(name)

            # Rule 3: Allow non-empty sequence in line 2 with neucleotides only (no space, no multiline)
            seq_raw = s.rstrip("\n")

            # No whitespace within sequence line
            if any(ch.isspace() for ch in seq_raw):
                raise ValueError(f"FASTQ error at line {line_no-2}: whitespaces are not allowed in sequence.")
            seq = seq_raw.upper()
            
            # Allow non-empty sequence line
            if not seq:
                raise ValueError(f"FASTQ error at line {line_no-2}: empty sequence line.")

            # Allow nucleotide characters only
            if any(ch not in FASTQ_ALLOWED for ch in seq):
                raise ValueError(f"FASTQ error at line {line_no-2}: invalid character in sequence.")

            # Rule 4: Line 3 starts with '+' (free text allowed after).
            plus_raw = p.rstrip("\n")
            if not plus_raw.startswith("+"):
                raise ValueError(f"FASTQ error at line {line_no-1}: third line must start with '+'.")

            # Rule 5: Line 4 is quality string — same length as seq, printable ASCII 33–126 only, no whitespace
            qual = q.rstrip("\n")
            if len(qual) != len(seq):
                raise ValueError(f"FASTQ error at line {line_no}: quality length {len(qual)} != sequence length {len(seq)}.")
            
            # No whitespace within quality line
            if any(ch.isspace() for ch in qual):
                raise ValueError(f"FASTQ error at line {line_no}: whitespaces are not allowed in quality.")
            
            # Phred+33 typical printable ASCII range '!' (33) to '~' (126)
            if any(not (33 <= ord(ch) <= 126) for ch in qual):
                raise ValueError(f"FASTQ error at line {line_no}: quality contains non-printable characters (must be ASCII 33–126).")
        
    if not seen_names:
        raise ValueError("FASTA error: no valid records found.")
    
#-------------------fasta parser---------------------
def read_fasta(path_or_file) -> Iterator[Tuple[str, str]]:
    """
    Yield (header, sequence) from a FASTA file.
    Header returned without leading '>'.
    Generator function for memory-efficient processing of large files.
    Automatically handles multi-line sequences by concatenating them.

    """
    fh = _open(path_or_file, "r")
    header = None
    seq_chunks = [] # Collect sequence lines for multi-line sequences
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
            seq_chunks.append(line)  # Accumulate sequence lines
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
    Generator function for memory-efficient processing of large FASTQ files.
    Reads exactly 4 lines per record to maintain FASTQ structure
    """
    fh = _open(path_or_file, "r")
    line_num = 0
    while True: # Read all 4 lines of FASTQ record
        h = fh.readline(); line_num += 1
        if not h:
            break
        seq = fh.readline(); line_num += 1
        plus = fh.readline(); line_num += 1
        qual = fh.readline(); line_num += 1 
        if not (seq and plus and qual): # Validate we got all 4 lines
            raise ValueError(f"Truncated FASTQ record near line {line_num}.")
        if not h.startswith("@") or not plus.startswith("+"):
            raise ValueError(f"Invalid FASTQ format near line {line_num-3}.")
        s = seq.strip(); q = qual.strip()
        if len(s) != len(q):    # Ensure sequence and quality have matching lengths
            raise ValueError(f"FASTQ sequence/quality length mismatch near line {line_num-3}.")
        yield h[1:].strip(), s, q

#---------------------write fastq sequences----------------------
def write_fastq(records: Iterable[Tuple[str, str, str]], path_or_file) -> None:
    fh = _open(path_or_file, "w")
    for header, seq, qual in records:
        fh.write(f"@{header}\n{seq}\n+\n{qual}\n")