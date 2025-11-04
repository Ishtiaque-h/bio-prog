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
                if in_seq:
                    # Error: We are inside a record, blank lines are not allowed here.
                    raise ValueError(
                        f"FASTA error at line {line_no}: blank lines are not permitted "
                        f"within a sequence record (after header '{current_name}')."
                    )
                else:
                    # We are between records (or at file start). This is fine.
                    continue

            if line.startswith(">"):  # header line
                # If we were in a record, ensure it had at least one sequence line
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
                seq_len_this = 0
                continue

            # Sequence line processing (multi-line allowed)
            # If we see a sequence line but weren't in a record, it's an error
            if not in_seq:
                raise ValueError(f"FASTA error at line {line_no}: expected header starting with '>' before sequence data.")

            # Rule 3: No whitespace within sequence lines
            if any(ch.isspace() for ch in line):
                raise ValueError(f"FASTA error at line {line_no}: whitespaces are not allowed within sequence lines.")

            # Allowed alphabet check (case-insensitive)
            su = line.upper()
            if any(ch not in FASTA_ALLOWED for ch in su):
                raise ValueError(f"FASTA error at line {line_no}: invalid character(s) found in sequence.")

            seq_len_this += len(su)     # We found valid sequence data

    # EOF: if file ended while inside a record, ensure we saw at least one sequence line
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

            # Find the next non-empty line (must be a header).
            # This loop skips blank lines at the file start and between records.
            while h and not h.rstrip('\n'):
                h = fh.readline(); line_no += 1

            if not h:
                break       # EOF cleanly on record boundary
            
            # Read the next 3 lines. They CANNOT be blank.
            s = fh.readline(); line_no += 1     
            p = fh.readline(); line_no += 1
            q = fh.readline(); line_no += 1

            # Check for blank lines within the record
            #    (We already know 'h' is not blank)
            if not s.rstrip('\n'):
                raise ValueError(f"FASTQ error at line {line_no-2}: blank line not allowed (after header).")
            if not p.rstrip('\n'):
                raise ValueError(f"FASTQ error at line {line_no-1}: blank line not allowed (after sequence).")
            if not q.rstrip('\n'):
                raise ValueError(f"FASTQ error at line {line_no}: blank line not allowed (after plus line).")

            # Rule 1: Records are 4-line groups. Handle truncated records.
            # Check for truncated file (EOF within a record)
            # Note: 'h' is guaranteed non-empty, so we just check s, p, q
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
        raise ValueError("FASTQ error: no valid records found.")
    
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
    """
    fh = _open(path_or_file, "r")

    while True:
        h = fh.readline()

        '''
        Find the next non-empty line (the header)
        This loop skips all blank lines between records.
        '''
        while h and not h.rstrip('\n'):
            h = fh.readline()

        # If no line was found, we're at the end of the file.
        if not h:
            break   # clean EOF
        '''
        Read the next three lines.
        We ASSUME these exist and are not blank,
        Because validate_fastq() already checked.
        '''
        seq = fh.readline()
        plus = fh.readline()
        qual = fh.readline()
    
        # Check only for a truncated file (catches unexpected EOF)
        if not (seq and plus and qual):
            raise ValueError(f"Incomplete FASTQ record (unexpected EOF).")
        '''
        Yield the data.
        We use .strip() to clean up header/seq/qual lines.
        The validator already checked for internal whitespace,
        so this is safe.
        '''
        yield h[1:].strip(), seq.strip(), qual.strip()

#---------------------write fastq sequences----------------------
def write_fastq(records: Iterable[Tuple[str, str, str]], path_or_file) -> None:
    fh = _open(path_or_file, "w")
    for header, seq, qual in records:
        fh.write(f"@{header}\n{seq}\n+\n{qual}\n")