#!/usr/bin/env python3
'''
------------------------------------------------
 Execution
------------------------------------------------
'''

#----------------import libraries-------------------------
import argparse
from pathlib import Path
from bio import io
from bio import stats as stats_mod
from tasks import ops


#----------------variables--------------------------
BANNER = "=" * 60

#----------------show outputs---------------------
def print_stats(d: dict) -> None:
    print(BANNER)
    print("File statistics")
    print(BANNER)
    print(f"Total sequences         : {d['count']}")
    print(f"Average length          : {d['avg_len']}")
    print(f"Largest length          : {d['largest_len']}")
    print("Largest sequence name(s): " + (", ".join(d['largest_names']) if d['largest_names'] else "-"))
    print(f"Smallest length         : {d['smallest_len']}")
    print("Smallest sequence name(s): " + (", ".join(d['smallest_names']) if d['smallest_names'] else "-"))
    print(f"Average GC-content (%)  : {d['avg_gc_percent']}")
    print(f"Average # of Ns/sequence: {d['avg_ns_per_seq']}")
    print(BANNER)

#------------------cmd interactive parser------------------------------------
def cmd_analyze(args: argparse.Namespace) -> None:
    input_path = Path(args.input)
    # Fails fast if the input path/file is empty.
    if not input_path.exists():
        raise SystemExit(f"Error: input file not found: {input_path}")
    if input_path.stat().st_size == 0:
        raise ValueError("Input file is empty.")

    try:
        fmt = io.check_format(input_path)  # 'fasta' or 'fastq'
    except Exception as e:
        raise SystemExit(f"Error detecting file format: {e}")
    
    """
    for invalid extension we are generating a soft warning, not restricting further.
    We will let the user to access a correct file with accidentally wrong extension
    """
    io.detect_extension(input_path, fmt) # soft warning 

    try:
        if fmt == "fasta":
            io.validate_fasta(input_path)
        else:
            io.validate_fastq(input_path)
    except Exception as e:
        raise SystemExit(f"Corrupted or invalid {fmt.upper()} file: {e}")

    try:
        if fmt == "fasta":
            summary = stats_mod.summarize_fasta(io.read_fasta(input_path))
        else:
            summary = stats_mod.summarize_fastq(io.read_fastq(input_path))
    except Exception as e:
        raise SystemExit(f"Error while reading/parsing sequences: {e}")

    print_stats(summary)

    # Interactive menu

    print("Choose an operation: 1 for extract | 2 for filter | 3 for convert")

    choice = int(input("Enter 1/2//3> ").strip())

    if choice == 1:

        out = args.output or default_out_name(input_path, fmt, "extract")

        ops.extract_random(input_path, fmt, out)

        print(f"Wrote random selection to: {out}")

    elif choice == 2:

        min_len = ask_min_len()

        out = args.output or default_out_name(input_path, fmt, f"filter_ge{min_len}")

        kept = ops.filter_by_min_len(input_path, fmt, out, min_len)

        print(f"Wrote {kept} sequences (len ≥ {min_len}) to: {out}")

    elif choice == 3:

        if fmt == "fasta":

            print("FASTA file can't be converted to FASTQ file.")

        else:

            out = args.output or input_path.with_suffix(".fasta")

            n = ops.convert_fastq_to_fasta(input_path, out)

            print(f"Converted {n} sequences FASTQ → FASTA: {out}")

    else:

        print("Unrecognized operation. Please choose: extract | filter | convert")