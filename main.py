#!/usr/bin/env python3
'''
------------------------------------------------
 Execution
  This is the main entry point for running the bioinformatics
 sequence analyzer program.

 Responsibilities:
   • Parse command-line arguments.
   • Detect and validate input files (FASTA / FASTQ).
   • Display sequence statistics.
   • Provide an interactive menu for further operations:
       - Extract random subset
       - Filter by sequence length
       - Convert FASTQ → FASTA
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
        raise SystemExit(f"Error: {input_path.name} is an empty file")

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

    # --- INTERACTIVE MENU ---
    while True:  # Loop until user makes a valid choice or exits
        print("\nChoose an operation:")
        print("  1: Extract random sequences")
        print("  2: Filter by minimum length")
        print("  3: Convert FASTQ to FASTA")
        print("  4: Exit menu")
        
        choice_str = input("Enter a number (1-4): ").strip()

        try:
            # Handle empty input
            if not choice_str:
                print("Error: No input provided. Please enter a number.")
                continue  # Go back to the start of the loop
            
            choice = int(choice_str)

            if choice == 1:
                out = args.output or default_out_name(input_path, fmt, "extract")
                ops.extract_random(input_path, fmt, out)
                print(f"Wrote random selection to: {out}")
                break  # Exit loop on success

            elif choice == 2:
                min_len = ask_min_len()

                # Handle when no sequence is filtered
                if min_len > summary["largest_len"]:
                    print(f"Provided length is greater than the largest seuence length {summary['largest_len']}. Enter a valid length")
                    continue
                
                # Write only if there is a filtered sequence
                out = args.output or default_out_name(input_path, fmt, f"filter_ge{min_len}")
                kept = ops.filter_by_min_len(input_path, fmt, out, min_len)
                print(f"Wrote {kept} sequences (len ≥ {min_len}) to: {out}")
                break  # Exit loop on success

            elif choice == 3:
                if fmt == "fasta":
                    print("FASTA file can't be converted to FASTQ file.")
                else:
                    out = args.output or input_path.with_suffix(".fasta")
                    n = ops.convert_fastq_to_fasta(input_path, out)
                    print(f"Converted {n} sequences FASTQ → FASTA: {out}")
                break  # Exit loop (even if conversion was not possible)

            elif choice == 4:
                print("Exiting interactive menu.")
                break # User chose to exit

            else:
                # Handle numbers out of range (e.g., 5, 0, -1)
                print(f"Error: Invalid choice '{choice}'. Please enter a number between 1 and 4.")
                # No break, loop continues

        except ValueError:
            # Handle non-integer input (e.g., "apple")
            print(f"Error: Unrecognized operation '{choice_str}'. Please enter a number (1-4).")
            # No break, loop continues
    
    # --- END OF INTERACTIVE MENU ---
        
#------------------minimum length checker------------------------------------
def ask_min_len() -> int:
    while True:
        s = input("Enter minimum sequence length (integer ≥ 0): ").strip()
        try:
            v = int(s)
            if v < 0:
                print("Length must be ≥ 0.")
                continue
            return v
        except ValueError:
            print("Please enter a valid integer.")

#------------------Preserving output file format------------------------------------
def default_out_name(inp: Path, fmt: str, tag: str) -> Path:
    # Preserve original format for extract/filter outputs
    if fmt == "fasta":
        return inp.with_suffix(f".{tag}.fasta")
    else:
        return inp.with_suffix(f".{tag}.fastq")

#------------------Building Parser------------------------------------
def build_parser() -> argparse.ArgumentParser:
    p = argparse.ArgumentParser(
        prog="bio-prog",
        description="Sequence statistics and operations (FASTA/FASTQ)."
    )
    sub = p.add_subparsers(dest="cmd", required=True)

    a = sub.add_parser("analyze", help="Print stats, then interactively choose an operation")
    a.add_argument("-i", "--input", required=True, help="Input FASTA or FASTQ")
    a.add_argument("-o", "--output", help="Optional output path (operation-dependent)")
    a.set_defaults(func=cmd_analyze)

    return p

#------------------Main Function------------------------------------
def main():
    args = build_parser().parse_args()
    args.func(args)

if __name__ == "__main__":
    main()
