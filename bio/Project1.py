def detect_format(filename):
    # Check file extension to know if it's FASTA or FASTQ
    if filename.endswith((".fasta", ".fa", ".fna", ".txt")):
        return "FASTA"
    elif filename.endswith(".fastq"):
        return "FASTQ"
    else:
        return None


def read_fasta(file):
    f = open(file, "r")
    lines = []
    for line in f:
        line = line.strip()
        if line != "":
            lines.append(line)
    f.close()

    if len(lines) == 0:
        print("Error: The file is empty.")
        return None

    if not lines[0].startswith(">"):
        print("Error: FASTA file must start with '>'.")
        return None

    seqs = []
    header = ""
    seq = ""

    for line in lines:
        if line.startswith(">"):
            if header != "":
                seqs.append((header, seq))
            header = line[1:]
            seq = ""
        else:
            seq = seq + line

    if header != "" and seq != "":
        seqs.append((header, seq))

    return seqs


def read_fastq(file):
    f = open(file, "r")
    lines = []
    for line in f:
        line = line.strip()
        if line != "":
            lines.append(line)
    f.close()

    if len(lines) == 0:
        print("Error: The file is empty.")
        return None

    if len(lines) % 4 != 0:
        print("Error: FASTQ files must have groups of 4 lines per record.")
        return None

    seqs = []
    i = 0
    while i < len(lines):
        header = lines[i]
        seq = lines[i + 1]
        plus = lines[i + 2]
        qual = lines[i + 3]

        if not header.startswith("@"):
            print("Error: FASTQ header must start with '@'.")
            return None
        if not plus.startswith("+"):
            print("Error: Third line must start with '+'.")
            return None
        if len(seq) != len(qual):
            print("Error: Sequence and quality must be same length.")
            return None

        seqs.append((header[1:], seq, qual))
        i = i + 4

    return seqs


def fasta_stats(data):
    total = len(data)
    total_len = 0
    for h, s in data:
        total_len += len(s)
    avg = total_len / total
    print("\nFASTA Stats:")
    print("Total sequences:", total)
    print("Total length:", total_len)
    print("Average length:", avg)


def fastq_stats(data):
    total = len(data)
    total_len = 0
    for h, s, q in data:
        total_len += len(s)
    avg = total_len / total
    print("\nFASTQ Stats:")
    print("Total reads:", total)
    print("Total bases:", total_len)
    print("Average read length:", avg)


def main():
    print("Welcome to the Sequence Analyzer Program!")
    filename = input("Enter the path to your sequence file: ").strip().strip('"')

    fmt = detect_format(filename)

    if fmt == "FASTA":
        data = read_fasta(filename)
        if data != None:
            fasta_stats(data)
    elif fmt == "FASTQ":
        data = read_fastq(filename)
        if data != None:
            fastq_stats(data)
    else:
        print("Error: Unknown file type. Must be FASTA or FASTQ.")


main()
