
---

## 🧬 Project 1 — Sequence File Statistics and Manipulation

### 🎯 Objective

A Python program that reads sequence files in **FASTA** or **FASTQ** format, computes essential statistics, and performs interactive operations such as extracting, filtering, or converting sequences.
The project reinforces skills in file handling, data validation, data analysis, and user interaction within a bioinformatics context.

---

### 📂 Repository Structure

```
project1_sequence_file_manipulation/
├── README.md
├── .gitignore
├── main.py
├── bio/
│   ├── io.py            # FASTA/FASTQ I/O + validators
│   ├── utils.py         # GC-content and helper functions
│   └── stats.py         # Statistical computations
├── tasks/
│   └── ops.py           # extract / filter / convert operations
├── data/
│   ├── example.fasta
│   └── example.fastq
└── tests/
```

---

### 🧠 Key Features & Requirements Implemented

| Requirement                    | Implementation                                                                                                                        |
| ------------------------------ | ------------------------------------------------------------------------------------------------------------------------------------- |
| **Support both FASTA & FASTQ** | Automatic format detection via first non-empty line (`>` or `@`).                                                                     |
| **Error handling**             | Handles missing, empty, and corrupted files with clear messages.                                                                      |
| **Validation rules**           | Enforces correct structure and header rules for both formats (see validators below).                                                  |
| **Statistics output**          | Prints total count, average length, largest/smallest sequences (up to 10 names), average GC%, and average number of N’s per sequence. |
| **Interactive operations**     | Prompts user for `extract`, `filter`, or `convert` (case-insensitive).                                                                |
| **Extract**                    | Randomly selects ≤ 25 sequences, writes in the same format.                                                                           |
| **Filter**                     | Keeps sequences ≥ user-defined length, writes in the same format.                                                                     |
| **Convert**                    | FASTQ → FASTA supported; FASTA → error with exact required message.                                                                   |
| **Soft warnings**              | For invalid extensions or unexpected alphabets (IUPAC non-standard).                                                                  |
| **Streaming & efficient**      | Linear runtime, minimal memory usage, clean modular design.                                                                           |

---

### ⚙️ Installation & Usage

```bash
# 1. Clone the repository
git clone http://beadle.dbi.udel.edu:9000/ishtiaq/project1_sequence_file_manipulation.git
cd project1_sequence_file_manipulation

# 2. Optional: create a virtual environment
python3 -m venv .venv && source .venv/bin/activate
```

#### ▶ Run the program

```bash
python main.py analyze -i data/example.fasta
```

Example output:

```
============================================================
File statistics
============================================================
Total sequences         : 120
Average length          : 512.36
Largest length          : 890
Largest sequence name(s): seq_0042, seq_0711
Smallest length         : 300
Smallest sequence name(s): seq_0098
Average GC-content (%)  : 48.73
Average # of Ns/sequence: 0.12
============================================================
Choose an operation (type keyword): extract | filter | convert
> extract
Wrote random selection to: data/example.extract.fasta
```

---

### 🧪 Validation Rules

#### FASTA

1. Header line begins with `>` immediately followed by sequence name (no space).
2. Next non-empty line after a header must be sequence (not another header).
3. Sequence names must be **unique**.
4. Each record ≥ 2 lines (header + sequence).
5. Sequence lines may include only DNA/RNA/Protein characters (IUPAC; case-insensitive).

#### FASTQ

1. Records are always **4 lines**.
2. Line 1: `@` immediately followed by name (no space).
3. Line 2: sequence (non-empty, DNA/RNA characters).
4. Line 3: starts with `+`.
5. Line 4: quality string same length as sequence.
6. Sequence names must be **unique**.

**Soft warnings:**

* Extension mismatch produces a warning, not an error.

---

### 🧮 Statistical Computation

* Total number of sequences
* Average sequence length (2 decimal places)
* Largest & smallest sequence names (up to 10 ties)
* Average GC-content (per-sequence GC%, averaged, Ns excluded)
* Average number of N’s per sequence

---

### ⚡ Operations

| Operation   | Description                                      | Output                                  |
| ----------- | ------------------------------------------------ | --------------------------------------- |
| **extract** | Randomly selects ≤ 25 sequences, retains format. | `input.extract.fasta` / `.fastq`        |
| **filter**  | Keeps sequences ≥ user-provided min length.      | `input.filter_ge<MIN>.fasta` / `.fastq` |
| **convert** | FASTQ → FASTA; FASTA → prints error message.     | `input.fasta`                           |

---

### 🧩 Development Notes & Testing Summary

*(compiled from internal progress notes)*

#### Part 1 – I/O & Validation

* **Format detection:** checked first non-empty line for `>` or `@`.
* **File-not-found:** tested with nonexistent file → error handled.
* **Empty document:** tested empty `.txt` → correct “empty file” error.
* **Corrupted content:**

  * FASTA or FASTQ validators enforced standard sequencing rules and exceptions were handles carefully.
  * FASTA or FASTQ file with gibberish after “>” or "@" correctly flagged as invalid.
  * Both file with unwanted empty lines and whitespaces were handled carefully.
  * Both with illegal characters in sequence raised a line-specific error.
* **Truncated file tests:** single `@` or `>` line correctly reported as “truncated record.”
* **Extension handling:** invalid extensions trigger a soft warning; program proceeds.

#### Part 2 – Statistics

* Header info ignored for calculations (some FASTA may lack headers).
* Used per-sequence stats for length, GC%, and N counts.

#### Part 3 – Interactive Operations

* Interactive menu uses simple number values as options to avoid typos.
* All possible exceptions were handled gracefully.
* `extract`: uses `random.sample()` to pick 25 sequences.
* `filter`: user-defined min length → sequences ≥ threshold written to output;
  if threshold > largest length, outputs “0 sequences” and doesn’t write a file.
* `convert`: FASTQ → FASTA implemented by keeping first 2 lines; reverse conversion disallowed with proper message.

All tests confirmed graceful failure handling — no unhandled exceptions or crashes.

---

### 🧠 Design Choices & Good Practices

* **Streaming validation:** one-pass, constant memory even for large files.
* **Strict structure + soft content validation:** only structural corruption halts execution.
* **Fail fast:** detect obvious corruption early.
* **Clear user feedback:** all errors include line numbers and context.
* **Git discipline:** modular commits and meaningful PR reviews (implementation vs. review roles maintained).

---

### 🧮 Complexity

| Operation           | Time | Space |
| ------------------- | ---- | ----- |
| FASTA/FASTQ parsing | O(n) | O(1)  |
| Validation          | O(n) | O(1)  |
| Stats aggregation   | O(n) | O(1)  |
| Extract (default)   | O(n) | O(n)  |
| Filter / Convert    | O(n) | O(1)  |

---

### 🧑‍🤝‍🧑 Team

| Role      | Name                | Responsibility                                  |
| --------- | ------------------- | ----------------------------------------------- |
| Lead      | Ishtiaque           | Design, Implementation & testing                |
| Developer | Zainab              | Implementation & testing                        |
| Developer | Eric                | Implementation & testing                        |
| Reviewer  | Ishtiaque           | Code review, documentation, validation strategy |

---

### 🧱 Extensibility Ideas

* Add reservoir sampling for very large FASTA/FASTQ datasets.
* Integrate unit tests with `pytest` for CI workflows.

---

### 👍 Acknowledgement

AI tools (ChatGPT & Gemini) was used to design, test and improve code; to preapre README file.  

---

### 🧾 License

Developed for **BINF690:Programming in Bioinformatics** graduate course.
© 2025 Team <ishtiaque, Zainab and Eric>. For academic use only.

---
