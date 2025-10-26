## What
Short summary of changes.

## Why
Tie to project requirement (quote the bullet from the spec).

## How to make a PR
Key design decisions (parsing, rounding, error handling).

```bash
git checkout -b feat/branch_name
# add files
git add .
git commit -m "chore: scaffold project structure"
git push -u origin feat/scaffold
```
Make sure to change the **"branch_name"** for every commits

## How to test PR
- [ ] Tiny FASTA happy-path
- [ ] Tiny FASTQ happy-path
- [ ] Edge case (describe)
- [ ] Manual run transcript

## How to review locally
- open PR → you review → merge
- loccally:

``'bash
# fresh clone
git clone <repo>
cd repo
python -m venv .venv && source .venv/bin/activate  # or your usual
python -m pip install -r requirements.txt          # likely empty/minimal

# check PR branch
git fetch origin
git checkout <branch-name>

# run minimal demos (PR description should list exact commands)
python main.py analyze -i data/tiny.fasta
python main.py analyze -i data/tiny.fastq

```

## Screenshots / Logs
Paste the exact console snippet of stats output (redact paths if needed).

## Checklist
- [ ] Reads both FASTA/FASTQ or is scoped to module
- [ ] Clean, simple, modular (no over-engineering)
- [ ] Docstrings/comments for non-obvious parts
- [ ] Handles invalid input gracefully
- [ ] Matches rounding & cap-on-names rules
