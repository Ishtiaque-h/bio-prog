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

