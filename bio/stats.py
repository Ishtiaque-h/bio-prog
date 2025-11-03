#!/usr/bin/env python3
'''
------------------------------------------------
 Statistics and data analytics
------------------------------------------------
'''

 #---------------imports---------------------------------------------------------------
from .io import _open                                # file utility from IO module
from typing import Tuple
from typing import List, Dict, Any, Iterable, Tuple
from statistics import mean                         # built-in mean for averages
from .utils import gc_and_counts                    # helper: computes GC and N counts
 
#------------rounding decimal----------------------------------------------------------
def _round2(x: float) -> float: # Round a float to two decimal places
    return round(x + 1e-12, 2) 

#-------------summarizing fasta-----------------------
def summarize_fasta(records: Iterable[Tuple[str, str]]) -> Dict[str, Any]:
    items: List[Dict[str, Any]] = []    # Process each sequence and collect statistics
    for h, s in records:
        gc, atgc, n = gc_and_counts(s)
        gc_frac = (gc / atgc) if atgc > 0 else 0.0  #Calculate GC fraction (avoid division by zero)
        items.append({"name": h, "len": len(s), "gc_frac": gc_frac, "n": n})

    if not items:
        raise ValueError("No sequences found.")
            
    lengths = [x["len"] for x in items] #Extract all lengths for statistical calculations
    avg_len = _round2(mean(lengths))
    max_len = max(lengths)
    min_len = min(lengths)

    largest = [x["name"] for x in items if x["len"] == max_len][:10]   # Find sequences with extreme lengths (up to 10 of each)
    smallest = [x["name"] for x in items if x["len"] == min_len][:10]

    avg_gc = _round2(mean(x["gc_frac"] for x in items) * 100.0)        # Calculate average GC percentage across all sequences
    avg_ns = _round2(mean(x["n"] for x in items))

    return {
        "count": len(items),         # total number of sequences
        "avg_len": avg_len,          # average sequence length
        "largest_names": largest,    # up to 10 names with longest length
        "largest_len": max_len,      # length of the longest sequence
        "smallest_names": smallest,  # up to 10 names with shortest length
        "smallest_len": min_len,     # length of the shortest sequence
        "avg_gc_percent": avg_gc,    # averaged per-seq GC%, Ns excluded from per-seq denominator
        "avg_ns_per_seq": avg_ns,    # average count of 'N' bases per sequence
    }

#----------------summarizing fastq-----------------------------
def summarize_fastq(records: Iterable[Tuple[str, str, str]]) -> Dict[str, Any]:
    # Reuses FASTA summary logic by converting FASTQ format
    # Quality scores are ignored for statistical analysis
    # Provides same metrics as summarize_fasta() for consistency
    
    def fastq_to_fasta_iter():
        for h, s, _q in records:
            yield h, s
    return summarize_fasta(fastq_to_fasta_iter())

