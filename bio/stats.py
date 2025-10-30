#!/usr/bin/env python3
'''
------------------------------------------------
 Statistics and data analytics
------------------------------------------------
'''

 #---------------imports------------------------------
from .io import _open
from typing import Tuple
from typing import List, Dict, Any, Iterable, Tuple
from statistics import mean
from .utils import gc_and_counts
 
#------------rounding decimal-------------------------
def _round2(x: float) -> float:
    return round(x + 1e-12, 2) 

#-------------summarizing fasta-----------------------
def summarize_fasta(records: Iterable[Tuple[str, str]]) -> Dict[str, Any]:
    items: List[Dict[str, Any]] = []
    for h, s in records:
        gc, atgc, n = gc_and_counts(s)
        gc_frac = (gc / atgc) if atgc > 0 else 0.0
        items.append({"name": h, "len": len(s), "gc_frac": gc_frac, "n": n})    if not items:
        raise ValueError("No sequences found.")    lengths = [x["len"] for x in items]
    avg_len = _round2(mean(lengths))
    max_len = max(lengths)
    min_len = min(lengths)    largest = [x["name"] for x in items if x["len"] == max_len][:10]
    smallest = [x["name"] for x in items if x["len"] == min_len][:10]    avg_gc = _round2(mean(x["gc_frac"] for x in items) * 100.0)
    avg_ns = _round2(mean(x["n"] for x in items))    return {
        "count": len(items),
        "avg_len": avg_len,
        "largest_names": largest,
        "largest_len": max_len,
        "smallest_names": smallest,
        "smallest_len": min_len,
        "avg_gc_percent": avg_gc,   # averaged per-seq GC%, Ns excluded from per-seq denominator
        "avg_ns_per_seq": avg_ns,
    } 

#----------------summarizing fastq-----------------------------
def summarize_fastq(records: Iterable[Tuple[str, str, str]]) -> Dict[str, Any]:
    # identical logic; ignore qualities for stats
    def fastq_to_fasta_iter():
        for h, s, _q in records:
            yield h, s
    return summarize_fasta(fastq_to_fasta_iter())
