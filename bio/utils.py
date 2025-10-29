#!/usr/bin/env python3
'''
------------------------------------------------
 Helper/ utility functions
------------------------------------------------
'''
#---------------imports------------------------------
from .io import _open
from typing import Tuple

#----------------gc and atgc estimator-------------------------
def gc_and_counts(seq: str) -> Tuple[int, int, int]:
    """
    Return (gc_count, atgc_count, n_count) for a sequence.
    Ns are counted separately; atgc_count excludes Ns.
    """
    s = seq.upper()
    g = s.count("G"); c = s.count("C")
    a = s.count("A"); t = s.count("T")
    n = s.count("N")
    gc = g + c
    atgc_total = a + t + g + c
    return gc, atgc_total, n
