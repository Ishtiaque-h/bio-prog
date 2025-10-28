#!/usr/bin/env python3
'''
------------------------------------------------
 I/O files (FASTA/ FASTQ readers & writers)
------------------------------------------------
'''
import os
from pathlib import Path

def _open(filename, mode: str):
    """Open a sequence file."""
    if not os.path.exists(filename):
        print('Error - Invalid filename.')
        return None
    return open(filename, mode)



