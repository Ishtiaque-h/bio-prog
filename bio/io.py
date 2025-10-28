#!/usr/bin/env python3
'''
------------------------------------------------
 I/O files (FASTA/ FASTQ readers & writers)
------------------------------------------------
'''
def open_file(filename):
    """Open a sequence file."""
    if not os.path.exists(filename):
        print('Error - Invalid filename.')
        return None

    
    return open(filename, 'r') as f



