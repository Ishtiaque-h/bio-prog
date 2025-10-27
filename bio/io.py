#!/usr/bin/env python3
'''
------------------------------------------------
 I/O files (FASTA/ FASTQ readers & writers)
------------------------------------------------
'''
import os

filename = input('Enter the path to the sequence file: ')
if not os.path.exists(filename):
  print('Error - Invalid filename.')

if filename.lower().endswith(('fasta','fa','fna','txt')):
  sequences = []
  with open(filename,'r') as f:
    header = None
    sequence = ''
    for line in f:
      line =line.strip()
      if line.startswith('>'):
        header = line[1:].split()[0]
        seq = []
      else:
        seq.append(line)
elif filename.lower().endswith(('fastq','fq')):
  pass
else:
  print('Error - Wrong File Format. Only fasta and fastq formats are allowed.')



