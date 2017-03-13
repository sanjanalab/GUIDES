# pool exome into single file
import os
import os.path

APP_STATIC = "/home/joshm/GUIDES/CRISPR-Library-Designer/static"

exons_path_hum = os.path.join(APP_STATIC, 'data/GRCh37_exons')

revcompl = lambda x: ''.join([{'N': 'N','A':'T','C':'G','G':'C','T':'A'}.get(B,'X') for B in x][::-1])

to_check = 'TTCGACCCCGACCTGCCAGG'

def pool(exon_path):
  for filename in os.listdir(exon_path):
    with open(os.path.join(exon_path, filename), 'r') as exon_file:
      seq = exon_file.read()
      rev_seq = revcompl(seq)
      if to_check in seq:
        print filename, 'fwd'
      if to_check in rev_seq:
        print filename, 'rev'

pool(exons_path_hum)
