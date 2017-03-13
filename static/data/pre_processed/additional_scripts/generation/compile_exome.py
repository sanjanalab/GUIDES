# pool exome into single file
import os
import os.path

APP_STATIC = "/home/joshm/GUIDES/CRISPR-Library-Designer/static"

exons_path_hum = os.path.join(APP_STATIC, 'data/GRCh37_exons')
exons_path_mus = os.path.join(APP_STATIC, 'data/GRCm38_exons')

output_path_hum = os.path.join(APP_STATIC, 'data', 'exome_hum.txt')
output_path_mus = os.path.join(APP_STATIC, 'data', 'exome_mus.txt')

revcompl = lambda x: ''.join([{'N': 'N','A':'T','C':'G','G':'C','T':'A'}[B] for B in x][::-1])

def pool(exon_path, filename_save):
  full_seq = ""
  for filename in os.listdir(exon_path):
    with open(os.path.join(exon_path, filename), 'r') as exon_file:
      seq = exon_file.read()
      full_seq += seq
      full_seq += "=" * 20
      try:
        full_seq += revcompl(seq)
        full_seq += "=" * 20
      except KeyError:
         print seq
  with open(filename_save, 'w') as output:
    output.write(full_seq)

pool(exons_path_hum, output_path_hum)
#pool(exons_path_mus, output_path_mus)
