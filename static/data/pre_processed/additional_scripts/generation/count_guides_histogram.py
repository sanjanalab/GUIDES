# Remove guides that have an exact match elsewhere in the genome
import os
import os.path
import msgpack

APP_STATIC = "/home/joshm/GUIDES/CRISPR-Library-Designer/static"

path1 = os.path.join(APP_STATIC, 'data/GRCh37_guides_msgpack_Azimuth_domains/')

exome_path_hum = os.path.join(APP_STATIC, 'data', 'exome_hum.txt')

mer_len = 10
seq_len = mer_len - 3

exome_mers = None

with open(exome_path_hum, 'r') as input:
  exome = input.read()
  exome_mers = {}
  for i in range(len(exome)):
    s = exome[i:i + mer_len]
    if exome[i + mer_len + 1 : i + mer_len + 3] != "GG":
      continue 
    if s in exome_mers:
      exome_mers[s] += 1
    else:
      exome_mers[s] = 1

multiple_matches = 0
for mer in exome_mers:
  if exome_mers[mer] > 1:
    multiple_matches += 1

print len(exome_mers), multiple_matches
