# Remove guides that have an exact match elsewhere in the genome
import os
import os.path
import msgpack

APP_STATIC = "/home/joshm/GUIDES/CRISPR-Library-Designer/static"

path1 = os.path.join(APP_STATIC, 'data/GRCh37_guides_msgpack_Azimuth_domains/')

exome_path_hum = os.path.join(APP_STATIC, 'data', 'exome_hum.txt')

mer_len = 23
seq_len = mer_len - 3

exome_mers = None

def valid(guide):
  hits = 0
  for middle in ["A", "G", "T", "C"]:
    guide_seq = guide['seq'][-seq_len:] + middle + "GG"
    if guide_seq in exome_mers:
      hits += exome_mers[guide_seq]
  return hits < 2

print "preparing hum kmers"
with open(exome_path_hum, 'r') as input:
  exome = input.read()
  exome_mers = {}
  for i in range(len(exome)):
    s = exome[i:i + mer_len]
    if s in exome_mers:
      exome_mers[s] += 1
    else:
      exome_mers[s] = 1

total_begin = 0
total_after = 0
print "pruning dup hum guides"
for file in os.listdir(path1):
  with open(os.path.join(path1, file), 'r') as datafile:
    gRNAs = msgpack.load(datafile)
  total_begin += len(gRNAs)
  gRNAs = [i for i in gRNAs if valid(i)]
  total_after += len(gRNAs)

print float(total_after) / total_begin * 100, total_begin, total_after
