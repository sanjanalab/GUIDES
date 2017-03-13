# Remove guides that have an exact match elsewhere in the genome
import os
import os.path
import msgpack

APP_STATIC = "/home/joshm/GUIDES/CRISPR-Library-Designer/static"

path1 = os.path.join(APP_STATIC, 'data/GRCh37_guides_msgpack_Azimuth/')
path2 = os.path.join(APP_STATIC, 'data/GRCm38_guides_msgpack_Azimuth/')

exome_path_hum = os.path.join(APP_STATIC, 'data', 'exome_hum.txt')
exome_path_mus = os.path.join(APP_STATIC, 'data', 'exome_mus.txt')

mer_len = 13

exome_mers = None

def valid(guide):
  hits = 0
  for middle in ["A", "G", "T", "C"]:
    guide_seq = guide['seq'] + middle + "GG"
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

print "pruning dup hum guides"
for file in os.listdir(path1):
  with open(os.path.join(path1, file), 'r') as datafile:
    gRNAs = msgpack.load(datafile)
  gRNAs = [i for i in gRNAs if valid(i)]
  with open(os.path.join(path1, file), 'w+') as datafile:
    msgpack.dump(gRNAs, datafile)

print "preparing mus kmers"
with open(exome_path_mus, 'r') as input:
  exome = input.read()
  exome_mers = {}
  for i in range(len(exome)):
    s = exome[i:i + mer_len]
    if s in exome_mers:
      exome_mers[s] += 1
    else:
      exome_mers[s] = 1

print "pruning dup mus guides"
for file in os.listdir(path2):
  with open(os.path.join(path2, file), 'r') as datafile:
    gRNAs = msgpack.load(datafile)
  gRNAs = [i for i in gRNAs if valid(i)]
  with open(os.path.join(path2, file), 'w+') as datafile:
    msgpack.dump(gRNAs, datafile)
