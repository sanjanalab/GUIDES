# Remove homopolymer (>4 A/C/Gs in a row) and U6 terminator (>3 Ts in a row)
import os
import os.path
import msgpack

APP_STATIC = "/home/joshm/GUIDES/CRISPR-Library-Designer/static"

path1 = os.path.join(APP_STATIC, 'data/GRCh37_guides_msgpack_Azimuth_domains/')
path2 = os.path.join(APP_STATIC, 'data/GRCm38_guides_msgpack_Azimuth_domains/')

def valid(guide):
  return not (('AAAAA' in guide['seq']) or ('CCCCC' in guide['seq']) or ('GGGGG' in guide['seq']) or ('TTTT' in guide['seq']))

for file in os.listdir(path1):
  with open(os.path.join(path1, file), 'r') as datafile:
    print os.path.join(path1, file)
    gRNAs = msgpack.load(datafile)
  gRNAs = [i for i in gRNAs if valid(i)]
  with open(os.path.join(path1, file), 'w+') as datafile:
    msgpack.dump(gRNAs, datafile)

for file in os.listdir(path2):
  with open(os.path.join(path2, file), 'r') as datafile:
    gRNAs = msgpack.load(datafile)
  gRNAs = [i for i in gRNAs if valid(i)]
  with open(os.path.join(path2, file), 'w+') as datafile:
    msgpack.dump(gRNAs, datafile)
