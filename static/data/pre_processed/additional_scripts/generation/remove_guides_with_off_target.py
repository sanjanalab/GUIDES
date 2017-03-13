# Remove guides that have an exact match elsewhere in the genome
import os
import os.path
import msgpack

APP_STATIC = "/home/joshm/GUIDES/CRISPR-Library-Designer/static"

path1 = os.path.join(APP_STATIC, 'data/GRCh37_guides_msgpack_Azimuth_domains/')
path2 = os.path.join(APP_STATIC, 'data/GRCm38_guides_msgpack_Azimuth_domains/')

def valid(guide, species, file_to_skip = None):
  if species == 'hum':
    exon_path = APP_STATIC + "/data/GRCh37_exons"
  elif species == 'mus':
    exon_path = APP_STATIC + "/data/GRCm38_exons"

  valid = True
  for filename in os.listdir(exon_path):
    if filename == file_to_skip:
      continue
    with open(os.path.join(exon_path, filename), 'r') as exon_file:
      exon_seq = exon_file.read()
      if guide['seq'] in exon_seq:
        valid = False
        break
  return valid

for file in os.listdir(path1):
  with open(os.path.join(path1, file), 'r') as datafile:
    print os.path.join(path1, file)
    gRNAs = msgpack.load(datafile)
  file_to_skip = file[:-2]
  gRNAs = [i for i in gRNAs if valid(i, 'hum', file_to_skip)]
  with open(os.path.join(path1, file), 'w+') as datafile:
    msgpack.dump(gRNAs, datafile)

for file in os.listdir(path2):
  with open(os.path.join(path2, file), 'r') as datafile:
    gRNAs = msgpack.load(datafile)
  file_to_skip = file[:-2]
  gRNAs = [i for i in gRNAs if valid(i, 'mus', file_to_skip)]
  with open(os.path.join(path2, file), 'w+') as datafile:
    msgpack.dump(gRNAs, datafile)
