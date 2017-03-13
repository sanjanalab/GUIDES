# pool exome into single file
import os
import pandas as pd
import gzip
from Bio import SeqIO

APP_STATIC = "/home/joshm/GUIDES/CRISPR-Library-Designer/static"

ccds_coords = os.path.join(APP_STATIC, 'data/pre_processed/CDS', 'CCDS_coords.csv')

chrom_seq_base = os.path.join(APP_STATIC, "data/GRCh37/Homo_sapiens.GRCh37.75.dna.chromosome.{}.fa.gz")

output_path = os.path.join(APP_STATIC, "data", "exome_hum_ccds.txt")

exome_output_str = ""

revcompl = lambda x: ''.join([{'N': 'N','A':'T','C':'G','G':'C','T':'A'}[B] for B in x][::-1])

def union_intervals(intervals):
  sorted_by_lower_bound = sorted(intervals, key=lambda tup: tup[0])
  merged = []
  for higher in sorted_by_lower_bound:
      if not merged:
          merged.append(higher)
      else:
          lower = merged[-1]
          # test for intersection between lower and higher:
          # we know via sorting that lower[0] <= higher[0]
          if higher[0] <= lower[1]:
              upper_bound = max(lower[1], higher[1])
              merged[-1] = (lower[0], upper_bound)  # replace by merged interval
          else:
              merged.append(higher)
  return merged

with open(ccds_coords, 'r') as ccds_coords_file:
  ccds_df = pd.read_csv(ccds_coords_file, sep="\t", header=0)

chroms = [str(i) for i in range(1, 23)]
chroms += ['X', 'Y']

total_len = 0

for chrom_sym in chroms:
  chrom = 'chr' + chrom_sym
  starts = ccds_df[ccds_df['chrom'] == chrom]['exonStarts']
  stops = ccds_df[ccds_df['chrom'] == chrom]['exonEnds']
  starts = ''.join(starts.tolist())
  stops = ''.join(stops.tolist())
  starts_nums = [int(num) for num in starts.split(',')[:-1]]
  stops_nums = [int(num) for num in stops.split(',')[:-1]]
  all_coords = union_intervals([(int(a),int(b)) for a, b in zip(starts_nums, stops_nums)])

  # get full chr sequence
  full_chr_seq_filename = chrom_seq_base.format(chrom_sym)
  handle = gzip.open(full_chr_seq_filename)
  seq = str(SeqIO.read(handle, "fasta").seq).upper()

  # add DNA between the coordinates to our output
  for coord_start, coord_stop in all_coords:
    to_output = seq[coord_start:coord_stop]
    total_len += 2 * len(to_output)
    exome_output_str += to_output
    exome_output_str += '=' * 20
    exome_output_str += revcompl(to_output)
    exome_output_str += '=' * 20

# save the new exome
with open(output_path, 'w') as output:
  output.write(exome_output_str)
