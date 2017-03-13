### Generate exon start and stop sites
### Mouse

import pandas as pd
import pickle
import time

start_time = time.time()

# Create ENSG -> cdsStart, cdsStop mapping
# this is then used inside exon_info generator.

# Find all the CCDS associated with some ENSG.

ccds_ensg_map = {}
with open('ENSMUSG-CCDS.txt', 'r') as ensg_ccds:
  for line in ensg_ccds:
    comps = line.strip('\n').split('\t')
    ensg = comps[0]
    ccds = comps[1] + '.'
    if len(ccds) > 1:
      if ensg not in ccds_ensg_map:
        ccds_ensg_map[ensg] = [ccds]
      else:
        ccds_ensg_map[ensg].append(ccds)

# CCDS -> range mapping (str |-> list)
with open('CCDS_coords.csv', 'r') as ccds_coords_file:
  df = pd.read_csv(ccds_coords_file, sep="\t", header=0)

  for i, row in df.iterrows():
    if i == 10: break
    starts_list = [int(num) for num in row['exonStarts'].split(',')[:-1]]
    ends_list = [int(num) for num in row['exonEnds'].split(',')[:-1]]

    # Expand to include intronic sequences (5 each side)
    for k in range(len(starts_list)):
      starts_list[k] -= 5

    for k in range(len(ends_list)):
      ends_list[k] += 5

    # recombine into string
    starts_list_str = [(','.join([str(n) for n in x]) + ',') for x in starts_list]
    ends_list_str = [(','.join([str(n) for n in x]) + ',') for x in ends_list]

    # reassign
    df.set_value(i, 'exonStarts', starts_list_str)
    df.set_value(i, 'exonEnds', ends_list_str)

    # if we have ccds info...
    if row['name'] in ccds_ensg_map:
      df.set_value(i, 'name', ccds_ensg_map[row['name']])
    else:
      print "could not find in map", row['name']

  # write exon_info
  exon_info = df[["name", "chrom", "strand", "exonCount", "exonStarts", "exonEnds"]]

  with open("results/exon_info_mus.p", "wb") as f:
    pickle.dump(exon_info, f)

  # write new refGene

  df.to_csv('results/refGene_mus.txt', sep="\t", index=False, header=False)

end_time = time.time()
hours, rem = divmod(end_time-start_time, 3600)
minutes, seconds = divmod(rem, 60)
print "time elapsed"
print("{:0>2}:{:0>2}:{:05.2f}".format(int(hours),int(minutes),seconds))
