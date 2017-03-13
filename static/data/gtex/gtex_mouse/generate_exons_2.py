'''
1. Generate CCDS -> ENSG mapping.
2. Get ENSG -> CCDS intervals mapping from CCDS_coords.
3. Go through CCDS_coords.csv and change all CCDS to ENSG.
4. Remove duplicates on name column.
5. Change ENSG exon positions to the intervals we found above.
6. Add 5 to either side.
7. Save
'''

import pandas as pd
import pickle
import time

t0 = time.time()

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

# 1. Generate CCDS -> ENSG mapping.
ccds_ensg_map = {}
with open('ENSMUSG-CCDS.txt', 'r') as ensg_ccds:
  for line in ensg_ccds:
    comps = line.strip('\n').split('\t')
    ensg = comps[0]
    ccds = comps[1]
    if len(ccds) > 0:
      if ccds not in ccds_ensg_map:
        ccds_ensg_map[ccds] = ensg

# 2. Get ENSG -> CCDS intervals mapping from CCDS_coords.

# set of starts/stops for an ensg
ensg_coords = {}

# CCDS -> range mapping (str |-> list)
with open('CCDS_coords.csv', 'r') as ccds_coords_file:
  ccds_df = pd.read_csv(ccds_coords_file, sep="\t", header=0)

  # go through each CCDS entry individually
  print "going through each CCDS entry for ensg", time.time() - t0
  for i, row in ccds_df.iterrows():

    if i % (len(ccds_df) / 100) == 0: print i, '/', len(ccds_df), time.time() - t0

    # find which ensembl gene it belongs to
    ccds_name = row['name'].split('.')[0]
    if ccds_name not in ccds_ensg_map:
      print ccds_name, "not in ccds_ensg_map"
      continue
    ensg = ccds_ensg_map[ccds_name]

    # Add the associated information to ensg_coords
    starts = row['exonStarts']
    stops = row['exonEnds']

    starts_nums = [int(num) for num in starts.split(',')[:-1]]
    stops_nums = [int(num) for num in stops.split(',')[:-1]]

    new_intervals = [(int(a),int(b)) for a, b in zip(starts_nums, stops_nums)]

    if ensg not in ensg_coords:
      ensg_coords[ensg] = {
        'df_data': row,
        'intervals': union_intervals(new_intervals)
      }
    else:
      all_intervals = ensg_coords[ensg]['intervals'] + new_intervals
      ensg_coords[ensg]['intervals'] = union_intervals(all_intervals)

  print "ready to go through results and move intervals into df", time.time() - t0
  # go through results and move intervals into df
  results_df = pd.DataFrame(columns=['','name','chrom','strand','txStart','txEnd','cdsStart','cdsEnd','exonCount','exonStarts','exonEnds','id','name2','cdsStartStat','cdsEndStat','exonFrames'])
  exon_info_df = pd.DataFrame(columns=['','name','chrom','strand','txStart','txEnd','cdsStart','cdsEnd','exonCount','exonStarts','exonEnds','id','name2','cdsStartStat','cdsEndStat','exonFrames'])
  for idx, ensg in enumerate(ensg_coords):

    if idx % (len(ensg_coords) / 100) == 0: print idx, '/', len(ensg_coords), time.time() - t0

    exonCount = len(ensg_coords[ensg]['intervals'])
    if exonCount == 0:
      continue

    # expand intervals to include intronic sequences (5 each side)
    starts_list = []
    ends_list = []
    for i in range(len(ensg_coords[ensg]['intervals'])):
      starts_list.append(ensg_coords[ensg]['intervals'][i][0] - 5)
      ends_list.append(ensg_coords[ensg]['intervals'][i][1] + 5)

    # recombine into string
    starts_list_str = (','.join([str(n) for n in starts_list])) + ','
    ends_list_str = (','.join([str(n) for n in ends_list])) + ','

    new_row = ensg_coords[ensg]['df_data']
    ccdsStart = ensg_coords[ensg]['intervals'][0][0]
    ccdsEnd = ensg_coords[ensg]['intervals'][-1][1]

    new_row['name'] = ensg
    new_row['txStart'] = ccdsStart
    new_row['cdsStart'] = ccdsStart
    new_row['txEnd'] = ccdsEnd
    new_row['cdsEnd'] = ccdsEnd
    new_row['exonCount'] = exonCount
    new_row['exonStarts'] = starts_list_str
    new_row['exonEnds'] = ends_list_str
    new_row['chrom'] = new_row['chrom'].split('chr')[1]

    results_df.loc[idx] = new_row

    new_row['exonStarts'] = starts_list
    new_row['exonEnds'] = ends_list

    exon_info_df.loc[idx] = new_row

print "writing results", time.time() - t0
# Write results
exon_info = exon_info_df[["name", "chrom", "strand", "exonCount", "exonStarts", "exonEnds"]]

with open("exon_info.p", "wb") as f:
  pickle.dump(exon_info, f)

results_df.to_csv('refGene_mouse.txt', sep="\t", index=False, header=False)

end_time = time.time()
hours, rem = divmod(end_time-t0, 3600)
minutes, seconds = divmod(rem, 60)
print "time elapsed"
print("{:0>2}:{:0>2}:{:05.2f}".format(int(hours),int(minutes),seconds))
