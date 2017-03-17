### Generate exon start and stop sites

import pandas as pd
import pickle
import time

t0 = time.time()

# Create ENSG -> cdsStart, cdsStop mapping
# this is then used inside exon_info generator.

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

def restrict_to_intervals(start, stop, intervals):
  for interval in intervals:
    if interval[0] >= start and interval[1] <= stop: # refGene should contain the CCDS. CCDS inside RefGene
      return interval
  return (start, start + 1)

print 'opening CCDS_coords.csv', time.time() - t0
# CCDS -> range mapping (str |-> list)
with open('../pre_processed/CDS/CCDS_coords.csv', 'r') as ccds_coords_file:
  ccds_df = pd.read_csv(ccds_coords_file, sep="\t", header=0)

# ENSG -> CCDS mapping --> list of intervals [(a,b), (c,d), ...]
print 'opening ENSG-CCDS_hum.txt', time.time() - t0
ensg_ccds_map = {}
with open('../pre_processed/CDS/ENSG-CCDS_hum.txt', 'r') as ensg_ccds:
  for line in ensg_ccds:
    comps = line.strip('\n').split('\t')
    ensg = comps[0] + '.' + comps[1]
    ccds = comps[2] + '.'
    if ensg == "ENSG00000137185.7":
      print ensg, ccds
    if len(ccds) > 0:
      starts = ccds_df[ccds_df['name'].str.startswith(ccds)]['exonStarts']
      stops = ccds_df[ccds_df['name'].str.startswith(ccds)]['exonEnds']
      if len(starts) > 0 and len(stops) > 0:
        starts = starts.tolist()[0]
        stops = stops.tolist()[0]
        starts_nums = [int(num) for num in starts.split(',')[:-1]]
        stops_nums = [int(num) for num in stops.split(',')[:-1]]

        if ensg not in ensg_ccds_map:
          ensg_ccds_map[ensg] = union_intervals([(int(a),int(b)) for a, b in zip(starts_nums, stops_nums)])
        else:
          new_intervals = [(int(a),int(b)) for a, b in zip(starts_nums, stops_nums)]
          all_intervals = ensg_ccds_map[ensg] + new_intervals
          ensg_ccds_map[ensg] = union_intervals(all_intervals)

print 'ensg_ccds_map["ENSG00000137185.7"]'
print ensg_ccds_map["ENSG00000137185.7"]

print 'starting main script', time.time() - t0
if __name__ == "__main__":
  filename = "refGene_base.txt"
  df = pd.read_csv(filename, sep="\t", header=None)
  df.columns=['','name','chrom','strand','txStart','txEnd','cdsStart','cdsEnd','exonCount','exonStarts','exonEnds','id','name2','cdsStartStat','cdsEndStat','exonFrames']
  # process the dataframe
  print len(df)
  results_df = pd.DataFrame(columns=['','name','chrom','strand','txStart','txEnd','cdsStart','cdsEnd','exonCount','exonStarts','exonEnds','id','name2','cdsStartStat','cdsEndStat','exonFrames'])
  res_idx = 0
  for i, row in df.iterrows():
    if i % len(df) / 1000 == 0: print i
    starts_list = [int(num) for num in row['exonStarts'].split(',')[:-1]]
    ends_list = [int(num) for num in row['exonEnds'].split(',')[:-1]]
    df.set_value(i, 'exonStarts', [0,0,0] + starts_list) # hacky - will force an error if 'exonStarts' never changed below.

    # if we have ccds info...
    if row['name'] in ensg_ccds_map:
      starts_list_processed = []
      ends_list_processed = []

      cds_intervals = ensg_ccds_map[row['name']]

      assert(len(starts_list) == len(ends_list))
      for j in range(len(starts_list)):
        start, end = starts_list[j], ends_list[j]
        new_interval = restrict_to_intervals(start, end, cds_intervals)

        if new_interval[1] == new_interval[0] + 1:
          print "removed exon " + str(j) + " of " + str(len(starts_list))

        starts_list_processed.append(new_interval[0])
        ends_list_processed.append(new_interval[1])

      # Expand to include intronic sequences (5 each side)
      for k in range(len(starts_list_processed)):
        starts_list_processed[k] -= 5

      for k in range(len(ends_list_processed)):
        ends_list_processed[k] += 5

      df.set_value(i, 'exonStarts', list(starts_list_processed))
      df.set_value(i, 'exonEnds', list(ends_list_processed))
      df.set_value(i, 'cdsStart', cds_intervals[0][0])
      df.set_value(i, 'cdsEnd', cds_intervals[-1][1])

      new_row = df.iloc[i]
      results_df.loc[res_idx] = new_row
      res_idx += 1
      
    #else: # we dont' have ccds... keep default
      # Expand to include intronic sequences (5 each side)
      #for k in range(len(starts_list)):
      #  starts_list[k] -= 5

      #for k in range(len(ends_list)):
      #  ends_list[k] += 5

      #df.set_value(i, 'exonStarts', starts_list)
      #df.set_value(i, 'exonEnds', ends_list)

  # write exon_info
  exon_info = results_df[["name", "chrom", "strand", "exonCount", "exonStarts", "exonEnds"]]

  with open("../pre_processed/exon_info.p", "wb") as f:
    pickle.dump(exon_info, f)

  # write new refGene
  results_df['exonStarts'] = results_df.apply(lambda x: (','.join([str(n) for n in x['exonStarts']]) + ','), axis=1)
  results_df['exonEnds'] = results_df.apply(lambda x: (','.join([str(n) for n in x['exonEnds']]) + ','), axis=1)

  results_df.to_csv('refGene.txt', sep="\t", index=False, header=False)

end_time = time.time()
hours, rem = divmod(end_time-t0, 3600)
minutes, seconds = divmod(rem, 60)
print "time elapsed"
print("{:0>2}:{:0>2}:{:05.2f}".format(int(hours),int(minutes),seconds))
