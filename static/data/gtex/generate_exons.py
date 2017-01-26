### Generate exon start and stop sites

import pandas as pd
import pickle

# Create ENSG -> cdsStart, cdsStop mapping
# this is then used inside exon_info generator.

# CCDS -> start stop mapping
ccds_coords_map = {}
with open('../pre_processed/CDS/ccds_coords_hum.txt', 'r') as ccds_coords:
  for line in ccds_coords:
    comps = line.strip('\n').split('\t')
    ccds = comps[0].split('.')[0]
    start = int(comps[1])
    stop = int(comps[2])
    ccds_coords_map[ccds] = (min(start, stop), max(start, stop))

# ENSG -> CCDS mapping
ensg_ccds_map = {}
with open('../pre_processed/CDS/ENSG-CCDS_hum.txt', 'r') as ensg_ccds:
  for line in ensg_ccds:
    comps = line.strip('\n').split('\t')
    ensg = comps[0] + '.' + comps[1]
    ccds = comps[2]
    if len(ccds) > 0:
      if ccds in ccds_coords_map:
        start, stop = ccds_coords_map[ccds]
        if ensg not in ensg_ccds_map:
          ensg_ccds_map[ensg] = (start, stop)
        else:
          prev_stop, prev_start = ensg_ccds_map[ensg]
          ensg_ccds_map[ensg] = (min(prev_start, start), max(prev_stop, stop))

if __name__ == "__main__":
  filename = "refGene_base.txt"
  df = pd.read_csv(filename, sep="\t", header=None)
  df.columns=['','name','chrom','strand','txStart','txEnd','cdsStart','cdsEnd','exonCount','exonStarts','exonEnds','id','name2','cdsStartStat','cdsEndStat','exonFrames']
  # process the dataframe
  for i, row in df.iterrows():
    starts_list = [int(num) for num in row['exonStarts'].split(',')[:-1]]
    ends_list = [int(num) for num in row['exonEnds'].split(',')[:-1]]
    df.set_value(i, 'exonStarts', [0,0,0] + starts_list) # hacky - will force an error if 'exonStarts' never changed below.

    # if we have ccds info...
    if row['name'] in ensg_ccds_map:
      starts_list_processed = []
      ends_list_processed = []

      cds_start, cds_stop = ensg_ccds_map[row['name']]

      assert(len(starts_list) == len(ends_list))
      for j in range(len(starts_list)):
        start, end = starts_list[j], ends_list[j]
        if end < cds_start or start > cds_stop: # whole exon outside cds
          starts_list_processed.append(start)
          ends_list_processed.append(start + 1)
          print "removed exon " + str(j) + " of " + str(len(starts_list))
          continue # don't add this exon
        starts_list_processed.append(max(start, cds_start))
        ends_list_processed.append(min(end, cds_stop))
      df.set_value(i, 'exonStarts', list(starts_list_processed))
      df.set_value(i, 'exonEnds', list(ends_list_processed))
      df.set_value(i, 'cdsStart', cds_start)
      df.set_value(i, 'cdsEnd', cds_stop)
    else: # we dont' have ccds... keep default
      df.set_value(i, 'exonStarts', starts_list)
      df.set_value(i, 'exonEnds', ends_list)

  # write exon_info
  exon_info = df[["name", "chrom", "strand", "exonCount", "exonStarts", "exonEnds"]]

  with open("../pre_processed/exon_info.p", "wb") as f:
    pickle.dump(exon_info, f)

  # write new refGene
  df['exonStarts'] = df.apply(lambda x: (','.join([str(n) for n in x['exonStarts']]) + ','), axis=1)
  df['exonEnds'] = df.apply(lambda x: (','.join([str(n) for n in x['exonEnds']]) + ','), axis=1)

  df.to_csv('refGene.txt', sep="\t", index=False, header=False)
