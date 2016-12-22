# Create ENSG -> cdsStart, cdsStop mapping
# this is then used inside exon_info generator.

# CCDS -> start stop mapping
ccds_coords_map = {}
with open('ccds_coords_hum.txt', 'r') as ccds_coords:
  for line in ccds_coords:
    comps = line.split('\t')
    ccds = comps[0].split('.')[0]
    start = comps[1]
    stop = comps[2]
    ccds_coords_map[ccds] = (start, stop)

# ENSG -> CCDS mapping
ensg_ccds_map = {}
with open('ENSG-CCDS_hum.txt', 'r') as ensg_ccds:
  for line in ensg_ccds:
    comps = line.split('\t')
    ensg = comps[0] + '.' + comps[1]
    ccds = comps[2]
    if len(ccds) > 0:
      ensg_ccds_map[ensg] = ccds


