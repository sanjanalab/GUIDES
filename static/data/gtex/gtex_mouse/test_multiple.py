# CCDS -> ENSG mapping --> gives a unique ENSMUSG for each CCDS.

ccds_ensg_map = {}
with open('ENSMUSG-CCDS.txt', 'r') as ensg_ccds:
  for line in ensg_ccds:
    comps = line.strip('\n').split('\t')
    ensg = comps[0]
    ccds = comps[1] + '.'
    if len(ccds) > 1:
      # check for conflicts
      if ccds in ccds_ensg_map:
        if ccds_ensg_map[ccds] != ensg:
          print "conflicting CCDS entries", ccds_ensg_map[ccds], ensg

      # assign
      ccds_ensg_map[ccds] = ensg

ccds_ensg_map = {}
with open('ENSMUSG-CCDS.txt', 'r') as ensg_ccds:
  for line in ensg_ccds:
    comps = line.strip('\n').split('\t')
    ensg = comps[0]
    ccds = comps[1] + '.'
    if len(ccds) > 1:
      # check for conflicts
      if ensg in ccds_ensg_map:
        if ccds_ensg_map[ensg] != ccds:
          print "conflicting CCDS entries", ccds_ensg_map[ensg], ccds

      # assign
      ccds_ensg_map[ensg] = ccds
