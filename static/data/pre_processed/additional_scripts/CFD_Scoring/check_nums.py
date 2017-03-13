import os

APP_STATIC = "/home/joshm/GUIDES/CRISPR-Library-Designer/static"
exome_path_hum = os.path.join(APP_STATIC, 'data', 'exome_hum.txt')
mer_len = 20

off_target_scores_count = 0
total_num_spacers = 0

off_target_scores = {}

print "preparing hum kmers"
with open(exome_path_hum, 'r') as input:
  exome = input.read()
  exome_mers = {}
  for i in range(len(exome) - mer_len - 3):
    s = exome[i:i + mer_len]
    if exome[i + mer_len + 1 : i + mer_len + 3] != "GG": # only PAMs
      continue
    if '=' in s or 'N' in s:
      continue
    total_num_spacers += 1
    if s in off_target_scores:
      off_target_scores_count -= 1
      off_target_scores[s] = 'inf'
      continue
    off_target_scores[s] = 0
    off_target_scores_count += 1
    if s in exome_mers:
      exome_mers[s] += 1
    else:
      exome_mers[s] = 1

print 'len(exome_mers)', len(exome_mers)
print 'off_target_scores_count', off_target_scores_count
print 'total_num_spacers', total_num_spacers
