'''
Requirements:
../CFD_Scoring/off_target_scores.p
tf_guides.csv
../GUIDES/CRISPR-Library-Designer/static/data/GRCh37_exons/

Writes results:
all_scores.pdf
selected_scores.pdf
'''

import numpy as np
import random
from matplotlib import pyplot as plt
import pickle
import glob
import csv

# parameters
mer_len = 20

# results
all_scores = []
selected_scores = []

# Load off-target scores
with open('../CFD_Scoring/G_scores_hum.p', 'rb') as inp:
  off_target_scores = pickle.load(inp)

# Grab list of ensgs
ensg_d = {}
with open('tf_guides_hum.csv', 'r') as tf_guides_file:
  reader = csv.DictReader(tf_guides_file)
  for row in reader:
    ensg_id = row['Ensembl ID']#.rstrip('\r\n')
    ensg_d[ensg_id] = True

ensg_list = ensg_d.keys()

# all scores
print 'grabbing all scores'
could_not_find_count = 0
for ensg_id in ensg_list:
  form = '../GUIDES/CRISPR-Library-Designer/static/data/GRCh37_exons/' + ensg_id + '_*'
  for file_name in glob.glob(form):
    with open(file_name, 'r') as seq_file:
      seq = seq_file.read()
      for i in range(len(seq) - mer_len - 3):
        s = seq[i:i + mer_len] # s is protospacer
        if seq[i + mer_len + 1 : i + mer_len + 3] != "GG": # only PAMs
          continue
        if '=' in s or 'N' in s:
          continue
        # I'm just going to assume that s is in off_target_scores...
        if s in off_target_scores:
          all_scores.append(off_target_scores[s])
        else:
          could_not_find_count += 1

print "{0} guides did not appear in off_target_scores".format(could_not_find_count)

# selected scores
print 'grabbing selected scores'
with open('tf_guides_hum.csv', 'r') as tf_guides_file:
  reader = csv.DictReader(tf_guides_file)
  for row in reader:
    ots = float(row['Off-target score'])
    selected_scores.append(ots)

# find max
print 'fixing maxes'
max_all_scores = 0
max_selected_scores = 0

for i in range(len(all_scores)):
  if all_scores[i] != 100000:
    max_all_scores = max(max_all_scores, all_scores[i])

for i in range(len(selected_scores)):
  if selected_scores[i] != 100000:
    max_selected_scores = max(max_selected_scores, selected_scores[i])

# Change all occurences of 100000
for i in range(len(all_scores)):
  if all_scores[i] == 100000:
    all_scores[i] = max_all_scores * 2

for i in range(len(selected_scores)):
  if selected_scores[i] == 100000:
    selected_scores[i] = max_selected_scores * 2

# Pickle results
with open('all_scores_hum.p', 'wb') as all_scores_f:
  pickle.dump(all_scores, all_scores_f)

with open('selected_scores_hum.p', 'wb') as selected_scores_f:
  pickle.dump(selected_scores, selected_scores_f)
