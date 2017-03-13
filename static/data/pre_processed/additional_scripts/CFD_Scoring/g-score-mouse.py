#Calculates the Cutting Frequency Determination score
#Requirements: 1. Pickle file with mismatch scores in working directory
#              2. Pickle file containing PAM scores in working directory
#Input: 1. 23mer WT sgRNA sequence
#       2. 23mer Off-target sgRNA sequence
#Output: CFD score
import pickle
import argparse
import re
import numpy as np
import os
from multiprocessing import Process, Manager
import time

t0 = time.time()

NUM_CORES = 16

# Did small experiment on Jan 11
# PAM does not matter - good!
# so we'll just use a TGG

APP_STATIC = "/home/joshm/GUIDES/CRISPR-Library-Designer/static"
exome_seq_path = os.path.join(APP_STATIC, 'data', 'GRCm38_exons')
mer_len = 20

def get_parser():
    parser = argparse.ArgumentParser(description='Calculates CFD score')
    parser.add_argument('--wt',
        type=str,
        help='WT 23mer sgRNA sequence')
    parser.add_argument('--off',
        type=str,
        help='Off-target 23mer sgRNA sequence')
    return parser

#Reverse complements a given string
def revcom(s):
    basecomp = {'A': 'T', 'C': 'G', 'G': 'C', 'T': 'A','U':'A'}
    letters = list(s[::-1])
    letters = [basecomp.get(base, 'N') for base in letters]
    return ''.join(letters)

#Unpickle mismatch scores and PAM scores
def get_mm_pam_scores():
    try:
        mm_scores = pickle.load(open('mismatch_score.pkl','rb'))
        pam_scores = pickle.load(open('pam_scores.pkl','rb'))
        return (mm_scores,pam_scores)
    except:
        raise Exception("Could not find file with mismatch scores or PAM scores")

#Calculates CFD score
def calc_cfd(wt,sg,pam):
    mm_scores,pam_scores = get_mm_pam_scores()
    score = 1
    sg = sg.replace('T','U')
    wt = wt.replace('T','U')
    s_list = list(sg)
    wt_list = list(wt)
    for i,sl in enumerate(s_list):
        if wt_list[i] == sl:
            score*=1
        else:
            key = 'r'+wt_list[i]+':d'+revcom(sl)+','+str(i+1)
            score*= mm_scores[key]
    score*=pam_scores[pam]
    return (score)

def get_pot_off_targets(seq):
    seq_list = list(seq)
    backup_seq_list = list(seq)
    nts = ['A','T','C','G']
    results = {}
    for a in range(len(seq)):
        for a_sym in nts:
            seq_list[a] = a_sym
            for b in range(a + 1, len(seq)):
                for b_sym in nts:
                    seq_list[b] = b_sym
                    for c in range(b + 1, len(seq)):
                        for c_sym in nts:
                            seq_list[c] = c_sym
                            new_seq = ''.join(seq_list)
                            results[new_seq] = True
                        seq_list[c] = backup_seq_list[c]
                seq_list[b] = backup_seq_list[b]
            seq_list[a] = backup_seq_list[a]
    if seq in results:
        del results[seq]
    return results.keys()

# parallel process container
manager = Manager()
off_target_scores = manager.dict()

# generate all 21-mers followed by GG
print "preparing mus kmers", time.time() - t0
exome_mers = {}
for file in os.listdir(exome_seq_path):
  file_loc = os.path.join(exome_seq_path, file)
  with open(file_loc, 'r') as file_data:
    fwdseq = file_data.read()
  revseq = revcom(fwdseq)

  for seq in [fwdseq, revseq]:
    for i in range(len(seq) - mer_len - 2):
      s = seq[i: i + mer_len]
      if seq[i + mer_len + 1 : i + mer_len + 3] != "GG": # only PAMs
        continue
      if 'N' in s:
        continue
      if s in off_target_scores:
        off_target_scores[s] = 'inf'
        continue
      off_target_scores[s] = 0
      if s in exome_mers:
        exome_mers[s] += 1
      else:
        exome_mers[s] = 1

# Parallelize
def process_core(pid, results, protospacers):
    i = 0
    for protospacer in protospacers:
        i += 1
        if i % 10000 == 0: print i, 'pid = ', pid
        score = 0
        off_targets = get_pot_off_targets(protospacer)
        for off_target in off_targets:
            if off_target in exome_mers:
                wt = protospacer + "CGG"
                sg = off_target
                pam = "GG"
                score += exome_mers[off_target] * calc_cfd(wt, sg, pam)
        off_target_scores[protospacer] = score

# throw onto cores
print 'throwing onto cores', time.time() - t0
processes = []
exome_mers_keys = exome_mers.keys()
unit = len(exome_mers_keys) / NUM_CORES + 1
print 'unit is', unit, time.time() - t0
for i in range(NUM_CORES):
  start = unit * i
  end = min(unit * (i + 1), len(exome_mers_keys))
  protospacers = exome_mers_keys[start:end]
  p = Process(target = process_core, args=(i, off_target_scores, protospacers,))
  processes.append(p)
for process in processes:
  process.start()
for process in processes:
  process.join()

print 'writing results', time.time() - t0
if __name__ == '__main__':
    with open("G_scores_mus.p", "wb") as output:
        pickle.dump(dict(off_target_scores), output)
print 'done', time.time() - t0
