# Python packages
import json
import os.path
from Queue import PriorityQueue
import random
import re
import time
import pickle
import numpy as np

start_time = time.time()

# CLD code from parent directories
import os
import inspect
currentdir = os.path.dirname(os.path.abspath(inspect.getfile(inspect.currentframe())))
rootdir = os.path.dirname(os.path.dirname(os.path.dirname(currentdir)))
os.sys.path.insert(0,rootdir)
import azimuth.model_comparison


class GuideRNA():
  """Holder of gRNA information"""
  def __init__(self, selected, start, seq, PAM, score, exon_ranking, ensembl_gene, gene_name, functional_domain, has_exome_repeat, off_target_score):
    self.start = start
    self.seq = seq
    self.PAM = PAM
    self.score = score
    self.exon_ranking = exon_ranking
    self.ensembl_gene = ensembl_gene
    self.gene_name = gene_name
    self.selected = selected
    self.functional_domain = functional_domain
    if functional_domain:
      self.has_functional_domain = True
    else:
      self.has_functional_domain = False
    self.has_exome_repeat = has_exome_repeat
    self.off_target_score = off_target_score
    if off_target_score == 'inf':
      self.off_target_score = 10000

  def serialize_for_display(self):
    """Serialize for the way we are returning json"""
    return {
      "score": self.score,
      "start": self.start,
      "seq": self.seq,
      "PAM": self.PAM,
      "selected": self.selected,
      "has_exome_repeat": self.has_exome_repeat,
      "off_target_score": self.off_target_score,
      "has_functional_domain": self.has_functional_domain
    }

  def cmp_scheme(self, g):
    return (-g.off_target_score, g.has_functional_domain, g.score)

  def __cmp__(self, other):
    return cmp(self.cmp_scheme(self), self.cmp_scheme(other))

params = {
  "PAM": "NGG",
  "protospacer_len": 20,
  "prime5": True,
  "use_Doench": True,
  "quantity": 100
}

# azimuth model
print "loading azimuth models", time.time() - start_time
azimuth_saved_model_dir = os.path.join(os.path.dirname(azimuth.__file__), 'saved_models')
model_name = 'V3_model_full.pickle'
azimuth_model_file = os.path.join(azimuth_saved_model_dir, model_name)
with open(azimuth_model_file, 'rb') as f:
  azimuth_model = pickle.load(f)

azimuth_scores = {}

def get_azimuth_score(mer30):
  if mer30 in azimuth_scores:
    return azimuth_scores[mer30]
  else:
    score = azimuth.model_comparison.predict(np.array([mer30]), aa_cut=None, percent_peptide=None, model=azimuth_model, model_file=azimuth_model_file)[0]
    print "generating Azimuth", mer30, score
    azimuth_scores[mer30] = score
    return score

modPAM = params["PAM"].upper()
modPAM = modPAM.replace('N', '[ATCG]')
params["modPAM"] = modPAM
params["PAM_len"] = len(params["PAM"])

revcompl = lambda x: ''.join([{'A':'T','C':'G','G':'C','T':'A','N':'N'}[B] for B in x][::-1])

def gene_exon_file(gene, exon):
  filename = gene + "_" + str(exon)
  seq_path = os.path.join('../GRCh37_exons/', filename)
  if os.path.isfile(seq_path):
    with open(seq_path) as infile:
      return infile.read()
  else:
    return None

def exome(gene):
  seq = ""
  exon = 0
  exon_seq = gene_exon_file(gene, exon)
  while exon_seq:
    seq += exon_seq
    exon += 1
    exon_seq = gene_exon_file(gene, exon)
  return seq

def random_guides(gene, library_size):
  seq = exome(gene["ensembl_id"])

  q = PriorityQueue()

  def process_guide(m, selected, max_queue_size, seq):
    PAM_start = m.start()
    # Doench score requires the 4 before and 6 after 20-mer (gives 30-mer)
    mer30 = seq[PAM_start-params["protospacer_len"]-4:PAM_start+params["PAM_len"]+3]
    if len(mer30) != 30:
      print "Error! The following guide is not long enough:", seq, mer30, gene["ensembl_id"], gene["name"]
    score = get_azimuth_score(mer30)
    protospacer = ""
    PAM = ""
    if params["prime5"]:
      protospacer = seq[PAM_start-params["protospacer_len"]:PAM_start]
      PAM = seq[PAM_start:PAM_start+params["PAM_len"]]
    else:
      protospacer = seq[PAM_start+params["PAM_len"]:PAM_start+params["PAM_len"]+params["protospacer_len"]]
      PAM = seq[PAM_start:PAM_start+params["PAM_len"]]
    potential_gRNA = GuideRNA(selected, PAM_start-params["protospacer_len"], protospacer, PAM, score, -1, gene["ensembl_id"], gene["name"])

    # If there's enough room, add it, no question.
    if q.qsize() < max_queue_size:
      q.put(potential_gRNA)

    # We remove the other case, since we want to add anything that is given to this function.
    #### Otherwise, take higher score
    ###else:
    ###  lowest_gRNA = q.get()
    ###  if potential_gRNA.score > lowest_gRNA.score:
    ###    q.put(potential_gRNA)
    ###  else:
    ###    q.put(lowest_gRNA)

  # Logic continues here, outside process_guide
  seq_rc = revcompl(seq)
  forward_matches = [('Forward', m) for m in re.finditer(params["modPAM"], seq)]
  reverse_matches = [('Reverse', m) for m in re.finditer(params["modPAM"], seq_rc)]
  all_matches = forward_matches + reverse_matches
  random.shuffle(all_matches)

  # generate library_size number of guides
  required = library_size
  for match in all_matches:
    if required == 0:
      break
    direction, m = match
    if params["prime5"] and (m.start() < params["protospacer_len"] + 4 or m.start() + params["PAM_len"] + 3 > len(seq)):
      continue
    elif not params["prime5"] and (m.start() + params["PAM_len"] + params["protospacer_len"] > len(seq)):
      continue

    if direction == 'Forward':
      process_guide(m, True, params["quantity"], seq)
    else:
      process_guide(m, True, params["quantity"], seq_rc)
    required -= 1

  # Pop gRNAs into our 'permament' storage, i.e. return them.
  gRNAs = []
  while not q.empty():
    gRNA = q.get()
    gRNAs.append(gRNA.serialize_for_display())
  return gRNAs
