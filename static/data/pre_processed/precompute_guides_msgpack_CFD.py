import msgpack
import json
import pickle
import os.path
from Queue import PriorityQueue
import re
import doench_score
import azimuth.model_comparison
import numpy as np
import pandas as pd
import csv
from intervaltree import IntervalTree
from multiprocessing import Process
import os
import time

start_time = time.time()

#Reverse complements a given string
def revcom(s):
    basecomp = {'A': 'T', 'C': 'G', 'G': 'C', 'T': 'A','U':'A', 'N':'N'}
    letters = list(s[::-1])
    letters = [basecomp[base] for base in letters]
    return ''.join(letters)

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
    serialization = {
      "score": self.score,
      "start": self.start,
      "seq": self.seq,
      "PAM": self.PAM,
      "selected": self.selected,
      "has_exome_repeat": self.has_exome_repeat,
      "off_target_score": self.off_target_score,
      "has_functional_domain": self.has_functional_domain
    }
    if self.functional_domain != None:
      serialization["functional_domain"] = self.functional_domain
    return serialization

  def cmp_scheme(self, g):
    return (-g.off_target_score, g.has_functional_domain, g.score)

  def __cmp__(self, other):
    return cmp(self.cmp_scheme(self), self.cmp_scheme(other))

params = {
  "PAM": "NGG",
  "protospacer_len": 20,
  "prime5": True,
  "scoring": "Azimuth",
  "quantity": 100,
  "functional_domains": False,
  "mer_len": 20
}

# azimuth model
print "loading azimuth models", time.time() - start_time
azimuth_saved_model_dir = os.path.join(os.path.dirname(azimuth.__file__), 'saved_models')
model_name = 'V3_model_full.pickle'
azimuth_model_file = os.path.join(azimuth_saved_model_dir, model_name)
with open(azimuth_model_file, 'rb') as f:
  azimuth_model = pickle.load(f)

azimuth_scores_file = 'azimuth_scores.p'
with open(azimuth_scores_file, 'rb') as inp:
  azimuth_scores = pickle.load(inp)

def get_azimuth_score(mer30):
  if mer30 in azimuth_scores:
    return azimuth_scores[mer30]
  else:
    score = azimuth.model_comparison.predict(np.array([mer30]), aa_cut=None, percent_peptide=None, model=azimuth_model, model_file=azimuth_model_file)[0]
    print "generating Azimuth", mer30, score
    azimuth_scores[mer30] = score
    return score

# load in exome
APP_STATIC = "/home/joshm/GUIDES/CRISPR-Library-Designer/static"
exome_seq_path = os.path.join(APP_STATIC, 'data', 'GRCh37_exons')
mer_len = params['mer_len']

# process kmers
# consider all kmers which are followed by NGG
print "preparing hum kmers", time.time() - start_time
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
      if s in exome_mers:
        exome_mers[s] += 1
      else:
        exome_mers[s] = 1
print 'len(exome_mers) = ', len(exome_mers), time.time() - start_time

# takes in guide OBJECT
# returns whether there is a duplicate in exome
def hasExomeRepeat(protospacer):
  guide_seq = protospacer[-mer_len:] # get PAM-proximal mer_len bases
  hits = exome_mers[guide_seq] # how many times does occur in genome followed by NGG?
  return hits >= 2

# loading CFD preprocessed

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

# load preprocessed info
with open("off_target_scores.p", "rb") as inp:
  off_target_scores = pickle.load(inp)

print 'len(off_target_scores) = ', len(off_target_scores), time.time() - start_time

def get_off_target_score(protospacer):
  if hasExomeRepeat(protospacer):
    return 100000
  if not protospacer in off_target_scores:
    score = 0
    off_targets = get_pot_off_targets(protospacer)
    for off_target in off_targets:
        if off_target in exome_mers:
            wt = protospacer + "CGG"
            sg = off_target
            pam = "GG"
            score += exome_mers[off_target] * calc_cfd(wt, sg, pam)
    off_target_scores[protospacer] = score
  return off_target_scores[protospacer]

# Create interval tree for functional domains
print "constructing interval tuples", time.time() - start_time
interval_tuples_dict = {}
ucsc_pfam_f = '../functional_domains/ucsc_pfam.txt'
with open(ucsc_pfam_f, 'r') as pfam_csv:
  csvreader = csv.reader(pfam_csv, delimiter='\t')
  next(csvreader) # skip header

  for row in csvreader:
    chrom = row[1]
    start = row[2]
    end = row[3]
    name = row[4]
    if chrom not in interval_tuples_dict:
      interval_tuples_dict[chrom] = []
    new_tuple = (int(start), int(end), name)
    interval_tuples_dict[chrom].append(new_tuple)

print "constructing interval trees", time.time() - start_time
interval_trees_dict = {}
for k, v in interval_tuples_dict.iteritems():
  interval_trees_dict[k] = IntervalTree.from_tuples(v)

modPAM = params["PAM"].upper()
modPAM = modPAM.replace('N', '[ATCG]')
params["modPAM"] = modPAM
params["PAM_len"] = len(params["PAM"])

revcompl = lambda x: ''.join([{'A':'T','C':'G','G':'C','T':'A','N':'N'}[B] for B in x][::-1])

print "constructing refGene", time.time() - start_time
refGeneFilename = '../gtex/refGene.txt'
refGene = pd.read_csv(refGeneFilename, sep="\t")
refGene.columns=['','name','chrom','strand','txStart','txEnd','cdsStart','cdsEnd','exonCount','exonStarts','exonEnds','id','name2','cdsStartStat','cdsEndStat','exonFrames']
refGene["exonStarts"] = refGene.apply(lambda x: x['exonStarts'].split(',')[:-1], axis=1)
refGene["exonEnds"] = refGene.apply(lambda x: x['exonEnds'].split(',')[:-1], axis=1)
refGene["exonFrames"] = refGene.apply(lambda x: x['exonFrames'].split(',')[:-1], axis=1)

def gene_exon_coords(gene, exon):
  try:
    start = list(refGene.loc[refGene['name'] == gene]['exonStarts'])[0][exon]
    end = list(refGene.loc[refGene['name'] == gene]['exonEnds'])[0][exon]
    chrom = list(refGene.loc[refGene['name'] == gene]['chrom'])[0]
    return {
      'start': int(start),
      'end': int(end),
      'chrom': str(chrom)
    }
  except IndexError:
    return None

def gene_exon_file(gene, exon):
  filename = gene + "_" + str(exon)
  seq_path = os.path.join('../GRCh37_exons/', filename)
  if os.path.isfile(seq_path):
    with open(seq_path) as infile:
      return infile.read()
  else:
    return None

# this is run on multiprocessing workflow
def run(genes_list):
  for gene in genes_list:
    exon = 0
    seq = gene_exon_file(gene["ensembl_id"], exon)
    coords = gene_exon_coords(gene["ensembl_id"], exon)
    while seq:
      # Check if we haven't done this in a previous run of the program
      outfile_name = gene["ensembl_id"] + "_" + str(exon) + ".p"
      folder = '../cfdGRCh37_guides_msgpack_' + params["scoring"] + '/'
      if params['functional_domains']:
        folder = '../cfdGRCh37_guides_msgpack_' + params['scoring'] + '_domains/'
      output_path = os.path.join(folder, outfile_name)

      if os.path.isfile(output_path):
        # prepare next exon
        exon += 1
        seq = gene_exon_file(gene["ensembl_id"], exon)
        coords = gene_exon_coords(gene["ensembl_id"], exon)
        continue

      q = PriorityQueue()
      def process_guide(m, selected, max_queue_size, seq, domain):
        if 'N' in seq:
          return
        PAM_start = m.start()
        score = 0
        if params["scoring"] == "Doench":
          # Doench score requires the 4 before and 6 after 20-mer (gives 30-mer)
          mer30 = seq[PAM_start-params["protospacer_len"]-4:PAM_start+params["PAM_len"]+3]
          if len(mer30) == 30:
            score = doench_score.calc_score(mer30)
        elif params["scoring"] == "Azimuth":
          # Azimuth requires the 4 before and 6 after 20-mer (gives 30-mer)
          mer30 = seq[PAM_start-params["protospacer_len"]-4:PAM_start+params["PAM_len"]+3]
          if len(mer30) == 30:
            score = get_azimuth_score(mer30)
        protospacer = ""
        PAM = ""
        if params["prime5"]:
          protospacer = seq[PAM_start-params["protospacer_len"]:PAM_start]
          PAM = seq[PAM_start:PAM_start+params["PAM_len"]]
        else:
          protospacer = seq[PAM_start+params["PAM_len"]:PAM_start+params["PAM_len"]+params["protospacer_len"]]
          PAM = seq[PAM_start:PAM_start+params["PAM_len"]]
        if protospacer not in exome_mers:
          print protospacer, 'NOT in exome_mers', gene["ensembl_id"], exon
          print 'PAM is', seq[PAM_start:PAM_start+params["PAM_len"]]
        has_exome_repeat = hasExomeRepeat(protospacer)
        off_target_score = get_off_target_score(protospacer)
        potential_gRNA = GuideRNA(selected, PAM_start-params["protospacer_len"], protospacer, PAM, score, exon, gene["ensembl_id"], gene["name"], domain, has_exome_repeat, off_target_score)

        # If there's enough room, add it, no question.
        if q.qsize() < max_queue_size:
          q.put(potential_gRNA)
        # Otherwise, take higher score
        else:
          lowest_gRNA = q.get()
          if cmp(potential_gRNA, lowest_gRNA) == 1: # if potential_gRNA > lowest_gRNA
            q.put(potential_gRNA)
          else:
            q.put(lowest_gRNA)

      for m in re.finditer(params["modPAM"], seq):
        if params["prime5"] and (m.start() < params["protospacer_len"] or m.start() + params["PAM_len"] > len(seq)):
          continue
        elif not params["prime5"] and (m.start() + params["PAM_len"] + params["protospacer_len"] > len(seq)):
          continue

        # Functional domains currently only supported for Cas9.
        # This needs to be modified for other genome editing proteins.
        domain = None
        if params["PAM"] == "NGG": # spCas9
          cut_site = coords['start'] + m.start() - 3
          chrom = 'chr' + coords['chrom']
          if chrom in interval_trees_dict:
            domain_matches = list(interval_trees_dict[chrom][cut_site])
            if len(domain_matches) > 0:
              domain = domain_matches[0].data
        process_guide(m, True, params["quantity"], seq, domain)

      seq_rc = revcompl(seq)

      for m in re.finditer(params["modPAM"], seq_rc):
        if params["prime5"] and (m.start() < params["protospacer_len"] or m.start() + params["PAM_len"] > len(seq)):
          continue
        elif not params["prime5"] and (m.start() + params["PAM_len"] + params["protospacer_len"] > len(seq)):
          continue

        # Functional domains currently only supported for Cas9.
        # This needs to be modified for other genome editing proteins.
        domain = None
        if params["PAM"] == "NGG": #spCas9
          cut_site = coords['end'] - m.start() + 3
          chrom = 'chr' + coords['chrom']
          if chrom in interval_trees_dict:
            domain_matches = list(interval_trees_dict[chrom][cut_site])
            if len(domain_matches) > 0:
              domain = domain_matches[0].data
        process_guide(m, True, params["quantity"], seq_rc, domain)

      # Pop gRNAs into our 'permanent' storage
      gRNAs = []
      while not q.empty():
        gRNA = q.get()
        gRNAs.append(gRNA.serialize_for_display())
      outfile_name = gene["ensembl_id"] + "_" + str(exon) + ".p"
      folder = '../cfdGRCh37_guides_msgpack_' + params['scoring'] + '/'
      if params['functional_domains']:
        folder = '../cfdGRCh37_guides_msgpack_' + params['scoring'] + '_domains/'
      output_path = os.path.join(folder, outfile_name)
      with open(output_path,  'w') as outfile:
        # Reverse gRNAs list.
        # Want highest on-target first.
        msgpack.dump(gRNAs[::-1], outfile)

      # prepare next exon
      exon += 1
      seq = gene_exon_file(gene["ensembl_id"], exon)
      coords = gene_exon_coords(gene["ensembl_id"], exon)

NUM_CORES = 16
print "beginning gene by gene processing", time.time() - start_time
with open('genes_list.json') as genes_list_file:
  full_genes_list = json.load(genes_list_file)
  # gene format: {"ensembl_id": "ENSG00000261122.2", "name": "5S_rRNA", "description": ""}
  processes = []
  unit = len(full_genes_list) / NUM_CORES + 1
  print 'unit is', unit, time.time() - start_time
  for i in range(NUM_CORES):
    start = unit * i
    end = min(unit * (i + 1), len(full_genes_list))
    genes_list = full_genes_list[start:end]
    p = Process(target = run, args=(genes_list,))
    processes.append(p)
  for process in processes:
    process.start()
  for process in processes:
    process.join()

with open('azimuth_scores.p', 'wb') as output:
  pickle.dump(azimuth_scores, output)

end_time = time.time()
hours, rem = divmod(end_time-start_time, 3600)
minutes, seconds = divmod(rem, 60)
print "time elapsed"
print("{:0>2}:{:0>2}:{:05.2f}".format(int(hours),int(minutes),seconds))
