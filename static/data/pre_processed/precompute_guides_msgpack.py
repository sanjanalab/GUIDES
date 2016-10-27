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

class GuideRNA():
  """Holder of gRNA information"""
  def __init__(self, selected, start, seq, PAM, score, exon_ranking, ensembl_gene, gene_name, functional_domain=None):
    self.start = start
    self.seq = seq
    self.PAM = PAM
    self.score = score
    self.exon_ranking = exon_ranking
    self.ensembl_gene = ensembl_gene
    self.gene_name = gene_name
    self.selected = selected
    self.functional_domain = functional_domain

  def serialize_for_display(self):
    """Serialize for the way we are returning json"""
    serialization = {
      "score": self.score,
      "start": self.start,
      "seq": self.seq,
      "PAM": self.PAM,
      "selected": self.selected,
    }
    if self.functional_domain != None:
      serialization["functional_domain"] = self.functional_domain
    return serialization

  def __cmp__(self, other):
    return cmp(self.score, other.score)

params = {
  "PAM": "NGG",
  "protospacer_len": 20,
  "prime5": True,
  "scoring": "Azimuth",
  "quantity": 100,
  "functional_domains": True
}

# azimuth mdoel
azimuth_saved_model_dir = os.path.join(os.path.dirname(azimuth.__file__), 'saved_models')
model_name = 'V3_model_full.pickle'
azimuth_model_file = os.path.join(azimuth_saved_model_dir, model_name)
with open(azimuth_model_file, 'rb') as f:
  azimuth_model = pickle.load(f)

# Create interval tree for functional domains
print "constructing interval tuples"
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

print "constructing interval trees"
interval_trees_dict = {}
for k, v in interval_tuples_dict.iteritems():
  interval_trees_dict[k] = IntervalTree.from_tuples(v)

modPAM = params["PAM"].upper()
modPAM = modPAM.replace('N', '[ATCG]')
params["modPAM"] = modPAM
params["PAM_len"] = len(params["PAM"])

revcompl = lambda x: ''.join([{'A':'T','C':'G','G':'C','T':'A','N':'N'}[B] for B in x][::-1])

print "constructing refGene"
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

with open('genes_list.json') as genes_list_file:
  genes_list = json.load(genes_list_file)
  # gene format: {"ensembl_id": "ENSG00000261122.2", "name": "5S_rRNA", "description": ""}
  for gene in genes_list:
    exon = 0
    seq = gene_exon_file(gene["ensembl_id"], exon)
    coords = gene_exon_coords(gene["ensembl_id"], exon)
    while seq:
      # Check if we haven't done this in a preivous run of the program
      outfile_name = gene["ensembl_id"] + "_" + str(exon) + ".p"
      folder = '../GRCh37_guides_msgpack_' + params["scoring"] + '/'
      if params['functional_domains']:
        folder = '../GRCh37_guides_msgpack_' + params['scoring'] + '_domains/'
      output_path = os.path.join(folder, outfile_name)

      if os.path.isfile(output_path):
        # prepare next exon
        exon += 1
        seq = gene_exon_file(gene["ensembl_id"], exon)
        coords = gene_exon_coords(gene["ensembl_id"], exon)
        continue

      q = PriorityQueue()
      domain_q = PriorityQueue()
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
            score = azimuth.model_comparison.predict(np.array([mer30]), aa_cut=None, percent_peptide=None, model=azimuth_model, model_file=azimuth_model_file)[0]
        protospacer = ""
        PAM = ""
        if params["prime5"]:
          protospacer = seq[PAM_start-params["protospacer_len"]:PAM_start]
          PAM = seq[PAM_start:PAM_start+params["PAM_len"]]
        else:
          protospacer = seq[PAM_start+params["PAM_len"]:PAM_start+params["PAM_len"]+params["protospacer_len"]]
          PAM = seq[PAM_start:PAM_start+params["PAM_len"]]
        potential_gRNA = GuideRNA(selected, PAM_start-params["protospacer_len"], protospacer, PAM, score, exon, gene["ensembl_id"], gene["name"], domain)

        if domain:
          domain_q.put(potential_gRNA)
        # If there's enough room, add it, no question.
        elif q.qsize() < max_queue_size:
          q.put(potential_gRNA)
        # Otherwise, take higher score
        else:
          lowest_gRNA = q.get()
          if potential_gRNA.score > lowest_gRNA.score:
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
      count = 0
      gRNAs = []
      while not q.empty() and count < params["quantity"]:
        gRNA = q.get()
        gRNAs.append(gRNA.serialize_for_display())
        count = count + 1
      while not domain_q.empty() and count < params["quantity"]:
        gRNA = domain_q.get()
        gRNAs.append(gRNA.serialize_for_display())
        count = count + 1
      domain_count = count
      outfile_name = gene["ensembl_id"] + "_" + str(exon) + ".p"
      if domain_count > 0:
        print "for {0} we had {1} domain and {2} ordinary guides.".format(outfile_name, domain_count, count - domain_count)
      folder = '../GRCh37_guides_msgpack_' + params['scoring'] + '/'
      if params['functional_domains']:
        folder = '../GRCh37_guides_msgpack_' + params['scoring'] + '_domains/'
      output_path = os.path.join(folder, outfile_name)
      with open(output_path,  'w') as outfile:
        # Reverse gRNAs list.
        # Want highest on-target first.
        msgpack.dump(gRNAs[::-1], outfile)

      # prepare next exon
      exon += 1
      seq = gene_exon_file(gene["ensembl_id"], exon)
      coords = gene_exon_coords(gene["ensembl_id"], exon)
