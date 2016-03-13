import doench_score # calc_score(seq)
from Queue import PriorityQueue
import re
import pandas as pd
import numpy as np
import pickle
import cPickle
import msgpack
from seq_generator import ExonError
import seq_generator
from settings import APP_STATIC
import os
import itertools
import json

df_normalized = pickle.load(open(os.path.join(APP_STATIC, 'data/pre_processed', 'pd_by_tissue_normalized.p'), "rb"))
__module__ = os.path.splitext(os.path.basename(__file__))[0]  ### look here ###

class GuideRNA():
  """Holder of gRNA information"""
  newid = itertools.count().next

  def __init__(self, selected, start, seq, PAM, score, exon_ranking, ensembl_gene, gene_name):
    self.start = start
    self.seq = seq
    self.PAM = PAM
    self.score = score
    self.exon_ranking = exon_ranking
    self.ensembl_gene = ensembl_gene
    self.gene_name = gene_name
    self.selected = selected

    self.uid = "customLibrary_guide" + str(GuideRNA.newid()).zfill(4)

  def serialize_for_display(self):
    """Serialize for the way we are returning json"""
    return {
      "score": self.score,
      "start": self.start,
      "seq": self.seq,
      "PAM": self.PAM,
      "selected": self.selected,
      "uid": self.uid
    }

  # Comparison technique - describes ranking DISCUSS (include exon_ranking?)
  # We assume higher deonch score is better -- otherwise, we should swap self and other scores.
  def __cmp__(self, other):
    return cmp(self.score, other.score)

# Species info is not yet incorporated
class Ranker():
  """Finds and ranks gRNAs from a given gene and species"""
  def __init__(self, genome, species, tissues, gtex_enabled, tissues_enabled, PAM='NGG', prime5=True, protospacer_len=20, scoring_alg={"Doench"}):
    self.genome = genome
    self.species = species
    self.tissues = tissues
    self.gtex_enabled = gtex_enabled
    self.tissues_enabled = tissues_enabled
    self.PAM_len = len(PAM)
    self.prime5 = prime5
    self.protospacer_len = protospacer_len
    self.scoring_alg = scoring_alg

    self.genes = []
    self.gRNAs = []
    self.countSelectedGuides = 0
    self.expression_values = {}

    # Prepare PAM
    PAM = PAM.upper()
    PAM = PAM.replace('N', '[ATCG]')
    self.PAM = "(?=" + PAM + ")"

    # Doench score
    self.use_Doench = "Doench" in scoring_alg

    # Load pre-processed GTEx data
    self.df_normalized = df_normalized

    self.guides_for_exons = {}

  def getGuides(self, gene_exon):
    try:
      filename = gene_exon + ".p"
      path = os.path.join('static/data/GRCh37_guides_msgpack/', filename)
      with open(path) as datafile:
        gRNAs = msgpack.load(datafile)
        return gRNAs
    except IOError:
      gene, exon = gene_exon.split('_')
      raise ExonError(gene, exon)

  # ensembl_gene - ensembl encocoding for the gene (grCH37)
  # gene_name - user-friendly name for the gene
  def rank(self, ensembl_gene, gene_name, quantity):
    # Get the revcompl of a sequence print revcompl("AGTCAGCAT")
    revcompl = lambda x: ''.join([{'A':'T','C':'G','G':'C','T':'A'}[B] for B in x][::-1])

    df_gene_filename = ensembl_gene + ".p"
    df_gene_path = os.path.join('static/data/df_gene_normalized_cpickle', df_gene_filename)
    df_gene_file = open(df_gene_path, 'r')
    df_gene = cPickle.load(df_gene_file)

    self.genes.append((ensembl_gene, gene_name))

    # Sort by exon number, removing first and last exon

    df_gene['median'] = df_gene[self.tissues].median(axis=1)
    expression_values = {}
    for index, row in df_gene.iterrows():
      expression_value = {
        'median': row['median'],
        'overall': row['overall'],
        'brain': row['Brain'],
        'heart': row['Heart'],
        'kidney': row['Kidney'],
        'liver': row['Liver'],
        'skin':  row['Skin']
      }
      expression_values[int(row['exon_num'])] = expression_value
    self.expression_values[gene_name] = expression_values

    total_exons = len(df_gene)

    if len(df_gene) > 4:
      df_gene = df_gene.sort(['exon_num'], ascending=True).iloc[1:-1]

    df_results =  df_gene[['Id', 'median', 'overall','exon_num']]
    df_results = df_results.sort(['median'], ascending=False)

    # For this gene, analyze top 4 exons, at most
    q = PriorityQueue()

    # Ensure we actually have 4 constitutive exons... if we don't, potentially use fewer.
    # NOTE: this is therefore capped at 4.
    constitutive_exon_count = 0
    for i in range(total_exons):
      exp_val = df_results.iloc[i]['median']
      if exp_val > 0.00000000001: # our epsilon value... below this is not real.
        constitutive_exon_count += 1
        if constitutive_exon_count == 4:
          break

    # Focus on top 4 exons if gtex_enabled...otherwise use all.
    exons_to_analyze = len(df_results)
    if self.gtex_enabled == True:
      exons_to_analyze = min(constitutive_exon_count, len(df_results))

    # for i in range(min(4, len(df_results))):
    i = 0
    exon_ends = [0] * total_exons # how many guides in the array do you want to take?

    storedGuides = {}

    while (i < total_exons) and (i < exons_to_analyze or q.qsize() < quantity):
      # Generate pseudogenome
      gtex_gene_exon = str(df_results.iloc[i]['Id'])
      gtex_exon_num = int(df_results.iloc[i]['exon_num'])

      # Get precomputed gRNAs from pickle
      storedGuides[gtex_exon_num] = self.getGuides(gtex_gene_exon)
      gRNAs = storedGuides[gtex_exon_num]

      # If there's enough room, add it, no question.
      add_more = True
      while q.qsize() < quantity:
        if exon_ends[gtex_exon_num] >= len(gRNAs):
          add_more = False
          break
        potential_gRNA = gRNAs[exon_ends[gtex_exon_num]]
        exon_ends[gtex_exon_num] += 1

        q.put((potential_gRNA["score"], gtex_exon_num))

      # Otherwise, keep going until we are too low.
      while add_more:
        if exon_ends[gtex_exon_num] >= len(gRNAs):
          add_more = False
          break
        lowest_gRNA_score, lowest_gRNA_exon = q.get()
        potential_gRNA = gRNAs[exon_ends[gtex_exon_num]]
        if potential_gRNA["score"] > lowest_gRNA_score:
          q.put((potential_gRNA["score"], gtex_exon_num))
          exon_ends[gtex_exon_num] += 1
          exon_ends[lowest_gRNA_exon] -= 1
        else:
          q.put((lowest_gRNA_score, lowest_gRNA_exon))
          add_more = False

      i += 1

    # Mark k = 10 after quantity as unselected
    # add those and beginning ones to guides_for_exons
    for exon_num in storedGuides.keys():
      gRNAs = storedGuides[exon_num]
      for guide in gRNAs[:exon_ends[exon_num]]:
        guide["selected"] = True
      for guide in gRNAs[exon_ends[exon_num]:exon_ends[exon_num]+10]:
        guide["selected"] = False
      human_gene_exon = gene_name + "+" + str(exon_num)
      self.guides_for_exons[human_gene_exon] = gRNAs[:exon_ends[exon_num]+10]
      self.countSelectedGuides += exon_ends[exon_num]

    df_gene_file.close()

  def get_guides_by_exon(self):
    # First, setup data as an associative array, to make guide insertion easier
    # Later, we'll transform to ordered array.
    gene_to_exon = {}
    for (ensembl_gene, gene_name) in self.genes:
      gene_info = self.genome.gene_info(ensembl_gene)

      gene_to_exon[gene_name] = {
        "name": gene_name,
        "ensembl_gene": ensembl_gene,
        "length": gene_info['txEnd'] - gene_info['txStart'],
        "start": gene_info['txStart'],
        "end": gene_info['txEnd'],
        "exons": [] # Gets filled in below
      }

      # Prepare each exon and add to gene_to_exon[gene_name].exons
      for i in range(gene_info['exonCount']):
        expression_value = self.expression_values[gene_name][i]
        expression_value_returned = dict(expression_value)
        if not self.tissues_enabled:
          del expression_value_returned["median"]
        exon = {
          "start": gene_info['exonStarts'][i] - gene_info['txStart'],
          "end": gene_info['exonEnds'][i] - gene_info['txStart'],
          "gRNAs": [], # Gets filled in below
          "expression": expression_value_returned
        }
        gene_to_exon[gene_name]['exons'].append(exon)

    # Iterate through guides and add to appropriate exon
    for k, v in self.guides_for_exons.iteritems():
      gene_name, exon = k.split('+')
      gene_to_exon[gene_name]['exons'][int(exon)]['gRNAs'] = v

    # for guide in self.gRNAs:
    #   gene_to_exon[guide["gene_name"]]['exons'][guide["exon_ranking"]]['gRNAs'].append(guide)

    # Convert associative array to ordered array
    guides_by_exon = gene_to_exon.values()
    return guides_by_exon

  def get_count_selected_guides(self):
    return self.countSelectedGuides

def run():
  suma = 0
  for i in range(5):
    suma += i
  return suma

# for testing/profiling issues
if __name__ == '__main__':
    # Request arguments
  genes = None
  species  = None
  quantity = None
  tissues  = None

  # Validations
  if genes == None:
    genes = [('ENSG00000115977.14', "AAK1"), ('ENSG00000142168.10', "SOD1")]
  if species == None:
    species = "human"
  if quantity == None:
    quantity = 60
  if tissues == None:
    tissues = ["Muscle", "Heart", "Brain"]

  quantity_per_gene = quantity/len(genes)

  # Setup ranker
  genome = {
    "human" : seq_generator.FastGenome()
  }
  ranker = Ranker(genome["human"], species, tissues)

  # Iterate over genes, finding guides for each
  for (ensembl_gene, gene_name) in genes:
    ranker.rank(ensembl_gene, gene_name, quantity_per_gene)

  guides_by_exon = ranker.get_guides_by_exon()
  run()
