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

with open(os.path.join(APP_STATIC, 'data/pre_processed', 'pd_by_tissue_normalized.p'), "rb") as infile:
  df_normalized = pickle.load(infile)

with open(os.path.join(APP_STATIC, 'data/pre_processed', 'strand_info.p'), "rb") as infile:
  gene_strand_mapping = cPickle.load(infile)

__module__ = os.path.splitext(os.path.basename(__file__))[0]  ### look here ###

# Species info is not yet incorporated
class Ranker():
  """Finds and ranks gRNAs from a given gene and species"""
  def __init__(self, genome, species, tissues, gtex_enabled, tissues_enabled, domains_enabled, PAM='NGG', prime5=True, protospacer_len=20, scoring_alg="Azimuth"):
    self.genome = genome
    self.species = species
    self.tissues = tissues
    self.gtex_enabled = gtex_enabled
    self.tissues_enabled = tissues_enabled
    self.domains_enabled = domains_enabled
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

    # Load pre-processed GTEx data
    self.df_normalized = df_normalized

    self.guides_for_exons = {}

  def getGuides(self, gene_exon):
    try:
      filename = gene_exon + ".p"
      if self.scoring_alg == "Doench":
        path = os.path.join(os.path.dirname(__file__), 'static/data/GRCh37_guides_msgpack_Doench/', filename)
      elif self.domains_enabled:
        path = os.path.join(os.path.dirname(__file__), 'static/data/GRCh37_guides_msgpack_Azimuth_domains/', filename)
      else:
        path = os.path.join(os.path.dirname(__file__), 'static/data/GRCh37_guides_msgpack_Azimuth/', filename)
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
    df_gene_path = os.path.join(os.path.dirname(__file__), 'static/data/df_gene_normalized_cpickle', df_gene_filename)
    df_gene_file = open(df_gene_path, 'r')
    df_gene = cPickle.load(df_gene_file)

    self.genes.append((ensembl_gene, gene_name))

    # Sort by exon number, removing first and last exon

    if self.tissues_enabled:
      df_gene['median'] = df_gene[self.tissues].median(axis=1)
    expression_values = {}
    constitutive_exon_count = 0
    for index, row in df_gene.iterrows():
      # Remember expression values
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

      # Count number of constitutive exons
      if row['median'] > 0.0000000001: # our epsilon value... below this is not real.
        constitutive_exon_count += 1

    self.expression_values[ensembl_gene] = expression_values

    total_exons = len(df_gene)

    # Only consider first and last exon if we don't have enough constitutive exons.
    if not (self.gtex_enabled and constitutive_exon_count < 4) and len(df_gene) > 2:
      df_gene = df_gene[(df_gene.exon_num != 0) & (df_gene.exon_num != len(df_gene) - 1)]
    df_results = df_gene[['Id', 'median', 'overall', 'exon_num']]

    if self.tissues_enabled:
      df_results = df_results.sort(['median'], ascending=False)

    # For this gene, analyze top 4 exons, at most
    q = PriorityQueue()

    # Focus on top 4 exons if gtex_enabled...otherwise use all.
    exon_entries = len(df_results) # we might have removed first and last
    exons_to_analyze = exon_entries
    if self.gtex_enabled == True:
      exons_to_analyze = min(4, constitutive_exon_count, exon_entries)

    # for i in range(min(4, len(df_results))):
    i = 0
    exon_ends = [0] * total_exons # how many guides in the array do you want to take?

    storedGuides = {}

    for i in range(total_exons):
      # Get precomputed gRNAs from pickle
      gtex_gene_exon = ensembl_gene + "_" + str(i)
      storedGuides[i] = self.getGuides(gtex_gene_exon)

    i = 0
    while (i < exon_entries) and (i < total_exons) and (i < exons_to_analyze or q.qsize() < quantity):
      gtex_gene_exon = str(df_results.iloc[i]['Id'])
      gtex_exon_num = int(df_results.iloc[i]['exon_num'])
      gRNAs = storedGuides[gtex_exon_num]

      # If there's enough room, add it, no question.
      add_more = True
      while q.qsize() < quantity:
        if exon_ends[gtex_exon_num] >= len(gRNAs):
          add_more = False
          break
        potential_gRNA = gRNAs[exon_ends[gtex_exon_num]]
        exon_ends[gtex_exon_num] += 1
        functional_presence = 0
        print potential_gRNA
        if "functional_domain" in potential_gRNA:
          print "yolo"
          functional_presence = 1
        q.put((functional_presence, potential_gRNA["score"], gtex_exon_num))

      # Otherwise, keep going until we are too low.
      while add_more:
        if exon_ends[gtex_exon_num] >= len(gRNAs):
          add_more = False
          break
        lowest_functional_presence, lowest_gRNA_score, lowest_gRNA_exon = q.get()
        potential_gRNA = gRNAs[exon_ends[gtex_exon_num]]
        functional_presence = 0
        if "functional_domain" in potential_gRNA:
          print "yolo"
          functional_presence = 1
        if (functional_presence and not lowest_functional_presence) or (functional_presence == lowest_functional_presence and potential_gRNA["score"] > lowest_gRNA_score):
          q.put((functional_presence, potential_gRNA["score"], gtex_exon_num))
          exon_ends[gtex_exon_num] += 1
          exon_ends[lowest_gRNA_exon] -= 1
        else:
          q.put((lowest_functional_presence, lowest_gRNA_score, lowest_gRNA_exon))
          add_more = False

      i += 1

    # Mark k = 10 after quantity as unselected
    # add those and beginning ones to guides_for_exons
    for exon_num in xrange(total_exons):
      gRNAs = storedGuides[exon_num]
      for guide in gRNAs[:exon_ends[exon_num]]:
        guide["selected"] = True
      for guide in gRNAs[exon_ends[exon_num]:exon_ends[exon_num]+10]:
        guide["selected"] = False
      human_gene_exon = ensembl_gene + "+" + str(exon_num)
      self.guides_for_exons[human_gene_exon] = gRNAs[:exon_ends[exon_num]+10]
      self.countSelectedGuides += exon_ends[exon_num]

    df_gene_file.close()

  def get_guides_by_exon(self):
    # First, setup data as an associative array, to make guide insertion easier
    # Later, we'll transform to ordered array.
    gene_to_exon = {}

    # hold onto exon_counts as they are computed,
    # so we don't have to look them up again later.
    exon_counts = {}

    for (ensembl_gene, gene_name) in self.genes:
      gene_info = self.genome.gene_info(ensembl_gene)
      exon_count = gene_info['exonCount']
      exon_counts[ensembl_gene] = exon_count

      gene_to_exon[ensembl_gene] = {
        "name": gene_name,
        "ensembl_gene": ensembl_gene,
        "length": gene_info['txEnd'] - gene_info['txStart'],
        "start": gene_info['txStart'],
        "end": gene_info['txEnd'],
        "exons": [] # Gets filled in below
      }

      # Prepare each exon and add to gene_to_exon[ensembl_gene].exons
      if gene_strand_mapping[ensembl_gene] == '+':
        gene_to_exon[ensembl_gene]['start'] = gene_info['txStart']
        gene_to_exon[ensembl_gene]['end'] = gene_info['txEnd']
        r = xrange(exon_count)
      else:
        gene_to_exon[ensembl_gene]['start'] = gene_info['txEnd']
        gene_to_exon[ensembl_gene]['end'] = gene_info['txStart']
        r = xrange(exon_count-1, -1, -1)

      for i in r:
        expression_value = self.expression_values[ensembl_gene][i]
        expression_value_returned = dict(expression_value)
        if (not self.tissues_enabled) and ("median" in expression_value_returned):
          del expression_value_returned["median"]
        exon = {
          "gRNAs": [], # Gets filled in below
          "expression": expression_value_returned
        }

        # for + strand, start < end
        # for - strand, end < start
        if gene_strand_mapping[ensembl_gene] == '+':
          exon['start'] = gene_info['exonStarts'][i] - gene_info['txStart']
          exon['end'] = gene_info['exonEnds'][i] - gene_info['txStart']
        else:
          exon['start'] = gene_info['exonEnds'][i] - gene_info['txStart']
          exon['end'] = gene_info['exonStarts'][i] - gene_info['txStart']

        gene_to_exon[ensembl_gene]['exons'].append(exon)

    # Iterate through guides and add to appropriate exon
    for k, v in self.guides_for_exons.iteritems():
      ensembl_gene, exon = k.split('+')
      if gene_strand_mapping[ensembl_gene] == '+':
        exon_num = int(exon)
      else:
        exon_num = exon_counts[ensembl_gene] - 1 - int(exon)
      gene_to_exon[ensembl_gene]['exons'][exon_num]['gRNAs'] = v

    # for guide in self.gRNAs:
    #   gene_to_exon[guide["gene_name"]]['exons'][guide["exon_ranking"]]['gRNAs'].append(guide)

    # Convert associative array to ordered array
    guides_by_exon = gene_to_exon.values()
    return guides_by_exon

  def get_count_selected_guides(self):
    return self.countSelectedGuides
