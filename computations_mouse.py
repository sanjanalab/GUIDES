#
# Rewriting Ranker for mouse genome
# We will call
# rank(self, gene_name, quantity)
# as usual.
#
# Thereafter, this needs ot still implement
#       get_guides_by_exon()
#

import doench_score # calc_score(seq)
from Queue import PriorityQueue
import re
import pandas as pd
import numpy as np
import pickle
import cPickle
import msgpack
from settings import APP_STATIC
import os
import itertools
import json

# Load gene_strand_mapping
with open(os.path.join(APP_STATIC, 'data/pre_processed', 'strand_info_GRCm38.p'), "rb") as infile:
  gene_strand_mapping = cPickle.load(infile)

# Load refGene
refGeneFilename = os.path.join(APP_STATIC, 'data/gtex', 'refGene_GRCm38.txt')
refGene = pd.read_csv(refGeneFilename, sep="\t")
refGene.columns=['','name','chrom','strand','txStart','txEnd','cdsStart','cdsEnd','exonCount','exonStarts','exonEnds','id','name2','cdsStartStat','cdsEndStat','exonFrames']
refGene["exonStarts"] = refGene.apply(lambda x: x['exonStarts'].split(',')[:-1], axis=1)
refGene["exonEnds"] = refGene.apply(lambda x: x['exonEnds'].split(',')[:-1], axis=1)
refGene["exonFrames"] = refGene.apply(lambda x: x['exonFrames'].split(',')[:-1], axis=1)

# Where to look for files
__module__ = os.path.splitext(os.path.basename(__file__))[0]  ### look here ###

class RankerMouse():
  """Finds and ranks gRNAs from a given gene and species"""
  def __init__(self, domains_enabled, PAM='NGG', prime5=True, protospacer_len=20, scoring_alg="Azimuth"):
    self.domains_enabled = domains_enabled
    self.PAM_len = len(PAM)
    self.prime5 = prime5
    self.protospacer_len = protospacer_len
    self.scoring_alg = scoring_alg

    # Local data storage for Ranker
    self.genes = []
    self.gRNAs = []
    self.countSelectedGuides = 0

    # Prepare PAM
    PAM = PAM.upper()
    PAM = PAM.replace('N', '[ATCG]')
    self.PAM = "(?=" + PAM + ")"

    # Ultimate return values
    self.guides_for_exons = {}

  # Fetch precomputed guides
  def getGuides(self, gene_exon):
    try:
      filename = gene_exon + ".p"
      if self.scoring_alg == "Doench":
        path = os.path.join(os.path.dirname(__file__), 'static/data/GRCm38_guides_msgpack_Doench/', filename)
      elif self.domains_enabled:
        path = os.path.join(os.path.dirname(__file__), 'static/data/GRCm38_guides_msgpack_Azimuth_domains/', filename)
      else:
        path = os.path.join(os.path.dirname(__file__), 'static/data/GRCm38_guides_msgpack_Azimuth/', filename)
      with open(path) as datafile:
        gRNAs = msgpack.load(datafile)
        return gRNAs
    except IOError:
      raise Exception("Could not find gene,exon: {0} with that combination".format(gene_exon))

  def rank(self, gene_name, quantity):
    # Get the revcompl of a sequence print revcompl("AGTCAGCAT")
    revcompl = lambda x: ''.join([{'A':'T','C':'G','G':'C','T':'A'}[B] for B in x][::-1])

    self.genes.append(gene_name)

    location = refGene.loc[refGene['name2'] == gene_name]
    exonCount = list(location['exonCount'])[-1]

    # Discard first and last exons, unless we won't have 4 otherwise
    rangeExons = None # declare
    if exonCount > 7:
      rangeExons = range(1, exonCount - 1) # skip first and last
    else:
      rangeExons = range(exonCount) # use all exons

    # store the results
    q = PriorityQueue()

    # Load guides into memory
    storedGuides = {}
    exon_ends = {}
    for exonNum in rangeExons:
      # Get precomputed gRNAs from msgpack
      gtex_gene_exon = gene_name + "_" + str(exonNum)
      storedGuides[exonNum] = self.getGuides(gtex_gene_exon)
      exon_ends[exonNum] = 0

      # If there's enough room, add all, no question
      add_more = True
      while q.qsize() < quantity:
        if exon_ends[exonNum] >= len(storedGuides[exonNum]):
          add_more = False
          break
        potential_gRNA = storedGuides[exonNum][exon_ends[exonNum]]
        exon_ends[exonNum] += 1

        q.put((potential_gRNA["score"], exonNum))

      # Otherwise, keep going until we are too low
      while add_more:
        if exon_ends[exonNum] >= len(storedGuides[exonNum]):
          add_more = False
          break
        lowest_gRNA_score, lowest_gRNA_exon = q.get()
        potential_gRNA = storedGuides[exonNum][exon_ends[exonNum]]
        if potential_gRNA["score"] > lowest_gRNA_score:
          q.put((potential_gRNA["score"], exonNum))
          exon_ends[exonNum] += 1
          exon_ends[lowest_gRNA_exon] -= 1
        else:
          q.put((lowest_gRNA_score, lowest_gRNA_exon))
          add_more = False

    # Mark k = 10 after quantity as unselected
    # add those and beginning ones to guides_for_exons
    for exonNum in rangeExons:
      gRNAs = storedGuides[exonNum]
      for guide in gRNAs[:exon_ends[exonNum]]:
        guide["selected"] = True
      for guide in gRNAs[exon_ends[exonNum]:exon_ends[exonNum]+10]:
        guide["selected"] = False
      human_gene_exon = gene_name + "+" + str(exonNum)
      self.guides_for_exons[human_gene_exon] = gRNAs[:exon_ends[exonNum]+10]
      self.countSelectedGuides += exon_ends[exonNum]

  def get_guides_by_exon(self):
    # First, setup data as an associative array, to make guide insertion easier
    # Later, we'll transform to ordered array.
    gene_to_exon = {}

    # hold onto exon_counts as they are computed,
    # so we don't have to look them up again later.
    exon_counts = {}

    for gene_name in self.genes:
      gene_info = refGene.loc[refGene['name2'] == gene_name]
      exon_count = list(gene_info['exonCount'])[-1]
      txStart = list(gene_info['txStart'])[-1]
      txEnd = list(gene_info['txEnd'])[-1]
      exonStarts = list(gene_info['exonStarts'])[-1]
      exonEnds = list(gene_info['exonEnds'])[-1]

      exon_counts[gene_name] = exon_count

      gene_to_exon[gene_name] = {
        "name": gene_name,
        "length": txEnd - txStart,
        "start": txStart,
        "end": txEnd,
        "exons": [] # Gets filled in below
      }

      # Prepare each exon and add to gene_to_exon[ensembl_gene].exons
      if gene_strand_mapping[gene_name] == '+':
        gene_to_exon[gene_name]['start'] = txStart
        gene_to_exon[gene_name]['end'] = txEnd
        r = xrange(exon_count)
      else:
        gene_to_exon[gene_name]['start'] = txEnd
        gene_to_exon[gene_name]['end'] = txStart
        r = xrange(exon_count-1, -1, -1)

      for i in r:
        exon = {
          "gRNAs": [] # Gets filled in below
        }

        # for + strand, start < end
        # for - strand, end < start
        if gene_strand_mapping[gene_name] == '+':
          exon['start'] = int(exonStarts[i]) - txStart
          exon['end'] = int(exonEnds[i]) - txStart
        else:
          exon['start'] = int(exonEnds[i]) - txStart
          exon['end'] = int(exonStarts[i]) - txStart

        gene_to_exon[gene_name]['exons'].append(exon)

    # Iterate through guides and add to appropriate exon
    for k, v in self.guides_for_exons.iteritems():
      gene_name, exon = k.split('+')
      if gene_strand_mapping[gene_name] == '+':
        exon_num = int(exon)
      else:
        exon_num = exon_counts[gene_name] - 1 - int(exon)
      gene_to_exon[gene_name]['exons'][exon_num]['gRNAs'] = v

    # Convert associative array to ordered array
    guides_by_exon = gene_to_exon.values()
    return guides_by_exon

  def get_count_selected_guides(self):
    return self.countSelectedGuides
