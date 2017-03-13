import doench_score # calc_score(seq)
from Queue import PriorityQueue
import re
import pandas as pd
import numpy as np
import pickle
import cPickle
import msgpack
import seq_generator_mus
from settings import APP_STATIC
import os
import itertools
import json
import operator

with open(os.path.join(APP_STATIC, 'data/pre_processed', 'strand_info_GRCm38.p'), "rb") as infile:
  gene_strand_mapping = cPickle.load(infile)

# Load refGene
refGeneFilename = os.path.join(APP_STATIC, 'data/gtex/gtex_mouse', 'refGene_mouse.txt')
refGene = pd.read_csv(refGeneFilename, sep="\t")
refGene.columns=['','name','chrom','strand','txStart','txEnd','cdsStart','cdsEnd','exonCount','exonStarts','exonEnds','id','name2','cdsStartStat','cdsEndStat','exonFrames']
refGene["exonStarts"] = refGene.apply(lambda x: x['exonStarts'].split(',')[:-1], axis=1)
refGene["exonEnds"] = refGene.apply(lambda x: x['exonEnds'].split(',')[:-1], axis=1)
refGene["exonFrames"] = refGene.apply(lambda x: x['exonFrames'].split(',')[:-1], axis=1)

__module__ = os.path.splitext(os.path.basename(__file__))[0]  ### look here ###

# Species info is not yet incorporated
class RankerMouse():
  """Finds and ranks gRNAs from a given gene and species"""
  def __init__(self, genome, species, domains_enabled, PAM='NGG', prime5=True, protospacer_len=20, scoring_alg="Azimuth"):
    self.genome = genome
    self.species = species
    self.domains_enabled = domains_enabled
    self.PAM_len = len(PAM)
    self.prime5 = prime5
    self.protospacer_len = protospacer_len
    self.scoring_alg = scoring_alg

    self.genes = []
    self.gRNAs = []
    self.countSelectedGuides = 0

    # Prepare PAM
    PAM = PAM.upper()
    PAM = PAM.replace('N', '[ATCG]')
    self.PAM = "(?=" + PAM + ")"

    self.guides_for_exons = {}

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
      gene, exon = gene_exon.split('_')
      print 'IOError in getGuides', gene, exon
      raise seq_generator_mus.ExonError(gene, exon)

  # ensembl_gene - ensembl encocoding for the gene (grCH37)
  # gene_name - user-friendly name for the gene
  def rank(self, ensembl_gene, gene_name, quantity):
    # Get the revcompl of a sequence print revcompl("AGTCAGCAT")
    revcompl = lambda x: ''.join([{'A':'T','C':'G','G':'C','T':'A'}[B] for B in x][::-1])

    self.genes.append((ensembl_gene, gene_name))
    location = refGene.loc[refGene['name'] == ensembl_gene]
    num_exons = int(list(location['exonCount'])[-1])

    # Sort by exon number, removing first and last exon
    if gene_strand_mapping[ensembl_gene] == '+':
      first = 0
      last = num_exons - 1
    else:
      first = num_exons - 1
      last = 0
    ordering = range(1, num_exons - 1) + [first, last]

    # preload all the guides
    storedGuides = {}
    for i in range(num_exons):
      # Get precomputed gRNAs from pickle
      gtex_gene_exon = ensembl_gene + "_" + str(i)
      storedGuides[i] = self.getGuides(gtex_gene_exon)

    # choose which guides we want.
    # this is represented by the exon_ends list.
    exon_ends = [0] * num_exons # how many guides in the array do you want to take?

    # maintain 3 lists of exons
    def range_with_guide_containing_exons(r):
      result = []
      for i in r:
        if len(storedGuides[i]) > 0:
          result.append(i)
      return result

    # primary --> where we would like to take exons from
    # secondary --> when primary is empty, choose from here.
    if num_exons >= 6:
      primary = range_with_guide_containing_exons(range(0,4)) # 0,1,2,3
      secondary = range_with_guide_containing_exons(range(4, num_exons)) # all but first 4
    else:
      primary = range_with_guide_containing_exons(range(0, num_exons - 2))
      secondary = range_with_guide_containing_exons(range(max(num_exons - 2, 0), num_exons)) # last two exons

    # now execute the data structure
    total_guides_selected = 0
    while total_guides_selected < quantity and (len(primary) > 0 or len(secondary) > 0):
      # get the list of guides
      potential_guides = []
      if len(primary) == 0:
        primary = [secondary[0]]
        secondary = secondary[1:]
      for i in primary:
        guides_for_exon = storedGuides[i]
        potential_guide = guides_for_exon[exon_ends[i]]
        potential_guides.append(potential_guide)

      # potential_guides is an array of the best guide from each exon in primary
      def cmp_key(x):
        _, guide = x
        if not 'has_functional_domain' in guide:
          guide['has_functional_domain'] = False
        return (-guide['off_target_score'], guide['has_functional_domain'], guide['score'])

      # find which guide is best
      max_index, max_value = max(enumerate(potential_guides), key=cmp_key)
      exon_chosen = primary[max_index] # which row in the dataframe

      # if our best possible guide has an off-target score of 100000,
      # then we should search secondary for a guide with lower off-target score.
      # If we find one, use that guide instead!
      # Then, go back to the beginning of the while loop (continue)
      # Otherwise, work with the info we've gotten above

      ###################### - should not affect any variables we need below
      got_new_guide = False
      if max_value['off_target_score'] == 100000:
        for idx,i in enumerate(secondary):
          gtex_exon_num = i
          guides_for_exon = storedGuides[gtex_exon_num]
          if exon_ends[gtex_exon_num] >= len(guides_for_exon): # we've used up all of these.
            del secondary[idx]
            continue
          potential_guide = guides_for_exon[exon_ends[gtex_exon_num]]
          if potential_guide['off_target_score'] < 100000:
            # use this guide!
            got_new_guide = True
            exon_ends[gtex_exon_num] += 1
            total_guides_selected += 1
            break
        # if we pulled a guide above, continue in the while loop (don't execute below)
        if got_new_guide:
          continue
      ######################

      # convert to actual exon num
      gtex_exon_num = exon_chosen
      exon_ends[gtex_exon_num] += 1
      total_guides_selected += 1

      # If we're out of guides, remove this
      guides_for_exon = storedGuides[gtex_exon_num]
      if exon_ends[gtex_exon_num] == len(guides_for_exon): # exon_ends is exclusive --> i.e., don't include the value listed in exon_ends
        del primary[max_index]

      # enforce sparsity: don't let a single exon account for more than 50% of guides
      # move the exon from primary to back of secondary
      elif exon_ends[gtex_exon_num] >= quantity / 2:
        secondary.append(primary[max_index]) # exon_chosen
        del primary[max_index]

    # Mark k = 10 after quantity as unselected
    # add those and beginning ones to guides_for_exons
    for exon_num in xrange(num_exons):
      gRNAs = storedGuides[exon_num]
      for guide in gRNAs[:exon_ends[exon_num]]:
        guide["selected"] = True
      for guide in gRNAs[exon_ends[exon_num]:exon_ends[exon_num]+10]:
        guide["selected"] = False
      human_gene_exon = ensembl_gene + "+" + str(exon_num)
      self.guides_for_exons[human_gene_exon] = gRNAs[:exon_ends[exon_num]+10]
      self.countSelectedGuides += exon_ends[exon_num]

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
        "end": gene_info['txEnd'],
        "exons": [] # Gets filled in below
      }

      # Prepare each exon and add to gene_to_exon[ensembl_gene].exons
      if gene_strand_mapping[ensembl_gene] == '+':
        r = xrange(exon_count)
      else:
        r = xrange(exon_count-1, -1, -1)

      for i in r:
        exon = {
          "gRNAs": [] # Gets filled in below
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
