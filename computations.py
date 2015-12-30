import doench_score # calc_score(seq)
from Queue import PriorityQueue
import re
import pandas as pd
import numpy as np
import pickle
from pyensembl import EnsemblRelease
import seq_generator
from settings import APP_STATIC
import os

class GuideRNA():
  """Holder of gRNA information"""
  def __init__(self, selected, start, seq, score, exon_ranking, ensembl_gene, gene_name):
    self.start = start
    self.seq = seq
    self.score = score
    self.exon_ranking = exon_ranking
    self.ensembl_gene = ensembl_gene
    self.gene_name = gene_name
    self.selected = selected

  def serialize_for_display(self):
    """Serialize for the way we are returning json"""
    return {
      "score": self.score,
      "start": self.start,
      "seq": self.seq,
      "selected": self.selected
    }

  # Comparison technique - describes ranking DISCUSS (include exon_ranking?)
  # We assume higher deonch score is better -- otherwise, we should swap self and other scores.
  def __cmp__(self, other):
    return cmp(self.score, other.score)

# Species info is not yet incorporated
class Ranker():
  """Finds and ranks gRNAs from a given gene and species"""
  def __init__(self, genome, species, tissues):
    self.genome = genome
    self.species = species
    self.tissues = tissues
    self.genes = []
    self.gRNAs = []
    self.countSelectedGuides = 0

    # Load pre-processed GTEx data
    self.df_normalized = pickle.load(open(os.path.join(APP_STATIC, 'data/pre_processed', 'pd_by_tissue_normalized.p'), "rb"))

  # ensembl_gene - ensembl encocoding for the gene (grCH37)
  # gene_name - user-friendly name for the gene
  def rank(self, ensembl_gene, gene_name, quantity):
    self.genes.append((ensembl_gene, gene_name))
    df_gene = self.df_normalized[self.df_normalized.Id.str.contains(ensembl_gene)]

    # Sort by exon number, removing first and last exon
    # Recall: Id are entered as ENSG000xxxxx.x_EXONNUM, e.g. ENSG00000000971.11_21
    df_gene['exon_num'] = df_gene['Id'].apply(lambda x: int(x.split('_')[1]))
    if len(df_gene) > 4:
      df_gene = df_gene.sort(['exon_num'], ascending=True).iloc[1:-1]

    # Get median expression for selected tissues
    df_gene['median'] = df_gene[self.tissues].median(axis=1)
    df_results =  df_gene[['Id', 'median', 'exon_num']]
    df_results = df_results.sort(['median'], ascending=False)

    # For this gene, analyze top 4 exons, at most
    q = PriorityQueue()

    exons_to_analyze = min(4, len(df_results))
    total_exons = len(df_results)

    #for i in range(min(4, len(df_results))):
    i = 0
    while (i < total_exons) and (i < exons_to_analyze or q.qsize() < quantity):
      # Generate pseudogenome
      gtex_gene_exon = str(df_results.iloc[i]['Id'])
      gtex_exon_num = int(df_results.iloc[i]['exon_num'])
      seq = str(self.genome.sequence_gtex_gene(gtex_gene_exon))

      # Search for the NGG PAM, beginning after the 21st base pair. DISCUSS (introns needed?)
      # Makes no sense to run if sequence is too short
      if len(seq) < 30:
        i += 1
        continue
      for m in re.finditer(r'(?=GG)', seq[21:]):
        end = m.start() + 21 # add 21 because are looking at seq[21:]
        assert seq[end] == seq[end+1] == 'G'

        # Doench score requires the 4 before and 6 after 20-mer (gives 30-mer)
        # Recall that end points to first G = 21st base pair 
        mer30 = seq[end-22-3:end+8-3]

        # For reference:
        # -- 20mer: seq[end-20 : end]
        # -- PAM: seq[end : end+3]

        # We might too far to the right
        if (len(mer30)) != 30:
          continue

        # Calculate the doench score
        score = doench_score.calc_score(mer30)
        potential_gRNA = GuideRNA(True, end-20, mer30, score, gtex_exon_num, ensembl_gene, gene_name)

        # If there's enough room, add it, no question.
        if q.qsize() < quantity:
          q.put(potential_gRNA)
        # Otherwise, take higher score
        else:
          lowest_gRNA = q.get()
          if potential_gRNA.score > lowest_gRNA.score:
            q.put(potential_gRNA)
          else:
            q.put(lowest_gRNA)

      i += 1

    # Pop gRNAs into our 'permanent' storage 
    while not q.empty():
      gRNA = q.get()
      if gRNA.selected:
        self.countSelectedGuides += 1
      self.gRNAs.append(gRNA)

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
        exon = {
          "start": gene_info['exonStarts'][i] - gene_info['txStart'],
          "end": gene_info['exonEnds'][i] - gene_info['txStart'],
          "gRNAs": [] # Gets filled in below
        }
        gene_to_exon[gene_name]['exons'].append(exon)

    # Iterate through guides and add to appropriate exon
    for guide in self.gRNAs:
      gene_to_exon[guide.gene_name]['exons'][guide.exon_ranking]['gRNAs'].append(guide.serialize_for_display())

    # Convert associative array to ordered array
    guides_by_exon = []
    for _, value in gene_to_exon.iteritems():
      guides_by_exon.append(value)

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
    "human" : seq_generator.Genome()
  }
  ranker = Ranker(genome["human"], species, tissues)

  # Iterate over genes, finding guides for each
  for (ensembl_gene, gene_name) in genes:
    ranker.rank(ensembl_gene, gene_name, quantity_per_gene)

  guides_by_exon = ranker.get_guides_by_exon()
  run()
