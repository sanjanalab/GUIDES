# Python packages
import json
import random
from random_guides import random_guides
import time

# CLD code from parent directories
import os
import inspect
currentdir = os.path.dirname(os.path.abspath(inspect.getfile(inspect.currentframe())))
rootdir = os.path.dirname(os.path.dirname(os.path.dirname(currentdir)))
os.sys.path.insert(0,rootdir)
import computations
import seq_generator

### Constant ranker parameters
genome = {
  "human" : seq_generator.FastGenome()
}
tissues = ['Thyroid', 'Testis', 'Cervix Uteri', 'Adipose Tissue', 'Breast', 'Vagina', 'Nerve', 'Pituitary', 'Stomach', 'Fallopian Tube', 'Bone Marrow', 'Bladder', 'Blood', 'Colon', 'Prostate', 'Pancreas', 'Blood Vessel', 'Liver', 'Spleen', 'Small Intestine', 'Uterus', 'Ovary', 'Muscle', 'Heart', 'Adrenal Gland', 'Brain', 'Salivary Gland', 'Lung', 'Skin', 'Esophagus', 'Kidney']

### Experimental parameters
library_size = 5
num_runs = 1000
num_genes = 500
scoring_alg = "Azimuth"

def compute_average_effeciency_with_on_target(genes, library_size):
  # Call into computations.py, using GTEx_enabled = False, and take results.
  ranker = computations.Ranker(genome["human"], "human", tissues, False, False, scoring_alg = scoring_alg)
  for gene in genes:
    ranker.rank(gene['ensembl_id'], gene['name'], library_size)
  guides_by_exon = ranker.get_guides_by_exon()
  total_score = 0.0
  for gene in guides_by_exon:
    for idx, exon in enumerate(gene['exons']):
      guide_count = 0
      for guide in exon['gRNAs']:
        if guide['selected']:
          total_score += guide['score']
  average_score = total_score / (len(genes) * library_size)
  return average_score

def compute_average_effeciency_without_on_target(genes, library_size):
  # Redo the code in precompute_guides_msgpack.py but don't do any ranking by on-target.
  # Basically, just regex and take from some random order.
  total_score = 0.0
  for gene in genes:
    guides = random_guides(gene, library_size)
    for guide in guides:
      total_score += guide['score']
  average_score = total_score / (len(genes) * library_size)
  return average_score

if __name__ == "__main__":
  print "Selecting guides with and without on-target effeciency. Computing average guide effeciency across these sets. Analyzing {0} genes, with {1} guides/gene. {2} runs will be performed.".format(num_genes, library_size, num_runs)
  t0 = time.time()

  # Stores average results
  results_Y = []
  results_N = []

  # open the list of genes
  with open('../pre_processed/genes_list.json', 'r') as genes_list_file:
    genes_list = json.load(genes_list_file)
    for i in xrange(num_runs):
      genes = random.sample(genes_list, num_genes)
      result_Y = compute_average_effeciency_with_on_target(genes, library_size)
      result_N = compute_average_effeciency_without_on_target(genes, library_size)

      results_Y.append(result_Y)
      results_N.append(result_N)

  filename = 'on_target_contribution.' + str(num_runs) + '.' + scoring_alg + '.results'
  with open(filename, 'w') as results:
    results.write("Average guide effeciency with and without on-target scores considered. Analyzing {0} genes, with {1} guides/gene.\nData from {2} runs is summarized below.\n\n".format(num_genes, library_size, num_runs))
    results.write("Guide selection with on-target considered:\n")
    for result in results_Y:
      results.write(str(result) + "\n")
    results.write("\n\nRandom guide selection:\n")
    for result in results_N:
      results.write(str(result) + "\n")

  print "Completed job in", time.time() - t0, "seconds."
