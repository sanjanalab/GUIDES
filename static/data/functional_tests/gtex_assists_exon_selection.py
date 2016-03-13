# Python packages
import cPickle
import json
import random
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

### User selected parameters
library_size = 5
num_runs = 1000
num_genes = 500

def overall_gene_expression(ensembl_gene):
  df_gene_filename = ensembl_gene + ".p"
  df_gene_path = os.path.join('../df_gene_normalized_cpickle', df_gene_filename)
  df_gene_file = open(df_gene_path, 'r')
  df_gene = cPickle.load(df_gene_file)

def gtex_effect(genes):
  sum_normalized_expression_Y = 0
  sum_normalized_expression_N = 0

  ranker_GTEX_Y = computations.Ranker(genome["human"], "human", tissues, True, False)
  ranker_GTEX_N = computations.Ranker(genome["human"], "human", tissues, False, False)

  for gene in genes:
    ranker_GTEX_Y.rank(gene['ensembl_id'], gene['name'], library_size)
    ranker_GTEX_N.rank(gene['ensembl_id'], gene['name'], library_size)

  guides_by_exon_GTEX_Y = ranker_GTEX_Y.get_guides_by_exon()
  guides_by_exon_GTEX_N = ranker_GTEX_N.get_guides_by_exon()

  for gene in guides_by_exon_GTEX_Y:
    for idx, exon in enumerate(gene['exons']):
      guide_count = 0
      for guide in exon['gRNAs']:
        if guide['selected']:
          guide_count += 1
      sum_normalized_expression_Y +=  guide_count * exon['expression']['overall']

  for gene in guides_by_exon_GTEX_N:
    for idx, exon in enumerate(gene['exons']):
      guide_count = 0
      for guide in exon['gRNAs']:
        if guide['selected']:
          guide_count += 1
      sum_normalized_expression_N +=  guide_count * exon['expression']['overall']

  sum_normalized_expression_Y /= len(genes)
  sum_normalized_expression_N /= len(genes)

  return sum_normalized_expression_Y, sum_normalized_expression_N

if __name__ == "__main__":
  print "Computing average of normalized exon expression across {0} genes, with {1} guides/gene. {2} runs will be performed.".format(num_genes, library_size, num_runs)
  t0 = time.time()
  # Stores the averages of normalized expressions across different sets of genes
  results_Y = []
  results_N = []

  # open the list of genes
  with open('../pre_processed/genes_list.json', 'r') as genes_list_file:
    genes_list = json.load(genes_list_file)

    for i in xrange(num_runs):
      genes = random.sample(genes_list, num_genes)
      result_Y, result_N = gtex_effect(genes)
      results_Y.append(result_Y)
      results_N.append(result_N)

  with open('gtex_assists_exon_selection.results', 'w') as results:
    results.write("Average of normalized exon expression across {0} genes, with {1} guides/gene.\nData from {2} runs is summarized below.\n\n".format(num_genes, library_size, num_runs))
    results.write("Guide selection with GTEX:\n")
    for result in results_Y:
      results.write(str(result) + "\n")
    results.write("\n\nGuide selection without GTEX:\n")
    for result in results_N:
      results.write(str(result) + "\n")

  print "Completed job in", time.time() - t0, "seconds."
