import cPickle
import os

def overall_gene_expression(ensembl_gene):
  df_gene_filename = ensembl_gene + ".p"
  df_gene_path = os.path.join('../df_gene_normalized_cpickle', df_gene_filename)
  df_gene_file = open(df_gene_path, 'r')
  df_gene = cPickle.load(df_gene_file)

def gtex_effect(gene_names):
  sum_normalized_expression_1 = 0
  sum_normalized_expression_2 = 0

  for gene_name in gene_names:
