# Python packages
import cPickle
import json
import random
import time
import msgpack
from os import listdir
import os.path

### Experimental parameters
num_runs = 10000
method = "tissues5" # tissues5 or overall
scoring_alg = "Doench"

def guidesFromFile(filename):
  with open(filename) as datafile:
    gRNAs = msgpack.load(datafile)
    return gRNAs

def getExpression(gene, exon):
  df_gene_filename = gene + ".p"
  df_gene_path = os.path.join(os.path.dirname(__file__), '../df_gene_normalized_cpickle', df_gene_filename)
  with open(df_gene_path, 'r') as df_gene_file:
    df_gene = cPickle.load(df_gene_file)
    df_exon = df_gene[df_gene.exon_num == exon]
    if method == "overall":
      result = df_exon['overall']
      return float(result)
    elif method == "tissues5":
      df_exon['median'] = df_exon[['Brain', 'Heart', 'Kidney', 'Liver', 'Skin']].median(axis=1)
      result = df_exon['median']
      return float(result)

if __name__ == "__main__":
  print "Computing on-target efficiency vs. the GTEx exon constitutiveness score for {0} guides.".format(num_runs)

  filename = 'expression_on_target_correlation.' + str(num_runs) + '.' + method + '.' + scoring_alg + '.results.csv'
  with open(filename, 'w') as results:
    results.write("score\texpression\n")
    files_all = os.listdir(os.path.join(os.path.dirname(__file__), '../GRCh37_guides_msgpack' + '_' + scoring_alg + '/'))
    files_test = [random.choice(files_all) for i in xrange(num_runs)]
    for filename in files_test:
      gene, exon = filename[:-2].split('_')
      path = os.path.join(os.path.dirname(__file__), '../GRCh37_guides_msgpack' + '_' + scoring_alg + '/', filename)
      guides = guidesFromFile(path)
      if len(guides) == 0:
        print "no guides for", filename
        continue
      guide = random.choice(guides)

      score = guide['score']
      if score == 0:
        continue
      expression = getExpression(gene, int(exon))
      results.write("{0}\t{1}\n".format(score, expression))
