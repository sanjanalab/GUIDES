# Goal: To precompute gene expression information for each exon using cPickle for easy access later.
# NO MORE PD DATAFRAMES. They are slow. Let's optimize on C.
import pickle
import msgpack
import json
import math

df_normalized = pickle.load(open('pd_by_tissue_normalized.p', "rb"))
tissues = df_normalized.columns[1:]

with open('genes_list.json') as genes_list_file:
  genes_list = json.load(genes_list_file)
   # gene format: {"ensembl_id": "ENSG00000261122.2", "name": "5S_rRNA", "description": ""}
  for gene in genes_list:
    exon = 0
    df_gene = df_normalized[df_normalized.Id.str.contains(ensembl_gene + '_' + exon)]
    # We're done here...
    if df_gene.shape[0] == 0:
      exon = 0
      continue

    # Make an associative array from tissues -> expression values
    expression = {}
    exon = df_gene.iloc[0]
    for t in tissues:
      expression_val = exon[t]
      if math.isnan(expression_val):
        expression_val = 0
      expression[t] = expression_val
      
    outfile_name = gene["ensembl_id"] + "_" + str(exon) + ".p"
    output_path = os.path.join('../exon_expression_msgpack/', outfile_name)
    with open(output_path,  'w') as outfile:
      msgpack.dump(gRNAs, outfile)
    



    seq = gene_exon_file(gene["ensembl_id"], exon)



import time

i = 0
def test_pandas_timing():
  t0 = time.time()
  with open('genes_list.json') as genes_list_file:
    genes_list = json.load(genes_list_file)
    for gene in genes_list[:100]:
      df_gene = df_normalized[df_normalized.Id.str.contains(gene['ensembl_id'] + '_0')]
  print time.time() - t0
  print (time.time() - t0) / 100

def test_guides_timing():
  t0 = time.time()
  with open('genes_list.json') as genes_list_file:
    genes_list = json.load(genes_list_file)
    for gene in genes_list[:100]:
      dfff_gene = getGuides(gene['ensembl_id'] + '_0')
  print time.time() - t0
  print (time.time() - t0) / 100


def getGuides(gene_exon):
  try:
    filename = gene_exon + ".p"
    path = os.path.join('../GRCh37_guides_cpickle/', filename)
    with open(path) as datafile:
      gRNAs = cPickle.load(datafile)
      return gRNAs
  except IOError:
    pass


