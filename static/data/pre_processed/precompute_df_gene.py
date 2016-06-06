# goal: precompute df_genes

import pickle
import cPickle
import json
import math
import os

df_normalized = pickle.load(open('pd_by_tissue_normalized.p', "rb"))

with open('genes_list.json') as genes_list_file:
  genes_list = json.load(genes_list_file)
   # gene format: {"ensembl_id": "ENSG00000261122.2", "name": "5S_rRNA", "description": ""}
  for gene in genes_list:
    df_gene = df_normalized[df_normalized.Id.str.contains(gene["ensembl_id"])]

    # Get median expression for all tissues
    df_gene['overall'] = df_gene.median(axis=1)
    df_gene['median'] = df_gene.median(axis=1)

    # Recall: Id are entered as ENSG000xxxxx.x_EXONNUM, e.g. ENSG00000000971.11_21
    total_exons = len(df_gene)
    df_gene['exon_num'] = df_gene['Id'].apply(lambda x: int(x.split('_')[1]))

    # Sort so that we don't have to for experiments where we use all tissues.
    df_gene = df_gene.sort(['median'], ascending=False)

    outfile_name = gene["ensembl_id"] + ".p"
    output_path = os.path.join('../df_gene_normalized_cpickle/', outfile_name)
    with open(output_path,  'w') as outfile:
      cPickle.dump(df_gene, outfile)
