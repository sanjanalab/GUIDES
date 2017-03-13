# produce list of genes in GRCm38
import pandas as pd
import json

# open refgene
refGeneFilename = '../gtex/gtex_mouse/refGene_mouse.txt'
refGene = pd.read_csv(refGeneFilename, sep="\t")
refGene.columns=['','name','chrom','strand','txStart','txEnd','cdsStart','cdsEnd','exonCount','exonStarts','exonEnds','id','name2','cdsStartStat','cdsEndStat','exonFrames']

# open biomart
biomartFilename = 'mart_export_mus_2.txt'
biomart = pd.read_csv(biomartFilename, sep="\t")

seen = {}
results = []
total_len = len(refGene)
for index, row in refGene.iterrows():
  ensembl_id = row['name']
  if ensembl_id not in seen:

    the_loc = biomart.loc[biomart['Gene ID'] == ensembl_id]
    gene_name = list(the_loc['Associated Gene Name'])[0]
    entrez = list(the_loc['EntrezGene ID'])[0]
    if pd.isnull(entrez):
      entrez = ''
      print ensembl_id, gene_name, 'has no entrez'
    else:
      entrez = str(int(entrez))

    if pd.isnull(gene_name):
      gene_name = ''
      print ensembl_id, 'has no gene_name'

    results.append({
      'name': gene_name,
      'ensembl_id': ensembl_id,
      'entrez_id': entrez,
      'description': ""
    })
    seen[ensembl_id] = True

with open('genes_list_GRCm38_processed.txt', 'w') as output:
  json.dump(results, output)

with open('genes_list_GRCm38.txt', 'w') as output:
  json.dump(results, output)
