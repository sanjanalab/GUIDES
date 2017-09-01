# produce list of genes in GRCm38
import pandas as pd
import json

# open refgene
refGeneFilename = '../gtex/refGene.txt'
refGene = pd.read_csv(refGeneFilename, sep="\t")
refGene.columns=['','name','chrom','strand','txStart','txEnd','cdsStart','cdsEnd','exonCount','exonStarts','exonEnds','id','name2','cdsStartStat','cdsEndStat','exonFrames']

# open biomart
biomartFilename = 'mart_export_hum_2.txt'
biomart = pd.read_csv(biomartFilename, sep="\t")

results = []
total_len = len(refGene)
for index, row in refGene.iterrows():
  full_ensembl_id = row['name']
  ensembl_id = full_ensembl_id.split('.')[0]
  the_loc = biomart.loc[biomart['Gene ID'] == ensembl_id]
  gene_name_l = list(the_loc['Associated Gene Name'])
  entrez_l = list(the_loc['EntrezGene ID'])

  if len(gene_name_l) != len(entrez_l):
    print "ASSERT: gene_name_l and entrez_l should have equal lengths"

  if len(gene_name_l) == 0 or len(entrez_l) == 0:
    print "no gene name and no entrez id"
    results.append({
      'name': '',
      'ensembl_id': full_ensembl_id,
      'entrez_id': '',
      'description': ""
    })

  for gene_name, entrez in zip(gene_name_l, entrez_l):
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
      'ensembl_id': full_ensembl_id,
      'entrez_id': entrez,
      'description': ""
    })

with open('genes_list_GRCh37_processed.txt', 'w') as output:
  json.dump(results, output)

with open('genes_list.json', 'w') as output:
  json.dump(results, output)
