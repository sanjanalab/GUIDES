# produce list of genes in GRCm38
import pandas as pd
import json

refGeneFilename = '../gtex/refGene_GRCm38.txt'
refGene = pd.read_csv(refGeneFilename, sep="\t")
refGene.columns=['','name','chrom','strand','txStart','txEnd','cdsStart','cdsEnd','exonCount','exonStarts','exonEnds','id','name2','cdsStartStat','cdsEndStat','exonFrames']

seen = {}
results = []
for index, row in refGene.iterrows():
  val = row['name2']
  ensembl_id = row['name']
  if val not in seen:
    results.append({
      'name': val,
      'ensembl_id': ensembl_id
    })
    seen[val] = True

with open('genes_list_GRCm38.txt', 'w') as output:
  json.dump(results, output)
