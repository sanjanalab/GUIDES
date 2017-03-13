import json

with open('genes_list_GRCm38.txt', 'r') as genes_list:
  data = json.load(genes_list)
  for i in range(len(data)):
    data[i]['mRNA'] = data[i]['ensembl_id']
    data[i]['ensembl_id'] = data[i]['ensembl_id_real']
    del data[i]['ensembl_id_real']

with open('genes_list_GRCm38_processed.txt', 'w') as outfile:
  json.dump(data, outfile)
