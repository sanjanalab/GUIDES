# Automate speed benchmarking
import json
import random
import requests

# Prepare random lists of genes
list_lens = [1,10,50,100,300,500,700,900,1000,3000,5000,7000,9000,10000,12000,14000,17000,20000]
quantities = [3,6,10,20]
n_trials = 3

# constant
all_tissues = ['Thyroid', 'Testis', 'Cervix Uteri', 'Adipose Tissue', 'Breast', 'Vagina', 'Nerve', 'Pituitary', 'Stomach', 'Fallopian Tube', 'Bone Marrow', 'Bladder', 'Blood', 'Colon', 'Prostate', 'Pancreas', 'Blood Vessel', 'Liver', 'Spleen', 'Small Intestine', 'Uterus', 'Ovary', 'Muscle', 'Heart', 'Adrenal Gland', 'Brain', 'Salivary Gland', 'Lung', 'Skin', 'Esophagus', 'Kidney']

if __name__ == "__main__":
  with open('../pre_processed/genes_list.json', 'r') as genes_list_file:
      genes_list = json.load(genes_list_file)
      for n in xrange(n_trials):
        for l in list_lens:
          genes = random.sample(genes_list, l)

          url = 'http://localhost:2500/generate'
          for quantity in quantities:
            data = {
              'genes': genes,
              'quantity': quantity,
              'tissues': all_tissues,
              'gtex_enabled': True
            }
            data_json = json.dumps(data)
            headers = {'Content-type': 'application/json'}
            response = requests.post(url, data=data_json, headers=headers) # blocking
