import json
import random

# Prepare random lists of genes
list_lens = [1,10,50,100,300,500,700,900,1000,3000,5000,7000,9000,10000,12000,14000,17000,20000]
n_trials = 3

if __name__ == "__main__":
  with open('../pre_processed/genes_list.json', 'r') as genes_list_file:
      genes_list = json.load(genes_list_file)
      for n in xrange(n_trials):
        for l in list_lens:
          genes = random.sample(genes_list, l)
          gene_names = [g['name'] for g in genes]
          filename = "speed_benchmarking/genes_{0}_{1}.csv".format(n, l)
          with open(filename, 'w') as results:
            results.write('\t'.join(gene_names))
