# Python packages
import json

# CLD code from parent directories
import os
import inspect
currentdir = os.path.dirname(os.path.abspath(inspect.getfile(inspect.currentframe())))
rootdir = os.path.dirname(os.path.dirname(os.path.dirname(currentdir)))
os.sys.path.insert(0,rootdir)
import computations
import seq_generator

### Constant ranker parameters
genome = {
  "human" : seq_generator.FastGenome()
}
tissues = ['Thyroid', 'Testis', 'Cervix Uteri', 'Adipose Tissue', 'Breast', 'Vagina', 'Nerve', 'Pituitary', 'Stomach', 'Fallopian Tube', 'Bone Marrow', 'Bladder', 'Blood', 'Colon', 'Prostate', 'Pancreas', 'Blood Vessel', 'Liver', 'Spleen', 'Small Intestine', 'Uterus', 'Ovary', 'Muscle', 'Heart', 'Adrenal Gland', 'Brain', 'Salivary Gland', 'Lung', 'Skin', 'Esophagus', 'Kidney']

#tissues = ['Lung', 'Liver', 'Skin', 'Heart', 'Brain']

### Experimental parameters
scoring_alg = "Doench"
library_size = 6

if __name__ == "__main__":
  print "Computing number of exons targeted for all genes, using {0} tissues, {1}".format(len(tissues), scoring_alg)

  filename = 'num_exons_{0}tissues_{1}_pGTEX.csv'.format(len(tissues), scoring_alg)
  with open(filename, 'w') as results:
    results.write("gene\tensembl\texons_targeted\n")

    with open('../pre_processed/genes_list.json', 'r') as genes_list_file:
      genes_list = json.load(genes_list_file)
      for gene in genes_list:
        ranker = computations.Ranker(genome["human"], "human", tissues, True, True, scoring_alg = scoring_alg)
        ranker.rank(gene['ensembl_id'], gene['name'], library_size)
        guides_by_exon = ranker.get_guides_by_exon()
        num_exons_targeted = 0
        for gene_guide in guides_by_exon:
          for idx, exon in enumerate(gene_guide['exons']):
            for guide in exon['gRNAs']:
              if guide['selected']:
                num_exons_targeted += 1
                break
        results.write("{0}\t{1}\t{2}\n".format(gene['name'], gene['ensembl_id'], num_exons_targeted))

