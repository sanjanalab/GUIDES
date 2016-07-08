# Goal: create a cPickle'd python dictionary
# containing the strand of each ensembl gene
import csv
import cPickle

if __name__ == "__main__":
  mapping = {}
  with open('../gtex/refGene.txt', 'r') as csvfile:
    reader = csv.reader(csvfile, delimiter='\t')
    for row in reader:
      ensembl_gene = row[1]
      strand = row[3]
      mapping[ensembl_gene] = strand

  cPickle.dump(mapping, open("strand_info.p", "wb"))
