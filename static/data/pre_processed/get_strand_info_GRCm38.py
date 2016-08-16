# Goal: create a cPickle'd python dictionary
# containing the strand of each ensembl gene

#m mouse
import csv
import cPickle

if __name__ == "__main__":
  mapping = {}
  with open('../gtex/refGene_GRCm38.txt', 'r') as csvfile:
    reader = csv.reader(csvfile, delimiter='\t')
    for row in reader:
      gene_name = row[12]
      print gene_name
      strand = row[3]
      mapping[gene_name] = strand

  cPickle.dump(mapping, open("strand_info_GRCm38.p", "wb"))
