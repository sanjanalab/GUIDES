# produce mouse genome exon by exon
import gzip
from Bio import SeqIO
import pickle
import os
import pandas as pd

chromosomes = {}
filename_base = os.path.join(os.path.dirname(__file__), "static/data/GRCm38/chr{}.fa.gz")

print("starting refGene load")
refGeneFilename = os.path.join(os.path.dirname(__file__), "static/data/gtex/refGene_GRCm38.txt")
refGene = pd.read_csv(refGeneFilename, sep="\t")
refGene.columns=['','name','chrom','strand','txStart','txEnd','cdsStart','cdsEnd','exonCount','exonStarts','exonEnds','id','name2','cdsStartStat','cdsEndStat','exonFrames']
refGene["exonStarts"] = refGene.apply(lambda x: x['exonStarts'].split(',')[:-1], axis=1)
refGene["exonEnds"] = refGene.apply(lambda x: x['exonEnds'].split(',')[:-1], axis=1)
refGene["exonFrames"] = refGene.apply(lambda x: x['exonFrames'].split(',')[:-1], axis=1)
print("completed refGene load")

chroms = ["1","2","3","4","5","6","7","8","9","10","11","12","12","13","14","15","15","16","17","18","19","M","X","Y"]

print("starting read in all chrom")
for c in chroms:
  filename = filename_base.format(c)
  handle = gzip.open(filename)
  chromosomes[c] = SeqIO.read(handle, "fasta")
print("finished read in all chrom")

print("processing each gene")
for index, row in refGene.iterrows():
  gene_name, exon_count = row["name2"], row["exonCount"]
  exonStarts = row["exonStarts"]
  exonEnds = row["exonEnds"]
  chrom = row["chrom"].split('chr')[1]
  if not chrom in chroms:
    continue
  full_sequence = chromosomes[chrom]
  for exon in range(int(exon_count)):
    start = int(exonStarts[exon])
    end = int(exonEnds[exon])
    sequence = full_sequence.seq[start:end]
    filename = "{0}_{1}".format(gene_name, exon)
    path = os.path.join('static/data/GRCm38_exons/', filename)
    print filename
    with open(path, 'w') as outfile:
      outfile.write(str(sequence))
print("finished generating GRCm38 exons")
