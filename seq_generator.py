import gzip
from Bio import SeqIO
import pickle
import os

class ExonError(Exception):
    """Error class for improper exon numbering"""
    def __init__(self, gene, exonNum):
        self.gene = gene
        self.exonNum = exonNum
    def __str__(self):
        return "Exon not found for exon %s at %s" % (self.exonNum, self.gene)

class Genome():
    """Tools for accessing sequences, coordinates, and exons."""
    def __init__(self, preload = False):
        self.chromosomes = {}
        self.filename_base = os.path.join(os.path.dirname(__file__), "static/data/GRCh37/Homo_sapiens.GRCh37.75.dna.chromosome.{}.fa.gz")
        with open(os.path.join(os.path.dirname(__file__), "static/data/pre_processed/exon_info.p"), "rb") as f:
            self.df = pickle.load(f)

        if preload:
            chroms = ["1","2","3","4","5","6","7","8","9","10","11","12","12","13","14","15","15","16","17","18","19","20","21","22","MT","X","Y"]
            for c in chroms:
                self.chrom_sequence(c)

    def chrom_sequence(self, c):
        if c not in self.chromosomes.keys():
            filename = self.filename_base.format(c)
            handle = gzip.open(filename)
            self.chromosomes[c] = SeqIO.read(handle, "fasta")
        return self.chromosomes[c]

    def count_exons_in_gene(self, gene):
        gene_data = self.df.loc[self.df['name'] == gene]
        return len(gene_data["exonStarts"].tolist()[0])

    def genes_exons(self):
        return [(row["name"], row["exonCount"]) for index, row in self.df.iterrows()]

    def sequence(self, gene, exon):
        try:
            gene_data = self.df.loc[self.df['name'] == gene]
            chrom = str(gene_data["chrom"].tolist()[0]).split('.0')[0]
            start = int(gene_data["exonStarts"].tolist()[0][exon])
            end = int(gene_data["exonEnds"].tolist()[0][exon])
            full_sequenece = self.chrom_sequence(chrom)
            full_sequenece_seq = self.chrom_sequence(chrom).seq
            return full_sequenece_seq[start:end].upper()
        except IndexError:
            raise ExonError(gene, exon)

    # gene must be in ensembl format
    def gene_info(self, gene):
        """Wrap information regarding a certain gene."""
        gene_data = self.df.loc[self.df['name'] == gene]
        gene_info = {
            'txStart': 0, #tx - transcription
            'txEnd': 0,
            'exonCount': int(gene_data['exonCount']),
            'exonStarts': map(int, gene_data['exonStarts'].tolist()[0]),
            'exonEnds': map(int, gene_data['exonEnds'].tolist()[0])
        }
        return gene_info

    def sequence_gtex_gene(self, gene):
        """Retrieve DNA sequence for gene in format """
        try:
            values = gene.split('_')
            return self.sequence(values[0], int(values[1]))
        except IndexError:
            raise ValueError('Gene string entered in wrong format.')

class FastGenome(Genome):
    def sequence(self, gene, exon):
        try:
            filename = "{0}_{1}".format(gene, exon)
            path = os.path.join(os.path.join(os.path.dirname(__file__),'static/data/GRCh37_exons/'), filename)
            with open(path) as infile:
                return infile.read()
        except IOError:
            raise ExonError(gene, exon)
