'''
Write the sequence of each exon to a seperate text file.
'''

from seq_generator_mus import *
import os

g = Genome()
for gene, exon_count in g.genes_exons():
	for exon in range(int(exon_count)):
		seq = g.sequence(gene, exon)
		filename = "{0}_{1}".format(gene, exon)
		path = os.path.join('static/data/GRCm38_exons/', filename)
		print "writing data to " + path
		with open(path, 'w') as outfile:
			outfile.write(str(seq))
