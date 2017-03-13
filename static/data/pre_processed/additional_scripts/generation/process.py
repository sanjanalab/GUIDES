# number of exons for gene
import msgpack
import json
import os
import csv

#APP_STATIC = "/Users/Josh/Documents/Research-Zhang-Lab/CRISPR_Library_Design_2/static"
APP_STATIC = "/home/joshm/GUIDES/CRISPR-Library-Designer/static"

class LargeScreen():
  def __init__(self):
    # human
    with open(os.path.join(APP_STATIC, 'data/pre_processed', 'genes_list.json'), "r") as infile:
      self.hum_genes = json.loads(infile.read())
    # self.hum_hash_map = {}
    # for gene in hum_genes:
    #   self.hum_hash_map[gene['name'].lower()] = gene

    # ensembl_id (ENSG), name, description

    # mouse
    with open(os.path.join(APP_STATIC, 'data/pre_processed', 'genes_list_GRCm38.txt'), "r") as infile:
      self.mus_genes = json.loads(infile.read())
      # self.mus_hash_map = {}
      # for gene in mus_genes:
      #   self.mus_hash_map[gene['name'].lower()] = gene

    # ensembl_id, name

    # read in Patt's list
    with open('genes_for_12.txt', 'r') as infile:
      self.patt_hash_map = {}
      for line in infile:
        self.patt_hash_map[line.strip()] = True

  def getGuides(self, gene, exon, species):
    gene_exon = gene + "_" + str(exon)
    try:
      filename = gene_exon + ".p"
      if species == 'hum':
        path = os.path.join(APP_STATIC, 'data/GRCh37_guides_msgpack_Azimuth_domains/', filename)
      elif species == 'mus':
        path = os.path.join(APP_STATIC, 'data/GRCm38_guides_msgpack_Azimuth_domains/', filename)
      with open(path) as datafile:
        gRNAs = msgpack.load(datafile)
        return gRNAs
    except IOError:
      raise IOError

  def getGuidesArray(self, gene, species):
    result = []
    i = 0
    while 1:
      try:
        result.append(self.getGuides(gene, i, species))
        i += 1
      except IOError:
        break

    if len(result) > 6:
      del result[0]
      del result[-1]

    return result

  def functional_yes(self, guide):
    if "functional_domain" in guide:
      return (1, guide["score"], guide)
    else:
      return (0, guide["score"], guide)

  def functional_no(self, guide):
    if "functional_domain" in guide:
      return (0, guide["score"], guide)
    else:
      return (1, guide["score"], guide)

  def getFinalGuides(self, gene_obj, gene, twelver, species):
    def pull_values(guide):
      if 'functional_domain' not in guide:
        guide['functional_domain'] = ''
      return [guide['seq'], guide['PAM'], guide['score'], guide['functional_domain'], guide['start'], gene_obj['name'], gene_obj['ensembl_id'], species, twelver]

    guides_array = self.getGuidesArray(gene, species)
    combined_guides_array = [item for sublist in guides_array for item in sublist]

    combined_guides_array.sort(key=self.functional_yes, reverse=True)
    if not twelver:
      return [pull_values(guide) for guide in combined_guides_array[:6]]
    else:
      result = combined_guides_array[:6]
      combined_guides_array.sort(key=self.functional_no, reverse=True)
      result += combined_guides_array[:6]

      seen_seqs = {}
      filtered_result = []
      for guide in result:
        if guide['seq'] in seen_seqs:
          continue
        filtered_result.append(guide)
        seen_seqs[guide['seq']] = True

      return [pull_values(guide) for guide in filtered_result]

  def process(self):
    guides = []
    for gene in self.hum_genes:
      if gene['name'].upper() in self.patt_hash_map:
        guides += self.getFinalGuides(gene, gene['ensembl_id'], True, 'hum')
      else:
        guides += self.getFinalGuides(gene, gene['ensembl_id'], False, 'hum')

    for gene in self.mus_genes:
      if gene['name'].upper() in self.patt_hash_map:
        guides += self.getFinalGuides(gene, gene['name'], True, 'mus')
      else:
        guides += self.getFinalGuides(gene, gene['name'], False, 'mus')

    return guides

if __name__ == "__main__":
  s = LargeScreen()
  guides = s.process()
  i = 0

  with open('large_patt_results.csv', 'w') as csvfile:
    writer = csv.writer(csvfile)
    writer.writerows(guides)
