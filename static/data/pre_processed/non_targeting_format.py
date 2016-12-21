import json

output = {'data': []}

def main(inp_filename, out_filename):
  with open(inp_filename) as inp:
    for idx,line in enumerate(inp):
      output['data'].append({
        'seq': line.strip(),
        'score': '',
        'gene': '',
        'exon': '',
        'functional_domain': '',
        'uid': 'GUIDES_nt_sg' + str(idx).zfill(4)
      })

  # write data
  with open(out_filename, 'w+') as out:
    out.write(json.dumps(output))

main('non_targeting_unformatted_hum.txt', 'non_targeting_hum.json')
main('non_targeting_unformatted_mus.txt', 'non_targeting_mus.json')
