import json

def main(inp_filename, out_filename, base_uid):
  output = {'data': []}

  with open(inp_filename) as inp:
    for idx,line in enumerate(inp):
      sequence = line.strip()
      output['data'].append({
        'seq': sequence,
        'score': '',
        'gene': '',
        'exon': '',
        'functional_domain': '',
        'uid': base_uid + str(idx + 1).zfill(4),
        'scaffold_standard': "GGAAAGGACGAAACACCG" + sequence + "GTTTTAGAGCTAGAAATAGCAAGTTAAAATAAGGC",
        'scaffold_EF': "GGAAAGGACGAAACACCG" + sequence + "GTTTAAGAGCTATGCTGGAAACAGC"
      })

  # write data
  with open(out_filename, 'w+') as out:
    out.write(json.dumps(output))

main('non_targeting_unformatted_hum.txt', 'non_targeting_hum.json', 'NonTargeting_Human_')
main('non_targeting_unformatted_mus.txt', 'non_targeting_mus.json', 'NonTargeting_Mouse_')
