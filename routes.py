from flask import session, request, render_template, jsonify, g, redirect
from flask.ext.basicauth import BasicAuth
from app import app, dev
from flask import jsonify, url_for
import pandas as pd
import numpy as np
import pickle
import urllib2
import seq_generator
import computations
from settings import APP_STATIC
import os
import json
import time

genome = {
  "human" : seq_generator.FastGenome()
}

# setup password protection
app.config['BASIC_AUTH_USERNAME'] = 'zhanglab'
app.config['BASIC_AUTH_PASSWORD'] = 'editMe23'
app.config['BASIC_AUTH_FORCE'] = True # protect entire site
# for individual routes, decorate with @basic_auth.required
app.config['PROPOGATE_EXCEPTIONS'] = True

basic_auth = BasicAuth(app)

### Processing
@app.before_request
def preprocess_request():
  pass

@app.after_request
def postprocess_request(response):
  return response

### pass routing to angular
@app.route('/')
@app.route('/designer')
@app.route('/about')
def index_view():
  return render_template('index.html')

### API
@app.route('/generate', methods=['POST'])
def generate():
  t0 = time.time()

  # Retrieve arguments
  post = request.get_json()

  genes = post.get('genes')
  species = 'human'
  quantity = post.get('quantity')
  tissues = post.get('tissues')
  gtex_enabled = post.get('gtex_enabled')

  # Validations
  if genes == None:
    genes = [{u'ensembl_id': u'ENSG00000186575.13', u'$$hashKey': u'object:50', u'name': u'NF2'}]
  if species == None:
    species = "human"
  if quantity == None:
    quantity = 6
  if tissues == None:
    tissues = ["Muscle", "Heart", "Brain"]

  if len(genes) == 0:
    return jsonify(gene_to_exon={}, guide_count=0)

  # Setup ranker
  tissues_enabled = False if len(tissues) == 31 else True # true, unless all tissues are selected (average all)
  ranker = computations.Ranker(genome["human"], species, tissues, gtex_enabled, tissues_enabled)

  # Iterate over genes, finding guides for each
  for g in genes:
    ranker.rank(g['ensembl_id'], g['name'], quantity)

  guides_by_exon = ranker.get_guides_by_exon()
  guide_count = ranker.get_count_selected_guides()
  result = jsonify(gene_to_exon = guides_by_exon, guide_count = guide_count)

  print "Spent {0} seconds generating {1} guides/gene for {2} genes.".format(time.time() - t0, quantity, len(genes))

  return result

def test_generate():
  # Generate data
  data = [
    {
      "name": "SOD1",
      "exons": [
          {
              "id": 0,
              "length": 121,
              "gRNAs": 1
          },
          {
              "id": 1,
              "length": 22,
              "gRNAs": 3
          },
          {
              "id": 2,
              "length": 22,
              "gRNAs": 3
          },
          {
              "id": 3,
              "length": 22,
              "gRNAs": 3
          },
          {
              "id": 4,
              "length": 22,
              "gRNAs": 3
          },
          {
              "id": 5,
              "length": 22,
              "gRNAs": 3
          },
          {
              "id": 6,
              "length": 22,
              "gRNAs": 3
          },
          {
              "id": 7,
              "length": 22,
              "gRNAs": 3
          },
          {
              "id": 8,
              "length": 22,
              "gRNAs": 3
          },
          {
              "id": 9,
              "length": 22,
              "gRNAs": 3
          },
          {
              "id": 10,
              "length": 22,
              "gRNAs": 3
          },
          {
              "id": 11,
              "length": 22,
              "gRNAs": 3
          },
          {
              "id": 12,
              "length": 22,
              "gRNAs": 3
          }
      ],
      "grnas": [
          {
              "id": 0,
              "exon": 5,
              "score": 0.925,
              "selected": True
          },
          {
              "id": 1,
              "exon": 9,
              "score": 0.654,
              "selected": True
          }
      ]
    },
    {
      "name": "AAK1",
      "exons": [
          {
              "id": 0,
              "length": 230,
              "gRNAs": 0
          },
          {
              "id": 1,
              "length": 50,
              "gRNAs": 1
          },
          {
              "id": 2,
              "length": 23,
              "gRNAs": 3
          }
      ],
      "grnas": [
          {
              "id": 0,
              "exon": 2,
              "score": 0.933,
              "selected": False
          },
          {
              "id": 1,
              "exon": 4,
              "score": 0.914,
              "selected": True
          },
          {
              "id": 2,
              "exon": 7,
              "score": 0.724,
              "selected": False
          }
      ]
    }
  ]
  gene_to_exon = [
      {
        "name": "AAK1",
        "length": 1300,
        "exons": [
          {
            "start": 150,
            "end": 300,
            "gRNAs": [
              {
                "score": 0.914,
                "start": 20, # from start of exon
                "seq": "ATGTCGATGAGATCGATGAGCGG"
              },
              {
                "score": 0.714,
                "start": 50, # from start of exon
                "seq": "TTGCCGAGGAGATCGATGAGCGG"
              }
            ]
          },
          {
            "start": 600,
            "end": 750,
            "gRNAs": [
              {
                "score": 0.972,
                "start": 10, # from start of exon
                "seq": "GGGTCGATCCCATCGATGAGCGG"
              },
              {
                "score": 0.814,
                "start": 90, # from start of exon
                "seq": "AAGCCGAGGAGATCGATGAGCGG"
              },
              {
                "score": 0.753,
                "start": 120,
                "seq": "TTTCCGAGGAGGTCGATGAGCGG"
              }
            ]
          },
          {
            "start": 800,
            "end": 900,
            "gRNAs": [
              {
                "score": 0.442,
                "start": 3, # from start of exon
                "seq": "GGGTAAATCCCATCGATGAGCGG"
              },
              {
                "score": 0.214,
                "start": 27, # from start of exon
                "seq": "AAGCCGAGGTTGTCGATGAGCGG"
              },
              {
                "score": 0.153,
                "start": 76,
                "seq": "GGTCCGAGGAGGTCGATGAGCGG"
              }
            ]
          }
        ]
      },
      {
        "name": "SOD1",
        "length": 2000,
        "exons": [
          {
            "start": 600,
            "end": 750,
            "gRNAs": [
              {
                "score": 0.972,
                "start": 10, # from start of exon
                "seq": "GGGTCGATCCCATCGATGAGCGG"
              },
              {
                "score": 0.814,
                "start": 90, # from start of exon
                "seq": "AAGCCGAGGAGATCGATGAGCGG"
              },
              {
                "score": 0.753,
                "start": 120,
                "seq": "TTTCCGAGGAGGTCGATGAGCGG"
              }
            ]
          },
          {
            "start": 800,
            "end": 900,
            "gRNAs": [
              {
                "score": 0.442,
                "start": 3, # from start of exon
                "seq": "GGGTAAATCCCATCGATGAGCGG"
              },
              {
                "score": 0.214,
                "start": 27, # from start of exon
                "seq": "AAGCCGAGGTTGTCGATGAGCGG"
              },
              {
                "score": 0.153,
                "start": 76,
                "seq": "GGTCCGAGGAGGTCGATGAGCGG"
              }
            ]
          }
        ]
      }
    ]
  # Send json to client
  guide_count = 6
  return jsonify(genes=data, gene_to_exon=gene_to_exon, guide_count=guide_count)

### Filters

@app.template_filter('currency')
def currency_filter(val):
  return "${:,.2f}".format(val)
