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
import computations_mouse
from settings import APP_STATIC
import os
import json
import time
from celery import Celery
from emailing import send_completed_run

genome = {
  "human" : seq_generator.FastGenome()
}

# setup password protection
# app.config['BASIC_AUTH_USERNAME'] = 'zhanglab'
# app.config['BASIC_AUTH_PASSWORD'] = 'editMe23'
# app.config['BASIC_AUTH_FORCE'] = True # protect entire site
# for individual routes, decorate with @basic_auth.required
app.config['PROPOGATE_EXCEPTIONS'] = True

# Celery configuration
app.config['CELERY_BROKER_URL'] = 'redis://localhost:6379/0'
app.config['CELERY_RESULT_BACKEND'] = 'redis://localhost:6379/0'

# Initialize Celery
celery = Celery(app.name, broker=app.config['CELERY_BROKER_URL'])
celery.conf.update(app.config)

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
@celery.task(bind=True)
def start_compute_mouse(self, params):
  with app.app_context():
    genes = params['genes']
    species = params['species']
    quantity = params['quantity']
    domains_enabled = params['domains_enabled']
    email_address = params['email_address']
    rejected_genes = params['rejected_genes']
    gene_statistics = params['gene_statistics']

    if len(genes) == 0:
      result = {'gene_to_exon': {}, 'guide_count': 0, 'genome': 'mouse'}
      return {'current': 0, 'total': 0, 'status': 'Task completed!', 'result': result, 'gene_statistics': gene_statistics}

    ranker = computations_mouse.RankerMouse(domains_enabled)

    # Iterate over genes finding guides for each
    total_gene_count = len(genes)
    for idx, g in enumerate(genes):
      ranker.rank(g['name'], quantity)
      gene_statistics['processed'] = idx
      self.update_state(state='PROGRESS', meta={'current': idx, 'total': total_gene_count, 'gene_statistics': gene_statistics})

    guides_by_exon = ranker.get_guides_by_exon()
    guide_count = ranker.get_count_selected_guides()

    result = {
      'gene_to_exon': guides_by_exon,
      'guide_count': guide_count,
      'genes': genes,
      'quantity': quantity,
      'domains_enabled': domains_enabled,
      'email_address': email_address,
      'rejected_genes': rejected_genes,
      'genome': 'mus',
      'gtex_enabled': False,
      'tissues_disabled': True,
      'tissues': []
    }

    # Send results to the user
    send_completed_run(email_address, start_compute.request.id)

    return {
      'current': total_gene_count,
      'total': total_gene_count,
      'status': 'Task completed!',
      'gene_statistics': gene_statistics,
      'result': result,
      'email_address': email_address,
      'genome': 'mus',
      'gtex_enabled': False,
      'tissues_disabled': True,
      'tissues': []
    }

@celery.task(bind=True)
def start_compute(self, params):
  with app.app_context():
    t0 = time.time()

    genes = params['genes']
    species = params['species']
    quantity = params['quantity']
    tissues = params['tissues']
    gtex_enabled = params['gtex_enabled']
    tissues_disabled = params['tissues_disabled']
    domains_enabled = params['domains_enabled']
    gene_statistics = params['gene_statistics']
    email_address = params['email_address']
    rejected_genes = params['rejected_genes']

    if len(genes) == 0:
      result = {'gene_to_exon': {}, 'guide_count': 0, 'genome': 'human'}
      return {'current': 0, 'total': 0, 'status': 'Task completed!', 'result': result, 'gene_statistics': gene_statistics}

    # Setup ranker
    tissues_enabled = False if len(tissues) == 31 else True # true, unless all tissues are selected (average all)
    ranker = computations.Ranker(genome["human"], species, tissues, gtex_enabled, tissues_enabled, domains_enabled)

    # Iterate over genes, finding guides for each
    total_gene_count = len(genes)
    for idx, g in enumerate(genes):
      ranker.rank(g['ensembl_id'], g['name'], quantity)
      gene_statistics['processed'] = idx
      self.update_state(state='PROGRESS', meta={'current': idx, 'total': total_gene_count, 'gene_statistics': gene_statistics})

    guides_by_exon = ranker.get_guides_by_exon()
    guide_count = ranker.get_count_selected_guides()
    result = {
      'gene_to_exon': guides_by_exon,
      'guide_count': guide_count,
      'genes': genes,
      'tissues': tissues,
      'quantity': quantity,
      'gtex_enabled': gtex_enabled,
      'tissues_disabled': tissues_disabled,
      'domains_enabled': domains_enabled,
      'email_address': email_address,
      'rejected_genes': rejected_genes,
      'genome': 'hum'
    }

    # print "Spent {0} seconds generating {1} guides/gene for {2} genes.".format(time.time() - t0, quantity, len(genes))

    # Send results to the user
    send_completed_run(email_address, start_compute.request.id)

    return {
      'current': total_gene_count,
      'total': total_gene_count,
      'status': 'Task completed!',
      'gene_statistics': gene_statistics,
      'result': result,
      'email_address': email_address,
      'genome': 'hum'
    }

@app.route('/generate', methods=['POST'])
def generate():
  # Retrieve arguments
  post = request.get_json()

  genes = post.get('genes')
  species = post.get('genome')
  quantity = post.get('quantity')
  domains_enabled = post.get('domains_enabled')
  email_address = post.get('email_address')
  rejected_genes = post.get('rejected_genes')
  gene_statistics = post.get('gene_statistics')

  print species
  if species == "mus":
    print "hi in mouse"
    print genes
    params = {
      'genes': genes,
      'species': species,
      'quantity': quantity,
      'domains_enabled': domains_enabled,
      'gene_statistics': gene_statistics,
      'email_address': email_address,
      'rejected_genes': rejected_genes
    }

    task = start_compute_mouse.apply_async(args=[params])
    return jsonify({'task_id': task.id})

  # Human setup
  else:
    tissues = post.get('tissues')
    gtex_enabled = post.get('gtex_enabled')
    tissues_disabled = post.get('tissues_disabled')

    # Validations
    if genes == None:
      genes = [{u'ensembl_id': u'ENSG00000186575.13', u'$$hashKey': u'object:50', u'name': u'NF2'}]
    if species == None:
      species = "human"
    if quantity == None:
      quantity = 6
    if tissues == None:
      tissues = ["Muscle", "Heart", "Brain"]

    params = {
      'genes': genes,
      'species': species,
      'quantity': quantity,
      'tissues': tissues,
      'gtex_enabled': gtex_enabled,
      'tissues_disabled': tissues_disabled,
      'domains_enabled': domains_enabled,
      'gene_statistics': gene_statistics,
      'email_address': email_address,
      'rejected_genes': rejected_genes
    }

    task = start_compute.apply_async(args=[params])
    return jsonify({'task_id': task.id})
    #return jsonify({}), 202, {'Location': url_for('taskstatus', task_id = task.id)}

@app.route('/status/<task_id>')
def taskstatus(task_id):
  task = start_compute.AsyncResult(task_id)
  try:
    if task.state == 'PENDING':
      response = {
        'state': task.state,
        'current': 0,
        'total': 1,
        'status': 'Pending...',
        'gene_statistics': task.info.get('gene_statistics', '')
      }
    elif task.state != 'FAILURE':
      response = {
        'state': task.state,
        'current': task.info.get('current', 0),
        'total': task.info.get('total', 1),
        'status': task.info.get('status', ''),
        'gene_statistics': task.info.get('gene_statistics', '')
      }
      if 'result' in task.info:
        response['result'] = task.info['result']
    else:
      # something went wrong in the background job
      response = {
        'state': task.state,
        'current': 1,
        'total': 1,
        'status': str(task.info), # exception raised
      }
    return jsonify(response)
  except AttributeError:
    return jsonify({
        'state': task.state,
        'current': 1,
        'total': 1,
        'status': "Unknown..."
      })
