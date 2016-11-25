FlaskStart.factory 'GuidesFactory', ['$http', '$q', '$filter', '$timeout', ($http, $q, $filter, $timeout) ->
  class GuidesFactory
    available: # can select
      'genes': [] # set later by $http in constructor
      'human_genes': []
      'mouse_genes': []
      'tissues': ['Thyroid', 'Testis', 'Cervix Uteri', 'Adipose Tissue', 'Breast', 'Vagina', 'Nerve', 'Pituitary', 'Stomach', 'Fallopian Tube', 'Bone Marrow', 'Bladder', 'Blood', 'Colon', 'Prostate', 'Pancreas', 'Blood Vessel', 'Liver', 'Spleen', 'Small Intestine', 'Uterus', 'Ovary', 'Muscle', 'Heart', 'Adrenal Gland', 'Brain', 'Salivary Gland', 'Lung', 'Skin', 'Esophagus', 'Kidney']

    data: # currently selected
      'loading': false
      'genes' : []
      'rejected_genes': []
      'genome': 'hum'
      'tissues' : []
      'quantity': 60
      'gtex_enabled': true
      'tissues_disabled': true # does the user want to consider *individualized* tissue expression?
      'domains_enabled': true
      'genesFromFile': []
      'emailAddress': ''
      'gene_statistics':
        'expected': 0
        'actual': 0
        'processed': 0
    # expected --> how many we would get if user's file is 100% correct.
    # actual   --> how much we could actually derived from the user's file.

    # Pre-computed results from the server
    waitTime: () ->
      idx = this.hasUploadedFile()
      gene_count = this.data.genes.length
      if idx != -1
        gene_count += this.data.genesFromFile.length
      return 0.13 * gene_count

    waitTimeText: () ->
      wait = this.waitTime()
      if wait > 60
        "#{Math.ceil(wait / 60)} mins"
      else
        "#{Math.ceil(wait)} secs"

    hasUploadedFile: () ->
      ret = -1
      angular.forEach this.data.genes, (gene, key) ->
        if gene.ensembl_id == "GENES_FROM_FILE"
          ret = key
      return ret

    noteFileUploaded: (name) ->
      idx = this.hasUploadedFile
      if idx != -1
        this.data.genes.splice(idx, 1)

      this.data.genes.push({
        'name': name
        'ensembl_id': "GENES_FROM_FILE"
      })

    constructor: () ->
      # Setup default available, and default selected
      # We are in /static/js/min/scripts.min.js
      # We want /static/data/pre_processed/genes_list.json
      if this.available.genes.length == 0 # implies we have not yet preloaded
        this_ = this
        $http.get('/static/data/pre_processed/genes_list.json').then (res) ->
          this_.available.human_genes = res.data
          this_.available.genes = this_.available.human_genes # only once on constructor
          this_.data.genes = []
        $http.get('/static/data/pre_processed/genes_list_GRCm38.txt').then (res) ->
          this_.available.mouse_genes = res.data
          this_.data.genes = []

    prepareGenesFromFile: () ->
      this_ = this
      $q (resolve, reject) ->
        expected_genes_num = this_.data.genes.length
        idx = this_.hasUploadedFile()
        if (idx != -1) and this_.data.genesFromFile
          expected_genes_num = expected_genes_num + this_.data.genesFromFile.length - 1 # subtract "GENES_FROM_FILE holder"
          this_.data.genes.splice(idx, 1) # remove the "GENES_FROM_FILE holder"
          seen = {}
          for geneText in this_.data.genesFromFile
            foundMatch = false
            for gene in this_.available.genes
              if ((gene.name.toUpperCase() == geneText.toUpperCase() or gene.ensembl_id == geneText or gene.ensembl_id.substring(0, geneText.length) == geneText) and not (gene.ensembl_id in seen))
                seen[gene.ensembl_id] = true
                this_.data.genes.push(gene)
                foundMatch = true
                break
            if not foundMatch
              this_.data.rejected_genes.push(geneText)

        this_.data.gene_statistics.expected = expected_genes_num
        this_.data.gene_statistics.actual = this_.data.genes.length
        resolve()

    generateGuides: () ->
      deferred = $q.defer()

      promise = this.prepareGenesFromFile()
      this_ = this
      promise.then () ->
        if this_.data.tissues.length == 0 or this_.data.tissues_disabled
          this_.data.tissues = this_.available.tissues
        $http {
          url: '/generate'
          method: 'POST'
          headers:
            'Content-Type': 'application/json'
          data: JSON.stringify(this_.data)
        }
        .success (data) ->
          task_id = data.task_id
          deferred.resolve task_id

      deferred.promise

    updateProgress: (task_id) ->
      deferred = $q.defer()

      this.data.loading = true

      this_ = this
      $http {
        url: '/status/' + task_id
        method: 'GET'
        headers:
          'Content-Type': 'application/json'
      }
      .success (data) ->
        this_.data.loading = false
        # check if we are finished
        if data.state == 'SUCCESS'
          deferred.resolve data
        else
          this_.data.gene_statistics = data.gene_statistics
          return $timeout (() ->
            this_.updateProgress(task_id).then (data) ->
              deferred.resolve data
            ), 2000, true, task_id

      deferred.promise

    getComputedGuides: (task_id) ->
      deferred = $q.defer()
      this_ = this

      this.updateProgress(task_id).then (data) ->
        this_.data.genome = data.result.genome
        this_.data.genes = data.result.genes
        this_.data.tissues = data.result.tissues
        this_.data.quantity = data.result.quantity
        this_.data.domains_enabled = data.result.domains_enabled
        this_.data.rejected_genes = data.result.rejected_genes
        this_.data.gene_statistics = data.gene_statistics
        this_.data.gtex_enabled = data.result.gtex_enabled
        this_.data.tissues_disabled = data.result.tissues_disabled

        # return the data itself
        deferred.resolve data.result

      deferred.promise

    setGenes: (genes) ->
      this.data.genes = genes

    setTissues: (tissues) ->
      this.data.tissues = tissues
]
