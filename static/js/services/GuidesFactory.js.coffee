FlaskStart.factory 'GuidesFactory', ['$http', '$q', '$filter', ($http, $q, $filter) ->
  class GuidesFactory
    available: # can select
      'genes': [] # set later by $http in constructor
      'tissues': ['Thyroid', 'Testis', 'Cervix Uteri', 'Adipose Tissue', 'Breast', 'Vagina', 'Nerve', 'Pituitary', 'Stomach', 'Fallopian Tube', 'Bone Marrow', 'Bladder', 'Blood', 'Colon', 'Prostate', 'Pancreas', 'Blood Vessel', 'Liver', 'Spleen', 'Small Intestine', 'Uterus', 'Ovary', 'Muscle', 'Heart', 'Adrenal Gland', 'Brain', 'Salivary Gland', 'Lung', 'Skin', 'Esophagus', 'Kidney']

    data: # currently selected
      'genes' : []
      'tissues' : []
      'quantity': 60
      'gtex_enabled': true
      'tissues_disabled': true # does the user want to consider *individualized* tissue expression?

    # paramters for getting genes from file
    genesFromFile: []

    # expected --> how many we would get if user's file is 100% correct.
    # actual   --> how much we could actually derived from the user's file.
    gene_statistics:
      expected: 0
      actual: 0

    hasUploadedFile: () ->
      ret = false
      angular.forEach this.data.genes, (gene, key) ->
        if gene.ensembl_id == "GENES_FROM_FILE"
          ret = key
      return ret

    noteFileUploaded: () ->
      this.data.genes.push({
        'name': "Genes from file..."
        'ensembl_id': "GENES_FROM_FILE"
      })

    constructor: () ->
      # Setup default available, and default selected
      # We are in /static/js/min/scripts.min.js
      # We want /static/data/pre_processed/genes_list.json
      if this.available.genes.length == 0 # implies we have not yet preloaded
        this_ = this
        $http.get('/static/data/pre_processed/genes_list.json').then (res) ->
          this_.available.genes = res.data
          this_.data.genes = []

    prepareGenesFromFile: () ->
      this_ = this
      $q (resolve, reject) ->
        expected_genes_num = this_.data.genes.length
        idx = this_.hasUploadedFile()
        if idx and this_.genesFromFile
          expected_genes_num = expected_genes_num + this_.genesFromFile.length - 1 # subtract "GENES_FROM_FILE holder"
          this_.data.genes.splice(idx, 1) # remove the "GENES_FROM_FILE holder"
          for geneText in this_.genesFromFile
            for gene in this_.available.genes
              if gene.name == geneText or gene.ensembl_id == geneText
                this_.data.genes.push(gene)
                break
            # replace (for gene in .... break) with below, if you don't need EXACT comparisons.
            # matchedGenes = $filter('propsFilter')(this.available.genes, {name: geneText, ensembl_id: geneText})
            # matchedGenes.sort()
            # if matchedGenes.length > 0
            #   this.data.genes.push(matchedGenes[0])
  
        this_.gene_statistics.expected = expected_genes_num
        this_.gene_statistics.actual = this_.data.genes.length
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
          deferred.resolve data

      deferred.promise

    setGenes: (genes) ->
      this.data.genes = genes

    setTissues: (tissues) ->
      this.data.tissues = tissues
]

# scope.postTest = function(){

#   var data = [obj1, obj2, obj3];
#   var jsonData=angular.toJson(data);
#   var objectToSerialize={'object':jsonData};

#   $http({
#     url: 'myURL',
#     method: "POST",
#     data: $.param(objectToSerialize),
#     headers: {
#          'Content-Type': 'application/x-www-form-urlencoded; charset=UTF-8'
#     }
#   }).success(function(data){
#     alert("done");
#   });
# }