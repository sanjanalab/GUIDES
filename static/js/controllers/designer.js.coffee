FlaskStart.controller 'DesignerCtrl', ['$scope', '$filter', 'GuidesFactory', ($scope, $filter, GuidesFactory) ->
  computeGuidesData = (gene_to_exon) ->
    # $scope.guidesData = guidesData["gene_to_exon"]
    # gene_to_exon = guidesData["gene_to_exon"]
    #guide_count = guidesData["guide_count"]
    $scope.gene_to_exon = gene_to_exon

    all_gRNAs = {}
    merged_gRNAs = []

    # Display attributes for data
    # p_ values are pixel values (as opposed ot sequencing data)
    # was using p_ approach before switching to a directive.
    pixel_width = 800
    countSelectedGuides = 0 #guide_count
    angular.forEach gene_to_exon, (gene, key1) ->
      all_gRNAs[gene.name] = []
      angular.forEach gene.exons, (exon, key2) ->
        exon.p_start = exon.start / gene.length * pixel_width
        exon.p_end = exon.end / gene.length * pixel_width
        angular.forEach exon.gRNAs, (guide, key3) ->
          # guide.selected = false # Change later to only include best guides -> might even come from server
          if guide.selected
            countSelectedGuides += 1
          guide.p_start = guide.start / gene.length * pixel_width
          guide.exon = key2 + 1
          guide.gene = gene.name
          all_gRNAs[gene.name].push(guide)
          merged_gRNAs.push(guide)
    $scope.countSelectedGuides = countSelectedGuides
    $scope.all_gRNAs = all_gRNAs
    $scope.merged_gRNAs = merged_gRNAs

  guidesFactory = new GuidesFactory()
  $scope.gene_statistics = guidesFactory.gene_statistics
  $scope.generateGuidesPromise = guidesFactory.generateGuides()

  $scope.generateGuidesPromise.then (guidesData) ->
    computeGuidesData(guidesData["gene_to_exon"])

    ## I think this is unnecessary, since we filter by order in the template.
    # angular.forEach all_gRNAs, (guides_for_gene, gene_name) ->
    #   all_gRNAs[gene_name] = $filter('orderBy')(guides_for_gene, 'score', true)

    # Server is now doing this, so this has been removed.
    # guide_count = $scope.countSelectedGuides
    # merged_gRNAs = $filter('orderBy')(merged_gRNAs, 'score', true)
    # angular.forEach merged_gRNAs, (guide, key) ->
    #   if guide_count > 0
    #     guide.selected = true
    #     guide_count -= 1
    #   else
    #     guide.selected = false

    # used for table column sorting
    $scope.orderByField = 'score'
    $scope.reverseSort = true

    #$scope.gene_to_exon = gene_to_exon
    #$scope.all_gRNAs = all_gRNAs
    
    $scope.gene = $scope.gene_to_exon[0]
    $scope.setGene = (idx) ->
      $scope.gene = $scope.gene_to_exon[idx]

    $scope.removeGene = (idx) ->
      $scope.gene_to_exon.splice(idx, 1)
      computeGuidesData($scope.gene_to_exon)

    # Searching
    $scope.geneTissueQuery = ""
    $scope.geneTissueSearch = () ->
      console.log "in search"
      for elt in $scope.geneTissueQuery.split(',')
        elt = elt.replace(/ /g,'')
        found = false
        console.log elt
        console.log guidesFactory.available.tissues
        for tissue in guidesFactory.available.tissues
          if tissue.toUpperCase() == elt.toUpperCase()
            console.log elt + "is tissue"
            guidesFactory.data.tissues.push(tissue)
            found = true
            break
        if found == false
          for gene in guidesFactory.available.genes
            if gene.name.toUpperCase() == elt.toUpperCase() or gene.ensembl_id.toUpperCase() == elt.toUpperCase()
              console.log elt + "is gene"
              guidesFactory.data.genes.push(gene)
              found = true
              break
      console.log guidesFactory
      console.log guidesFactory.data
      $scope.generateGuidesPromise = guidesFactory.generateGuides().then (guidesData) ->
        computeGuidesData(guidesData["gene_to_exon"])


    $scope.guideSelected = (guide) ->
      if guide.selected == false
        $scope.countSelectedGuides -= 1
      else
        $scope.countSelectedGuides += 1

    $scope.getGuidesCSV = ->
      guidesCSV = $filter('filter')($scope.merged_gRNAs, {selected:true}, true)
      guidesCSV = $filter('orderBy')(guidesCSV, 'score', true)
      guidesCSV
]