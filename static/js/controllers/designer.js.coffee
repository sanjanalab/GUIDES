FlaskStart.controller 'DesignerCtrl', ['$scope', '$filter', 'GuidesFactory', ($scope, $filter, GuidesFactory) ->
  guidesFactory = new GuidesFactory()

  guidesFactory.generateGuides("AAk1").then (guidesData) ->
    console.log guidesData
    $scope.guidesData = guidesData["data"]["gene_to_exon"]
    gene_to_exon = guidesData["data"]["gene_to_exon"]
    guide_count = guidesData["data"]["guide_count"]
    all_gRNAs = {}
    merged_gRNAs = []

    # Display attributes for data
    # p_ values are pixel values (as opposed ot sequencing data)
    pixel_width = 800
    $scope.countSelectedGuides = guide_count
    angular.forEach gene_to_exon, (gene, key1) ->
      all_gRNAs[gene.name] = []
      angular.forEach gene.exons, (exon, key2) ->
        exon.p_start =  exon.start / gene.length * pixel_width
        exon.p_end = exon.end / gene.length * pixel_width
        angular.forEach exon.gRNAs, (guide, key3) ->
          # guide.selected = false # Change later to only include best guides -> might even come from server
          # if guide.selected
          #   $scope.countSelectedGuides += 1
          guide.p_start = guide.start / gene.length * pixel_width
          guide.exon = key2 + 1
          guide.gene = gene.name
          all_gRNAs[gene.name].push(guide)
          merged_gRNAs.push(guide)

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

    $scope.gene_to_exon = gene_to_exon
    $scope.all_gRNAs = all_gRNAs
    
    $scope.gene = gene_to_exon[0]
    $scope.setGene = (idx) ->
      $scope.gene = gene_to_exon[idx]
      console.log $scope.gene

    console.log $scope.gene

    $scope.guideSelected = (guide) ->
      console.log $scope.svgUnit
      if guide.selected == false
        $scope.countSelectedGuides -= 1
      else
        $scope.countSelectedGuides += 1

    $scope.getGuidesCSV = ->
      guidesCSV = $filter('filter')(merged_gRNAs, {selected:true}, true)
      guidesCSV = $filter('orderBy')(guidesCSV, 'score', true)
      guidesCSV

    # console.log "in filter"
    # console.log gene_length
    # console.log exon
    # pixel_width = 900 # length of displayed gene
    # {
    #   "start": exon.start / gene_length * pixel_width,
    #   "end": exon.end / gene_length * pixel_width
    # }

# exon.start = exon.start / gene_length * pixel_width

# exon.end = exon.end / gene_length * pixel_width

# exon

  # $scope.normalizeExons = (gene_length) -> (exon) ->
  # 	console.log exon
  # 	pixel_width = 800 #length of displayed gene
  # 	exon.start = exon.start / gene_length * pixel_width
  # 	exon.end = exon.end / gene_length * pixel_width
  # 	exon

  # getGuides = () ->
  # 	guidesFactory.generateGuides("AAk1").success((guidesData) ->
  # 		guidesData
  # 	)

  # console.log getGuides()
  # $scope.getGuides = (gene) ->
  # 	console.log gene
  # 	#console.log guidesFactory.generateGuides(gene)
  # 	# guidesFactory.generateGuides(gene).then((guidesData) ->
  # 	# 	console.log guidesData
  # 	# 	guidesData
  # 	# )
]