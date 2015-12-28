FlaskStart.controller 'DesignerCtrl', ['$scope', '$filter', 'GuidesFactory', ($scope, $filter, GuidesFactory) ->
  $scope.data = {}
  $scope.posts = ["hi", "hi2", "hi3"]
  #$scope.guidesData = ['aa1']
  guidesFactory = new GuidesFactory()

  guidesFactory.generateGuides("AAk1").then (guidesData) ->
    console.log guidesData
    $scope.guidesData = guidesData["data"]["gene_to_exon"]
    gene_to_exon = guidesData["data"]["gene_to_exon"]
    all_gRNAs = {}
    merged_gRNAs = []


    # Display attributes for data
    # p_ values are pixel values (as opposed ot sequencing data)
    pixel_width = 800
    angular.forEach gene_to_exon, (gene, key1) ->
      all_gRNAs[gene.name] = []
      angular.forEach gene.exons, (exon, key2) ->
        exon.p_start =  exon.start / gene.length * pixel_width
        exon.p_end = exon.end / gene.length * pixel_width
        angular.forEach exon.gRNAs, (guide, key3) ->
          guide.selected = false # Change later to only include best guides -> might even come from server
          guide.p_start = guide.start / gene.length * pixel_width
          guide.exon = key2
          guide.gene = gene.name
          all_gRNAs[gene.name].push(guide)
          merged_gRNAs.push(guide)

    ## I think this is unnecessary, since we filter by order in the template.
    # angular.forEach all_gRNAs, (guides_for_gene, gene_name) ->
    #   all_gRNAs[gene_name] = $filter('orderBy')(guides_for_gene, 'score', true)

    guide_count = guidesData["data"]["guide_count"]
    $scope.countSelectedGuides = guide_count

    merged_gRNAs = $filter('orderBy')(merged_gRNAs, 'score', true)
    angular.forEach merged_gRNAs, (guide, key) ->
      if guide_count > 0
        guide.selected = true
        guide_count -= 1
      else
        guide.selected = false

    $scope.gene_to_exon = gene_to_exon
    $scope.all_gRNAs = all_gRNAs
    

    $scope.guideSelected = (guide) ->
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