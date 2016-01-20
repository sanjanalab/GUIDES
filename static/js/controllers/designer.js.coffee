FlaskStart.controller 'DesignerCtrl', ['$scope', '$filter', 'GuidesFactory', ($scope, $filter, GuidesFactory) ->
  guidesFactory = new GuidesFactory()
  $scope.generateGuidesPromise = guidesFactory.generateGuides()
  $scope.gene_statistics = guidesFactory.gene_statistics

  # intitalize the svg_unit. It will be modified later by the drawIndividualExon directive. 
  $scope.modifySvgUnit = (unit) ->
    $scope.svg_unit_global = unit / 2
    $scope.exon_chart_options = {
      animation: false,
      barValueSpacing : $scope.svg_unit_global,
      responsive: true,
      maintainAspectRatio: false,
      scaleShowHorizontalLines: false,
      scaleIntegersOnly: true,
      scaleBeginAtZero: true,
      scaleShowGridLines : false,
      legendTemplate : "<ul class=\"<%=name.toLowerCase()%>-legend\"><% for (var i=0; i<datasets.length; i++){%><li><span style=\"background-color:<%=datasets[i].fillColor%>\"></span><%if(datasets[i].label){%><%=datasets[i].label%><%}%></li><%}%></ul>"
      scaleOverride: true,
      scaleSteps : 1,
      scaleStepWidth: 1,
      scaleStartValue: 0,
      barShowStroke: false
    }

  # Initialize
  $scope.modifySvgUnit(15)

  $scope.series = ['Mean All Tissues', 'Mean Selected Tissues']

  $scope.chartcolors = [{
    "fillColor": "#2B333B" # dark color from sidebar
  },
  {
    "fillColor": "#51D2B7" # light gray from exon background
  }
  {
    "fillColor": "rgba(224, 108, 112, 1)",
    "strokeColor": "rgba(207,100,103,1)",
    "pointColor": "rgba(220,220,220,1)",
    "pointStrokeColor": "#fff",
    "pointHighlightFill": "#fff",
    "pointHighlightStroke": "rgba(151,187,205,0.8)"
  }]

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

    # simulate setting up the first gene
    $scope.setGene(0)


  $scope.generateGuidesPromise.then (guidesData) ->
    computeGuidesData(guidesData["gene_to_exon"])
    $scope.gene = $scope.gene_to_exon[0]
    $scope.guidesReady = true

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
  
  $scope.setGene = (idx) ->
    $scope.gene = $scope.gene_to_exon[idx]
    labels = []
    data1 = []
    data2 = []
    # find normalizing constant
    # Normalize by max(max(expression_overalls), max(expression_median))
    max_expression = 0
    angular.forEach $scope.gene.exons, (exon, key) ->
      if exon.expression_overall > max_expression
        max_expression = exon.expression_overall
      if exon.expression_median > max_expression
        max_expression = exon.expression_median
    angular.forEach $scope.gene.exons, (exon, key) ->
      data1.push exon.expression_overall / max_expression
      data2.push exon.expression_median / max_expression
      labels.push('Exon ' + (key+1))
    $scope.labels = labels
    if guidesFactory.data.tissues_enabled
      $scope.data = [data1,data2]
    else
      $scope.data = [data1]

  # returns the actual exons
  $scope.exonsUtilized = (gene) ->
    exons = []
    if not gene or not gene.hasOwnProperty('exons')
      return 0
    for exon in gene.exons
      for guide in exon.gRNAs
        if guide.selected
          exons.push(exon)
          break
    exons

  $scope.selectedGuides = (gene_name) ->
    guides = []
    for guide in $scope.all_gRNAs[gene_name]
      if guide.selected
        guides.push(guide)
    guides

  $scope.removeGene = (idx) ->
    guidesFactory.data.genes.splice(idx, 1)
    $scope.gene_to_exon.splice(idx, 1)
    computeGuidesData($scope.gene_to_exon)

  # Searching
  $scope.geneTissueQuery = ""
  $scope.geneTissueSearch = () ->
    for elt in $scope.geneTissueQuery.split(',')
      elt = elt.replace(/ /g,'')
      found = false
      for tissue in guidesFactory.available.tissues
        if tissue.toUpperCase() == elt.toUpperCase()
          guidesFactory.data.tissues.push(tissue)
          guidesFactory.data.tissues_enabled = true
          containsTissue = true
          break
      if found == false
        for gene in guidesFactory.available.genes
          if gene.name.toUpperCase() == elt.toUpperCase() or gene.ensembl_id.toUpperCase() == elt.toUpperCase()
            guidesFactory.data.genes.push(gene)
            found = true
            break

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