FlaskStart.controller 'DesignerCtrl', ['$scope', '$filter', '$location', '$routeParams', 'GuidesFactory', 'Analytics', ($scope, $filter, $location, $routeParams, GuidesFactory, Analytics) ->
  $scope.guidesReady = false
  guidesFactory = new GuidesFactory()

  # Check for task_id
  if $routeParams.task_id?
    $scope.getGuidesPromise = guidesFactory.getComputedGuides($routeParams.task_id)

  else
    $scope.generateGuidesPromise = guidesFactory.generateGuides()

  # Track Analytics
  Analytics.trackEvent('designer', 'begin', 'genes', guidesFactory.data.genes.length, true, { genes: guidesFactory.data.genes })

  # Bar chart setup
  base_options = {
    animation: true,
    barDatasetSpacing : 1,
    responsive: true,
    maintainAspectRatio: false,
    scaleShowHorizontalLines: false,
    scaleIntegersOnly: true,
    scaleBeginAtZero: true,
    scaleShowGridLines : false,
    #legendTemplate : "<ul class=\"<%=name.toLowerCase()%>-legend\"><% for (var i=0; i<datasets.length; i++){%><li><span style=\"background-color:<%=datasets[i].fillColor%>\"></span><%if(datasets[i].label){%><%=datasets[i].label%><%}%></li><%}%></ul>"
    legendTemplate: "",
    scaleOverride: true,
    scaleSteps : 1,
    scaleStepWidth: 1,
    scaleStartValue: 0,
    barShowStroke: false
  }

  expression_colors = [
      {
        "fillColor": "#2B333B"
      },
      {
        "fillColor": "#EDD6F9"
      },
      {
        "fillColor": "#51D2B7"
      },
      {
        "fillColor": "#FED373"
      },
      {
        "fillColor": "#E7454E"
      },
      {
        "fillColor": "#a0a0a0"
      },
      {
        "fillColor": "#f2eaea"
      }
    ]

  expression_options = base_options
  guide_options = base_options
  guide_options.scaleOverride = false

  $scope.chart_config = {
    'expression': {
      'data':    [[]] # set later
      'labels':  [] # set later
      'series':  ['All Tissues', 'Selected Tissues', 'Brain', 'Heart', 'Kidney', 'Liver', 'Skin']
      'options': expression_options
      'colors': expression_colors
    },
    'guides': {
      'data':    [[]] # set later
      'labels':  [] # set later
      'series':  ['Guides per Exon']
      'options': guide_options
      'colors': [
        {
          "fillColor": "#2B333B"
        },
        {
          "fillColor": "#51D2B7" # light gray from exon background
        },
        {
          "fillColor": "rgba(224, 108, 112, 1)",
          "strokeColor": "rgba(207,100,103,1)",
          "pointColor": "rgba(220,220,220,1)",
          "pointStrokeColor": "#fff",
          "pointHighlightFill": "#fff",
          "pointHighlightStroke": "rgba(151,187,205,0.8)"
        }
      ]
    }
  }

  console.log $scope.chart_config.guides.options
  console.log $scope.chart_config.expression.options

  # intitalize the svg_unit. It will be modified later by the drawIndividualExon directive.
  $scope.modifySvgUnit = (unit) ->
    $scope.svg_unit_global = unit / 2
    $scope.chart_config.expression.options.barValueSpacing = $scope.svg_unit_global
    $scope.chart_config.guides.options.barValueSpacing = $scope.svg_unit_global

    if angular.element('#exon_graph').length > 0
      chartwidth = angular.element('#exon_graph')[0].getBoundingClientRect().width

      angular.element('.chart-container').css('width', chartwidth - unit - 22)
      angular.element('.chart-container').css('margin-left', unit / 2 + 10)
      return

  # Initialize
  $scope.modifySvgUnit(15)

  # For highlighting exons on overlay
  $scope.exonHovered = -1

  $scope.setExonHovered = (val) ->
    $scope.exonHovered = val

  # restructure the received json
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

  if $routeParams.task_id?
    $scope.getGuidesPromise.then (guidesData) ->
      computeGuidesData(guidesData["gene_to_exon"])
      $scope.gene = $scope.gene_to_exon[0]
      $scope.guidesReady = true

      $scope.guidesFactoryData = guidesFactory.data
      $scope.tissues = guidesFactory.data.tissues
      $scope.tissues_enabled = not guidesFactory.data.tissues_disabled

      # Change series if we are not going to display median
      if guidesFactory.data.tissues_disabled
        $scope.chart_config.expression.series = ['All Tissues', 'Brain', 'Heart', 'Kidney', 'Liver', 'Skin']
        $scope.chart_config.expression.colors.splice(1,1)

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
  else
    $scope.generateGuidesPromise.then (task_id) ->
      console.log 'generateguidepromise finished with task_id = ' + task_id
      $location.path('/designer/' + task_id)

  $scope.setGene = (idx) ->
    $scope.gene = $scope.gene_to_exon[idx]
    expression_labels = []
    guide_labels = []
    expression_data1 = []
    expression_data2 = []
    expression_data = []
    guides_data = []

    # Setup the correct number of arrays to store our data later.
    # Assumption: Every exon has same number of expression values.
    # ...but this is required anyway by our chart library.
    angular.forEach $scope.gene.exons[0].expression, (tissue_value, key) ->
      expression_data.push []
    # find normalizing constant
    # Normalize by max(max(expression_overalls), max(expression_median))
    max_expression = 0
    angular.forEach $scope.gene.exons, (exon, key) ->
      angular.forEach exon.expression, (tissue_value, key) ->
        if tissue_value > max_expression
          max_expression = tissue_value
      guides_count = ($filter('filter')(exon.gRNAs, {selected:true}, true)).length
      guides_data.push guides_count
    angular.forEach $scope.gene.exons, (exon, key) ->
      if guidesFactory.data.tissues_disabled # don't use median
        expression_data[0].push (exon.expression.overall / max_expression).toFixed(2)
        expression_data[1].push (exon.expression.brain / max_expression).toFixed(2)
        expression_data[2].push (exon.expression.heart / max_expression).toFixed(2)
        expression_data[3].push (exon.expression.kidney / max_expression).toFixed(2)
        expression_data[4].push (exon.expression.liver / max_expression).toFixed(2)
        expression_data[5].push (exon.expression.skin / max_expression).toFixed(2)
      else
        expression_data[0].push (exon.expression.overall / max_expression).toFixed(2)
        expression_data[1].push (exon.expression.median / max_expression).toFixed(2)
        expression_data[2].push (exon.expression.brain / max_expression).toFixed(2)
        expression_data[3].push (exon.expression.heart / max_expression).toFixed(2)
        expression_data[4].push (exon.expression.kidney / max_expression).toFixed(2)
        expression_data[5].push (exon.expression.liver / max_expression).toFixed(2)
        expression_data[6].push (exon.expression.skin / max_expression).toFixed(2)
      expression_labels.push('Exon ' + (key+1))
      guide_labels.push('') # empty labels

    $scope.chart_config.expression.labels = expression_labels
    $scope.chart_config.guides.labels = guide_labels
    $scope.chart_config.guides.data = [guides_data]
    $scope.chart_config.expression.data = expression_data

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

  $scope.removeRejectedGene = (idx) ->
    guidesFactory.data.rejected_genes.splice(idx, 1)

  $scope.removeTissue = (idx) ->
    $scope.guidesReady = false
    guidesFactory.data.tissues.splice(idx, 1)
    $scope.generateGuidesPromise = guidesFactory.generateGuides()
    $scope.generateGuidesPromise.then (task_id) ->
      $location.path('/designer/' + task_id)

  # Searching
  $scope.geneTissueQuery = ""
  $scope.geneTissueSearch = () ->
    for elt in $scope.geneTissueQuery.split(',')
      elt = elt.replace(/ /g,'')
      found = false
      for tissue in guidesFactory.available.tissues
        if tissue.toUpperCase() == elt.toUpperCase()
          guidesFactory.data.tissues.push(tissue)
          guidesFactory.data.tissues_disabled = false
          containsTissue = true
          break
      if found == false
        for gene in guidesFactory.available.genes
          if gene.name.toUpperCase() == elt.toUpperCase() or gene.ensembl_id.toUpperCase() == elt.toUpperCase()
            guidesFactory.data.genes.push(gene)
            found = true
            break

    $scope.guidesReady = false
    $scope.generateGuidesPromise = guidesFactory.generateGuides().then (task_id) ->
      $location.path('/designer/' + task_id)

  # individual guide selection
  $scope.show_different_guides = false

  # UI-Select setup
  $scope.guidesFactory = guidesFactory # only used by ui-select!
  $scope.selectedGene = (gene) ->
    guidesFactory.data.genes.push gene
    $scope.guidesReady = false
    $scope.generateGuidesPromise = guidesFactory.generateGuides().then (task_id) ->
      $location.path('/designer/' + task_id)

  $scope.selectedTissue = (tissue) ->
    guidesFactory.data.tissues.push tissue
    $scope.guidesReady = false
    $scope.generateGuidesPromise = guidesFactory.generateGuides().then (task_id) ->
      $location.path('/designer/' + task_id)

  # $scope.additionalGene = "" # what the person is typing.
  # $scope.$watch 'additionalGene', () ->
  #   console.log "hi"
  #   console.log $scope.additionalGene
  #   guidesFactory.data.genes.push $scope.additionalGene
  #   $scope.guidesReady = false
  #   $scope.generateGuidesPromise = guidesFactory.generateGuides().then (guidesData) ->
  #     computeGuidesData(guidesData["gene_to_exon"])
  #     $scope.guidesReady = true
  # , true


  $scope.guideSelected = (guide) ->
    exon_key = guide.exon - 1 # dynamically update chart
    if guide.selected == false
      $scope.countSelectedGuides -= 1
      $scope.chart_config.guides.data[0][exon_key] -= 1
    else
      $scope.countSelectedGuides += 1
      $scope.chart_config.guides.data[0][exon_key] += 1

  # not in usage right now
  # colorBar = (exon, color) ->
  #   $scope.chart_config.guides.colors[0]["fillColor"][exon] = color

  pad = (n, width, z) ->
    z = z or '0'
    n = n + ''
    if n.length >= width then n else new Array(width - n.length + 1).join(z) + n

  $scope.getGuidesCSV = ->
    guidesCSV = $filter('filter')($scope.merged_gRNAs, {selected:true}, true)
    guidesCSV = $filter('orderBy')(guidesCSV, ['gene','score'], true)
    angular.forEach guidesCSV, (guide, idx) ->
      guide.uid = "customLibrary_guide" + pad(idx, 4)
    guidesCSV

]
