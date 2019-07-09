FlaskStart.controller 'DesignerCtrl', ['$scope', '$filter', '$location', '$window', '$routeParams', '$http', 'GuidesFactory', 'Analytics', ($scope, $filter, $location, $window, $routeParams, $http, GuidesFactory, Analytics) ->
  $scope.guidesReady = false
  guidesFactory = new GuidesFactory()

  # expose metadata to view
  $scope.guidesFactoryData = guidesFactory.data

  $scope.permalink = $window.location.href

  # non-targeting guides info
  $scope.non_targeting_guides_info = {
    'count': 500
    'default': 500
    'min': 0
    'max': 1000
    'text': "Default: Make non-targeting sgRNAs 5% of final library size"
  }

  # synthesis ready compounds
  $scope.synthesis_ready = {
    'enabled': false
    'type': 'standard' # or EF
    'csv_header_normal': ['GUIDE ID', 'Gene', 'Ensembl ID', 'Sequence', 'Guide with flanking sequence (10bp before and after the guide + PAM)','Chromosome', 'cut position in chromosome', 'On-target efficiency', 'Exon', 'Targets last exon', '10bp off-target match', 'Off-target score', 'Protein domain']
    'csv_column_order_normal': ['uid', 'gene', 'ensembl_gene', 'seq', 'full_seq', 'chrom', 'cut_pos', 'score', 'exon', 'targets_last_exon', 'has_exome_repeat', 'off_target_score', 'functional_domain']
    'csv_header_scaffold_standard': ['GUIDE ID', 'Gene', 'Ensembl ID', 'Sequence', 'Guide with flanking sequence (10bp before and after the guide + PAM)','Full-length sgRNA scaffold oligo for synthesis','Chromosome', 'cut position in chromosome', 'On-target efficiency', 'Exon', 'Targets last exon', '10bp off-target match', 'Off-target score', 'Protein domain']
    'csv_column_order_scaffold_standard':['uid', 'gene', 'ensembl_gene', 'seq', 'full_seq', 'scaffold_standard', 'chrom', 'cut_pos', 'score', 'exon', 'targets_last_exon', 'has_exome_repeat', 'off_target_score', 'functional_domain']
    'csv_header_scaffold_EF': ['GUIDE ID', 'Gene', 'Ensembl ID', 'Sequence', 'Guide with flanking sequence (10bp before and after the guide + PAM)','E+F modified sgRNA scaffold oligo for synthesis','Chromosome', 'cut position in chromosome', 'On-target efficiency', 'Exon', 'Targets last exon', '10bp off-target match', 'Off-target score', 'Protein domain']
    'csv_column_order_scaffold_EF':['uid', 'gene', 'ensembl_gene', 'seq', 'full_seq', 'scaffold_EF', 'chrom', 'cut_pos', 'score', 'exon', 'targets_last_exon', 'has_exome_repeat', 'off_target_score', 'functional_domain']
  }
  $scope.synthesis_ready.csv_header = $scope.synthesis_ready.csv_header_normal
  $scope.synthesis_ready.csv_column_order = $scope.synthesis_ready.csv_column_order_normal


  # going backwards
  $scope.$on '$routeChangeStart', (scope, next, current) ->
    if $scope.guidesReady and next.$$route.controller == 'DesignerCtrl'
      event.preventDefault()
      $location.path('index.html')

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
        "fillColor": "#EDD6F9"
      },
      {
        "fillColor": "#2B333B"
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
          if key2 == gene.exons.length - 1
            guide.targets_last_exon = true
          else
            guide.targets_last_exon = false
          guide.p_start = guide.start / gene.length * pixel_width
          guide.exon = key2 + 1
          guide.gene = gene.name
          guide.ensembl_gene = gene.ensembl_gene
          guide.full_seq = guide.seq_before + guide.seq + guide.PAM + guide.seq_after
          guide.scaffold_standard = "GGAAAGGACGAAACACCG" + guide.seq + "GTTTTAGAGCTAGAAATAGCAAGTTAAAATAAGGC"
          guide.scaffold_EF = "GGAAAGGACGAAACACCG" + guide.seq + "GTTTAAGAGCTATGCTGGAAACAGC"
          all_gRNAs[gene.name].push(guide)
          merged_gRNAs.push(guide)
    $scope.countSelectedGuides = countSelectedGuides
    $scope.all_gRNAs = all_gRNAs
    $scope.merged_gRNAs = merged_gRNAs

    # non-targeting guides count
    $scope.non_targeting_guides_info.count = Math.min(Math.ceil($scope.countSelectedGuides / 20), $scope.non_targeting_guides_info.max)
    if $scope.non_targeting_guides_info.count < 10
      $scope.non_targeting_guides_info.count = 10
      $scope.non_targeting_guides_info.default = 10
      $scope.non_targeting_guides_info.text = "Default: Make non-targeting sgRNAs 5% of final library size (min 10 sgRNAs)"
    else
      $scope.non_targeting_guides_info.default = $scope.non_targeting_guides_info.count

    # simulate setting up the first gene
    $scope.setGene(0)

  if $routeParams.task_id?
    $scope.getGuidesPromise.then (guidesData) ->
      $scope.guidesFactoryData = guidesFactory.data
      $scope.tissues = guidesFactory.data.tissues
      $scope.tissues_enabled = not guidesFactory.data.tissues_disabled

      # Change series if we are not going to display median
      if guidesFactory.data.tissues_disabled
        $scope.chart_config.expression.series = ['All Tissues', 'Brain', 'Heart', 'Kidney', 'Liver']
        $scope.chart_config.expression.colors.splice(1,1)

      else
        $scope.chart_config.expression.series = []
        if guidesFactory.data.tissues.length > 1
          $scope.chart_config.expression.series = ['All Tissues', 'Selected Tissues']

        tissue_index = 0
        while $scope.chart_config.expression.series.length < 6 and tissue_index < guidesFactory.data.tissues.length
          tissue_name = guidesFactory.data.tissues[tissue_index]
          $scope.chart_config.expression.series.push tissue_name
          tissue_index += 1

      $scope.chart_config.expression.colors = $scope.chart_config.expression.colors[...$scope.chart_config.expression.series.length]

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
  else
    $scope.generateGuidesPromise.then (task_id) ->
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
    angular.forEach $scope.chart_config.expression.series, (tissue_value, key) ->
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
      if exon.expression
        cur_index = 0

        if guidesFactory.data.tissues.length > 1
          expression_data[cur_index].push (exon.expression.overall / max_expression).toFixed(2)
          cur_index += 1

          if not guidesFactory.data.tissues_disabled # include median
            expression_data[cur_index].push (exon.expression.median / max_expression).toFixed(2)
            cur_index += 1

        tissue_index = 0
        while cur_index < $scope.chart_config.expression.series.length and tissue_index < guidesFactory.data.tissues.length
          tissue_name = guidesFactory.data.tissues[tissue_index]

          expression_data[cur_index].push (exon.expression[tissue_name] / max_expression).toFixed(2)
          tissue_index += 1
          cur_index += 1

      expression_labels.push('Exon ' + (key+1))
      guide_labels.push('') # empty labels

    $scope.chart_config.expression.labels = expression_labels
    $scope.chart_config.expression.data = expression_data

    $scope.chart_config.guides.labels = guide_labels
    $scope.chart_config.guides.data = [guides_data]

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
    non_targeting_guides_href = '/static/data/pre_processed/non_targeting_hum.json'
    if guidesFactory.data.genome == 'mus'
      non_targeting_guides_href = '/static/data/pre_processed/non_targeting_mus.json'

    $http.get(non_targeting_guides_href).then (res) ->
      non_targeting = res.data.data
      guidesCSV = $filter('filter')($scope.merged_gRNAs, {selected:true}, true)
      guidesCSV = $filter('orderBy')(guidesCSV, ['-gene','score'], true)
      padding = Math.floor(Math.log(guidesCSV.length) / Math.log(10)) + 1
      angular.forEach guidesCSV, (guide, idx) ->
        guide.uid = "GUIDES_sg" + pad(idx + 1, padding)
      stop_pos = $scope.non_targeting_guides_info.count - 1
      if stop_pos > 0
        guidesCSV.concat non_targeting[0..stop_pos]
      else
        guidesCSV
]
