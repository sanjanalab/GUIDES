FlaskStart.controller 'IndexCtrl', ['$scope', '$http', '$timeout', 'GuidesFactory', 'Upload', ($scope, $http, $timeout, GuidesFactory, Upload) ->

  # Color of step 2 arrow
  $scope.designer_loading = false

  # Color change on hover (controlled from directive)
  $scope.learn_more_fill = "#4A4A4A"
  $scope.learn_more_to_black = () ->
    $scope.learn_more_fill = "#000"

  $scope.learn_more_to_white = () ->
    $scope.learn_more_fill = "#4A4A4A"

  # tissuesAvailable and genesAvailable is prepared here
  # selected results is placed directly back into factory
  # this way, we can just push to next controller and everything is ready :)
  $scope.guidesFactory = new GuidesFactory()
  $scope.guidesFactory.data.gtex_enabled = true

  # Validations
  $scope.genesWarning = false
  $scope.tissuesWarning = false
  $scope.emailWarning = false

  emailPattern = /// ^ # beginning of line
   ([\w.-]+)           # one or more letters, numbers, _ . or -
   @                   # followed by an @ sign
   ([\w.-]+)           # then one or more letters, numbers, _ . or -
   \.                  # followed by a period
   ([a-zA-Z.]{2,6})    # followed by 2 to 6 letters or periods
   $ ///i              # end of line and ignore case

  # Should the next page load?
  $scope.$on '$locationChangeStart', (event, newUrl) ->
    urlComps = newUrl.split('/')
    if urlComps[urlComps.length - 1] == "designer"
      if $scope.guidesFactory.data.genes.length < 1
        $scope.genesWarning = true
        event.preventDefault()
      else
        $scope.genesWarning = false
      if not $scope.guidesFactory.data.tissues_disabled and $scope.guidesFactory.data.tissues.length < 1
        $scope.tissuesWarning = true
        event.preventDefault()
      else
        $scope.tissuesWarning = false
      if ($scope.guidesFactory.data.email_address) and ($scope.guidesFactory.data.email_address != '') and (not $scope.guidesFactory.data.email_address.match emailPattern)
        $scope.emailWarning = true
        event.preventDefault()
      else
        $scope.emailWarning = false
      if not ($scope.genesWarning or $scope.tissuesWarning or $scope.emailWarning)
        $scope.designer_loading = true

  # file upload
  $scope.$watch 'file', () ->
    if $scope.file
      if not $scope.file.$error
        # read file and give contents to factory for parsing later
        reader = new FileReader()
        reader.onload = (e) ->
          text = reader.result
          text = text.replace(/[\"\'=>]+/g,'')
          newGenesFromFile = text.split(/[\s\t,;#\| ]+/)
          $scope.guidesFactory.data.genesFromFile = newGenesFromFile

          # change the text for file uploaded
          # add "Genes from file..." button if necessary
          displayTxt = "Genes from file..."
          if newGenesFromFile.length > 2
            displayTxt = newGenesFromFile[0] + ", " + newGenesFromFile[1] + ", ..."
          else if newGenesFromFile.length == 2
            displayTxt = newGenesFromFile[0] + ", " + newGenesFromFile[1]
          else if newGenesFromFile.length == 1
            displayTxt = newGenesFromFile[0]
          $scope.guidesFactory.noteFileUploaded(displayTxt)

        reader.readAsText($scope.file) # assume UTF-8 encoding
      else
        console.log $scope.file.$error

  # Refresh on genome change
  $scope.$watch (() -> $scope.guidesFactory.data.genome), (current, original) ->
    console.log "changed genome to " + current
    $scope.guidesFactory.data.genes = []
    $scope.guidesFactory.data.rejected_genes = []
    if current == 'hum' # human
      $scope.guidesFactory.available.genes = $scope.guidesFactory.available.human_genes
    else # mouse
      $scope.guidesFactory.available.genes = $scope.guidesFactory.available.mouse_genes

  # Slider changes
  $scope.md_slider_quantity = 20
  slider_vals = [1,  2,  3,  4,  5,  6,  7,  8,  9, 10, 11, 12, 15, 20, 25, 30, 35, 40, 45, 50, 60, 70, 80, 90, 100]
  $scope.$watch 'md_slider_quantity', (value) ->
    if value < 100
      $scope.guidesFactory.data.quantity = slider_vals[Math.floor(value/4)]
    else
      $scope.guidesFactory.data.quantity = 100
]
