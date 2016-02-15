FlaskStart.controller 'IndexCtrl', ['$scope', '$http', '$timeout', 'GuidesFactory', 'Upload', ($scope, $http, $timeout, GuidesFactory, Upload) ->

  # Color change on hover (controlled from directive)
  $scope.learn_more_fill = "#fff"
  $scope.learn_more_to_black = () ->
    $scope.learn_more_fill = "#000"

  $scope.learn_more_to_white = () ->
    $scope.learn_more_fill = "#fff"

  # tissuesAvailable and genesAvailable is prepared here
  # selected results is placed directly back into factory
  # this way, we can just push to next controller and everything is ready :)
  $scope.guidesFactory = new GuidesFactory()
  $scope.guidesFactory.data.gtex_enabled
  $scope.gtex_enabled = true

  # Validations
  $scope.genesWarning = false
  $scope.tissuesWarning = false

  $scope.$on '$locationChangeStart', (event, newUrl) ->
    urlComps = newUrl.split('/')
    if urlComps[urlComps.length - 1] == "designer"
      if $scope.guidesFactory.data.genes.length < 1
        $scope.genesWarning = true
        event.preventDefault()
      if not $scope.guidesFactory.data.tissues_disabled and $scope.guidesFactory.data.tissues.length < 1
        $scope.tissuesWarning = true
        event.preventDefault()

  # file upload
  $scope.$watch 'file', () ->
    if $scope.file
      if not $scope.file.$error
        # add "Genes form file..." button if necessary
        if not $scope.guidesFactory.hasUploadedFile()
          $scope.guidesFactory.noteFileUploaded()

        # read file and give contents to factory for parsing later
        reader = new FileReader()
        reader.onload = (e) ->
          text = reader.result
          newGenesFromFile = text.split(/[\s\t,#| ]+/)
          $scope.guidesFactory.genesFromFile = $scope.guidesFactory.genesFromFile.push.apply($scope.guidesFactory.genesFromFile, newGenesFromFile)
        reader.readAsText($scope.file) # assume UTF-8 encoding
      else
        console.log $scope.file.$error

  # Slider changes
  $scope.md_slider_quantity = 20
  slider_vals = [1,  2,  3,  4,  5,  6,  7,  8,  9, 10, 11, 12, 15, 20, 25, 30, 35, 40, 45, 50, 60, 70, 80, 90, 100]
  $scope.$watch 'md_slider_quantity', (value) ->
    if value < 100
      $scope.guidesFactory.data.quantity = slider_vals[Math.floor(value/4)]
    else
      $scope.guidesFactory.data.quantity = 100

  # From pre-calculated experimental results
  $scope.waitTime = () ->
    3 + Math.max(30,($scope.guidesFactory.data.quantity - 100)) // 8 * $scope.guidesFactory.data.genes.length
]