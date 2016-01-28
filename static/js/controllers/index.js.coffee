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

  $scope.$on '$locationChangeStart', (event) ->
    if $scope.guidesFactory.data.genes.length < 1
      $scope.genesWarning = true
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
          newGenesFromFile = text.split(',')
          $scope.guidesFactory.genesFromFile = $scope.guidesFactory.genesFromFile.push.apply($scope.guidesFactory.genesFromFile, newGenesFromFile)
          console.log "Just set genesFromFile!"
          console.log $scope.guidesFactory.genesFromFile 
        reader.readAsText($scope.file) # assume UTF-8 encoding
      else
        console.log $scope.file.$error

  # From pre-calculated experimental results
  $scope.waitTime = () ->
    3 + Math.max(30,($scope.guidesFactory.data.quantity - 100)) // 8 * $scope.guidesFactory.data.genes.length
]