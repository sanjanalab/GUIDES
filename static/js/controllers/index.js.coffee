FlaskStart.controller 'IndexCtrl', ['$scope', '$http', '$timeout', 'GuidesFactory', ($scope, $http, $timeout, GuidesFactory) ->

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
]