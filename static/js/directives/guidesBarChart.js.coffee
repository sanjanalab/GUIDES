Directives.directive 'guidesBarChart', ['$timeout', '$window', ($timeout, $window) ->
  replace: false
  restrict: 'E'
  scope: true # prototypically inherit scope, but create a new one
  templateUrl: 'static/partials/guides-bar-chart.html'
  link: (scope, elem, attrs) ->
    # resizeFunction = () ->
    #   scope.chartwidth = $('#exon_graph')[0].getBoundingClientRect().width - 200

    # angular.element($window).bind 'resize', () ->
    #   console.log scope.chartwidth
    #   resizeFunction()
    #   scope.$apply()

    # $timeout () ->
    #   resizeFunction()
]
