Directives.directive 'drawindividualexon', () ->
  replace: false
  restrict: 'A'
  scope: {
    index: '@'
    exonsLength: '@'
    exon: '='
  }
  #template: '<div></div>'
  templateUrl: 'static/partials/individual-exon.html'
  link: (scope, elem, attrs) ->
    # x - 2x - x - 2x - ... - 2x - x
    # Do the math - you'll find we get 3kx + x = x(3k+1) = total width
    # x = (total width) / (3k+1)
    # total width comes from parent of parent, which is the svg.
    scope.svgUnit = elem.parent().parent()[0].getBoundingClientRect().width / (3 * scope.exonsLength + 1)
    scope.toInt = (num) ->
      parseInt(num, 10)
