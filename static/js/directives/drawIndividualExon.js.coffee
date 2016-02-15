Directives.directive 'drawindividualexon', ['$timeout', '$window', ($timeout, $window) ->
  replace: false
  restrict: 'A'
  scope: {
    index: '@'
    exonsLength: '@'
    y: '@'
    exon: '='
    modifySvgUnit: '&' # hacky way to get this to change svgUnit
    exonHovered: '='
  }
  #template: '<div></div>'
  templateUrl: 'static/partials/individual-exon.html'
  link: (scope, elem, attrs) ->
    # x - 2x - x - 2x - ... - 2x - x
    # Do the math - you'll find we get 3kx + x = x(3k+1) = total width
    # x = (total width) / (3k+1)
    # total width comes from parent of parent, which is the svg.
    scope.toInt = (num) ->
      parseInt(num, 10)

    scope.hoveringExon = false
    # scope.hoveringExonFunc = () ->
    #   if scope.exonHovered == scope.toInt(scope.index)
    #     true
    #   else
    #     false
    #scope.hoveringExonFunc = (scope.exonHovered == scope.toInt(scope.index))

    scope.gRNA_count_y = scope.toInt(scope.y) + 23
    gRNA_count_y_backup = scope.gRNA_count_y

    resizeFunction = () ->
      scope.svgUnit = ($('#exon_graph')[0].getBoundingClientRect().width - 20 - 20) / (3 * scope.exonsLength + 1)
      scope.modifySvgUnit({unit: scope.svgUnit})
      if 2 * scope.svgUnit <= 41
        scope.gRNA_count_y = 31
      else
        scope.gRNA_count_y = gRNA_count_y_backup

    angular.element($window).bind 'resize', () ->
      resizeFunction()
      scope.$apply()

    # scope.$watch 'exonHovered', (val) ->
    #   if val == scope.toInt(scope.index) + 1
    #     scope.hoveringExon = true
    #   else
    #     scope.hoveringExon = false

    # Kick things off
    resizeFunction()
]