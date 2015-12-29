Directives = angular.module 'FlaskStart.directives', []

Directives.directive 'growshrinkarrow', ['$interval', ($interval) ->
  restrict: 'A'
  link: (scope, elem, attrs) ->
    elem.mouseenter ->
      elem.animate({
        height: '85px',
        width: '85px'
      }, {
        duration: 50
      })

    elem.mouseleave ->
      elem.animate({
        height: '75px',
        width: '75px'
      }, {
        duration: 50
      })
]