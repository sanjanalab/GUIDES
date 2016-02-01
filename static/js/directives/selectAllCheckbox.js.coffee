Directives.directive 'selectAllCheckbox', () ->
  replace: true
  restrict: 'E'
  scope: {
    checkboxes: '='
    counter: '='
    chartdata: '='
    guidechange: '&'
  }
  template: '<input type="checkbox" ng-model="master" ng-change="masterChange()">'
  controller: ($scope, $element) ->
    $scope.masterChange = () ->
      if $scope.master
        angular.forEach $scope.checkboxes, (cb, idx) ->
          cb.selected = true
          $scope.guidechange({guide:cb})
      else
        angular.forEach $scope.checkboxes, (cb, idx) ->
          cb.selected = false
          $scope.guidechange({guide:cb})

    $scope.$watch 'checkboxes', (() ->
      allSet = true
      allClear = true
      angular.forEach $scope.checkboxes, (cb, idx) ->
        if cb.selected
          allClear = false
        else
          allSet = false
      $element.prop('indeterminate', false)
      if allSet
        $scope.master = true
      else if allClear
        $scope.master = false
      else
        $scope.master = false
        $element.prop('indeterminate', true)
    ), true