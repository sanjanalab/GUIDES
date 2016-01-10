Directives.directive 'selectAllCheckbox', () ->
  replace: true
  restrict: 'E'
  scope: {
    checkboxes: '='
    counter: '='
  }
  template: '<input type="checkbox" ng-model="master" ng-change="masterChange()">'
  controller: ($scope, $element) ->
    $scope.masterChange = () ->
      if $scope.master
        angular.forEach $scope.checkboxes, (cb, idx) ->
          if cb.selected == false
            $scope.counter += 1
          cb.selected = true
      else
        angular.forEach $scope.checkboxes, (cb, idx) ->
          if cb.selected == true
            $scope.counter -= 1
          cb.selected = false

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