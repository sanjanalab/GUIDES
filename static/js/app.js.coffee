FlaskStart = angular.module 'FlaskStart', ['ui.bootstrap', 'FlaskStart.filters', 'FlaskStart.directives', 'ngSanitize', 'ngCsv', 'ui.select', 'ngMaterial', 'ngMessages', 'ngRoute', 'cgBusy', 'ngFileUpload', 'uiSwitch', 'chart.js']

FlaskStart.controller 'RootCtrl', ['$scope', ($scope) ->
]

FlaskStart.config ['$routeProvider', '$locationProvider', ($routeProvider, $locationProvider) ->
  factoryConfigured = ($location, $q, GuidesFactory) ->
    guidesFactory = new GuidesFactory
    deferred = $q.defer()
    if guidesFactory.data.genes.length > 0
      deferred.resolve()
    else
      deferred.reject()
      $location.url('/')
    deferred.promise

  $routeProvider.when "/designer", {
    templateUrl: "static/partials/designer.html", 
    controller: 'DesignerCtrl', 
    resolve: {
      factoryConfigured: factoryConfigured,
    }
  }
  .when "/", {
    templateUrl: "static/partials/home.html", 
    controller: 'IndexCtrl'
  }
  .otherwise {
    redirectTo: "/"
  } 
]

Filters = angular.module 'FlaskStart.filters', []
Directives = angular.module 'FlaskStart.directives', []