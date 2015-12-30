FlaskStart = angular.module 'FlaskStart', ['ui.bootstrap', 'FlaskStart.filters', 'FlaskStart.directives', 'ngSanitize', 'ngCsv', 'ui.select', 'ngMaterial', 'ngMessages', 'ngRoute', 'cgBusy']

FlaskStart.controller 'RootCtrl', ['$scope', ($scope) ->
]

FlaskStart.config ['$routeProvider', '$locationProvider', ($routeProvider, $locationProvider) ->
  $routeProvider.when("/designer", {templateUrl: "static/partials/designer.html", controller: 'DesignerCtrl'})
  .when("/", {templateUrl: "static/partials/home.html", controller: 'IndexCtrl'}) 
  .otherwise({redirectTo: "/s"})
  # $locationProvider.html5Mode(true);
]

Filters = angular.module 'FlaskStart.filters', []
Directives = angular.module 'FlaskStart.directives', []