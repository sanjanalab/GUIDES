FlaskStart = angular.module 'FlaskStart', ['ui.bootstrap', 'FlaskStart.filters', 'FlaskStart.directives', 'ngSanitize', 'ngCsv', 'ui.select', 'ngMaterial', 'ngMessages', 'ngRoute', 'cgBusy', 'ngFileUpload', 'uiSwitch', 'chart.js', 'angular-google-analytics']

FlaskStart.controller 'RootCtrl', ['$scope', ($scope) ->
]

FlaskStart.config ['$routeProvider', '$locationProvider', 'AnalyticsProvider', ($routeProvider, $locationProvider, AnalyticsProvider) ->
  AnalyticsProvider
    .logAllCalls(true)
    .setAccount('UA-72726433-1')
    .setDomainName('none')
    .useDisplayFeatures(true)

  factoryConfigured = ($location, $q, GuidesFactory) ->
    guidesFactory = new GuidesFactory
    deferred = $q.defer()
    if guidesFactory.data.genes.length > 0 or $location.path().split('?').length > 1
      deferred.resolve()
    else
      deferred.reject()
      $location.url('/')
    deferred.promise

  $routeProvider.when "/about", {
    templateUrl: "static/partials/about.html",
    controller: 'AboutCtrl'
  }
  .when "/designer/:task_id", {
    templateUrl: "static/partials/designer.html",
    controller: 'DesignerCtrl'
  }
  .when "/designer", {
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

# FlaskStart.run ['$rootScope', '$location', '$window', ($rootScope, $location, $window) ->
#   $rootScope.$on '$stateChangeSuccess', (event) ->
#     if not $window.ga
#       return
#     $window.ga('send', 'pageview', { page: $location.path() })
# ]

Filters = angular.module 'FlaskStart.filters', []
Directives = angular.module 'FlaskStart.directives', []
