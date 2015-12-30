FlaskStart.controller 'IndexCtrl', ['$scope', '$http', '$timeout', 'GuidesFactory', ($scope, $http, $timeout, GuidesFactory) ->

  # tissuesAvailable and genesAvailable is prepared here
  # selected results is placed directly back into factory
  # this way, we can just push to next controller and everything is ready :)
  guidesFactory = new GuidesFactory()
  $scope.guidesFactory = guidesFactory

  $scope.tissuesAvailable = ['Thyroid', 'Testis', 'Cervix Uteri', 'Adipose Tissue', 'Breast', 'Vagina', 'Nerve', 'Pituitary', 'Stomach', 'Fallopian Tube', 'Bone Marrow', 'Bladder', 'Blood', 'Colon', 'Prostate', 'Pancreas', 'Blood Vessel', 'Liver', 'Spleen', 'Small Intestine', 'Uterus', 'Ovary', 'Muscle', 'Heart', 'Adrenal Gland', 'Brain', 'Salivary Gland', 'Lung', 'Skin', 'Esophagus', 'Kidney']
  guidesFactory.data.tissues = [
    $scope.tissuesAvailable[8]
  ]

  $scope.genesAvailable = []
  # We are in /static/js/min/scripts.min.js
  # We want /static/data/pre_processed/genes_list.json
  $http.get('/static/data/pre_processed/genes_list.json').then (res) ->
    $scope.genesAvailable = res.data
    guidesFactory.data.genes = [
      $scope.genesAvailable[5]
      $scope.genesAvailable[4]
    ]

  $scope.guidesFactory.data.quantity = 60
]