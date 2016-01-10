FlaskStart.factory 'GuidesFactory', ['$http', '$q', ($http, $q) ->
    class GuidesFactory
        available: # can select
            'genes': [] # set later by $http in constructor
            'tissues': ['Thyroid', 'Testis', 'Cervix Uteri', 'Adipose Tissue', 'Breast', 'Vagina', 'Nerve', 'Pituitary', 'Stomach', 'Fallopian Tube', 'Bone Marrow', 'Bladder', 'Blood', 'Colon', 'Prostate', 'Pancreas', 'Blood Vessel', 'Liver', 'Spleen', 'Small Intestine', 'Uterus', 'Ovary', 'Muscle', 'Heart', 'Adrenal Gland', 'Brain', 'Salivary Gland', 'Lung', 'Skin', 'Esophagus', 'Kidney']

        data: # currently selected
            'genes' : []
            'tissues' : []
            'quantity': 60

        constructor: () ->
            # Setup default available, and default selected
            # We are in /static/js/min/scripts.min.js
            # We want /static/data/pre_processed/genes_list.json
            this_ = this
            $http.get('/static/data/pre_processed/genes_list.json').then (res) ->
                this_.available.genes = res.data
                this_.data.genes = [
                    this_.available.genes[28284]
                    this_.available.genes[494]
                ]

        generateGuides: () ->
            if this.data.tissues.length == 0
                this.data.tissues = this.available.tissues
            $http {
                url: '/generate'
                method: 'POST'
                headers:
                    'Content-Type': 'application/json'
                data: JSON.stringify(this.data)
            }
            .success (data) ->
                console.log(data)

        generateExons: () ->
            $http.get('/generate')

        setGenes: (genes) ->
            this.data.genes = genes

        setTissues: (tissues) ->
            this.data.tissues = tissues

]

# scope.postTest = function(){

#     var data = [obj1, obj2, obj3];
#     var jsonData=angular.toJson(data);
#     var objectToSerialize={'object':jsonData};

#     $http({
#         url: 'myURL',
#         method: "POST",
#         data: $.param(objectToSerialize),
#         headers: {
#                  'Content-Type': 'application/x-www-form-urlencoded; charset=UTF-8'
#         }
#     }).success(function(data){
#         alert("done");
#     });
# }