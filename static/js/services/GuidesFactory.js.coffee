FlaskStart.factory 'GuidesFactory', ['$http', '$q', ($http, $q) ->
	class GuidesFactory
		data: 
			'genes' : []
			'tissues' : []
			'quantity': 25

		constructor: () ->
			# set genes, guides, etc.

		generateGuides: () ->
			# use that information to grab from server
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