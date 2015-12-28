FlaskStart.factory 'GuidesFactory', ['$http', '$q', ($http, $q) ->
	class GuidesFactory
		constructor: () ->
			# set genes, guides, etc.

		generateGuides: (gene) ->
			# use that information to grab from server
			$http.get('/generate')

		generateExons: () ->
			$http.get('/generate')

]