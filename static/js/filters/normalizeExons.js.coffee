Filters = angular.module 'FlaskStart.filters', []

Filters.filter 'normalizeExonsFilter', ->
  (exons, gene_length) ->
    pixel_width = 800 #length of displayed gene

    normalizedExons = []

    angular.forEach exons, (value, key) ->
      newExon = {}
      angular.copy(value, newExon)

      newExon.start = value.start / gene_length * pixel_width
      newExon.end = value.end / gene_length * pixel_width
      normalizedExons.push(newExon)

      # newExon = {}
      # newExon["start"] = value.start / gene_length * pixel_width
      # newExon["end"] = value.end / gene_length * pixel_width

      # newGRNAs = []

      # angular.forEach value["gRNAs"], (value_, key_) ->


      # newExon["gRNAs"] = value["gRNAs"]

      # normalizedExons.push ({
      #   "start": value.start / gene_length * pixel_width
      #   "end": value.end / gene_length * pixel_width
      # })
    normalizedExons

    # {
    #   "start": exon.start / gene_length * pixel_width,
    #   "end": exon.end / gene_length * pixel_width
    # }
    	# exon.start = exon.start / gene_length * pixel_width
    	# exon.end = exon.end / gene_length * pixel_width
    	# exon