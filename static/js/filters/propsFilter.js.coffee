Filters.filter 'propsFilter', ->
  (items, props) ->
    out = []
    if angular.isArray(items)
      items.forEach (item) ->
        itemMatches = false
        keys = Object.keys(props)
        i = 0
        while i < keys.length
          prop = keys[i]
          text = props[prop].toLowerCase()
          # Check if the property starts with the query
          if item[prop].toString().toLowerCase().indexOf(text) == 0
            itemMatches = true
            break
          i++
        if itemMatches
          out.push item
        return
    else
      # Let the output be the input untouched
      out = items
    
    out

# Used for our Angular select library