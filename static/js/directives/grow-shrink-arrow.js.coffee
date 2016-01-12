Directives.directive 'growshrinkarrow', ['$interval', ($interval) ->
  restrict: 'A'
  link: (scope, elem, attrs) ->
    elem.height('75px')
    elem.width('75px')

    elem.mouseenter ->
      elem.animate({
        height: '85px',
        width: '85px'
      }, {
        duration: 50
      })

    elem.mouseleave ->
      elem.animate({
        height: '75px',
        width: '75px'
      }, {
        duration: 50
      })
]


    #           <svg height="250" id="exon_graph" width="100%" x="20">
    #             <g id="background_exon_graph">
    #               <rect fill="#eee" height="250" width="100%" x="0" y="0"></rect>
    #             </g>
    #             <g id="all_exons">
    #               <rect class="connector" fill="#000" height="3" stroke="#000" width="100%" x="0" y="100" id="connector_rect"></rect>

    #               <g ng-repeat="exon in gene.exons">
    #                 <rect class="exon" fill="#fff" height="50" ng-attr-width="{{2 * svgUnit()}}" stroke="#aaa" ng-attr-x="{{ svgUnit() + $index*3*svgUnit() }}" y="75"></rect>

    #                 <text class="gRNA_count" fill="#111" font-size="16" ng-attr-x="{{svgUnit()*1.8  + $index*3*svgUnit() }}" y="108">
    #                   {{ (exon.gRNAs | filter:{selected:true}).length }}
    #                 </text>

    #                 <text class="exon_label" fill="#111" font-size="8" ng-attr-x="{{svgUnit()*1.6  + $index*3*svgUnit() }}" y="145">
    #                   Exon {{ $index + 1 }}
    #                 </text>
    #               </g>
    #             </g>
    #           </svg>



    # $scope.svgUnit = () -> 
    #   $element.find('#connector_rect')[0].getBBox().width / (3 * $scope.gene.exons.length + 1)
