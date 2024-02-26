table.on( 'row-reordered', function ( e, diff, edit ) {
    for ( var i=0, ien=diff.length ; i<ien ; i++ ) {
        $(diff[i].node).addClass("reordered");
    }
    console.log("test");
} );
