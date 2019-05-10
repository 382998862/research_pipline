jQuery(window).load(function() {

    // apply supersubs and superfish for navigation menu
    jQuery('ul.topnav').supersubs({minWidth: 12, maxWidth: 27, extraWidth: 1}).superfish();
    jQuery('ul.topnav').superfish({
        delay: 1000,
        animation: {opacity:'false',height:'show'},
        speed: 'normal',
        autoArrows: false,
        dropShadows: false
    });

    // nivo slider
    jQuery('#slider').nivoSlider({
        effect: 'fade',
        slices:15,
        boxCols:8,
        boxRows:8,
        animSpeed:500,
        pauseTime:5000,
        directionNav:false,
        directionNavHide:false,
        controlNav:true});
});
