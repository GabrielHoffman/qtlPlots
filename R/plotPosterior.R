# Gabriel Hoffman
# July 22, 2021
# based on mmQTL_plots code 


#' Zoom in manhattan plot of posteriors
#'
#' Zoom in manhattan plot of posteriors
#'
#' @param grObj GRanges object with scores as the posterior inclusion probability
#' @param wh genome interval to plot
#' @param size font size
#'
#' @import ggplot2 GenomicRanges 
#' @importFrom scales comma
#' @export
plotPosterior = function( grObj, wh, size=8 ){

	# extact scores *within* the region wh
    gr_sub = grObj[grObj %within% wh]
    df_sub = data.frame(Position=start(gr_sub), score=score(gr_sub))

	fig = ggplot(df_sub, aes(Position, score)) + theme_classic(size)

	fig + geom_point(size=.4, color="red") + ylab("Posterior") + scale_y_continuous(expand=c(0,0.01), limits=c(0, 1)) + scale_x_continuous(expand=c(0,0), limits=c(start(wh), end(wh)), label=comma) + theme(legend.position="none") 
}	
