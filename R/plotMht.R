# Gabriel Hoffman
# July 22, 2021
# based on mmQTL_plots code 

#' Zoom in manhattan plot
#'
#' Zoom in manhattan plot
#'
#' @param grObj GRanges object with scores as the -log10 p-value
#' @param wh genome interval to plot
#' @param size font size
#'
#' @import ggplot2 GenomicRanges GenomeInfoDb
#' @importFrom scales comma
#' @importFrom stats approx
#' @export
plotMht = function( grObj, wh, size=8 ){

	# load recombination rate
	dir <- system.file("data", package="qtlPlots")
	# dir = "/Users/gabrielhoffman/workspace/repos/qtlPlots/data/"
    gr_recomb <- readRDS(file.path(dir,"gr_recomb_hg19.RDS"))

    # keep only entries on the same chromosome as grObj 
    # gr_recomb = gr_recomb[seqnames(gr_recomb) %in% seqnames(wh)] 
    # gr_recomb = gr_recomb[S4Vectors::match(seqnames(gr_recomb), seqnames(wh))]
    # gr_recomb = subset(gr_recomb, array(seqnames(gr_recomb)==unique(seqnames(wh))))
    gr_recomb = gr_recomb[as.character(seqnames(gr_recomb)) %in% as.character(seqnames(wh))]

    # extact scores *within* the region wh
    gr_sub = grObj[grObj %within% wh]
    df_sub = data.frame(Position=start(gr_sub), score=score(gr_sub), inCandidateSet=gr_sub$inCandidateSet)

    # recombination rate
	# expand window until a recombination measurement is found 
	wh_query = wh
	n = 0
	while( n < 1 ){
		# get intersection
		idx = which(gr_recomb %within% wh)

		# count number
		n = length(idx)

		# expand interval
		wh_query = flank(wh_query, width=width(wh_query)*2, both=TRUE)
	}

	# then expand one more step in each direction to avoid edge effects
	# but avoid steping over the last element
	idx_augment = c(min(idx)-1, idx, max(idx)+1)
	idx_augment = pmin(pmax(1, idx_augment), length(gr_recomb))
	gr_recomb_sub = gr_recomb[idx_augment] 

	# Add recombination rate to first plot
	# linear interpolation evaluate at points in resComposite
	app = approx( start(gr_recomb_sub), gr_recomb_sub$recomb_rate, df_sub$Position )
	df_sub$rate = pmin(app$y, 100) # cap local rate at 100
	# df_window$log10.p.value = pmin(300,df_window$log10.p.value)

	ymax = max(df_sub$score)
	ymax.rate = ymax

	# if( (ord == 1) || ! showConditional ){
		fig = ggplot(df_sub, aes(Position, score, color=inCandidateSet)) + scale_color_manual(values = c("black", "red")) + theme_classic(size)
	# }else{		
	# 	df_window[,cat := paste(inCandSet, eQTL_order < ord)]
	# 	df_window$cat = factor(df_window$cat, c("no TRUE", "no FALSE", "yes FALSE"))
	# 	setorder(df_window, eQTL_order, cat)
	# 	fig_eqtl = ggplot(df_window, aes(Position, log10.p.value, color=cat)) + scale_color_manual(values = c("grey", "black", "red"))
	# }

	fig = fig + geom_point(size=.4) + ylab(bquote(-log[10]~P)) + scale_x_continuous(expand=c(0,0), label=comma, limits=c(start(wh), end(wh))) + scale_y_continuous(expand=c(0,0), limits=c(0, ymax*1.05), sec.axis = sec_axis(~./(ymax.rate/100), name = "Recombination\nrate [cM/Mb]")) + geom_line(aes(y=rate*(ymax.rate/100)), color="dodgerblue", size=.5) + theme(legend.position="none")

	fig
}