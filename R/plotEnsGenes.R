# Gabriel Hoffman
#
# October 2, 2020
# Adapt decorate::plotEnsGenes to plot using ggplot2 instead of grid

# adapted from decorate::get_exon_coords(), but return gene_id here
#' @importFrom GenomicFeatures exonsByOverlaps
get_exon_coords = function (ensdb, query, biotypes = c("protein_coding")){
    if (!is(ensdb, "EnsDb")) {
        stop("ensdb must be an ENSEMBL databsed of type EnsDb")
    }
    gr_exon = exonsByOverlaps(ensdb, query, columns = c("exon_seq_start", 
        "exon_seq_end", "symbol", "gene_biotype", "tx_seq_start", 
        "tx_seq_end", "tx_cds_seq_start", "tx_cds_seq_end", 'gene_id'))
    if (!is.na(biotypes)) {
        gr_exon = gr_exon[gr_exon$gene_biotype %in% biotypes]
    }
    gr_exon
}


#' Plot ENSEMBL genes
#' 
#' Plot ENSEMBL genes in region
#' 
#' @param ensdb ENSEMBL database object like EnsDb.Hsapiens.v86 
#' @param wh genome interval to plot
# @param plot_lines_distance veritcal distance between genes
# @param vp viewport
#' @param splice_variants if TRUE, show multiple transcripts from the same gene
#' @param non_coding if TRUE, also show non-coding genes
#' @param arrow.size arrow size for genes
#' @param ensGene show specific non-coding locus
#' @param size font size
#'
#' @return ggplot2 of genome region
#' 
#' @examples
#' library(EnsDb.Hsapiens.v86)
#' library(GenomicRanges) 
# library(grid) 
#' 
#' # gene database
#' ensdb = EnsDb.Hsapiens.v86
#' 
#' # interval
#' query = GRanges("20", IRanges(62045027,62164563))
#'
#' # plot genes
#' plotEnsGenes( ensdb, query)
#' 
#' @export
#' @import GenomicRanges ggplot2
#' @importFrom scales comma
# @import grid
#' @importFrom data.table data.table
# @importFrom grDevices pdf dev.off
plotEnsGenes = function(ensdb, wh, splice_variants = FALSE, non_coding = FALSE, arrow.size=0.05, ensGene = NULL, size=8){

    minRange = start(wh)
    maxRange = end(wh)
    chromosome = seqnames(wh)

	gr = GRanges(gsub( "^chr", "", chromosome), IRanges(minRange, maxRange))

	if( non_coding ){
		gr_exons = get_exon_coords( ensdb, gr, NA )
	}else{
		gr_exons = decorate::get_exon_coords( ensdb, gr, 'protein_coding' )

		# get locus for ensGene (even if it is non-coding)
		if( !is.null(ensGene) ){
			# get all loci
			gr_exons_focus = get_exon_coords( ensdb, gr, NA )

			# add locus of focus to other loci
			gr_exons_focus = gr_exons_focus[gr_exons_focus$gene_id == ensGene]
			gr_exons = unique(c(gr_exons, gr_exons_focus))
		}
	}

    # get gene coordinates
    if( non_coding ){
        biotype = NA
    }else{
        biotype = c("protein_coding")
    }

    if( length(gr_exons) > 0){
   
    	df = data.table(data.frame(gr_exons))

        gene_biotype = 0
        symbol = 0
        tx_cds_seq_end = 0
        tx_cds_seq_start = 0
        tx_seq_end = 0
        tx_seq_start = 0

    	if( !splice_variants ){
    		# single body per gene
    		suppressWarnings(df_wide <- df[,data.frame(
    			gene_name = unique(symbol),
    			chrom = unique(seqnames),
    			strand = unique(strand), 
    			exonStarts = paste(start, collapse=','), 
    			exonEnds = paste(end, collapse=','),
    			exonCount = length(start),
    			biotype = unique(gene_biotype),
    			txStart = min(tx_seq_start, na.rm=TRUE),
    			txEnd = max(tx_seq_end, na.rm=TRUE),
    			cdsStart = as.integer(min(tx_cds_seq_start, na.rm=TRUE)),
    			cdsEnd = as.integer(max(tx_cds_seq_end, na.rm=TRUE)),
    			stringsAsFactors=FALSE),by=c("symbol")])
    	}else{
    		stop("splice_variants = TRUE not currently supported ")
    		# multiple transcripts per gene		
            suppressWarnings(df_wide <- df[,data.frame(
    			gene_name = unique(symbol),
    			chrom = unique(seqnames),
    			strand = unique(strand), 
    			exonStarts = paste(start, collapse=','), 
    			exonEnds = paste(end, collapse=','),
    			exonCount = length(start),
    			biotype = unique(gene_biotype),
    			txStart = tx_seq_start,
    			txEnd = tx_seq_end,
                cdsStart = as.integer(min(tx_cds_seq_start, na.rm=TRUE)),
                cdsEnd = as.integer(max(tx_cds_seq_end, na.rm=TRUE)),
    			stringsAsFactors=FALSE),by=c("symbol", 'tx_seq_start', 'tx_seq_end')])
    	}
        t = data.frame(df_wide, stringsAsFactors=FALSE)
        t$plot_line <- 0
        t$plot_line[1] <- 1
    }else{
        t = data.frame()
    }
	
	# get min and max values based on strand
	t$txStart = pmax(t$txStart, minRange)
	t$txEnd = pmin(t$txEnd, maxRange)

	t$txMin = ifelse( t$strand == "+", t$txStart, t$txEnd )
	t$txMax = ifelse( t$strand == "+", t$txEnd, t$txStart )

	plot_lines_no <- 1

    if (dim(t)[1] > 1) {
        for(i in 2:dim(t)[1]){
            # gene_name_width <- Range/map_len * convertWidth(grobWidth(textGrob(paste(t[1,"gene_name"], "   "), gp = gpar(fontsize = 7))),
                # "npc", valueOnly = TRUE)
            gene_name_width = width(gr) / 10
            for (j in 1:plot_lines_no) {
                if (max(t[1:(i - 1), ][t[, "plot_line"] == j, "txEnd"], na.rm = TRUE) < t[i, "txStart"] - gene_name_width*2) {
                  t[i, "plot_line"] <- j
                  break
                }
            }
            if (!t[i, "plot_line"])
                t[i, "plot_line"] <- plot_lines_no <- plot_lines_no + 1
        }

       	t$plot_line = max(t$plot_line) + 1 - t$plot_line
    }   

    if( nrow(t) > 0 ){

		# genes using ggplot2
		ylim = c(min(t$plot_line), max(t$plot_line))
		d = ylim[2] - ylim[1]

		fig = ggplot(t, aes(x=txMin, xend=txMax, y=plot_line, yend=plot_line, color = ifelse(biotype=='protein_coding', '1', '2'), label=paste0('  ',symbol, '  '), hjust=ifelse(strand=='+', 1, 0))) + geom_segment( arrow = arrow(length = unit(arrow.size, "npc"), end="last",type = "open")) + ylim(ylim[1] - 0.3*d, ylim[2] + 0.1*d) + geom_text(size=2.2) + scale_x_continuous(label=comma, expand=c(0,0), limits=c(minRange, maxRange)) + scale_color_manual(values=c("navy", "grey40")) 
	}else{
		fig = ggplot() + scale_x_continuous(label=comma, expand=c(0,0), limits=c(minRange, maxRange))
	}

    fig + theme_bw() + theme(legend.position = "none", axis.title=element_blank(), axis.ticks.y=element_blank(), axis.line.y=element_blank(), panel.grid=element_blank(), axis.text.y=element_blank(),panel.background = element_blank(), strip.background = element_blank(), rect = element_rect(fill="white", linetype=0))
}

# chr20 47223127-48223127
# chromosome = "9"
# minRange = 19000000
# maxRange = 19300000


# ,panel.border = element_blank(), panel.grid.major = element_blank(),
# 		panel.grid.minor = element_blank(), axis.line = element_blank(),
# 		axis.title=element_blank(),
#         axis.text=element_blank(),
#         axis.ticks=element_blank())


