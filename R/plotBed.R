
# Gabriel Hoffman
# July 22, 2021
# based on mmQTL_plots code 

#' Plot BED intervals
#'
#' Plot BED intervals
#'
#' @param grObj GRanges object with scores as the posterior inclusion probability
#' @param wh genome interval to plot
#' @param size font size
#'
#' @import ggplot2 GenomicRanges 
#' @importFrom scales comma
#' @export
plotBed = function( grObj, wh ){

    gr_sub = grObj[grObj %within% wh]

    # convert to data.frame
    df_loc = as.data.frame(gr_sub)

    if( ! ("status" %in% colnames(df_loc)) ){
        df_loc$status = 0
    }

    df_loc$status = factor(df_loc$status)

    col = c()
    if( any(df_loc$status==0)){
        col = append(col, "black")
    }
    if( any(df_loc$status==1)){
        col = append(col, "red")
    }

    df_loc$start = pmax(df_loc$start, start(wh))
    df_loc$end = pmin(df_loc$end, end(wh))

    ggplot(df_loc) + geom_segment(aes(x=start, y=1, xend = end, yend=1, color=status), size=4) + scale_x_continuous(expand=c(0,0), limits=c(start(wh), end(wh)), label=comma) + scale_y_continuous(expand=c(0,0), limits=c(0,1)) + scale_color_manual(values = col) + theme_bw() + theme(legend.position = "none", axis.title=element_blank(), axis.ticks.y=element_blank(), axis.line.y=element_blank(), panel.grid=element_blank(), axis.text.y=element_blank(),panel.background = element_blank(), strip.background = element_blank(), rect = element_rect(fill="white", linetype=0))
}



