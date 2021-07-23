# qtlPlots

### Install
```r
devtools::install_github("GabrielHoffman/qtlPlots")
```

### Description
Plotting functions take in a `GRanges` object where `scores` defines the value to be plotted on the y-axis:
 
 * `plotMht()` Manhattan plot for local region
 * `plotPosterior()` Posterior probabilities in local region

 In addition, there are 2 other functions:
 
 * `plotEnsGenes()` Plot ENSBMLE genes
 * `plotBed()` Plot genome intervals


Here is a simple example for simulated data:
 
```r
library(GenomicRanges)
library(qtlPlots)
 
# simple example of QTL results
df = data.frame(Chr         = 1,
                Position    = 1:100, 
                Variant     = paste0('SNP', 1:100),
                p.value     = runif(100))

# Convert to GRanges
gr = with(df,   GRanges(Chr, 
                IRanges(start   = Position, 
                        width   = 1, 
                        Variant = Variant, 
                        score   = -log10(p.value)))
 
# where: genome interval to show
wh = GRanges(1, IRanges(1, 100))			
 						
# example plot
plotMht(gr, wh, recombRate=FALSE) 
```
 
# Example on real caQTL data 

You can also add any other plots to this for eQTL and GWAS

```r
# test combined plotting
suppressPackageStartupMessages({
library(qtlPlots)
library(ggplot2)
library(GenomicRanges) 
})

ensdb = EnsDb.Hsapiens.v75::EnsDb.Hsapiens.v75 # genes on hg19
# ensdb = EnsDb.Hsapiens.v86::EnsDb.Hsapiens.v86 # genes on hg38

# read data
dir <- system.file("data", package="qtlPlots")	
gr_caQTL = readRDS(paste0(dir, '/', "gr_caQTL.RDS"))
gr_peaks = readRDS(paste0(dir, '/', "gr_peaks.RDS"))
gr_finemap = readRDS(paste0(dir, '/', "gr_finemap.RDS"))

# plot results for this peak
targetPeak = 'Peak_82668'

# where: get interval to plot
wh = range(gr_caQTL[gr_caQTL$Gene == targetPeak])

# only keep results for target peak
gr_caQTL_sub = gr_caQTL[which(gr_caQTL$Gene == targetPeak)]
gr_finemap = gr_finemap[gr_finemap$feature == targetPeak]

# Manhattan plot
################

# identify variants in the candidate set
gr_caQTL_sub$inCandidateSet = FALSE
gr_caQTL_sub$inCandidateSet[gr_caQTL_sub$Variant %in% gr_finemap$Variant] = TRUE

# recomb rate is hg19
fig_mht = plotMht( gr_caQTL_sub, wh )

# Posteriors
#############

fig_post = plotPosterior( gr_finemap, wh )

# Plot genes
############

fig_gene = plotEnsGenes( ensdb, wh)

# ATAC-seq peaks
#################

# highlight the target peak
gr_peaks$status = 0
gr_peaks$status[names(gr_peaks) == targetPeak] = 1

fig_bed = plotBed( gr_peaks, wh )

# Combine into multi-track plot
###############################

fig_track = ggbio::tracks(
     "caQTL"        = fig_mht, 
    "caQTL Fine map"= fig_post,
    'Gene'          = fig_gene,
    'ATAC'          = fig_bed,
    xlim = wh,
    padding = unit(-.65, "lines"), 
    label.bg.fill="navy", 
    label.text.color="white",
    heights=c(1, .6, .2, .1),
    title = targetPeak) 
fig_track@mutable['Genes'] = FALSE
fig_track
```

 