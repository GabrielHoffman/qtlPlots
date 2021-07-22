# Gabriel Hoffman
#
# October 2, 2020
# Intersect mmQTL with CAUSALdb

# library(RColorBrewer)
# library(EnsDb.Hsapiens.v75)

# Plot results
##############

make_plot = function( ensGene, window = 5e5, ord=1, non_coding=FALSE, showConditional=FALSE ){

	# geneSymbol = tryCatch(
	# 		mapIds(org.Hs.eg.db,
	#                      keys=ensGene,
	#                      column="SYMBOL",
	#                      keytype="ENSEMBL")
	# 		, error = function(e) '')

	# get gene name from ENSEMBL id
	geneInfo = select(EnsDb.Hsapiens.v75, keys=ensGene, keytype="GENEID", column=c('GENENAME', 'GENEBIOTYPE'))
	geneSymbol = geneInfo$GENENAME

	# get location of top fine-mapped eQTL SNP	
	df_fm = df_finemap[(Gene == ensGene) & (eQTL_order == ord),]
	pos_center = df_fm[which.max(PIP),Position]
	wh = GRanges(paste0("chr", df_fm$Chr[1]), IRanges(pos_center - window, pos_center + window))

	# eQTL
	if( (ord == 1) || ! showConditional ){		
		df_window = df_eqtl[(Gene == ensGene) & (eQTL_order == ord),]
	}else{		
		df_window = df_eqtl[(Gene == ensGene) & (eQTL_order <= ord),]
	}

	# get candidate causal set
	df_window = df_window[(Variant %in% df_fm$Variant) & (eQTL_order == ord),inCandSet := "yes"]
	df_window$inCandSet[is.na(df_window$inCandSet)]  = 'no'

	# recombination rate
	# expand window until a recombination measurement is found 
	scl = 1
	n = 0
	while( n ==0 ){
		df_recomb_sub = df_recomb[(position > start(wh)/scl) & (position < end(wh)*scl) & (chromosome == df_window$Chr[1]),]
		n = nrow(df_recomb_sub)
		scl = scl * 3
	}

	# Add recombination rate to first plot
	# linear interpolation evaluate at points in resComposite
	app = approx( df_recomb_sub$position, df_recomb_sub$recomb_rate, df_window$Position )
	df_window$rate = pmin(app$y, 100) # cap local rate at 100
	df_window$log10.p.value = pmin(300,df_window$log10.p.value)

	ymax = max(df_window$log10.p.value)
	ymax.rate = ymax
	
	if( (ord == 1) || ! showConditional ){
		fig_eqtl = ggplot(df_window, aes(Position, log10.p.value, color=inCandSet)) + scale_color_manual(values = c("black", "red"))
	}else{		
		df_window[,cat := paste(inCandSet, eQTL_order < ord)]
		df_window$cat = factor(df_window$cat, c("no TRUE", "no FALSE", "yes FALSE"))
		setorder(df_window, eQTL_order, cat)
		fig_eqtl = ggplot(df_window, aes(Position, log10.p.value, color=cat)) + scale_color_manual(values = c("grey", "black", "red"))
	}

	fig_eqtl = fig_eqtl + geom_point(size=.4) + ggtitle("FINEMAP - gene") + ylab(bquote(-log[10]~P)) + scale_x_continuous(expand=c(0,0)) + scale_y_continuous(expand=c(0,0), limits=c(0, ymax*1.05), sec.axis = sec_axis(~./(ymax.rate/100), name = "Recombination\nrate [cM/Mb]")) + geom_line(aes(y=rate*(ymax.rate/100)), color="dodgerblue", size=.5) + theme(legend.position="none")

	# eQTL fine mapping
	# df_fm = df_finemap[(Gene == ensGene) & (eQTL_order == ord),]
	#, color=as.character(eQTL_order)
	fig_finemap_gene = ggplot(df_fm, aes(Position, PIP, size=PIP^2)) + geom_point(color="red") + ggtitle("Finemap - gene") + ylab("Posterior") + ylim(0,1) + scale_x_continuous(expand=c(0,0)) + scale_size_continuous(limits = c(0, 1), range(0.1, 6)) 

	# CAUSALdb
	# c('OT294', 'OT344', 'OT311', 'GA553', 'OT310')
	df_pip = get_credible_set( causalDB, unique(df_show$meta_id), wh)
	setorder(df_pip, Trait, -FINEMAP)

	# intersect finemapping from gene (df_fm) and GWAS (df_pip)
	df_coloc = merge(df_fm, df_pip, by.x="Position", by.y="BP")
	df_select = df_coloc[, data.frame(max=max(FINEMAP*PIP)),by='meta_id']

	# only keep GWAS that shares a causal variant with eQTl at least 0.01 
	df_coloc = df_coloc[meta_id %in% df_select[max>0.01,meta_id],]

	# keep GWAS fine mapping for these traits
	df_pip = df_pip[(meta_id %in% df_coloc$meta_id) | (Trait %in% df_coloc$Trait),]
	df_pip$Trait = factor(df_pip$Trait)

	if( nrow(df_pip) > 0){

		if( nlevels(df_pip$Trait) < 8 ){
			cols = brewer.pal(n = nlevels(df_pip$Trait), name = "Dark2")
		}else{
			cols = rainbow( nlevels(df_pip$Trait) )
		}

		fig_causal = ggplot(df_pip, aes(BP, FINEMAP, color=Trait, size=FINEMAP^2)) + geom_point() + ylim(0,1) + scale_x_continuous(expand=c(0,0)) + ylab("Posterior") + scale_size_continuous(limits = c(0, 1), range(0.1, 6)) + scale_color_manual(values = cols) # + geom_text_repel(data=subset(df_pip, !duplicated(Trait)), aes(label=Trait), size=3, force=100, nudge_x = .1, nudge_x=window/2)
	}else{
		fig_causal = ggplot() + scale_x_continuous(expand=c(0,0), limits=c(start(wh), end(wh))) + ylim(0,1) + ylab("Posterior")
	}

	# plot coloc
	# df_coloc = merge(df_fm, df_pip, by.x="Position", by.y="BP")
	df_coloc = df_coloc[FINEMAP*PIP>0,]
	if( nrow(df_coloc) > 0){
		df_coloc$Trait = factor(df_coloc$Trait, levels(df_pip$Trait))
		df_coloc[,prob.prod := FINEMAP*PIP,]
		setorder(df_coloc, Trait, -prob.prod)
		ymax = max(df_coloc[,prob.prod])
		fig_coloc = ggplot(df_coloc[prob.prod>0,], aes(Position, FINEMAP*PIP, color=Trait, size=prob.prod^2)) + geom_point() + scale_x_continuous(expand=c(0,0)) + ylab("CLPP") + scale_size_continuous(limits = c(0, 1), range(0.1, 6)) + ylim(0,1) + scale_color_manual(values = cols) + geom_text_repel(data=subset(df_coloc[prob.prod>0,], !duplicated(Trait)), aes(label=Trait), size=3, force=1, nudge_y = .1, nudge_x=window/2, hjust=0)
	}else{	
		fig_coloc = ggplot() + scale_x_continuous(expand=c(0,0), limits=c(start(wh), end(wh))) + ylim(0,1) + ylab("CLPP")
	}

	# Gene models     
	fig_genebody = plotEnsGenes_gg( EnsDb.Hsapiens.v75, start(wh), end(wh), seqnames(wh), splice_variants=FALSE, non_coding=non_coding, ensGene=ensGene) 
	
	# Combine plots
	fig_track = ggbio::tracks( "eQTL"		= fig_eqtl, 
		"eQTL\nFine map"	= fig_finemap_gene,
		"GWAS\nFine map" 	= fig_causal,
		"Shared\ncandidates"= fig_coloc,
		"Genes"				= fig_genebody,
		xlim = wh,
		padding = unit(-.65, "lines"),
		label.bg.fill="navy", label.text.color="white",
		heights=c(.9,.5,.5,.5,.4),
		theme = theme_bw(8) + theme(legend.position="none", panel.grid.minor = element_blank()),
			title=paste0( geneSymbol, ' (',ensGene,')') )
	fig_track@mutable['Genes'] = FALSE
	fig_track
}


  # tracks(fig_eqtl, fig_finemap_gene, xlim=wh)









