library(ggplot2)
library(gridExtra)

options(stringsAsFactors = FALSE)
args = commandArgs(trailingOnly=TRUE)

celltypes <- args[1] ## Tab delimited file with celltype in first column and number of replicates in the second. Should include DNA/plasmid
count_table <- args[2] ## Barcode level count table from MPRAcount pipeline
id_out <- args[3]
floc <- args[4]

`%notin%` <- Negate(`%in%`)

cell_reps <- read.delim(celltypes, header=F, stringsAsFactors=F)
dataCount <- read.delim(count_table, header=T, stringsAsFactors=F)

colnames(cell_reps) <- c("file_loc","replicate","celltype","material")
cell_reps <- unique(cell_reps[,2:3])

dataCond <- as.data.frame(cell_reps$celltype)
rownames(dataCond) <- cell_reps$replicate
colnames(dataCond) <- "condition"

write.table(dataCond, paste0(id_out, "_condition.txt"), quote=F, sep="\t", row.names=T, col.names=F)

dataCount <- dataCount[,colnames(dataCount) %notin% c("Error","CIGAR","MD","cs","Aln_Start.Stop")]

message(paste0(colnames(dataCount), collapse="\t"))

agg_rep_counts <- list()
agg_rep_bcs <- list()
indv_rep_bcs <- list()
indv_rep_counts <- list()
for(celltype in unique(cell_reps$celltype)){
	message(celltype)
	indv_rep_bcs[[celltype]] <- list()
	indv_rep_counts[[celltype]] <- list()
	reps <- rownames(dataCond)[which(dataCond$condition==celltype)]
	message(paste0(reps, collapse="\t"))
	cell_sub <- dataCount[,reps]
	agg_rep_bc <- data.frame(table(dataCount$Oligo[which(rowSums(cell_sub)>0)]))
	agg_gt10 <- nrow(agg_rep_bc[which(agg_rep_bc$Freq>10),])
	message("aggregating counts")
	agg_count <- aggregate(. ~Oligo, data=dataCount[,-1], FUN=sum)
	agg_count$means <- rowMeans(agg_count[,reps])
	mean_bound <- as.numeric(quantile(agg_count$means, seq(0,1,0.01))[91])
	tot_bound <- 0.5*max(agg_count$means, na.rm=T)
	if(celltype=="DNA"){
		xup <- max(mean_bound,tot_bound)
	}
	if(celltype!="DNA"){
		xup <- min(mean_bound,tot_bound)
	}
	message(mean_bound)
	message(tot_bound)
	agg_ct_gt20 <- nrow(agg_count[which(agg_count$means > 20),])
	total_oligos <- nrow(agg_count)

	agg_rep_bcs[[celltype]] <- ggplot(agg_rep_bc, aes(x=Freq)) + geom_histogram(bins=200) + geom_vline(xintercept=10, col="red") + xlab("Barcodes per aggregated Oligo") + ggtitle(paste0("Aggregated Barcode Count\n",agg_gt10," Oligos with > 10 Barcodes (",round(100*(agg_gt10/total_oligos),digits=1),")")) + theme_light()
	agg_rep_counts[[celltype]] <- ggplot(agg_count, aes(x=means)) + geom_histogram(bins=300) + geom_vline(xintercept=20, col="red") + xlim(0,xup) + xlab("Mean Count per  aggregated Oligo") + ggtitle(paste0("Mean Oligo Counts\n", agg_ct_gt20," Oligos with Mean Count > 20 (",round(100*(agg_ct_gt20/total_oligos),digits=1),")")) + theme_light()

	for(rep in reps){
		### Individual replicate barcode histograms
		message(rep)
		indv_hist <- data.frame(table(dataCount$Oligo[which(dataCount[,rep] > 0)]))
		indv_gt10 <- nrow(indv_hist[which(indv_hist$Freq>10),])
		indv_rep_bcs[[celltype]][[rep]] <- ggplot(indv_hist, aes(x=Freq)) + geom_histogram(bins=200) + geom_vline(xintercept=10, col="red") + xlab("Barcodes per Oligo") + ggtitle(paste0("Barcode Count ", rep, "\n", indv_gt10, " Oligos with > 10 Barcodes")) + theme_light()
		indv_agg <- aggregate(. ~Oligo, data=dataCount[,c("Oligo",rep)], FUN=sum)
		colnames(indv_agg) <- c("Oligo","counts")
		counts_bound <- as.numeric(quantile(indv_agg$counts, seq(0,1,0.01))[81])
		indv_ct_gt20 <- nrow(indv_agg[which(indv_agg$counts>20),])
		max_ct <- max(indv_agg$counts)
		indv_rep_counts[[celltype]][[rep]] <- ggplot(indv_agg, aes(x=counts)) + geom_histogram(bins=200) + geom_vline(xintercept=10, col="red") + xlim(0,counts_bound) + xlab("Counts per Oligo") + ggtitle(paste0("Oligo Count ", rep, "\n", indv_ct_gt20, " Oligos with Count > 20, max count = ", max_ct)) + theme_light()
	}
	message("plotting")
	grid.title <- paste0("Individual Replicate Barcode Histograms - ", celltype, "\nTotal Oligos: ", total_oligos)
	pdf(paste0(floc,"/",id_out,"_",celltype,"_indv_rep_barcode_QC.pdf"), width=10, height=10)
	#grid.arrange(grobs=list(indv_rep_bcs[[celltype]][1:length(reps)]), top=grid.title)
	do.call(grid.arrange, c(indv_rep_bcs[[celltype]], top=grid.title))
	dev.off()

	message("plotting")
	grid.count.title <- paste0("Individual Replicate Count Histograms - ", celltype, "\nTotal Oligos: ", total_oligos)
	pdf(paste0(floc,"/",id_out,"_",celltype,"_indv_rep_counts_QC.pdf"), width=10, height=10)
	#grid.arrange(grobs=list(indv_rep_bcs[[celltype]][1:length(reps)]), top=grid.title)
	do.call(grid.arrange, c(indv_rep_counts[[celltype]], top=grid.count.title))
	dev.off()

	agg.grid.title <- paste0("Aggregated Histograms - ", celltype, "\nTotal Oligos: ", total_oligos)
	pdf(paste0(floc,"/",id_out,"_",celltype,"_agg_count_QC.pdf"), width=10, height=10)
	grid.arrange(grobs=list(agg_rep_bcs[[celltype]], agg_rep_counts[[celltype]]), top=agg.grid.title)
	dev.off()
}
