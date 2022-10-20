library(ggplot2)
options(stringsAsFactors = FALSE)
args = commandArgs(trailingOnly=TRUE)

stats_out_file <- args[1]
acc_file <- args[2]
id_out <- args[3]
out_dir <- args[4]

stats_out <- read.delim(stats_out_file, stringsAsFactors = F, header=F, sep = "\t")
acc <- read.delim(acc_file, stringsAsFactors = F, header=F)
colnames(stats_out) <- c("rep","proj","barcodes","good_reads","total_reads","per_good")

stats_out$cell <- NA
for(rep in unique(stats_out$rep)){
  message(rep)
  stats_out$cell[which(stats_out$rep==rep)] <- acc$V3[which(acc$V2==rep)]
}

a <- ggplot(stats_out, aes(x=rep, y=barcodes, fill=cell)) + geom_bar(stat="identity") + coord_flip() + theme_light()
b <- ggplot(stats_out, aes(x=rep, fill=cell)) + geom_bar(stat="identity", aes(y=total_reads), alpha=0.5) + geom_bar(stat="identity", aes(y=good_reads)) + coord_flip() + theme_light()


pdf(paste0(out_dir,"/",id_out,"_read_stats.pdf"))
print(a)
print(b)
dev.off()