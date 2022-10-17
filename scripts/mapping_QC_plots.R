library(reshape2)
library(ggplot2)
library(gridExtra)

options(stringsAsFactors = FALSE)
args = commandArgs(trailingOnly=TRUE)

parsed_file <- args[1]
hist_file <- args[2]
preseq_out <- args[3]
preseq_in <- args[4]
fasta_file<-args[5]
id_out <- args[6]


count.hist<-read.delim(hist_file,header=FALSE)
flags<-read.delim(parsed_file,header=FALSE)[,c(5,7)]
fasta<-read.table(fasta_file,header=FALSE)
preseq.out <- read.delim(preseq_out, header = T)
preseq.in <- read.delim(preseq_in, header = F)

message(paste0(count.hist[1,], collapse = "\t"))
message(paste0(preseq.out[1,], collapse = "\t"))
message(paste0(preseq.in[1,], collapse = "\t"))

message(max(flags$V7))
message(min(flags$V7))

fasta<-data.frame("ID"=fasta[seq(1,nrow(fasta),2),],"seq"=fasta[seq(2,nrow(fasta),2),])
fasta[,"ID"]<-gsub("^>","",fasta[,"ID"])

parsed<-count.hist[which(count.hist$V2<quantile(count.hist$V2,.99)),]
maxb<-max(count.hist$V2)

plotA<-ggplot(parsed, aes(x=V2)) +
  geom_histogram(bins=200)+
  theme(legend.position="top") +
  xlab("Barcodes per Oligo") + ggtitle(paste0("Barcode Count - truncated, max: ",maxb)) +
  theme_bw() + theme(panel.grid.major = element_line(size = .25,colour = rgb(0,0,0,75,maxColorValue=255)), panel.grid.minor = element_blank()) 

xlim<-sum(quantile(count.hist$V3,0.99))
mean<-mean(count.hist$V3)
maxo<-max(count.hist$V3)

plotB<-ggplot(count.hist, aes(V3)) +
  stat_ecdf(geom = "step") +
  theme(legend.position="top") +
  coord_cartesian(xlim=c(0,xlim)) +
  geom_vline(xintercept=mean, linetype="solid", color = "red", size=0.5) +
  geom_vline(xintercept=mean*5, linetype="dashed", color = "red", size=0.5) +
  geom_vline(xintercept=mean/5, linetype="dashed", color = "red", size=0.5) +
  xlab("Per Oligo Seq Coverage") + ggtitle(paste0("Oligo Cov CDF - truncated, max: ",maxo)) +
  theme_bw() + theme(panel.grid.major = element_line(size = .25,colour = rgb(0,0,0,75,maxColorValue=255)), panel.grid.minor = element_blank()) 




parsed<-flags[flags$V5==0,]
parsed$V7<-as.numeric(parsed$V7)
plotC<-ggplot(parsed, aes(x=V7)) +
  geom_histogram()+
  theme(legend.position="top") +
  xlab("Error Rate for Passing Barcodes") + ggtitle("Oligo Error Rate") +
  theme_bw() + theme(panel.grid.major = element_line(size = .25,colour = rgb(0,0,0,75,maxColorValue=255)), panel.grid.minor = element_blank()) 


flag.ct<-data.frame(table(flags$V5))
colnames(flag.ct)<-c("Flag","Freq")
flag.ct$Flag<-as.character(flag.ct$Flag)
flag.ct[flag.ct$Flag==0,]$Flag<-"Passing"
flag.ct[flag.ct$Flag==2,]$Flag<-"Failed or No Mapping"
flag.ct[flag.ct$Flag==1,]$Flag<-"Conflict"

flag.ct$percent <- (flag.ct$Freq/sum(flag.ct$Freq)) * 100
flag.ct$row<-1

plotD<-ggplot(flag.ct, aes(x = row,y = Freq, fill = Flag)) +
  geom_bar(stat="identity") + 
  theme_bw() + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
  ggtitle("Sequence Mapping") +
  geom_text(aes(label = percent), position = position_stack(vjust = 0.5))


### Here total distinct actually needs to be what is total found and total found needs to be the sum of product of the two columns
plotE <- ggplot(preseq.out, aes(x=TOTAL_READS, y=EXPECTED_DISTINCT)) + geom_point() + geom_line() + geom_ribbon(aes(ymin=LOWER_0.95CI,ymax=UPPER_0.95CI, fill='prediction'), fill="black", alpha=0.4) + 
  geom_hline(yintercept = sum(preseq.in$V2), col="red") + geom_vline(xintercept = sum(preseq.in$V2*preseq.in$V1), col="red") +
  xlab("Total Reads") + ylab("Expected Distinct Reads") + 
  ggtitle("Expected Dinstinct per Total Reads", subtitle = paste0("total found: ", sum(preseq.in$V2*preseq.in$V1)," total distinct: ", sum(preseq.in$V2))) +
  theme_light()

seen<-nrow(count.hist)
total<-nrow(fasta)
per<-round(seen/total,4)*100

grid.title<-paste0(id_out," - ",per,"% captured - ",seen,"/",total)

pdf(paste0(id_out,"_barcode_qc.pdf"), width = 10, height = 10) # Open a new pdf file
grid.arrange(plotA,plotB,plotC,plotD, plotE,top=grid.title)
dev.off()