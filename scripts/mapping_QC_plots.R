library(reshape2)
library(ggplot2)
library(gridExtra)

args = commandArgs(trailingOnly=TRUE)

parsed_file <- args[1]
id_out <- args[2]

dict <- read.delim(parsed_file, header=F,stringsAsFactors = F)
flags_dict <- dict[,c(5,7)]
colnames(flags_dict) <- c("flag_code","error_rate")
dict <- dict[,1:3]
colnames(dict) <- c("Barcode","Oligo","Seen")
# dict$Oligo <- as.factor(dict$Oligo)
dict$Seen <- as.numeric(dict$Seen)

libP_oligos <- colsplit(dict$Oligo,",",c("keep","reject"))
libP_reject <- dict[libP_oligos$keep=="*"|libP_oligos$reject!="",]
libP_keep <- dict[libP_oligos$keep!="*"&libP_oligos$reject=="",]

oligo_freq_P <- table(libP_keep$Oligo)
oligo_freq_P <- as.data.frame(oligo_freq_P)
colnames(oligo_freq_P) <- c("oligo","number_of_barcodes")

libP_keep_agg <- aggregate(Seen ~ Oligo, data = libP_keep, FUN=sum)
colnames(oligo_freq_P) <- c("Oligo","Barcodes")
count_hist_P <- merge(oligo_freq_P,libP_keep_agg,by="Oligo",all = T)

xlimP<-sum(quantile(count_hist_P$Seen,0.99))
meanP<-mean(count_hist_P$Seen)
maxoP<-max(count_hist_P$Seen)

flags_parsedP <- flags_dict[flags_dict$flag_code==0,]
flags_parsedP$error_rate <- as.numeric(flags_parsedP$error_rate)

flagP_ct<-data.frame(table(flags_dict$flag_code))
colnames(flagP_ct)<-c("Flag","Freq")
message(paste0(flagP_ct$Flag, collapse = "\t"))
message(paste0(flagP_ct$Freq, collapse = "\t"))
flagP_ct$Flag<-as.character(flagP_ct$Flag)
flagP_ct[flagP_ct$Flag==0,]$Flag<-"Passing"
flagP_ct[flagP_ct$Flag==2,]$Flag<-"Failed or No Mapping"
flagP_ct[flagP_ct$Flag==1,]$Flag<-"Conflict"

flagP_ct$percent <- (flagP_ct$Freq/sum(flagP_ct$Freq)) * 100
flagP_ct$row<-1



a <- ggplot(oligo_freq_P,aes(x=Barcodes)) + geom_histogram(bins=300) + labs(title=paste0("Number of Barcode Frequency - ", id_out)) + theme_light()
d <- ggplot(count_hist_P, aes(Seen)) +
  stat_ecdf(geom = "step") +
  theme(legend.position="top") +
  coord_cartesian(xlim=c(0,xlimP)) +
  geom_vline(xintercept=meanP, linetype="solid", color = "red", size=0.5) +
  geom_vline(xintercept=meanP*5, linetype="dashed", color = "red", size=0.5) +
  geom_vline(xintercept=meanP/5, linetype="dashed", color = "red", size=0.5) +
  xlab("Per Oligo Seq Coverage") + ggtitle(paste0("Oligo Cov CDF - truncated, max: ",maxoP,"\n",id_out)) +
  theme_bw() + theme(panel.grid.major = element_line(size = .25,colour = rgb(0,0,0,75,maxColorValue=255)), panel.grid.minor = element_blank()) 
f <- ggplot(flags_parsedP, aes(x=error_rate)) +
  geom_histogram()+
  theme(legend.position="top") +
  xlab("Error Rate for Passing Barcodes") + ggtitle(paste0("Oligo Error Rate - ", id_out)) +
  theme_bw() + theme(panel.grid.major = element_line(size = .25,colour = rgb(0,0,0,75,maxColorValue=255)), panel.grid.minor = element_blank()) 
h <- ggplot(flagP_ct, aes(x = row,y = Freq, fill = Flag)) +
  geom_bar(stat="identity") + 
  theme_bw() + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
  ggtitle("Sequence Mapping") +
  geom_text(aes(label = percent), position = position_stack(),vjust = 0.5)


pdf(paste0(id_out,"_barcode_qc.pdf"), width = 10, height = 10) # Open a new pdf file
grid.arrange(a,d,f,h)
dev.off()
