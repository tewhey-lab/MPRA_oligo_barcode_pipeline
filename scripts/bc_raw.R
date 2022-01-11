# Write a celltype specific version of the count table with plasmid and each celltype

options(stringsAsFactors = FALSE)
args = commandArgs(trailingOnly=TRUE)

cond.data <- args[1] ## Tab delimited file with celltype in first column and number of replicates in the second. Should include DNA/plasmid
count.table <- args[2] ## Barcode level count table from MPRAcount pipeline
id_out <- args[3]
floc <- args[4]

bcRawOut <- function(countsData, conditionData, file_prefix, floc){
  conditionData <- conditionStandard(conditionData)
  for(celltype in levels(conditionData$condition)){
    if(celltype=="DNA") next
    reps <- rownames(conditionData)[which(conditionData$condition=="DNA" | conditionData$condition==celltype)]
    count_temp <- countsData[,c("Barcode","Oligo",reps)]
    write.table(count_temp, paste0(floc, "/", file_prefix,"_",celltype,".counts"), quote=F, sep = "\t")
  }
}

cond_data <- read.delim(cond.data, row.names=1, header=F, stringsAsFactors=F)
colnames(cond_data) <- "condition"

count_table <- read.delim(count.table, stringsAsFactors=F, header=T)

bcRawOut(count_table, cond_data, id_out, floc)
