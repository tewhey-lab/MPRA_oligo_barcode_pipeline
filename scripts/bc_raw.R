# Write a celltype specific version of the count table with plasmid and each celltype

options(stringsAsFactors = FALSE)
args = commandArgs(trailingOnly=TRUE)

cond.data <- args[1] ## Tab delimited file with celltype in first column and number of replicates in the second. Should include DNA/plasmid
count.table <- args[2] ## Barcode level count table from MPRAcount pipeline
id_out <- args[3]
floc <- args[4]

conditionStandard <- function(conditionData){
  cond_data <- as.data.frame(conditionData)
  colnames(cond_data)[1] <- "condition"
  cond_data[,1] <- factor(cond_data[,1])
  cond_data$condition <- relevel(cond_data$condition, "DNA")

  return(cond_data)
}

bcRawOut <- function(countsData, conditionData, file_prefix, floc){
  conditionData <- conditionStandard(conditionData)
  for(celltype in levels(conditionData$condition)){
    if(celltype=="DNA") next
    reps <- rownames(conditionData)[which(conditionData$condition=="DNA" | conditionData$condition==celltype)]
    count_temp <- countsData[,c("Barcode","Oligo","cs",reps)]
    write.table(count_temp, paste0(floc, "/", file_prefix,"_",celltype,".counts"), quote=F, sep = "\t", row.names=F)
  }
}

cond_data <- read.delim(cond.data, row.names=1, header=F, stringsAsFactors=F)
colnames(cond_data) <- "condition"

count_table <- read.delim(count.table, stringsAsFactors=F, header=T)

bcRawOut(count_table, cond_data, id_out, floc)
