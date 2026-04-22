args = commandArgs(trailingOnly=TRUE)

# test if there are all the necessary arguments
if (length(args)<6) {
  stop("at least 6 arguments must be supplied: counts_file_path, metadata_file_path, up_value, outdir_path, targeted_gene_id, confounders flag", call.=FALSE)
}

require(dplyr)
require(stringr)
require(data.table)
require(DESeq2)
require(apeglm)

count_file_path = args[1]
metadata_file_path = args[2]
up_value = args[3]
outdir_path = args[4]
targeted_gene_id = args[5]
confounders = args[6]
lfcThresholdVal = args[7] # this allows to specify custom value for null hypothesis of Wald Test
if (is.na(lfcThresholdVal)) {
  lfcThresholdVal = 0 # by default we test log2FC against zero
}
apeglm = args[8] # wether to use apeglm shrinkage of log2FC values
if (is.na(apeglm)) {
  apeglm = 'FALSE' # by default we test log2FC against zero
}

# helper function to parse input paths
trim_quotes_from_path<-function(path_string){
path_string = gsub("'", '', path_string)
path_string = gsub('"', '', path_string)
return(path_string)
}

# loading count matrix
cts <- read.csv(trim_quotes_from_path(count_file_path), sep='\t',header = TRUE,stringsAsFactors=FALSE,check.names = FALSE)
id_column <- colnames(cts)[1]
cts <- cts[!duplicated(cts[ , c(id_column)]),]

# loading metadata
metadata <- read.csv(trim_quotes_from_path(metadata_file_path), sep='\t',header = TRUE,stringsAsFactors=FALSE,check.names = FALSE)
colnames(metadata) <- c('sample','condition','confounder_1','confounder_2','confounder_3')
rownames(metadata) <- metadata$sample
metadata$condition <- factor(metadata$condition==up_value)

# prepare count matrix for DESEQ format
rownames(cts) <- cts[,id_column]
cts <- select(cts, rownames(metadata))

# create a DESeqDataSet
if(confounders=='three_confounders'){
    dds <- DESeqDataSetFromMatrix(countData = cts, colData = metadata, design = ~condition+confounder_1+confounder_2+confounder_3)
}
if(confounders=='two_confounders'){
    dds <- DESeqDataSetFromMatrix(countData = cts, colData = metadata, design = ~condition+confounder_1+confounder_2)
}
if(confounders=='one_confounder'){
    dds <- DESeqDataSetFromMatrix(countData = cts, colData = metadata, design = ~condition+confounder_1)
}
if(confounders=='no_confounders'){
    dds <- DESeqDataSetFromMatrix(countData = cts, colData = metadata, design = ~condition)
}

# VST COUNT TRANSFORMATION
# vsd <- vst(dds, blind = FALSE)
# vsd_table = assay(vsd)
# vsd_table = data.frame(vsd_table)
# cols <- colnames(vsd_table)
# vsd_table <- cbind(rownames(vsd_table), vsd_table)
# colnames(vsd_table) <- c(id_column,cols)
# # write the VST transformed counts
# write.table(vsd_table,file=paste(trim_quotes_from_path(outdir_path),"/VST_counts.tsv",sep = ""),row.names = FALSE, quote=FALSE, sep='\t')

# do not run Deseq2 if there is only one replicate per condition
if(nrow(metadata)==2){
    res <- vsd_table
    a = metadata[metadata$condition==TRUE,]['index'][[1]]
    b = metadata[metadata$condition==FALSE,]['index'][[1]]
    res$log2FoldChange <- log2(res[,a])-log2(res[,b])
    res$pvalue <- 1
    res$padj <- 1
    res$gene_name = rownames(res)
    res <- res[,c('gene_name','log2FoldChange','pvalue','padj')] # To be changed to stick to same output format
}else{
    dds <- DESeq(dds)
    print(resultsNames(dds))
    if(apeglm=='TRUE'){
        res <- lfcShrink(dds, coef="condition_TRUE_vs_FALSE", type="apeglm", lfcThreshold=lfcThresholdVal)
    }else{
        res <- results(dds, contrast=c("condition",TRUE,FALSE), lfcThreshold=lfcThresholdVal) # two-sided p-value    
    }
    # res <- lfcShrink(dds, coef="condition_TRUE_vs_FALSE", type="apeglm", lfcThreshold=lfcThresholdVal, altHypothesis="greaterAbs")
    # res <- lfcShrink(dds, coef="conditionTRUE", apeAdapt=FALSE,type="apeglm", lfcThreshold=lfcThresholdVal, altHypothesis="greaterAbs")
    # res <- lfcShrink(dds, coef="conditionTRUE", type="normal", lfcThreshold=lfcThresholdVal, altHypothesis="greaterAbs")
    res = data.frame(res[order(res$padj),])
    res$gene_name = rownames(res)
    res <- res[,c('gene_name','baseMean','log2FoldChange','lfcSE','pvalue','padj')]
}

# write the log2FCs
write.table(res,file=paste(trim_quotes_from_path(outdir_path),"/DE_table.tsv",sep = ""),row.names = FALSE, quote=FALSE, sep='\t')
