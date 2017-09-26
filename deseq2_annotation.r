#!/usr/bin/env Rscript

##########################
# gene annotation module #
##########################
args=commandArgs(TRUE)
library("getopt")


spec <- matrix(c(
  'input', 'i', 1, "character", "output from Deseq2 (required)",
  'output','o', 2, "character", "output file (optional)"
  
),ncol=5,byrow=T)

opt = getopt(spec)

if ( is.null(opt$input)) {
  cat(paste(getopt(spec, usage=T)))
  q();
}

# define default value
if ( is.null(opt$output ) ) { opt$output = "output.txt" }

library("data.table")

# loading annotation
annotation=readRDS("dbNSFP3.2_gene_lite.rds")
deg=fread(opt$input,header=T,sep="\t")
genelist=paste(deg[abs(log2FoldChange)>=1&pvalue<=0.05,]$gene,collapse=",")
write.table(genelist,"output/genelist.txt",quote=F,col.names=F,row.names=F)
output=merge(deg,annotation,all.x=T,by.x="gene",by.y="Gene_name")

# merge gene with annotation
write.table(output,opt$output,quote=F,sep="\t",row.names=F)
#system(paste0("/share/apps/lancer2/tabulate_lncRNA_RPB_r1.R -q ",deg[abs(log2FoldChange)>=1&pvalue<=0.05,]$gene[1]))
#system(paste0("/share/apps/lancer2/correlation_20161206.r"," -n output/norm_readstable.txt -f count_matrix_conditon.txt -a gencodev25 -s symbol -m lnc -q ",genelist))
