#!/usr/bin/env Rscript
#####!/share/apps/R/bin/Rscript
args=commandArgs(TRUE)
library("getopt")

spec <- matrix(c(
  'input', 'i', 1, "character", "output from correlation.r (required)",
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

annotation=readRDS("dbNSFP3.2_gene_lite.rds")
coexp=fread(opt$input,header=T,sep="\t")
output=merge(coexp,annotation,all.x=T,by.x="co_exp_gene",by.y="Gene_name")
write.table(output,opt$output,quote=F,sep="\t",row.names=F)
