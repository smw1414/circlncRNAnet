#!/usr/bin/env Rscript
args=commandArgs(TRUE)
library("getopt")

spec <- matrix(c(
  'rawreads', 'r', 1, "character", "gene raw reads table (required) 
  \ngene\tSamples1  Sample2\nGene1\t23\t234\nGene2\t565\t23\n",
  'factorlist'     , 'f', 1, "character", "factor list, no header required, at least 3 reaplicates for each condition, place the reference condition to the top of list (required)
  \nSample1\tNormal <--- reference condition\nSample2\tNormal\nSample3\tTumor\nSample4\tTumor\n",
  'core'     , 'c', 2, "numeric","Number of cores for computation",
  'annotation'    , 'a', 2, "character", "-a gencodev25, -a gencodev19 (required)",
  'circmatrix'     , 'i', 2, "character","circRNA raw reads table",
  'mode'    , 'm', 2, "character", "lnc for lncRNA correlation, circ for circRNA correlation (required)"
  
),ncol=5,byrow=T)

opt = getopt(spec);

if ( is.null(opt$mode)) {
  cat(paste(getopt(spec, usage=T),"\nExample ./deseq2.r -r count_matrix.txt -f count_matrix_conditon.txt -m lnc -a gencodev25
            Example ./deseq2.r -r count_matrix_circ.txt -f factor_OSCC_Rseq_circrna.txt -m lnc -a gencodev19  
            Example ./deseq2.r -r OSCC_circRNA_raw2.txt -f factor_OSCC_Rseq_circrna.txt -m circ -a gencodev19
DEMO DATA CMD CIRCRNA ./deseq2.r -r encode_example_Gene_raw_read_count_casted.txt -i encode_example_circRNA_raw_read_count_casted.txt -f encode_example_circRNA_condition.txt -m circ -a gencodev25
DEMO DATA CMD LNCRNA ./deseq2.r -r TCGA_COADREAD_GENCODEV25_raw_read_count.txt -f TCGA_COADREAD_GENCODEV25_condition.txt -m lnc -a gencodev25 -c 6
DEMO DATA CMD TCGA  ./deseq2.r  -m TCGA-COAD \n"));
  q();
}

# define default value
if ( is.null(opt$core ) ) { opt$core = 4 }

library("data.table")
#library("DESeq2")
library("S4Vectors")
library("factoextra")
#library("BiocParallel")
#library("reshape2")

#read table and factor

ann<-opt$annotation[1]
mode<-opt$mode[1]

TCGA_folder<-""
TCGA_df<-data.frame(rds=system(paste0("ls ",TCGA_folder,"deseq_list*TCGA*rds"),intern = T))
TCGA_df$cancer<-gsub(".*(TCGA.*)\\.rds","\\1",TCGA_df$rds)

if (!(mode %like% "TCGA-")){
  if (ann=="gencodev25" ) {
    gencodev25_sym_id<-readRDS("gencodev25_sym_id.rds")
    symbol_id<-gencodev25_sym_id[["symbol_id"]]
    lncRNA_symbol_id<-gencodev25_sym_id[["lncRNA_symbol_id"]]
    coding_symbol_id<-gencodev25_sym_id[["coding_symbol_id"]]
  }
  if (ann=="gencodev19" ) {
    gencodev25_sym_id<-readRDS("gencodev19_sym_id.rds")
    symbol_id<-gencodev25_sym_id[["symbol_id"]]
    lncRNA_symbol_id<-gencodev25_sym_id[["lncRNA_symbol_id"]]
    coding_symbol_id<-gencodev25_sym_id[["coding_symbol_id"]]
  }
}


get_symbol_id<-function(gtf){
  v25lncRNAgtf<-data.table::fread(gtf,stringsAsFactors = F)
  lncRNA_symbol<-strsplit(unlist(unique(subset(v25lncRNAgtf,V3== "exon",select = V9))),";")
  symbol<-gsub(" gene_name \\\"","",sapply(lncRNA_symbol,"[[",5))
  id<-gsub("gene_id \\\"","",sapply(lncRNA_symbol,"[[",1))
  symbol<-(gsub("\\\"","",symbol))
  id<-(gsub("\\\"","",id))
  symbol_id<-unique(data.frame(symbol=symbol,id=id))
  return(symbol_id)
}

get_transcript_symbol_id<-function(gtf){
  v25lncRNAgtf<-data.table::fread(gtf,stringsAsFactors = F)
  lncRNA_symbol<-strsplit(unlist(unique(subset(v25lncRNAgtf,V3== "exon",select = V9))),";")
  symbol<-gsub(" gene_name \\\"","",sapply(lncRNA_symbol,"[[",5))
  id<-gsub(" transcript_id \\\"","",sapply(lncRNA_symbol,"[[",2))
  symbol<-(gsub("\\\"","",symbol))
  id<-(gsub("\\\"","",id))
  symbol_id<-unique(data.frame(symbol=symbol,id=id))
  return(symbol_id)
}

get_gene_transcript_symbol_id<-function(gtf){
  v25lncRNAgtf<-data.table::fread(gtf,stringsAsFactors = F)
  lncRNA_symbol<-strsplit(unlist(unique(subset(v25lncRNAgtf,V3== "exon",select = V9))),";")
  gene_id<-gsub("gene_id \\\"","",sapply(lncRNA_symbol,"[[",1))
  transcript_id<-gsub(" transcript_id \\\"","",sapply(lncRNA_symbol,"[[",2))
  gene_id<-as.character(gsub("\\\"","",gene_id))
  transcript_id<-as.character(gsub("\\\"","",transcript_id))
  symbol_id<-unique(data.frame(gene_id=gene_id,transcript_id=transcript_id))
  return(symbol_id)
}

get_symbol_id_coding<-function(gtf){
  v25lncRNAgtf<-data.table::fread(gtf,stringsAsFactors = F)
  lncRNA_symbol<-strsplit(unlist(unique(subset(v25lncRNAgtf,V3== "exon",select = V9))),";")
  #lncRNA_symbol<-lncRNA_symbol[grep(" gene_type \"protein_coding\"",lncRNA_symbol)]
  symbol<-gsub(" gene_name \\\"","",sapply(lncRNA_symbol,"[[",5))
  id<-gsub("gene_id \\\"","",sapply(lncRNA_symbol,"[[",1))
  type<-gsub("gene_id \\\"","",sapply(lncRNA_symbol,"[[",3))
  symbol<-(gsub("\\\"","",symbol))
  id<-(gsub("\\\"","",id))
  symbol_id<-unique(data.frame(symbol=symbol,id=id,type=type))
  symbol_id<-symbol_id[symbol_id$type %in% " gene_type \"protein_coding\"",]
  return(symbol_id[,1:2])
}

pca<-function(input,condition,output){
  
  gene_profile=data.table(input)
  cols<-colnames(gene_profile)[2:ncol(gene_profile)]
  gene_profile[,(cols):=lapply(.SD, function(x) {log2(as.numeric(x)+1)}), .SDcols = cols]
  gene_profile=dcast.data.table(melt(gene_profile, id.vars = "gene",variable.name = "sample"), sample ~ gene)
  sample_condition=condition
  colnames(sample_condition)=c("sample","group")
  gene_profile=merge(gene_profile,sample_condition,all.x=T,by="sample")
  gene_profile_pca= prcomp(gene_profile[,-c(1,ncol(gene_profile)),with=F],scale=TRUE)
  
  png("output/pca.png")
  p<-fviz_pca_ind(gene_profile_pca, geom = "point",
               habillage=gene_profile$group, addEllipses=TRUE,
               ellipse.level= 0.95)+scale_color_manual(values=c('#0C4B8E','#BF382A'))+scale_fill_manual(values=c("#3BBDF9","#F9713B"))+theme_minimal()
  print(p)
  dev.off()
  
  pca=as.data.frame(gene_profile_pca$x)
  pca$group=gene_profile$group
  pca$sample=gene_profile$sample
  write.table(pca,output,quote=F,row.names=F,sep="\t")
}


if (mode=="lnc") {
  rawreads <- as.data.frame(data.table::fread(opt$rawreads[1],stringsAsFactors = F))
  factor_list <- as.data.frame(data.table::fread(opt$factorlist[1], header=FALSE))
  factor_list$V2<-as.factor(factor_list$V2 )
  # table check
  if ( length(setdiff(colnames(rawreads[,2:ncol(rawreads)]),factor_list$V1))!=0 ){  
    system("echo Unmacthed sample names "); q();}
  
  if ( length(unique(factor_list$V2)) != 2 ){  
    system("echo Only support 2 levels factor "); q(); }
  
  if (  length(factor_list$V2[grep(levels(factor_list$V2)[1],factor_list$V2)]) <3 &   
        length(factor_list$V2[grep(levels(factor_list$V2)[2],factor_list$V2)]) <3  ) {  
    system("echo Not enought replicates "); q(); }

  symbol_id$id<-gsub("\\..*","",symbol_id$id)
  rawreads[,1]<-gsub("\\..*","",rawreads[,1])
  rawreads<-merge(symbol_id,rawreads,by.x="id",by.y=colnames(rawreads)[1])
  rawreads$id<-NULL
  
  rawreads_melt <- data.table::data.table(reshape2::melt(rawreads, id = "symbol"))
  rawreads<-data.table::dcast.data.table(rawreads_melt, symbol ~ variable, sum)
  
  countsTable <- as.data.frame(rawreads)
  names(countsTable)[1]<-"Geneid"
  rownames(countsTable)<-countsTable$Geneid
  countsTable$Geneid<-NULL
  countsTable<-round(as.matrix(countsTable))
  
  if (!is.numeric(countsTable)){ system("echo Data matrix contains non-numeric "); q(); }
  
  deseq_factor<-factor(factor_list[order(match(factor_list$V1,colnames(countsTable))), ]$V2)
  deseq_factor<-factor(deseq_factor, levels=c(levels(deseq_factor)[grep(as.character(factor_list$V2[1]),levels(deseq_factor))], # Tvs N or NvsT
                                              levels(deseq_factor)[-grep(as.character(factor_list$V2[1]),levels(deseq_factor))]))
  colData<-data.frame(row.names=colnames(countsTable),condition=deseq_factor)
  dds<-DESeq2::DESeqDataSetFromMatrix(countsTable,colData, ~condition)
  dds <- dds[ rowSums(DESeq2::counts(dds)) > 1, ]
  BiocParallel::register( BiocParallel::MulticoreParam(opt$core))
  minreplicates<-min(length(factor_list[factor_list$V2 %in% levels(factor_list$V2)[1],]$V1),
                     length(factor_list[factor_list$V2 %in% levels(factor_list$V2)[2],]$V1))
  dds<-DESeq2::DESeq(dds, parallel = T,minReplicatesForReplace= minreplicates)
  res <- DESeq2::results(dds, parallel = T) 
  res_tab<-as.data.frame(res)
  
  # condition mean reads
  system("echo condition mean reads")
  deseq_norm_reads<-as.data.frame(DESeq2::counts(dds,normalized=TRUE))
  level_1<-apply(deseq_norm_reads[,colnames(deseq_norm_reads) %in% 
                                    factor_list[factor_list$V2 %in% levels(factor_list$V2)[1],]$V1],1,mean)
  level_2<-apply(deseq_norm_reads[,colnames(deseq_norm_reads) %in% 
                                    factor_list[factor_list$V2 %in% levels(factor_list$V2)[2],]$V1],1,mean)
  mean_reads<-data.frame(level_1,level_2)
  colnames(mean_reads)<-c(paste0(levels(factor_list$V2)[1]," BaseMean"),paste0(levels(factor_list$V2)[2]," BaseMean"))
  res_tab<-merge(res_tab,mean_reads, by='row.names')[,c(1,3,6,7,2,8,9)]
  colnames(res_tab)[1]<-"gene"
  write.table(res_tab[res_tab$gene %in% lncRNA_symbol_id$symbol,], file="output/DEGs_lncRNA.txt",sep = '\t', row.names = F,quote = F)
  write.table(res_tab, file="output/DEGs.txt",sep = '\t', row.names = F,quote = F)
  deseq_norm_reads$gene<-rownames(deseq_norm_reads)
  # normalized all gene table
  deseq_norm_reads<-deseq_norm_reads[,c(ncol(deseq_norm_reads),1:ncol(deseq_norm_reads)-1)]
  write.table(deseq_norm_reads, file="output/norm_readstable.txt",sep = '\t', row.names = F,quote = F)
  # normalized lnc gene table
  deseq_norm_reads_lnc<-deseq_norm_reads[rownames(deseq_norm_reads) %in% lncRNA_symbol_id$symbol,]
  write.table(deseq_norm_reads_lnc, file="output/norm_readstable_lncRNA.txt",sep = '\t', row.names = F,quote = F)
  # coding gene
  deseq_norm_reads_coding<-deseq_norm_reads[rownames(deseq_norm_reads) %in% coding_symbol_id$symbol,]
  write.table(deseq_norm_reads_coding, file="output/norm_readstable_coding_gene.txt",sep = '\t', row.names = F,quote = F)
  
  pca(deseq_norm_reads,factor_list,"output/pca.txt")
}


if (mode=="circ") {
  factor_list <- as.data.frame(data.table::fread(opt$factorlist[1], header=FALSE))
  factor_list$V2<-as.factor(factor_list$V2 )
  
  rawreads <- as.data.frame(data.table::fread(opt$rawreads[1],stringsAsFactors = F))
  symbol_id$id<-gsub("\\..*","",symbol_id$id)
  rawreads[,1]<-gsub("\\..*","",rawreads[,1])
  rawreads<-merge(symbol_id,rawreads,by.x="id",by.y=colnames(rawreads)[1])
  rawreads$id<-NULL
  rawreads_melt <- data.table::data.table(reshape2::melt(rawreads, id = "symbol"))
  rawreads<-as.data.frame(data.table::dcast.data.table(rawreads_melt, symbol ~ variable, sum))
  rawreads<-reshape2::melt(rawreads)
  
  circmatrix <- as.data.frame(data.table::fread(opt$circmatrix[1],stringsAsFactors = F))
  circname<-circmatrix[,1] # load circRNA id from first col
  circmatrix<-reshape2::melt(circmatrix)
  colnames(circmatrix)[1]<-colnames(rawreads)[1]
  
  rawreads<-rbind(rawreads,circmatrix)
  
  rawreads<-reshape2::dcast(rawreads,symbol~variable,sum)
  
  # table check
  if ( length(setdiff(colnames(rawreads[,2:ncol(rawreads)]),factor_list$V1))!=0 ){  
    system("echo Unmacthed sample names "); q();}
  
  if ( length(unique(factor_list$V2)) != 2 ){  
    system("echo Only support 2 levels factor "); q(); }
  
  if (  length(factor_list$V2[grep(levels(factor_list$V2)[1],factor_list$V2)]) <3 &   
        length(factor_list$V2[grep(levels(factor_list$V2)[2],factor_list$V2)]) <3  ) {  
    system("echo Not enought replicates "); q(); }
  
  countsTable <- as.data.frame(rawreads)
  names(countsTable)[1]<-"Geneid"
  rownames(countsTable)<-countsTable$Geneid
  countsTable$Geneid<-NULL
  countsTable<-round(as.matrix(countsTable))
  
  if (!is.numeric(countsTable)){ system("echo Data matrix contains non-numeric "); q(); }
  
  deseq_factor<-factor(factor_list[order(match(factor_list$V1,colnames(countsTable))), ]$V2)
  deseq_factor<-factor(deseq_factor, levels=c(levels(deseq_factor)[grep(as.character(factor_list$V2[1]),levels(deseq_factor))], # Tvs N or NvsT
                                              levels(deseq_factor)[-grep(as.character(factor_list$V2[1]),levels(deseq_factor))]))
  colData<-data.frame(row.names=colnames(countsTable),condition=deseq_factor)
  dds<-DESeq2::DESeqDataSetFromMatrix(countsTable,colData, ~condition)
  dds <- dds[ rowSums(DESeq2::counts(dds)) > 1, ]
  BiocParallel::register(BiocParallel::MulticoreParam(opt$core))
  minreplicates<-min(length(factor_list[factor_list$V2 %in% levels(factor_list$V2)[1],]$V1),
                     length(factor_list[factor_list$V2 %in% levels(factor_list$V2)[2],]$V1))
  dds<-DESeq2::DESeq(dds, parallel = T,minReplicatesForReplace= minreplicates)
  res <- DESeq2::results(dds, parallel = T) 
  res_tab<-as.data.frame(res)
  
  
  # condition mean reads
  system("echo condition mean reads")
  deseq_norm_reads<-as.data.frame(DESeq2::counts(dds,normalized=TRUE))
  level_1<-apply(deseq_norm_reads[,colnames(deseq_norm_reads) %in% 
                                    factor_list[factor_list$V2 %in% levels(factor_list$V2)[1],]$V1],1,mean)
  level_2<-apply(deseq_norm_reads[,colnames(deseq_norm_reads) %in% 
                                    factor_list[factor_list$V2 %in% levels(factor_list$V2)[2],]$V1],1,mean)
  mean_reads<-data.frame(level_1,level_2)
  colnames(mean_reads)<-c(paste0(levels(factor_list$V2)[1]," BaseMean"),paste0(levels(factor_list$V2)[2]," BaseMean"))
  res_tab<-merge(res_tab,mean_reads, by='row.names')[,c(1,3,6,7,2,8,9)]
  colnames(res_tab)[1]<-"gene"
  #write.table(res_tab[res_tab$gene %in% lncRNA_symbol_id$symbol,], file="DEGs_circRNA.txt",sep = '\t', row.names = F,quote = F)
  write.table(res_tab[ res_tab$gene %in% circname,], file="output/DEGs_circRNA.txt",sep = '\t', row.names = F,quote = F)
  deseq_norm_reads$gene<-rownames(deseq_norm_reads)
  deseq_norm_reads<-deseq_norm_reads[,c(ncol(deseq_norm_reads),1:ncol(deseq_norm_reads)-1)]
  write.table(deseq_norm_reads[deseq_norm_reads$gene %in% circname,], file="output/norm_readstable_circRNA.txt",sep = '\t', row.names = F,quote = F)
  write.table(deseq_norm_reads[!deseq_norm_reads$gene %in% circname,], file="output/norm_readstable.txt",sep = '\t', row.names = F,quote = F)

  # normalized lnc gene table
  deseq_norm_reads_lnc<-deseq_norm_reads[rownames(deseq_norm_reads) %in% lncRNA_symbol_id$symbol,]
  write.table(deseq_norm_reads_lnc, file="output/norm_readstable_lncRNA.txt",sep = '\t', row.names = F,quote = F)
  # coding gene
  deseq_norm_reads_coding<-deseq_norm_reads[rownames(deseq_norm_reads) %in% coding_symbol_id$symbol,]
  write.table(deseq_norm_reads_coding, file="output/norm_readstable_coding_gene.txt",sep = '\t', row.names = F,quote = F)
  
  pca(deseq_norm_reads[deseq_norm_reads$gene %in% circname,],factor_list,"output/pca.txt")
  }



if (mode %like% "TCGA-") {
  message("opening TCGA deseq2 result")
  tcga.rds.path<-as.character(TCGA_df[ TCGA_df$cancer == mode,]$rds)
  tcga.deseq.res<-readRDS(tcga.rds.path)
  message("writing TCGA deseq2 result")
  factor_list<-data.frame(V1=names(tcga.deseq.res[[ "norm_readstable" ]][2:ncol(tcga.deseq.res[[ "norm_readstable" ]])]))
  factor_list$V2<-ifelse(substr(factor_list$V1,14,14) >=1,"N","T")
  write.table(factor_list, file="count_matrix_condition.txt",sep = '\t', row.names = F,quote = F,col.names = F)
  write.table(tcga.deseq.res[[ "DEGs_lncRNA" ]], file="output/DEGs_lncRNA.txt",sep = '\t', row.names = F,quote = F)
  write.table(tcga.deseq.res[[ "DEGs" ]], file="output/DEGs.txt",sep = '\t', row.names = F,quote = F)
  write.table(tcga.deseq.res[[ "norm_readstable" ]], file="output/norm_readstable.txt",sep = '\t', row.names = F,quote = F)
  write.table(tcga.deseq.res[[ "norm_readstable_lncRNA" ]], file="output/norm_readstable_lncRNA.txt",sep = '\t', row.names = F,quote = F)
  write.table(tcga.deseq.res[[ "norm_readstable_coding_gene" ]], file="output/norm_readstable_coding_gene.txt",sep = '\t', row.names = F,quote = F)
  message("processing PCA")
  pca(tcga.deseq.res[[ "norm_readstable" ]],factor_list,"output/pca.txt")
}







