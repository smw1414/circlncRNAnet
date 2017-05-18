#!/usr/bin/env Rscript
args=commandArgs(TRUE)

library(getopt)

# 1=required argument;2=optional argument
spec <- matrix(c(
  'normreads', 'n', 2, "character", "normalized reads table (required)",
  'factorlist'     , 'f', 2, "character", "factor list (required)",
  'circRNA'     , 'c', 2, "character", "circRNA normalized reas table (optional)",
  'query'     , 'q', 1, "character", "Query genes, separate by comma (required)",
  'mode'    , 'm', 2, "character", "lnc for lncRNA correlation, circ for circRNA correlation (required)",
  'annotation'    , 'a', 2, "character", "-a gencodev25, -a gencodev19 (required)",
  'symbolid'    , 's', 2, "character", "use gene symbol or ensembl id. -s symbol, -s id ",
  'demo'    , 'd', 2, "character", "-d COADREAD "
),ncol=5,byrow=T)

opt = getopt(spec);

if (is.null(opt$query)) {
  cat(paste(getopt(spec, usage=T),"\n","Example:  ./correlation.r -n CRC_Rseq_partekflow_rpkm_norm.txt -f factor_CRC_Rseq_partekflow_norm.txt -q FAM83H-AS1,ELFN1-AS1,LINC00346,SNHG15,AC021218.2 -a gencodev25 -s symbol -m lnc
 Example:  ./correlation.r -d COADREAD -q FAM83H-AS1,ELFN1-AS1,LINC00346,SNHG15,AC021218.2,DDX11-AS1
 Example:  ./correlation.r -n CRC_Rseq_hg38_id_deseq.txt -f factor_ID.txt -q ENSG00000203499,ENSG00000236081,ENSG00000255874,ENSG00000232956,ENSG00000204876 -a gencodev25 -s id -m lnc 
 Example:  ./correlation.r -n OSCC_Rseq.txt -f factor_OSCC_Rseq_circrna.txt -q chrX_47755339_47705503_fwd,chr2_191537878_191523883_fwd -a gencodev19 -s symbol -m circ -c OSCC_circRNA2.txt
DEMO DATA CMD CIRCRNA ./correlation.r -n output/norm_readstable.txt -f encode_example_circRNA_condition.txt -q chr11_35204640_35201082_fwd,chr10_97437191_97438703_rev,chr9_128515639_128508876_fwd,chr12_68836749_68828771_fwd -a gencodev25 -s symbol -m circ -c norm_readstable_circRNA.txt
DEMO DATA CMD LNCRNA ./correlation.r -n output/norm_readstable.txt -f TCGA_COADREAD_GENCODEV25_condition.txt -q CCAT1,ELFN1-AS1,LINC00346,SNHG15,AC021218.2 -a gencodev25 -s symbol -m lnc
DEMO DATA CMD TCGA ./correlation.r -n output/norm_readstable.txt -f count_matrix_condition.txt -q CCAT1,ELFN1-AS1,LINC00346,SNHG15,AC021218.2 -a gencodev25 -s symbol -m lnc\n"));
  q();
}

#if (is.null(opt$circRNA) & opt$mode == "circ") {
##  cat(paste("please attach a circRNA table\n"))
#  q();
#}
run_time_message<-function(msg){
  message(paste0(Sys.time()," : ",msg))
}
# #idlnc
# idlnc<-function(){
#   opt$query = "ENSG00000203499,ENSG00000236081,ENSG00000255874,ENSG00000232956,ENSG00000204876"
#   opt$factorlist = "factor_ID.txt"
#   opt$normreads= "CRC_Rseq_hg38_id_deseq.txt"
#   opt$annotation = "gencodev25"
#   opt$symbolid = "id"
#   opt$mode = "lnc"
# }
# symbollnc<-function(){
#  # opt$query = "MAGEA6,KLK6,MAGEA3,SLC26A9,NOTUM,TCN1,MMP7,COL10A1,PAH,ONECUT3,CST1,KRT23,OTOP2,WNT2,IGF2BP1,SIM2,HABP2,INHBA,ESM1,FOXQ1"
#   opt$query = "FAM83H-AS1,ELFN1-AS1,LINC00346,SNHG15,AC021218.2,CCAT1"
#   opt$factorlist = "factor_CRC_Rseq_partekflow_norm.txt"
#   opt$normreads= "CRC_Rseq_partekflow_norm.txt"
#   opt$annotation = "gencodev25"
#   opt$symbolid = "symbol"
#   opt$mode = "lnc"
# }
# 
# symboldemo<-function(){
#   # opt$query = "MAGEA6,KLK6,MAGEA3,SLC26A9,NOTUM,TCN1,MMP7,COL10A1,PAH,ONECUT3,CST1,KRT23,OTOP2,WNT2,IGF2BP1,SIM2,HABP2,INHBA,ESM1,FOXQ1"
#   opt$query = "FAM83H-AS1,ELFN1-AS1,LINC00346,SNHG15,AC021218.2,DDX11-AS1"
#   opt$factorlist = "NO"
#   opt$normreads= "NO"
#   opt$annotation = "gencodev25"
#   opt$symbolid = "symbol"
#   opt$mode = "lnc"
#   opt$demo="COADREAD"
# }
# 
# symcirc<-function(){
#  
#   opt$query = "chr11_35204640_35201082_fwd,chr10_97437191_97438703_rev,chr9_128515639_128508876_fwd,chr12_68836749_68828771_fwd"
#   opt$factorlist = "encode_example_circRNA_condition.txt"
#   opt$normreads= "norm_readstable.txt"
#   opt$annotation = "gencodev25"
#   opt$symbolid = "symbol"
#   opt$circRNA = "norm_readstable_circRNA.txt"
#   opt$mode = "circ"
# }
# 
# 
# symcirc<-function(){
#   opt$query = "chrX_47755339_47705503_fwd,chr2_191537878_191523883_fwd"
#   opt$factorlist = "factor_OSCC_Rseq_circrna.txt"
#   opt$normreads= "OSCC_Rseq.txt"
#   opt$annotation = "gencodev19"
#   opt$symbolid = "symbol"
#   opt$circRNA = "OSCC_circRNA2.txt"
#   opt$mode = "circ"
# }
# 
# 
 symcirc<-function(){
   #OSCC all
  opt$query = "chr11_102871594_102871952_rev,chr5_38523419_38530666_rev,chr20_17955365_17957037_rev,chr3_58081776_58077046_fwd,chr10_115129535_115120185_fwd,chr9_16435555_16437524_rev,chr17_69128619_69134742_rev"
  opt$factorlist = "OSCC/oscc_allsample_condition.txt"
  opt$normreads= "norm_readstable.txt"
  opt$annotation = "gencodev25"
  opt$symbolid = "symbol"
  opt$circRNA = "norm_readstable_circRNA.txt"
  opt$mode = "circ"
}
symbollnc<-function(){
 # opt$query = "MAGEA6,KLK6,MAGEA3,SLC26A9,NOTUM,TCN1,MMP7,COL10A1,PAH,ONECUT3,CST1,KRT23,OTOP2,WNT2,IGF2BP1,SIM2,HABP2,INHBA,ESM1,FOXQ1"
  opt$query = "C5orf66,CDR1-AS,DSG1-AS1,HOXA10-AS,RP11-1223D19.3,RP11-132A1.4,RP11-146D12.2,RP11-161M6.2,RP11-681B3.4,RP11-774O3.3"
  opt$factorlist = "OSCC/oscc_allsample_condition.txt"
  opt$normreads= "norm_readstable.txt"
  opt$annotation = "gencodev25"
  opt$symbolid = "symbol"
  opt$mode = "lnc"
}

symbollnc<-function(){
  # opt$query = "MAGEA6,KLK6,MAGEA3,SLC26A9,NOTUM,TCN1,MMP7,COL10A1,PAH,ONECUT3,CST1,KRT23,OTOP2,WNT2,IGF2BP1,SIM2,HABP2,INHBA,ESM1,FOXQ1"
  opt$query = "C5orf66,CDR1-AS,DSG1-AS1,HOXA10-AS,RP11-1223D19.3,RP11-132A1.4,RP11-146D12.2,RP11-161M6.2,RP11-681B3.4,RP11-774O3.3"
  opt$factorlist = "output/TCGA_factor_list.txt"
  opt$normreads= "output/norm_readstable.txt"
  opt$annotation = "gencodev25"
  opt$symbolid = "symbol"
  opt$mode = "lnc"
}


#library("GenomicFeatures")
#library("biomaRt")
library("data.table")
library("ggplot2")
library("feather")
#library("WGCNA")
#system("export ALLOW_WGCNA_THREADS=6")
#allowWGCNAThreads()
#### reads table ########
run_time_message("reads table")
if(length(opt$demo)>0) {
  if (opt$demo[1]=="COADREAD"){
    
    normreads <- readRDS("COADREAD_fpkm_demo.rds")[[1]]
    factor_list <- readRDS("COADREAD_fpkm_demo.rds")[[2]]
    colnames(factor_list)<-c("V1","V2")
    normreads<-as.data.frame(normreads)
    row.names(normreads)<-as.character(normreads[,1])
    normreads[,1]<-NULL
    # opt$query = "FAM83H-AS1,ELFN1-AS1,LINC00346,SNHG15,AC021218.2"
    opt$factorlist = "readRDS(\"COADREAD_fpkm_demo.rds\")[[2]]"
    opt$normreads= "readRDS(\"COADREAD_fpkm_demo.rds\")[[1]]"
    opt$annotation = "gencodev25"
    opt$symbolid = "symbol"
    opt$mode = "lnc"
  } 
} else {
  normreads <- as.data.frame(fread(opt$normreads[1],stringsAsFactors = F))
  factor_list <- read.delim(opt$factorlist[1], header=FALSE)
  if (!is.null(opt$circRNA[1]) & opt$mode[1]=="circ" ){
    circrna<-as.data.frame(fread(opt$circRNA[1],stringsAsFactors = F))
    circrna<-as.data.frame(circrna)
    row.names(circrna)<-as.character(circrna[,1])
    circrna[,1]<-NULL
  }
  normreads<-as.data.frame(normreads)
  row.names(normreads)<-as.character(normreads[,1])
  normreads[,1]<-NULL
  
}


#
#system("rm -rf output ")
#system("mkdir -p output ")

###### loading db files #####
run_time_message("loading db files")
 if (opt$annotation[1]=="gencodev19"){
  chrsize <- read.delim("hg19.chrom.sizes", header=FALSE)
  lnc_gene <-readRDS("v19_lncrna_gene.rds")
  codeinggenelist<-readRDS("v19_protein_coding.rds")
  lncpedia<-read.delim("lncipedia_4_hc_hg19_id.txt", header=FALSE)
  exonidx<-"hg19.gencode.v19.annotation.gtf.exon.list.with.count.idx.gz"
  circdb<-readRDS("cirrnadb_hg19.rds")
  gene_coordinate<-readRDS("v19_gene_coordinate.rds")

} else if (opt$annotation[1]=="gencodev25") {
  chrsize <- read.delim("hg38.chrom.sizes", header=FALSE)
  lnc_gene <- readRDS("v25_lncrna_gene.rds")
  codeinggenelist<-readRDS("v25_protein_coding.rds")
  all_gene<-rbind(lnc_gene,codeinggenelist)
  lncpedia<-read.delim("lncipedia_4_hc_hg38_id.txt", header=FALSE)
  exonidx<-"hg38.gencode.v25.annotation.gtf.exon.list.with.count.idx.gz"
  circdb<-readRDS("cirrnadb_hg38.rds")
  gene_coordinate<-readRDS("v25_gene_coordinate.rds")
  
  
} else { cat(paste(getopt(spec, usage=T)));q(); }

rownames(chrsize)<-chrsize$V1
chrsize$V1<-NULL
chrsize<-t(chrsize)
names(chrsize)<-colnames(chrsize)
chrsize<-chrsize[1,]  
gene_coordinate<-unique(gene_coordinate[,c(1:8)]) # remove entrzid  #add 

df_opt<-as.data.frame(do.call("rbind", opt),stringsAsFactors = F)
saveRDS(df_opt,'output/opt.rds')
###### loading db files END ##### 



for_lnc_mode<-function(){
  run_time_message("lncRNA analysis start")
  ##### start table check and queries ##########
  if ( sum(length(setdiff(colnames(normreads [,1:ncol(normreads )]),factor_list$V1)),
           length(setdiff(factor_list$V1,colnames(normreads [,1:ncol(normreads )]))))!=0 ){  
    system("echo Unmacthed sample names > output/errormsg.txt"); q();}
  
  if ( length(unique(factor_list$V2)) != 2 ){  
    system("echo Only support 2 levels factor > output/errormsg.txt"); q(); }
  
  if (  length(factor_list$V2[grep(levels(factor_list$V2)[1],factor_list$V2)]) <3 |   
        length(factor_list$V2[grep(levels(factor_list$V2)[2],factor_list$V2)]) <3  ) { system("echo Not enought replicates > output/errormsg.txt"); q(); }
  #call known lncRNA id or names
  if (opt$symbolid[1] == "symbol") {geneidcolumn =all_gene$gencode_gene_symbol }
  if (opt$symbolid[1] == "id") {geneidcolumn  = all_gene$ensembl_gene_id }
  
  # add functin to check if its really id or gene name
  
  genequery<-unique(unlist(strsplit(opt$query[1],",")[1]))
  
  
  if (opt$symbolid[1] == "id"){  rownames(normreads) <-gsub("\\..*","",rownames(normreads)) }
  if (opt$symbolid[1] == "id"){  genequery <-gsub("\\..*","",genequery) }
  
  queriesindb<-genequery[genequery %in% geneidcolumn] 
  queriesnotindb <- setdiff(genequery,queriesindb)
  run_time_message(paste0("NOT in the gene/circ column: \n" ,paste0(paste(queriesnotindb, sep="", collapse="\n"))))
  run_time_message(paste0("In the gene/circ column: \n",paste0(paste(queriesindb, sep="", collapse="\n"))))
  
  q_cor_summary<-data.frame()  #generate q_cor_summary table 
  #for (i in queriesindb) { system(paste0("echo ",i," >> output/genequery_check.txt"))  }
  for (i in queriesnotindb) { 
    system(paste0("echo \'",i," no found \' >> output/genequery_check.txt"))  ;
    q_cor_summary<-rbind(q_cor_summary,data.frame(query=paste0(i,"(no found)"))) # write not found query
  }
  
  ####### check end ############
  
  
  ###generte linc_coexp_pairs ######
  run_time_message("generate linc_coexp_pairs")
  
  if (opt$symbolid[1] == "symbol") {geneidcolumn = gene_coordinate$gencode_gene_symbol }
  if (opt$symbolid[1] == "id") {geneidcolumn  = gene_coordinate$ensembl_gene_id }
  
  allquery<-data.frame()
  for (q_gene in queriesindb){
    #q_gene="LINC00346"
    query_coordinate<-gene_coordinate[geneidcolumn %in% q_gene,]
    #check if the query gene is unique
    if(nrow(query_coordinate)>1) { system(paste0("echo \'",q_gene," (multiple ensembl id) \' >> output/genequery_check.txt")) ;}
    if(nrow(query_coordinate)>=1) { 
      system(paste0("echo \'",q_gene,"  \' >> output/genequery_check.txt"))  
      q_cor_summary<-rbind(q_cor_summary,data.frame(query=paste0(q_gene))) # for correlation summary
    }
    query_coordinate<-query_coordinate[rep(1, each=nrow(gene_coordinate)),][,c(1,2,3,6,4,8)]
    query_coordinate<-cbind(query_coordinate,gene_coordinate[,c(1,2,3,6,4,8)])
    allquery<-rbind(allquery,query_coordinate)
  }
  
  
  
  colnames(allquery)<-c("chr","lncRNA_start","lncRNA_end","lncRNA_strand","lncRNA_id","lncRNA","co_exp_gene_chr","co_exp_gene_start",
                        "co_exp_gene_end","co_exp_gene_strand","co_exp_gene_id","co_exp_gene")
  linc_coexp_pairs<-allquery 
  #####generte linc_coexp_pairs  END######
  
  
  
  ##### chek min read count#####
  run_time_message("chek min read count")
  if (min(normreads) <= 0) { # add offset value
  run_time_message("Add 1 to normreads if min <=0 ")
    normreads<-normreads+1
  }
  if (min(normreads) < 0) { #stop run if it's log format
    system('echo "Examine the table if contain negative values "')
    system('echo "Dose not support matrix with negative values" > output/errormsg.txt')
    stop("Dose not support logarithmic format")
  }
  ##### chek min read count end####
  
  
  ##### density plot################
  # for test
  #coding_name_id<-codeinggenelist$gencode_gene_symbol
  #lnc_name_id<-lnc_gene$gencode_gene_symbol
 
  run_time_message("density plot")
  coding_gene_name<-codeinggenelist$gencode_gene_symbol
  lncRNA_gene_name<-lnc_gene$gencode_gene_symbol
  
  coding_gene_id<-codeinggenelist$ensembl_gene_id
  lncRNA_gene_id<-lnc_gene$ensembl_gene_id
  
  generate_density<-function(coding_name_id,lnc_name_id){
    run_time_message("Start Distribution of lncRNAs and coding genes")
    normreads_sd<-data.table(normreads[apply(normreads,1,mean)!=min(normreads) ,]) # rm row sum=0)
    #normreads_sd<-melt.data.table(data.table(merge(t(normreads),factor_list,by.x ="row.names",by.y = "V1")))
    normreads_sd<-melt((merge(t(normreads),factor_list,by.x ="row.names",by.y = "V1")))
    colnames(normreads_sd)<-c("samples","condition","gene","normreads")
    normreads_sd<-data.table(normreads_sd)
    normreads_sd[,aaa:=mean(normreads),by=list(gene,condition)]
    normreads_sd<-as.data.table(subset(normreads_sd,select=c(condition,gene,aaa)))
    normreads_sd<-unique(normreads_sd)
    colnames(normreads_sd)[3]<-"normreads"
    # normreads_sd<-aggregate(normreads ~ condition + gene, data = normreads_sd[,2:4], mean)
    #smalldat[, aggGroup1 := mean(x), by = group1]
    meltcodinggene<-normreads_sd[normreads_sd$gene %in% unique(coding_name_id),]
    meltcodinggene$group<-paste0("coding genes ",meltcodinggene$condition)
    meltlncgene<-normreads_sd[normreads_sd$gene %in% unique(lnc_name_id),]
    meltlncgene$group<-paste0("lncRNAs ",meltlncgene$condition)
    mergedmelt<-rbind(meltlncgene,meltcodinggene)
    mergedmelt$normreads<-log2(mergedmelt$normreads)
    
    ggplot(mergedmelt, aes(x=group,y=normreads,fill =group))+
      geom_boxplot(alpha = 0.5)+
      #geom_violin(alpha = 0.5,lwd=0.25,draw_quantiles=c(0.25,0.5,0.75))+
      theme_bw(base_size = 25)+ #15
      coord_flip()+
      theme(
        aspect.ratio = 0.75,
        legend.position="none"
      )+
      xlab(bquote("")) +
      ylab(bquote("Normalized Reads("~log[2]~")"))
    ggsave(paste0("output/boxplot.png"), dpi=300, width = 8, height = 6)
  }
  
  if (opt$symbolid[1]=="symbol"  ){
    generate_density(coding_gene_name,lncRNA_gene_name)
    
  } else if (opt$symbolid[1]=="id" ){
    generate_density(coding_gene_id,lncRNA_gene_id)
  }
  
  ##### density plot end################
  
  ##### Add distance information########

  run_time_message("Add distance information")
  lncRNAAA<-linc_coexp_pairs$chr
  coexpp<-linc_coexp_pairs$co_exp_gene_chr
  disttt<-linc_coexp_pairs$distance
  lncS<-linc_coexp_pairs$lncRNA_start
  lncE<-linc_coexp_pairs$lncRNA_end
  coS<-linc_coexp_pairs$co_exp_gene_start
  coE<-linc_coexp_pairs$co_exp_gene_end
  
  linc_coexp_pairs$distance<-ifelse(linc_coexp_pairs$lncRNA==linc_coexp_pairs$co_exp_gene,0,
                                    ifelse(linc_coexp_pairs$chr!=linc_coexp_pairs$co_exp_gene_chr,"Different chr",
                                    ifelse(  lncS >= coE & lncE >= coE , paste0("Upstream ",(lncS-coE)/1000," kbp(s)"),
                                             ifelse(  lncE <= coS & lncE <= coE , paste0("Downstream ",(coS-lncE)/1000," kbp(s)"),
                                                      ifelse(  lncS >= coS & lncE >= coE & lncS <= coE & lncE >= coS ,paste0("Upstream overlapped"),
                                                               ifelse(  lncS <= coS & lncE <= coE & lncS <= coE & lncE >= coS ,paste0("Downstream overlapped") ,
                                                                        ifelse(  lncS <= coE & lncS <= coS & lncE >= coS & lncE >= coE,paste0("co-express gene residents query gene"),
                                                                                 ifelse(  lncS >= coS & lncS <= coE & lncE >= coS & lncE <= coE , paste0("query gene residents co-express gene"),"NA"))))))))
  linc_coexp_pairs$query_Lncipedia_HC<-ifelse(linc_coexp_pairs$lncRNA_id %in% lncpedia$V1,"Yes","No")
  linc_coexp_pairs$co_exp_gene_Lncipedia_HC<-ifelse(linc_coexp_pairs$co_exp_gene_id %in% lncpedia$V1,"Yes","No")
  
  #### add distance end #########

  run_time_message("Calculating Pearson correlation")
  if (opt$symbolid[1]=="symbol"  ){
    sub_normreads<-subset(normreads,rownames(normreads) %in% unique(c(linc_coexp_pairs[,"lncRNA"],linc_coexp_pairs[,"co_exp_gene"])))
    
  } else if (opt$symbolid[1]=="id" ){
    sub_normreads<-subset(normreads,rownames(normreads) %in% unique(c(linc_coexp_pairs[,"lncRNA_id"],linc_coexp_pairs[,"co_exp_gene_id"])))
  }
  
  ###### Pearson correlation   #######
  run_time_message("All samples Pearson correlation")
  
  # change the selection of column accroading the input gene symbol or id
  if (opt$symbolid[1]=="symbol" ){
    lncRNAcolumn<-linc_coexp_pairs[,"lncRNA"]
    coexpcolumn<-linc_coexp_pairs[,"co_exp_gene"]
    coexpcolumn_str<-"co_exp_gene"
    lncRNAcolumn_str<-"lncRNA"
  } else if (opt$symbolid[1]=="id"){    
    lncRNAcolumn<-linc_coexp_pairs[,"lncRNA_id"]
    coexpcolumn<-linc_coexp_pairs[,"co_exp_gene_id"]
    coexpcolumn_str<-"co_exp_gene_id"
    lncRNAcolumn_str<-"lncRNA_id"
  }
  
  corlist<-WGCNA::corAndPvalue(t(log2(sub_normreads)[unique(lncRNAcolumn),]),t(log2(sub_normreads)),nThreads = 6 ,method = "pearson", verbose=1) # calcualte pearson correlation
  
  for (i in  1:length(unique(lncRNAcolumn))) {
    run_time_message(paste0("No. of query :", i))
     cor_mat<-melt(corlist$cor) 
    colnames(cor_mat)<-c("query","co_exp_gene","cor")
    cor_mat[,1:2]<-sapply(cor_mat[,1:2],as.character)
    #cor_mat<-cor_mat[ !(cor_mat$co_exp_gene %in% unique(cor_mat$query)),]
    cor_p<-melt(corlist$p) #extract correlation pvalue
    colnames(cor_p)<-c("query","co_exp_gene","cor_p")
    cor_p[,1:2]<-sapply(cor_p[,1:2],as.character)
    #cor_p<-cor_p[ !(cor_p$co_exp_gene %in% unique(cor_p$query)),]
    x<-merge(linc_coexp_pairs[lncRNAcolumn==unique(lncRNAcolumn)[i],1:12], 
             cor_mat[cor_mat$query ==unique(lncRNAcolumn)[i] & cor_mat$co_exp_gene %in% unique(coexpcolumn),][,2:3],by.x=coexpcolumn_str, by.y="co_exp_gene",all.x=T,sort=F)
    linc_coexp_pairs[lncRNAcolumn==unique(lncRNAcolumn)[i],"cor"]<-x[order(match(x$co_exp_gene_id,linc_coexp_pairs[lncRNAcolumn==unique(lncRNAcolumn)[i],]$co_exp_gene_id)),]$cor
    
    x<-merge(linc_coexp_pairs[lncRNAcolumn==unique(lncRNAcolumn)[i],1:12], 
             cor_p[cor_p$query ==unique(lncRNAcolumn)[i] & cor_p$co_exp_gene %in% unique(coexpcolumn),][,2:3],by.x=coexpcolumn_str, by.y="co_exp_gene",all.x=T,sort=F)
    linc_coexp_pairs[lncRNAcolumn==unique(lncRNAcolumn)[i],"cor_p"]<-x[order(match(x$co_exp_gene_id,linc_coexp_pairs[lncRNAcolumn==unique(lncRNAcolumn)[i],]$co_exp_gene_id)),]$cor_p
    
    
    # for ( factor_levels in 1:length(levels(factor_list$V2))){ # calculation condition correlation
    #   # calculation of factor pearson
    #  
    #   run_time_message(paste0("condition of sample :", paste0(levels(factor_list$V2)[factor_levels])))
    #   
    #   corlist<-WGCNA::corAndPvalue(t(log2(sub_normreads)[unique(lncRNAcolumn),colnames(log2(sub_normreads)) %in% as.character(factor_list[factor_list$V2==levels(factor_list$V2)[factor_levels],1])]),
    #                         t(log2(sub_normreads)[,colnames(log2(sub_normreads)) %in% as.character(factor_list[factor_list$V2==levels(factor_list$V2)[factor_levels],1])]),nThreads = 6 ,method = "pearson", verbose=1)
    #   cor_mat<-melt(corlist$cor)
    #   colnames(cor_mat)<-c("query","co_exp_gene",paste0(levels(factor_list$V2)[factor_levels],"_cor"))
    #   cor_mat[,1:2]<-sapply(cor_mat[,1:2],as.character)
    #  # cor_mat<-cor_mat[ !(cor_mat$co_exp_gene %in% unique(cor_mat$query)),]
    #   cor_p<-melt(corlist$p)
    #   colnames(cor_p)<-c("query","co_exp_gene",paste0(levels(factor_list$V2)[factor_levels],"_cor_pvalue"))
    #   
    #   x<-merge(linc_coexp_pairs[lncRNAcolumn==unique(lncRNAcolumn)[i],1:12], 
    #            cor_mat[cor_mat$query ==unique(lncRNAcolumn)[i] & cor_mat$co_exp_gene %in% unique(coexpcolumn),][,2:3],by.x=coexpcolumn_str, by.y="co_exp_gene",all.x=T,sort=F)
    #   linc_coexp_pairs[lncRNAcolumn==unique(lncRNAcolumn)[i],paste0(levels(factor_list$V2)[factor_levels],"_cor")]<-
    #     x[order(match(x$co_exp_gene_id,linc_coexp_pairs[lncRNAcolumn==unique(lncRNAcolumn)[i],]$co_exp_gene_id)),][,paste0(levels(factor_list$V2)[factor_levels],"_cor")]
    #   
    #   x<-merge(linc_coexp_pairs[lncRNAcolumn==unique(lncRNAcolumn)[i],1:12], 
    #            cor_p[cor_p$query ==unique(lncRNAcolumn)[i] & cor_p$co_exp_gene %in% unique(coexpcolumn),][,2:3],by.x=coexpcolumn_str, by.y="co_exp_gene",all.x=T,sort=F)
    #   linc_coexp_pairs[lncRNAcolumn==unique(lncRNAcolumn)[i],paste0(levels(factor_list$V2)[factor_levels],"_cor_pvalue")]<-
    #     x[order(match(x$co_exp_gene_id,linc_coexp_pairs[lncRNAcolumn==unique(lncRNAcolumn)[i],]$co_exp_gene_id)),][,paste0(levels(factor_list$V2)[factor_levels],"_cor_pvalue")]
    #   
    #   
    # }
    # add cor summary to q_cor_summary
    qq<-unique(lncRNAcolumn)[i]  # easily    indorduce bug
    x<-subset(linc_coexp_pairs[lncRNAcolumn==qq ,], abs(cor) > 0.5)
    q_cor_summary[q_cor_summary$query==qq,"Fraction of absolute cor > 0.5"]<-round(nrow(x)/nrow(linc_coexp_pairs[lncRNAcolumn==qq ,]),2)
    x<-subset(linc_coexp_pairs[lncRNAcolumn==unique(lncRNAcolumn)[i] ,], cor_p < 0.05 )
    q_cor_summary[q_cor_summary$query==qq,"Fraction of cor pvalue < 0.05"]<-round(nrow(x)/nrow(linc_coexp_pairs[lncRNAcolumn==unique(lncRNAcolumn)[i] ,]),2)  
  }
  
  for ( factor_levels in 1:length(levels(factor_list$V2))){ # calculation condition correlation
    # calculation of factor pearson
    
    run_time_message(paste0("condition of sample :", paste0(levels(factor_list$V2)[factor_levels])))

    corlist<-WGCNA::corAndPvalue(t(log2(sub_normreads)[unique(lncRNAcolumn),colnames(log2(sub_normreads)) %in% as.character(factor_list[factor_list$V2==levels(factor_list$V2)[factor_levels],1])]),
                          t(log2(sub_normreads)[,colnames(log2(sub_normreads)) %in% as.character(factor_list[factor_list$V2==levels(factor_list$V2)[factor_levels],1])]),nThreads = 6 ,method = "pearson", verbose=1)
    for (i in  1:length(unique(lncRNAcolumn))) {
    run_time_message(paste0("No. of query :", i))
    
    cor_mat<-melt(corlist$cor)
    colnames(cor_mat)<-c("query","co_exp_gene",paste0(levels(factor_list$V2)[factor_levels],"_cor"))
    cor_mat[,1:2]<-sapply(cor_mat[,1:2],as.character)
   # cor_mat<-cor_mat[ !(cor_mat$co_exp_gene %in% unique(cor_mat$query)),]
    cor_p<-melt(corlist$p)
    colnames(cor_p)<-c("query","co_exp_gene",paste0(levels(factor_list$V2)[factor_levels],"_cor_pvalue"))

    x<-merge(linc_coexp_pairs[lncRNAcolumn==unique(lncRNAcolumn)[i],1:12],
             cor_mat[cor_mat$query ==unique(lncRNAcolumn)[i] & cor_mat$co_exp_gene %in% unique(coexpcolumn),][,2:3],by.x=coexpcolumn_str, by.y="co_exp_gene",all.x=T,sort=F)
    linc_coexp_pairs[lncRNAcolumn==unique(lncRNAcolumn)[i],paste0(levels(factor_list$V2)[factor_levels],"_cor")]<-
      x[order(match(x$co_exp_gene_id,linc_coexp_pairs[lncRNAcolumn==unique(lncRNAcolumn)[i],]$co_exp_gene_id)),][,paste0(levels(factor_list$V2)[factor_levels],"_cor")]

    x<-merge(linc_coexp_pairs[lncRNAcolumn==unique(lncRNAcolumn)[i],1:12],
             cor_p[cor_p$query ==unique(lncRNAcolumn)[i] & cor_p$co_exp_gene %in% unique(coexpcolumn),][,2:3],by.x=coexpcolumn_str, by.y="co_exp_gene",all.x=T,sort=F)
    linc_coexp_pairs[lncRNAcolumn==unique(lncRNAcolumn)[i],paste0(levels(factor_list$V2)[factor_levels],"_cor_pvalue")]<-
      x[order(match(x$co_exp_gene_id,linc_coexp_pairs[lncRNAcolumn==unique(lncRNAcolumn)[i],]$co_exp_gene_id)),][,paste0(levels(factor_list$V2)[factor_levels],"_cor_pvalue")]

    }
  }
  
  ###### END of Pearson correlation   #######
  
  
  
  ###### add reads mean information of lncRNA and co express genes ######
  
  run_time_message("calculating mean reads for all samples")
  xx<-as.data.frame(apply(sub_normreads[,],1,mean))
  colnames(xx)<-"query_normalized_mean"
  linc_coexp_pairs<-merge(linc_coexp_pairs,xx,by.x = lncRNAcolumn_str, by.y="row.names",all.x=T)
  colnames(xx)<-"co_exp_gene_normalized_mean"
  linc_coexp_pairs<-merge(linc_coexp_pairs,xx,by.x = coexpcolumn_str, by.y="row.names",all.x=T)
  
  
  ## add conditional mean reads
  run_time_message("calculating mean reads for conditional samples")
  # add reads mean information of lncRNA and co express genes accroding factors
  for ( factor_levels in 1:length(levels(factor_list$V2))){
    xx<-as.data.frame(apply(sub_normreads[,colnames(sub_normreads) %in% as.character(factor_list[factor_list$V2==levels(factor_list$V2)[factor_levels],1])],1,mean))
    colnames(xx)<-paste0(levels(factor_list$V2)[factor_levels],"_query_normalized_mean")
    linc_coexp_pairs<-merge(linc_coexp_pairs,xx,by.x = lncRNAcolumn_str, by.y="row.names",all.x=T)
    colnames(xx)<-paste0(levels(factor_list$V2)[factor_levels],"_co_exp_gene_normalized_mean")
    linc_coexp_pairs<-merge(linc_coexp_pairs,xx,by.x = coexpcolumn_str, by.y="row.names",all.x=T)
  }
  ##### END of add reads mean information of lncRNA and co express genes ######
  
  linc_coexp_pairs_write<-linc_coexp_pairs[,c(2,3,4,5,6,7,1,seq(8,ncol(linc_coexp_pairs)))]
  linc_coexp_pairs<-linc_coexp_pairs[,c(2,3,4,5,6,7,1,seq(8,ncol(linc_coexp_pairs)))]
  linc_coexp_pairs<-linc_coexp_pairs[order(linc_coexp_pairs$lncRNA),]
  
  run_time_message("Write linc_coexp_pairs.rds")
  saveRDS(linc_coexp_pairs, paste0("output/linc_coexp_pairs.rds"))
  write_feather(linc_coexp_pairs, paste0("output/linc_coexp_pairs.feather"))
  # saveRDS(normreads, paste0("output/normreads.rds"))
  saveRDS(q_cor_summary, paste0("output/q_cor_summary.rds"))
  
  write.table(linc_coexp_pairs,"output/linc_coexp_pairs.txt",sep = '\t',row.names = F,quote = F)
  write.table(q_cor_summary,"output/q_cor_summary.txt",sep = '\t',row.names = F,quote = F)
  
  circ_gene_merged<-normreads
  run_time_message("Write circ_gene_merged_not.rds")
  saveRDS(circ_gene_merged, paste0("output/circ_gene_merged_not.rds"))
  circ_gene_merged<-as.data.frame(t(circ_gene_merged))
  circ_gene_merged<-merge(factor_list,circ_gene_merged, by.x= "V1" , by.y="row.names")
  row.names(circ_gene_merged)<-circ_gene_merged$V1
  circ_gene_merged$V1<-NULL
  colnames(circ_gene_merged)[1]<-"attr"
  run_time_message("Write circ_gene_merged.rds")
  saveRDS(circ_gene_merged, paste0("output/circ_gene_merged.rds"))
  run_time_message("finish")
 # write_feather(circ_gene_merged, paste0("output/circ_gene_merged.feather"))
}

#for_lnc_mode()

for_circ_mode<-function(){
  ####### circ RNA ##########
  ####### check table start ############
  run_time_message("circRNA analysis start")
  if ( sum(length(setdiff(colnames(normreads [,1:ncol(normreads )]),factor_list$V1)),
           length(setdiff(factor_list$V1,colnames(normreads [,1:ncol(normreads )]))))!=0 ){  
    system("echo Unmacthed sample names, factor list and normalized reads table > output/errormsg.txt"); q();}
  
  if ( sum(length(setdiff(colnames(circrna [,1:ncol(circrna )]),factor_list$V1)),
           length(setdiff(factor_list$V1,colnames(circrna [,1:ncol(circrna )]))))!=0 ){  
    system("echo Unmacthed sample names, factor list and circRNA read table > output/errormsg.txt"); q();}
  
  if ( length(unique(factor_list$V2)) != 2 ){  
    system("echo Only support 2 levels factor > output/errormsg.txt"); q(); }
  
  if (  length(factor_list$V2[grep(levels(factor_list$V2)[1],factor_list$V2)]) <3 &   
        length(factor_list$V2[grep(levels(factor_list$V2)[2],factor_list$V2)]) <3  ) { system("echo Not enought replicates > output/errormsg.txt"); q(); }
  
  # add functin to check if its really id or gene name
  genequery<-unique(unlist(strsplit(opt$query[1],",")[1]))
  
  for (i in genequery){
  } 
  if (opt$symbolid[1] == "id"){  rownames(normreads) <-gsub("\\..*","",rownames(normreads)) }
  if (opt$symbolid[1] == "id"){  rownames(circrna) <-gsub("\\..*","",rownames(circrna)) } # for test only 
  
  if (opt$symbolid[1] == "id"){  genequery <-gsub("\\..*","",genequery) } # for test only 
  
  #check if query in the circ rna table
  queriesindb<-genequery[genequery %in% rownames(circrna)] 
  queriesindb<-intersect(genequery,rownames(circrna))
  queriesnotindb <- setdiff(genequery,queriesindb)
  #queriesindb<-genequery
 
  run_time_message(paste0("NOT in the gene/circ column: \n" ,paste0(paste(queriesnotindb, sep="", collapse="\n"))))
  run_time_message(paste0("In the gene/circ column: \n",paste0(paste(queriesindb, sep="", collapse="\n"))))
  
  system(paste0("echo ",length(genequery)," in genequery\n")) 
  
  q_cor_summary<-data.frame()  #generate q_cor_summary table 
  #for (i in queriesindb) { system(paste0("echo ",i," >> output/genequery_check.txt"))  }
  for (i in queriesnotindb) { 
    system(paste0("echo ",i," \\(no found\\) >> output/genequery_check.txt"))  ;
    q_cor_summary<-rbind(q_cor_summary,data.frame(query=paste0(i,"(no found)"))) # write not found query
  }
  
  
  ####### check end ############
  
  
  ###generte linc_coexp_pairs ######
  
  
  if (opt$symbolid[1] == "symbol") {geneidcolumn = gene_coordinate$gencode_gene_symbol }
  if (opt$symbolid[1] == "id") {geneidcolumn  = gene_coordinate$ensembl_gene_id }
  
  allquery<-data.frame()
  for (q_gene in queriesindb){
    #q_gene="LINC00346"
    query_coordinate<-data.frame(circRNA=q_gene,stringsAsFactors=F) #<-gene_coordinate[geneidcolumn %in% q_gene,]
    
    #check if the query gene is unique
    if(nrow(query_coordinate)>1) { system(paste0("echo \'",q_gene," (multiple ensembl id) \' >> output/genequery_check.txt")) ;}
    if(nrow(query_coordinate)>=1) { 
      system(paste0("echo \'",q_gene,"\' >> output/genequery_check.txt")) ;
      q_cor_summary<-rbind(q_cor_summary,data.frame(query=paste0(q_gene))  ) # for correlation summary
    }
    query_coordinate<-data.frame(circRNA=query_coordinate[rep(1, each=nrow(gene_coordinate)),],stringsAsFactors=F) #will be
    query_coordinate<-cbind(query_coordinate,gene_coordinate[,c(1,2,3,6,4,8)])
    allquery<-rbind(allquery,query_coordinate)
  }
  
  run_time_message("add allquery column names")
  #colnames(allquery)<-c("chr","lncRNA_start","lncRNA_end","lncRNA_strand","lncRNA_id","lncRNA","co_exp_gene_chr","co_exp_gene_start",
  #                      "co_exp_gene_end","co_exp_gene_strand","co_exp_gene_id","co_exp_gene")
  colnames(allquery)<-c("circRNA","co_exp_gene_chr","co_exp_gene_start",
                        "co_exp_gene_end","co_exp_gene_strand","co_exp_gene_id","co_exp_gene")
  
  linc_coexp_pairs<-allquery 
  #####generte linc_coexp_pairs  END######
  
  
  ##### chek min read count#####
  if (min(normreads) <= 0) { # add offset value
    run_time_message("Add 1 to normreads if min <=0")
    normreads<-normreads+1
  }
  if (min(normreads) < 0) { #stop run if it's log format
    run_time_message("Examine the table if contain negative values")
    system('echo "Dose not support matrix with negative values" > output/errormsg.txt')
    stop("Dose not support logarithmic format")
  }
  
  if (min(circrna) <= 0) { # add offset value
    run_time_message("Add 1 to circrnas if min <=0")
    circrna<-circrna+1
  }
  if (min(circrna) < 0) { #stop run if it's log format
    run_time_message("Examine the table if contain negative values")
    system('echo "Dose not support matrix with negative values" > output/errormsg.txt')
    stop("Dose not support logarithmic format")
  }
  
  ##### chek min read count end####
  
 #####add gene abundance dstribution ##### 
  run_time_message("density plot")
  coding_gene_name<-codeinggenelist$gencode_gene_symbol
  lncRNA_gene_name<-lnc_gene$gencode_gene_symbol
  
  coding_gene_id<-codeinggenelist$ensembl_gene_id
  lncRNA_gene_id<-lnc_gene$ensembl_gene_id
  
  generate_density<-function(coding_name_id,lnc_name_id){
    system('echo "Distribution of lncRNAs and coding genes"')
    normreads_sd<-data.table(normreads[apply(normreads,1,mean)!=min(normreads) ,]) # rm row sum=0)
    #normreads_sd<-melt.data.table(data.table(merge(t(normreads),factor_list,by.x ="row.names",by.y = "V1")))
    normreads_sd<-melt((merge(t(normreads),factor_list,by.x ="row.names",by.y = "V1")))
    colnames(normreads_sd)<-c("samples","condition","gene","normreads")
    normreads_sd<-data.table(normreads_sd)
    normreads_sd[,aaa:=mean(normreads),by=list(gene,condition)]
    normreads_sd<-as.data.table(subset(normreads_sd,select=c(condition,gene,aaa)))
    normreads_sd<-unique(normreads_sd)
    colnames(normreads_sd)[3]<-"normreads"
    # normreads_sd<-aggregate(normreads ~ condition + gene, data = normreads_sd[,2:4], mean)
    #smalldat[, aggGroup1 := mean(x), by = group1]
    meltcodinggene<-normreads_sd[normreads_sd$gene %in% unique(coding_name_id),]
    meltcodinggene$group<-paste0("coding genes ",meltcodinggene$condition)
    meltlncgene<-normreads_sd[normreads_sd$gene %in% unique(lnc_name_id),]
    meltlncgene$group<-paste0("lncRNAs ",meltlncgene$condition)
    mergedmelt<-rbind(meltlncgene,meltcodinggene)
   
    
    circrna_sd<-melt((merge(t(circrna),factor_list,by.x ="row.names",by.y = "V1")))
    colnames(circrna_sd)<-c("samples","condition","gene","normreads")
    circrna_sd<-data.table(circrna_sd)
    circrna_sd[,aaa:=mean(normreads),by=list(gene,condition)]
    circrna_sd<-as.data.table(subset(circrna_sd,select=c(condition,gene,aaa)))
    circrna_sd<-unique(normreads_sd)
    colnames(circrna_sd)[3]<-"normreads"
    circrna_sd$group<-paste0("circRNAs ",circrna_sd$condition)
    mergedmelt<-rbind(circrna_sd,mergedmelt)
    
    mergedmelt$normreads<-log2(mergedmelt$normreads)
    
    ggplot(mergedmelt, aes(x=group,y=normreads,fill =group))+
      geom_boxplot(alpha = 0.5)+
      #geom_violin(alpha = 0.5,lwd=0.25,draw_quantiles=c(0.25,0.5,0.75))+
      theme_bw(base_size = 25)+ #15
      coord_flip()+
      theme(
        aspect.ratio = 0.75,
        legend.position="none"
      )+
      xlab(bquote("")) +
      ylab(bquote("Normalized Reads("~log[2]~")"))
    ggsave(paste0("output/boxplot.png"), dpi=300, width = 8, height = 6)
  }
  
  if (opt$symbolid[1]=="symbol"  ){
    generate_density(coding_gene_name,lncRNA_gene_name)
    
  } else if (opt$symbolid[1]=="id" ){
    generate_density(coding_gene_id,lncRNA_gene_id)
  }
  #####add gene abundance dstribution  END##### 
  
  ##### Add distance information########

  asda<-function(){lncRNAAA<-linc_coexp_pairs$chr
  coexpp<-linc_coexp_pairs$co_exp_gene_chr
  disttt<-linc_coexp_pairs$distance
  lncS<-linc_coexp_pairs$lncRNA_start
  lncE<-linc_coexp_pairs$lncRNA_end
  coS<-linc_coexp_pairs$co_exp_gene_start
  coE<-linc_coexp_pairs$co_exp_gene_end
  
  linc_coexp_pairs$distance<-ifelse(linc_coexp_pairs$chr!=linc_coexp_pairs$co_exp_gene_chr,"Different chr",
                                    ifelse(  lncS >= coE & lncE >= coE , paste0("Upstream ",(lncS-coE)/1000," kbp(s)"),
                                             ifelse(  lncE <= coS & lncE <= coE , paste0("Downstream ",(coS-lncE)/1000," kbp(s)"),
                                                      ifelse(  lncS >= coS & lncE >= coE & lncS <= coE & lncE >= coS ,paste0("Upstream overlapped"),
                                                               ifelse(  lncS <= coS & lncE <= coE & lncS <= coE & lncE >= coS ,paste0("Downstream overlapped") ,
                                                                        ifelse(  lncS <= coE & lncS <= coS & lncE >= coS & lncE >= coE,paste0("co-express gene residents query gene"),
                                                                                 ifelse(  lncS >= coS & lncS <= coE & lncE >= coS & lncE <= coE , paste0("query gene residents co-express gene"),"NA")))))))
  linc_coexp_pairs$query_Lncipedia_HC<-ifelse(linc_coexp_pairs$lncRNA_id %in% lncpedia$V1,"Yes","No")
  linc_coexp_pairs$co_exp_gene_Lncipedia_HC<-ifelse(linc_coexp_pairs$co_exp_gene_id %in% lncpedia$V1,"Yes","No")}
  
  #### add distance end #########

  run_time_message("Calculating Pearson correlation")
  if (opt$symbolid[1]=="symbol"  ){
    sub_normreads<-subset(normreads,rownames(normreads) %in% unique(linc_coexp_pairs[,"co_exp_gene"]))
    
  } else if (opt$symbolid[1]=="id" ){
    sub_normreads<-subset(normreads,rownames(normreads) %in% unique(linc_coexp_pairs[,"co_exp_gene_id"]))
  }
  
  ###### Pearson correlation   #######
  
  if (opt$symbolid[1]=="symbol" ){
    querycolumn<-linc_coexp_pairs[,"circRNA"]
    coexpcolumn<-linc_coexp_pairs[,"co_exp_gene"]
    coexpcolumn_str<-"co_exp_gene"
    querycolumn_str<-"circRNA"
  } else if (opt$symbolid[1]=="id"){
    querycolumn<-linc_coexp_pairs[,"circRNA"]
    coexpcolumn<-linc_coexp_pairs[,"co_exp_gene_id"]
    coexpcolumn_str<-"co_exp_gene_id"
    querycolumn_str<-"circRNA"
  }
  
  circrna<-circrna[,colnames(sub_normreads)] # sort circrna column by sub_normreads coulmn
  
  for (i in  1:length(unique(querycolumn))) {
    run_time_message(paste0("No. of query :", i))
    #linc_coexp_pairs<-linc_coexp_pairs[,1:7] for test
    corlist<-WGCNA::corAndPvalue(t(log2(circrna)[unique(querycolumn),]),t(log2(sub_normreads)),nThreads = 6 ,method = "pearson", verbose=1)
    cor_mat<-melt(corlist$cor)
    colnames(cor_mat)<-c("query","co_exp_gene","cor")
    cor_mat[,1:2]<-sapply(cor_mat[,1:2],as.character)
    #cor_mat<-cor_mat[ !(cor_mat$co_exp_gene %in% unique(cor_mat$query)),]
    cor_p<-melt(corlist$p)
    colnames(cor_p)<-c("query","co_exp_gene","cor_p")
    cor_p[,1:2]<-sapply(cor_p[,1:2],as.character)
    #cor_p<-cor_p[ !(cor_p$co_exp_gene %in% unique(cor_p$query)),]
    x<-merge(linc_coexp_pairs[querycolumn==unique(querycolumn)[i],1:7],  # merge correlation to linc_coexp_pairs
             cor_mat[cor_mat$query ==unique(querycolumn)[i] & cor_mat$co_exp_gene %in% unique(coexpcolumn),][,2:3],by.x=coexpcolumn_str, by.y="co_exp_gene",all.x=T,sort=F)
    linc_coexp_pairs[querycolumn==unique(querycolumn)[i],"cor"]<-x[order(match(x$co_exp_gene_id,linc_coexp_pairs[querycolumn==unique(querycolumn)[i],]$co_exp_gene_id)),]$cor
    
    x<-merge(linc_coexp_pairs[querycolumn==unique(querycolumn)[i],1:7], 
             cor_p[cor_p$query ==unique(querycolumn)[i] & cor_p$co_exp_gene %in% unique(coexpcolumn),][,2:3],by.x=coexpcolumn_str, by.y="co_exp_gene",all.x=T,sort=F)
    linc_coexp_pairs[querycolumn==unique(querycolumn)[i],"cor_p"]<-x[order(match(x$co_exp_gene_id,linc_coexp_pairs[querycolumn==unique(querycolumn)[i],]$co_exp_gene_id)),]$cor_p
    
    
    for ( factor_levels in 1:length(levels(factor_list$V2))){ # calculation condition correlation
      # calculation of factor pearson
      run_time_message(paste0("condition of sample :", paste0(levels(factor_list$V2)[factor_levels])))
      corlist<-WGCNA::corAndPvalue(t(log2(circrna)[unique(querycolumn),colnames(log2(circrna)) %in% as.character(factor_list[factor_list$V2==levels(factor_list$V2)[factor_levels],1])]),
                            t(log2(sub_normreads)[,colnames(log2(sub_normreads)) %in% as.character(factor_list[factor_list$V2==levels(factor_list$V2)[factor_levels],1])]),nThreads = 6 ,method = "pearson", verbose=1)
      cor_mat<-melt(corlist$cor)
      colnames(cor_mat)<-c("query","co_exp_gene",paste0(levels(factor_list$V2)[factor_levels],"_cor"))
      cor_mat[,1:2]<-sapply(cor_mat[,1:2],as.character)
      #cor_mat<-cor_mat[ !(cor_mat$co_exp_gene %in% unique(cor_mat$query)),]
      cor_p<-melt(corlist$p)
      colnames(cor_p)<-c("query","co_exp_gene",paste0(levels(factor_list$V2)[factor_levels],"_cor_pvalue"))
      
      #x<-merge(linc_coexp_pairs[linc_coexp_pairs$lncRNA==unique(linc_coexp_pairs$lncRNA)[i],1:ncol(linc_coexp_pairs)], 
      #         subset(cor_mat,query ==unique(linc_coexp_pairs$lncRNA)[i] & co_exp_gene %in% unique(linc_coexp_pairs$co_exp_gene) )[,2:3],by.x="co_exp_gene", by.y="co_exp_gene",all.x=T,sort=F)
      #linc_coexp_pairs[linc_coexp_pairs$lncRNA==unique(linc_coexp_pairs$lncRNA)[i],paste0(levels(factor_list$V2)[factor_levels],"_cor")]<-
      #  x[order(match(x$co_exp_gene_id,subset(linc_coexp_pairs,lncRNA==unique(linc_coexp_pairs$lncRNA)[i])$co_exp_gene_id)),][,paste0(levels(factor_list$V2)[factor_levels],"_cor")]
      
      # x<-merge(linc_coexp_pairs[linc_coexp_pairs$lncRNA==unique(linc_coexp_pairs$lncRNA)[i],1:ncol(linc_coexp_pairs)], 
      #           subset(cor_p,query ==unique(linc_coexp_pairs$lncRNA)[i] & co_exp_gene %in% unique(linc_coexp_pairs$co_exp_gene))[,2:3],by.x="co_exp_gene", by.y="co_exp_gene",all.x=T,sort=F)
      #  linc_coexp_pairs[linc_coexp_pairs$lncRNA==unique(linc_coexp_pairs$lncRNA)[i],paste0(levels(factor_list$V2)[factor_levels],"_cor_pvalue")]<-
      #    x[order(match(x$co_exp_gene_id,subset(linc_coexp_pairs,lncRNA==unique(linc_coexp_pairs$lncRNA)[i])$co_exp_gene_id)),][,paste0(levels(factor_list$V2)[factor_levels],"_cor_pvalue")]
      
      
      x<-merge(linc_coexp_pairs[querycolumn==unique(querycolumn)[i],1:7], 
               cor_mat[cor_mat$query ==unique(querycolumn)[i] & cor_mat$co_exp_gene %in% unique(coexpcolumn),][,2:3],by.x=coexpcolumn_str, by.y="co_exp_gene",all.x=T,sort=F)
      linc_coexp_pairs[querycolumn==unique(querycolumn)[i],paste0(levels(factor_list$V2)[factor_levels],"_cor")]<-
        x[order(match(x$co_exp_gene_id,linc_coexp_pairs[querycolumn==unique(querycolumn)[i],]$co_exp_gene_id)),][,paste0(levels(factor_list$V2)[factor_levels],"_cor")]
      
      x<-merge(linc_coexp_pairs[querycolumn==unique(querycolumn)[i],1:7], 
               cor_p[cor_p$query ==unique(querycolumn)[i] & cor_p$co_exp_gene %in% unique(coexpcolumn),][,2:3],by.x=coexpcolumn_str, by.y="co_exp_gene",all.x=T,sort=F)
      linc_coexp_pairs[querycolumn==unique(querycolumn)[i],paste0(levels(factor_list$V2)[factor_levels],"_cor_pvalue")]<-
        x[order(match(x$co_exp_gene_id,linc_coexp_pairs[querycolumn==unique(querycolumn)[i],]$co_exp_gene_id)),][,paste0(levels(factor_list$V2)[factor_levels],"_cor_pvalue")]
      
     
      
    }
    
    qq<-unique(linc_coexp_pairs[querycolumn==unique(querycolumn)[i] ,]$circRNA)
    x<-subset(linc_coexp_pairs[querycolumn==unique(querycolumn)[i] ,], abs(cor) > 0.5)
    q_cor_summary[q_cor_summary$query==qq,"Fraction of absolute cor > 0.5"]<-round(nrow(x)/nrow(linc_coexp_pairs[querycolumn==unique(querycolumn)[i] ,]),2)
    x<-subset(linc_coexp_pairs[querycolumn==unique(querycolumn)[i] ,], cor_p < 0.05 )
    q_cor_summary[q_cor_summary$query==qq,"Fraction of cor pvalue < 0.05"]<-round(nrow(x)/nrow(linc_coexp_pairs[querycolumn==unique(querycolumn)[i] ,]),2)
  }
  
  ###### END of Pearson correlation   #######
  
  
  
  ###### add reads mean information of lncRNA and co express genes ######
  
  run_time_message("calculating mean reads for all samples")
  xx<-as.data.frame(apply(circrna[,],1,mean)) #circRNA table
  colnames(xx)<-"query_normalized_mean"
  linc_coexp_pairs<-merge(linc_coexp_pairs,xx,by.x = querycolumn_str, by.y="row.names",all.x=T)
  xx<-as.data.frame(apply(sub_normreads[,],1,mean))
  colnames(xx)<-"co_exp_gene_normalized_mean"
  linc_coexp_pairs<-merge(linc_coexp_pairs,xx,by.x = coexpcolumn_str, by.y="row.names",all.x=T)
  
  
  ## add conditional mean reads
  run_time_message("calculating mean reads for conditional samples")
  # add reads mean information of lncRNA and co express genes accroding factors
  for ( factor_levels in 1:length(levels(factor_list$V2))){
    xx<-as.data.frame(apply(circrna[,colnames(circrna) %in% as.character(factor_list[factor_list$V2==levels(factor_list$V2)[factor_levels],1])],1,mean)) #circRNA table
    colnames(xx)<-paste0(levels(factor_list$V2)[factor_levels],"_query_normalized_mean")
    linc_coexp_pairs<-merge(linc_coexp_pairs,xx,by.x = querycolumn_str, by.y="row.names",all.x=T)
    xx<-as.data.frame(apply(sub_normreads[,colnames(sub_normreads) %in% as.character(factor_list[factor_list$V2==levels(factor_list$V2)[factor_levels],1])],1,mean))
    colnames(xx)<-paste0(levels(factor_list$V2)[factor_levels],"_co_exp_gene_normalized_mean")
    linc_coexp_pairs<-merge(linc_coexp_pairs,xx,by.x = coexpcolumn_str, by.y="row.names",all.x=T)
  }
  #linc_coexp_pairs_write<-linc_coexp_pairs[,c(2,3,4,5,6,7,1,seq(8,ncol(linc_coexp_pairs)))]
  #linc_coexp_pairs<-linc_coexp_pairs[,c(2,3,4,5,6,7,1,seq(8,ncol(linc_coexp_pairs)))]
  linc_coexp_pairs<-linc_coexp_pairs[order(linc_coexp_pairs[,querycolumn_str]),]
  ##### END of add reads mean information of lncRNA and co express genes ######
  
  ###### add circRNA position ########
  run_time_message("add circRNA position")
  linc_coexp_pairs$circRNA_chr<-sapply(strsplit(linc_coexp_pairs$circRNA,"_"),"[[",1)
  
  circstartend<-data.frame(V1=linc_coexp_pairs$circRNA_5exon<-sapply(strsplit(linc_coexp_pairs$circRNA,"_"),"[[",3),
             V2=linc_coexp_pairs$circRNA_5exon<-sapply(strsplit(linc_coexp_pairs$circRNA,"_"),"[[",2))
  linc_coexp_pairs$circRNA_5exon <-apply(circstartend,1,min)
  linc_coexp_pairs$circRNA_3exon <-apply(circstartend,1,max)
  linc_coexp_pairs$circRNA_strand<-ifelse(sapply(strsplit(linc_coexp_pairs$circRNA,"_"),"[[",4)=="fwd","+",
                                          ifelse(sapply(strsplit(linc_coexp_pairs$circRNA,"_"),"[[",4)=="rev","-","NA"))
  
  ######### add exon informaion #####################
  circ_bed<-unique(linc_coexp_pairs[,c("circRNA","circRNA_chr","circRNA_5exon","circRNA_3exon","circRNA_strand")])
  circ_bed$circRNA_5exon<-as.numeric(circ_bed$circRNA_5exon)
  circ_bed$circRNA_3exon<-as.numeric(circ_bed$circRNA_3exon)
  exonidx_df<-as.data.frame(data.table::fread(paste0("gzip -dc ",exonidx)))  #loading exon info db
  exonidx_df$V5<-gsub("\\.[0-9]*","",exonidx_df$V5) 
  exonidx_df$V6<-gsub("\\.[0-9]*","",exonidx_df$V6)
  colnames(exonidx_df)<-c("chr","start","end","strand","ensembl_geneid","ensembl_transciptid","gene_symbol","transcipt_symbol","exon_sum","exon_number")
  exonidx_df$exon_info<-paste0(exonidx_df$exon_number,"/",exonidx_df$exon_sum)
  
  new_linc_coexp_pairs<-data.frame()
  linc_coexp_pairs_sub<-data.frame()
  for (i in seq(1:nrow(circ_bed))){   
    ## extract 5exon info
    linc_coexp_pairs_sub<-linc_coexp_pairs[linc_coexp_pairs$circRNA %in% circ_bed[i,1],]
    tolerance<-2
    #find nearest exon
    #if close to start find the exons near to the start
    #if close to end find the exons near to the end
    # test if circ within the exon and rbind, unique would be better
    if(min(abs(circ_bed[i,3]-exonidx_df$start)) >= min(abs(exonidx_df$end-circ_bed[i,3]))){ 
      # creating subtracting list , find which start or end smaller or equal to the  min end or start
      # only find either start or end  
      exon5_ann<-exonidx_df[which(abs(exonidx_df$end-circ_bed[i,3])<= min(abs(exonidx_df$end-circ_bed[i,3]))+tolerance),] # 
      exon5_ann<-exon5_ann[exon5_ann$chr %in% circ_bed[i,2] & exon5_ann$strand %in% circ_bed[i,5],]
    } else{
      exon5_ann<-exonidx_df[which(abs(exonidx_df$start-circ_bed[i,3])<=min(abs(exonidx_df$start-circ_bed[i,3]))+tolerance),]
      exon5_ann<-exon5_ann[exon5_ann$chr %in% circ_bed[i,2] & exon5_ann$strand %in% circ_bed[i,5],]
    }
    ## extract 3exon info
    if(min(abs(circ_bed[i,4]-exonidx_df$start)) >= min(abs(exonidx_df$end-circ_bed[i,4]))){
      exon3_ann<-exonidx_df[which(abs(exonidx_df$end-circ_bed[i,4])<= min(abs(exonidx_df$end-circ_bed[i,4]))+tolerance),]
      exon3_ann<-exon3_ann[exon3_ann$chr %in% circ_bed[i,2] & exon3_ann$strand %in% circ_bed[i,5],]
    } else{
      exon3_ann<-exonidx_df[which(abs(exonidx_df$start-circ_bed[i,4])<=min(abs(exonidx_df$start-circ_bed[i,4]))+tolerance),]
      exon3_ann<-exon3_ann[exon3_ann$chr %in% circ_bed[i,2] & exon3_ann$strand %in% circ_bed[i,5],]
    }
    
    #ov_txnid<-intersect(exon5_ann$ensembl_transciptid ,exon3_ann$ensembl_transciptid )
    #cated_exonann<-rbind(exon3_ann,exon5_ann)
   # ov_geneid_paired<-cated_exonann[cated_exonann$ensembl_transciptid %in% ov_txnid,]$ensembl_geneid
    #unpaired_txnid<-cated_exonann[cated_exonann$ensembl_geneid %in% !ov_geneid,]$ensembl_transciptid
    #ov_txnid<-c(ov_txnid,unpaired_txnid)
    
   # if (length(ov_txnid) > 0){
   # exon3_ann<-exon3_ann[exon3_ann$ensembl_transciptid %in% ov_txnid, ] # 
   # exon5_ann<-exon5_ann[exon5_ann$ensembl_transciptid %in% ov_txnid, ] # 
   # } 
    
    colnames(exon5_ann)<-paste0("5exon_",colnames(exon5_ann))
    colnames(exon3_ann)<-paste0("3exon_",colnames(exon3_ann))
    #### exon5 #####
    if (nrow(exon5_ann) ==0 ) {exon5_ann[1,]<-"NA"}
    
    if (nrow(exon5_ann) ==1 ){
      exon5_ann<-cbind(circ_bed[rep(i, each=nrow(exon5_ann)),],exon5_ann)
      linc_coexp_pairs_sub<-merge(linc_coexp_pairs_sub,exon5_ann[,c(1,10,11,12,16)], all.x =T , by="circRNA" )
    } 
    if (nrow(exon5_ann) > 1 ){
      #exon5_ann[2,]<-exon5_ann[1,] for testing row more than2 
      exon5_ann<-cbind(circ_bed[rep(i, each=nrow(exon5_ann)),],exon5_ann)
      for (a in seq(2,ncol(exon5_ann))){
      exon5_ann[1,a]<-paste(as.character(exon5_ann[,a]),collapse=",",sep="")
      }
      linc_coexp_pairs_sub<-merge(linc_coexp_pairs_sub,exon5_ann[1,c(1,10,11,12,16)], all.x =T , by="circRNA" )
    }  
    #####  exon3 ######
    if (nrow(exon3_ann) ==0 ) {exon3_ann[1,]<-"NA"}
    if (nrow(exon3_ann) ==1 ){
      exon3_ann<-cbind(circ_bed[rep(i, each=nrow(exon3_ann)),],exon3_ann)
      linc_coexp_pairs_sub<-merge(linc_coexp_pairs_sub,exon3_ann[,c(1,10,11,12,16)], all.x =T , by="circRNA" )
    } 
    if (nrow(exon3_ann) > 1 ){
      #exon3_ann[2,]<-exon3_ann[1,] for testing row more than2 
      exon3_ann<-cbind(circ_bed[rep(i, each=nrow(exon3_ann)),],exon3_ann)
      for (a in seq(2,ncol(exon3_ann))){
        exon3_ann[1,a]<-paste(as.character(exon3_ann[,a]),collapse=",",sep="")
      }
      linc_coexp_pairs_sub<-merge(linc_coexp_pairs_sub,exon3_ann[1,c(1,10,11,12,16)], all.x =T , by="circRNA" )
    }  
    
    new_linc_coexp_pairs<-rbind(new_linc_coexp_pairs,linc_coexp_pairs_sub)
  }
  ######### add exon informaion END #####################
  
  
  # reorder column names
  linc_coexp_pairs<-new_linc_coexp_pairs[,c(1,20:ncol(new_linc_coexp_pairs),2:19)]
  
  ###### add circdb cirRNA id #############
#x<-c("chrY", "15435434" ,"15448215"    ,    "-" )
  #x<-as.character(circ_bed[1,2:5])
  circdb<-as.data.table(circdb)
  circrna_tolerance<-2
  search_circdb<-function(x){
    search_circdb_run<-circdb[chr == x[1] & strand== x[4]]
    search_circdb_run<-search_circdb_run[ !((start >  as.numeric(x[2]) &  start > as.numeric(x[3])) |  (end <  as.numeric(x[2]) &  end < as.numeric(x[3])))]
    search_circdb_run<-search_circdb_run[order(abs(as.numeric(x[2])-as.numeric(search_circdb_run$start))+abs(as.numeric(x[3])-as.numeric(search_circdb_run$end)))[1],]
    search_circdb_run<-search_circdb_run[!is.na(search_circdb_run$chr),]
    if (nrow( search_circdb_run) ==0 ) {
      search_circdb_run<-NULL
    } else if ((abs(search_circdb_run$start[1]-as.numeric(x[2])) > circrna_tolerance) | ( abs(search_circdb_run$end[1]-as.numeric(x[3])) > circrna_tolerance )){
       search_circdb_run<-NULL
    }
     paste0(search_circdb_run$circRNA_ID,collapse = ",")
  }
  
  ### for test ###
  circirna_test<-function(){  
  testt<-data.frame(circRNA=row.names(circrna))
  testt$chr<-sapply(strsplit(as.character(testt$circRNA),"_"),"[",1)
  testt$circRNA_exon1<-sapply(strsplit(as.character(testt$circRNA),"_"),"[",3)
  testt$circRNA_exon2<-sapply(strsplit(as.character(testt$circRNA),"_"),"[",2)
  testt$start <-apply(testt[,3:4],1,min)
  testt$end <-apply(testt[,3:4],1,max)
  testt$strand<- ifelse(sapply(strsplit(as.character(testt$circRNA),"_"),"[",4)=="fwd","+",
         ifelse(sapply(strsplit(as.character(testt$circRNA),"_"),"[",4)=="rev","-","NA"))
  testt<-testt[,c(1,2,5,6,7)]
  testt$idd<-apply(testt[,c(2,3,4,5)],1, try(search_circdb))
  nrow(testt[testt$idd=="",])
  }
  
  circ_bed$best_circRNA_id<-apply(circ_bed[,c(2,3,4,5)],1, try(search_circdb))
  
  linc_coexp_pairs$best_circRNA_id<-""
  for (i in seq(1,nrow(circ_bed))){
    linc_coexp_pairs[linc_coexp_pairs$circRNA %in% circ_bed[i,]$circRNA ,]$best_circRNA_id <- circ_bed[i,]$best_circRNA_id
  }
  
  linc_coexp_pairs<- linc_coexp_pairs[,c(1,ncol(linc_coexp_pairs),2:(ncol(linc_coexp_pairs)-1))]
  
  ###### add circdb cirRNA id END#############
  
  saveRDS(linc_coexp_pairs, paste0("output/linc_coexp_pairs.rds"))
  saveRDS(normreads, paste0("output/normreads.rds"))
  saveRDS(q_cor_summary, paste0("output/q_cor_summary.rds"))
  saveRDS(circrna, paste0("output/circrna_reads.rds"))
  write_feather(linc_coexp_pairs, paste0("output/linc_coexp_pairs.feather"))
  
  write.table(circrna,"output/circRNA.txt",sep = '\t',quote = F)
  write.table(linc_coexp_pairs,"output/linc_coexp_pairs.txt",sep = '\t',row.names = F,quote = F)
  write.table(q_cor_summary,"output/q_cor_summary.txt",sep = '\t',row.names = F,quote = F)
  
  #### save merged circrna gene table
  circ_gene_merged<-rbind(normreads,circrna[,colnames(normreads)])
  saveRDS(circ_gene_merged, paste0("output/circ_gene_merged_not.rds")) # no transposed
  circ_gene_merged<-as.data.frame(t(circ_gene_merged))
  circ_gene_merged<-merge(factor_list,circ_gene_merged, by.x= "V1" , by.y="row.names")
  row.names(circ_gene_merged)<-circ_gene_merged$V1
  circ_gene_merged$V1<-NULL
  colnames(circ_gene_merged)[1]<-"attr"
  saveRDS(circ_gene_merged, paste0("output/circ_gene_merged.rds"))
}

if( opt$mode[1]== "circ" ) {
  for_circ_mode() 
  
}else if ( opt$mode[1]== "lnc"){
  for_lnc_mode() 
  
 
}











