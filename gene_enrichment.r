#!/usr/bin/env Rscript
args=commandArgs(TRUE)
library("getopt")


spec <- matrix(c(
  'query_gene'     , 'g', 1, "character","gene",
  'cutoff_pos'     , 'p', 1, "numeric"," positive correlation cutoff, -p 0.5 ",
  'cutoff_neg'     , 'n', 1, "numeric"," inverse correlation cutoff, -n -0.5 ",
  'cutoff_reg'     , 'r', 1, "character"," interval for excleude(ex) or include(in), -r ex, -c in ",
  'enrich_type'     , 'e', 1, "character"," enrichment geneset , -e bp ,GO BP. -e cc ,GO cc. -e mf ,GO mf.  -e kegg ,kegg.  -e hm, msigdb hallmark. -e tf, transcription factor. -e encodetf, ENCODE transcription factor   "
),ncol=5,byrow=T)

opt = getopt(spec);

if ( is.null(opt$enrich_type) |is.null(opt$query_gene  )) {
  cat(paste(getopt(spec, usage=T),"\nExample (linux): ./gene_enrichment.r -p 0.5 -n -0.5 -r ex -g chrX_47755339_47705503_fwd -e bp \n"));
  q();
}

#https://www.r-bloggers.com/kegg-enrichment-analysis-with-latest-online-data-using-clusterprofiler/
#library(goseq)
library("clusterProfiler")
library(org.Hs.eg.db)
library(ggplot2)
library(visNetwork)
source("universal_function.R")

query_gene<-opt$query_gene[1]
cutoff_pos<-opt$cutoff_pos[1]
cutoff_neg<-opt$cutoff_neg[1]
cutoff_reg<-opt$cutoff_reg[1]

enrich_type<-opt$enrich_type[1]


library(feather)
run_time_message("loading linc_coexp_pairs")
linc_coexp_pairs<-as.data.frame(read_feather(paste0("output/linc_coexp_pairs.feather")))
df_opt<-readRDS(paste0("output/opt.rds"))

if ("gencodev19" %in% df_opt["annotation",1] ){
  gene_coordinate<-readRDS("v19_gene_coordinate.rds")
} else if ("gencodev25" %in% df_opt["annotation",1] ) {
  gene_coordinate<-readRDS("v25_gene_coordinate.rds")
}

db_path<-list(gencodev19=data.frame(msigdb_hm="msigdb/h.all.v5.2.entrez.gmt",
                                    msigdb_tf="msigdb/c3.tft.v5.2.entrez.gmt",
                                    encode_tf="enchriment_ENCODE_hg19_ChIP_3000_1000_q0.01.feather"),
              gencodev25=data.frame(msigdb_hm="msigdb/h.all.v5.2.entrez.gmt",
                                    msigdb_tf="msigdb/c3.tft.v5.2.entrez.gmt",
                                    encode_tf="enchriment_ENCODE_hg38_ChIP_3000_1000_q0.01.feather")
)

msigdb_hm<-as.character(db_path[[df_opt["annotation",1]]][1,"msigdb_hm"])
msigdb_tf<-as.character(db_path[[df_opt["annotation",1]]][1,"msigdb_tf"])
encode_tf<-as.character(db_path[[df_opt["annotation",1]]][1,"encode_tf"])

run_time_message(paste0("Analysis started\n",
                        "Query: ",query_gene,"\n",
                        "cutoff pos: ",cutoff_pos,"\n",
                        "cutoff neg: ",cutoff_neg,"\n",
                        "cutoff reg: ",cutoff_reg,"\n",
                        "enrich_type: ",enrich_type,"\n"
                        ))

################################################################################
rna_idx<-ifelse( length(which("lncRNA" %in% colnames(linc_coexp_pairs))) == 1,which("lncRNA" %in% colnames(linc_coexp_pairs)),
                 ifelse(length(which("circRNA" %in% colnames(linc_coexp_pairs))) == 1,which("circRNA" %in% colnames(linc_coexp_pairs)),"NA"))
if (rna_idx =="NA"){q();}
rnacolumn<-linc_coexp_pairs[,colnames(linc_coexp_pairs)[rna_idx]]

linc_coexp_pairs_filtered<-linc_coexp_pairs[ rnacolumn %in%  query_gene , ]

linc_coexp_pairs_filtered<-linc_coexp_pairs_filtered[!is.na(linc_coexp_pairs_filtered$cor),]

if( cutoff_reg == "ex") {
  linc_coexp_pairs_filtered<-linc_coexp_pairs_filtered[which(linc_coexp_pairs_filtered$cor > cutoff_pos |  linc_coexp_pairs_filtered$cor < cutoff_neg), ]
} 
if(cutoff_reg == "in" ) { 
  linc_coexp_pairs_filtered<-linc_coexp_pairs_filtered[which(!(linc_coexp_pairs_filtered$cor > cutoff_pos &  linc_coexp_pairs_filtered$cor < cutoff_neg)), ]
}


################################################################################## 
write_error_fun<-function(comment){ 
  run_time_message(paste0("Error: ",comment))
  # error tab
  write.table(data.frame(Error=paste0(comment)),
              paste0("output/enrichment_res_",enrich_type,"_",query_gene,"_",cutoff_reg,"_",cutoff_pos,cutoff_neg,"_.txt")
              ,sep = '\t',row.names = F,quote = F)
  # bar plot
  error_df<-data.frame(Description=paste0(comment),pvalue=0.05)
  error_df$Description<-gsub(",",",\n",error_df$Description)
  width=4.5
  height=2
  ggplot(error_df,  aes(x=Description, y= -log10(pvalue)) )+ 
    geom_bar(stat="identity") +
    theme_bw(base_size = 11) +
    coord_flip()+
    xlab("")+
    ylab("pvalue (-log10)")+
    #coord_fixed(ra1)+
    scale_x_discrete(limits=paste(error_df$Description[length(error_df$Description):1]))+  # order by original df
    theme(
      aspect.ratio = 2,
      panel.grid.minor = element_blank(),
      plot.title = element_text(hjust = 1)
    )  + 
    scale_y_continuous( expand = c(0.02,0.02) )+
    ggtitle("No enrichment result available")
  ggsave(paste0("output/enrichment_res_",enrich_type,"_",query_gene,"_",cutoff_reg,"_",cutoff_pos,cutoff_neg,"_.png"),height = height,width =width,dpi=300 )
  
  # error node and link
  node_df<-data.frame(label=c("No network result available",paste0(comment)))
  node_df<-cbind(data.frame(id=rep(1:nrow(node_df)),stringsAsFactors=F),node_df) # generate ID
  link_df<-data.frame(from=1,to=2)
  write.table(link_df,paste0("output/pathway_network_link_",enrich_type,"_",query_gene,"_",cutoff_reg,"_",cutoff_pos,cutoff_neg,".txt"),sep = "\t",quote=F,row.names = F)
  write.table(node_df,paste0("output/pathway_network_node_",enrich_type,"_",query_gene,"_",cutoff_reg,"_",cutoff_pos,cutoff_neg,".txt"),sep = "\t",quote=F,row.names = F)
  
}
gene_symbol<-linc_coexp_pairs_filtered$co_exp_gene
run_time_message(paste0("number of gene before entrenz id conversion: ", length(gene_symbol[!gene_symbol %in% query_gene])))
if(length(gene_symbol[!gene_symbol %in% query_gene])<=0) {write_error_fun("No co-expressed gene, please lower the correlation cutoff");q()  } 

gene_coordinate_noNA<-gene_coordinate[!is.na(gene_coordinate$entrezgene),]

# entrenz ID conversion
gene.df1 <- try(bitr(unique(linc_coexp_pairs_filtered$co_exp_gene), fromType = "SYMBOL", toType = c("ENTREZID", "ENSEMBL"),OrgDb = "org.Hs.eg.db"))
gene.df3 <- try(bitr(unique(linc_coexp_pairs$co_exp_gene), fromType="SYMBOL", toType= c("ENTREZID", "ENSEMBL"), OrgDb="org.Hs.eg.db"))

if(class(gene.df1)=="try-error") { write_error_fun("entrenz ID conversion error, please contact authors");q()  } 

gene<-gene.df1$ENTREZID

# 
gene_bd<-unique(Reduce(union, list(gene.df3$ENTREZID)))
gene_bd_symbol<-linc_coexp_pairs$co_exp_gene

run_time_message(paste0("number of gene after entrenz id conversion: ", length(gene)))
if(length(gene)<=0) {write_error_fun("No co-expressed gene after entrenz ID conversion, please lower the correlation cutoff");q()  } 

run_kegg<-function(){
  enrich_res <- enrichKEGG(gene, organism="hsa",
                           pvalueCutoff = 0.05,
                           keyType = "kegg",
                           minGSSize=5,
                           maxGSSize=2000,
                           pAdjustMethod="BH",
                           universe=gene_bd)
  enrich_res<-as.data.frame(enrich_res) 
  return(enrich_res)
}

run_mf<-function(){
  enrich_res<-enrichGO(gene,
                       'org.Hs.eg.db',
                       pvalueCutoff = 0.05,
                       ont = "MF", 
                       pAdjustMethod = "BH",
                       minGSSize=5,
                       maxGSSize=2000,
                       universe=gene_bd  )
  enrich_res<-as.data.frame(enrich_res) 
  return(enrich_res)
}
run_bp<-function(){
  enrich_res<-enrichGO(gene,
                       'org.Hs.eg.db',
                       ont = "BP", 
                       minGSSize=5,
                       maxGSSize=2000,
                       pvalueCutoff = 0.05,
                       pAdjustMethod = "BH",
                       universe=gene_bd  )
  enrich_res<-as.data.frame(enrich_res) 
  return(enrich_res)
}
run_cc<-function(){
  enrich_res<-enrichGO(gene,
                       'org.Hs.eg.db',
                       ont = "CC", 
                       minGSSize=5,
                       maxGSSize=2000,
                       pvalueCutoff = 0.05,
                       pAdjustMethod = "BH",
                       universe=gene_bd  )
  #head(summary(goCC))
  enrich_res<-as.data.frame(enrich_res) 
  return(enrich_res)
}

run_hm<-function(){
  h <- read.gmt(msigdb_hm)
  enrich_res<- enricher(gene, TERM2GENE=h,
                        universe=gene_bd,
                        pvalueCutoff = 0.05,
                        minGSSize=5,
                        maxGSSize=2000,
                        pAdjustMethod = "BH")
  enrich_res<-as.data.frame(enrich_res) 
  enrich_res$Description<-gsub("HALLMARK_","",enrich_res$Description)
  enrich_res$Description<-gsub("_"," ",enrich_res$Description)
  return(enrich_res)
}


run_tf<-function(){
  h <- read.gmt(msigdb_tf)
  enrich_res<- enricher(gene, TERM2GENE=h,
                        universe=gene_bd,
                        minGSSize=5,
                        maxGSSize=2000,
                        pvalueCutoff = 0.05,
                        pAdjustMethod = "BH")
  enrich_res<-as.data.frame(enrich_res) 
  return(enrich_res)
}

run_encodetf<-function(){
  library(feather)
  
  encode_chipp<-as.data.frame(read_feather(encode_tf))
  encode_chipp<-encode_chipp[,c("ont","gene")]

  enrich_res<- enricher(gene_symbol, TERM2GENE=encode_chipp,
                        universe=gene_bd_symbol, 
                        minGSSize=5,
                        maxGSSize=40000,
                        pvalueCutoff = 0.05,
                        pAdjustMethod = "BH")
  enrich_res<-as.data.frame(enrich_res) 
  return(enrich_res)
}

if (enrich_type == "kegg") {
  enrich_res <- try(run_kegg())
}
if (enrich_type == "cc") {
  enrich_res <- try(run_cc())
}
if (enrich_type == "mf") {
  enrich_res <- try(run_mf())
}
if (enrich_type == "bp") {
  enrich_res <- try(run_bp())
}
if (enrich_type == "hm") {
  enrich_res <- try(run_hm())
}
if (enrich_type == "tf") {
  enrich_res <- try(run_tf())
}
if (enrich_type == "encodetf") {
  enrich_res <- try(run_encodetf())
}
if(class(enrich_res)=="try-error") { write_error_fun(paste0("Enrichment error, number of input gene after entrenz ID conversion: ",length(gene)));q()  } 


egmt_df<-as.data.frame(enrich_res) 
egmt_df<-egmt_df[!is.na(egmt_df$pvalue),]
run_time_message(paste0("number of enriched terms: ", nrow(egmt_df)))

if(is.null(nrow(egmt_df)) )  { write_error_fun(paste0("Enrichment error, number of input gene after entrenz ID conversion: ",length(gene)));q()  } 
if(nrow(egmt_df)<=0)   { write_error_fun(paste0("No term was enriched, number of input gene after entrenz ID conversion: ",length(gene)));q()  } 


if (enrich_type!= "encodetf"){
  term_id_dup<-data.frame() # conver ID to symbol
  for (j in seq(1:nrow(egmt_df[,]))){ 
    term_id_dup_tmp<-egmt_df[rep(j, each=length(strsplit(egmt_df[j,]$geneID,"/")[[1]])),]
    term_id_dup_tmp$geneID<-strsplit(as.character(egmt_df[j,]$geneID),"/")[[1]]
    symbol_entrenz<-gene.df1[gene.df1$ENTREZID %in% unique(term_id_dup_tmp$geneID),]
    term_id_dup_tmp<-merge(term_id_dup_tmp,symbol_entrenz,by.x= "geneID", by.y= "ENTREZID",all.x=T)
    egmt_df[j,]$geneID<-paste(unique(term_id_dup_tmp$SYMBOL),collapse = "/", sep = "")
    term_id_dup<-rbind(term_id_dup,term_id_dup_tmp)
  }
}


egmt_df<-data.frame(egmt_df[order(egmt_df$pvalue,decreasing = F),])
egmt_df<-egmt_df[,c(2,5,6,7,8,9)]
colnames(egmt_df)[5]<-"gene"
write.table(egmt_df,
            paste0("output/enrichment_res_",enrich_type,"_",query_gene,"_",cutoff_reg,"_",cutoff_pos,cutoff_neg,"_.txt"),sep = '\t',row.names = F,quote = F)
run_time_message("writing enrichment table")

# remove pvalue=0
egmt_df<-egmt_df[egmt_df$pvalue !=0,]

if(nrow(egmt_df)> 15 ){
  egmt_df_15<-egmt_df[1:15,]
} else {
  egmt_df_15<-egmt_df
}
if(nrow(egmt_df)<=0)   { write_error_fun(paste0("No term was enriched, number of input gene after entrenz ID conversion: ",length(gene)));q()  } 
# break the terms that have too long length 
# for (stt in seq(1,nrow(egmt_df_15))) {  
#   x<-gregexpr(" ",egmt_df_15$Description[stt])[[1]][7]
#   if (!is.na(x)) {substring(egmt_df_15$Description[stt], x, x) <-"\n"  
#   } else {next}
# }

width=(max(nchar(egmt_df_15$Description))/16)+4
height=(nrow(egmt_df_15)/10)+3
ggplot(egmt_df_15,  aes(x=Description, y= -log10(pvalue)) )+ 
  geom_bar(stat="identity") +
  theme_bw(base_size = 11) +
  coord_flip()+
  xlab("")+
  ylab("pvalue (-log10)")+
  #coord_fixed(ra1)+
  scale_x_discrete(limits=paste(egmt_df_15$Description[length(egmt_df_15$Description):1]))+  # order by original df
  theme(
    #axis.text.x = element_text(size = 18),
    #axis.text.y = element_text(size =13),
    #legend.key = element_rect(fill = "navy"),
    #plot.title = element_text(size=18, hjust =0),
    aspect.ratio = 2,
    #legend.background = element_rect(fill = "white"),
    panel.grid.minor = element_blank()
    #axis.text.x = element_text(angle = 0, hjust = 1, vjust =1),
    #aspect.ratio = 2,
    #panel.grid.major.x = element_blank()+
    #plot.margin = unit( c(0,0,0,0) , "in" )
  )  + 
  scale_y_continuous( expand = c(0.02,0.02) )
ggsave(paste0("output/enrichment_res_",enrich_type,"_",query_gene,"_",cutoff_reg,"_",cutoff_pos,cutoff_neg,"_.png"),height = height,width =width,dpi=300 )
run_time_message("output barplot, finished\n")

# plot network
#prep node 
maxrow<-ifelse(nrow(egmt_df) >=3,3,nrow(egmt_df))
egmt_df_maxrow<-egmt_df[ 1:maxrow ,]

node_df<-data.frame(label=unique(c(egmt_df_maxrow$Description,unique(unlist(strsplit(egmt_df_maxrow$gene,"/"))))))
node_df$group<-ifelse( node_df$label %in% egmt_df_maxrow$Description,"terms","genes") 
node_df<-cbind(data.frame(id=rep(1:nrow(node_df)),stringsAsFactors=F),node_df) # generate ID
node_df$value<-""
node_df[node_df$group %in% "terms",]$value<-5
node_df[node_df$group %in% "genes",]$value<-2
node_df$title<-paste0(node_df$label) #add tooltip

#prep link 
link_df<-data.table::as.data.table(egmt_df_maxrow ,key="Description")
link_df<- link_df[, list(gene = unlist(strsplit(gene, "/"))), by=Description]
names(link_df)<-paste0("link_",names(link_df))
#merge  link table
link_df<-merge(link_df,node_df[node_df$group %in% "terms",1:2] ,by.x="link_Description" , by.y="label")
colnames(link_df)[grep("^id$",colnames(link_df))]<-"from"  # from_id == link_Description

#merge co_express_gene_id to link table
link_df<-merge(link_df,node_df[!node_df$group %in% "gene",1:2],by.x="link_gene", by.y="label" )
colnames(link_df)[grep("^id$",colnames(link_df))]<-"to"

write.table(link_df,paste0("output/pathway_network_link_",enrich_type,"_",query_gene,"_",cutoff_reg,"_",cutoff_pos,cutoff_neg,".txt"),sep = "\t",quote=F,row.names = F)
write.table(node_df,paste0("output/pathway_network_node_",enrich_type,"_",query_gene,"_",cutoff_reg,"_",cutoff_pos,cutoff_neg,".txt"),sep = "\t",quote=F,row.names = F)

visNetwork(node_df, link_df, height = "1000px", width = "100%") %>%
   visOptions(selectedBy = "group", 
              highlightNearest = TRUE, 
              nodesIdSelection = TRUE) %>%
  visPhysics(maxVelocity=1,stabilization = T) %>%
  visLayout(randomSeed = 123) %>%
 visEdges(smooth =  list(enabled = TRUE, type = 'continuous',roundness=0),physics=F)

asdasdfas<-function(){

  
  #enrich_type<-"hm"
  #pathway_network_node_hm_chr17_69128619_69134742_rev_ex_0.5-0.5
  
  query_gene<-"CCAT1"
 # query_gene<-"chr10_115129535_115120185_fwd"
  cutoff_pos<-0.15
  cutoff_neg<--0.15
  cutoff_reg<-"ex"
  
  enrich_type<-"kegg"
}


