#!/usr/bin/env Rscript
args=commandArgs(TRUE)
library("getopt")

spec <- matrix(c(
  'query'     , 'q', 1, "character", "Query genes, separate by comma (required)",
  'cutoff_pos'     , 'p', 1, "numeric"," positive correlation cutoff, -p 0.5 ",
  'cutoff_neg'     , 'n', 1, "numeric"," inverse correlation cutoff, -n -0.5 ",
  'cutoff_reg'     , 'r', 1, "character"," interval for excleude(ex) or include(in), -r ex, -c in ",
  'display_supported_data'     , 'd', 1, "numeric"," number of supported data source to display  -d 2 -d 3 ",
  'RBP'     , 'b', 1, "character"," RBP -b HNRNPK "
),ncol=5,byrow=T)

opt = getopt(spec);

if ( is.null(opt$query  )) {
  cat(paste(getopt(spec, usage=T),"\nExample (linux): ./triple_network.R -p 0.5 -n -0.5 -r ex -q gene_name -d 2 -b RBP_name \n"));
  q();
}
source("universal_function.R")

gencode_files<-list(gencodev19=c("v19_gene_coordinate.rds",
                                 "ENCODE_eCLIP_all_RBP_V19.feather",
                                 "fimo_all_v19.feather",
                                 "all_hg19_pipseq.txt"),
                    gencodev25=c("v25_gene_coordinate.rds",
                                 "ENCODE_eCLIP_all_RBP_V25.feather",
                                 "fimo_all_v25.feather",
                                 "all_hg38_pipseq.txt")
)

library(data.table)
library(plyr)
#library(ggplot2)
library(visNetwork)

query_symbol_input<-opt$query[1]
cutoff_pos<-opt$cutoff_pos[1]
cutoff_neg<-opt$cutoff_neg[1]
cutoff_reg<-opt$cutoff_reg[1]
DSD<-opt$display_supported_data[1]
acq<-opt$RBP[1] # asscociated component query

asas<-function(){
  query_gene<<-"AC021218.2"
  query_symbol_input<<-"AC021218.2"
  cutoff_pos<<-0.5
  cutoff_neg<<--0.5
  cutoff_reg<<-"in"
  acq<-"PTBP1"
}

qry<-"query"
ac_col<-"RBP"
ac<-"RBP"
coexp<-"co-expressed_gene"
# visnetwork assign the color by the order of vector,  

read_tri_network_nodelink("RBP",query_symbol_input)

if(!is.null(acq) & (file.exists(out_txt_node) & file.exists(out_txt_link) )){
  
  testnode<-fread(out_txt_node)
  testlink<-fread(out_txt_link)

  #level2 link
  ids<-testnode[label %in% c(acq,query_symbol_input)]$id
  testlink_level2<-unique(testlink[ to  %in% ids | from %in% ids])
  
  #level2 node
  testnode_level2 <-
    unique( testnode[testnode$id %in% c(testlink_level2$to, testlink_level2$from) , ])
  
  # new added
  ##########################################################################################
  co_exp_lim<-200
  if (!(coexp %in% unique(testnode_level2$group))) {
    write_error_triple_tabnode2_fun("RBP",query_symbol_input,"No co-expressed gene sharing the same RBP")
    q()
  }
  if (nrow(testnode_level2[group %in% coexp ]) > co_exp_lim ) {
    write_error_triple_tabnode2_fun("RBP",query_symbol_input,paste0("Number of co-expressed genes > ",co_exp_lim,", no co-expressed network will be shown for the concern of browser performance "))
    q()
  }
  
  ###########################################################################################
  # end of new added
  
  trinetwork_nodelinkout_level2("RBP",query_symbol_input)
  
  run_time_message("Terminated by reading node and link file")
  q()
}

# load correlation result files 
trinetwork_loading_cor_res()
run_time_message("filter co-express genes")
################################# getting co-expressed gene  #####################################

rna_idx<-ifelse( length(which(colnames(linc_coexp_pairs) %in% "lncRNA")) == 1,which(colnames(linc_coexp_pairs) %in% "lncRNA"),
                 ifelse(length(which( colnames(linc_coexp_pairs) %in% "circRNA")) == 1,which( colnames(linc_coexp_pairs) %in% "circRNA"),"NA"))

rnacolumn<-linc_coexp_pairs[,colnames(linc_coexp_pairs)[rna_idx]]

linc_coexp_pairs_filtered<-linc_coexp_pairs[ rnacolumn %in%  query_symbol_input , ]

linc_coexp_pairs_filtered<-linc_coexp_pairs_filtered[!is.na(linc_coexp_pairs_filtered$cor),]

if( cutoff_reg == "ex") {
  linc_coexp_pairs_filtered<-linc_coexp_pairs_filtered[which(linc_coexp_pairs_filtered$cor > cutoff_pos |  linc_coexp_pairs_filtered$cor < cutoff_neg), ]
} 
if(cutoff_reg == "in" ) { 
  linc_coexp_pairs_filtered<-linc_coexp_pairs_filtered[which(!(linc_coexp_pairs_filtered$cor > cutoff_pos &  linc_coexp_pairs_filtered$cor < cutoff_neg)), ]
}

################################################################################## 
run_time_message(paste0("number of coexpressed gene: ", 
             nrow(linc_coexp_pairs_filtered[! linc_coexp_pairs_filtered[,"co_exp_gene"]  %in%  query_symbol_input , ])))

# if (nrow(linc_coexp_pairs_filtered[! linc_coexp_pairs_filtered[,"co_exp_gene"]  %in%  query_symbol_input , ])==0) {
#   system('echo RBP binding site not found > output/RBP_netowrk_error.txt')
#   q("Please lower the cutoff for co-expression network")
#   }

# write_error_triple_tabnode0_fun(word,comment,query)
# write_error_triple_tabnode2_fun(word,comment,query)

query_symbol_input<-query_symbol_input
gene_coordinate<-readRDS(gencode_files[[df_opt["annotation","V1"]]][1])
run_time_message(paste0("loading PIP-seq data "))
pipseq<-fread(gencode_files[[df_opt["annotation","V1"]]][4])
run_time_message("loading ENCODE data")
exp_lncRNA_RBP<-read_feather_dt(gencode_files[[df_opt["annotation","V1"]]][2])
run_time_message("loading FIMO data ")
com_lncRNA_RBP<-read_feather_dt(gencode_files[[df_opt["annotation","V1"]]][3])

run_time_message("filter RBP db ")
com_lncRNA_RBP<-com_lncRNA_RBP[query_symbol %in% c(query_symbol_input,linc_coexp_pairs_filtered$co_exp_gene) ]
exp_lncRNA_RBP<-exp_lncRNA_RBP[query_symbol %in% c(query_symbol_input,linc_coexp_pairs_filtered$co_exp_gene) ]
pipseq<-pipseq[query_symbol %in% c(query_symbol_input,linc_coexp_pairs_filtered$co_exp_gene)]

run_rbp<-unique(rbind(com_lncRNA_RBP[ ,1:4,with=F],exp_lncRNA_RBP[ ,1:4,with=F]))

notin<-unique( c(which(pipseq$query_id %in% setdiff(unique(pipseq$query_id),unique(run_rbp$query_id ))), #extract RBP only in PIPSEQ
                 which(pipseq$query_symbol %in% setdiff(unique(pipseq$query_symbol),unique(run_rbp$query_symbol )) )))
pipseq2<-pipseq[ notin,]
pipseq2$pvalue<-""
pipseq2<-pipseq2[,c(1,2,3,4,5,7,6),with=F]
run_rbp<-rbind(run_rbp,pipseq2[,1:4,with=F])
run_rbp<-run_rbp[query_symbol %in% unique(c(query_symbol_input,linc_coexp_pairs_filtered$co_exp_gene)),1:4,with=F]
run_rbp<-run_rbp[run_rbp$RBP %in% run_rbp[run_rbp$query_symbol %in% query_symbol_input, ]$RBP,]  # filter query lncRNA  associated components
run_rbp<-as.data.table(data.frame(lapply(run_rbp, as.character), stringsAsFactors=FALSE)) # conver all coulmn to factor
####
run_time_message(paste0("lncRNA RBP dependent lncRNA/co-expressed gene-associated compment link count: ",nrow(run_rbp)))

#output error if no RBP site found 
if(length(run_rbp[run_rbp$query_symbol %in% query_symbol_input, ]$RBP)==0) {
  write_error_triple_tabnode0_fun("RBP",query_symbol_input,"RBP binding site not found")
  q("RBP binding site not found")
  }

CISBP<-com_lncRNA_RBP[grep("CISBP",com_lncRNA_RBP$data)]
CISBP$CISBP<-paste0(CISBP$count,"(",round(CISBP$score,2),")")
CISBP<-CISBP[,c( "query_symbol","query_id","RBP","RBP_id","CISBP"),with=F]

Ray2013<-com_lncRNA_RBP[grep("Ray2013",com_lncRNA_RBP$data)]
Ray2013$Ray2013<-paste0(Ray2013$count,"(",round(Ray2013$score,2),")")
Ray2013<-Ray2013[,c( "query_symbol","query_id","RBP","RBP_id","Ray2013"),with=F]

exp_lncRNA_RBP$cellline<-gsub(".*\\(|\\)","",exp_lncRNA_RBP$data)

K562<-exp_lncRNA_RBP[exp_lncRNA_RBP$cellline %in% "ENCODE_eCLIP_K562"]
K562$K562<-paste0(K562$count,"(",round(K562$score,2),")")
K562<-K562[,c( "query_symbol","query_id","RBP","RBP_id","K562"),with=F]

HepG2<-exp_lncRNA_RBP[exp_lncRNA_RBP$cellline %in% "ENCODE_eCLIP_HepG2"]
HepG2$HepG2<-paste0(HepG2$count,"(",round(HepG2$score,2),")")
HepG2<-HepG2[,c( "query_symbol","query_id","RBP","RBP_id","HepG2"),with=F]

pipseq$PIPseq<-paste0(pipseq$count,"(",round(pipseq$score,2),")")
pipseq<- pipseq[,c( "query_symbol","query_id","PIPseq"),with=F]

run_rbp<-merge(run_rbp,CISBP,by=c("query_symbol","query_id","RBP","RBP_id" ),all.x=T)
run_rbp<-merge(run_rbp,Ray2013,by=c("query_symbol","query_id","RBP","RBP_id" ),all.x=T)
run_rbp<-merge(run_rbp,K562,by=c("query_symbol","query_id","RBP","RBP_id" ),all.x=T)
run_rbp<-merge(run_rbp,HepG2,by=c("query_symbol","query_id","RBP","RBP_id" ),all.x=T)
run_rbp<-merge(run_rbp,pipseq,by=c("query_symbol","query_id" ),all.x=T)

run_rbp[is.na(run_rbp)]<-""

run_rbp$com_count<-ifelse(run_rbp$CISBP != "" | run_rbp$Ray2013 != "",1,0)
run_rbp$exp_count<-ifelse(run_rbp$K562 != "" | run_rbp$HepG2 != "",1,0)
run_rbp$pipseq_count<-ifelse(run_rbp$PIPseq != "",1,0)

#### generate df for table output
run_rbp$support_sources_count<-apply(run_rbp[,c("com_count","exp_count","pipseq_count")],1,sum)
run_rbp<-run_rbp[run_rbp$RBP %in% run_rbp[run_rbp$query_symbol %in% query_symbol_input, ]$RBP,]  # filter  mutual query lncRNA  associated components
# trinetwork_tabout(run_rbp,"RBP",query_symbol_input)

#### generate df for network use
run_rbp_network<-run_rbp[run_rbp$support_sources_count >= DSD,]
run_rbp_network<-run_rbp_network[run_rbp_network$RBP %in% run_rbp_network[run_rbp_network$query_symbol %in% query_symbol_input, ]$RBP,] # filter  mutual query lncRNA  associated components


 #####parse nod and link########
run_time_message(paste0("RBP count: ", length(run_rbp_network[run_rbp_network$query_symbol %in% query_symbol_input, ]$RBP)))
run_time_message(paste0("query and co-expressed gene: ",length(unique(run_rbp_network$query_symbol))))
run_time_message(paste0("lncRNA RBP dependent lncRNA/co-expressed gene-associated compment link count with 2 or more supported data: ",nrow(run_rbp_network)))


#output error if no RBP site found  with DSD setting
if(length( unique(run_rbp_network[ run_rbp_network$query_symbol %in% query_symbol_input, ]$RBP))==0){
  #system('echo RBP binding site not found, after DSD filter > output/RBP_netowrk_error.txt' )
  write_error_triple_tabnode0_fun("RBP",query_symbol_input,"No RBP binding site found after supported source filter")
  q("RBP binding site not found, after DSD filter")
  }
#test if no co-expressed gene after the mutual associated component filter
# if(length( unique(run_rbp_network[ !run_rbp_network$query_symbol %in% query_symbol_input, ]$query_symbol))==0)
# {system('echo "No co-expressed genes after mutual associated component filter" > output/RBP_netowrk_error.txt' );q("No co-expressed genes after mutual associated component filter")}

trinetwork_tabout(run_rbp_network,"RBP",query_symbol_input)

#if(nrow(run_rbp_network) > 2000) {
#  write.table(data.frame(Error="too many network links, please lower the cutoff for co-expression net"),
#              paste0("output/enrichment_res_",query_symbol_input,".txt"),sep = '\t',row.names = F,quote = F)
#  q("too many network links")
# } 

#group name table
groupnames<-data.frame(query="lncRNA",coexp="co-expressed_gene",associated_compoment="RBP",stringsAsFactors=F)


###  plot network
#creating nodes
# testnode<-data.frame(label=unique(c(run_rbp_network[[1]])),stringsAsFactors=F) # 1=query_symbol , 3=RBP
# testnode$group<-ifelse( testnode$label %in% query_symbol_input,groupnames$query[1],groupnames$coexp[1])
# testnode<-rbind(data.frame(label=unique(run_rbp_network[[3]])[unique(run_rbp_network[[3]]) %in% non_zero_gene ] #add RBP non-zero filter
#                            ,group=groupnames$associated_compoment[1],stringsAsFactors=F),testnode)
# #testnode<-testnode[order(testnode$label),]
# testnode<-cbind(data.frame(id=rep(1:nrow(testnode)),stringsAsFactors=F),testnode)
# testnode$value<-""
# testnode[testnode$group %in% groupnames$query[1],]$value<-5
# testnode[testnode$group %in% groupnames$coexp[1],]$value<-2
# testnode[testnode$group %in% groupnames$associated_compoment[1],]$value<-3
# testnode$title<-paste0(testnode$label) #add tooltip
# # id=id , lable= gene_name , value=size of node , title=tooltip
# 
# #creating links
# testlink<-run_rbp_network[,c(groupnames$associated_compoment[1],"query_symbol","support_sources_count"),with=F]
# 
# #merge RBP_id to link table
# testlink<-merge(testlink,testnode[testnode$group %in% groupnames$associated_compoment[1],1:3],by.x=groupnames$associated_compoment[1] , by.y="label")
# colnames(testlink)[grep("^id$",colnames(testlink))]<-"from"  # from_id == RBP
# 
# #merge co_express_gene_id to link table
# testlink<-merge(testlink,testnode[!testnode$group %in% groupnames$associated_compoment[1],1:3],by.x="query_symbol", by.y="label" )
# colnames(testlink)[grep("^id$",colnames(testlink))]<-"to"  # to_id == lncRNA or co-expressed gene
# #testlink$label<-paste0(testlink$group.x,"_",testlink$group.y)
# testlink$title<-paste0(testlink[,groupnames$associated_compoment[1]]," to ",testlink$query_symbol)
# testlink<- testlink[,c("from","to","support_sources_count","title")] #,"label"
# colnames(testlink)[grep("support_sources_count",colnames(testlink))]<-"color"
# testlink$color<-ifelse(testlink$color==3,"red","lightblue")
# testlink$color<-NULL

# creating network
creating_nodelink_RBP()

# write out all node and link
trinetwork_nodelinkout( "RBP",query_symbol_input)

#level1 node
testnode_level1<-testnode[testnode$group %in% c(ac, qry),]
#level1 link
testlink_level1<-testlink[testlink$to %in% testnode_level1[testnode_level1$label %in% query_symbol_input ,]$id,]

trinetwork_nodelinkout_level1("RBP",query_symbol_input)

asdas<-function(){
  
  visNetwork(testnode_level1, testlink_level1, height = "1000px", width = "100%") %>%
    visOptions(selectedBy = "group", 
               highlightNearest = TRUE, 
               nodesIdSelection = TRUE) %>%
    visPhysics(maxVelocity=50,stabilization = T) %>%
    visLayout(randomSeed = 123) %>%
    visEdges(smooth =  list(enabled = TRUE, type = 'continuous',roundness=0),physics=F)# %>%
  #visIgraphLayout(randomSeed=123)
  
  visNetwork(testnode_level2, testlink_level2, height = "1000px", width = "100%") %>%
    visOptions(selectedBy = "group", 
               highlightNearest = TRUE, 
               nodesIdSelection = TRUE) %>%
    visPhysics(maxVelocity=50,stabilization = T) %>%
    visLayout(randomSeed = 123) %>%
    visEdges(smooth =  list(enabled = TRUE, type = 'continuous',roundness=0),physics=F)# %>%
  
  
  pipseq<-"all_hg38_pipseq.txt"
  encode<-"ENCODE_eCLIP_all_RBP_V25.rds"
  fimo<-"fimo_all_v25.rds"
  query_symbol_input<-"SNHG15"
  query_symbol_input<-"CCAT1" 
  #query_symbol_input<-"CCAT1"

  #acq<-"FMR1"
  cutoff_pos<-0.7
  cutoff_neg<--0.7
  cutoff_reg<-"ex"
   DSD<-2
}














