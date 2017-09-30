#!/usr/bin/env Rscript

############################################################################
# circRNA-RBP triple network                                               #
# identify the circRNA/co-expressed genes RBP binding and plot the network #
############################################################################

args=commandArgs(TRUE)
library("getopt")


spec <- matrix(c(
  'query'     , 'q', 1, "character", "Query genes, separate by comma (required)",
  'cutoff_pos'     , 'p', 1, "numeric"," positive correlation cutoff, -p 0.5 ",
  'cutoff_neg'     , 'n', 1, "numeric"," inverse correlation cutoff, -n -0.5 ",
  'cutoff_reg'     , 'r', 1, "character"," interval for excleude(ex) or include(in), -r ex, -c in ",
  'display_supported_data'     , 'd', 1, "numeric"," number of supported data source to display  -d 2 -d 3 ",
  'RBP'     , 'b', 1, "character"," RBP -b HNRNPK "
  # 'color_supported_data'     , 'd', 1, "character"," number of supported data to be colored  -d 2 -d 3 "
),ncol=5,byrow=T)

opt = getopt(spec);

if ( is.null(opt$query  )) {
  cat(paste(getopt(spec, usage=T),"\nExample (linux): ./triple_network.R -p 0.95 -n -0.95 -r ex -q chr11_35204640_35201082_fwd -d 2 
Example (linux): ./triple_network.R -p 0.95 -n -0.95 -r ex -q chr11_35204640_35201082_fwd -d 2 -b KHDRBS2\n"));
  q();
}

source("universal_function.R")

gencode_files<-list(gencodev19=c("v19_gene_coordinate.rds",
                                 "ENCODE_eCLIP_all_RBP_V19.feather",
                                 "ENCODE_eCLIP_circ_RBP_V19.feather",
                                 "fimo_all_v19.feather",
                                 "fimo_v19_circ.feather",
                                 "all_hg19_pipseq.txt",
                                 "cirrnadb_hg19_pipseq.txt"),
                    gencodev25=c("v25_gene_coordinate.rds",
                                 "ENCODE_eCLIP_all_RBP_V25.feather",
                                 "ENCODE_eCLIP_circ_RBP_V25.feather",
                                 "fimo_all_v25.feather",
                                 "fimo_v25_circ.feather",
                                 "all_hg38_pipseq.txt",
                                 "cirrnadb_hg38_pipseq.txt")
)

#for (i in c("ENCODE_eCLIP_all_RBP_V19.rds","ENCODE_eCLIP_circ_RBP_V19.rds","ENCODE_eCLIP_all_RBP_V25.rds","ENCODE_eCLIP_circ_RBP_V25.rds")){
#  write_feather(readRDS(i),paste0(gsub("rds$","",i),"feather"))
#}


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

# assign column name
qry<-"query"
ac_col<-"RBP"
ac<-"RBP"
coexp<-"co-expressed_gene"

# read the path of linknode file 
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
  if (!("co-expressed_gene" %in% unique(testnode_level2$group))) {
    write_error_triple_tabnode2_fun("RBP",query_symbol_input,"No co-expressed gene sharing the same RBP")
    q()
  }
  if (nrow(testnode_level2[group %in% "co-expressed_gene"]) > co_exp_lim ) {
    write_error_triple_tabnode2_fun("RBP",query_symbol_input,paste0("Number of co-expressed genes > ",co_exp_lim,", no co-expressed network will be shown for the concern of browser performance "))
    q()
  }
  
  ###########################################################################################
  # end of new added
  
  trinetwork_nodelinkout_level2("RBP",query_symbol_input)
  
  run_time_message("Terminated by reading node and link file")
  q()
}


#loading files
trinetwork_loading_cor_res()
run_time_message("filter co-express genes")
################################# getting co-expressed gene  #####################################
# identify whether the linc_coexp_pairs  are circRNA mode 

  query_id_input<-unique(linc_coexp_pairs[linc_coexp_pairs$circRNA %in% query_symbol_input,]$best_circRNA_id)
  # query_id_input circdb_id 
  # query_symbol_input circ query symbol e.g. chr12_68836749_68828771_fwd

if(as.numeric(nchar(query_id_input)) < 1) {system('echo circBase ID not found > output/sponge_netowrk_error.txt' );q("circBase ID not found")}

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
run_time_message(paste0("number of coexpressed gene: ",nrow(linc_coexp_pairs_filtered)))

# if (nrow(linc_coexp_pairs_filtered)==0) {
#   system('echo RBP binding site not found > output/RBP_netowrk_error.txt')
#   q("Please lower the cutoff for co-expression network")
# }

#merge_com_exp_rbp_queries_circ<-function(encode,encode_circ,fimo,fimo_circ,pipseq,pipseq_circ){

# loading ENCODE FIMO and PIPseq data
gene_coordinate<-readRDS(gencode_files[[df_opt["annotation","V1"]]][1])
run_time_message(paste0("loading PIP-seq data "))
pipseq<-unique(rbind(fread(gencode_files[[df_opt["annotation","V1"]]][6]),fread(gencode_files[[df_opt["annotation","V1"]]][7])))
run_time_message("loading ENCODE data")
exp_lncRNA_RBP<-rbind(read_feather_dt(gencode_files[[df_opt["annotation","V1"]]][2]),read_feather_dt(gencode_files[[df_opt["annotation","V1"]]][3]))
run_time_message("loading FIMO data ")
com_lncRNA_RBP<-rbind(read_feather_dt(gencode_files[[df_opt["annotation","V1"]]][4]),read_feather_dt(gencode_files[[df_opt["annotation","V1"]]][5]))
run_time_message("filter RBP db ")
com_lncRNA_RBP<-com_lncRNA_RBP[query_symbol %in% c(query_id_input,linc_coexp_pairs_filtered$co_exp_gene) ]
exp_lncRNA_RBP<-exp_lncRNA_RBP[query_symbol %in% c(query_id_input,linc_coexp_pairs_filtered$co_exp_gene) ]
pipseq<-pipseq[query_symbol %in% c(query_id_input,linc_coexp_pairs_filtered$co_exp_gene)]

# formated data, ready for filter
run_rbp<-unique(rbind(com_lncRNA_RBP[ ,1:4,with=F],exp_lncRNA_RBP[ ,1:4,with=F]))
 
# extract RBP only in PIPSEQ
notin<-unique( c(which(pipseq$query_id %in% setdiff(unique(pipseq$query_id),unique(run_rbp$query_id ))),
                 which(pipseq$query_symbol %in% setdiff(unique(pipseq$query_symbol),unique(run_rbp$query_symbol )) )))
pipseq2<-pipseq[ notin,]
pipseq2$pvalue<-""
pipseq2<-pipseq2[,c(1,2,3,4,5,7,6),with=F]

# combine with pipseq data
run_rbp<-rbind(run_rbp,pipseq2[,1:4,with=F])

# filter the table with query and co-expressed gene
run_rbp<-run_rbp[query_symbol %in% unique(c(query_id_input,linc_coexp_pairs_filtered$co_exp_gene)),1:4,with=F]

# filter query circRNA  associated components
run_rbp<-run_rbp[run_rbp$RBP %in% run_rbp[run_rbp$query_symbol %in% query_id_input, ]$RBP,]  

# conver all coulmn to factor
run_rbp<-data.table(data.frame(lapply(run_rbp, as.character), stringsAsFactors=FALSE))


run_time_message(paste0("circRNA RBP dependent lncRNA/co-expressed gene-associated compment link count: ",nrow(run_rbp)))
#output error if no RBP site found 
if(length(run_rbp[run_rbp$query_symbol %in% query_id_input, ]$RBP)==0) {
  write_error_triple_tabnode0_fun("RBP",query_symbol_input,"RBP binding site not found")
  system('echo RBP binding site not found > output/RBP_netowrk_error.txt' );q("RBP binding site not found")
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
run_rbp<-as.data.table(run_rbp[run_rbp$RBP %in% run_rbp[run_rbp$query_symbol %in% query_symbol, ]$RBP,] ) # filter mutual query circRNA  associated components
run_rbp[run_rbp$query_id %in% query_id_input,]$query_symbol<-query_symbol_input #replace circdbid to circsymbol
# trinetwork_tabout(run_rbp,"RBP",query_symbol_input)

#### generate df for network use
run_rbp_network<-run_rbp[run_rbp$support_sources_count >= DSD,]
run_rbp_network<-as.data.table(run_rbp_network[run_rbp_network$RBP %in% run_rbp_network[run_rbp_network$query_symbol %in% query_symbol_input, ]$RBP,] ) # filter mutual query circRNA  associated components


#####parse nod and link #####
run_time_message(paste0("RBP count: ", length(run_rbp_network[run_rbp_network$query_symbol %in% query_symbol_input, ]$RBP)))
run_time_message(paste0("query and co-expressed gene: ",length(unique(run_rbp_network$query_symbol))))
run_time_message(paste0("circRNA RBP dependent circRNA/co-expressed gene-associated compment link count with 2 or more supported data: ",nrow(run_rbp_network)))

#output error if no RBP site found  with DSD setting
if(length( run_rbp_network[ run_rbp_network$query_symbol %in% query_symbol_input, ]$RBP)==0){
  write_error_triple_tabnode0_fun("RBP",query_symbol_input,"No RBP binding site found after supported source filter")
  system('echo RBP binding site not found, after DSD filter > output/RBP_netowrk_error.txt' );q("RBP binding site not found, after DSD filter")
  }
#test if no co-expressed gene after the mutual associated component filter
# if(length( unique(run_rbp_network[ !run_rbp_network$query_symbol %in% query_symbol_input, ]$query_symbol))==0)
# {system('echo "No co-expressed genes after mutual associated component filter" > output/RBP_netowrk_error.txt' );q("No co-expressed genes after mutual associated component filter")}
trinetwork_tabout(run_rbp_network,"RBP",query_symbol_input)

#group name table
groupnames<-data.frame(query="query",coexp="co-expressed_gene",associated_compoment="RBP",stringsAsFactors=F)


# creating network
creating_nodelink_RBP()

# write out all node and link
trinetwork_nodelinkout( "RBP",query_symbol_input)


#level1 node
testnode_level1<-testnode[testnode$group %in% c(ac, qry),]
#level1 link
testlink_level1<-testlink[testlink$to %in% testnode_level1[testnode_level1$label %in% query_symbol_input ,]$id,]

trinetwork_nodelinkout_level1("RBP",query_symbol_input)













