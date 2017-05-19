#!/usr/bin/env Rscript
args=commandArgs(TRUE)
library("getopt")


spec <- matrix(c(
  'query'     , 'q', 1, "character", "Query genes, separate by comma (required)",
  'cutoff_pos'     , 'p', 1, "numeric"," positive correlation cutoff, -p 0.5 ",
  'cutoff_neg'     , 'n', 1, "numeric"," inverse correlation cutoff, -n -0.5 ",
  'cutoff_reg'     , 'r', 1, "character"," interval for excleude(ex) or include(in), -r ex, -c in ",
  'display_supported_data'     , 'd', 1, "numeric"," number of supported data source to display  -d 2 -d 3 ",
  'miRNA'     , 'b', 1, "character"," miRNA -b hsa-mir-365b-3p "
  #'color_supported_data'     , 'd', 1, "character"," number of supported data to be colored  -d 2 -d 3 "
),ncol=5,byrow=T)

opt = getopt(spec);

if ( is.null(opt$query  )) {
  cat(paste(getopt(spec, usage=T),"\nExample (linux): ./triple_network.R -p 0.5 -n -0.5 -r ex -q gene_name -d 2 -b miRNA_id \n"));
  q();
}

source("universal_function.R")

# load correlation result files 


gencode_files<-list(gencodev19=c("v19_gene_coordinate.rds",
                                 "all_v19_sponge.feather",
                                 "circrna_hg19_sponge.feather"),
                    gencodev25=c("v25_gene_coordinate.rds",
                                 "all_v25_sponge.feather",
                                 "circrna_hg38_sponge.feather")
)

library(data.table)
library(plyr)
#library(ggplot2)
library(visNetwork)

query_symbol_input <- opt$query[1]
cutoff_pos <- opt$cutoff_pos[1]
cutoff_neg <- opt$cutoff_neg[1]
cutoff_reg <- opt$cutoff_reg[1]
DSD <- opt$display_supported_data[1]
acq <- opt$miRNA[1] # asscociated component query


qry<-"query"
ac<-"miRNA"
ac_col<-"miRNA"
coexp<-"co-expressed_gene"


# read the path of linknode file 
read_tri_network_nodelink("sponge",query_symbol_input)
run_time_message("filter co-express genes")
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
    write_error_triple_tabnode2_fun("sponge",query_symbol_input,"No co-expressed gene sharing the same miRNA targetting site")
    q()
  }
  if (nrow(testnode_level2[group %in% "co-expressed_gene" ]) > co_exp_lim) {
    write_error_triple_tabnode2_fun("sponge",query_symbol_input,paste0("Number of co-expressed genes > ",co_exp_lim,", no co-expressed network will be shown for the concern of browser performance"))
    q()
  }
  ###########################################################################################
  # end of new added
  
  trinetwork_nodelinkout_level2("sponge",query_symbol_input)

  run_time_message("Terminated by reading node and link file")
  q()
  }

# loading files 
trinetwork_loading_cor_res()

################################# getting co-expressed gene  #####################################
# identify whether the linc_coexp_pairs  are circRNA mode 

 query_id_input<-unique(linc_coexp_pairs[linc_coexp_pairs$circRNA %in% query_symbol_input,]$best_circRNA_id)
# query_sid circdb_id 
# symbol_qqgene circ coexp symbol e.g. chr12_68836749_68828771_fwd

# output error if no circBase ID found
if(as.numeric(nchar(query_id_input)) < 1) {
  write_error_triple_tabnode0_fun("sponge",query_symbol_input,"circBase ID not found")
  #system('echo circBase ID not found > output/sponge_netowrk_error.txt' );q("circBase ID not found")
  }

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

# 
# if (nrow(linc_coexp_pairs_filtered)==0) {
#   system('echo sponge not found > output/sponge_netowrk_error.txt')
#   q("Please lower the cutoff for co-expression network")
# }

run_time_message("loading spong db")
gene_coordinate<-readRDS(gencode_files[[df_opt["annotation","V1"]]][1])
sponge_co<-read_feather_dt(gencode_files[[df_opt["annotation","V1"]]][2])
sponge_q<-read_feather_dt(gencode_files[[df_opt["annotation","V1"]]][3])

run_time_message("filter co-expressed gene miRNA targetting site")
sponge_q<-unique(sponge_q[sponge_q$query_id %in%  query_id_input ])
sponge_co<-sponge_co[sponge_co$miRNA %in% sponge_q[sponge_q$query_id %in%  query_id_input, ]$miRNA,]
sponge_co<-unique(sponge_co[query_symbol %in% c(linc_coexp_pairs_filtered$co_exp_gene) ])
sponge<-rbind(sponge_q,sponge_co)

sponge$data2<-paste0(sponge$count,"(",round(sponge$score,2),")")

sponge[query_id %in%  query_id_input ]$query_symbol<-query_symbol_input # replace id to symbol
#miranda<-sponge[software=="miranda"]
#RNAHybrid<-sponge[software=="RNAHybrid"]
#TarpMir<-sponge[software=="TarpMir"]
run_time_message("spliting data subset")
for (i in c("miranda","RNAHybrid","TarpMir")) {
  colnames(sponge)[8]<-i
  assign(i,sponge[software==i])
}

# output error if no ncRNA sponge detected
# if(nrow(sponge[query_symbol %in% query_symbol_input,])==0) {
#   run_time_message("number of circRNA miRNA targetting site: ",paste0(nrow(sponge[query_symbol %in% query_symbol_input,])))
#   system('echo sponge not found > output/sponge_netowrk_error.txt' )
#   q()}

# output error if no ncRNA sponge detected
if(nrow(sponge[query_symbol %in% query_symbol_input,])==0) {
  write_error_triple_tabnode0_fun("sponge",query_symbol_input,"miRNA targeting site not found")
  q("miRNA targeting site not found")}

run_time_message(paste0("circRNA sponge dependent circRNA/co-expressed gene-associated compment link count: ",nrow(sponge)))

#output error if no sponge site found 
if(length(sponge[sponge$query_symbol %in% query_symbol_input, ]$miRNA)==0) {system('echo sponge not found > output/sponge_netowrk_error.txt' );q("spoge not found")}

run_time_message("obtain unique mirna")
sponge<-unique(sponge[,1:3 ,with=F]) 

sponge_list<-list(miranda=miranda,RNAHybrid=RNAHybrid,TarpMir=TarpMir)

run_time_message("merge three software ")
#merge three software 
for (i in seq(1:length(sponge_list))){
  sponge<-merge(sponge,sponge_list[[i]][,c(1,2,3,8),with=F], #col= c(query_symbol,query_id,miRNA,software_col)
                by=c("query_symbol","query_id","miRNA" ),all.x=T)
}

# count the for each software
sponge[is.na(sponge)]<-""

for (i in c("miranda","RNAHybrid","TarpMir")) {
  sponge[,paste0(i,"_count")]<-ifelse(sponge[[i]] != "",1,0)
}

#### generate df for table output
sponge$support_sources_count<-apply(sponge[,c("miranda_count","RNAHybrid_count","TarpMir_count")],1,sum)
sponge<-as.data.table(sponge[sponge$miRNA %in% sponge[sponge$query_id %in%  query_id_input, ]$miRNA,]) 
# trinetwork_tabout(sponge, "sponge", query_symbol_input)

sponge_network<-sponge[sponge$support_sources_count >= DSD,]
sponge_network<-as.data.table(sponge_network[sponge_network$miRNA %in% sponge_network[sponge_network$query_id %in%  query_id_input, ]$miRNA,])  # filter  mutual co-exp circRNA  associated components


#####parse nod and link########
run_time_message(paste0("miRNA count: ", length(sponge_network[sponge_network$query_symbol %in% query_symbol_input, ]$miRNA)))
run_time_message(paste0("query and co-expressed gene count: ",length(unique(sponge_network$query_symbol))))
run_time_message(paste0("circRNA sponge dependent circRNA/co-expressed gene-associated compment link count with 2 or more supported data: ",nrow(sponge_network)))

#output error if no sponge site found  with DSD setting
# if(length( unique(sponge_network[ sponge_network$query_id %in%  query_id_input, ]$miRNA))==0)
# {system('echo sponge not found, after DSD filter > output/sponge_netowrk_error.txt' );q("miRNA sponge not found, after DSD filter")}

if(length( unique(sponge_network[ sponge_network$query_id %in% query_id_input, ]$miRNA))==0){
  #system('echo sponge not found, after DSD filter > output/sponge_netowrk_error.txt' );q("miRNA sponge not found, after DSD filter")
  write_error_triple_tabnode0_fun("sponge",query_symbol_input,"No miRNA targeting site found after supported source filter")
  q("No miRNA targeting site found after supported source filter")
}

#test if no co-expressed gene after the mutual associated component filter
# if(length( unique(sponge_network[ !sponge_network$query_symbol %in% query_id_input, ]$query_symbol))==0)
# {system('echo "No co-expressed genes after mutual associated component filter" > output/sponge_netowrk_error.txt' );q("No co-expressed genes after mutual associated component filter")}
trinetwork_tabout(sponge_network, "sponge", query_symbol_input)

#group name table
groupnames<-data.frame(query="query",coexp="co-expressed_gene",associated_compoment="miRNA",stringsAsFactors=F)

###  plot network 
#creating nodes
# testnode<-data.frame(label=unique(as.character(sponge_network[[1]])),stringsAsFactors=F) # 1=query_symbol ,
# # assign gene belongs to query,co-express gene, miRNA 
# testnode$group<-ifelse( testnode$label %in% query_symbol_input,groupnames$query[1],groupnames$coexp[1])
# testnode<-rbind(data.frame(label=unique(as.character(sponge_network[[3]])) 
#                            ,group=groupnames$associated_compoment[1],stringsAsFactors=F),testnode)
# 
# testnode<-cbind(data.frame(id=rep(1:nrow(testnode)),stringsAsFactors=F),testnode)
# testnode$value<-""
# #assign node size value
# testnode[testnode$group %in% groupnames$query[1],]$value<-5
# try(testnode[testnode$group %in% groupnames$coexp[1],]$value<-2)
# testnode[testnode$group %in% groupnames$associated_compoment[1],]$value<-3
# testnode$title<-paste0(testnode$label) #add tooltip
# # id=id , lable= gene_name , value=size of node , title=tooltip
# 
# # creating links
# testlink<-sponge_network[,c(groupnames$associated_compoment[1],"query_symbol","support_sources_count"),with=F]
# 
# # merge miRNA_id to link table
# testlink<-merge(testlink,testnode[testnode$group %in% groupnames$associated_compoment[1],1:3],by.x=associated_compoment[1] , by.y="label")
# colnames(testlink)[grep("^id$",colnames(testlink))]<-"from"  # from_id == miRNA
# 
# # merge co_express_gene_id to link table
# testlink<-merge(testlink,testnode[!testnode$group %in% groupnames$associated_compoment[1],1:3],by.x="query_symbol", by.y="label" )
# colnames(testlink)[grep("^id$",colnames(testlink))]<-"to"  # to_id == circRNA or co-expressed gene
# #testlink$label<-paste0(testlink$group.x,"_",testlink$group.y)
# testlink$title<-paste0(testlink$miRNA," to ",testlink$query_symbol)
# testlink<- testlink[,c("from","to","support_sources_count","title")] #,"label"
# colnames(testlink)[grep("support_sources_count",colnames(testlink))]<-"color"
# testlink$color<-ifelse(testlink$color==3,"red","lightblue")
# testlink$color<-NULL

# creating node and link 
creating_nodelink_sponge()

# write out all node and link
trinetwork_nodelinkout( "sponge",query_symbol_input)


# level1 node
testnode_level1<-testnode[testnode$group %in% c(ac, qry),]
# level1 link
testlink_level1<-testlink[testlink$to %in% testnode_level1[testnode_level1$label %in% query_symbol_input ,]$id,]


trinetwork_nodelinkout_level1("sponge",query_symbol_input )



asdas<-function(){
  
  visNetwork(testnode_level1, testlink_level1, height = "1000px", width = "100%") %>%
    visOptions(selectedBy = "group", 
               highlightNearest = TRUE, 
               nodesIdSelection = TRUE) %>%
    visPhysics(maxVelocity=50,stabilization = T) %>%
    visLayout(randomSeed = 123) %>%
    visEdges(smooth =  list(enabled = TRUE, type = 'continuous',roundness=0),physics=F)# %>%
  #visIgraphLayout(randomSeed=123)
  
  visNetwork(testnode_level2, testlink_level2, height = "1000px", width = "150%") %>%
    visOptions(selectedBy = "group", 
               highlightNearest = TRUE, 
               nodesIdSelection = TRUE) %>%
    visPhysics(maxVelocity=50,stabilization = T) %>%
    visLayout(randomSeed = 123) %>%
    visEdges(smooth =  list(enabled = TRUE, type = 'continuous',roundness=0),physics=F)# %>%
  
  #query_symbol_input<-"chr11_35204640_35201082_fwd"
  acq<-"hsa-mir-1260a"
  query_symbol_input<-"chr20_17955365_17957037_rev"
  query_symbol_input<-"chr10_97437191_97438703_rev"
 
  #acq<-"hsa-mir-574-5p"
  cutoff_pos<-0.5
  cutoff_neg<--0.5
  cutoff_reg<-"ex"
  DSD<-2
}














