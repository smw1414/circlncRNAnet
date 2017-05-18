#!/usr/bin/env Rscript

run_time_message<-function(msg){
  message(paste0(Sys.time()," : ",msg))
}


run_time_prom<-function(){
 h <- taskCallbackManager()
 h$add(function(expr, value, ok, visible) { 
       options("prompt"=format(Sys.time(), "%H:%M:%S > ")); 
               return(TRUE) }, 
         name = "simpleHandler")
}

trinetwork_loading_cor_res<-function(){
  library(feather)
  run_time_message("loading linc_coexp_pairs")
  # linc_coexp_pairs<<-readRDS(paste0("output/linc_coexp_pairs.rds"))
  linc_coexp_pairs<<-as.data.frame(read_feather(paste0("output/linc_coexp_pairs.feather")))
  df_opt<<-readRDS(paste0("output/opt.rds"))
  run_time_message("loading circ_gene_merged.rds")
  readsetable<-readRDS(paste0("output/circ_gene_merged_not.rds"))
  non_zero_gene<<-rownames(readsetable)
  
}

read_feather_dt<-function(x) {
  dt<-data.table::as.data.table(feather::read_feather(x))  
  return(dt)
}

#trinetwork_tabout(sponge_network,"sponge",query_symbol_input)

trinetwork_tabout<-function(pre_tab,words,query){
  run_time_message("writing table")
  parameters<-paste0(query,"_",cutoff_reg,"_",cutoff_pos,cutoff_neg,"_",DSD)
  base<-paste0("output/trinetwork_tab_",words,"_",parameters)
  
  co_expressed_along<-unique(pre_tab[!(query_symbol %in% query),!colnames(pre_tab)[ncol(pre_tab)-1:3],with=F])
  co_expressed_along<-merge(co_expressed_along,unique(linc_coexp_pairs_filtered[,c("co_exp_gene","cor")]),by.x = "query_symbol",by.y="co_exp_gene")
  
  write.table(pre_tab[,!colnames(pre_tab)[ncol(pre_tab)-1:3],with=F],paste0(base,"_query_co_expressed_merged.txt"),sep = '\t',quote = F,row.names = F)
  write.table(pre_tab[query_symbol %in% query,!colnames(pre_tab)[ncol(pre_tab)-1:3],with=F],paste0(base,"_query_along.txt"),sep = '\t',quote = F,row.names = F)
  write.table(co_expressed_along,paste0(base,"_co_expressed_along.txt"),sep = '\t',quote = F,row.names = F)
  # old
  write.table(pre_tab[query_symbol %in% query,!colnames(pre_tab)[ncol(pre_tab)-1:3],with=F],paste0("output/",words,"_co_expressed_merged.txt"),sep = '\t',quote = F,row.names = F)
  
}



trinetwork_nodelinkout<-function(words,query){
  run_time_message("writing node and link file level 0")
  parameters<-paste0(query,"_",cutoff_reg,"_",cutoff_pos,cutoff_neg,"_",DSD)
  
  out_txt_node=paste0("output/trinetwork_",words,"_node_",parameters,".txt")
  out_txt_link=paste0("output/trinetwork_",words,"_link_",parameters,".txt")
  
  write.table(testlink,out_txt_link,sep = "\t",quote=F,row.names = F)
  write.table(testnode,out_txt_node,sep = "\t",quote=F,row.names = F)
}



trinetwork_nodelinkout_level1<-function(words,query){
  run_time_message("writing node and link file level 1")
  parameters1<-paste0(query,"_",cutoff_reg,"_",cutoff_pos,cutoff_neg,"_",DSD)
  
  out_txt_nodel1=paste0("output/trinetwork_",words,"_node_level1_",parameters1,".txt")
  out_txt_linkl1=paste0("output/trinetwork_",words,"_link_level1_",parameters1,".txt")
  
  write.table(testlink_level1,out_txt_linkl1,sep = "\t",quote=F,row.names = F)
  write.table(testnode_level1,out_txt_nodel1,sep = "\t",quote=F,row.names = F)
  run_time_message("finished\n\n")
}



trinetwork_nodelinkout_level2<<-function(words,query){
  run_time_message("query and corresponding file found, writing node and link file level 2")
  parameters2<-paste0(query,"_",cutoff_reg,"_",cutoff_pos,cutoff_neg,"_",DSD,"_",acq)
  
  out_txt_nodel2=paste0("output/trinetwork_",words,"_node_level2_",parameters2,".txt")
  out_txt_linkl2=paste0("output/trinetwork_",words,"_link_level2_",parameters2,".txt")
  testnode_level2<-testnode_level2[order(group,decreasing = T),]
  
  write.table(testlink_level2,out_txt_linkl2,sep = "\t",quote=F,row.names = F)
  write.table(testnode_level2,out_txt_nodel2,sep = "\t",quote=F,row.names = F)
 
   run_time_message("finished\n\n")
}

read_tri_network_nodelink<-function(words,query){
  run_time_message(paste0("Analysis started\n",
                 "Query: ",query,"\n",
                 "cutoff pos: ",cutoff_pos,"\n",
                 "cutoff neg: ",cutoff_neg,"\n",
                 "cutoff reg: ",cutoff_reg,"\n",
                 "DSD: ",DSD,"\n",
                 "acq: ",acq,""))
  parameters<-paste0(query,"_",cutoff_reg,"_",cutoff_pos,cutoff_neg,"_",DSD)
  out_txt_node<<-paste0("output/trinetwork_",words,"_node_",parameters,".txt")
  out_txt_link<<-paste0("output/trinetwork_",words,"_link_",parameters,".txt")
}


creating_nodelink_RBP<-function(){
  run_time_message("generating node and link file ")
  ###  plot network
  #creating nodes
  testnode<-data.frame(label=unique(c(run_rbp_network[[1]])),stringsAsFactors=F) # 1=query_symbol , 3=RBP
  testnode$group<-ifelse( testnode$label %in% query_symbol_input,qry,coexp)
  testnode<-rbind(data.frame(label=unique(run_rbp_network[[3]])[unique(run_rbp_network[[3]]) %in% non_zero_gene ] #add RBP non-zero filter
                             ,group=ac,stringsAsFactors=F),testnode)
  #testnode<-testnode[order(testnode$label),]
  testnode<-cbind(data.frame(id=rep(1:nrow(testnode)),stringsAsFactors=F),testnode)
  testnode$value<-""
  testnode[testnode$group %in% qry,]$value<-5
  
  if (length(testnode[testnode$group %in% coexp,]$value)>0) {
    testnode[testnode$group %in% coexp,]$value<-2
  }
  if(length(testnode[testnode$group %in% ac,]$value)>0){
    testnode[testnode$group %in% ac,]$value<-3
  }
  
  testnode$title<-paste0(testnode$label) #add tooltip
  testnode<<-testnode
  # id=id , lable= gene_name , value=size of node , title=tooltip
  
  #creating links
  testlink<-run_rbp_network[,c(ac_col,"query_symbol","support_sources_count"),with=F]
  
  #merge RBP_id to link table
  testlink<-merge(testlink,testnode[testnode$group %in% ac,1:3],by.x=ac_col , by.y="label")
  colnames(testlink)[grep("^id$",colnames(testlink))]<-"from"  # from_id == RBP
  
  #merge co_express_gene_id to link table
  testlink<-merge(testlink,testnode[!testnode$group %in% ac,1:3],by.x="query_symbol", by.y="label" )
  colnames(testlink)[grep("^id$",colnames(testlink))]<-"to"  # to_id == lncRNA or co-expressed gene
  #testlink$label<-paste0(testlink$group.x,"_",testlink$group.y)
  testlink$title<-paste0(testlink[[ac_col]], " to ",testlink$query_symbol)
  testlink<- testlink[,c("from","to","support_sources_count","title")] #,"label"
  colnames(testlink)[grep("support_sources_count",colnames(testlink))]<-"color"
  testlink$color<-ifelse(testlink$color==3,"red","lightblue")
  testlink$color<-NULL
  testlink<<-testlink
}

creating_nodelink_sponge<-function(){
  run_time_message("generating node and link file ")
  ###  plot network 
  #creating nodes
  testnode<-data.frame(label=unique(as.character(sponge_network[[1]])),stringsAsFactors=F) # 1=query_symbol ,
  # assign gene belongs to query,co-express gene, miRNA 
  testnode$group<-ifelse( testnode$label %in% query_symbol_input,qry,coexp)
  
  #add miRNA non-zero filter, but no miRNA read count availavle
  testnode<-rbind(data.frame(label=unique(as.character(sponge_network[[3]])) 
                             ,group=ac,stringsAsFactors=F),testnode)
  
  # 
  testnode<-cbind(data.frame(id=rep(1:nrow(testnode)),stringsAsFactors=F),testnode)
  testnode$value<-""
  #assign node size value
  testnode[testnode$group %in% qry,]$value<-5
  if (length(testnode[testnode$group %in% coexp,]$value)>0) {
    testnode[testnode$group %in% coexp,]$value<-2
  }
  if(length(testnode[testnode$group %in% ac,]$value)>0){
    testnode[testnode$group %in% ac,]$value<-3
  }
  testnode$title<-paste0(testnode$label) #add tooltip
  testnode<<-testnode
  # id=id , lable= gene_name , value=size of node , title=tooltip
  
  #creating links
  testlink<-sponge_network[,c(ac_col,"query_symbol","support_sources_count"),with=F]
  
  #merge miRNA_id to link table
  testlink<-merge(testlink,testnode[testnode$group %in% ac,1:3],by.x=ac_col , by.y="label")
  colnames(testlink)[grep("^id$",colnames(testlink))]<-"from"  # from_id == miRNA
  
  #merge co_express_gene_id to link table
  testlink<-merge(testlink,testnode[!testnode$group %in% ac,1:3],by.x="query_symbol", by.y="label" )
  colnames(testlink)[grep("^id$",colnames(testlink))]<-"to"  # to_id == lncRNA or co-expressed gene
  #testlink$label<-paste0(testlink$group.x,"_",testlink$group.y)
  testlink$title<-paste0(testlink[[ac_col]], " to ",testlink$query_symbol)
  testlink<- testlink[,c("from","to","support_sources_count","title"),with=F] #,"label"
  colnames(testlink)[grep("support_sources_count",colnames(testlink))]<-"color"
  testlink$color<-ifelse(testlink$color==3,"red","lightblue")
  testlink$color<-NULL
  testlink<<-testlink
}

write_error_triple_tabnode0_fun<-function(words,query,comment){ # words= RBP, sponge
  # for level0 and table
  run_time_message(paste0("Error: ",comment))
  run_time_message("writing table")
  parameters<-paste0(query,"_",cutoff_reg,"_",cutoff_pos,cutoff_neg,"_",DSD)
  base<-paste0("output/trinetwork_tab_",words,"_",parameters)
  
  pre_tab<-data.frame(Error=paste0(comment))
  
  # tab out
  write.table(pre_tab,paste0(base,"_query_co_expressed_merged.txt"),sep = '\t',quote = F,row.names = F)
  write.table(pre_tab,paste0(base,"_query_along.txt"),sep = '\t',quote = F,row.names = F)
  write.table(pre_tab,paste0(base,"_co_expressed_along.txt"),sep = '\t',quote = F,row.names = F)
  write.table(pre_tab,paste0("output/",words,"_co_expressed_merged.txt"),sep = '\t',quote = F,row.names = F)
  
  # node level 0 error 
  node_df<-data.frame(label=c("No network result available",paste0(comment)),Error="Error")
  node_df<-cbind(data.frame(id=rep(1:nrow(node_df)),stringsAsFactors=F),node_df) # generate ID
  link_df<-data.frame(from=1,to=2)
  
  out_txt_node=paste0("output/trinetwork_",words,"_node_",parameters,".txt")
  out_txt_link=paste0("output/trinetwork_",words,"_link_",parameters,".txt")
  
  write.table(link_df,out_txt_link,sep = "\t",quote=F,row.names = F)
  write.table(node_df,out_txt_node,sep = "\t",quote=F,row.names = F)
  
  # node level 1 error 
  run_time_message("writing node and link file level 1")
  parameters1<-paste0(query,"_",cutoff_reg,"_",cutoff_pos,cutoff_neg,"_",DSD)
  
  out_txt_nodel1=paste0("output/trinetwork_",words,"_node_level1_",parameters1,".txt")
  out_txt_linkl1=paste0("output/trinetwork_",words,"_link_level1_",parameters1,".txt")
  
  write.table(link_df,out_txt_linkl1,sep = "\t",quote=F,row.names = F)
  write.table(node_df,out_txt_nodel1,sep = "\t",quote=F,row.names = F)
  
  # node level 2 error
  run_time_message("query and corresponding file found, writing node and link file level 2")
  parameters2<-paste0(query,"_",cutoff_reg,"_",cutoff_pos,cutoff_neg,"_",DSD,"_",acq)
  
  out_txt_nodel2=paste0("output/trinetwork_",words,"_node_level2_",parameters2,".txt")
  out_txt_linkl2=paste0("output/trinetwork_",words,"_link_level2_",parameters2,".txt")
  
  write.table(link_df,out_txt_linkl2,sep = "\t",quote=F,row.names = F)
  write.table(node_df,out_txt_nodel2,sep = "\t",quote=F,row.names = F)
  
}


write_error_triple_tabnode2_fun<-function(words,query,comment){ # words= RBP, sponge
  # for level0 and table
  run_time_message(paste0("Error: ",comment))
  run_time_message("writing table")

  # node level 2 error
  run_time_message("No co-expressed genes, writing node and link file level 2")
  parameters2<-paste0(query,"_",cutoff_reg,"_",cutoff_pos,cutoff_neg,"_",DSD,"_",acq)
  
  node_df<-data.frame(label=c("No network result available",paste0(comment)),Error="Error")
  node_df<-cbind(data.frame(id=rep(1:nrow(node_df)),stringsAsFactors=F),node_df) # generate ID
  link_df<-data.frame(from=1,to=2)
  
  out_txt_nodel2=paste0("output/trinetwork_",words,"_node_level2_",parameters2,".txt")
  out_txt_linkl2=paste0("output/trinetwork_",words,"_link_level2_",parameters2,".txt")
  
  write.table(link_df,out_txt_linkl2,sep = "\t",quote=F,row.names = F)
  write.table(node_df,out_txt_nodel2,sep = "\t",quote=F,row.names = F)
  
}





