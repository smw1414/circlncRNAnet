#!/usr/bin/env Rscript
#######!/share/apps/R/bin/Rscript
args=commandArgs(TRUE)
library("getopt")


spec <- matrix(c(
  'query_gene'     , 'g', 1, "character","gene",
  'cutoff_pos'     , 'p', 1, "numeric"," positive correlation cutoff, -p 0.5 ",
  'cutoff_neg'     , 'n', 1, "numeric"," inverse correlation cutoff, -n -0.5 ",
  'cutoff_reg'     , 'r', 1, "character"," interval for excleude(ex) or include(in), -r ex, -c in "

),ncol=5,byrow=T)

opt = getopt(spec);

if (is.null(opt$query_gene  )) {
  cat(paste(getopt(spec, usage=T),"\nExample 
./heatmap.R -p 0.5 -n -0.5 -r ex -g chrX_47755339_47705503_fwd
./heatmap.R -p 0.5 -n -0.5 -r ex -g ELFN1-AS1 
./heatmap.R -p 0.5 -n -0.5 -r ex -g LINC00346 \n"));
  q();
}

library(pheatmap)
library(ggplot2)
library(data.table)

linc_coexp_pairs<-readRDS(paste0("output/linc_coexp_pairs.rds"))
df_opt<-readRDS(paste0("output/opt.rds"))
circ_gene_merged<-readRDS("output/circ_gene_merged.rds")
circ_gene_merged_t<-t(circ_gene_merged[,2:ncol(circ_gene_merged)])
# factor_list<-fread(df_opt["factorlist","V1"],header = F)
factor_list<-fread(paste0(getwd(),"/",df_opt["factorlist","V1"]),header = F)
factor_list[,V2:=as.factor(V2)]

query_gene<-opt$query_gene[1]
cutoff_pos<-opt$cutoff_pos[1]
cutoff_neg<-opt$cutoff_neg[1]
cutoff_reg<-opt$cutoff_reg[1]

asas<-function(){
  query_gene<<-"chr10_97437191_97438703_rev"
  cutoff_pos<<-0.5
  cutoff_neg<<-0.5
  cutoff_reg<<-"ex"
}

# asas()
# if ("gencodev19" %in% df_opt$V1 ){
#   gene_coordinate<-readRDS("v19_gene_coordinate.rds")
# } else if ("gencodev25" %in% df_opt$V1 ) {
#   gene_coordinate<-readRDS("v25_gene_coordinate.rds")
# } 

################################# filtering  #####################################
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
#filter list if gene lusist to large
gene_limit<-5000
 if (nrow(linc_coexp_pairs_filtered) > gene_limit ) {
   rank<-rank(abs(linc_coexp_pairs_filtered$cor))
   linc_coexp_pairs_filtered<-linc_coexp_pairs_filtered[which(rank> (nrow(linc_coexp_pairs_filtered)-gene_limit)),]

 }
print(nrow(linc_coexp_pairs_filtered
           ))
#######################  filtering END  #########################################
 
  
  circ_gene_merged_t<-as.data.table(melt(circ_gene_merged_t))
  circ_gene_merged_t<-as.data.frame(dcast(circ_gene_merged_t,Var1~Var2,sum))
  row.names(circ_gene_merged_t)<-circ_gene_merged_t$Var1
  circ_gene_merged_t$Var1<-NULL
  # order correlation by cor
  coexpgene_sym<-unique(linc_coexp_pairs_filtered[order(linc_coexp_pairs_filtered[,"cor"],decreasing = T),]$co_exp_gene)
  heatmaprow<-c(query_gene,  setdiff(coexpgene_sym, query_gene))
 
  tmp<-circ_gene_merged_t[heatmaprow,]
  sortedtab<-tmp
  sortedtab$samples<-rownames(sortedtab)
  sortedtab<-subset(sortedtab,select = samples)
  for ( factor_levels in 1:length(levels(factor_list[[2]]))){
    tmp2<-as.data.frame(tmp[heatmaprow,])
    tmp2<-tmp2[,colnames(tmp2) %in% as.character(factor_list[factor_list$V2 %in% levels(factor_list$V2)[factor_levels],1][[1]])]
    #for (revorder in length(rownames(tmp2)):1){
    #  tmp2<-tmp2[,order(tmp2[revorder,],decreasing = T)]
      #tmp2<-tmp2[,order(colnames(tmp2),decreasing = F)]
   # }
    sortedtab<-cbind(sortedtab,tmp2)
  }
  heatmapann<-as.data.frame(factor_list)
  rownames(heatmapann)<-heatmapann$V1
  heatmapann$V1<-NULL
  colnames(heatmapann)<-c("type")
  sortedtab$samples<-NULL
  # ###########################
  # # transform to max and min
  # sortedtab2<-log2(sortedtab)
  #data<-sortedtab2
  #  maxminnorm<-data.frame()
  #  for (row in 1:nrow(data)){
  #    for (col in 1:ncol(data)){
  #    maxminnorm[row,col]<-(10*(data[row,col]-min(na.omit(data[row,]))))/(max(na.omit(data[row,]))-min(na.omit(data[row,])))+(-5)
  #       #maxminnorm[row,col]<- (data[row,col]-median(as.numeric(data[row,]),na.rm=T))/mad(as.numeric(data[row,]),na.rm=T)
  #    }
  #  }
  #  row.names(maxminnorm)<-rownames(sortedtab)
  #  colnames(maxminnorm)<-colnames(sortedtab)
  #  sortedtab2<-maxminnorm
  ####################
  # annotation processing
  #roder factorlevel
  
  factor_list$V2<- factor(factor_list$V2, levels = c(as.character(factor_list[[2]][1]),
                                                     as.character(factor_list$V2[-grep(as.character(factor_list[[2]][1]),factor_list$V2)][1])))
  color_type = c(  N="blue",T="red")
  ann_colors = list(type= color_type)
  
  pdfwidth<-length(colnames(sortedtab))*0.05+5
  pdfwidth<-ifelse( pdfwidth > 45,45, pdfwidth)
  
  pdfheight <-(length(rownames(sortedtab))/length(colnames(sortedtab)))* pdfwidth*1
  pdfheight<-ifelse( pdfheight > 45,45, pdfheight)
  
  xx<-max(abs(min(scale(log2(t(sortedtab))),na.rm=T)),abs(max(scale(log2(t(sortedtab))),na.rm=T)))+0.5
  xx2<-t(scale((t(log2(sortedtab)))))
  qq<-max(abs(quantile(xx2, probs = 0.001,na.rm = T)),abs(quantile(xx2, probs = 0.999,na.rm = T)))
  maxx<-qq
  minn<--qq
  xx2[xx2 < minn] <- minn # replace 
  xx2[xx2 > maxx] <- maxx
  xx<-max(abs(min(xx2,na.rm=T)),abs(max(xx2,na.rm=T)))
  breaksList = seq(-xx-1, xx+1, by = 0.1)
  
  xx2<-melt(xx2)
  
  ggplot(xx2, aes(  Var2, Var1)) + 
    geom_raster(aes(fill = value)) + 
   scale_fill_gradientn(colours=c('blue', 'black', 'gold'),
                         name="Z-score")+
    theme_bw(base_size = 8)+
    coord_fixed(ratio=1)+
    theme(
      axis.text.x = element_text(angle =90, vjust = 0.5, 
                                 hjust = 1),
      axis.title.x = element_blank(),
      axis.title.y = element_blank(),
      panel.grid.major = element_blank(),
      panel.border = element_blank(),
      panel.background = element_blank(),
      axis.ticks = element_blank()
    )
  par(mar=c(0,0,0,0))
  ggsave(paste("output/heatmap",query_gene,cutoff_pos,cutoff_neg,cutoff_reg,".png",sep="_"), 
         width =pdfwidth,
         height = pdfheight)
  
  
  p<-ggplot(xx2, aes(  Var2, Var1)) + 
    geom_raster(aes(fill = value)) + 
    scale_fill_gradientn(colours=c('blue', 'black', 'gold'),
                         name="Z-score")+
    theme_bw(base_size = 8)+
    #coord_fixed(ratio=1)+
    theme (
      axis.text.x = element_text(
        angle = 90,
        vjust = 0.5,
        hjust = 1
      ) ,
      axis.title.x = element_blank(),
      axis.title.y = element_blank(),
      panel.grid.major = element_blank(),
      panel.border = element_blank(),
      panel.background = element_blank(),
      axis.ticks = element_blank()
    )
saveRDS(p,paste("output/heatmap",query_gene,cutoff_pos,cutoff_neg,cutoff_reg,".rds",sep="_"))
  
  
  
  


