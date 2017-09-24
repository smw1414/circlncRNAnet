#!/usr/bin/env Rscript
######!/share/apps/R/bin/Rscript
###### ! /usr/bin/env Rscript
args=commandArgs(TRUE)
library("getopt")

spec <- matrix(c(
  'query_gene', 'g', 1, "character", "gene for randomization of pearson correlation"#,
  # 'cutoff_pos'     , 'p', 1, "numeric"," positive correlation cutoff, -p 0.5 ",
  # 'cutoff_neg'     , 'n', 1, "numeric"," inverse correlation cutoff, -n -0.5 ",
  # 'cutoff_reg'     , 'r', 1, "character"," interval for excleude(ex) or include(in), -r ex, -c in "
  
),ncol=5,byrow=T)

opt = getopt(spec)

if ( is.null(opt$query_gene)) {
  cat(paste(getopt(spec, usage=T)))
  q();
}

query_gene<-opt$query_gene[1]
# cutoff_pos<-opt$cutoff_pos[1]
# cutoff_neg<-opt$cutoff_neg[1]
# cutoff_reg<-opt$cutoff_reg[1]
# query_gene <- "CCAT1"
# query_gene <- "chr10_97437191_97438703_rev"
# cutoff_pos <- 0.5
# cutoff_neg <- -0.5
# cutoff_reg <- "ex"
library(data.table)
library(WGCNA)
library(reshape2)
library(feather)
library(plotly)
norm_tbl<-fread(paste0("output/norm_readstable.txt"))
linc_coexp_pairs<-as.data.frame(read_feather(paste0("output/linc_coexp_pairs.feather")))
df_opt<-readRDS(paste0("output/opt.rds"))


################################################################################
rna_idx<-ifelse( length(which("lncRNA" %in% colnames(linc_coexp_pairs))) == 1,which("lncRNA" %in% colnames(linc_coexp_pairs)),
                 ifelse(length(which("circRNA" %in% colnames(linc_coexp_pairs))) == 1,which("circRNA" %in% colnames(linc_coexp_pairs)),"NA"))
if (rna_idx =="NA"){q();}
rnacolumn<-linc_coexp_pairs[,colnames(linc_coexp_pairs)[rna_idx]]

linc_coexp_pairs_filtered<-linc_coexp_pairs[ rnacolumn %in%  query_gene , ]

# linc_coexp_pairs_filtered<-linc_coexp_pairs_filtered[!is.na(linc_coexp_pairs_filtered$cor),]
# 
# if( cutoff_reg == "ex") {
#   linc_coexp_pairs_filtered<-linc_coexp_pairs_filtered[which(linc_coexp_pairs_filtered$cor > cutoff_pos |  linc_coexp_pairs_filtered$cor < cutoff_neg), ]
# } 
# if(cutoff_reg == "in" ) { 
#   linc_coexp_pairs_filtered<-linc_coexp_pairs_filtered[which(!(linc_coexp_pairs_filtered$cor > cutoff_pos &  linc_coexp_pairs_filtered$cor < cutoff_neg)), ]
# }
##################################################################################

# extract query gene expresion value
if(df_opt[grep("mode",row.names(df_opt)),]=="circ") {
  norm_tbl_circ<-fread(paste0("output/norm_readstable_circRNA.txt"))
  query_gene_vec<-data.table(t(norm_tbl_circ[gene == query_gene,2:ncol(norm_tbl),with=F]))[[1]]
} else {
  query_gene_vec<-data.table(t(norm_tbl[gene == query_gene,2:ncol(norm_tbl),with=F]))[[1]]
}



# bg_gene subset
number_subset_genes<-5000
bg_genes<-intersect(linc_coexp_pairs$co_exp_gene,norm_tbl$gene)
bg_genes<-bg_genes[!bg_genes %in% query_gene]
bg_genes<-sample(bg_genes,number_subset_genes)
#co_exp_genes<-bg_genes
rand_r<-function(reptimes,co_exp_genes,titles){
  
  rest_tbl<-t(norm_tbl[gene != query_gene  & gene %in% co_exp_genes ,2:ncol(norm_tbl),with=F])
  colnames(rest_tbl)<-norm_tbl[gene != query_gene & gene %in% co_exp_genes][["gene"]]
  ##### all gene END
  randt<-data.table(replicate(reptimes,sample(query_gene_vec)))
  corlist<-WGCNA::corAndPvalue(log2(randt+1),log2(rest_tbl+1),nThreads =6 ,method = "pearson", verbose=1)
  cortbl<-data.table(corlist$cor)
  cortbl<-cbind(data.table(V1=seq(1,reptimes)),cortbl)
  cortbl<-melt(cortbl, id=1)
  cortbl$variable<-as.character(cortbl$variable)
  cortbl<-cortbl[order(variable,decreasing = T)]
  cortbl$V1<-NULL
  
  #exped_pairs<-unique(linc_coexp_pairs[ rnacolumn %in%  query_gene & linc_coexp_pairs$co_exp_gene %in% norm_tbl$gene , c("co_exp_gene","cor") ])
  
  exped_pairs<-data.table(unique(linc_coexp_pairs[ rnacolumn %in%  query_gene & linc_coexp_pairs$co_exp_gene %in% co_exp_genes,c("co_exp_gene","cor") ]))
  exped_pairs<-exped_pairs[order(co_exp_gene,decreasing = T)]
  exped_pairs<-exped_pairs[rep(1:.N,each=reptimes)]
  names(exped_pairs)[2]<-"corr"
  xx<-cbind(cortbl,exped_pairs)
  xx<-xx[variable==co_exp_gene]
  xx$count<-ifelse(xx$corr < 0 & xx$corr < xx$value,0, ifelse(xx$corr > 0 & xx$corr > xx$value,0,1))
  
  pvalue<-sum(xx$count)/nrow(xx)
  pvalue<-paste0("p-value =",format(pvalue,digits =3))
  
  return(list(pvalue=pvalue,cortbl=cortbl))
  
  # ggplot(cortbl, aes(x=value)) + 
  #   geom_histogram(binwidth =0.01)+
  #   scale_x_continuous(breaks = seq(-1, 1, by = 0.1),limits = c(-1,1))+
  #   theme_bw(base_size = 18)+
  #   theme(
  #     axis.text.x = element_text(angle = 45, hjust = 1)
  #   )+
  #   ylab("Counts")+
  #   xlab("R value")+ 
  #   ggtitle(paste0(titles," ", pvalue))
  # annotate("text", x = 0, y = 0, label = pvalue,
  #          size = 20)#,hjust =2, vjust =25)
  # ggsave(paste0(titles,".png"),dpi = 300)
}

reptimes<-500

bg_genes_randam<-rand_r(500,bg_genes,"asas")

bg_genes_randam$cortbl$type<-"Rand" 



# generation of signifcantly correlated table
p_cutoff<-1
# rna_idx<-ifelse( length(which("lncRNA" %in% colnames(linc_coexp_pairs))) == 1,which("lncRNA" %in% colnames(linc_coexp_pairs)),
#                  ifelse(length(which("circRNA" %in% colnames(linc_coexp_pairs))) == 1,which("circRNA" %in% colnames(linc_coexp_pairs)),"NA"))
# if (rna_idx =="NA"){q();}
# rnacolumn<-linc_coexp_pairs[,colnames(linc_coexp_pairs)[rna_idx]]

# linc_coexp_pairs_filtered<-linc_coexp_pairs[ rnacolumn %in%  query_gene , ]

co_exp_tbl<-unique(linc_coexp_pairs_filtered[,c("co_exp_gene","cor","cor_p")])
co_exp_tbl<-co_exp_tbl[co_exp_tbl$co_exp_gene != query_gene,]
co_exp_tbl$type<-"Obs"
co_exp_tbl<-co_exp_tbl[!is.na(co_exp_tbl$cor),]
co_exp_tbl<-co_exp_tbl[co_exp_tbl$cor_p < p_cutoff,]
co_exp_tbl$cor_p<-NULL

names(bg_genes_randam$cortbl)<-names(co_exp_tbl)

co_exp_tbl<-rbind(bg_genes_randam$cortbl,co_exp_tbl)

# plot
p<-ggplot(co_exp_tbl[,], aes(cor,  fill = type, colour = type)) +
  geom_density(alpha=0.3)+
  # geom_histogram(aes(y=c(..count..[..group..==1]/sum(..count..[..group..==1]),
  #                        ..count..[..group..==2]/sum(..count..[..group..==2]))),
  #                position='dodge',alpha=0.1,bins = 100,size=0.2)+
  theme_linedraw(base_size =12)+
  xlim(-1, 1)+
  theme_linedraw(base_size =12)+
  xlab("correlation")#+
 # ylab("Percentage")#+
 # scale_y_continuous(labels = scales::percent)


# ggsave(paste0(paste(query_gene,cutoff_pos,cutoff_neg,cutoff_reg,"random",reptimes,"rep",number_subset_genes,"genes",sep = "_"),".png"),
#         width = 5,
#        height = 4)#,dpi = 300)

#ggplotly(p)

saveRDS(p,paste0("output/random_",query_gene,".rds"))
