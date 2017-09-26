#!/usr/bin/env Rscript
#####!/share/apps/R/bin/Rscript

#######################
# scatter plot module #
#######################

args=commandArgs(TRUE)
library("getopt")


spec <- matrix(c(
  'query'     , 'q', 1, "character","X = query gene",
  'coexp'     , 'c', 1, "character","Y = co expression gene"
),ncol=5,byrow=T)

opt = getopt(spec);

if ( is.null(opt$query) | is.null(opt$coexp)) {
  cat(paste(getopt(spec, usage=T),"\nExample (linux): Rscript scatterplot.r -q ENSG00000232956 -c ENSG00000106603  
Example (linux): Rscript scatterplot.r -q chr19_21216990_21216261_fwd -c ENSG00000106603 \n"));
  q();
}

#reads table
circ_gene_merged<-readRDS("output/circ_gene_merged.rds") 
#linc_coexp_pairs<-readRDS("output/linc_coexp_pairs.rds")
df_opt<-readRDS("output/opt.rds")

# assign X and y for scatter plot
xx<-opt$query[1]
yy<-opt$coexp[1]


#lm_scatter = lm(log2(circ_gene_merged[,yy])~log2(circ_gene_merged[,xx]), data = circ_gene_merged)
#sum_lm_scatter <- summary(lm_scatter )
#r2_lm_scatter  = sum_lm_scatter $adj.r.squared
#pvalue_lm_scatter  = sum_lm_scatter $coefficients[2,4]

# calculation of correlation and pvalue
label_r2_lm_scatter  <- paste("Pearson's r = ", format(as.numeric(cor(log2(circ_gene_merged[,yy]),log2(circ_gene_merged[,xx]))) , digits = 3))
label_pvalue_lm_scatter  <- paste("p-value = ", format(cor.test(log2(circ_gene_merged[,yy]),log2(circ_gene_merged[,xx]))$p.value , digits = 3))


# plot
library("ggplot2")
ggplot((circ_gene_merged), aes(x= log2(circ_gene_merged[,xx]), y= log2(circ_gene_merged[,yy]) )) + 
  geom_point(aes(color=factor(attr)),size = 2,alpha=0.5)+
  theme_linedraw(base_size =12)+
  geom_smooth(method = "lm", se = FALSE, size =0.5 , color = "black")+
  scale_color_manual(values = c("blue", "red"),name="condition")+
  xlab(paste0(xx," (log2)"))+
  ylab(paste0(yy," (log2)"))+
  annotate("text", label = paste0(label_r2_lm_scatter,"\n",label_pvalue_lm_scatter) ,parse=F ,
           x = Inf, y = -Inf, color = "black",size = 3.5,hjust = 1.1, vjust = -1 ) 
ggsave(file=paste0("output/",xx,"_",yy,"_scatterplot.png"), width = 5, height = 4)

p<-ggplot((circ_gene_merged), aes(x= log2(circ_gene_merged[,xx]), y= log2(circ_gene_merged[,yy]) )) +
            geom_point(aes(color=factor(attr)),size = 2,alpha=0.5)+
            theme_linedraw(base_size =12)+
            geom_smooth(method = "lm", se = FALSE, size =0.5 , color = "black")+
            scale_color_manual(values = c("blue", "red"),name="condition")+
            xlab(paste0(xx," (log2)"))+
            ylab(paste0(yy," (log2)"))+
            annotate("text", label = paste0(label_r2_lm_scatter,"\n",label_pvalue_lm_scatter) ,parse=F ,
                     x = Inf, y = -Inf, color = "black",size = 3.5,hjust = 1.1, vjust = -1 )



