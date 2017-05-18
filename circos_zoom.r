#!/usr/bin/env Rscript
######!/share/apps/R/bin/Rscript
args=commandArgs(TRUE)


library(getopt)

# 1=required argument;2=optional argument
spec <- matrix(c(
  'query'     , 'q', 1, "character", "-q SNHG15, -q chr2_191537878_191523883_fwd",
  'top'       , 't', 2, "numeric", "top ranked of coexpressed genes, default:25",
  'folder'    , 'f', 1, "character", "-f output/ "
),ncol=5,byrow=T)

opt2 = getopt(spec);

if (is.null(opt2$query) || is.null(opt2$folder)) {
  cat(paste(getopt(spec, usage=T),"\n","Example:  ./circ_zoom.r -q chr2_191537878_191523883_fwd -f output
 Example: ./circ_zoom.r -q SNHG15 -f output \n"));
  q();
}



##### for test ######
#  opt2$query<-"chr2_191537878_191523883_fwd"
#  opt2$query<-"chrX_47755339_47705503_fwd"
#  opt2$query<-"CTD-2396E7.11"
# opt2$query<-"LINC00483"
# opt2$query<-"MIR4435-2HG"
# opt2$query<-"RP11-395B7.2"
# opt2$query<-"SNHG15"
# opt2$query<-"XXbac-B461K10.4"
# opt2$query<-"XXbac-B476C20.9"
# opt2$query<-"XXbac-B476C20.9"
# opt2$query<-"XXbac-B476C20.9"
#  opt2$folder<-"output"
# define default value
if ( is.null(opt2$top ) ) { opt2$top = 25 }
ncoexp<-opt2$top
# read cytoband info
library("circlize")


df_opt<-readRDS(paste0(opt2$folder[1],'/opt.rds'))
if (df_opt["annotation","V1"]=="gencodev25"){
  cytoBand = read.cytoband(species = "hg38")

}
if (df_opt["annotation","V1"]=="gencodev19"){
  cytoBand = read.cytoband(species = "hg19")
}

cytoband = cytoBand$df
chromosome = cytoBand$chromosome
chr.len = cytoBand$chr.len
chr.len<-chr.len[nchar(names(chr.len)) < 7 ]
chr.len<-chr.len[!names(chr.len) %in% "chrM"  ]
cytoband$V1<-as.character(cytoband$V1)
cytoband<-cytoband[nchar(cytoband$V1) < 7,]
cytoband<-cytoband[!cytoband$V1 %in% "chrM",]


#bed processing  
#linc_coexp_pairs<-readRDS("output/")   ### input table
linc_coexp_pairs<-readRDS(paste0(opt2$folder[1],'/linc_coexp_pairs.rds'))  
#  ### gene 

if (df_opt["mode","V1"]=="lnc"){
  bed<-linc_coexp_pairs[linc_coexp_pairs$lncRNA ==opt2$query[1],]
  bed2<-unique(bed[,c("chr","lncRNA_start","lncRNA_start","lncRNA")])  #
}
if (df_opt["mode","V1"]=="circ"){
  bed<-linc_coexp_pairs[linc_coexp_pairs$circRNA ==opt2$query[1],]
  bed2<-unique(bed[,c("circRNA_chr","circRNA_5exon","circRNA_5exon","circRNA")])  #
  
}


bed2$value<-1
bed2<-bed2[,c(1,2,3,5,4)]
colnames(bed2)<-c("chr","start","end","value","symbol")
bed2[,c(2:3)]<-sapply(bed2[,c(2:3)],as.numeric)

bed<-bed[,c("co_exp_gene_chr","co_exp_gene_start","cor","co_exp_gene")]
bed[,c(2:3)]<-sapply(bed[,c(2:3)],as.numeric)
bed$end<-bed$co_exp_gene_start
bed<-bed[,c(1,2,5,3,4)]
colnames(bed)<-c("chr","start","end","value","symbol")
bed<-bed[complete.cases(bed),]
bed<-bed[with(bed, order(chr, start)), ]

cytoband_zoom = cytoband[cytoband[[1]] %in% c(as.character(bed2$chr)), ]
cytoband<-cytoband[-grep(as.character(bed2$chr),cytoband$V1),]
cytoband_zoom[[1]] = paste0(cytoband_zoom[[1]])
cytoband = rbind(cytoband, cytoband_zoom)

bed_zoom = bed[bed[[1]] %in% c(as.character(bed2$chr)), ]
bed<-bed[-grep(as.character(bed2$chr),bed$chr),]
bed_zoom[[1]] = paste0(bed_zoom[[1]])
bed = rbind(bed, bed_zoom)


region2_inc<-bed[bed$value < 0,]# inverse correlation
region2_inc<-bed[order(bed$value) ,][1:ncoexp,] 
region2_inc$value<-0
region2_dec<-bed[bed$value > 0,]# positive correlation
region2_dec<-bed[order(bed$value,decreasing = T) ,][1:ncoexp,]  
region2_dec$value<-0
region1<-bed2[rep(1, each=nrow(region2_inc)),]

chr.len<-chr.len[-grep(as.character(bed2$chr),names(chr.len))]

png(paste0("output/",opt2$query,"_",opt2$top,".png"),units="px", width=800, height=800, res=150)
par(mar = c(1, 1, 1, 1))
circos.par(start.degree = 135)
circos.initializeWithIdeogram(cytoband, sort.chr = FALSE, sector.width = c(chr.len/sum(chr.len),0.4),plotType = c("axis", "labels"))

textcex<-1.2
if (nchar(opt2$query[1])> 10){textcex<-0.8}
if (nchar(opt2$query[1])> 15){textcex<-0.6}
circos.genomicTrackPlotRegion(bed2, ylim = c(0, 1), panel.fun = function(region, value, ...) {
  circos.genomicText(region, value, y = 0, adj = c(0.5, 0), labels = opt2$query[1] ,cex=textcex,facing = "inside")#, #niceFacing = TRUE,
                     #posTransform = posTransform.default)
}, bg.border = NA,track.height=0.05)

#circos.genomicTrackPlotRegion(region2_dec,ylim = c(0, 1), panel.fun = function(region, value, ...) {
#  circos.genomicText(region, value, y = 0, adj = c(0.5, 0),labels.column = 2,niceFacing = TRUE,
#                     posTransform = posTransform.default,cex=0.5,facing = "reverse.clockwise", ...)
#},bg.border = NA,track.height=0.1)

circos.genomicPosTransformLines(bed2,  direction = "outside",track.height = 0.05)#,posTransform = posTransform.default)


circos.genomicTrackPlotRegion(bed, panel.fun = function(region, value, ...) {
  col =  ifelse(value[[1]] > 0.75, "darkred", 
                ifelse(value[[1]] > 0.5,"red",
                       ifelse(value[[1]] <  -0.75,"darkblue",
                              ifelse(value[[1]] <  -0.55,"blue","lightgrey"))))
  cex = ifelse(value[[1]] > 0.75, 0.8, 
               ifelse(value[[1]] > 0.5, 0.2,
                      ifelse(value[[1]] <  -0.75,0.8,
                             ifelse(value[[1]] <  -0.5,0.2,0.1))))
  circos.genomicPoints(region, value, cex = cex, pch = 16, col = col, ...)
},track.height=0.35)


circos.genomicLink(region1, region2_inc, border = NA,lwd=0.25,col="blue")
circos.genomicLink(region1, region2_dec, border = NA,lwd=0.25,col="red")

circos.clear()       

dev.off()

