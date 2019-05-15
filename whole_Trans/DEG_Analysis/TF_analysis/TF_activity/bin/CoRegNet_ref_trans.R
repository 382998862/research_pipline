library('getopt')

spec = matrix(c(
  'help', 'h', '0', "logical",
  'fpkm', 'e', '1', "character",
  'TF', 'f', '1', "character",
  'threads', 't', '2', "numeric",
  'out.prefix', 'o', '1', "character",
  'filter', 'l', '2', "numeric",
  'minTarg','m','2',"numeric",
  'fdr', 'd', '2', "numeric"
), byrow=TRUE, ncol=4)
opt = getopt(spec)

print_usage <- function(spec=NULL){
  cat(getopt(spec, usage=TRUE));
  cat("Usage example: \n");
  cat("
      Usage example: 
      1. Rscript CoRegNet_ref_trans.r --fpkm CIT_BLCA_EXP.txt --TF HumanTF.txt --out testres 
      
      Options: 
      --help		-h	NULL	get this help
      --fpkm		numeric;	the expression file [forced]
      --TF		character;       the list of transcription factors [forced]
      --threads   numeric;  Default set to 2. [optional]
      --out.prefix	character	the prefix of output file [forced]
      --filter	a threshold is defined for telling a gene belong to which class of three value categorical : over-expression (+1), under-expression (-1) and no change (0).
      --minTarg  the minimum number of targets for a regulator to be considered for actvity prediction. Default set to 10. 
      --fdr   The threshold to consider a pair of co-regulator significant (after pvalue correction). Default is 0.01.

\n")
	q(status=1);
}

if ( !is.null(opt$help) ) { print_usage(spec) }
if ( is.null(opt$fpkm) )	{ print_usage(spec) }
if ( is.null(opt$TF) )	{ print_usage(spec) }
if ( is.null(opt$threads) )	{ opt$threads = 4}
if ( is.null(opt$minTarg) )	{ opt$minTarg = 10}
if ( is.null(opt$fdr) )	{ opt$fdr = 0.01}
if ( is.null(opt$out.prefix) )	{ print_usage(spec) }


# load library
require(CoRegNet)
require(gplots)
require(RColorBrewer)
require(pheatmap)

#raw data
fpkm_data <- as.matrix(read.table(opt$fpkm, row.names = 1, header=TRUE,sep = "\t"))
fpkm_data <- log2(fpkm_data+1)
tfs <- as.matrix(read.table(opt$TF,header = F))


#Discretization of expression data and Reconstruct a large-scale regulatory network from gene expression data

options("mc.cores"=opt$threads)
if(!is.null(opt$filter)){
  if (is.na(opt$filter)) {
    disc=discretizeExpressionData(fpkm_data)
    grn = hLICORN(disc, TFlist=tfs)
  }else{
    disc=discretizeExpressionData(fpkm_data,threshold=opt$filter)
    grn = hLICORN(disc, TFlist=tfs)
  }
} else {
  grn = hLICORN(fpkm_data, TFlist=tfs)
}

#Calculate influence or transcriptional activity of TFs
CITinf=regulatorInfluence(grn,fpkm_data,minTarg = opt$minTarg)
corNet <- coregulators(grn,alpha = opt$fdr)
TFs <- as.matrix(rownames(CITinf))
colnames(TFs) <- "TFs"
influence <- cbind(TFs,CITinf)
write.table(influence,paste(opt$out.prefix,"_influence.txt",sep=""),sep="\t",quote=FALSE,row.names = F)
write.table(corNet,paste(opt$out.prefix,"_cornet.txt",sep=""),sep="\t",quote=FALSE,row.names = F)
write.table(grn@GRN,paste(opt$out.prefix,"_activity_grn.txt",sep=""),sep="\t",quote=FALSE,row.names = F)


######acquire cytoscape coregulate TFs
if(is.null(nrow(corNet) )| nrow(corNet) ==0){
  stop("No co-regulators in the provided network. If it was inferred with the hLICORN function, try a lower minCoregSupport.")
}
links = corNet[c(1:3, 5:7)]

#################################################################
### Drawing plot for coregulators
####
color <- colorRampPalette(c('#436eee','white','#EE0000'))(100)
if (length(colnames(CITinf))-1 > 15) {
        pdf(paste(opt$out.prefix,"_influences_heatmap.pdf",sep=""))
        pheatmap(CITinf, color = color, show_colnames=F, main="the heatmap of significant TFs' influence")
        dev.off()
}else{
        pdf(paste(opt$out.prefix,"_influences_heatmap.pdf",sep=""))
        pheatmap(CITinf, color = color, show_colnames=T, main="the heatmap of significant TFs' influence")
        dev.off()
}
######
require(igraph)
####
net <- graph.data.frame(links, directed = F)
V(net)$color <- "tomato"
V(net)$shape="circle"
V(net)$label.color <- "black"
E(net)$width <- log10(E(net)$nGRN)*0.1
E(net)$arrow.size <- .2   #设置箭头及边的大小
E(net)$edge.color <- "gray80"   #设置边的颜色
ll <- layout_nicely(net)
pdffile = paste(opt$out.prefix,"_network.pdf",sep="")
pdf(pdffile)
plot(net, layout=ll, vertex.label.cex=.7)
dev.off()

