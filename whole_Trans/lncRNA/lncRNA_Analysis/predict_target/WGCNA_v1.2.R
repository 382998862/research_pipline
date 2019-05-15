time1=Sys.time()
print(paste("Start time is",time1))

#!/share/nas2/genome/biosoft/R/2.15.1/lib64/R/bin/Rscript

# load library
library('getopt');

#-----------------------------------------------------------------
# getting parameters
#-----------------------------------------------------------------
#get options, using the spec as defined by the enclosed list.
#we read the options from the default: commandArgs(TRUE).
spec = matrix(c(
  'help' ,          'h', 0, "logical",
  'indir' ,         'i', 1, "character",
  'outdir' ,        'o', 1, "character",
  'is.log' ,        'J', 0, "logical",
  'meanFPKM' ,      'm', 1, "integer",
  'cormethod' ,     'c', 1, "character",
  'clustmethod' ,   'l', 1, "character",#'color' , 'r', 1, "character", #������ɫ
  'samplenum',      'n',1,"integer",
  'sampletype' ,    't', 1, "character",#��������
  'minModuleSize',  's' ,1, "integer",
  'fold' ,          'f' ,1, "integer",
  'ntop' ,          'p', 1, "integer"), 
  byrow=TRUE, ncol=4);
opt = getopt(spec);


# define usage function
print_usage <- function(spec=NULL){
  cat(getopt(spec, usage=TRUE));
  cat("Usage example: \n")
  cat("
      Usage example: 
      1) Rscript bmk_WGCNA_2.R --indir /share/nas13/ref/project/BMK141210-J66_Rice/analysis/Basic_Analysis/geneExpression/ \\
	  --outdir /share/nas21/chengxc/testing/linhj/WGCNA2/ --meanFPKM 0.5 --cormethod spearman \\
	  --samplenum 6 --sampletype control,test --ntop 150
      
      Options: 
      --help	-h 	NULL 		get this help
      --indir	character	the input file directory [forced]
	  --outdir	character	the output file directory [optional, default: ./WGCNA_Result/ ]
	  --meanFPKM	integer	mean FPKM number [optional, default: 0.5]
	  --cormethod	character	correlation method [optional, default: spearman]
	  --clustmethod	character	hclust method [optional, default: ward]
	  --samplenum	integer	sample number, default: 6]

	  --sampletype	character	type of sample [optional, default: c(control,test)]
	  --minModuleSize	integer	min Module gene number [optional, default: 50]
	  --fold	integer	cut module gene fold number [optional, default: 0.23]
	  --ntop	integer	top n genes to output [optional, default: 150]
	  --is.log	logical	log(data) [optional, default: FALSE]
	  \n")
  q(status=1);
}



# if help was asked for print a friendly message
# and exit with a non-zero error code
if ( !is.null(opt$help) ) { print_usage(spec) }


# check non-null args
if ( is.null(opt$indir) )	{ print_usage(spec) }

if ( is.null(opt$outdir) )	{ dir.create("./WGCNA_Result/");
opt$outdir = "./WGCNA_Result/" } else { dir.create(opt$outdir) }

#WGCNA_v1.2.R --indir $in --outdir $od/ --meanFPKM 0.5 -f 0.2 -n $count_sample
#set some reasonable defaults for the options that are needed,
#but were not specified.
if ( is.null(opt$meanFPKM ) )		{ opt$meanFPKM = 0.5 }
if ( is.null(opt$cormethod ) )		{ opt$cormethod = c('spearman') }
if ( is.null(opt$clustmethod ) )	{ opt$clustmethod = c('ward') }
#if ( is.null(opt$color ) )			{ opt$color = c(rep("blue",3),rep("red",3)) }
if ( is.null(opt$samplenum ) )		{ opt$samplenum = 6 }
if ( is.null(opt$sampletype ) )		{ opt$sampletype = c("control","test") }
if ( is.null(opt$minModuleSize ) )	{ opt$minModuleSize = 50 }
if ( is.null(opt$fold ) )			{ opt$fold = 0.23 }
if ( is.null(opt$ntop ) )			{ opt$ntop = 150 }
if ( is.null(opt$is.log) )			{ opt$is.log = FALSE }



#####--------------------------------------#####
#####--------------- get fpkm -------------#####
#####--------------------------------------#####
library(WGCNA)
options(stringsAsFactors = FALSE)

###get gene fpkm
a = list.files(opt$indir)
expression_file = a[grep("All_gene_fpkm.list$",a)]
merge.data = read.table(file = paste(opt$indir,"/",expression_file,sep=""),sep="\t",header=T,comm="",check=F)
n<-ncol(merge.data)-1
#get lncRNA fpkm
lncRNA_expression_file = a[grep("All_lncRNA_fpkm.list$",a)]
lncRNA_fpkm = read.table(file = paste(opt$indir,"/",lncRNA_expression_file,sep=""),sep="\t",header=T,comm="",check=F)

merge.data.final=rbind(merge.data,lncRNA_fpkm)
rownames(merge.data.final)=merge.data.final[,1]
write.table(merge.data.final, file=paste0(opt$outdir,"/","FPKM_merge.xls"),row.names=F, col.names=T, quote=FALSE,sep="\t")


#---------------------------------------------------------------------------#
#------------------------------ data mining  -------------------------------#
#---------------------------------------------------------------------------#
aa=merge.data.final[,-1]
aa[,n+1]=apply(aa[,c(1:length(aa))],1,mean)
fpkm=aa[aa[,n+1] > opt$meanFPKM ,1:n]
write.table(fpkm, file=paste0(opt$outdir,"/","FPKM_filter.xls"),row.names=F, col.names=T, quote=FALSE,sep="\t")

if( opt$is.log ){
	# check value < 0
	if( sum(fpkm<0)/sum(fpkm) > 0.5 ) {
		cat("Final Error: there exist half values in data < 0, option is.log ERROR\n")
		print_usage(spec)
	}
	# log
	data <- log2(fpkm+0.000001)
}



#####---------------------------------------------#####
#####-------------1. sample cluster --------------#####
#####---------------------------------------------#####


#sizeGrWindow(12, 9)
#png(filename = "Plots/1.1.Sample Clustering by (spearman).png",width = 600, height = 480)
c=cor(as.matrix(fpkm),method = opt$cormethod );
d <- as.dist(1-c); 
hr <- hclust(d, method = opt$clustmethod, members=NULL)
pdf(file=paste0(opt$outdir,"/","Fig_1_hclust.pdf"))
plot(as.dendrogram(hr),
     edgePar=list(col=1, lwd=2),
     horiz=F,
     main=paste("Sample Clustering by",opt$cormethod),
     cex=0.5)
dev.off()
png(filename = paste0(opt$outdir,"/","Fig_1_hclust.png"),width= 6000,height = 4800,res=800)      ###################################����png�����############
plot(as.dendrogram(hr),
     edgePar=list(col=1, lwd=2),
     horiz=F,
     main=paste("Sample Clustering by",opt$cormethod),
     cex=0.5)
dev.off()
#color=c(rep("blue",3),rep("red",3))

#library(ape)
#plot(as.phylo(hr), type = "fan", tip.color = labels2colors(c), label.offset = 0.01, color=labels2colors(c),
#     cex = log(mtcars$mpg, 10), col = "red")
#pdf(file="Fig.1.hclust_phylo.pdf",width=14,height = 14)
#plot(as.phylo(hr), type = "fan" , use.edge.length = TRUE,
#     node.pos = NULL, show.tip.label = TRUE, show.node.label = 1,
#     edge.color = "black", edge.width = 1.6, edge.lty = 1, font = 2,
#     cex = 1.5,#1/sqrt(hr$height)*0.9, 
#     adj = NULL, srt = 0, no.margin = F,
#     root.edge = T, label.offset = 0.01, underscore = T,
#     x.lim = NULL, y.lim = NULL, direction = "rightwards",
#     lab4ut = NULL, tip.color = color, plot = TRUE,
#     rotate.tree = 5, open.angle = 0, node.depth = 1)
#dev.off()


#-------------------------------------------------#
#------------- compare hclust result -------------#
#-------------------------------------------------#

#two random trees
#tree1<-as.phylo(hr)
#rtree(20) #random tree with 40 leaves

#c2=cor(as.matrix(fpkm.0.1),method="pearson");
#d2 <- as.dist(1-c2); 
#hr2 <- hclust(d2, method = "complete", members=NULL)

#tree2<-as.phylo(hr2)
#rtree(20)

#tree1$tip.label=paste(substr(tree1$tip.label,1,1),1:124,sep="-")
#tree2$tip.label=paste(substr(tree1$tip.label,1,1),1:124,sep="-")
#tree3$tip.label=paste(substr(tree1$tip.label,1,1),1:124,sep="-")

#creation of the association matrix
#association<-matrix(ncol=2, nrow=6)
#association[,1]<-association[,2]<-tree2$tip.label

#plot two tree
#pdf(file="Fig.2.compare.hclust_result.spearman_pearson.pdf",width=6,height = 8)
#cophyloplot(tree1, tree2, assoc=association, length.line=0, space=200, gap=35,
#            col = color,tip.color = color,edge.lty = 1, font = 2,edge.width = 15)
#dev.off()

#save.image(file=paste0(opt$outdir,"exp_data_hclust_analysis.RData"))


#-------------------------------------------------#
#----------------------- PCA result ---------------#
#-------------------------------------------------#

pca = prcomp(t(fpkm))


pdf(file=paste0(opt$outdir,"/","Fig_2_PCA_result_2D.pdf"))
#plot(pca$x[,c(1,2)],pch=16,col=rep(rainbow(length(opt$sampletype)),each=opt$samplenum/2),cex=1.5,main = "PCA map")
plot(pca$x[,c(1,2)],pch=16,col=rep(rainbow(opt$samplenum),each=1),cex=1.5,main = "PCA map")
text(pca$x[,c(1,2)],row.names(pca$x),col="black",pos=3)
#legend("topleft",legend=opt$sampletype,
#       pch=16,cex=1.2,col=rainbow(length(opt$sampletype)),bty="p")
legend("topleft",legend=row.names(pca$x),
       pch=16,cex=1.2,col=rainbow(opt$samplenum),bty="p")
dev.off()
png(filename = paste0(opt$outdir,"/","Fig_2_PCA_result_2D.png"),width= 6000,height = 4800,res=800)      ###################################����png�����############
#plot(pca$x[,c(1,2)],pch=16,col=rep(rainbow(length(opt$sampletype)),each=opt$samplenum/2),cex=1.5,main = "PCA map")
plot(pca$x[,c(1,2)],pch=16,col=rep(rainbow(opt$samplenum),each=1),cex=1.5,main = "PCA map")
text(pca$x[,c(1,2)],row.names(pca$x),col="black",pos=3)
#legend("topleft",legend=opt$sampletype,
#       pch=16,cex=1.2,col=rainbow(length(opt$sampletype)),bty="p")
legend("topleft",legend=row.names(pca$x),
       pch=16,cex=1.2,col=rainbow(opt$samplenum),bty="p")
dev.off()

################################# �ж��������� #######################################################
if(opt$samplenum >= 3){
	library(scatterplot3d)
	pdf(file=paste0(opt$outdir,"/","Fig_3_PCA_result_3D.pdf"))
	#scatterplot3d(pca$x[,1:3], highlight.3d=F, col.axis="black",color = rep(rainbow(length(opt$sampletype)),each=opt$samplenum/2),cex.symbols=1.5,cex.lab=1.5,cex.axis=1.4,
	#	      col.grid="lightblue", main="PCA map", pch=16)
	#legend("topright",legend = opt$sampletype ,
	#       pch=16,cex=1.2,col=rainbow(length(opt$sampletype)), ncol = 1) #"bottomright"
	scatterplot3d(pca$x[,1:3], highlight.3d=F, col.axis="black",color = rep(rainbow(opt$samplenum),each=1),cex.symbols=1.5,cex.lab=1.5,cex.axis=1.4,
              col.grid="lightblue", main="PCA map", pch=16)
	legend("topright",legend = row.names(pca$x) ,pch=16,cex=1.2,col=rainbow(opt$samplenum), ncol = 1) #"bottomright"
	dev.off()
	png(filename = paste0(opt$outdir,"/Fig_3_PCA_result_3D.png"),width= 6000,height = 4800,res=800)      ###################################����png�����############
	#scatterplot3d(pca$x[,1:3], highlight.3d=F, col.axis="black",color = rep(rainbow(length(opt$sampletype)),each=opt$samplenum/2),cex.symbols=1.5,cex.lab=1.5,cex.axis=1.4,
	#	      col.grid="lightblue", main="PCA map", pch=16)
	#legend("topright",legend = opt$sampletype ,
	#       pch=16,cex=1.2,col=rainbow(length(opt$sampletype)), ncol = 1) #"bottomright"
	scatterplot3d(pca$x[,1:3], highlight.3d=F, col.axis="black",color = rep(rainbow(opt$samplenum),each=1),cex.symbols=1.5,cex.lab=1.5,cex.axis=1.4,
              col.grid="lightblue", main="PCA map", pch=16)
	legend("topright",legend = row.names(pca$x) ,pch=16,cex=1.2,col=rainbow(opt$samplenum), ncol = 1) #"bottomright"
	dev.off()
	write.table(pca$x, file=paste0(opt$outdir,"/pca_result.xls"),row.names=T, col.names=NA, quote=FALSE,sep="\t")   ########################pca_result.xls##########
	#write.table()
}


#----------------------------------------------------------------------#
#-----------------ת�ñ���׾�����֤������Ϊ0�����----------------#
#----------------------------------------------------------------------#

datExpr = as.data.frame(t(fpkm[,]))#extract expression data as profile,NOTICE:t()
#sample(1:length(fpkm[,1]),5000)
gsg = goodSamplesGenes(datExpr, verbose = 3)
#If the last statement returns TRUE, all genes have passed the cuts.
gsg$allOK

if (!gsg$allOK)
{
  # Optionally, print the gene and sample names that were removed:
  if (sum(!gsg$goodGenes)>0)
    printFlush(paste("Removing genes:", paste(names(datExpr)[!gsg$goodGenes], collapse = ", ")));
  if (sum(!gsg$goodSamples)>0)
    printFlush(paste("Removing samples:", paste(rownames(datExpr)[!gsg$goodSamples], collapse = ", ")));
  # Remove the offending genes and samples from the data:
  datExpr = datExpr[gsg$goodSamples, gsg$goodGenes]
}
write.table(names(datExpr)[!gsg$goodGenes], file=paste0(opt$outdir,"/removeGene.xls"), row.names=FALSE, col.names=FALSE, quote=FALSE)
write.table(names(datExpr)[!gsg$goodSamples], file=paste0(opt$outdir,"/removeSample.xls"), row.names=FALSE, col.names=FALSE, quote=FALSE)


#--------------------------------------------------------------------------------------#
#-----------------------------------Phetmap correlation--------------------------------#
#--------------------------------------------------------------------------------------#

#Correlation in sample with sample
library(pheatmap)

#annotation = data.frame(Sample=factor(c(rep("Control",3),rep("Test",3))))
#annotation$Sample = factor(annotation$Sample, levels = c("Control","Test"))
#rownames(annotation) = colnames(fpkm.0.1)

#ann_colors = unique(color)

pdf(file = paste0(opt$outdir,"/Fig_4_sample_correlation_pheatmap.pdf"));
pheatmap(c,color=colorRampPalette(c("green","black","red"))(100),border_color = NA,
         main = "Sample Correlation")
dev.off()
png(filename = paste0(opt$outdir,"/Fig_4_sample_correlation_pheatmap.png"),width= 6000,height = 4800,res=800)      ###################################����png�����############
pheatmap(c,color=colorRampPalette(c("green","black","red"))(100),border_color = NA,
         main = "Sample Correlation")
dev.off()

#-----------------------------------------------------------------------------------#
#----------------------------Loading Trait data-------------------------------------#
#-----------------------------------------------------------------------------------#
if(opt$samplenum >= 4){
	trait_matrix=as.data.frame(diag(1,n))
	colnames(trait_matrix)=colnames(fpkm)
	row.names(trait_matrix)=colnames(fpkm)
	write.table(trait_matrix,file = paste0(opt$outdir,"/trait_data.xls"),col.names =NA,row.names = T,quote = F,sep="\t")

	allTraits=trait_matrix
	#traitData=read.table("trait_data_6_sample.txt",header=T,sep="\t",row.names = 1)# Load trait information
	#allTraits = traitData#[1:25,] #-c(59:61)];# Remove information we don't need
	#allTraits

	Samples = rownames(datExpr)#get sample names from expression data
	datTraits = allTraits;
	#rownames(datTraits) = allTraits[traitRows, 1];
	collectGarbage();# Performs garbage collection until free memory idicators show no change.
	#We now have the expression data in the variable datExpr, and the corresponding clinical traits in the variable datTraits.
	#Before we continue with network construction and module detection,
	#we visualize how traits relate to the sample dendrogram.
	# Re-cluster samples

	sampleTree2 = hclust(dist(datExpr,method="euclidean"), method = opt$clustmethod)
	# Convert traits to a color representation: white means low, red means high, grey means missing entry
	traitColors = numbers2colors(datTraits, signed = FALSE);
	# Plot the sample dendrogram and the colors underneath.
	#par(mar = c(1,12,2,1))
	#sizeGrWindow(12, 9)
	pdf(file = paste0(opt$outdir,"Fig_5_Sample dendrogram and trait heatmap.pdf"))      ########################����pdf���
	plotDendroAndColors(sampleTree2, traitColors,
			    groupLabels = row.names(datTraits),
			    marAll = c(1, 8, 3, 1),
			    cex.rowText = 0.8,
			    main = "Sample dendrogram and trait heatmap")
	dev.off()
	png(filename = paste0(opt$outdir,"/Fig_5_Sample_dendrogram_and_trait_heatmap.png"),width= 6000,height = 4800,res=800)  #�޸Ĵ�С������ʹ��ʾ��������
	plotDendroAndColors(sampleTree2, traitColors,
			    groupLabels = row.names(datTraits),
			    marAll = c(1, 8, 3, 1),
			    cex.rowText = 0.8,
			    main = "Sample dendrogram and trait heatmap")
	dev.off()

	#save.image(file = paste0(opt$outdir,"exp_data_hclust_analysis_trait.RData"))



#---------------------------------------------------------------------#
#---------------------- Choose a set of soft-thresholding powers------#
#---------------------------------------------------------------------#

	powers = c(1:30)#c(c(1:10), seq(from = 12, to=60, by=2))
	# Call the network topology analysis function
	sft = pickSoftThreshold(datExpr, powerVector = powers, verbose = 5)
	# Plot the results:
	sizeGrWindow(9, 5)# in Rstudio, this argument is not available
	png(filename = paste0(opt$outdir,"Fig_6_Scale_Free_Topology_and_mean.png"),width = 6000, height = 4800,res=800);
	par(mfrow = c(1,2));
	par(mar=c(3,5,3,3))
	cex1 = 0.9;
	# Scale-free topology fit index as a function of the soft-thresholding power
	plot(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
	     xlab="Soft Threshold (power)",ylab="Scale Free Topology Model Fit,signed R^2",type="n",
	     main = paste("Scale independence"),ylim=c(0,1));
	text(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
	     labels=powers,cex=cex1,col="red");
	# this line corresponds to using an R^2 cut-off of h
	abline(h=0.90,col="red")# threshold=0.9, this could be alternate to 0.8 or anything you like.
	# Mean connectivity as a function of the soft-thresholding power
	plot(sft$fitIndices[,1], sft$fitIndices[,5],
	     xlab="Soft Threshold (power)",ylab="Mean Connectivity", type="n",
	     main = paste("Mean connectivity"))
	text(sft$fitIndices[,1], sft$fitIndices[,5], labels=powers, cex=cex1,col="red")
	dev.off()
	pdf(file = paste0(opt$outdir,"Fig_6_Scale Free Topology and mean.pdf"));         #############################����һ��PDF���###################                        
	par(mfrow = c(1,2));
	par(mar=c(3,5,3,3))
	cex1 = 0.9;
	# Scale-free topology fit index as a function of the soft-thresholding power
	plot(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
	     xlab="Soft Threshold (power)",ylab="Scale Free Topology Model Fit,signed R^2",type="n",
	     main = paste("Scale independence"),ylim=c(0,1));
	text(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
	     labels=powers,cex=cex1,col="red");
	# this line corresponds to using an R^2 cut-off of h
	abline(h=0.90,col="red")# threshold=0.9, this could be alternate to 0.8 or anything you like.
	# Mean connectivity as a function of the soft-thresholding power
	plot(sft$fitIndices[,1], sft$fitIndices[,5],
	     xlab="Soft Threshold (power)",ylab="Mean Connectivity", type="n",
	     main = paste("Mean connectivity"))
	text(sft$fitIndices[,1], sft$fitIndices[,5], labels=powers, cex=cex1,col="red")
	dev.off()

#---------------------------------------------------------------------
#---------------------------- Co-expression similarity and adjacency
#---------------------------------------------------------------------
	soft <- -sign(sft$fitIndices[,3])*sft$fitIndices[,2]
	softPower<-0
	maxpower<-0
	power<-0
	for (i in powers)
	{
	  if (maxpower < soft[i])
	  {maxpower <- soft[i]
	   power <- i
	  }
	  
	  if (soft[i]>=0.9)
	  {
	    softPower<-i #ѡ�����˽ṹ������ֵ
	    maxpower<-soft[i]
	    break
	  }
	  
	}
	if (softPower==0) softPower<-power
	if (maxpower<0.9) softPower<- 15


	softPower = softPower;#Based on pickSFT result, 15 is default.
	adjacency = adjacency(datExpr, type="signed",power = softPower);

#---------------------------------------------------------------------
#---------------------------- Topological Overlap Matrix (TOM)
#---------------------------------------------------------------------
	# Turn adjacency into topological overlap
	TOM = TOMsimilarity(adjacency,TOMType = "signed",TOMDenom = "mean");
	dissTOM = 1-TOM

#---------------------------------------------------------------------
#---------------------------- Clustering using TOM
#---------------------------------------------------------------------
	# Call the hierarchical clustering function
	# flashClust is faster than rountine hclust, you can use hclust instead.
	geneTree = hclust(as.dist(dissTOM), method = opt$clustmethod);
	#sizeGrWindow(12,9)
	png(filename = paste0(opt$outdir,"Fig_7_1_Gene clustering on TOM-based dissimilarity.png"),width = 6000,height = 4800,res=800)
	plot(geneTree, xlab="", sub="", main = "Gene clustering on TOM-based dissimilarity",
	     labels = FALSE, hang = 0.04,cex=2)
	dev.off()
	pdf(file = paste0(opt$outdir,"Fig_7_1_Gene clustering on TOM-based dissimilarity.pdf"))         #################����һ��pdf���##########
	plot(geneTree, xlab="", sub="", main = "Gene clustering on TOM-based dissimilarity",
	     labels = FALSE, hang = 0.04,cex=2)
	dev.off()
#---------------------------------------------------------------------
#---------------------------- Cut tree to identificate modules
#---------------------------------------------------------------------
	# There are several methods for branch cutting;
	# our standard method is the Dynamic Tree Cut from the package dynamicTreeCut.
	# We like large modules, so we set the minimum module size relatively high:
	minModuleSize = opt$minModuleSize;
	# Module identification using dynamic tree cut:
	dynamicMods = cutreeDynamic(dendro = geneTree, distM = dissTOM,
				    deepSplit = 2, pamRespectsDendro = FALSE,
				    minClusterSize = minModuleSize);
	#table(dynamicMods)

#---------------------------------------------------------------------
#---------------------------- Convert Module To Color
#---------------------------------------------------------------------
	#The function returned 48 modules labeled 1???48 largest to smallest.
	# Label 0 is reserved for unassigned genes.
	# The above command lists the sizes of the modules.
	# We now plot the module assignment under the gene dendrogram
	# Convert numeric lables into colors
	dynamicColors = labels2colors(dynamicMods)
	#table(dynamicColors)
	# Plot the dendrogram and colors underneath
	#sizeGrWindow(8,6)
	png(filename = paste0(opt$outdir,"Fig_7_2_Gene dendrogram and module colors.png"),width = 6000,height = 4800,res=800)
	plotDendroAndColors(geneTree, dynamicColors, "Dynamic Tree Cut",
			    dendroLabels = FALSE, hang = 0.03,
			    addGuide = TRUE, guideHang = 0.05,
			    main = "Gene dendrogram and module colors")
	dev.off()
	pdf(file = paste0(opt$outdir,"Fig_7_2_Gene dendrogram and module colors.pdf"))        ###############����һ��PDF���################
	plotDendroAndColors(geneTree, dynamicColors, "Dynamic Tree Cut",
			    dendroLabels = FALSE, hang = 0.03,
			    addGuide = TRUE, guideHang = 0.05,
			    main = "Gene dendrogram and module colors")
	dev.off()
#---------------------------------------------------------------------
#---------------------------- merging of similar modules
#---------------------------------------------------------------------
#The Dynamic Tree Cut may identify modules whose expression profiles are very similar.
#It may be prudent to merge such modules since their genes are highly co-expressed.
#To quantify co-expression similarity of entire modules,
#we calculate their eigengenes and cluster them on their correlation
# Calculate eigengenes
	MEList = moduleEigengenes(datExpr, colors = dynamicColors)
	MEs = MEList$eigengenes
	# Calculate dissimilarity of module eigengenes
	MEDiss = 1-cor(MEs);
	# Cluster module eigengenes
	METree = hclust(as.dist(MEDiss), method = opt$clustmethod);

	tmp<-0
	dimModule<-dim(MEDiss)[1]
	for (i in 1:dimModule)
	{
	  for (j in 1:dimModule)
	  {
	    tmp<-c(tmp,MEDiss[i,j]) 
	  }
	}
	MEDissThres = opt$fold 

	# Plot the result
	par(mfrow = c(1,1));
	pdf(file = paste0(opt$outdir,"Fig_8_Clustering of module eigengenes.pdf"))
	plot(METree, main = "Clustering of module eigengenes",xlab = "", sub = "")
	abline(h=MEDissThres, col = "red")
	dev.off()
	png(filename = paste0(opt$outdir,"Fig_8_Clustering of module eigengenes.png"),width = 6000,height = 4800,res=800)  ############PNG#################
	plot(METree, main = "Clustering of module eigengenes",xlab = "", sub = "")
	abline(h=MEDissThres, col = "red")
	dev.off()
	# Call an automatic merging function
	merge = mergeCloseModules(datExpr, dynamicColors, cutHeight = MEDissThres, verbose = 3)
	# The merged module colors
	mergedColors = merge$colors;
	# Eigengenes of the new merged modules:
	mergedMEs = merge$newMEs;

	#To see what the merging did to our module colors,
	#we plot the gene dendrogram again,
	#with the original and merged module colors underneath
	png(filename = paste0(opt$outdir,"Fig_7_3_Dynamic Tree Cut.png"),width = 6000,height = 4800,res=800)
	plotDendroAndColors(geneTree, cbind(dynamicColors, mergedColors),
			    c("Dynamic Tree Cut", "Merged dynamic"),
			    dendroLabels = FALSE, hang = 0.03,
			    addGuide = TRUE, guideHang = 0.05)
	dev.off()
	pdf(file= paste0(opt$outdir,"Fig_7_3_Dynamic Tree Cut.pdf"))       ##################PDF###############
	plotDendroAndColors(geneTree, cbind(dynamicColors, mergedColors),
			    c("Dynamic Tree Cut", "Merged dynamic"),
			    dendroLabels = FALSE, hang = 0.03,
			    addGuide = TRUE, guideHang = 0.05)
	dev.off()
	# In the subsequent analysis, we will use the merged module colors in mergedColors.
	# Rename to moduleColors
	moduleColors = mergedColors
	# Construct numerical labels corresponding to the colors
	colorOrder = c("grey", standardColors(50));
	moduleLabels = match(moduleColors, colorOrder)-1;
	MEs = mergedMEs;
	# Save module colors and labels for use in subsequent parts
	#save(MEs, moduleLabels, moduleColors, geneTree, file = "blastomere-moduleIdentification.RData")


#---------------------------------------------------------------------
#----------------------------TOMplot for selected genes
#---------------------------------------------------------------------
	nGenes = ncol(datExpr);
	nSamples = nrow(datExpr);

	dissTOM2 = 1 - TOMsimilarityFromExpr(datExpr, power = softPower)  
	nSelect = opt$ntop * 10  
	set.seed(10)  
	select = sample(nGenes, size = nSelect,replace=TRUE)  
	selectTOM = dissTOM2[select, select]  
	selectTree = hclust(as.dist(selectTOM), method = opt$clustmethod)  
	selectColors = moduleColors[select]  
	plotDiss = selectTOM^7  
	diag(plotDiss) = NA  
	pdf(file = paste0(opt$outdir,"Fig_8_1_networkHeatmap.pdf"))
	TOMplot(plotDiss, selectTree, selectColors, main = "Network heatmap plot")  
	dev.off()
	png(filename = paste0(opt$outdir,"Fig_8_1_networkHeatmap.png"),width = 6000,height = 4800,res=800)      ###############png###############
	TOMplot(plotDiss, selectTree, selectColors, main = "Network heatmap plot")  
	dev.off()

#---------------------------------------------------------------------------#
#----------------------------Quantifying module trait associations ---------#
#---------------------------------------------------------------------------#
#correlate eigengenes with external traits
# Define numbers of genes and samples


	# Recalculate MEs with color labels
	MEs0 = moduleEigengenes(datExpr, moduleColors)$eigengenes
	MEs = orderMEs(MEs0)
	moduleTraitCor = cor(MEs, datTraits, use = "p");
	moduleTraitPvalue = corPvalueStudent(moduleTraitCor, nSamples);
	# Will display correlations and their p-values
	textMatrix = paste(signif(moduleTraitCor, 2), "\n(",
			   signif(moduleTraitPvalue, 1), ")", sep = "");
	dim(textMatrix) = dim(moduleTraitCor)

	# Display the correlation values within a heatmap plot
	#sizeGrWindow(3, 2)
	#png(filename = "Plots/Dynamic Tree Cut.png",width = 1280, height = 960)
	pdf(file = paste0(opt$outdir,"Fig_9_Module-Trait Correlation.pdf"))#, width = 36, height = 5.5)
	par(mar = c(6, 10, 2, 1.5))#set margin
	labeledHeatmap(Matrix = moduleTraitCor,
		       xLabels = colnames(datTraits),
		       yLabels = names(MEs),
		       ySymbols = names(MEs),
		       colorLabels = FALSE,
		       colors = blueWhiteRed(50),#colorRampPalette(c("green","black","red"))(100),
		       textMatrix = textMatrix,
		       setStdMargins = FALSE,
		       cex.text = 0.45,
		       zlim = c(-1,1),
		       yColorWidth=0.02,
		       xColorWidth = 0.05,
		       main = paste("Module-trait relationships"))
	dev.off()
	png(filename = paste0(opt$outdir,"Fig_9_Module-Trait Correlation.png"),width = 6000,height = 4800,res=800)#, width = 36, height = 5.5)
	par(mar = c(6, 10, 2, 1.5))#set margin
	labeledHeatmap(Matrix = moduleTraitCor,
		       xLabels = colnames(datTraits),
		       yLabels = names(MEs),
		       ySymbols = names(MEs),
		       colorLabels = FALSE,
		       colors = blueWhiteRed(50),#colorRampPalette(c("green","black","red"))(100),
		       textMatrix = textMatrix,
		       setStdMargins = FALSE,
		       cex.text = 0.45,
		       zlim = c(-1,1),
		       yColorWidth=0.02,
		       xColorWidth = 0.05,
		       main = paste("Module-trait relationships"))
	dev.off()


	#save.image(file=paste0(opt$outdir,"exp_data_hclust_analysis_trait_TOM_module.RData"))

#---------------------------------------------------------------------
#----------------------------Gene Significance and Module Membership
#---------------------------------------------------------------------
# Define variable time containing the time column of datTrait
#time = as.data.frame(datTraits$time);
#names(time) = "time"

	# names (colors) of the modules
	modNames = substring(names(MEs), 3)
	# MM value= cor(gene, ME)
	geneModuleMembership = as.data.frame(signedKME(datExpr, MEs));
	MMPvalue = as.data.frame(corPvalueStudent(as.matrix(geneModuleMembership), nSamples));
	names(geneModuleMembership) = paste("MM", modNames, sep="");
	names(MMPvalue) = paste("p.MM", modNames, sep="");
	# GS value = cor(gene, Trait)
	geneTraitSignificance = as.data.frame(cor(datExpr, datTraits, use = "p"));
	GSPvalue = as.data.frame(corPvalueStudent(as.matrix(geneTraitSignificance), nSamples));
	names(geneTraitSignificance) = paste("GS.", colnames(datTraits), sep="");
	names(GSPvalue) = paste("GS.", colnames(datTraits), sep="");
	# Now we have GS of each gene, MM of every member in a moduled
	# Then we can identify highly trait-relative gene and Modules.
	# The belowing Code illustrating that
	# genes highly significantly associated with a trait
	# are often also the most important (central) elements of modules associated
	# with the trait.
#---------------------------------------------------------------------
#----------------------------Visualizing the network of eigengenes(correlations between each module)
#---------------------------------------------------------------------

	#Plot network of eigengenes(modules)
	pdf(file = paste0(opt$outdir,"Fig_10_1_meta-module hclust and heatmap.pdf"));
	#Plot the relationships among the eigengenes and the trait
	plotEigengeneNetworks(MEs, "", marDendro = c(0,4,1,5.5),
			      heatmapColors=colorRampPalette(c("green","black","red"))(100),
			      marHeatmap = c(4,4,1,2), cex.lab = 0.7, xLabelsAngle = 90)
	dev.off()
	png(filename = paste0(opt$outdir,"Fig_10_1_meta-module hclust and heatmap.png"),width = 6000,height = 4800,res=800);  ####################png###############
	#Plot the relationships among the eigengenes and the trait
	plotEigengeneNetworks(MEs, "", marDendro = c(0,4,1,5.5),
			      heatmapColors=colorRampPalette(c("green","black","red"))(100),
			      marHeatmap = c(4,4,1,2), cex.lab = 0.7, xLabelsAngle = 90)
	dev.off()

#---------------------------------------------------------------------
#----------------------------Gene dendragram and trait
#---------------------------------------------------------------------
	geneTraitColor=as.data.frame(numbers2colors(geneTraitSignificance,signed=TRUE,colors = colorRampPalette(c("green","black","red"))(100)))
	geneTraitColor=geneTraitColor;
	names(geneTraitColor)= colnames(datTraits)
	pdf(file = paste0(opt$outdir,"Fig_11_1_gene_dendrogram_with_trait.pdf"));
	par(mar = c(3.5, 7, 2, 1))
	plotDendroAndColors(geneTree, cbind(mergedColors, geneTraitColor),
			    c("Module",colnames(datTraits)),dendroLabels = FALSE, hang = 0.03,
			    addGuide = TRUE, guideHang = 0.05)
	dev.off()
	png(filename = paste0(opt$outdir,"Fig_11_1_gene_dendrogram_with_trait.png"),width = 6000,height = 4800,res=800);  ##########################png###############
	par(mar = c(3.5, 7, 2, 1))
	plotDendroAndColors(geneTree, cbind(mergedColors, geneTraitColor),
			    c("Module",colnames(datTraits)),dendroLabels = FALSE, hang = 0.03,
			    addGuide = TRUE, guideHang = 0.05)
	dev.off()


#---------------------------------------------------
#----------------export MM and p.MM-----------------
#---------------------------------------------------

	MMlist=data.frame(colnames(datExpr),merge$color)
	names(MMlist)=c("ID","module")
	for (module in modNames){
	  oldname=names(MMlist)
	  MMlist=data.frame(MMlist, geneModuleMembership[,paste("MM",module,sep="")],MMPvalue[,paste("p.MM",module,sep="")]);
	  names(MMlist)=c(oldname,paste(paste("MM",module,sep="-")),paste("p.MM",module,sep="-"))}
	write.table(MMlist,file=paste0(opt$outdir,"MMlist.xls"),sep="\t",quote=F,row.names=F,col.names=T)


#---------------------------------------------------------------------
#----------------------------Visualizing the Cor of GS and MM in a given module
#---------------------------------------------------------------------

	# Interesting module colors
	for (module in modNames)
	{
	  column = match(module, modNames); # col number of interesting modules
	  moduleGenes = moduleColors==module;
	  #sizeGrWindow(7, 7);
	  pdf(file = paste(opt$outdir,"Fig_12_",module,"_Module membership vs gene significance.pdf",sep=""), width = 7, height = 7);
	  par(mfrow = c(1,1));
	  verboseScatterplot(abs(geneModuleMembership[moduleGenes, column]),
			     abs(geneTraitSignificance[moduleGenes, 1]),
			     xlab = paste("Module Membership in", module, "module"),
			     ylab = "Gene significance for body weight",
			     main = paste("Module membership vs. gene significance\n"),
			     cex.main = 1.2, cex.lab = 1.2, pch=19,cex.axis = 1.2, col = module)
	  dev.off()
	  png(filename = paste(opt$outdir,"Fig_12_",module,"_Module membership vs gene significance.png",sep=""),width = 6000,height = 4800,res=800); ##################png############
	  par(mfrow = c(1,1));
	  verboseScatterplot(abs(geneModuleMembership[moduleGenes, column]),
			     abs(geneTraitSignificance[moduleGenes, 1]),
			     xlab = paste("Module Membership in", module, "module"),
			     ylab = "Gene significance for body weight",
			     main = paste("Module membership vs. gene significance\n"),
			     cex.main = 1.2, cex.lab = 1.2, pch=19,cex.axis = 1.2, col = module)
	  dev.off()
	  
	}


#---------------------------------------------------------------------
#----------------------------Export All Genes all module
#---------------------------------------------------------------------
	intModules = modNames#c("turquoise","black","brown","blue","red","yellow") #,"pink","brown","black","purple","turquoise"
	for (module in intModules)  
	{
	  # Select module probes
	  modGenes = (moduleColors==module)
	  # Get their entrez ID codes
	  modLLIDs = names(datExpr)[modGenes];
	  assign(paste0("IDs_",module),modLLIDs)
	  # Write them into a file
	  fileName = paste(opt$outdir,"IDs-", module, ".txt", sep="");
	  write.table(as.data.frame(modLLIDs), file = fileName,
		      quote = FALSE,row.names = FALSE, col.names = FALSE)
	}

#---------------------------------------------------------------------
#----------------------------Export ntop Genes all module
#---------------------------------------------------------------------

	for (module in intModules) 
	{
	  # Select module probes
	  modGenes = (moduleColors==module)
	  # Get their entrez ID codes
	  modLLIDs_top = MMlist[modGenes,][,match(c("ID",paste("MM",module,sep="-")),colnames(MMlist))];
	  order_modLL_top = as.data.frame(modLLIDs_top[order(modLLIDs_top[,2]),])[1:opt$ntop,1]
	  assign(paste0("IDs_ntop_",module),order_modLL_top)
	  # Write them into a file
	  fileName = paste(opt$outdir,"IDs-top-",opt$ntop, module,".txt", sep="");
	  write.table(order_modLL_top, file = fileName,
		      quote = FALSE,row.names = FALSE, col.names = FALSE)
	}

	#save.image(file=paste0(opt$outdir,"exp_data_hclust_analysis_trait_interest_module_TOM.RData"))


#---------------------------------------------------------------------
#----------------------------Export All Genes In module to Visant
#---------------------------------------------------------------------
# Select module
	intModules = intModules #c("turquoise","black","brown","blue","red","yellow")
	for (module in intModules)
	{
	  # Select module probes were,probe = id. for microarray, probe is just probe.
	  probes = names(datExpr)
	  inModule = (moduleColors==module);
	  modProbes = probes[inModule];
	  # Select the corresponding Topological Overlap
	  modTOM = TOM[inModule, inModule];
	  dimnames(modTOM) = list(modProbes, modProbes)
	  # Export the network into an edge list file VisANT can read
	  vis = exportNetworkToVisANT(modTOM,
				      file = paste(opt$outdir,"VisANTInput-", module, "_all.txt", sep=""),
				      weighted = TRUE,
				      threshold = 0,
				      probeToGene = NULL)
	}


#---------------------------------------------------------------------
#----------------------------Export Top N Genes In module to Visant(use KME)
#---------------------------------------------------------------------

	intModules = intModules
	for (module in intModules)
	{
	  nTop = opt$ntop
	  topid=rownames(as.matrix(geneModuleMembership[order(geneModuleMembership[,paste("MM",module,sep="")],
							      decreasing=TRUE),]))#
	  topnumber=match(head(topid,n=nTop),colnames(datExpr))
	  #top300=(order(geneModuleMembership$kMEturquoise,decreasing=TRUE)<=6)
	  topProbes = head(topid,n=nTop);
	  # Select the corresponding Topological Overlap
	  testTOM = TOM[topnumber, topnumber];
	  dimnames(testTOM) = list(topProbes, topProbes)
	  vis = exportNetworkToVisANT(testTOM,
				      file = paste(opt$outdir,"VisANTInput-",module,"-top",nTop,".txt", sep=""),
				      weighted = TRUE,
				      threshold = 0,
				      probeToGene = NULL)
	}


#---------------------------------------------------------------------
#----------------------------Export Top N Genes In module to Cytoscape (use KME)
#---------------------------------------------------------------------

	intModules =intModules#c("turquoise","black","brown","blue","red","yellow")
	for (module in intModules)
	{
	  nTop=opt$ntop
	  topid=rownames(as.matrix(geneModuleMembership[order(geneModuleMembership[,paste("MM",module,sep="")],
							      decreasing=TRUE),]))#
	  topnumber=match(head(topid,n=nTop),colnames(datExpr))
	  #top300=(order(geneModuleMembership$kMEturquoise,decreasing=TRUE)<=6)
	  topProbes = head(topid,n=nTop);
	  inModule = (moduleColors==module);
	  # Select the corresponding Topological Overlap
	  testTOM = TOM[topnumber, topnumber];
	  dimnames(testTOM) = list(topProbes, topProbes)
	  cyt = exportNetworkToCytoscape(testTOM,
					 edgeFile = paste(opt$outdir,"Cytoscape_Input-edges-top",nTop,paste(module, collapse="-"), ".txt", sep=""),
					 nodeFile = paste(opt$outdir,"Cytoscape_Input-nodes-top",nTop,paste(module, collapse="-"), ".txt", sep=""),
					 weighted = TRUE,
					 threshold = 0,
					 nodeNames = topProbes,
					 nodeAttr = moduleColors[inModule][1:300]);
	}

#---------------------------------------------------------------------
#----------------------------Module Heatmap
#---------------------------------------------------------------------
	intModules = intModules#c("turquoise","black","brown","blue","red","yellow")
	for (module in intModules)
	{
	  #jpeg(file=paste("Fig.13.1.module_heatmap_",module,".jpg",sep=""),width = 3000, height = 3000,res = 300)
	  pdf(file=paste(opt$outdir,"Fig_13_1_module_heatmap_",module,".pdf",sep=""),width = 9, height = 9)
	  #tiff(file=paste(module," heatmap-0.2.tiff",sep=""),width = 3.5, height = 4,units = "in",res=96)
	  #library(devEMF);emf(file = paste(module," heatmap-0.2.emf",sep=""),width = 3.5, height = 4)
	  ME=MEs[,paste("ME",module,sep="")]
	  par(mfrow=c(3,1))
	  par(mar=c(0.3, 5, 3.5, 2.8))
	  plotMat(t(scale(datExpr[,moduleColors==module])),labels=NA,
		  nrgcols=30,rlabels=F,rcols=module,clabels=F,main=module,cex.main=2,axis(3,labels=F,tick=F))
	  axis(1,at=seq(1:134),labels=F,tick=T)
	  par(mar=c(5, 3, 1, 0.7))
	  barplot(ME, col=module, main="", cex.main=3,border = "grey50",
		  names.arg=row.names(datExpr),las=2,offset =0,cex.names=1)
	  box()
	  axis(1,at=seq(1:15),labels=F,tick=T,ylab="eigengene expression",
	       xlab="array sample")
	  par(mar=c(0, 0, 2.5, 0))
	 pie(c((length(names(datExpr[,moduleColors==module]))-length(grep("newGene",names(datExpr[,moduleColors==module])))-length(grep("MSTRG",names(datExpr[,moduleColors==module])))),
	    length(grep("newGene",names(datExpr[,moduleColors==module]))),length(grep("MSTRG",names(datExpr[,moduleColors==module])))),invert=T,
	  labels=paste(c("Known_gene","new_gene","new_lncRNA"),
		       c((length(names(datExpr[,moduleColors==module]))-length(grep("newGene",names(datExpr[,moduleColors==module])))-length(grep("MSTRG",names(datExpr[,moduleColors==module])))),
			 length(grep("newGene",names(datExpr[,moduleColors==module]))),length(grep("MSTRG",names(datExpr[,moduleColors==module])))),
		       sep="\n"),
	  col=c("gray","gray4","black"),
	  cex=1.3)
	  dev.off()
	  
	  
	  png(filename=paste(opt$outdir,"Fig_13_1_module_heatmap_",module,".png",sep=""),width = 6000,height = 4800,res=800)  ############## png #############
	  #tiff(file=paste(module," heatmap-0.2.tiff",sep=""),width = 3.5, height = 4,units = "in",res=96)
	  #library(devEMF);emf(file = paste(module," heatmap-0.2.emf",sep=""),width = 3.5, height = 4)
	  ME=MEs[,paste("ME",module,sep="")]
	  par(mfrow=c(3,1))
	  par(mar=c(0.3, 5, 3.5, 2.8))
	  plotMat(t(scale(datExpr[,moduleColors==module])),labels=NA,
		  nrgcols=30,rlabels=F,rcols=module,clabels=F,main=module,cex.main=2,axis(3,labels=F,tick=F))
	  axis(1,at=seq(1:134),labels=F,tick=T)
	  par(mar=c(5, 3, 1, 0.7))
	  barplot(ME, col=module, main="", cex.main=3,border = "grey50",
		  names.arg=row.names(datExpr),las=2,offset =0,cex.names=1)
	  box()
	  axis(1,at=seq(1:15),labels=F,tick=T,ylab="eigengene expression",
	       xlab="array sample")
	  par(mar=c(0, 0, 2.5, 0))
	#  pie(c(length(names(datExpr[,moduleColors==module]))-length(grep("newGene",names(datExpr[,moduleColors==module])))),
	#        length(grep("newGene",names(datExpr[,moduleColors==module])))),invert=T,
	#      labels=paste(c("Known_gene","new_gene"),
	#                   c(length(names(datExpr[,moduleColors==module]))-length(grep("newGene",names(datExpr[,moduleColors==module])))),
	#                     length(grep("newGene",names(datExpr[,moduleColors==module])))),sep="\n"),
	#      col=c("gray","black"),cex=1.3)
	  pie(c((length(names(datExpr[,moduleColors==module]))-length(grep("newGene",names(datExpr[,moduleColors==module])))-length(grep("MSTRG",names(datExpr[,moduleColors==module])))),
	    length(grep("newGene",names(datExpr[,moduleColors==module]))),length(grep("MSTRG",names(datExpr[,moduleColors==module])))),invert=T,
	    labels=paste(c("Known_gene","new_gene","new_lncRNA"),
		       c((length(names(datExpr[,moduleColors==module]))-length(grep("newGene",names(datExpr[,moduleColors==module])))-length(grep("MSTRG",names(datExpr[,moduleColors==module])))),
			 length(grep("newGene",names(datExpr[,moduleColors==module]))),length(grep("MSTRG",names(datExpr[,moduleColors==module])))),
		       sep="\n"),
	  col=c("gray","gray4","black"),
	  cex=1.3)
	  dev.off()
	}
	save.image(file=paste0(opt$outdir,"all_data.RData"))


	time2=Sys.time()
	print(paste("Start time is",time1))
	print(paste("End time is",time2))
	print(time2 - time1)
}else{
	print("Your sample number is not valid for WGCNA")
}
