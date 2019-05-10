#!/share/nas2/genome/biosoft/R/3.1.1/lib64/R/bin/Rscript

# usage function
usage <- function() {
	print("-------------------------------------------------------------------------------")
	print("Usage: Rscript fpkm_cor.r read_fpkm.txt outdir")
        print("1) read_fpkm.txt: the read fpkm file for RNA-SEQ")
	print("2) outdir: the dir for output")
	print("3) pearson or spearman default pearson")
	print("-------------------------------------------------------------------------------")
}


# get args
args <-commandArgs(TRUE)

# check args length
if( length(args) != 2 && length(args)!=3) {
	print(args)
	usage()
	stop("the length of args != 3")
}
if(length(args)==2){args[3]="pearson"}

library(Hmisc)
library(corrplot)
library(pheatmap)
##read fpkm
fpkm_data <- read.delim(args[1], row.names = 1, header=TRUE,check.names =F)
fpkm<-as.matrix(fpkm_data)

##filter low exp genes
which(rowSums(fpkm)==0)->low_exp
if(length(low_exp)>0){
	fpkm<-fpkm[-low_exp,]
}
####col numbers
sample_num <- dim(fpkm)[2]
print(sample_num)

cor<-rcorr(fpkm,type=args[3])

write.table(cor$P,paste(args[2],"/sample_pvalue.txt",sep=""),row.names=T,col.names=T,sep="\t",quote=F)
write.table(cor$r,paste(args[2],"/sample_coefficient.txt",sep=""),row.names=T,col.names=T,sep="\t",quote=F)

#########plot
pngfile <- paste(args[2], "/sample_corelation.png", sep="")
pdffile <- paste(args[2], "/sample_corelation.pdf", sep="")
if (sample_num > 25) {
	png(filename=pngfile)
	pheatmap(cor$r,display_numbers=F,color=colorRampPalette(c("orchid1","white","cyan"))(100), show_rownames=F, show_colnames=F)
	dev.off()

	pheatmap(cor$r,display_numbers=F,color=colorRampPalette(c("orchid1","white","cyan"))(100), filename=pdffile,show_rownames=F, show_colnames=F)
	dev.off()
	plotree <- pheatmap(cor$r,scale = "row")
	######
	newOrder=cor$r[plotree$tree_row$order, plotree$tree_col$order]
	write.table(newOrder,paste(args[2],"/cluster_tree.txt",sep=""),row.names=T,col.names=T,sep="\t",quote=F)
	file <- paste(args[2], "/cluster_tree.png", sep="")
	png(filename=file,height=3000, width=8000, res=500)
	plot(plotree$tree_row)
	dev.off()
}else{
	fontsize=10
	if(sample_num < 5)      fontsize=20
	if (sample_num <= 15){
		png(filename=pngfile)
	}else{
		png(filename=pngfile,height=5000, width=5000, res=500)
	}
        pheatmap(cor$r,display_numbers=T,color=colorRampPalette(c("orchid1","white","cyan"))(100),fontsize_number=fontsize)
        dev.off()

	pheatmap(cor$r,display_numbers=T,color=colorRampPalette(c("orchid1","white","cyan"))(100) ,filename=pdffile)
        dev.off()	
}
#corrplot(cor$r, type="upper", order="hclust", tl.col="black", tl.srt=45)
