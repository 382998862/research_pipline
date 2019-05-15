#!/share/nas2/genome/biosoft/R/current/bin/Rscript

# usage function
usage <- function() {
  print("-------------------------------------------------------------------------------")
  print("Usage: Rscript rpkm_density_plot_func.r read_count.txt outdir")
  print("1) ciriRNA.outfile")
  print("2) outdir: the dir for output")
  print("-------------------------------------------------------------------------------")
}


# get args
args <-commandArgs(TRUE)

# check args length
if( length(args) != 2 ) {
  print(args)
  usage()
  stop("the length of args != 2")
}


# load library
#library(plotrix)


# read count data
print("read count data ...")
# All_gene_counts.list
#test<-"e:/All_gene_counts.list"
count_data <- read.delim(args[1],header=TRUE,sep="\t")
head(count_data)
print("read count data is over")


tmp_data <-count_data

tmp1<- count_data[,-ncol(count_data)]
tmp2 <- as.data.frame(tmp1[,-1])
colnames(tmp2)<- colnames(tmp1)[-1]
Total_reads <- colSums(tmp2)

unique_count <-function(x) sum(x!=0)
Total_circRNA <-apply(tmp2,2,unique_count)
data <- data.frame(Total_reads,Total_circRNA)
write.table(data.frame(sample=rownames(data),data), file = args[2], row.names =F,quote=F,sep="\t")
### plot 
# T <- venn.diagram(list(A=blast,B=blastFast),filename = NULL,lwd=1,lty=2,col=c("red","green")
#                   ,cat.col=c("red","green"),rotation.degree=90)
# 
# grid.draw(T)
# 
# draw.pairwise.venn(area1=length_blast,area2=length_blastFast,cross.area=length_blast_blastFast
#                    ,category=c('blast','blastFast'),lwd=rep(1,1),lty=rep(2,2)
#                    ,col=c('red','green'),fill=c('red','green')
#                    ,cat.col=c('red','green')
#                    ,rotation.degree=10)
# 
# mytable <- table(slices)
# lbls <- paste(names(mytable), " : ", mytable, sep="")
# #colours <- rainbow(10)
# colours <- c("#8dd3c7", "#ffffb3", "#bebada", "#80b1d3", "#fb8072", "#fdb462", "#b3de69", "#d9d9d9", "#fccde5")
# # output
# file <- paste(args[2], "/circRNA", ".type.pie.png", sep="")
# png(filename=file, height = 3000, width = 3400, res = 500, units = "px")
# pie(mytable, labels = lbls, main="Pie Chart of circRNA type", col = colours)
# dev.off()







