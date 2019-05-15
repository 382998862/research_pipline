#!/share/nas2/genome/biosoft/R/3.1.1/lib64/R/bin/Rscript
library('getopt');
# usage function
spec = matrix(c(
  'help' , 'h', 0, "logical",
  'infile' , 'i', 1, "character",
  'outfile','o',1,"character",
  'label','l',1,"character"
), byrow=TRUE, ncol=4);
opt = getopt(spec);

#遗传图与基因组的共线性分析
# define usage function
print_usage <- function(spec=NULL){
  cat(getopt(spec, usage=TRUE));
  cat("Usage example: \n")
  cat("
      Usage example: 
      1) Rscript fpkm_density_plot_func.r --infile a.list,b.list,c.list --label miRNA,circRNA,lncRNA --outfile box.png



      
      Options: 
      --help          -h  NULL        get this help
      --infile        -i  character   the tab delimited input file saving numeric matrix of the values to be plotted.[forced],split by commas
      --outfile       -o  character   file path where to save the picture. Filetype is decided by the extension in the path. [optional,heatmap in current working directory]
      --label        -l      character RNA type,split by commas
      \n")
  q(status=1);
}
if (!is.null(opt$help)) { print_usage(spec) }

# check non-null args
if ( is.null(opt$infile) )  { print_usage(spec) }
if( is.null(opt$outfile))opt$outfile="fpkm _box"
# load library
require(edgeR)
require(ggplot2)


# read count data
print("read count data ...")
file <- unlist(strsplit(opt$infile, ","))
label <- unlist(strsplit(opt$label, ","))
count_data<-NULL
for(i in 1:length(file))
{
	# check geneLength key
	tmp <- read.delim(file[i],  header=TRUE,check.names = FALSE)
	if( !("geneLength" %in% colnames(tmp))) {
		stop("geneLength is not in read count file")
	}
	head(tmp)
	print("read circRNA count data is over")
	tmp[,1] <- label[i]
	if(i==1)
	{
		count_data<-tmp	
	}
	else
	{
		colnames(count_data)<-colnames(tmp)
		count_data <- rbind(count_data,tmp)
		colnames(count_data)[1]<- "ID"
	}
}
# calculate FPKM
sam_name <- colnames(count_data)[ !(colnames(count_data)%in%c("geneLength")) & !(colnames(count_data)%in%c("ID"))]
fpkm <- count_data[ ,sam_name ]
log10_fpkm <- fpkm
for ( i in 1:(dim(fpkm)[2]) ){
	fpkm[,i] <- 10^9 * fpkm[,i] / sum(fpkm[,i]) / count_data[,"geneLength"]
	log10_fpkm[,i] <- log10(fpkm[,i])
}
print("calculate FPKM is over")


# init 
all <- NULL
all_sam <- NULL
type <- NULL
# iter plot fpkm density
for( i in 1:length(sam_name) ){
	# fpkm
	c <- count_data[,sam_name[i]]
	keep <- c > 0
	r <- 10^9 * c[keep] / sum(c[keep]) / count_data[,"geneLength"][keep]
	log10fpkm <- data.frame(log10fpkm=log10(r))

	# update all
	all <- c(all, log10(r))
	### rep
	all_sam <- c(all_sam, rep(sam_name[i], length(r)))
	type <- c(type,count_data[,"ID"][keep])
}


# plot fpkm density for all
# create data.frame
log10fpkm <- data.frame(log10fpkm=all, sample=all_sam,type=type)
log10fpkm$sample <- factor(log10fpkm$sample)
log10fpkm$type <- factor(log10fpkm$type)
# plot fpkm box for all
# plot
m <- ggplot(log10fpkm, aes(x=sample, y=log10fpkm))
p <- theme_set(theme_bw())
p <- m + geom_boxplot(aes(fill=type)) + xlab("Sample") + ylab("log10(FPKM)")
p <- p + theme(axis.title.x = element_text(face="bold", size=14),
	axis.text.x  = element_text(face="bold", size=12),
	axis.title.y = element_text(face="bold", size=14),
	axis.text.y  = element_text(face="bold", size=12) )
p <- p + theme(legend.title = element_text(face="bold", size=14),
	legend.text = element_text(size=12) )
p <- p + theme(panel.grid.major =element_blank(), panel.grid.minor = element_blank(),panel.background = element_blank())
# output
file <- paste(opt$outfile, "/all", ".fpkm_box.png", sep="")
png(filename=file, height = 3000, width = 3400, res = 500, units = "px")
print(p)
dev.off()





