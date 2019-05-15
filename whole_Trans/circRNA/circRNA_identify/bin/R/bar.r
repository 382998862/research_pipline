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
library(ggplot2)


# read count data
print("read count data ...")
count_data <- read.delim(args[1],header=TRUE,sep="\t")
head(count_data)
print("read count data is over")


#chr <- count_data[,"chr"]
#junction <- count_data[,5]

#mydat <- data.frame(a = chr, b = junction)
#chrsta <- tapply(mydat$b,mydat$a, function(t) sum(t))

data <- table(count_data[,"chr"])

lbls <- paste("chr",names(data))
colours <- rainbow(10)



#colours <- c("#8dd3c7", "#ffffb3", "#bebada", "#80b1d3", "#fb8072", "#fdb462", "#b3de69", "#d9d9d9", "#fccde5")
# output
file <- paste(args[2], "/chr", ".distrution.bar.png", sep="")
png(filename=file, height = 3000, width = 3400, res = 500, units = "px")
bar<-barplot(data,                            #条形图
        main="chr distrution",
        xlab="chr",ylab="cirRNA number",
        col=colours,
        names.arg=lbls,   #设定横坐标名称
        border=NA,                      #条形框不设置边界线
        font.main=4,                
        font.lab=2,
        beside=TRUE,xaxt="n")
text(x=apply(bar,2,mean),y=par('usr')[3],labels =rownames(data),xpd=T, srt=32,adj=c(1,0.5))
dev.off()






