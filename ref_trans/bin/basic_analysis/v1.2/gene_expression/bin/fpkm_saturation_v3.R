

#!/share/nas2/genome/biosoft/R/3.1.1/lib64/R/bin/Rscript

#******************************************************************part 1 usage ***************************************

library('getopt');
library(ggplot2)
library(reshape2)
require(RColorBrewer)
spec = matrix(c(
  'help'    , 'h', 0, "logical",
  'i'  , 'i', 1, "character",
  'o'  , 'o', 1, "character",
  'pre', 'p', 1, "character"
), byrow=TRUE, ncol=4);

opt = getopt(spec);


# define usage function
print_usage <- function(spec=NULL){
  cat(getopt(spec, usage=TRUE));
  cat("Usage example: \n")
  cat("
      Usage example: 
      Rscript  fpkm_saturation_v3.R --i T01.geneExpression.xls  --o T01_saturation.png
      
      Options: 
      --help     -h 	NULL 		get this help
      --i 		character 	    the geneExpression file  [forced]
      --o 		character 	    the outfile png    [forced]
      \n")
  q(status=1);
}
################################################################################################
if ( !is.null(opt$help) ) { print_usage(spec) }
if ( is.null(opt$i) )  { print_usage(spec) }
if ( is.null(opt$o) )	{ print_usage(spec) }
################################################################

# opt=list()
# opt$i="T01.geneExpression.xls"
# opt$o="T01.png"

d=read.table(opt$i,sep="\t",header=T,comm="",check=F,row=1)

Total=sum(d[,6])

Len=d[,6]/d[,2]/Total*1e9
Len[is.na(Len)]=1
prob=d[,6]/Total
count=rep(0,nrow(d))
names(count)=rownames(d)

mysample=function(rate){
  s=sample(rownames(d),size=ceiling(Total*rate),re=T,prob=prob)
  ct=table(s)
  count[names(ct)]=ct[names(ct)]
  count/sum(count)/Len*1e9
}

d=cbind(sapply(seq(0.1,0.9,0.1),mysample),d[,2])

#------------------------------------------------------------------------
d=d[d[,10]>0.1,]

#--------------------------------------------------------------
#calculate the absolution difference percent out of final fpkm
#--------------------------------------------------------------
ch=apply(d,1,function(x){
  abs((x-x[10])/x[10])
})
ch=t(ch)

#------------------------------------------------------------------
#the gene number with  diff percnet <0.15 in ervery interval fpkm
#-----------------------------------------------------------------
fpkm=d[,10]
interval=apply(ch,2,function(x){
  #fm0=fpkm[x<0.20]   
  fm0=fpkm[x<0.15] 
  c(length(which(fm0<=1)),
    length(which(fm0<=5&fm0>1)),
    length(which(fm0<=10&fm0>5)),
    length(which(fm0<=50&fm0>10)),
    length(which(fm0>50))
  )
  
})

#-------------------------------------------------------
#the percent gene out of final value
#-------------------------------------------------------
interval=apply(interval,1,function(x){x/x[10]})
colnames(interval)=c("0-1","1-5","5-10","10-50",">50")
data=data.frame(precent=seq(10,100,10),interval,check.names =F)
data=melt(data,id=1)

#----------------------------------------------------------------
#ggplot 
#---------------------------------------------------------------
width=140
height=140
size_ratio=width/210
mytheme=theme_get()
mytheme$axis.title.x$size=14*size_ratio
mytheme$axis.title.y$size=14*size_ratio
mytheme$axis.text.x$size=12*size_ratio
mytheme$axis.text.y$size=12*size_ratio
mytheme$legend.title$size=11*size_ratio
mytheme$legend.text$size=9*size_ratio
theme_set(mytheme)

p=ggplot(data,aes(x=precent,y=value,colour=variable))+
  geom_line()+
  geom_point(size=2*width/210)+
  ylab("Percent of genes within 15% of final fpkm")+
  xlab("Mapped Reads(%)")+
  scale_colour_manual(values =brewer.pal(5,"Dark2"),name="fpkm interval")
# remove background, horizontal & vertical grid lines 
p <- p + theme_bw() + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())

ggsave(filename=opt$o,width=width,height=height,unit="mm",dpi=300)









