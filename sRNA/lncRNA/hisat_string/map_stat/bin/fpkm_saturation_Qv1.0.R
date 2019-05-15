

#!/share/nas2/genome/biosoft/R/3.1.1/lib64/R/bin/Rscript

#******************************************************************part 1 usage ***************************************

library('getopt');
library(ggplot2)
library(reshape2)
require(RColorBrewer)
spec = matrix(c(
  'help'    , 'h', 0, "logical",
  'count'  , 'count', 1, "character",
  'fpkm'   , 'fpkm' , 1, "character",
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
      Rscript  fpkm_saturation_v3.R --count All_gene_count.list --fpkm All_gene_fpkm.list  --o ./*Saturtion.png
      
      Options: 
      --help     -h 	NULL 		get this help
      --count 		character 	    the count file  [forced]
      --fpkm		character	    the fpkm file 
      --o 		character 	    the outfile png    [forced]
      \n")
  q(status=1);
}
################################################################################################
if ( !is.null(opt$help) ) { print_usage(spec) }
if ( is.null(opt$count) )  { print_usage(spec) }
if ( is.null(opt$fpkm) ) {print_usage(spec) }
if ( is.null(opt$o) )	{ print_usage(spec) }
################################################################
##ID     L01     L02     L04     L05
#ENSG00000086015 2.124531        1.797813        1.969637        2.2617570000000002

# opt=list()
# opt$i="T01.geneExpression.xls"
# opt$o="T01.png"
#GeneID	Length	FPKM	Locus	Strand	Count	Normalized
#ENSG00000279928	718	0.0919267	1:182393-184158	+	1	0.857301

count=read.table(opt$count,sep="\t",header=T,comm="",check=F,row=1)
print ("read count ")
sam_name <- colnames(count)[ !(colnames(count)%in%c("geneLength")) ]
print ("get sample name")
count <- as.matrix(count[ ,sam_name ])
print ("matrix count")

fpkm=read.table(opt$fpkm,sep="\t",header=T,comm="",check=F,row=1)
print ("read fpkm")
fpkm <- as.matrix(fpkm[ ,sam_name ])
print ("matrix fpkm")

sam_num = length(sam_name)
print (sam_num)
for (i in 1:sam_num ){
	Total=sum(count[,sam_name[i]])	#计算样品i的Total count
	Len=count[,sam_name[i]]/fpkm[,sam_name[i]]/Total*1e9	#计算每个基因的count
	Len[is.na(Len)]=1	#空值赋值为1
	prob=count[,sam_name[i]]/Total	#计算每个基因count占所有count的比值
	count1=rep(0,nrow(count))	#获取矩阵行数后将0复制n行
	names(count1)=rownames(count)	#将基因名称赋给每行的0

	mysample=function(rate){
	  s=sample(rownames(count),size=ceiling(Total*rate),re=T,prob=prob)	
	  ct=table(s)
	  count1[names(ct)]=ct[names(ct)]
	  count1/sum(count1)/Len*1e9
	}
	d=cbind(sapply(seq(0.1,0.9,0.1),mysample),fpkm[,sam_name[i]])

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
fpkm_data=d[,10]
interval=apply(ch,2,function(x){
  #fm0=fpkm[x<0.20]   
  fm0=fpkm_data[x<0.15] 
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

result = paste(opt$o,"/",sam_name[i],".Saturation.png",sep = "")
ggsave(filename=result,width=width,height=height,unit="mm",dpi=300)
print (i)
}
