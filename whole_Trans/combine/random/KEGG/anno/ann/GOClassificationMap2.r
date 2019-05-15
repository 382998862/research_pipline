#!/share/nas2/genome/biosoft/R/2.15.1/lib64/R/bin/Rscript

# load library
library('getopt');

#opt<-data.frame(infile="E:/R_workplace/20150605GObarplot/test.GO_enrichment.stat.xls",rotation=75,txt.size=0.6)

#-----------------------------------------------------------------
# getting parameters
#-----------------------------------------------------------------
#get options, using the spec as defined by the enclosed list.
#we read the options from the default: commandArgs(TRUE).
spec = matrix(c(
  'help' , 'h', 0, "logical",
  'infile' , 'i', 1, "character",
  'rotation','r',2,"integer",
  'txt.size' , 't', 2, "double",
  "outpath",'o',2,"character",
  "key",'k',2,"character"
), byrow=TRUE, ncol=4);
opt = getopt(spec);


#遗传图与基因组的共线性分析
# define usage function
print_usage <- function(spec=NULL){
  cat(getopt(spec, usage=TRUE));
  cat("Usage example: \n")
  cat("
      Usage example: 
      1) Rscript GOClassificationMap.R --infile in.dat --outpath /share/nas1/map --txt.size 0.6
      2) Rscript GOClassificationMap.R --infile in.dat --rotation 0.5 --outpath /share/nas1
      
      
      
      
      
      Options: 
      --help          -h  NULL        get this help
      --infile        -i  character   the tab delimited input file saving functoin class,function,gene number,the line begin with # will be skipped[forced]
      --rotation      -r  integer     A numerical value specifying (in degrees) how the x axis annotation should be rotated[optional, default:75]
      --outpath       -o  character   the path for output file [optional,current working directory]
      --txt.size      -t  double      the size of the text in the map [optional, default: 0.5]
      --key       -k  character   the key words for output file [optional,default:'demo'] 
 \n")
  q(status=1);
}


# if help was asked for print a friendly message
# and exit with a non-zero error code
if ( !is.null(opt$help) ) { print_usage(spec) }


# check non-null args
if ( is.null(opt$infile) )  { print_usage(spec) }else {opt$infile<-gsub("\\\\",replacement = "/",x = opt$infile)}
infile<-as.vector(opt$infile)
#set some reasonable defaults for the options that are needed,
#but were not specified.
if ( is.null(opt$outpath) )  { opt$outpath =getwd() }
outpath<-gsub("\\\\",replacement = "/",x = opt$outpath)
if(!file.exists(outpath)){dir.create(outpath,recursive = T)}
if(is.null(opt$rotation)){opt$rotation=75}
if(is.null(opt$txt.size)){opt$txt.size=1}
if(is.null(opt$key)){opt$key='demo'}
rot<-as.numeric(opt$rotation)
txt.size<-as.numeric(opt$txt.size)
#
ppi=500
group<-c("molecular function","cellular component","biological process")
col7<-c("#A6CEE3", "#1F78B4" ,"#B2DF8A", "#33A02C", "#FB9A99", "#E31A1C")

#-----------------------------------------------------------------
# reading data
#-----------------------------------------------------------------
# reading data
#import in.dat, check data file formate
if(file.exists(infile)){
  indat<-read.table(as.vector(infile),sep="\t",comment.char = "#",header = F) 
  sample<-paste("dat",1:(ncol(indat)-2),sep="_")
  colnames(indat)<-c("class","funct",sample) 
  bardat<-as.matrix(indat[,sample])
  line1st<-as.matrix(read.table(infile,sep="\t",header = F,comment.char = "",nrows= 2,colClasses = "character"))  
  if(sum(is.element("#Total_gene",line1st))>0){
    genes.total<-as.numeric(line1st[which(line1st[,1]=="#Total_gene",arr.ind = T),-c(1:2)])
  }else{
    genes.total<-apply(bardat,2,max)
  }
  sam.lab=paste("Group",1:length(sample),sep="")
  if(sum(is.element("#GO_classify1",line1st[,1]))>0){
    sam.lab<-line1st[which(line1st[,1]=="#GO_classify1",arr.ind = T),-c(1:2)]
  }    
}else{
  stop("open infile error:infile does not exist!")
}

#sort bardat by the max gene.tota
class<-as.vector(indat[,1])
rk<-as.numeric(indat[,which.max(genes.total)+2])
indat<-indat[order(class,-rk),]
class<-as.vector(indat[,1])
func<-as.vector(indat[,2])
bardat<-as.matrix(indat[,-c(1:2)])
if(sum(!is.element(class,set = group))>0){
  stop(paste(c("the first collumn in input file contains Unknown classification:",setdiff(class,group)),collapse = "\t"))
}

#set labels colour
if(length(sample)==1)col=col1 else col=col7
names(col)<-paste(rep(group,each=length(sample)),sample,sep="_")
class<-paste(rep(class,each=length(sample)),sample,sep="_")
text.len<-max(nchar(as.matrix(indat)[,2]))*1.2
temp<-100*bardat/rep(genes.total,each=nrow(bardat))
width =(7+length(bardat))*0.1*ppi

height =1.75*text.len*0.09*sin(rot*pi/180)*txt.size*ppi
if(length(sample)==1){
  space=0.4
}else if(length(sample)==2){
  space=c(0,0.6)
}else if(length(sample)==3){
  space=c(0,0,0.6)
}


#par(op)
#对百分比进行线性转换
percent<-log10(temp)+1
percent[percent<0]<-0
png(paste(outpath,paste(opt$key,"GO.png",sep="."),sep = "/"),width =width,height =height,res = ppi,units = "px")
par(cex=1,cex.lab =txt.size,cex.axis = txt.size,font.lab = 2,font.axis = 2,mar = c(txt.size*text.len*0.09*sin(rot*pi/180)/par("csi"),3,2,ceiling(nchar(max(genes.total))/3)+2),mgp=c(1.6,0.5,0))
bar=barplot(t(percent),col=col[class],beside = T,border =NA,axes=F,names.arg = rep(NA,length(percent)),ylim=c(0,3), ylab="Percentage of genes",space =space)
mtext(side = 4,line =ceiling(nchar(max(genes.total))/3)+0.5,text = "Number of genes" ,cex=txt.size,font=2)
box(bty="u")
bar<-apply(bar,2,mean)
axis(side = 1,at = bar-(bar[2]-bar[1])*0.5,labels =F,tck=-0.02)
axis(side = 2,at=c(0,1,2,3),lab=c(0.1,1,10,100),cex=txt.size,tck=-0.02,font=2,las=2)
text(x=bar, y=-0.02,labels=func, srt=rot, adj=1, xpd=TRUE,tcl=-0.02,cex=txt.size,font = 2)
yl1<-c(0.1,1,10,100)*genes.total[1]/100
yl2<-c(0.1,1,10,100)*genes.total[2]/100
yl1[yl1<1]<-1
yl2[yl2<1]<-1
axis(side = 4,at=c(0,1,2,3),tck=-0.02,labels = F)
axis(side = 4,at=c(0,1,2,3)-0.1,lab=floor(yl1),cex=txt.size,tck=0,las=1,xpd=T,font = 2,col.axis="black")
axis(side = 4,at=c(0,1,2,3)+0.1,lab=floor(yl2),cex=txt.size,tck=0,col.axis="blue",las=1,xpd=T,font = 2)
le.x<-par("usr")[2]-1.2*(strwidth(c("a ","a a ","a a a "),units = "user")+max(strwidth(sam.lab,units = "user")))*txt.size
points(x =le.x,y=rep(3-0.1,3),pch=15,col=c("#A6CEE3","#B2DF8A","#FB9A99"),cex=2*txt.size,xpd=T)
points(x =le.x,y=rep(3+0.1,3),pch=15,col=c("#1F78B4", "#33A02C","#E31A1C"),cex=2*txt.size,xpd=T)
text(x=max(le.x),y=3-0.1,labels = sam.lab[1],cex=txt.size,pos=4,font = 2,col="black")
text(x=max(le.x),y=3+0.1,labels = sam.lab[2],cex=txt.size,pos=4,xpd=T,font = 2,col="blue")

xdis<-strwidth(indat[which.max(nchar(func)),2])*cos(rot*pi/180)*txt.size*1.2
yf<-tan(rot*pi/180)/strwidth("a")*strwidth("a","in")/strheight("a","in")*strheight("a")
ydis<-xdis*yf+0.05
for(i in 1:length(group)){
  lab<-indat[,1]==group[i]
  l<-func[lab]
  x1<-bar[lab][1]-strwidth(l[1])*cos(rot*pi/180)*txt.size*1.2
  y1<--yf*strwidth(l[1])*cos(rot*pi/180)*txt.size*1.2
  x2<-tail(bar[lab],1)-strwidth(tail(l,1))*cos(rot*pi/180)*txt.size*1.2
  y2<--yf*strwidth(tail(l,1))*cos(rot*pi/180)*txt.size*1.2
  segments(x0 = x1,y0 = y1-0.06,x1 =bar[lab][1]-xdis ,y1=-ydis-0.06,xpd=T)
  segments(x0=tail(bar[lab],1)-xdis,y0=-ydis-0.06,x1 =bar[lab][1]-xdis ,y1=-ydis-0.06,xpd=T)
  segments(x0=tail(bar[lab],1)-xdis,y0=-ydis-0.06,x1 =x2 ,y1=y2-0.06,xpd=T)    
  text(labels = group[i],x =mean(range(bar[lab]))-xdis,y=-ydis-0.1,xpd=T,cex=txt.size,adj =c(0.5,1),font = 2)
}
dev.off()
