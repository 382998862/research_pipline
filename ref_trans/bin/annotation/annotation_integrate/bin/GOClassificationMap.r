#!/share/nas2/genome/biosoft/R/2.15.1/lib64/R/bin/Rscript

# load library
library('getopt');

#opt<-data.frame(infile="E:/R_workplace/20150605GObarplot/in.dat",rotation=75,txt.size=0.5,method="log",txt.size=0.6)

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
  "fileName",'f',2,"character",
  "method",'m',2,"character"
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
      2) Rscript GOClassificationMap.R --infile in.dat --rotation 0.5 --outpath /share/nas1 --method log
      
      
      
      
      
      Options: 
      --help          -h  NULL        get this help
      --infile        -i  character   the tab delimited input file saving functoin class,function,gene number,the line begin with # will be skipped[forced]
      --fileName      -f  character   the name of the saved file[optional,default:GO_classificatoin_map]
      --rotation      -r  integer     A numerical value specifying (in degrees) how the x axis annotation should be rotated[optional, default:75]
      --method        -m  character   the method used to trans data [optional,log]
      --outpath       -o  character   the path for output file [optional,current working directory]
      --txt.size      -t  double      the size of the text in the map [optional, default: 0.5]
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
if(is.null(opt$fileName)){opt$fileName ="GO_classificatoin_map"}
if(is.null(opt$rotation)){opt$rotation=75}
if(is.null(opt$txt.size)){opt$txt.size=0.5}
if(is.null(opt$method)){opt$method="log"}
method<-as.vector(opt$method)
rot<-as.numeric(opt$rotation)
txt.size<-as.numeric(opt$txt.size)

#
ppi=500
group<-c("molecular function","cellular component","biological process")
col1=c("#E41A1C","#377EB8","#4DAF4A")
col2<-c("#377EB8","#F781BF","#FF7F00")
col3<-c("#E41A1C","#FF7F00","#377EB8")
col4<-c("red","yellow","blue")
col5<-c("#1B9E77","#D95F02","#7570B3")
col6<-c("#7FC97F","#BEAED4","#FDC086")

#-----------------------------------------------------------------
# reading data
#-----------------------------------------------------------------
# reading data
#import in.dat, check data file formate
if(file.exists(infile)){
  indat<-as.matrix(read.table(as.vector(infile),sep="\t",comment.char = "#",header = F,colClasses = "character"))  
  colnames(indat)<-c("class","funct","dat")
  class<-indat[,"class"]
  bardat<-as.numeric(indat[,"dat"])
  
  line1st<-as.matrix(read.table(infile,sep="\t",header = F,comment.char = "",nrows= 2,colClasses = "character"))  
  if(sum(is.element("#Total_gene",line1st))>0){
    genes.total<-as.numeric(line1st[which(line1st[,1]=="#Total_gene",arr.ind = T),-c(1:2)])
  }else{
    genes.total<-apply(bardat,2,max)
  }
}else{
  stop("open infile error:infile does not exist!")
}

if(sum(!is.element(class,set = group))>0){
  stop(paste(c("the first collumn in input file contains Unknown classification:",setdiff(class,group)),collapse = "\t"))
}

#set labels colour
col<-col1
names(col)<-group
text.len<-max(nchar(indat[,2]))


temp<-100*bardat/genes.total
percent<-temp
width =(7+length(bardat))*0.1*ppi
height =((10+2)*0.2+text.len*0.09*sin(rot*pi/180)*txt.size)*ppi

if(is.na(opt[["method"]])){
  #par(op)
  #分区段进行线性变换
  percent[temp<=1&temp>=0.1]=temp[temp<=1&temp>=0.1]*33.3/0.9+0.1-33.3/9
  percent[temp<=10&temp>=1]=temp[temp<=10&temp>=1]*33.3/9+33.4-33.3/9
  percent[temp<=100&temp>=10]=temp[temp<=100&temp>=10]*33.3/90+66.7-33.9/9  
  
  png(paste(c(outpath,"/",opt$fileName,".png"),collapse =  ""),width =width,height =height,res = ppi,units = "px")
  par(cex=1,cex.lab =txt.size,cex.axis = txt.size,mar = c(txt.size*text.len*0.09*sin(rot*pi/180)/par("csi")+0.5,3,1.5,3),mgp=c(1.6,0.5,0))
  bar=barplot(percent,col=col[indat[,"class"]],border =NA,axes=F,ylim=c(0.1,100),names.arg =NA, ylab="Percentage of genes",space =0.4,)
  mtext(side = 4,line =1.6,text = "Number of genes" ,cex=txt.size,font=2)
  box(bty="n")
  axis(side = 1,at = bar-(bar[2]-bar[1])*0.5,labels =F,tck=-0.02)
  axis(side = 2,at=c(0.1,33.4,66.7,100),lab=c(0.1,1,10,100),cex=txt.size,tck=-0.02,xpd=T,las=1)
  text(x=bar, y=-1.5,labels=indat[,2], srt=rot, adj=1, xpd=TRUE,tcl=-0.02,cex=txt.size)
  yl<-c(0.1,1,10,100)*genes.total/100
  yl[yl<1]<-1
  yl<-round(yl)
  axis(side = 4,at=c(0.1,33.4,66.7,100),lab=yl,cex=txt.size,tck=-0.02,xpd=T,las=1)
  
  xdis<-strwidth(indat[which.max(nchar(indat[,2])),2])*cos(rot*pi/180)*txt.size
  yf<-tan(rot*pi/180)/strwidth("a")*strwidth("a","in")/strheight("a","in")*strheight("a")
  ydis<-xdis*yf+2.5
  for(i in 1:length(group)){
    lab<-class==group[i]
    l<-indat[lab,2]
    x1<-bar[lab][1]-strwidth(l[1])*cos(rot*pi/180)*txt.size
    y1<--yf*strwidth(l[1])*cos(rot*pi/180)*txt.size
    x2<-tail(bar[lab],1)-strwidth(tail(l,1))*cos(rot*pi/180)*txt.size
    y2<--yf*strwidth(tail(l,1))*cos(rot*pi/180)*txt.size
    segments(x0 = x1,y0 = y1-3,x1 =bar[lab][1]-xdis ,y1=-ydis-3,xpd=T)
    segments(x0=tail(bar[lab],1)-xdis,y0=-ydis-3,x1 =bar[lab][1]-xdis ,y1=-ydis-3,xpd=T)
    segments(x0=tail(bar[lab],1)-xdis,y0=-ydis-3,x1 =x2 ,y1=y2-3,xpd=T)    
    text(labels = group[i],x =mean(range(bar[lab]))-xdis,y=-ydis-5,xpd=T,cex=txt.size,adj =c(0.5,1) )
  }
  dev.off()
  
  pdf(paste(c(outpath,"/",opt$fileName,".pdf"),collapse =  ""),width =width/ppi,height =height/ppi)
  par(cex=1,cex.lab =txt.size,cex.axis = txt.size,mar = c(txt.size*text.len*0.09*sin(rot*pi/180)/par("csi")+0.5,3,1.5,3),mgp=c(1.6,0.5,0))
  bar=barplot(percent,col=col[indat[,"class"]],border =NA,axes=F,ylim=c(0.1,100),names.arg =NA, ylab="Percentage of genes",space =0.4,)
  mtext(side = 4,line =1.6,text = "Number of genes" ,cex=txt.size,font=2)
  box(bty="n")
  axis(side = 1,at = bar-(bar[2]-bar[1])*0.5,labels =F,tck=-0.02)
  axis(side = 2,at=c(0.1,33.4,66.7,100),lab=c(0.1,1,10,100),cex=txt.size,tck=-0.02,xpd=T,las=1)
  text(x=bar, y=-1.5,labels=indat[,2], srt=rot, adj=1, xpd=TRUE,tcl=-0.02,cex=txt.size)
  yl<-c(0.1,1,10,100)*genes.total/100
  yl[yl<1]<-1
  yl<-round(yl)
  axis(side = 4,at=c(0.1,33.4,66.7,100),lab=yl,cex=txt.size,tck=-0.02,xpd=T,las=1)
  
  xdis<-strwidth(indat[which.max(nchar(indat[,2])),2])*cos(rot*pi/180)*txt.size
  yf<-tan(rot*pi/180)/strwidth("a")*strwidth("a","in")/strheight("a","in")*strheight("a")
  ydis<-xdis*yf+2.5
  for(i in 1:length(group)){
    lab<-class==group[i]
    l<-indat[lab,2]
    x1<-bar[lab][1]-strwidth(l[1])*cos(rot*pi/180)*txt.size
    y1<--yf*strwidth(l[1])*cos(rot*pi/180)*txt.size
    x2<-tail(bar[lab],1)-strwidth(tail(l,1))*cos(rot*pi/180)*txt.size
    y2<--yf*strwidth(tail(l,1))*cos(rot*pi/180)*txt.size
    segments(x0 = x1,y0 = y1-3,x1 =bar[lab][1]-xdis ,y1=-ydis-3,xpd=T)
    segments(x0=tail(bar[lab],1)-xdis,y0=-ydis-3,x1 =bar[lab][1]-xdis ,y1=-ydis-3,xpd=T)
    segments(x0=tail(bar[lab],1)-xdis,y0=-ydis-3,x1 =x2 ,y1=y2-3,xpd=T)    
    text(labels = group[i],x =mean(range(bar[lab]))-xdis,y=-ydis-5,xpd=T,cex=txt.size,adj =c(0.5,1) )
  }
  dev.off()
  
  
}else if(opt[["method"]]=="log"){
  #par(op)
  #对百分比进行线性转换
  percent<-log10(temp)+1
  percent[percent<0]<-0
  png(paste(c(outpath,"/",opt$fileName,".png"),collapse =  ""),width =width,height =height,res = ppi,units = "px")
  par(cex=1,cex.lab =txt.size,cex.axis = txt.size,mar = c(txt.size*text.len*0.09*sin(rot*pi/180)/par("csi")+0.5,3,1.5,3),mgp=c(1.6,0.5,0))
  bar=barplot(percent,col=col[indat[,"class"]],border =NA,axes=F,ylim=c(0,3),names.arg =NA, ylab="Percentage of genes",space =0.4,)
  mtext(side = 4,line =1.6,text = "Number of genes" ,cex=txt.size,font=2)
  box(bty="n")
  axis(side = 1,at = bar-(bar[2]-bar[1])*0.5,labels =F,tck=-0.02)
  axis(side = 2,at=c(0,1,2,3),lab=c(0.1,1,10,100),cex=txt.size,tck=-0.02,xpd=T,las=1)
  text(x=bar, y=-0.02,labels=indat[,2], srt=rot, adj=1, xpd=TRUE,tcl=-0.02,cex=txt.size)
  yl<-c(0.1,1,10,100)*genes.total/100
  yl[yl<1]<-1
  yl<-round(yl)
  axis(side = 4,at=c(0,1,2,3),lab=yl,cex=txt.size,tck=-0.02,xpd=T,las=1)
  
  xdis<-strwidth(indat[which.max(nchar(indat[,2])),2])*cos(rot*pi/180)*txt.size
  yf<-tan(rot*pi/180)/strwidth("a")*strwidth("a","in")/strheight("a","in")*strheight("a")
  ydis<-xdis*yf+0.05
  for(i in 1:length(group)){
    lab<-class==group[i]
    l<-indat[lab,2]
    x1<-bar[lab][1]-strwidth(l[1])*cos(rot*pi/180)*txt.size
    y1<--yf*strwidth(l[1])*cos(rot*pi/180)*txt.size
    x2<-tail(bar[lab],1)-strwidth(tail(l,1))*cos(rot*pi/180)*txt.size
    y2<--yf*strwidth(tail(l,1))*cos(rot*pi/180)*txt.size
    segments(x0 = x1,y0 = y1-0.06,x1 =bar[lab][1]-xdis ,y1=-ydis-0.06,xpd=T)
    segments(x0=tail(bar[lab],1)-xdis,y0=-ydis-0.06,x1 =bar[lab][1]-xdis ,y1=-ydis-0.06,xpd=T)
    segments(x0=tail(bar[lab],1)-xdis,y0=-ydis-0.06,x1 =x2 ,y1=y2-0.06,xpd=T)    
    text(labels = group[i],x =mean(range(bar[lab]))-xdis,y=-ydis-0.1,xpd=T,cex=txt.size,adj =c(0.5,1) )
  }
  dev.off()
  
  
  pdf(paste(c(outpath,"/",opt$fileName,".pdf"),collapse =  ""),width =width/ppi,height =height/ppi)
  par(cex=1,cex.lab =txt.size,cex.axis = txt.size,mar = c(txt.size*text.len*0.09*sin(rot*pi/180)/par("csi")+0.5,3,1.5,3),mgp=c(1.6,0.5,0))
  bar=barplot(percent,col=col[indat[,"class"]],border =NA,axes=F,ylim=c(0,3),names.arg =NA, ylab="Percentage of genes",space =0.4,)
  mtext(side = 4,line =1.6,text = "Number of genes" ,cex=txt.size,font=2)
  box(bty="n")
  axis(side = 1,at = bar-(bar[2]-bar[1])*0.5,labels =F,tck=-0.02)
  axis(side = 2,at=c(0,1,2,3),lab=c(0.1,1,10,100),cex=txt.size,tck=-0.02,xpd=T,las=1)
  text(x=bar, y=-0.02,labels=indat[,2], srt=rot, adj=1, xpd=TRUE,tcl=-0.02,cex=txt.size)
  yl<-c(0.1,1,10,100)*genes.total/100
  yl[yl<1]<-1
  yl<-round(yl)
  axis(side = 4,at=c(0,1,2,3),lab=yl,cex=txt.size,tck=-0.02,xpd=T,las=1)
  
  xdis<-strwidth(indat[which.max(nchar(indat[,2])),2])*cos(rot*pi/180)*txt.size
  yf<-tan(rot*pi/180)/strwidth("a")*strwidth("a","in")/strheight("a","in")*strheight("a")
  ydis<-xdis*yf+0.05
  for(i in 1:length(group)){
    lab<-class==group[i]
    l<-indat[lab,2]
    x1<-bar[lab][1]-strwidth(l[1])*cos(rot*pi/180)*txt.size
    y1<--yf*strwidth(l[1])*cos(rot*pi/180)*txt.size
    x2<-tail(bar[lab],1)-strwidth(tail(l,1))*cos(rot*pi/180)*txt.size
    y2<--yf*strwidth(tail(l,1))*cos(rot*pi/180)*txt.size
    segments(x0 = x1,y0 = y1-0.06,x1 =bar[lab][1]-xdis ,y1=-ydis-0.06,xpd=T)
    segments(x0=tail(bar[lab],1)-xdis,y0=-ydis-0.06,x1 =bar[lab][1]-xdis ,y1=-ydis-0.06,xpd=T)
    segments(x0=tail(bar[lab],1)-xdis,y0=-ydis-0.06,x1 =x2 ,y1=y2-0.06,xpd=T)    
    text(labels = group[i],x =mean(range(bar[lab]))-xdis,y=-ydis-0.1,xpd=T,cex=txt.size,adj =c(0.5,1) )
  }
  dev.off()
  
}else{
 stop("data trans method is wrong!")
}
