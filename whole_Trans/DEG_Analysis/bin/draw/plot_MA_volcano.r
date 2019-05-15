### load library
library(getopt)
spec=matrix(c(
        'help','h',0,"logical",
        'all','i',1,"character",
        'FDR','F',1,"character",
        'PValue','P',1,"character",
        'FC','f',1,"character",
	'od','o',1,"character",
	'type','t',1,"character"
        ),byrow=TRUE,ncol=4
)
opt=getopt(spec)
print_usage<-function(spec=NULL){
        cat(getopt(spec,usage=TRUE));
        cat("Usage example: \n")
        cat("
                Description: 
                Contact: Yanhua Wen <wenyh@biomarker>
                Usage example:
                
		Options:
                --help          NULL            get this help
                --all		character       DEG.all
                --FDR		character       FDR cutoff
                --PValue	character       PValue cutoff
						One of FDR and PValue must be existed!
						If given both parameter, FDR would be used!
                --FC		character       FC, default 2
		--od		character	output prefix
		--type		character	FPKM/TPM
                \n")

        q(status=1);
}
if(is.null(opt$all))	print_usage(spec)
if(is.null(opt$od))    print_usage(spec)
if(is.null(opt$FC))	opt$FC=2
if(is.null(opt$type))	opt$type="FPKM"

flag="FDR"
cutoff=1
opt$FC=as.numeric(as.character(opt$FC))
if(is.null(opt$FDR) && is.null(opt$PValue)){
        print_usage(spec)
}
if(!is.null(opt$PValue)){
        flag="PValue"
        cutoff=as.numeric(as.character(opt$PValue))
}
if(!is.null(opt$FDR)){
        flag="FDR"
        cutoff=as.numeric(as.character(opt$FDR))
}

require(ggplot2)
#################################################################
### Drawing MA plot for differential expression analysis.
### For each gene: the log2(fold change) between the two samples is plotted (y axis)
### against the gene's log2(average expression) in the two samples(x axis)
### NOTE: here using FPKM stand for expression
#################################################################
plot_MA <- function(log2exp=NULL, log2FC=NULL, FDR=NULL, Significant=NULL,
xlab="log2(FPKM)", ylab="log2(FC)", main="MA plot") {
	# check args
	# check null
	if( is.null(log2FC) ) stop("log2FC is NULL")
	if( is.null(FDR) ) stop("FDR is NULL")
	if( is.null(Significant) ) stop("Significant is NULL")
	# check length
	len <- c(length(log2exp), length(log2FC), length(FDR))
	if( len[1] != len[2] ) 
		stop(paste("length(log2exp) != length(log2FC): ", len[1], " != ", len[2], sep=""))
	if( len[2] != len[3] ) 
		stop(paste("length(log2FC) != length(FDR): ", len[2], " != ", len[3], sep=""))
	if( len[1] == 0 ) stop("length(log2FC) == 0")

	# plot
	Significant<-factor(Significant,levels=c("Up","Down","Normal"))
	p <- qplot(log2exp, log2FC, xlab=xlab, ylab=ylab, main=main, size=I(0.7), colour=Significant)
	p <- p+ scale_color_manual(values = c("Up"="red","Normal"="black","Down"="green"))+ theme_bw()+
		theme( panel.background = element_rect(colour="black", size=1, fill="white"), panel.grid = element_blank())
	# return
	return(p)
}


#################################################################
### Drawing Volcano plot for differential expression analysis.
### For each gene: the log2(fold change) between the two samples is plotted (x axis)
### against the gene's -log10(FDR) in the two samples(y axis)
#################################################################
plot_Volcano <- function(log2FC=NULL, FDR=NULL, Significant=NULL, 
xlab="log2(FC)", ylab="-log10(FDR)", main="Volcano plot",yline=NULL,xline=NULL) {
	# check args
	# check null
	Significant=as.vector(Significant)
        Significant=factor(Significant,levels=unique(Significant))
	if( is.null(log2FC) ) stop("log2FC is NULL")
	if( is.null(FDR) ) stop("FDR is NULL")
	if( is.null(Significant) ) stop("Significant is NULL")
	if( is.null(yline) ) stop("yline is NULL")
	if( is.null(xline) ) stop("xline is NULL")
	# check length
	len <- c(length(log2FC), length(FDR))
	if( len[1] != len[2] ) 
		stop(paste("length(log2FC) != length(log2FDR): ", len[1], " != ", len[2], sep=""))
	if( len[1] == 0 ) stop("length(log2FC) == 0")

	# plot
	Significant<-factor(Significant,levels=c("Up","Down","Normal"))
	p <- qplot(log2FC, -log10(FDR), xlab=xlab, ylab=ylab, main=main, size=I(0.7), colour=Significant)
	p <- p+ scale_color_manual(values = c("Up"="red","Normal"="black","Down"="green"))
	p <- p+geom_vline(xintercept=xline,lty=2,size=I(0.2),colour="grey11")
	p <- p+geom_hline(yintercept=yline,lty=2,size=I(0.2),colour="grey11")
	p <- p+theme_bw() + theme( panel.background = element_rect(colour="black", size=1, fill="white"),
		panel.grid = element_blank())
	# return
	return(p)
}


read.csv(opt$all,header=T,sep="\t",row.names=1)->all
exp_column<-grep(pattern=opt$type,colnames(all))
if(length(exp_column)==0){
	print(paste(opt$all," does not have columns contained ",opt$type,sep=""))
	q(status=1)	
}
exp<-as.matrix(all[,exp_column])
exp<-log(rowMeans(exp),2)


log2FC<-all$log2FC
logfc<-log(opt$FC,2)
which(colnames(all)==flag)->s
FDR<-all[,s]

#"Up"="red","Normal"="black","Down"
Significant<-rep("Normal",length(FDR))

which(FDR<=cutoff)->s1
which(log2FC>=logfc)->fc1
which(log2FC<= -logfc)->fc2

up<-intersect(s1,fc1)
down<-intersect(s1,fc2)

Significant[up]="Up"
Significant[down]="Down"
####################################plot MA

ma_xlab<-paste("log2(",opt$type,")",sep="")
png(paste(opt$od,"_MA.png",sep=""),width=3000,height=3000,res=500)
plot_MA(log2exp=exp, log2FC=log2FC, FDR=FDR, Significant=Significant,xlab=ma_xlab, ylab="log2(FC)", main="MA plot")
dev.off()
pdf(paste(opt$od,"_MA.pdf",sep=""))
plot_MA(log2exp=exp, log2FC=log2FC, FDR=FDR, Significant=Significant,xlab=ma_xlab, ylab="log2(FC)", main="MA plot")
dev.off()

####################################plot volcano
vol_ylab=paste("-log10(",flag,")",sep="")
png(paste(opt$od,"_Volcano.png",sep=""),width=3000,height=3000,res=500)
plot_Volcano(log2FC=log2FC, FDR=FDR, Significant=Significant,xlab="log2(FC)", ylab=vol_ylab, main="Volcano plot",xline=c(-logfc,logfc),yline=-log(cutoff,10))
dev.off()
pdf(paste(opt$od,"_Volcano.pdf",sep=""))
plot_Volcano(log2FC=log2FC, FDR=FDR, Significant=Significant,xlab="log2(FC)", ylab=vol_ylab, main="Volcano plot",yline=-log(cutoff,10),xline=c(-logfc,logfc))
dev.off()
