args<-commandArgs(TRUE)
if(length(args)!=2){
	print("1. input file")
	print("2. output prefix")
}
############## main codes ###############
data<-read.csv(args[1],sep="\t",header=T,row.names=1) 
data<-as.matrix(data)

rnatype<-colnames(data)
chr<-rownames(data)
sum<-0
for(i in 1:length(chr)){sum=sum + nchar(chr[i])}
mean<-sum/length(chr)

colors <- c("green","yellow","purple","blue")
colors<-colors[1:length(rnatype)]

file<-paste(args[2],".png",sep="")
png(file,width=1440*dim(data)[1]/30,height=480)
position<-barplot(t(data), beside = TRUE,main="Different RNA expression Ratio in different chromosomes",ylim=c(0,max(data)*1.3),ylab="Ratio",col=colors,axisnames=F)
legend("topright", colnames(data), cex=1.0, fill=colors,bty="n")
x <- (position[1,] + position[2,])/2+2
y <- par("usr")[3]-0.05
if(mean>3){
	text(x, y, font = 10,pos=2,labels=rownames(data), adj=1, srt=45, xpd=TRUE)
}
if(mean<=3){
	text(x, y, font = 10,pos=2,labels=rownames(data), adj=1, xpd=TRUE)
}
box()
dev.off()



file<-paste(args[2],".pdf",sep="")
pdf(file,width=21*dim(data)[1]/30,height=7)
position<-barplot(t(data), beside = TRUE,main="Different RNA expression Ratio in different chromosomes",ylim=c(0,max(data)*1.3),ylab="Ratio",col=colors,axisnames=F)
legend("topright", colnames(data), cex=1.0, fill=colors,bty="n")
x <- (position[1,] + position[2,])/2+2
y <- par("usr")[3]-0.05
if(mean>3){
        text(x, y, font = 10,pos=2,labels=rownames(data), adj=1, srt=45, xpd=TRUE)
}
if(mean<=3){
        text(x, y, font = 10,pos=2,labels=rownames(data), adj=1, xpd=TRUE)
}
box()
dev.off()

