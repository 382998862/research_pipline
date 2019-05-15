#!/share/nas2/genome/biosoft/R/2.15.1/lib64/R/bin/Rscript 
library('getopt');


#-----------------------------------------------------------------
# getting parameters
#-----------------------------------------------------------------
#get options, using the spec as defined by the enclosed list.
#we read the options from the default: commandArgs(TRUE).
spec = matrix(c(
				'help' , 'h', 0, "logical",
				'infile' , 'i', 1, "character",
				'outfile' , 'o', 1, "character",
				'outdir', 'd', 1, "character",
				'term.col' , 't', 1, "integer",
				'ontology.col' , 'O', 1, "integer",
				'items_num.col' , 'I', 1, "integer",
				'pvalue.col' , 'p', 1, "integer",
				'height' , 'H', 1, "integer",
				'width' , 'W', 1, "integer",
				'toplines' , 'l', 1, "integer",
				'x.lab' , 'X', 1, "character",
				'y.lab' , 'Y', 1, "character",
				'title.lab' , 'T', 1, "character",
				'p_threhold' , 'P', 1, "double",
				'header' , 'k', 0, "logical",
				'items_num_size','s',1,"double",
				'lab.size' , 'L', 1, "integer",
				'axis.size' , 'a', 1, "integer",
				'legend.size' , 'D', 1, "integer"
		), byrow=TRUE, ncol=4);
opt = getopt(spec);


# define usage function
print_usage <- function(spec=NULL){
	cat(getopt(spec, usage=TRUE));
	cat("\n\n")
	cat("
				data format:
					Data format  like below: 
					columns split by '\\t',one column must include one type data,four columns must be telling this script 
					which represent term,ontology,Items and pvalue, respectively;
					Note:if '#' in the front of line, this line will be ignored;
						column of Items could be integer or type of '12 out of 440'

					GO_ID	GO_Term	Ontology	Conserved_Items	Background_Items	P-value	Adjust_P-value
					GO:0003743	translation initiation factor activity	Molecular Function	28	66	5.45E-13	2.95E-10
					GO:0016592	mediator complex	Cellular Component	17	23	7.12E-11	2.15E-08

					
				Usage example: 
					1) Rscript GO_bar.R --infile Go.txt --outfile bbb --term.col 2 --ontology.col 3 --items_num.col 5 --pvalue.col 6 --header --width 6000
					
					
				Options: 
					--help		-h 	NULL 		get this help
					--infile 	-i 	character 	the input file [forced]
					--outfile 	-o 	character 	the filename prefix for output graph [forced]
					--term.col 	-t 	integer 	the term column to draw for y value[forced]
					--ontology.col 	-O 	integer 	the ontology column  [forced]
					--items_num.col 	-y 	integer 	the items_num column for bar label [forced]
					--pvalue.col	-p	the pvalue column to draw for x value[forced]
					--outdir	-d	output file dir [optional, default: current dir]
					--p_threhold	-p	the threhold of pvalue to draw [optional, default: null]
					--height 	-H 	integer 	the height of graph [optional, default: 3000]
					--width 	-W 	integer 	the width of graph [optional, default: 6000]
					--toplines 	-l 	integer 	the top enrich lines to draw [ optional,default of 20]
					--header	-k	NULL	infile whether has a header [[optional, default: NULL]]
					--x.lab 	-X 	character 	the lab for x [optional,default of '-log10(P-value)']
					--y.lab 	-Y 	character 	the lab for y [optional,default of 'GO term']
					--title.lab 	-T 	character 	the lab for title [optional, default: 'The Most enriched GO Terms']
					--lab.size 	-l 	integer 	the font size of lab [optional, default: 14]
					--items_num_size -s	double	the font size of num on the bar [optional, default: 3] 
					--lab.size 	-L 	integer 	the font size of lab [optional, default: 14]
					--axis.size 	-a 	integer 	the font size of text for axis [optional, default: 14]
					--legend.size 	-D 	integer 	the font size of text for legend [optional, default: 12]
				Author:
					huangls
				Version:
					1.1;2014-10-17
					\n")
	q(status=1);
}



# if help was asked for print a friendly message
# and exit with a non-zero error code
if ( !is.null(opt$help) ) { print_usage(spec) }
if ( is.null(opt$infile) )	{ print_usage(spec) }
if ( is.null(opt$outfile) )	{ print_usage(spec) }
if ( is.null(opt$term.col) )	{ print_usage(spec) }
if ( is.null(opt$ontology.col) )	{ print_usage(spec) }
if ( is.null(opt$items_num.col) )	{ print_usage(spec) }
if ( is.null(opt$pvalue.col) )	{ print_usage(spec) }



#set some reasonable defaults for the options that are needed,
#but were not specified.
if ( is.null(opt$y.lab) )	{ opt$y.lab='GO term' }
if ( is.null(opt$x.lab) )	{ opt$x.lab=paste('-log10(','P-value',')',sep="") }else{opt$x.lab=paste('-log10(',opt$x.lab,')',sep="")}
if ( is.null(opt$height ) )		{ opt$height = 3000 }
if ( is.null(opt$width ) )		{ opt$width = 6000 }
if ( is.null(opt$lab.size ) )		{ opt$lab.size = 14 }
if ( is.null(opt$axis.size ) )		{ opt$axis.size = 9 }
if ( is.null(opt$title.lab) )		{ opt$title.lab = 'The Most enriched GO Terms' }
if ( is.null(opt$items_num_size) )	{ opt$items_num_size = 3 }

#if ( is.null(opt$is.point) )		{ opt$is.point = FALSE }
#if ( is.null(opt$is.shape) )		{ opt$is.shape = FALSE }
#if ( is.null(opt$is.line) )		{ opt$is.line = FALSE }
#if ( is.null(opt$is.linetype) )		{ opt$is.linetype = FALSE }
#if ( is.null(opt$legend.xpos ) )	{ opt$legend.xpos = NULL }
#if ( is.null(opt$legend.ypos ) )	{ opt$legend.ypos = NULL }
#if ( is.null(opt$legend.col ) )		{ opt$legend.col = NULL }
if ( is.null(opt$legend.size ) )	{ opt$legend.size = 12 }
#if ( is.null(opt$line.size) )		{ opt$line.size = 0.8 }
#if ( is.null(opt$point.size) )		{ opt$point.size = 2 }

# check combination non-null args

# load library

require(ggplot2)

# check legend position
if( (!is.null(opt$legend.xpos)) && (!is.null(opt$legend.ypos)) ){
	# check xpos
	if( (opt$legend.xpos < 0) || (opt$legend.xpos > 1) ){
		cat("Final Error: legend.xpos must be in [0.0, 1.0]")
		print_usage(spec)
	}
	# check ypos
	if( (opt$legend.ypos < 0) || (opt$legend.ypos > 1) ){
		cat("Final Error: legend.ypos must be in [0.0, 1.0]")
		print_usage(spec)
	}
}
get_gene_number<-function(d){
	len=length(d)
	new_name<-c(rep=0,times=len)
	for (i in 1:len){
		
		new_name[i]<-as.numeric(unlist(strsplit(as.character(d[i])," ",perl=T))[1])
	}
	return(new_name)
}
get_threhold<-function(d,fraction){
	
	num<-floor(nrow(d)*fraction)
	d.sorted<-d[order(d$pvalue),]
	d.sorted<-d.sorted$pvalue
	res<-c(d.sorted[num],d.sorted[nrow(d)-num])
	return(res)
}
change_name<-function(old_name){
	len<-length(old_name)
	new_name<-c(rep=0,times=len)
	for (i in 1:len){
		tmp<-paste(' ',old_name[i],sep="")
		new_name[i]<-tmp
	}
	return(new_name)
}
#-----------------------------------------------------------------
# reading data
#-----------------------------------------------------------------
# reading data
if (!is.null(opt$header)){
	data <- read.table(opt$infile, header=T,sep='\t',quote='')
}else{
	data <- read.table(opt$infile, header=F,sep='\t',quote='')
}


df <- data.frame(term=data[,opt$term.col], ontology=data[,opt$ontology.col], items_num=data[,opt$items_num.col],pvalue=data[,opt$pvalue])
# order data
df[,4]<-as.double(df[,4])
df <- df[order(df[,4]),]
if (!is.null(opt$toplines) ){
	df<-head(df,opt$toplines)
}else{
	if (!is.null(opt$p_threhold)){
		df<-subset(df,pvalue<opt$p_threhold)
	}else{
		df<-head(df,20)
	}
	
}


df

class_ontology<-unique(df$ontology)

df_classed<-data.frame()

for (i in class_ontology){
	
	tmp<-subset(df,df$ontology==i)
	df_classed<-rbind(df_classed,tmp) 
	
}

colnames(df_classed)<-c('term','ontology','items_num','pvalue')

#df_classed$term<-change_name(df_classed$term)
#keep the order of  occurrence in the data
df_classed$term <- factor(df_classed$term, levels = rev(df_classed$term))
#df_classed$ontology <- factor(df_classed$ontology, levels = df_classed$term)
#class_ontology<-factor(class_ontology,levels=rev(class_ontology))

#-----------------------------------------------------------------
# plot
#-----------------------------------------------------------------
# mian plot
#write.table(df_classed,file='go_out.txt',sep='\t',quote=F,row.names=F)
#rownames(df_classed)<-seq(from=1,to=nrow(df_classed))

#df_classed[df_classed$pvalue==0,4]<-as.double(1e-100)

if (!is.numeric(df_classed$items_num)){df_classed$items_num<-get_gene_number(df_classed$items_num)}


df_classed
p<-ggplot(data=df_classed, aes(x=term, y=-log10(pvalue),fill=ontology))+
				geom_bar(stat="identity")+coord_flip()+
				xlab(opt$y.lab) + ylab(opt$x.lab) + 
				labs(title=opt$title.lab)+
				guides(col = guide_legend(reverse = TRUE))+
				scale_fill_brewer(palette="Set2",guide = guide_legend(title = NULL))+
				geom_text(aes(label=items_num),color='black',size=opt$items_num_size,hjust=-0.25)

#p <- p + theme(title = element_text(face="bold", size=opt$lab.size), axis.text = element_text(face="bold", size=opt$axis.size))
#p <- p + theme(legend.title = element_text(face="bold", size=opt$legend.size),legend.text = element_text(size=opt$legend.size) )
p <- p + theme(title = element_text(size=opt$lab.size), axis.text = element_text(size=opt$axis.size))
p <- p + theme(legend.title = element_text(size=opt$legend.size),legend.text = element_text(size=opt$legend.size) )

#-----------------------------------------------------------------
# output plot
#-----------------------------------------------------------------



if ( is.null(opt$outdir))	{
	#summary(p)
	
	pdf(file=paste(opt$outfile,".pdf",sep=""), height=opt$height*2/1000, width=opt$width*2/1000)
	print(p)
	dev.off()
	png(filename=paste(opt$outfile,".png",sep=""), height=opt$height, width=opt$width, res=500, units="px")
	print(p)
	dev.off()
	tiff(filename=paste(opt$outfile,".tiff",sep=""), height=opt$height, width=opt$width, res=500, units="px")
	print(p)
	dev.off()
	write.table(as.matrix(df_classed),file=paste(opt$outfile,".list",sep=""),quote=F,sep='\t',row.names=F,col.names=T)

	
}else{
	#summary(p)
	if(substr(opt$outdir,nchar(opt$outdir),nchar(opt$outdir))=="/"){
		opt$outdir<-substr(opt$outdir,1,nchar(opt$outdir)-1)
	}
	
	#setwd(opt$outdir)
	pdf(file=paste(opt$outdir,"/",opt$outfile,".pdf",sep=""), height=opt$height*2/1000, width=opt$width*2/1000)
	print(p)
	dev.off()
	png(filename=paste(opt$outdir,"/",opt$outfile,".png",sep=""), height=opt$height, width=opt$width, res=500, units="px")
	print(p)
	dev.off()
	tiff(filename=paste(opt$outdir,"/",opt$outfile,".tiff",sep=""), height=opt$height, width=opt$width, res=500, units="px")
	print(p)
	dev.off()
	
	write.table(as.matrix(df_classed),file=paste(opt$outdir,"/",opt$outfile,".list",sep=""),quote=F,sep='\t',row.names=F,col.names=T)
	#setwd(curr_dir)
}
q(save='no')





