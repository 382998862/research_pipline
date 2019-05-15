#!/share/nas2/genome/biosoft/R/2.15.1/lib64/R/bin/Rscript 
# load library
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
				'pathway_term.col' , 't', 1, "integer",
				'rich_factor.col' , 'O', 1, "integer",
				'qvalue.col' , 'I', 1, "integer",
				'gene_number.col' , 'p', 1, "integer",
				'outdir', 'd', 1, "character",
				'height' , 'H', 1, "integer",
				'width' , 'W', 1, "integer",
				'toplines' , 'l', 1, "integer",
				'x.lab' , 'X', 1, "character",
				'y.lab' , 'Y', 1, "character",
				'title.lab' , 'T', 1, "character",
				'q_threhold' , 'q', 1, "double",
				'header' , 'k', 0, "logical",
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
					which represent pathway	gene_number	enrichment_factor	and correct_p, respectively;
					Note:if '#' in the front of line, this line will be ignored;
						column of gene_nubher could be integer or type of '12 out of 440'

					pathway	gene_number	enrichment_factor	correct_p
					Starch and sucrose metabolism	1 out of 330	0.33	1.24E-07
					Galactose metabolism	2 out of 440	0.23	1.05E-05
					
					
					
					Usage example: 
					1) Rscript KEGG_dot.R --infile Kegg_pathway.txt --outfile qqq --pathway_term.col 1 --rich_factor.col 3 --qvalue.col 4 --gene_number.col 2 --header
					
					
					Options: 
					--help		-h 	NULL 		get this help
					--infile 	-i 	character 	the input file [forced]
					--outfile 	-o 	character 	the filename prefix for output graph [forced]
					--pathway_term.col 	-t 	integer 	the term column to draw for y axis[forced]
					--rich_factor.col 	-O 	integer 	the rich factor column  to draw for x axis [forced]
					--qvalue.col 	-y 	integer 	the qvalue(corrected pvalue) column  [forced]
					--gene_number.col	-p	the gene number column [forced]
					--q_threhold	-q	the threhold of pvalue to draw [optional, default: null]
					--outdir	-d	output file dir [optional, default: current dir]
					--height 	-H 	integer 	the height of graph [optional, default: 3000]
					--width 	-W 	integer 	the width of graph [optional, default: 4000]
					--toplines 	-l 	integer 	the top enrich lines to draw [ optional,default of 20]
					--header	-k	NULL	infile whether has a header [[optional, default: NULL]]
					--x.lab 	-X 	character 	the lab for x [optional,default of 'Rich factor']
					--y.lab 	-Y 	character 	the lab for y [optional,default of '']
					--title.lab 	-T 	character 	the lab for title [optional, default: 'Statistics of Pathway Enrichment']
					--lab.size 	-l 	integer 	the font size of lab [optional, default: 14]
					--axis.size 	-a 	integer 	the font size of text for axis [optional, default: 14]
					--legend.size 	-D 	integer 	the font size of text for legend [optional, default: 12]
					Author:
						huangls
					Version:
						1.0;2014-10-17
					\n")
	q(status=1);
}
#srcdir <- dirname(get_Rscript_filename()) ## get_Rscript_filename() is a method in package getopt
#source(paste(srcdir, "/", "BSseq.util.r", sep=""))
#source(paste(srcdir, "/", "methylplot.util.r", sep=""))



# if help was asked for print a friendly message
# and exit with a non-zero error code
if ( !is.null(opt$help) ) { print_usage(spec) }
if ( is.null(opt$infile) )	{ print_usage(spec) }
if ( is.null(opt$outfile) )	{ print_usage(spec) }
if ( is.null(opt$pathway_term.col) )	{ print_usage(spec) }
if ( is.null(opt$rich_factor.col) )	{ print_usage(spec) }
if ( is.null(opt$qvalue.col) )	{ print_usage(spec) }
if ( is.null(opt$gene_number.col) )	{ print_usage(spec) }



#set some reasonable defaults for the options that are needed,
#but were not specified.
if ( is.null(opt$y.lab) )	{ opt$y.lab="" }
if ( is.null(opt$x.lab) )	{ opt$x.lab='Rich factor' }
if ( is.null(opt$height ) )		{ opt$height = 3000 }
if ( is.null(opt$width ) )		{ opt$width = 4000 }
if ( is.null(opt$lab.size ) )		{ opt$lab.size = 14 }
if ( is.null(opt$axis.size ) )		{ opt$axis.size = 9 }
if ( is.null(opt$title.lab) )		{ opt$title.lab = "Statistics of Pathway Enrichment" }
if ( is.null(opt$legend.size ) )	{ opt$legend.size = 12 }


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

get_threhold<-function(d,fraction){
	
	num<-floor(nrow(d)*fraction)
	d.sorted<-d[order(d$qvalue),]
	d.sorted<-d.sorted$qvalue
	res<-c(d.sorted[num],d.sorted[nrow(d)-num])
	return(res)
}
get_gene_number<-function(d){
	len=length(d)
	new_name<-c(rep=0)
	for (i in 1:len){
		
		new_name[i]<-as.numeric(unlist(strsplit(as.character(d[i])," ",perl=T))[1])
	}
	return(new_name)
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
# load library
require(ggplot2, quietly=TRUE)

# reading data
read.csv(opt$infile,header=T,sep="\t")->df
#print OUT "ID\tDescription\tGeneRatio\tBgRatio\tenrich_factor\tpvalue\tqvalue\tgeneID\tContained\n";

colnames(df)<-c("ID","pathway_term","geneRatio","bgRatio","rich_factor","pvalue","qvalue","geneID","gene_number","Diff")
#pathway_term    rich_factor     qvalue  gene_number     shape   group
if (!is.numeric(df$gene_number)){df$gene_number<-get_gene_number(df$gene_number)}
#-----------------------------------------------------------------
# plot
#-----------------------------------------------------------------
# mian plot
#kegg pathway plot

#p<-ggplot(df, aes(rich_factor,pathway_term))
p<-ggplot(df, aes(rich_factor,y=reorder(pathway_term, qvalue) ) )
p<-p+geom_point(aes(colour=qvalue,size=gene_number,shape=Diff))+scale_colour_gradientn(colours=rainbow(4),guide = "colourbar") +expand_limits(color=seq(0, 1, by=0.25))
p<-p+ggtitle(opt$title.lab) + xlab(opt$x.lab) +ylab(opt$y.lab)+theme_bw()
p<-p+theme(panel.border=element_rect(colour = "black"))
p<-p+theme(plot.title=element_text(vjust=1), legend.key=element_blank())

p <- p + theme(title = element_text(face="bold", size=opt$lab.size), axis.text.x = element_text(face="bold", size=opt$axis.size), axis.text.y=element_text(face="bold", size=opt$axis.size))
p <- p + theme(legend.title = element_text(face="bold", size=opt$legend.size),legend.text = element_text(size=opt$legend.size))

#-----------------------------------------------------------------
# output plot
#-----------------------------------------------------------------


if ( is.null(opt$outdir))	{
	png(filename=paste(opt$outfile,".png",sep=""), height=opt$height, width=opt$width, res=opt$height/6, units="px")
	print(p)
	dev.off()
	tiff(filename=paste(opt$outfile,".tiff",sep=""), height=opt$height, width=opt$width, res=opt$height/6, units="px")
	print(p)
	dev.off()
	write.table(as.matrix(df),file=paste(opt$outfile,".list",sep=""),quote=F,sep='\t',row.names=F,col.names=T)
}else{
    if(!file.exists(opt$outdir)){
        if(!dir.create(opt$outdir,showWarnings = FALSE, recursive = TRUE)){
            stop(paste("dir.create failed: outdir=",opt$outdir,sep=""))
        }
    }

	if(substr(opt$outdir,nchar(opt$outdir),nchar(opt$outdir))=="/"){
		opt$outdir<-substr(opt$outdir,1,nchar(opt$outdir)-1)
	}
	
	png(filename=paste(opt$outdir,"/",opt$outfile,".png",sep=""), height=opt$height, width=opt$width, res=opt$height/6, units="px")
	print(p)
	dev.off()
	tiff(filename=paste(opt$outdir,"/",opt$outfile,".tiff",sep=""), height=opt$height, width=opt$width, res=opt$height/6, units="px")
	print(p)
	dev.off()
}
q(save='no')





