#!/share/nas2/genome/biosoft/R/3.1.1/lib64/R/bin/Rscript

library('getopt');

#get options, using the spec as defined by the enclosed list.
#we read the options from the default: commandArgs(TRUE).

spec = matrix(c(
	'help' , 'h', 0, "logical",
	'infile' , 'i', 1, "character",
	'outdir' , 'o', 1, "character",
	'method' , 'r', 1, "character",
	'height' , 'H', 1, "integer",
	'width' , 'W', 1, "integer",
	'show_rownames' , 'a', 0, "logical",
	'show_colnames' , 'b', 0, "logical",
	'cluster_rows' , 'c', 0, "logical",
	'cluster_cols' , 'd', 0, "logical",
	'treeheight_row' , 'e', 1, "integer",
	'treeheight_col' , 'g', 1, "integer",
	'legend' , 'j', 0, "logical",
	'fontsize' , 'k', 1, "integer",
	'fontsize_row' , 'm', 1, "integer",
	'fontsize_col' , 'n', 1, "integer",
	'cellwidth' , 'p', 1, "integer",
	'cellheight' , 'q', 1, "integer",
	'display_number', 't', '1', "integer",
	'scale' , 's', 1, "character"
	),byrow=TRUE,ncol=4);
opt = getopt(spec);

# define usage function
print_usage <- function(spec=NULL){
	cat(getopt(spec, usage=TRUE));
	cat("Usage example: \n")
	cat("
1) Rscript fpkm_cor.r --infile All_gene_fpkm.list --outdir density/

Options:
--help          -h      NULL            get this help
--infile                character       the input file [forced]
--outdir		character       the filename for output graph [forced]
--method		character	cor method [optional(pearson/spearman), default pearson]
--height		integer         the height of graph [optional, default: 3000]
--width			integer         the width of graph [optional, default: 3000]
--show_rownames         logical         show rownames [optional, default: FALSE]
--show_colnames         logical         show colnames [optional, default: FALSE]
--cluster_rows          logical         show row cluster tree [optional, default: FALSE]
--cluster_cols          logical         show col cluster tree [optional, default: FALSE]
--treeheight_row        integer         the height of row cluster tree [optional, default: 50]
--treeheight_col        integer         the height of col cluster tree [optional, default: 50]
--legend                logical         show legend [optional, default: FALSE]
--fontsize              integer         the size of font [optional, default: NULL]
--fontsize_row          integer         the size of font for rowname [optional, default: NULL]
--fontsize_col          integer         the size of font for colname [optional, default: NULL]
--cellwidth             integer         the width of cell [optional, default: NA]
--cellheight            integer         the height of cell [optional, default: NA]
--display_number	integer		display number on plot [optional, default: T]
--scale                 character       row or column or none [optional, default: none]
\n")
	q(status=1);
}

# if help was asked for print a friendly message
# and exit with a non-zero error code
if ( !is.null(opt$help) ) { print_usage(spec) }
# check non-null args
if ( is.null(opt$infile) )      { print_usage(spec) }
if ( is.null(opt$outdir) )      { print_usage(spec) }

#set some reasonable defaults for the options that are needed,
#but were not specified.
if ( is.null(opt$method ) )		{ opt$method = "pearson" }
if ( is.null(opt$height ) )             { opt$height = 3000 }
if ( is.null(opt$width ) )              { opt$width = 3000 }
if ( is.null(opt$show_rownames ) )      { opt$show_rownames = FALSE }
if ( is.null(opt$show_colnames ) )      { opt$show_colnames = FALSE }
if ( is.null(opt$cluster_rows ) )       { opt$cluster_rows = FALSE }
if ( is.null(opt$cluster_cols ) )       { opt$cluster_cols = FALSE }
if ( is.null(opt$treeheight_row ) )     { opt$treeheight_row = ifelse(opt$cluster_rows, 50, 0) }
if ( is.null(opt$treeheight_col ) )     { opt$treeheight_col = ifelse(opt$cluster_cols, 50, 0) }
if ( is.null(opt$legend) )              { opt$legend = FALSE }
if ( is.null(opt$fontsize ) )           { opt$fontsize = 10 }
if ( is.null(opt$fontsize_row ) )       { opt$fontsize_row = NULL }
if ( is.null(opt$fontsize_col ) )       { opt$fontsize_col = NULL }
if ( is.null(opt$cellwidth ) )          { opt$cellwidth = NA }
if ( is.null(opt$cellheight ) )         { opt$cellheight = NA }
if ( is.null(opt$display_number ) )	{ opt$display_number = TRUE }
if ( is.null(opt$scale) )               { opt$scale = "none" }

# check scale
if( (opt$scale != "column") && (opt$scale != "row") && (opt$scale != "none") ){
        cat("Final Error: the scale must be column or row or none\n")
        print_usage(spec)
}

library(Hmisc)
library(corrplot)
library(pheatmap)
##read fpkm
fpkm_data <- read.delim(opt$infile, row.names = 1, header=TRUE,check.names =F)
fpkm<-as.matrix(fpkm_data)

##filter low exp genes
which(rowSums(fpkm)==0)->low_exp
if(length(low_exp)>0){
	fpkm<-fpkm[-low_exp,]
}

cor<-rcorr(fpkm,type=opt$method)

write.table(cor$P,paste(opt$outdir,"/sample_pvalue.txt",sep=""),row.names=T,col.names=T,sep="\t",quote=F)
write.table(cor$r,paste(opt$outdir,"/sample_coefficient.txt",sep=""),row.names=T,col.names=T,sep="\t",quote=F)

#########font size
sample_num <- dim(fpkm)[2]
print(sample_num)
percent =0.8
if (sample_num<=15){
	percent = 0.8
}else if(sample_num>15 && sample_num<=20){
	percent = 0.7
}else if(sample_num>20 && sample_num<=30){
	percent = 0.5
}else if(sample_num>30 && sample_num<=40){
	percent = 0.2
}else if(sample_num>40 && sample_num<=50){
	opt$fontsize= 8
	percent = 0.2
}else if(sample_num>50 && sample_num<=100){
	opt$fontsize = 2
}else{
	opt$fontsize = 2
	opt$display_number = FALSE
}

#########plot
file <- paste(opt$outdir, "/sample_corelation.png", sep="")
png(filename=file, height=opt$height, width=opt$width, res=500, units="px")
#corrplot(cor$r, type="upper", order="hclust", tl.col="black", tl.srt=45)
pheatmap(cor$r,color=colorRampPalette(c("orchid1","white","cyan"))(100),scale=opt$scale,margins=c(10,10),
	display_numbers=opt$display_number,fontsize_number=percent*opt$fontsize,number_format = "%.2f",number_color = "grey30",
	show_rownames=opt$show_rownames, show_colnames=opt$show_colnames,
	cluster_rows=opt$cluster_rows, cluster_cols=opt$cluster_cols,
	treeheight_row=opt$treeheight_row, treeheight_col=opt$treeheight_col,
	legend=opt$legend, fontsize=opt$fontsize,
	fontsize_row=opt$fontsize_row, fontsize_col=opt$fontsize_col,
	cellwidth=opt$cellwidth, cellheight=opt$cellheight )
dev.off()

file <- paste(opt$outdir, "/sample_corelation.pdf", sep="")
pheatmap(filename = file,cor$r,color=colorRampPalette(c("orchid1","white","cyan"))(100),scale=opt$scale,margins=c(10,10),
	display_numbers=opt$display_num,fontsize_number=percent*opt$fontsize,number_format = "%.2f",number_color = "grey30",
	show_rownames=opt$show_rownames, show_colnames=opt$show_colnames,
	cluster_rows=opt$cluster_rows, cluster_cols=opt$cluster_cols,
	treeheight_row=opt$treeheight_row, treeheight_col=opt$treeheight_col,
	legend=opt$legend, fontsize=opt$fontsize,
	fontsize_row=opt$fontsize_row, fontsize_col=opt$fontsize_col,
	cellwidth=opt$cellwidth, cellheight=opt$cellheight )
dev.off()
#corrplot(cor$r, type="upper", order="hclust", tl.col="black", tl.srt=45)


