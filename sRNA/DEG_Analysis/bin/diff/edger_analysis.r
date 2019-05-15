library('getopt')

spec = matrix(c(
	'help', 'h', '0', "logical",
	'count', 'i', '1', "character",
	'config', 'c', '1', "character",
	'out.prefix', 'o', '1', "character",
	'filter', 'f', '1', "character",
	'fpkm', 'e', '2', "character",
	'type', 't', '2', "character",
	'disper', 'd', '2', "character"
), byrow=TRUE, ncol=4)

opt = getopt(spec)

print_usage <- function(spec=NULL){
        cat(getopt(spec, usage=TRUE));
        cat("Usage example: \n");
        cat("
Usage example: 
1. Rscript edger_analysis.r --count All_gene_counts.list --config L01_L02_L06_vs_L03_L04_L05.de.config --out.prefix L01_L02_L06_vs_L03_L04_L05 --filter cmp --fpkm All_gene_expression.list --type SRPBM
2. Rscript edger_analysis.r --count All_gene_counts.list --config L01_L02_L06_vs_L03_L04_L05.de.config --out.prefix L01_L02_L06_vs_L03_L04_L05 --filter cmp 
3. Rscript edger_analysis.r --count All_gene_counts.list --config L01_vs_L03.de.config --out.prefix L01_vs_L03 --filter cmp --fpkm All_gene_expression.list --type SRPBM --disper 0.16
4. Rscript edger_analysis.r --count All_gene_counts.list --config L01_vs_L03.de.config --out.prefix L01_vs_L03 --filter cmp --disper 0.01

Options: 
--help		-h	NULL	get this help
--count		character       the read count file for RNA-SEQ [forced]
--config	character	the multiple factor matrix file [forced]
--out.prefix	character	the prefix of output file [forced]
--filter	character	the method for filtering low gene expression,cmp/count/no can be choosed [forced, default: cpm]
--fpkm		character	the expression file [optional]
--type		character	the type of expression,FPKM/TPM/SRPBM if define --fpkm [optional, default: FPKM]
--disper	character	dispersion：human(0.16)/mouse(0.01), because it's genetically identical model organisms,rat and other species default is 0.01(no duplication is must) [optional, default: 0.01]

\n")
	q(status=1);
}

if ( !is.null(opt$help) ) { print_usage(spec) }
if ( is.null(opt$count) )	{ print_usage(spec) }
if ( is.null(opt$config) )	{ print_usage(spec) }
if ( is.null(opt$out.prefix) )	{ print_usage(spec) }
if ( is.null(opt$filter) )	{ opt$filter = "cpm"}

require(edgeR)
f_mat <- read.table(opt$config, header=TRUE,check.names =F)  ##read de.config
ori_data <- read.delim(opt$count, row.names = 1, header=TRUE,check.names =F)  ##read count file
colnames(ori_data)<-read.delim(opt$count, row.names = 1, header=F,check.names =F,stringsAsFactors=F,nrows=1)
count_data <- ori_data[ , as.vector(f_mat[,1]) ]

condition <- f_mat[,2]
con <- as.factor(condition)

filter=opt$filter   ###cpm /count
if(filter=="cpm"){ # filter low expression gene
	keep <- rowSums(cpm(count_data)>1) >= min( table(con) )
}else if(filter=="count"){
	keep <- rowSums(count_data)>=2
}else if(filter=="no"){
	keep <- rowSums(count_data)>0
}

filter_count<-count_data[keep,]

##############begin deg analysis
y <- DGEList(counts=filter_count, group=condition)
y <- calcNormFactors(y)
if(dim(filter_count)[2]>2){
	y <- estimateCommonDisp(y)
	y <- estimateTagwiseDisp(y)
	et <- exactTest(y)
}
if(dim(filter_count)[2]==2){
	if( is.null(opt$disper) ){
		opt$disper = "0.01"
		dis =as.numeric(opt$disper)
	}else{
		dis = as.numeric(opt$disper)
	}
	print(dis)
	et <- exactTest(y,dispersion=dis)
}
res <- as.data.frame(topTags(et,dim(et)[1]))

res<-as.matrix(res)
res<-res[rownames=rownames(filter_count),]
res<-as.data.frame(res)

if(!is.null(opt$fpkm)){
	#####add fpkm into count file
	fpkm <- read.delim(opt$fpkm, row.names = 1, header=TRUE,check.names =F)  ##read count file
	colnames(fpkm)<-read.delim(opt$fpkm, row.names = 1, header=F,check.names =F,stringsAsFactors=F,nrows=1)
	fpkm_data <- fpkm[ , as.vector(f_mat[,1]) ]
	fpkm_data <- as.matrix(fpkm_data)
	fpkm_data <- fpkm_data[rownames=rownames(filter_count),]

	type="FPKM"
	if(!is.null(opt$type)){type=opt$type}
	deg_all	<-data.frame(filter_count,fpkm_data,res$PValue,res$FDR,res$logFC)
	header<-c("#ID",paste(colnames(y),"_Count",sep=""),paste(colnames(y),"_",type,sep=""),"PValue","FDR","log2FC")
}else{
	deg_all<-data.frame(filter_count,res$PValue,res$FDR,res$logFC)
	header<-c("#ID",paste(colnames(y),"_Count",sep=""),"PValue","FDR","log2FC")
}

write.table(t(header),paste(opt$out.prefix,".all",sep=""),sep="\t",quote=FALSE,row.names=F,col.names=F)
write.table(deg_all,paste(opt$out.prefix,".all",sep=""),sep="\t",quote=FALSE,row.names=T,col.names=F,append=T)

