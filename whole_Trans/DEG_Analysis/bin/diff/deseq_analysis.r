args<-commandArgs(TRUE)
usage <- function() {
        print("-------------------------------------------------------------------------------")
        print("Usage: Rscript deseq_analysis.r read_count.txt multi_factor_matrix.txt 0.01 0 out.de")
        print("1) read_count.txt: the read count file for RNA-SEQ")
        print("2) multi_factor_matrix.txt: the multiple factor matrix")
        print("3) out.de: the output of DE")
	print("4) filter: the method used for filter low gene expression(cpm/count/no) no means only filter the genes of which all count equal 0")
        print("5) expression.list: not must ")
        print("6) expression.type:(FPKM or TPM) :not must if givn expression.list default FPKM")
        print("-------------------------------------------------------------------------------")
}
if(length(args)!=4 && length(args)!=5 && length(args)!=6){
	usage()
	stop("Parameter length Wrong!")
}
filter=args[4]

require(edgeR)
require(DESeq)

f_mat <- read.table(args[2], header=TRUE,check.names =F)  ##read de.config
ori_data <- read.delim(args[1], row.names = 1, header=TRUE,check.names =F)  ##read count file
colnames(ori_data)<-read.delim(args[1], row.names = 1, header=F,check.names =F,stringsAsFactors=F,nrows=1)
count_data <- ori_data[ , as.vector(f_mat[,1]) ]

condition <- f_mat[,2]
con <- as.factor(condition)
if(filter=="cpm"){ # filter low expression gene
	keep <- rowSums(cpm(count_data)>1) >= min( table(con) )
}else if(filter=="count"){
	keep <- rowSums(count_data)>=2
}else if(filter=="no"){
	keep <- rowSums(count_data)>0
}

filter_count<-as.matrix(count_data[keep,])

cds <- newCountDataSet(countData=filter_count, conditions = condition)
cds <- estimateSizeFactors(cds)
sizeFactors(cds)
cds <- estimateDispersions(cds, fitType="local")        # can for miRNA or small gene number

## Differential Expression
# Identify those features that are direntially expressed in the two groups.
print("Differential Expression ... ")
c1 <- names( table(condition) )[1]
c2 <- names( table(condition) )[2]
res <- nbinomTest(cds, c1, c2)
print("Differential Expression is over ")

if(length(args)==5 || length(args)==6 ){
        #####add fpkm into count file
        fpkm <- read.delim(args[5], row.names = 1, header=TRUE,check.names =F)  ##read count file
        colnames(fpkm)<-read.delim(args[5], row.names = 1, header=F,check.names =F,stringsAsFactors=F,nrows=1)
        fpkm_data <- fpkm[ , as.vector(f_mat[,1]) ]
        fpkm_data <- as.matrix(fpkm_data)
        fpkm_data <- fpkm_data[rownames=rownames(filter_count),]

        type="FPKM"
        if(length(args)==6){type=args[6]}
        deg_all <-data.frame(filter_count,fpkm_data,res$pval,res$padj,res$log2FoldChange)
        header<-c("#ID",paste(colnames(filter_count),"_Count",sep=""),paste(colnames(filter_count),"_",type,sep=""),"PValue","FDR","log2FC")
}else{

	deg_all<-data.frame(filter_count,res$pval,res$padj,res$log2FoldChange)
	header<-c("#ID",paste(colnames(filter_count),"_Count",sep=""),"PValue","FDR","log2FC")
}

write.table(t(header),paste(args[3],".all",sep=""),sep="\t",quote=FALSE,row.names=F,col.names=F)
write.table(deg_all,paste(args[3],".all",sep=""),sep="\t",quote=FALSE,row.names=T,col.names=F,append=T)






