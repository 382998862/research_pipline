args<-commandArgs(TRUE)
usage <- function() {
  print("-------------------------------------------------------------------------------")
  print("Usage: Rscript TFBStools.r spe_id input strand outpath")
  print("1) spe_id: the taxid of every species, which can be found in Taxonomy of NCBI")
  print("2) input: one gene's upstream sequence (1 kb), xxx.fa")
  print("3) strand: strand of gene")
  print("4) min.score: the min score of pridict TFBS")
  print("5) outpath: output path")
  print("-------------------------------------------------------------------------------")
  q(status=1);
}

if(length(args)!=5){
  usage()
}

#######require R-package 
library(JASPAR2018)
library(TFBSTools)

#####Load needed species dataset from JASPAR 
opts <- list()
opts[["species"]] <- args[1] ####(the spe_id in detail.cfg file;  Human:9606;Mouse:10090;Rat:10116)
opts[["all_versions"]] <- TRUE
PFMatrixList <- getMatrixSet(JASPAR2018, opts)
PWMatrixList <- toPWM(PFMatrixList, pseudocounts=0.8)

########import sequende of one gene
library(Biostrings)
gene_fa <- readDNAStringSet(args[2])
gene_id <- names(gene_fa)
subject <- DNAString(gene_fa[[1]])
sitesetList <- searchSeq(PWMatrixList, subject, seqname=gene_id, min.score=args[4], strand=args[3])
output <- as.matrix(writeGFF3(sitesetList,scoreType="relative"))

if(length(output) != 0){

pval <- pvalues(sitesetList, type="TFMPvalue")
pvalMatrix <- as.matrix(unlist(pval))
colnames(pvalMatrix)="Pvalue"

Model_id <- as.matrix(rownames(output))
colnames(Model_id)="Model_id"
out <- gsub(" ","",output[,1:8])
Attributes <- as.matrix(output[,9,drop=F])
colnames(Attributes)="Attributes"
rawRes <- cbind(Model_id, out, Attributes)

##############merge results 
if (nrow(rawRes) == nrow(pvalMatrix)) {
	Res <- cbind(rawRes[,-c(3:4)], pvalMatrix)
	Res <- as.matrix(Res)
}else{
	num <- matrix(NA,nrow=nrow(rawRes),ncol=1)
	for (i in 1:nrow(rawRes)) {
		model <- strsplit(as.character(rawRes[i,1]),".",fixed=TRUE)
		modelID <- paste(as.character(model[[1]][1]),as.character(model[[1]][2]),sep=".")
		#model <- substr(as.character(rawRes[i,1]),1,8)
		tmp <- strsplit(rawRes[i,10],"=")
		tf <- strsplit(as.character(tmp[[1]][2]),";")[[1]][1]
		seq <- as.character(tmp[[1]][4])
		id <- paste(modelID,rawRes[i,5],rawRes[i,6],rawRes[i,7],tf,seq,sep="_")
		num[i,1] <- id
	}	
	uni_num <- unique(num)
	newRes <- cbind(rawRes,num)
	
	newOut <- NULL
	for (i in 1:nrow(uni_num)) {
		id <- uni_num[i,1]
		tmp <- newRes[which(newRes[,11] == id),,drop=F]
		if(length(tmp[,1]) !=1) {
			newOut <- rbind(newOut,tmp[1,])
		}else{
			newOut <- rbind(newOut,tmp)
		}
	}
	Res <- cbind(newOut[,-c(3:4,11)], pvalMatrix)
	Res <- as.matrix(Res)
}

write.table(Res,file = paste(args[5],"/",gene_id,"_TFBS_predictRes.txt",sep=""), col.names = T,row.names = F, quote = F,sep = "\t")
}else{
write.table(output,file = paste(args[5],"/",gene_id,"_TFBS_predictRes.txt",sep=""), col.names = T,row.names = F, quote = F,sep = "\t")
}
