library(getopt)
#library(Vennerable)
#author:huangls
get_data_list<-function(str,names,out_file){
	input_names<-unlist(strsplit(names,','))
	data_paths<-unlist(strsplit(str,','))
	
	if (length(data_paths) != length(input_names)){
		cat("data name not eq data Num\n\n")
		h(0)
	}
	data_list<-list()
	uniq_all_id<-c()
	i<-1
	for(infile in data_paths){
		#cat(infile,"\n")
		mydata<-read.table(infile,header=F,sep="\t")
		
		uniq_all_id<-union(as.character(mydata[,1]),uniq_all_id)
		data_list[[i]]<-mydata[,1]
		i<-i+1
	}
        #aa=matrix(data = "NA", nrow = length(uniq_all_id), ncol = length(data_list)+1, byrow = FALSE,dimnames = NULL)
        #aa<-as.data.frame(aa)

        #aa[,1]=as.character(uniq_all_id)
        #print(head(aa))
        names(data_list)<-input_names
        aa<-data.frame(all_ID=uniq_all_id)
        for(i in 1:length(data_list)){

                tmp=uniq_all_id %in% as.character(data_list[[i]])
                aa<-cbind(aa,tmp)




        }
        colnames(aa)=c("all_ID",input_names)
        write.table(aa,file=out_file,quote=F,sep='\t',row.names=F,col.names=T)
        return(data_list)

}
h<-function(x){
	cat(getopt(spec, usage=TRUE));
	cat("\n\n")
	cat("
		Usage: Rscript  Vennerable_v1.0.R -d data1,data2    -n name1,name2 [ -o Prefix_of_output_files (output)] [-h]
		
		Note:
			multi input files must split by ',', and  labels order of '-n' option must correspond to input files;
			each input files should include one column data and there isn't a header;  


");
	
	
	quit(save="no",status=x)
}
spec<-matrix(c(
				'help','h', 0, "logical","for help",
				'mydata', 'd', 1, "character","datafile,must split by ',' ,max files was 9",
				'dataname','n',1, "character","sample titile in Venn Pic",
				'outdir','o',1,"character","outdir,optional",
				'doWeight','W' , 0, "logical","weights to associate with each possible combination of Set intersections",
				'prefix', 'p', 1, "character","output file prefix"
		), ncol=5, byrow=TRUE)
opt <- getopt(spec);

if (! is.null(opt$help)) { h(0) }
if (is.null(opt$dataname)) { h(0) }
if (is.null(opt$mydata)) { h(0) }
if (is.null(opt$doWeight)) { opt$doWeight=FALSE }
file.prefix <- ifelse(is.null(opt$prefix), "output", opt$prefix)
#cat(data_paths,"\n")
library(Vennerable)


##########################

if ( is.null(opt$outdir))	{
	pic_name_pdf<-paste(getwd(),"/",file.prefix,"_","venn",".pdf",sep="")
	pic_name_png<-paste(getwd(),"/",file.prefix,"_","venn",".png",sep="")
	out_file<-paste(getwd(),"/",file.prefix,"_","venn",".id.list",sep="")
}else{
	if( !file.exists(opt$outdir) ){
		if( !dir.create(opt$outdir, showWarnings = FALSE, recursive = TRUE) ){
			stop(paste("dir.create failed: outdir=",opt$outdir,sep=""))
		}
	}

	if(substr(opt$outdir,nchar(opt$outdir),nchar(opt$outdir))=="/"){
		opt$outdir<-substr(opt$outdir,1,nchar(opt$outdir)-1)
	}
	pic_name_pdf<-paste(opt$outdir,"/",file.prefix,"_","venn",".pdf",sep="")
	pic_name_png<-paste(opt$outdir,"/",file.prefix,"_","venn",".png",sep="")
	out_file<-paste(opt$outdir,"/",file.prefix,"_","venn",".list",sep="")
}




#data_list_up<-list()
#if (!is.null(opt$mydata_aa)){
#	data_list_up<-get_data_list(opt$mydata_aa,opt$dataname)
#}

data_list<-list()
data_list<-get_data_list(opt$mydata,opt$dataname,out_file)

data<-Venn(data_list)
infile_num=length(data_list)
isWeight=opt$doWeight
if(infile_num>=6){
	isWeight=FALSE
}
pdf(file=pic_name_pdf, height=10, width=10)
plot(data,doWeight=isWeight)
dev.off()
png(file=pic_name_png, height=800, width=800)
plot(data,doWeight=isWeight)
dev.off()
