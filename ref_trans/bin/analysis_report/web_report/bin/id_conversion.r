#!/share/nas2/genome/biosoft/R/2.15.1/lib64/R/bin/Rscript

# load library
library('getopt');


#listMarts()                  #The function listMarts will display all available BioMart web services
#ensembl=useMart("ensembl")   #list the datasets
#useMart("ENSEMBL_MART_ENSEMBL",dataset="sscrofa_gene_ensembl", host="www.ensembl.org") ##2015-11-10
#listDatasets(ensembl)

#-----------------------------------------------------------------
# getting parameters
#-----------------------------------------------------------------
#get options, using the spec as defined by the enclosed list.
#we read the options from the default: commandArgs(TRUE).
spec = matrix(c(
  'help' , 'h', 0, "logical",
  'infile' , 'i', 1, "character",
  'outfile' , 'o', 1, "character",
  'id_name' , 'n',1,"character",
  'id.col' , 'x', 1, "integer",
  'skip' , 'k', 1, "integer",
  'database','d',1,"character",
  'dataset','s',1,"character",
  'host', 't',1,"character"
), byrow=TRUE, ncol=4);
opt = getopt(spec);

#-----------------------------------------------------------------
# define usage function
#-----------------------------------------------------------------
print_usage <- function(spec=NULL){
  cat(getopt(spec, usage=TRUE));
  cat("Usage example: \n")
  cat("
      Usage example:
      Rscript id_conversion.r --infile in.xls --outfile out_filename  --id_name ensembl_gene_id --id.col 1 \\
      --database ENSEMBL_MART_ENSEMBL --dataset sscrofa_gene_ensembl --host www.ensembl.org \\1
      Options:
      --help		NULL 		    get this help
      --infile 	character 	the input file [forced]
      --outfile	character 	the filename for output file[forced]
      --id_name	character	the class of id [forced]
      --id.col 	integer 	the col for id [forced]
      --host	character	get the name of host you use[forced ,we list all host in cfg file]
      --database  character	get the name of database you use[forced ,if default ,we list all databases]
      --dataset	character	get the name of dataset you use[forced ,if default ,we list all datasets
                                of the database you choose before]
      --skip 	integer 	the number of line for skipping [optional, default: 0]
      \n")
  q(status=1);
}


# if help was asked for print a friendly message
if ( !is.null(opt$help) ) { print_usage(spec) }

# check non-null args
if ( is.null(opt$infile))	  { print_usage(spec) }
if ( is.null(opt$outfile))	{ print_usage(spec) }
if ( is.null(opt$id_name))  { print_usage(spec) }
if ( is.null(opt$id.col) )	{ print_usage(spec) }
if ( is.null(opt$database))	{
  cat("Please choose one of the following databases:\n")
  listMarts()
  print_usage(spec)
  }
if (is.null(opt$dataset) && !is.null(opt$database)){
  cat("Please choose one of the following datasets:\n")
  ensembl=useMart(opt$database)
  listDatasets(ensembl)
  print_usage(spec)
}



library("biomaRt")

#set some reasonable defaults for the options that are needed,
#but were not specified.
if ( is.null(opt$skip ) )		{ opt$skip = 0 }

#-----------------------------------------------------------------
# reading data
#-----------------------------------------------------------------
data <- read.table(opt$infile, skip=opt$skip,header = T,sep = "\t",fill=TRUE,comment.char="",check.names =F)
# check dim
data.dim <- dim(data)
if ( is.null(data.dim) ){
  cat("Final Error: the format of infile is error, dim(data) is NULL\n")
  print_usage(spec)
}
# check col size
if ( data.dim[2] < opt$id.col ){
  cat("Final Error: id.col > the col of infile\n")
  print_usage(spec)
}


#set the database
#ensembl = useMart(opt$database,dataset=opt$dataset)   #调用数据库
ensembl = useMart(opt$database,dataset=opt$dataset,host=opt$host)   #调用数据库2015-11-10
#useMart("ENSEMBL_MART_ENSEMBL",dataset="sscrofa_gene_ensembl", host="www.ensembl.org")
# create ids
ids<-as.character(data[,opt$id.col])

#########连接不上服务器的数据库###############
#ensembl = useMart("ensembl", dataset="ensembl", host="192.168.1.207:22",
#          path="/share/nas1/shumy/Rattus_norvegicus", port=22,archive=FALSE)
#useMart("ENSEMBL_MART_ENSEMBL",dataset="sscrofa_gene_ensembl", host="www.ensembl.org")
#############################################



#-----------------------------------------------------------------
# id conversion
#-----------------------------------------------------------------
#filters = listFilters(ensembl)         #[input]The listFilters function shows you all available lters in the selected dataset.
#attributes = listAttributes(ensembl)   #[output]The listAttributes function displays all available attributes in the selected dataset.

#newdata<-getBM(attributes=c(opt$id_name, 'external_gene_name','description'),
#      filters = opt$id_name, values = ids, mart = ensembl)

if(opt$host=="www.ensembl.org"){
  newdata<-getBM(attributes=c(opt$id_name, 'external_gene_name','description'),
                 filters = opt$id_name, values = ids, mart = ensembl)
}else{
  newdata<-getBM(attributes=c(opt$id_name, 'external_gene_id','description'),
                 filters = opt$id_name, values = ids, mart = ensembl)
}
row.num <- length(newdata[,1])
col.num <- length(newdata[1,])
for (i in 1:row.num){
  for (j in 2:col.num){
    if(newdata[i,j]==""){
      newdata[i,j]="--"
    }
  }
}


#set new col names
colnames(newdata)<-c(colnames(data)[opt$id.col],"gene symbol","gene name")

#-----------------------------------------------------------------
# output
#-----------------------------------------------------------------
write.table(newdata,file = paste(opt$outfile,".xls",sep=""),row.names = F,col.names = T,quote = F,sep="\t",na = "NA")
