library(getopt)
spec=matrix(c(
	'help','h',0,"logical",
	'inpath','i',1,"character",
	'outpath','o',1,"character"
	),byrow=TRUE,ncol=4
)
opt=getopt(spec)
print_usage<-function(spec=NULL){
	cat(getopt(spec,usage=TRUE));
	cat("Usage example: \n")
	cat("
		Description: extract info from multiple KEGG pathways
		Contact: Yanhua Wen <wenyh@biomarker>	
		Usage example:
			
		Options:
		--help		NULL		get this help
		--inpath	character	inpath contaned all kegg xml/kgml
		--outpath	character	out put path
		\n")

	q(status=1);
}
if(is.null(opt$inpath))		print_usage(spec)
if(is.null(opt$outpath))	print_usage(spec)


library(XML)

c("gene","compound","ortholog","group","map","enzyme","reaction")->type
list.files(opt$inpath,full.names=T)->file
out_relation<-paste(opt$outpath,"relation.txt",sep="/")
out_entry<-paste(opt$outpath,"entry.txt",sep="/")
out_group<-paste(opt$outpath,"group.txt",sep="/")

data.frame("id","name","type","path_name")->en
write.table(en,out_entry,sep="\t",row.names=F,col.names=F,eol="\n",quote=F) 
  
data.frame("path_name","type","entry1","entry2","subtype")->re
write.table(re,out_relation,sep="\t",row.names=F,col.names=F,eol="\n",quote=F)

data.frame("path_name","group_id","component")->gr
write.table(gr,out_group,sep="\t",row.names=F,col.names=F,eol="\n",quote=F)
 
for(i in 1:length(file)){
  doc<-xmlParse(file[i])
  path_title<-as.character(getNodeSet(doc,'//pathway/@title'))
  path_name<-as.character(getNodeSet(doc,'//pathway/@name'))
#  chartr("path:", "", path_name)->path_name
  path_name=strsplit(path_name,split=":")[[1]][2]
  as.character(getNodeSet(doc,'//pathway/entry/@id'))->eid
  as.character(getNodeSet(doc,'//pathway/entry/@name'))->ename 
  as.character(getNodeSet(doc,'//pathway/entry/@type'))->etype
  data.frame(eid,ename,etype,path_name)->oo
  write.table(oo,out_entry,sep="\t",append=T,row.names=F,col.names=F,eol="\n",quote=F)##提取所有entry信息
  ##获取relation信息
  as.character(getNodeSet(doc,'//pathway/relation/@type'))	->rtype
  as.character(getNodeSet(doc,'//pathway/relation/@entry1'))	->rentry1
  as.character(getNodeSet(doc,'//pathway/relation/@entry2'))	->rentry2
  getNodeSet(doc,'//pathway/relation')->relation_node
  if(length(relation_node)!=0){
	  for(k in 1:length(relation_node)){	
	    xmlSApply(relation_node[[k]],xmlAttrs)->subtype
	    if(subtype=="NULL"){subtype="--"}else{subtype<-paste(subtype[1,],collapse=",")}
	    ss<-data.frame(path_name,rtype[k],rentry1[k],rentry2[k],subtype)
	    write.table(ss,out_relation,sep="\t",append=T,row.names=F,col.names=F,eol="\n",quote=F)
	 }
   }
   ##获取group id信息
   getNodeSet(doc,"//entry[@type='group']")->node
   as.character(getNodeSet(doc,"//entry[@type='group']/@id"))->node_id
   if(length(node)!=0){
	for(k in 1:length(node)){
		xmlSApply(node[[k]],xmlAttrs)->component
		id_num<-length(component)
		component<-as.character(component[2:id_num])
		ss<-data.frame(path_name,node_id[k],paste(component,collapse=","))
		write.table(ss,out_group,sep="\t",append=T,row.names=F,col.names=F,eol="\n",quote=F)
	}    
   }      
}



######################################把group拆分
read.csv(out_group,header=T,sep="\t")->group
if(dim(group)[1]==0){
	q(status=1);
}
read.csv(out_relation,header=T,sep="\t")->relation
retmp1<-paste(out_relation,".tmp1",sep="")
retmp2<-paste(out_relation,".tmp2",sep="")
write.table(re,retmp1,sep="\t",row.names=F,col.names=F,eol="\n",quote=F)
write.table(re,retmp2,sep="\t",row.names=F,col.names=F,eol="\n",quote=F)

for(i in 1:dim(relation)[1]){
   id<-as.character(relation[i,3])   #########看relation的第一个元素是否有group 
   path=as.character(relation[i,1])
   which(group[,1]==path & group[,2]==id)->s
   if(length(s)==0)
       write.table(relation[i,],retmp1,sep="\t",row.names=F,col.names=F,eol="\n",append=T,quote=F)
   else{
       as.character(group[s,3])->component
       strsplit(component,split=",")[[1]]->component
       data.frame(path,relation[i,2],component,relation[i,4],relation[i,5])->out
       write.table(out,retmp1,sep="\t",row.names=F,col.names=F,eol="\n",append=T,quote=F)
   }
}
 
read.csv(retmp1,sep="\t",header=T)->relation
for(i in 1:dim(relation)[1]){
   id<-as.character(relation[i,4])     #########看relation的第二个元素是否有group 
   path=as.character(relation[i,1])
   which(group[,1]==path&group[,2]==id)->s
   if(length(s)==0)
	write.table(relation[i,],retmp2,sep="\t",row.names=F,col.names=F,eol="\n",append=T,quote=F)
   else{
       as.character(group[s,3])->component
       strsplit(component,split=",")[[1]]->component
       data.frame(path,relation[i,2],relation[i,3],component,relation[i,5])->out
       write.table(out,retmp2,sep="\t",row.names=F,col.names=F,eol="\n",append=T,quote=F)
   }
}                                                              ###############以上两步是把group替换成单一的元素 

file.remove(out_relation)
file.remove(retmp1)
file.rename(retmp2, out_relation)


