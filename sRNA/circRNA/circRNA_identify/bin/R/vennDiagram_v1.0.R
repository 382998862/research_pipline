library(VennDiagram)
library(getopt)
library(scales)
library(RColorBrewer)
draw_venn_5<-function(data_list_down,data_list_up,file_name,main){
	if (length(data_list_up)==0){
		png(file_name, height = 3000, width = 3000, res = 500, units = "px")
	
		draw.venn<-venn.diagram(
			x =data_list_down,
			main=main,
			main.cex=3,
			filename = NULL,
			col = c("dodgerblue", "goldenrod1", "darkorange1", "seagreen3", "orchid3"),
			fill = c("dodgerblue", "goldenrod1", "darkorange1", "seagreen3", "orchid3"),
			alpha = 0.50,
			cex = c(1.5, 1.5, 1.5, 1.5, 1.5, 1, 0.8, 1, 0.8, 1, 0.8, 1, 0.8,
				1, 0.8, 1, 0.55, 1, 0.55, 1, 0.55, 1, 0.55, 1, 0.55, 1, 1, 1, 1, 1, 1.5),
			cat.col = c("dodgerblue", "goldenrod1", "darkorange1", "seagreen3", "orchid3"),
			cat.cex = 1.5,
			cat.fontface = "bold",
			cat.fontfamily = "serif",
			#cat.default.pos="outer",
			#cat.dist = c(-0.08, -0.09, -0.08,-0.08,-0.08),
			margin = 0.05
		);
	grid.draw(draw.venn);
	dev.off()
	
	} else {
		png(file_name, height = 3000, width = 3000, res = 500, units = "px")
		draw.venn.down<-venn.diagram(
				x =data_list_up,
				#main=main,
				main.cex=4,
				main.pos = c(0.5, 1.05),
				filename = NULL,
				label.col="red",
				col = c("dodgerblue", "goldenrod1", "darkorange1", "seagreen3", "orchid3"),
				fill = c("dodgerblue", "goldenrod1", "darkorange1", "seagreen3", "orchid3"),
				alpha = 0.50,
				cex = c(1.5, 1.5, 1.5, 1.5, 1.5, 1, 0.8, 1, 0.8, 1, 0.8, 1, 0.8,
						1, 0.8, 1, 0.55, 1, 0.55, 1, 0.55, 1, 0.55, 1, 0.55, 1, 1, 1, 1, 1, 1.5),
				cat.col = c("dodgerblue", "goldenrod1", "darkorange1", "seagreen3", "orchid3"),
				cat.cex = 1.5,
				cat.fontface = "bold",
				cat.fontfamily = "serif",
				cat.default.pos="outer",
				#cat.dist = c(0.2, 0.2, 0.2,0.2,0.2),
				margin = 0.05
		);
		draw.venn.up<-venn.diagram(
				x =data_list_down,
				lty=0,
				main="",
				label.col="blue",
				category.names="",
				main.cex=1.5,
				main.pos=c(0.5,1),
				filename = NULL,
				col = "white",
				fill = "white",
				alpha = 0,
				cex = c(1.5, 1.5, 1.5, 1.5, 1.5, 1, 0.8, 1, 0.8, 1, 0.8, 1, 0.8,
						1, 0.8, 1, 0.55, 1, 0.55, 1, 0.55, 1, 0.55, 1, 0.55, 1, 1, 1, 1, 1, 1.5),
				cat.col = "white",
				cat.cex = 1,
				cat.fontface = "bold",
				cat.fontfamily = "serif",
				#cat.default.pos="outer",
				#cat.dist = c(1.08, 1.09, 1.08,1.08,1.08),
				#cat.dist = c(0.31, 0.31, 0.31,0.31,0.31),
				margin = 0.05
		);
		
		grid.draw(draw.venn.down);
		grid.draw(draw.venn.up);
		dev.off()
	}
	
}


draw_venn_4<-function(data_list_down,data_list_up,file_name,main){

	png(file_name, height = 3000, width = 3000, res = 500, units = "px")

	draw.venn<-venn.diagram(
		x =data_list_down,
		main=main,
		main.cex=3,
		
		col = c("dodgerblue", "goldenrod1", "darkorange1", "seagreen3"),
		
		lty = 1,
		lwd = 2,
		fill = c("dodgerblue", "goldenrod1", "darkorange1", "seagreen3"),
		alpha = 0.40,
		#label.col = c("orange", "white", "darkorchid4", "white", "white", "white", "white", "white", "darkblue", "white", "white", "white", "white", "darkgreen", "white"),
		#label.col=c("dodgerblue", "goldenrod1", "darkorange1", "seagreen3"),	
		cex = 2.3,
		fontfamily = "serif",
		#fontface = "bold",
		cat.col = c("dodgerblue", "goldenrod1", "darkorange1", "seagreen3"),
		cat.cex = 2,
		cat.fontfamily = "serif",
		cat.fontface = "bold",
		filename=NULL
		);

	grid.draw(draw.venn);
	dev.off()
	if (length(data_list_up)!=0){
		png(file_name, height = 3000, width = 3000, res = 500, units = "px")
		draw.venn.down<-venn.diagram(
				x =data_list_down,
				main=main,
				main.cex=4,
				main.pos = c(0.5, 0.98),
				col = c("dodgerblue", "goldenrod1", "darkorange1", "seagreen3"),
				label.col="blue",
				lty = 1,
				lwd = 2,
				fill = c("dodgerblue", "goldenrod1", "darkorange1", "seagreen3"),
				alpha = 0.40,
				#label.col = c("orange", "white", "darkorchid4", "white", "white", "white", "white", "white", "darkblue", "white", "white", "white", "white", "darkgreen", "white"),
				#label.col=c("dodgerblue", "goldenrod1", "darkorange1", "seagreen3"),	
				cex = 1.7,
				fontfamily = "serif",
				#fontface = "bold",
				cat.col = c("dodgerblue", "goldenrod1", "darkorange1", "seagreen3"),
				cat.cex = 2,
				cat.fontfamily = "serif",
				cat.fontface = "bold",
				filename=NULL
		);
		draw.venn.up<-venn.diagram(
				x =data_list_up,
				
				main.cex=1.5,
				category.names="",
				col = "white",
				label.col="red",
				lty = 0,
				lwd = 2,
				fill = "white",
				alpha = 0,
				#label.col = c("orange", "white", "darkorchid4", "white", "white", "white", "white", "white", "darkblue", "white", "white", "white", "white", "darkgreen", "white"),
				#label.col=c("dodgerblue", "goldenrod1", "darkorange1", "seagreen3"),	
				cex = 1.7,
				fontfamily = "serif",
				#fontface = "bold",
				cat.col = "white",
				cat.cex = 2,
				cat.fontfamily = "serif",
				cat.fontface = "bold",
				filename=NULL
		);
		grid.draw(draw.venn.down);
		grid.draw(draw.venn.up);

		dev.off()
	}
	
}

draw_venn_3<-function(data_list_down,data_list_up,file_name,main){
	png(file_name, height = 3000, width = 3000, res = 500, units = "px")
	draw.venn<-venn.diagram(
		x =data_list_down,
		filename = NULL,
		main=main,
		main.cex=3,
		lty=1,
		lwd=2,
		col = c("dodgerblue", "goldenrod1", "orchid3"),
		fill = c("dodgerblue", "goldenrod1", "orchid3"),
		alpha = 0.5,
		#label.col = c("dodgerblue", "goldenrod1", "orchid3"),
		cex = 2.5,
		fontfamily = "serif",
		#fontface = "bold",
		#cat.default.pos = "text",
		cat.col = c("dodgerblue", "goldenrod1", "orchid3"),
		cat.cex = 2,
		cat.fontfamily = "serif",
		cat.fontface = "bold",
		cat.dist = c(0.04, 0.04, 0.04),
		#cat.pos = c(0,0,6)
	);
	grid.draw(draw.venn);
	dev.off()
	if (length(data_list_up)!=0){
		png(file_name, height = 3000, width = 3000, res = 500, units = "px")
		draw.venn.down<-venn.diagram(
				x =data_list_down,
				filename = NULL,
				#main.fontfamily="",
				main.cex=5,
				main.pos = c(0.5, 1.07),
				main=main,
				lty=1,
				lwd=2,
				col = c("dodgerblue", "goldenrod1", "orchid3"),
				fill = c("dodgerblue", "goldenrod1", "orchid3"),
				alpha = 0.5,
				label.col="blue",
				#label.col = c("dodgerblue", "goldenrod1", "orchid3"),
				cex = 2,
				fontfamily = "serif",
				#fontface = "bold",
				#cat.default.pos = "text",
				cat.col = c("dodgerblue", "goldenrod1", "orchid3"),
				cat.cex = 1.5,
				cat.fontfamily = "serif",
				cat.fontface = "bold",
				cat.dist = c(0.05, 0.05, 0.05),
		#cat.pos = c(0,0,6)
		);
		draw.venn.up<-venn.diagram(
				x =data_list_up,
				filename = NULL,
				
				main.cex=2,
				lty=0,
				lwd=2,
				category.names="",
				col = "white",
				fill = "white",
				alpha = 0,
				#label.col = c("dodgerblue", "goldenrod1", "orchid3"),
				label.col="red",
				cex = 2,
				fontfamily = "serif",
				#fontface = "bold",
				#cat.default.pos = "text",
				cat.col ="white",
				cat.cex = 0,
				cat.fontface = "bold",
				cat.fontfamily = "serif",
				cat.dist = c(0.05, 0.05, 0.05),
		#cat.pos = c(0,0,6)
		);
		grid.draw(draw.venn.down);
		grid.draw(draw.venn.up);
		
		dev.off()
		
	}
				
	
}
draw_venn_2<-function(data_list_down,data_list_up,file_name,main){
	png(file_name, height = 800, width = 800, res = 200, units = "px")
	draw.venn.down<-venn.diagram(
		x =data_list_down,
		lty=1,
		
		main=main,
		main.cex=2.5,
		filename = NULL,
		lwd = 2,
		label.col = "black",
		col=c("dodgerblue", "goldenrod1"),
		fill = c("dodgerblue", "goldenrod1"),
		alpha = 0.5,
		#label.col = c("dodgerblue", "goldenrod1"),
		cex = 2.5,
		fontfamily = "serif",
		#fontface = "bold",
		cat.col = c("dodgerblue", "goldenrod1"),
		cat.cex = 2,
		cat.fontfamily = "serif",
		cat.fontface = "bold",
		cat.default.pos="outer",
		cat.dist = c(0.03, 0.03),
		cat.pos = 6,
		scaled=FALSE,
	);
	
	
	
	grid.draw(draw.venn.down);
	dev.off()
	
	if (length(data_list_up)!=0){
		png(file_name, height = 3000, width = 3000, res = 500, units = "px")
		draw.venn.down<-venn.diagram(
				x =data_list_down,
				lty=1,
				
				main=main,
				main.cex=3,
				filename = NULL,
				lwd = 2,
				label.col = "blue",
				col=c("dodgerblue", "goldenrod1"),
				fill = c("dodgerblue", "goldenrod1"),
				alpha = 0.5,
				#label.col = c("dodgerblue", "goldenrod1"),
				cex = 2.5,
				fontfamily = "serif",
				#fontface = "bold",
				cat.col = c("dodgerblue", "goldenrod1"),
				cat.cex = 2,
				cat.fontfamily = "serif",
				cat.fontface = "bold",
				cat.default.pos="outer",
				cat.dist = c(0.03, 0.03),
				cat.pos = 6
		);
		draw.venn.up<-venn.diagram(
				x =data_list_up,
				lty=0,
				
				category.names="",
				main.cex=2,
				filename = NULL,
				lwd = 2,
				col="white",
				fill = "white",
				alpha = 0,
				label.col = "red",
				cex = 2.5,
				fontfamily = "serif",
				#fontface = "bold",
				cat.col = "white",
				cat.cex = 2,
				cat.fontfamily = "serif",
				cat.fontface = "bold",
				cat.default.pos="outer",
				cat.dist = 0.03,
				cat.pos = 6
		);
		
		grid.draw(draw.venn.down);
		grid.draw(draw.venn.up);
		dev.off()
	}
	
	
	
	
}
get_union_gene<-function(data_list,out_file_name){
	out_file<-paste(out_file_name,'_common_id.list',sep="")	
	
	n<-length(data_list)
	dd<-data_list[[1]]
	dd
	for (i in 2:n){
		dd<-intersect(dd,data_list[[i]])
	}
	if (length(dd)!=0){
		write.table(as.matrix(dd),file=out_file,quote=F,sep='\t',row.names=F,col.names=F)
	}else{
		cat("No common id!!!!\n")	
	}
	
	
	
	
}



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
		mydata<-read.table(infile,header=F)
		
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
		Usage: Rscript  draw_venn.R -d data1,data2   -p main_titile -n name1,name2 [-i data2,data3 -o Prefix_of_output_files (output)] [-h]
		
		Note:
			multi input files must split by ',', and  labels order of '-n' option must correspond to input files;
			each input files should include one column data and there isn't a header;  
		
		Example usage:
			(1)Rscript draw_venn_5.R -p aaa -d T1_vs_T3_up.id,T1_vs_T2_down.id,T2_vs_T3_down.id,T1_vs_T4_up.id,T3_vs_T4_up.id 
				-i T1_vs_T3_up.id,T1_vs_T2_down.id,T2_vs_T3_down.id,T1_vs_T4_up.id,T3_vs_T4_up.id -n T1_vs_T3,T1_vs_T2,T2_vs_T3,T1-T4,T3-T4 -o test1
			(2)Rscript draw_venn_5.R -p aaa -d T1_vs_T3_up.id,T1_vs_T2_down.id,T2_vs_T3_down.id,T1_vs_T4_up.id,T3_vs_T4_up.id 
				 -n T1_vs_T3,T1_vs_T2,T2_vs_T3,T1-T4,T3-T4 -o test1
		Options:	
			-d datafile,must split by ',' ,max files was 5\n
			-p main titile in Venn Pic
			-i up gene list,optional
			-o outfile name prefix,optional
			-O outdir,optional 
			-n sample titile in Venn Pic
");
	
	
	quit(save="no",status=x)
}
spec<-matrix(c(
				'help','h', 0, "logical",
				'prefix', 'o', 1, "character",
				'mydata', 'd', 1, "character",
				'dataname','n',1, "character",
				'datatype','p',1,"character",
				'mydata_aa','i',1,"character",
				'outdir','O',1,"character"
		), ncol=4, byrow=TRUE)
opt <- getopt(spec);

if (! is.null(opt$help)) { h(0) }
if (is.null(opt$dataname)) { h(0) }
if (is.null(opt$datatype)) { main="venn" }
if (is.null(opt$mydata)) { h(0) }
file.prefix <- ifelse(is.null(opt$prefix), "output", opt$prefix)
#cat(data_paths,"\n")



##########################

if ( is.null(opt$outdir))	{
	pic_name<-paste(file.prefix,"_","venn",".png",sep="")
	opt$outdir=getwd()
}else{
	if(substr(opt$outdir,nchar(opt$outdir),nchar(opt$outdir))=="/"){
		opt$outdir<-substr(opt$outdir,1,nchar(opt$outdir)-1)
	}
	pic_name<-paste(opt$outdir,"/",file.prefix,"_","venn",".png",sep="")
}


main<-paste(opt$datatype,"venn",sep=" ")
#alldata<-read.table(opt$mydata,header=T)

data_list_up<-list()
if (!is.null(opt$mydata_aa)){
	data_list_up<-get_data_list(opt$mydata_aa,opt$dataname)
}






data_list<-get_data_list(opt$mydata,opt$dataname,paste(paste(opt$outdir,"/",file.prefix,"_","ID",".xls",sep="")))
infile_num=length(data_list)
if (infile_num == 2){draw_venn_2(data_list_down=data_list,data_list_up=data_list_up,file_name=pic_name,main=main);}
if (infile_num == 3){draw_venn_3(data_list_down=data_list,data_list_up=data_list_up,file_name=pic_name,main=main);}
if (infile_num == 4){draw_venn_4(data_list_down=data_list,data_list_up=data_list_up,file_name=pic_name,main=main);}
if (infile_num == 5){draw_venn_5(data_list_down=data_list,data_list_up=data_list_up,file_name=pic_name,main=main);}



