library(getopt)
spec=matrix(c(
	'help','h',0,"logical",
	'edge','e',1,"character",
	'odir','o',1,"character",

	'node','n',2,"character",
	'node.type','T',2,"integer",
	'node.size','S',2,"integer",
	'node.color','C',2,"integer",
	'node.shape','s',2,"integer",

	'from','f',2,"integer",
	'to','t',2,"integer",
	'direction','d',2,"character",
	'edge.lty','y',2,"integer",
	'edge.label','l',2,"integer",
	'edge.width','w',2,"integer",
	'edge.color','r',2,"integer",
	'edge.labcol','b',2,"character",
	'arrow.mode','a',2,"integer",

#	'layout','L',2,"character",
	'key','k',2,"character",
	'set.seed','B',2,"integer"
	),byrow=TRUE,ncol=4
)
opt=getopt(spec)
print_usage<-function(spec=NULL){
	cat(getopt(spec,usage=TRUE));
	cat("Usage example: \n")
	cat("
		Description: network plot
		Contact: Yanhua Wen <wenyh@biomarker>	
		Usage example:
			Rscript iGraph_network.r --edge edge.txt --node node.txt --odir outpath

		Note: The node in node.txt must be equal to the node in edge.txt
	Forced:
		 --help		NULL	get this help
		 --edge		str	edge file, must contain source/target node,need header
		 --odir		str	output path

	Opitional:

		Plot node parameters:
		 --node		str	node file, not must, need header
					if exists, must contail all nodes in edge, the first column is node

		 --node.type	int	in node.txt which column set for node type (lncRNA/circRNA/miRNA/mRNA)
					will define node shape based on RNA type
		 --node.size	int	in node.txt which column set for node size, if none the node size ranged with degree
		 --node.color	int	in node.txt which column set for node color(red/blue)
		 --node.shape	int	in node.txt which column set for node shape(circle/rectangle/crectangle/sphere)
		 

		Plot edge parameters:
		 --from		int	in edge.txt which column set for the source node,default 1
		 --to		int	in edge.txt which column set for the target node,default 2
		 --direction	str	whether edge in network have direction, FALSE/TRUE, default FALSE
		 --edge.lty		int	in edge.txt which column set for edge type
					0 for no edges, 1 for solid lines,2 for dashed, 3 for dotted, 4 for dotdash, 5 for longdash, 6 for twodash. 
		 --edge.width	int	in edge.txt which column set for edge width 
		 --edge.label	int	in edge.txt which column set for edge label 
		 --edge.color	int	in edge.txt which column set for edge color,same with edge.labcol
		 --edge.labcol	str	label color,default blue
		 --arrow.mode	int	in edge.txt which column set for arrow mode
					0 means no arrows(-), 1 means backward arrows(<-), 2 is for forward arrows(->) and 3 for both(<->)

		Network parameters:
		 --key		str	outputfile key, default NA
#		 --layout	str	layout method, default 	layout.auto
					layout、layout.auto、layout.bipartite、layout.circle、layout.drl、layout.fruchterman.reingold、layout.fruchterman.reingold.grid、layout.graphopt、layout.grid、layout.grid.3d、layout.kamada.kawai、layout.lgl、layout.mds、layout.merge、layout.norm、layout.random、layout.reingold.tilford、layout.sphere、layout.spring、layout.star、layout.sugiyama、layout.svd

		 --set.seed	int	default 200, meaning the figure will not be changed 
		 			if set as 0, the figure will be changed

Formula:	node.txt
BRCA1	mRNA	circle	
PTEN	lncRNA	circle	blue

Formula:	edge.txt
BRCA1	PTEN	25	target

		Rscript iGraph_network.r --edge edge.txt --from 1 --to 2 --direction TRUE --edge.width 4 --edge.label 3  --edge.label.color red --node node.txt --node.type 2 --key test --set.seed 200 --odir ./
		\n")

	q(status=1);
}

if(is.null(opt$edge)) print_usage(spec)
if(is.null(opt$odir)) print_usage(spec)


if(is.null(opt$direction))	opt$direction="FALSE"
if(is.null(opt$from))		opt$from=1
if(is.null(opt$to))		opt$to=2
if(is.null(opt$edge.labcol))	opt$edge.labcol="blue"

if(is.null(opt$set.seed))	opt$set.seed=200
if(is.null(opt$key))		opt$key=""


range<-function(values,size){
	diff<-max(values)-min(values)
	if(diff==0){
		return(size)
	}
	if(diff!=0){
		newvalues=size+(values-min(values))/diff*size
		return(newvalues)
	}
}


library(igraph)
read.csv(opt$edge,sep="\t",header=T)->net
edge<-data.frame(net[,opt$from],net[,opt$to])


shapes<-c("circle","rectangle","crectangle","sphere")
col=rainbow(4)
# mRNA, circRNA, lncRNA, miRNA
if(!is.null(opt$node)){
	read.csv(opt$node,sep="\t",header=T)->node
	g <- graph.data.frame(edge, directed = opt$direction, vertices=as.character(node[,1]))

	if(!is.null(opt$node.type)){	
		V(g)$type=as.character(node[,opt$node.type])
	
		V(g)[type=="mRNA"]$shape	=shapes[1]
		V(g)[type=="lncRNA"]$shape	=shapes[2]
		V(g)[type=="miRNA"]$shape	=shapes[3]
		V(g)[type=="circRNA"]$shape	=shapes[4]
	}

	if(!is.null(opt$node.color)){
		V(g)$color=as.character(node[,opt$node.color])
	}	

	if(!is.null(opt$node.shape)){
		V(g)$shape=as.character(node[,opt$node.shape])
	}
}

#############node 
if(is.null(opt$node)){
	g <- graph.data.frame(edge, directed = opt$direction)
	V(g)$color="pink"
	V(g)$shape="circle"
}


############edge
	##edge.lty
E(g)$lty=1
if(!is.null(opt$edge.lty)){
	E(g)$lty=as.numeric(as.character(net[,opt$edge.lty]))
}
	##edge.label
if(!is.null(opt$edge.label)){
	E(g)$label=as.character(net[,opt$edge.label])
}
	##edge.color
E(g)$color="gray"
if(!is.null(opt$edge.color)){
	E(g)$color=as.character(net[,opt$edge.color])

	unique(as.character(net[,opt$edge.color]))->color
	intersect(unique(substr(color,1,1)),"#")->Fcolor
	if(length(intersect(color,colors()))==0){
		if(length(Fcolor)==0){
			newcol=rainbow(length(color))
			for(i in 1:length(color)){
				which(E(g)$color==color[i])->s
				E(g)$color[s]=newcol[i]
			}
		}
	}
}

	###edge.width
E(g)$width=1.5
if(!is.null(opt$edge.width)){
	score=as.numeric(as.character(net[,opt$edge.width]))
	E(g)$width=range(values=score,size=1.5)
}

	####arrow.mode
E(g)$arrow.mode=2
if(!is.null(opt$arrow.mode)){
	E(g)$arrow.mode=as.numeric(as.character(net[,opt$arrow.mode]))
}




###################
V(g)$degree=degree(g)
V(g)$size=range(V(g)$degree,2)
V(g)$labcex=range(V(g)$degree,0.2)
if(!is.null(opt$node.size)){
	nodesize=as.numeric(as.character(net[,opt$node.size]))
	V(g)$size=range(nodesize,10)
	V(g)$labcex=range(nodesize,0.5)
}

if(opt$set.seed>0){
	set.seed(opt$set.seed)
}

pngfile<-paste(opt$odir,"/",opt$key,"_network.png",sep="")
pdffile<-paste(opt$odir,"/",opt$key,"_network.pdf",sep="")
#layout、layout.auto、layout.bipartite、layout.circle、layout.drl、layout.fruchterman.reingold、layout.fruchterman.reingold.grid、layout.graphopt、layout.grid、layout.grid.3d、layout.kamada.kawai、layout.lgl、layout.mds、layout.merge、layout.norm、layout.random、layout.reingold.tilford、layout.sphere、layout.spring、layout.star、layout.sugiyama、layout.svd
png(pngfile,height=3000,width=3000,res=300)
plot(	g,	
	edge.arrow.size=0.2,edge.arrow.mode=E(g)$arrow.mode,
	layout=layout.fruchterman.reingold.grid,#layout.kamada.kawai,
	vertex.size=V(g)$size, vertex.color=V(g)$color, vertex.shape=V(g)$shape,vertex.frame.color="gray",
	vertex.label=V(g)$name,vertex.label.cex=V(g)$labcex, vertex.label.color='black',
	edge.label=E(g)$label,edge.label.cex=0.2, edge.label.color=opt$edge.labcol,edge.width=E(g)$width,edge.lty=E(g)$lty
)
dev.off()

pdf(pdffile)
plot(	g,	
	edge.arrow.size=0.2,edge.arrow.mode=E(g)$arrow.mode,
	layout=layout.fruchterman.reingold.grid,#layout.kamada.kawai,
	vertex.size=V(g)$size, vertex.color=V(g)$color, vertex.shape=V(g)$shape,vertex.frame.color="gray",
	vertex.label=V(g)$name,vertex.label.cex=V(g)$labcex, vertex.label.color='black',
	edge.label=E(g)$label,edge.label.cex=0.2, edge.label.color=opt$edge.labcol,edge.width=E(g)$width,edge.lty=E(g)$lty
)
dev.off()


