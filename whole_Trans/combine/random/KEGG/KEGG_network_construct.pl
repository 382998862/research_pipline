#!/usr/bin/perl -w
use Getopt::Long;
use Getopt::Std;
use Config::General;
use FindBin qw($Bin $Script);
use Cwd qw(abs_path getcwd);
my ($gene,$ko,$od,$symbol,$fdr,$min,$key,$prefix);
GetOptions(
        "h|?"           =>\&USAGE,
	"gene:s"	=>\$gene,
	"ko:s"		=>\$ko,
	"od:s"		=>\$od,
	"symbol:s"	=>\$symbol,
	"fdr:s"		=>\$fdr,
	"min:s"		=>\$min,
	"max:s"		=>\$max,
	"key:s"		=>\$key,
	"prefix:s"	=>\$prefix,
)or &USAGE;
&USAGE unless ( $ko and $gene );

$fdr||=0.05;
$max||=10;
$min||=2;

`mkdir -p $od`	unless(-d $od);

####################get kegg pathways to intergrate
my @paths=();
if(-e $ko){
	open(KO,$ko)||die $!;
	while(<KO>){
		chomp;
		next if($_=~/^#/);
		my @tmp=split(/\t/,$_);
		if($tmp[5]<$fdr && @paths<$max){
			push @paths,$tmp[1];
		}
		if($tmp[5]>=$fdr && @paths<$min){
			push @paths,$tmp[1];
		}
	}
	close(KO);
}else{
	@paths=split(",",$ko);
}
die "KEGG pathway shoud exist at least one!\n"	if(@paths<1);

###################get gene and KO relation
my %symbol=();
if(defined $symbol){
	open(SYM,$symbol)||die $!;
	while(<SYM>){
		chomp;
		next if($_=~/^#/);
		my ($g,$s)=split(/\t/,$_);	
		$symbol{$g}=$s;	##$symbol{ENSG}=symbol
	}
	close(SYM);
}

my %K=();
open(GENE,$gene)||die $!;
while(<GENE>){
	chomp;
	next if($_=~/^#/);
	my ($id,$info)=split(/\t/,$_);
	my $kid=(split(/\|/,$info))[0];
	$K{$kid}=$id;	##$K{K50581}=ENSG
	##if given symbol, the ko id would be changed to gene symbol rather than Ensembl id
	$K{$kid}=$symbol{$id}	if(exists $symbol{$id});  ##$K{K505874}=symbol
}
close(GENE);
###########################key RNA 
my %keys=();
if(defined $key){
	open(KEY,"$key")||die $!;
	while(<KEY>){
		chomp;
		my $id=(split(/\s+/,$_))[0];
		next if($id eq "--");
		$keys{$id}++;
		$keys{$symbol{$id}}++	if(exists $symbol{$id});
	}
	close(KEY);
}
##########################get relation from kegg pathways
mkdir "$od/KEGG"	unless(-d "$od/KEGG");
foreach my $p(@paths){
	`cp /share/nas2/database/KEGG/201703/kegg_png/$p.kgml $od/KEGG`;
}

mkdir "$od/relation"        unless(-d "$od/relation");
print "Rscript KEGG_network_extract.r -i $od/KEGG -o $od/relation\n";
`Rscript $Bin/KEGG_network_extract.r -i $od/KEGG -o $od/relation`;


##########################get info from relation
my %entry=();	##$entry{path}{id}
open(ENTRY,"$od/relation/entry.txt") or die $! ;
while(<ENTRY>){
	chomp;
	next if($_=~/^id/);
	my @tmp=split(/\t/,$_);	##id	ko:K11262	ortholog	path
	$entry{$tmp[3]}{$tmp[0]}{type}=$tmp[2];
	my @kos=split(/\s+/,$tmp[1]);
	my $id="$tmp[3].$tmp[0]";
	foreach my $k(@kos){
		push @{$id},(split(/:/,$k))[1];
	}	##@{path.103}
}
close(ENTRY);

open(RELATION,"$od/relation/relation.txt")||die $!;
open(OUT,">$od/network_tmp.txt")||die $!;
open(NODE,">$od/node_tmp.txt")||die $!;

while(<RELATION>){
	chomp;
	next if($_=~/^path_name/);
	my ($path,$type,$en1,$en2,$subtype)=split(/\t/,$_);	
	my $type1=$entry{$path}{$en1}{type};
	my $type2=$entry{$path}{$en2}{type};
	next if($type1 eq "compound"	|| $type2 eq "compound");
	if($type1 eq "group" ||$type2 eq "group"){	## KEGG_network_extract.r will take apart group into component
		print $_,"\n";
		die "Something must be wrong in advanced process!\n";
	}
	#############
	my @entry1=();my @entry2=();
	if($type1 ne "map"){
		@entry1=&convertKtoGene(\@{"$path.$en1"});
	}else{	@entry1=@{"$path,$en1"};}
	if($type2 ne "map"){
		@entry2=&convertKtoGene(\@{"$path.$en2"});
	}else{	@entry2=@{"$path.$en2"};	
	}
	#############
	next if(@entry1==0 ||@entry2==0);
	my $entry_id1=join(",",@entry1);
	my $entry_id2=join(",",@entry2);
#	next if($entry_id1 eq $entry_id2);	##remove 自成环的边
	my ($lty,$width,$arrow,$lab)=&defineEdge($subtype);
	foreach my $e1(@entry1){
		if($type1 eq "map"){
			print NODE "$e1\trectangle\tblue\n";
		}else{	
			print NODE "$e1\tcircle\t";
			if(exists $keys{$e1}){	print NODE "red\n";}else{print NODE "gray\n";}
			
		}
		foreach my $e2(@entry2){
			if($type2 eq "map"){
				print NODE "$e2\trectangle\tblue\n";
			}else{	
				print NODE "$e2\tcircle\t";
				if(exists $keys{$e2}){  print NODE "red\n";}else{print NODE "gray\n";}
			}
			print OUT "$path\t$e1\t$e2\t$type\t$subtype\t$lty\t$width\t$arrow\t$lab\n";
		}
	}
}
close(OUT);
close(NODE);
close(RELATION);

`sort $od/network_tmp.txt|uniq >$od/${prefix}_network.xls && sort $od/node_tmp.txt |uniq>$od/${prefix}_node.xls`;
`sed -i "1ipath\tentry1\tentry2\ttype\tsubtype\tlty\twidth\tarrow\tlab" $od/${prefix}_network.xls`;
`sed -i "1iGene\tshape\tcolor" $od/${prefix}_node.xls`;

`rm $od/network_tmp.txt && rm $od/node_tmp.txt`;
print "Rscript $Bin/../bin/iGraph_network.r --edge $od/${prefix}_network.xls --odir $od --node $od/${prefix}_node.xls --node.color 3 --node.shape 2 --from 2 --to 3  --direction TRUE --edge.lty 6 --edge.width 7 --arrow.mode 8 --edge.label 9 --edge.color 1 --set.seed 200 --key $prefix \n";
`Rscript $Bin/../bin/iGraph_network.r --edge $od/${prefix}_network.xls --odir $od --node $od/${prefix}_node.xls --node.color 3 --node.shape 2 --from 2 --to 3  --direction TRUE --edge.lty 6 --edge.width 7 --arrow.mode 8 --edge.label 9 --edge.color 1 --set.seed 200 --key $prefix`;

#######################################
sub convertKtoGene{
	my $s=shift;
	my @set=@{$s};
	my @new=();
	foreach my $s(@set){
		if(exists $K{$s}){
			push @new,$K{$s};
		}
	}
	return @new;
}

################根据subtype定义边的类型，粗细以及箭头方向,以及边上的字
#0 for no edges, 1 for solid lines,2 for dashed, 3 for dotted, 4 for dotdash, 5 for longdash, 6 for twodash. 
#0 means no arrows(-), 1 means backward arrows(<-), 2 is for forward arrows(->) and 3 for both(<->)
sub defineEdge{
	my $subtype=shift;
	my ($lty,$width,$arrow,$label)=(1,2,2,"");
	if($subtype =~/compound/){
		$lty=1;
		$arrow=0;
	}
	$width=1	if($subtype =~/hidden|dissociation|missing/);
	$lty=2		if($subtype=~/activation|expression/);
	$lty=5		if($subtype=~/inhibition|repression/);
	$lty=3		if($subtype=~/indirect effect/);
	$lty=4          if($subtype=~/indirect effect/ && $subtype=~/activation/);
	if($subtype=~/state change/){
		$lty=3;$arrow=0;
	}
	if($subtype=~/binding\/association/){
		$lty=2;$arrow=0;
	}
	$label="+p"	if($subtype=~/phosphorylation/);
	$label="-p"     if($subtype=~/dephosphorylation/);
        $label="+g"     if($subtype=~/glycosylation/);
        $label="+u"     if($subtype=~/ubiquitination/);
        $label="+m"     if($subtype=~/methylation/);
	return ($lty,$width,$arrow,$label);
}


sub readConfig{
	my $configFile=shift;
	my $d=Config::General->new(-ConfigFile => "$configFile");
	my %config=$d->getall;	
	return %config;
}
sub USAGE{
	my $usage=<<"USAGE";
Program: $0
Version: 1.0.0
Contact: wenyh\@biomarker.com.cn
Usage:
	Options:
	-gene	<file>	gene and KO id relation (two column)
			Format:ENSG000004545874 K021547
	-symbol	<file>	gene and symbol relation (two column)
			Format:ENSG000004545874 BRAT1
	-od	<path>	output path
	-ko	<file>	kegg list file or ko00040,ko00042,ko00054
	####if provide kegg list file  ##pathwayname ko_id Cluster_fre Genome_frequency P-value fdr rich_factor
	-fdr	<float>	default 0.05
	-min	<int>	min pathway should be 2
	-max	<int>	max pathway should be 10

	-key	<file>	key RNA outstanding in the network by red color	



#####description
rectangle/blue	stand for other KEGG pathway
circle/gray	stand for gene
circle/red	stand for key gene

	
	-h	Help

Example:

USAGE
	print $usage;
	exit;
}


