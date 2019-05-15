use Getopt::Long;
use Getopt::Std;
use Config::General;
use FindBin qw($Bin $Script);
use Cwd qw(abs_path getcwd);
my ($net,$seed,$od,$node,$type);
GetOptions(
        "h|?"           =>\&USAGE,
	"net:s"		=>\$net,
	"node:s"	=>\$node,	
	"seed:s"	=>\$seed,	
	"od:s"		=>\$od,
	"key:s"		=>\$key,
	"type:s"	=>\$type,
)or &USAGE;
&USAGE unless ($net and $seed and $od);

$key||="sub_network";
$type||="onestep";
$net=abs_path($net);
$seed=abs_path($seed);
`mkdir -p $od`	unless(-d $od);
$od=abs_path($od);

my %seeds=();
open(SEED,$seed)||die $!;
while(<SEED>){
	chomp;
	next if($_=~/^#/);
	my @tmp=split(/\s+/,$_);
	$seeds{$tmp[0]}++;
}
close(SEED);

open(NET,$net)||die $!;
open(OUT,">$od/$key.txt")||die $!;
my %newnet=();
while(<NET>){
	chomp;
	my @tmp=split(/\t/,$_);
	if($type eq "onestep"){
		if(exists $seeds{$tmp[0]} || exists $seeds{$tmp[1]}){
			print OUT "$tmp[0]\t$tmp[1]\n";	
			$newnet{$tmp[0]}++;
			$newnet{$tmp[1]}++;
		}
	}elsif($type eq "only"){
		if(exists $seeds{$tmp[0]} && exists $seeds{$tmp[1]}){
			print OUT "$tmp[0]\t$tmp[1]\n";
			$newnet{$tmp[0]}++;
			$newnet{$tmp[1]}++;
		}
	}	
}
close(OUT);
close(NET);
`sed -i "1inode1\tnode2\n" $od/$key.txt`;

if(defined $node){
	my %nodes=();
	open(NODE,$node)||die $!;
	open(ONODE,">$od/$key.node")||die $!;
	while(<NODE>){
		chomp;
		my $point=(split(/\s+/,$_))[0];
		$nodes{$point}=$_;
	}
	foreach my $n(keys %nodes){
		if(exists $newnet{$n}){
			print ONODE "$nodes{$n}\t";
			if(exists $seeds{$n}){print ONODE "red\n";}else{print ONODE "gray\n";}
		}
	}
	
	close(ONODE);
	close(NODE);
	`sed -i "1inode\ttype\tcolor" $od/$key.node`;
}


my $net_line=`wc -l $od/$key.txt`;
chomp($net_line);
if($net_line>2){
	my $plot="Rscript $Bin/iGraph_network.r --edge $od/$key.txt --odir $od --set.seed 200 --key $key --direction FLASE";
	$plot .= " --node $od/$key.node --node.type 2 --node.color 3 " if(defined $node);
	print "$plot\n";
	system($plot);
}else{
	my $plot="Rscript $Bin/iGraph_network.r --edge $net --odir $od --set.seed 200 --key $key --direction FLASE";
	print "$plot";
	system($plot);	
}


sub USAGE{
	my $usage=<<"USAGE";
Program: $0
Version: 1.0.0
Contact: wenyh\@biomarker.com.cn
Usage:
	Options:
	-seed	<file>	forced		seed node to extract from network based on one step neighbour
	-net	<file>	forced		backgroud network, the first two column is source/target node
	-od	<file>	forced		output path,contained one setp neighbour network
	-node	<file>	not forced	node file, if provide will generate new node file
	-key	<str>	not forced	key word	

	-typ	default onestep
		onestep:this parameter represents onestep network
		only:	this parameter represents sub network source/target node must both come frome seed rna

	-h	Help

Example:

USAGE
	print $usage;
	exit;
}


