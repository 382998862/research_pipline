use Getopt::Long;
use Getopt::Std;
use Config::General;
use FindBin qw($Bin $Script);
use File::Basename qw(basename dirname);
use Cwd qw(abs_path getcwd);
my ($input, $od, $anno);
GetOptions(
        "h|?"           =>\&USAGE,
	"anno:s"	=>\$anno,
	"input:s"     	=>\$input,
	"od:s"		=>\$od,
)or &USAGE;
&USAGE unless ($od and $input and $anno);

`mkdir -p $od`	unless(-d $od);
$od=abs_path($od);
$input=abs_path($input);
$anno=abs_path($anno);

open(ANNO,"$anno/Integrated_Function.annotation.xls")||die $!;
my $header=<ANNO>;chomp($header);my @header=split(/\t/,$header);
shift @header;
my %annos=();
while(<ANNO>){
	chomp;
	my @tmp=split(/\t/,$_);
	my $id=shift @tmp;
	$annos{$id}=join("\t",@tmp);
}
close(ANNO);


my $go_list=(glob("$anno/*.GO.list.txt"))[0];
open(GO,$go_list)||die $!;
my %gos=();
while(<GO>){
	chomp;
	my ($id,$go)=split(/\t/,$_,2);
	$gos{$id}=$go;
}
close(GO);


my $vs=basename $input;
$vs=(split(/\./,$vs))[0];
open(OUT ,">$od/$vs.annotation.xls")||die $!;
open(OUT_GO,">$od/$vs.GO.list.txt")||die $!;
open(IN,$input)||die $!;
my $head=<IN>;chomp($head);
print OUT "$head\t",join("\t",@header),"\n";
while(<IN>){
	chomp;
	my @tmp=split(/\t/,$_);	
	if(exists $annos{$tmp[0]}){
		print OUT join("\t",@tmp),"\t$annos{$tmp[0]}\n";
	}
	if(exists $gos{$tmp[0]}){
		print OUT_GO "$tmp[0]\t$gos{$tmp[0]}\n";
	}
}
close(IN);
close(OUT);
close(OUT_GO);


print "perl $Bin/anno/gene_ontology_graph_v1.2.pl -i $go_list -i $od/$vs.GO.list.txt -mark All -mark DE -o $od -k $vs\n";
system "perl $Bin/anno/gene_ontology_graph_v1.2.pl -i $go_list -i $od/$vs.GO.list.txt -mark All -mark DE -o $od -k $vs";

sub USAGE{
	my $usage=<<"USAGE";
Program: $0
Version: 1.0.0
Contact: wenyh\@biomarker.com.cn
Usage:
	Options:
	-input	<file>	input DEG file,forced
	-od	<dir>	output dir,forced
	-anno	<dir>	anno dir, Result 

	-h	Help

Example:

USAGE
	print $usage;
	exit;
}


