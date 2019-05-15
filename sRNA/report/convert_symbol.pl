use Getopt::Long;
use Getopt::Std;
use Config::General;
use FindBin qw($Bin $Script);
use Cwd qw(abs_path getcwd);
my ($input,$id,$key);
GetOptions(
        "h|?"           =>\&USAGE,
	"i:s"		=>\$input,
	"id:s"		=>\$id,
	"key:s"		=>\$key,
)or &USAGE;
&USAGE unless ($input and $id);

my @in=split(/,/,$input);
$key||="ENS";
print "key is:$key:\n";
open(REL,$id)||die $!;
my %relate=();
while(<REL>){
	chomp;
	next if($_=~/^Gene|^#/);
	my @tmp=split(/\s+/,$_);
	$relate{$tmp[0]}=$tmp[1];
	$relate{"gene:$tmp[0]"}=$tmp[1];
	$relate{"gene_$tmp[0]"}=$tmp[1];
	$relate{"lncRNA:$tmp[0]"}="lncRNA:$tmp[1]";
}
close(REL);


foreach my $file(@in){
	$file=abs_path($file);
	my $ok = "gene:$key";
	my $lok = "lncRNA:$key";
	my $line=`grep -E '$key|$ok|$lok' $file|head -n 1`;
	chomp($line);
	next if($line eq "");

	print $file,"\n";
	my @temp=split(/\t/,$line);
	my $flag;my @colnumn;
	for (my $i=0;$i<@temp;$i++){
		if($temp[$i]=~/$key/){
			$flag=$i;
			push @colnumn,$flag;
			print "Flag is :$flag:\n";
			last;
		}
	}
	open(IN,$file)||die $!;
	my $header=<IN>;
	chomp($header);
	my @headers=split(/\t/,$header); 	
	open(OUT,">$file.tmp")||die $!;
#	foreach my $col(@colnumn){
#		$col+=1;
	open(IN,$file)||die $!;
	my $header=<IN>;
	chomp($header);
	my @headers=split(/\t/,$header);
	open(OUT,">$file.tmp")||die $!;	
	if($headers[$flag]=~/$key/){
		my @symbols=&tranSymbols($headers[$flag]);	
		$headers[$flag] .="\t".join(";",@symbols);
	}else{
		$headers[$flag] .="\tSymbol";
	}

	print OUT join("\t",@headers),"\n";
	while(<IN>){
		chomp;
		my @tmp=split(/\t/,$_);
		my $ensgs=$tmp[$flag];
		my @symbols=&tranSymbols($ensgs);
		$tmp[$flag] .="\t".join(";",@symbols);
		print OUT join("\t",@tmp),"\n";
	}	
	close(IN);
	close(OUT);
	`mv $file.tmp $file`;
#	}
}

close(LOG);

`rm $log`	if(-z $log);
sub tranSymbols{
	my $ensgs=shift;
	$ensgs =~ s/;|\//,/g;
	my @genes=split(/,/,$ensgs);
	my @symbols=();
	foreach my $g (@genes){
		if(exists $relate{$g}){
			push @symbols,$relate{$g};
		}else{
			push @symbols,$g;
		}
	}
	return @symbols;
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
	-i	<file>	input file list, separated by comma			forced
	-log	<file>	log file, record the files which convert ensembl id to symbol, not must	
	-id	<file>	ensembl and symbol id relationship file
			eg:/share/nas2/database/genome/Homo_sapiens/GRCh37/id_name.list
	-key	<str>	ENS or MSTRG
	-h	Help

Example:

USAGE
	print $usage;
	exit;
}


