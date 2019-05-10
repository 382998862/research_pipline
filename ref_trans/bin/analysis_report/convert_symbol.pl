use Getopt::Long;
use Getopt::Std;
use Config::General;
use FindBin qw($Bin $Script);
use Cwd qw(abs_path getcwd);
my ($input,$db);
GetOptions(
        "h|?"           =>\&USAGE,
	"i:s"		=>\$input,
	"db:s"		=>\$db,
	"id:s"		=>\$id,
)or &USAGE;
&USAGE unless ($input);
my $key="ENS";

my @in=();
if(-e $input){
	open(LIST,$input)||die $!;
	while(<LIST>){
		chomp;	
		push @in,$_	if(-e $_);
	}
		close(LIST);
}else{
	@in=split(/,/,$input);
}
if(!defined $id){
	if(defined $db){
		my %db_key=(
			"GRCh37"	=>"ENSG",
			"GRCh38"	=>"ENSG",
			"GRCm38"	=>"ENSMUSG",
			"mm9"		=>"ENSMUSG",
			"Rnor_5.0"	=>"ENSRNOG",
			"Rnor_6.0"      =>"ENSRNOG",
		);
		$key=$db_key{$db};
#		my $database="$Bin/database.config";
#		my %dataconfig=&readConfig($database);
#		my $symbolPath=$dataconfig{symbolPath};
		my $symbolPath="/share/nas1/wenyh/develop/lncRNA/medical_lncRNA_v1.2/anno/bin/symbol";
		$id="$symbolPath/$db.txt";
	}else{
		print "One of parameter db or id must exist!\n";
		die;
	}
}
open(REL,$id)||die $!;
my %relate=();
while(<REL>){
	chomp;
	next if($_=~/^Gene|^#/);
	my @tmp=split(/\t/,$_);
	$relate{$tmp[0]}=$tmp[1];
	$relate{"gene:$tmp[0]"}=$tmp[1];
	$relate{"gene_$tmp[0]"}=$tmp[1];
}
close(REL);


foreach my $file(@in){
	$file=abs_path($file);
	my $line=`grep $key $file|head -n 1`;
	chomp($line);
	next if($line eq "");

	print $file,"\n";
	my @temp=split(/\t/,$line);
	my $flag;
	for (my $i=0;$i<@temp;$i++){
		if($temp[$i]=~/$key/){
			$flag=$i;
			last;
		}
	}
	open(IN,$file)||die $!;
	my $header=<IN>;
	chomp($header);
	my @headers=split(/\t/,$header); 	
	open(OUT,">$file.tmp")||die $!;
	if($headers[$flag]=~/$key/){
		my @symbols=&tranSymbols($headers[$flag]);	
		$headers[$flag] .="\t".join(";",@symbols);
	}else{
		$headers[$flag] .="\tgene_symbol";
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
}

close(LOG);

`rm $log`	if(-z $log);
sub tranSymbols{
	my $ensgs=shift;
	$ensgs =~ s/;/,/g;
	my @genes=split(/,|\;|\//,$ensgs);
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
			if too much files, you could give a list contained all files, one line one file.

	-db	<str>	ref genome: GRCh37,GRCh38,GRCm38,mm9,Rnor_5.0,Rnor_6.0	
			one of db and id must exist
	-log	<file>	log file, record the files which convert ensembl id to symbol, not must	
	-id	<file>	ensembl and symbol id relationship file
			eg:/share/nas2/database/genome/Homo_sapiens/GRCh37/id_name.list

	-h	Help

Example:

USAGE
	print $usage;
	exit;
}


