use Getopt::Long;
use Getopt::Std;
use Config::General;
use FindBin qw($Bin $Script);
use newPerlBase;
use File::Basename;
use Cwd qw(abs_path getcwd);
my ($deg, $outpath, $db);
GetOptions(
        "h|?"           =>\&USAGE,
	"deg:s"		=>\$deg,
	"db:s"		=>\$db,
	"o:s"     	=>\$outpath,
)or &USAGE;
&USAGE unless ($outpath and $deg and $db);

$db ||="GRCh38";
my %config=&readConfig($input);
my $basename=basename $deg;
my $key=(split(/\./,$basename))[0];

`mkdir -p $outpath`	if(!-d $outpath);
$outpath=abs_path($outpath);

####################gene symbol ###########################

open(SYMBOL,"$Bin/bin/symbol/$db.txt")||die $!;
<SYMBOL>;
my %relation=();
while(<SYMBOL>){
	chomp;
	my @tmp=split;
	$relation{$tmp[0]}=$tmp[1];
	$relation{$tmp[1]}=$tmp[0];
}
close(SYMBOL);

#########################  DEG  ###########################

open(DEG,$deg)||die $!;
my $deghead=<DEG>;
chomp($deghead);
my @tmp=split(/\s+/,$deghead);
my $length=@tmp;
# my $samplenum=($length-4)/2;  #适用于有参
my $samplenum=$length-4;

my $head="ENSG_ID\tgene symbol\tlog2FC\tregulated";
#for my $i(1..$samplenum) {$head .="\t$tmp[$i]";}
#$head .="\tlog2FC";

my %degen=();
while(<DEG>){
	chomp;
	my @tmp=split;
	$degen{$tmp[0]}  =$relation{$tmp[0]};
=head
	for(my $i=0;$i<$samplenum;$i++){
		$degen{$tmp[0]} .="\t$tmp[$i+1]";
	}
=cut
	$degen{$tmp[0]} .="\t$tmp[-2]\t$tmp[-1]";
}
close(DEG);

##################### cosmic anno (only do for human) ######

if($db eq "GRCh37" || $db eq "GRCh38"){
	open(COSMIC,"$Bin/bin/$db.cancer_gene_census.txt")||die $!;
	my $coshead=<COSMIC>;;
	my %cosmic=();
	while (<COSMIC>) {
		chomp;
		my @tmp=split(/\t/,$_);
		$cosmic{$tmp[0]}="$tmp[1]\t$tmp[7]\t$tmp[8]";
	}
	close(COSMIC);
	open(OUT,">$outpath/$key.DEG_cosmic.xls")||die $!;
	print OUT "$head\tDescription\tTumour_Types(Somatic)\tTumour_Types(Germline)\n";
	foreach my $k(keys %degen){
		if(exists $relation{$k} && exists $cosmic{$relation{$k}}){
			print OUT "$k\t$degen{$k}\t$cosmic{$relation{$k}}\n";
		}
	}
	close(OUT);
}
######################### TF anno ##########################

open(OUT,">$outpath/$key.DEG_TF.xls")||die $!;
open(TF,"$Bin/bin/animalTFDB3.0/Animal_TFDB3.0.db")||die $!;
my $tfhead=<TF>;
print OUT "$head\tFamily\n";
while(<TF>){
	chomp;
	my @tmp=split(/\t/,$_);
	print OUT "$tmp[0]\t$degen{$tmp[0]}\t$tmp[3]\n"	if(exists $degen{$tmp[0]});
}
close(TF);
close(OUT);


#################################################
#	read config
#################################################

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
Date: 2016-11-18
Contact: wenyh\@biomarker.com.cn
Usage: Do DEG gene cosmic (only for human) and animalTFBS annotation.

	Options:
	-deg	<file>	input DEG file,forced## T04_T05_T06_vs_T01_T02_T03.DEG.final.xls
	-o	<str>	output path,forced
	-db	<str>	optional GRCh37,GRCh38,GRCm38,mm9,Rn_Celera,Rnor_5.0,Rnor_6.0
			default GRCh37

	-h	Help

Example:

USAGE
	print $usage;
	exit;
}


