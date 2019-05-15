use Getopt::Long;
use Getopt::Std;
use Config::General;
use FindBin qw($Bin $Script);
use newPerlBase;
use Cwd qw(abs_path getcwd);
my ($input,$output,$FDR,$PValue,$FC);
GetOptions(
        "h|?"           =>\&USAGE,
	"input:s"	=>\$input,
	"FDR:s"		=>\$FDR,
	"PValue:s"	=>\$PValue,
	"FC:s"		=>\$FC,
	"output:s"	=>\$output,
)or &USAGE;
&USAGE unless ($input);

if(!defined $PValue && !defined $FDR){
	$FDR||=0.05;
}

$FC||=2;
my $logfc=log($FC)/log(2);

$input=abs_path($input);
$output||="$input.DEG_final.xls";
$output=~s/\.all//;

open(IN,$input)||die $!;
open(OUT,">$output")||die $!;
my $head=<IN>;chomp($head);
my @header=split(/\t/,$head);
my %index=();
for(my $i=0;$i<@header;$i++){
	$index{FDR}=$i	if($header[$i] eq "FDR");
	$index{PValue}=$i  if($header[$i] eq "PValue");
	$index{log2FC}=$i  if($header[$i] eq "log2FC");
}

my $del;
if(!defined $PValue && exists $index{PValue}){
	splice @header,$index{PValue},1;
	$del=$index{PValue};
}
if(!defined $FDR && exists $index{FDR}){
        splice @header,$index{FDR},1;
	$del=$index{FDR};
}

if(defined $PValue && !exists $index{PValue}){
	print "Parameter PValue is given! While the header of $input does not have the PValue column!\n";
}

print OUT join("\t",@header),"\tregulated\n";
while(<IN>){
	chomp;
	my @tmp=split(/\t/,$_);
	next if(abs($tmp[$index{log2FC}])<$logfc);
	if(defined $FDR){
		next if($tmp[$index{FDR}] > $FDR);
	}
	if(defined $PValue){
		next if($tmp[$index{PValue}] >$PValue);
	}
	my $dir="up";
	$dir="down"	if($tmp[$index{log2FC}] <0);
	if(defined $del){
		splice @tmp,$del,1;
	}

	print OUT join("\t",@tmp),"\t$dir\n";
}
close(IN);
close(OUT);


#################
sub readConfig{
	my $configFile=shift;
	my $d=Config::General->new(-ConfigFile => "$configFile");
	my %config=$d->getall;	
	return %config;
}
sub qsub()
{
        my $shfile= shift;
        my $cmd = "sh /share/nas2/genome/bmksoft/tool/qsub_sge_plus/v1.0/qsub_sge.plus.sh --reqsub $shfile --independent";
        &run_or_die($cmd);              
        return ;
}
sub run_or_die()
{
        my ($cmd) = @_ ;
        &show_log($cmd);
        my $flag = system($cmd) ;
        if ($flag != 0){
                &show_log("Error: command fail: $cmd");
                exit(1);
        }
        &show_log("done.");
        return ;
}
sub show_log()
{
        my ($txt) = @_ ;
        my $time = time();
        my($sec, $min, $hour, $day, $mon, $year, $wday, $yday, $isdst) = localtime($time);
        $wday = $yday = $isdst = 0;
        my $Time=sprintf("%4d-%02d-%02d %02d:%02d:%02d", $year+1900, $mon+1, $day, $hour, $min, $sec);
        print "$Time:\t$txt\n" ;
        return ($time) ;
}

sub USAGE{
	my $usage=<<"USAGE";
Program: $0
Version: 1.0.0
Contact: wenyh\@biomarker.com.cn
Usage:
	
	Options:
	-h	Help
	-input	file	forced		input deg.all file
	-output	file	forced		output deg_final.xls
	-FDR	double	not must	default not used
	-PValue	double	not must	default not used
					Note: if both FDR and PValue not give, will default us FDR 0.05
	-FC	double	not must	default 2
Example:

USAGE
	print $usage;
	exit;
}


