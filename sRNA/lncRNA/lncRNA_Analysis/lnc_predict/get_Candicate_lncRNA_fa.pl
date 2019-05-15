use Getopt::Long;
use Getopt::Std;
use Config::General;
use FindBin qw($Bin $Script);
use newPerlBase;
use Cwd qw(abs_path getcwd);
my ($input, $output);
GetOptions(
        "h|?"           =>\&USAGE,
	"gtf:s"		=>\$gtf,
	"fa:s"     	=>\$fa,
	"fpkm:s"	=>\$fpkm,
	"o:s"		=>\$lnc_fa,
)or &USAGE;
&USAGE unless ($gtf and $fa and $fpkm and $lnc_fa);

my %rnas=();############
open(FPKM,$fpkm)||die $!;
while(<FPKM>){
	chomp;next if($_=~/^#/);
	my @tmp=split(/\s+/,$_);
	$rnas{$tmp[0]}++;
}
close(FPKM);
print "perl $Bin/get_lncRNA_fa.pl $gtf $fa $lnc_fa.tmp\n";
`perl $Bin/get_lncRNA_fa.pl $gtf $fa $lnc_fa.tmp`;

$/=">";
my %ref=();
open(FA,"$lnc_fa.tmp")||die $!;
open(OUT,">$lnc_fa")||die $!;
while(<FA>){
	chop;
	next if($_=~/^$/);	
	my ($id,$seq)=split(/\n/,$_,2);
	$seq=~s/\n//g;
	$id=(split(/\s+/,$id))[0];
	print OUT ">$id\n$seq\n"	if(exists $rnas{$id});
}
close(FA);


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
	-gtf	<file>	candicate lncRNA gtf file(Some was filtered)
	-fpkm	<file>	candicate lncRNA fpkm
	-fa	<file>	genome fa file
	-o	<file>	output candicate lncRNA.fa

	-h	Help

Example:

USAGE
	print $usage;
	exit;
}


