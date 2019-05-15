#!/usr/bin/perl -w
use strict;
use Getopt::Long;
use Data::Dumper;
use FindBin qw($Bin $Script);
use File::Basename qw(basename dirname);
use Cwd qw(abs_path);
use newPerlBase;
my $BEGIN_TIME=time();
my $version="1.0";

my($idir,$sample,$genome,$od);
GetOptions(
	"help|?"		=>\&USAGE,
	"idir:s"		=>\$idir,
	"sample:s"		=>\$sample,
#	"genome:s"		=>\$genome,
	"od:s"			=>\$od,
) or &USAGE;
&USAGE unless ($idir and $sample);

my %CFG=%{readconf("$Bin/../CFG")};

$idir=abs_path($idir);
$od ||="$idir/isomiR";
mkdir $od unless(-d $od);
$od=abs_path($od);

#my $genome_name=basename ($genome);
my @fas=glob("$idir/sRNA_Alignment/Ref_Database/*.fa");
$genome=$fas[0]	if($fas[0] !~ /mature\.fa/);
$genome=$fas[1] if($fas[1] !~ /mature\.fa/);

my @samples=split /,/,$sample;
my @sort_samples = sort @samples;

chdir $od;
&log_current_time("Change the bases begin:");
&cmd_call("perl $Bin/bin/upper.pl -i $idir/sRNA_Analysis/miRDeep2/miRNA_Quantify/All_miRNA.expressed.fa -o $od/All_miRNA.expressed.change.fa");
&cmd_call("perl $Bin/bin/upper.pl -i $idir/sRNA_Analysis/miRDeep2/miRNA_Quantify/All_miRNA_Pre.expressed.fa -o $od/All_miRNA_Pre.expressed.change.fa");
&log_current_time("Change the bases done!");

open(OUT,">$od/isomiR.cfg") or die $!;
print OUT "bowtie_path: $CFG{bowtie} \n";
print OUT "bowtie-build_path: $CFG{'bowtie_build'} \n";
for(my $i=0;$i<@sort_samples;$i++){
	print OUT "lib: $idir/sRNA_Alignment/$sort_samples[$i]/$sort_samples[$i].clean.fa $sort_samples[$i] fa \n";
}
print OUT "main_ref: $od/All_miRNA_Pre.expressed.change.fa yes \n";
print OUT "filter_ref: $genome  no \n";
print OUT "known_miRNAs:  $od/All_miRNA.expressed.change.fa \n";
print OUT "M3: yes 1 \n";
print OUT "M5: yes 1 \n";
print OUT "RangeSize: 18 26 \n";
print OUT "cutoff: 20 \n";
close OUT;

stepStart(1,"isomiR");
chdir $od;
if(-e "$genome.1.ebwtl"){
	&cmd_call("$CFG{python} $Bin/bin/isomiR_large.py $od/isomiR.cfg -f");
}elsif(-e "$genome.1.ebwt"){
	&cmd_call("$CFG{python} $Bin/bin/isomiR.py $od/isomiR.cfg -f");
}
stepTime(1);


sub cmd_call {
	print "@_\n";
	system(@_) == 0 or die "system @_ failed: $?";
}


sub log_current_time{
	# get parameter
     my ($info) = @_;

     # get current time with string
     my $curr_time = date_time_format(localtime(time()));

     # print info with time
     print "[$curr_time] $info\n";
}


sub date_time_format{
	my ($sec, $min, $hour, $day, $mon, $year, $wday, $yday, $isdst)=localtime(time());
    return sprintf("%4d-%02d-%02d %02d:%02d:%02d", $year+1900, $mon+1, $day, $hour, $min, $sec);
}


sub USAGE {#
	my $usage=<<"USAGE";
Program: $Script
Version: $version
Usage:
  Options:
	-idir	input directory(contained sRNA_Alignment and sRNA_Analysis dir)

	-sample	sample list

	-od	output directory

	-h	Help

USAGE
	print $usage;
	exit;
}
