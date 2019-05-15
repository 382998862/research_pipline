#!/usr/bin/perl -w
use strict;
use Getopt::Long;
use Cwd qw(abs_path);
use FindBin qw($Bin $Script);
use File::Basename qw(basename dirname);
my $version="1.0";

my ($indir);
GetOptions(
	"help|?" =>\&USAGE,
	"indir:s"=>\$indir,
) or &USAGE;
&USAGE unless ($indir);

$indir=abs_path($indir);

my $cmd;
my $name_list   = "$indir/Result/chname.list"; # yaling at 17-12-14
my $known_mrd	= (glob "$indir/miRNA_Quantify/expression_analyses/*/miRBase.mrd")[0];
my $csv			= (glob "$indir/miRNA_Quantify/miRNAs_expressed_all_samples_*.csv")[0];
my $i			= (glob "$indir/miRNA_Quantify/expression_analyses/expression_analyses_*/All_miRNA_mapped.arf")[0];
my $name=basename($csv);
$name=~/miRNAs_expressed_all_samples_(\S+).csv/;
my $time=$1;
print "$time\n";

my $knownmiR=(glob "$indir/miRNA_Quantify/Known_mature_expressed.fa")[0];
my $line=`less $knownmiR|wc -l`;
chomp $line;
if($line!=0){

	my %KNOWN;
	open(IN,$csv) or die $!;
	while(<IN>){
		chomp;
		next if (/^#/ || /novel/);
		my ($id)=(split /\s+/,$_)[2];
		$KNOWN{$id}=1;
	}
	close IN;

	open(KN,">$indir/miRNA_Quantify/known.mrd") or die $!;
	open(MRD,"$known_mrd") or die $!;
	$/=">";
	<MRD>;
	while(<MRD>){
		chomp;
		my ($mi)=(split /\s+/,$_)[0];
		if(exists $KNOWN{$mi}){
			print KN ">$_\n";
		}
	}
	$/="\n";
	close MRD;
	close KN;


	$cmd="cd $indir/miRNA_Quantify &&";
	$cmd .="perl $Bin/make_html2.pl -q $indir/miRNA_Quantify/known.mrd -k $indir/miRNA_Quantify/All_miRNA.expressed.fa -i $i -M $csv  -o  -l -y $time \n";
	print "$cmd\n";
	system ($cmd);
}

my %name;
open(LIST, $name_list) || die $!;
while (<LIST>) {
	next if /^\s+$/;
	my @a = split;
	$name{$a[0]} = $a[1];
}
close LIST;

open (R, "$indir/Result/novel_conservative.mrd") || die $!;
open (W, ">$indir/Result/novel_conservative.mrd.new") || die $!;
while (<R>) {
	s/^>(\S+)/>$name{$1}/ if /^>/;
	print W $_;
}
close R;
close W;
open (R, "$indir/Result/novel_unconservative.mrd") || die $!;
open (W, ">$indir/Result/novel_unconservative.mrd.new") || die $!;
while (<R>) {
	s/^>(\S+)/>$name{$1}/ if /^>/;
	print W $_;
}
close R;
close W;

`cat $indir/Result/novel_conservative.mrd.new $indir/Result/novel_unconservative.mrd.new >$indir/Result/novel.mrd.new`;

$cmd="cd $indir/miRNA_Quantify &&";
$cmd .="perl $Bin/make_html.pl -f $indir/Result/novel.mrd.new -s $csv -c -e -v -10 -y $time \n";

print "$cmd\n";
system ($cmd);

#`cp -r $indir/miRDeep2/pdf_$time/* $indir/Diff_Analysis/pdf_$time`;






sub USAGE {#
	my $usage=<<"USAGE";
Program:	$Script
Version:	$version
Usage:
  Options:
	-indir	input directory
	-h	Help
USAGE
	print $usage;
	exit;
}


