#!/usr/bin/perl -w
use strict;
use warnings;
use Getopt::Long;
use Data::Dumper;
use threads;
use Cwd qw(abs_path);
use FindBin qw($Bin $Script);
use File::Basename qw(basename dirname);
use newPerlBase;
my $Title="LncRNA";  
my $version="v2.0.3"; 
my %config=%{readconf("$Bin/../../lncRNA_pip.cfg")};
my ($in,$odir);
GetOptions(
	        "help|?"=>\&USAGE,
			"o:s"=>\$odir,
			"i:s"=>\$in,
)or &USAGE;
&USAGE unless ($odir and $in );

system "mkdir -p $odir" unless (-d $odir);
$odir = abs_path($odir);
&log_current_time("$Script start...");

#stat class code
my %class;
open IN ,"$in";
while (<IN>){
	chomp;
	next if /^$/ || /^#/;
	my @tmp=split /\t/,$_;
	next if ($tmp[2]!~/transcript/);
	$tmp[8]=~/transcript_id \"(.+?)\".*?class_code \"(.+?)\"/;
#	$ID=~/gene_id \"(.+?)\".*?transcript_id \"(.+?)\"/;
	my $transid=$1;
	my $classcode=$2;
	if (exists $class{$transid}){
		next;
	}else {
		$class{$transid}=$classcode;
	}
}
close IN;
open OUT ,">$odir/lnc_class_code.xls";
for my $id (sort keys %class){
	print OUT "$id\t$class{$id}\n";
}
close OUT;
my %code;
#id list
open IN ,"$odir/lnc_class_code.xls";
while(<IN>){
	chomp;
	next if /^$/ || /^#/;
	my $c=(split /\t/,$_)[1];
	if (exists $code{$c}){
		$code{$c}++;
	}else{
		$code{$c}=1
	}
}
close IN;
#stat
open OUT,">$odir/lnc.plot.class_code.stat";
print OUT "#lncRNA type\tlncRNA number\n";
if (exists $code{"u"}){
	print OUT "lincRNA\t$code{u}\n";
}else {
	print OUT "lincRNA\t0\n"
}
if (exists $code{"x"}){
        print OUT "Antisense-lncRNA\t$code{x}\n";
}else {
	print OUT "Antisense-lncRNA\t0\n"
}
if (exists $code{"i"}){
        print OUT "Intronic-lncRNA\t$code{i}\n";
}else {
	print OUT "Intronic-lncRNA\t0\n";

}
my $num;
if (exists $code{"j"}){
	$num=$num+$code{j};
        #print OUT "novel_isoform_lncRNA\t$code{j}\n";
}
if (exists $code{"e"}){
        $num=$num+$code{e};
}
if (exists $code{"o"}){
        $num=$num+$code{o};
}
if (exists $code{"j"} || exists $code{"e"}  || exists $code{"o"} ){
	print OUT "sense_lncRNA\t$num\n";
	}
####
#"j" => "novel_isoform_lncRNA",
#        "e" => "pre-mRNA_lncRNA",
#                "o" => "exonic_overlap_lncRNA",
#foreach my $key (sort keys %code){
#	print OUT "$key\t$code{$key}\n";
#}
close OUT;
`$config{Rscript} $Bin/simpleBar.r --infile $odir/lnc.plot.class_code.stat --outfile $odir/lncRNA_filter_pie.stat.png --x.col 1 --y.col 2 --x.lab type --y.lab Number `;

#################################
sub log_current_time{
	my ($info) = @_;
	my $curr_time = date_time_format(localtime(time()));
	print "[$curr_time] $info\n";
}
sub date_time_format {
	my ($sec, $min, $hour, $day, $mon, $year, $wday, $yday, $isdst)=localtime(time());
	return sprintf("%4d-%02d-%02d %02d:%02d:%02d", $year+1900, $mon+1, $day, $hour, $min, $sec);
}

###############################
sub USAGE {
	my $usage=<<"USAGE";
------------------------------------------------------------------------------------------------
	Program: $Script
	Contact: renhd\@biomarker.com.cn
	   Date:
	  Usage:


		-i input file , filter_final.gtf 
		-o output diretory
  		-h      help documents
------------------------------------------------------------------------------------------------
USAGE
	print $usage;
	exit();
}
