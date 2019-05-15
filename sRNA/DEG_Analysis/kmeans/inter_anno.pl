#!/usr/bin/perl -w
use strict;
use warnings;
use Getopt::Long;
use Data::Dumper;
use FindBin qw($Bin $Script);
use File::Basename qw(basename dirname);
use Cwd qw(abs_path getcwd);
my $BEGIN_TIME=time();
my $version="1.0.0";
#######################################################################################
#
#######################################################################################
# ------------------------------------------------------------------
# GetOptions
# ------------------------------------------------------------------
my ($anno,$cluster,$out);

GetOptions(
				"help|?" =>\&USAGE,
				"out|o:s"=>\$out,
				"anno|i:s"=>\$anno,
				"cluster|clu:s"=>\$cluster,
				) or &USAGE;
&USAGE unless ($anno and $cluster and $out);
# ------------------------------------------------------------------
# 
# ------------------------------------------------------------------

$anno = &ABSOLUTE_DIR($anno);
$cluster = &ABSOLUTE_DIR($cluster);
$out = &ABSOLUTE_DIR($out);

open (IN,"$anno") or die $!;
open (IN1,"$cluster") or die $!;
open (OUT,">$out") or die $!;

my %anno;my $head;my $str;
while (<IN>) {
	chomp;
	my ($id,$info) = split /\t/,$_,2;
	if (/\#/){
		$head = $info;
		my @tmp = split /\t/,$head;
		$str = "--" x scalar(@tmp);
		next;
	}
	$anno{$id}=$info;
}

while(<IN1>){
	chomp;
	if(/^id/){
		print OUT "$_\t$head\n";next;
	}
	my ($id,$other)=split /\t/,$_,2;
	if(exists $anno{$id}){
		print OUT "$id\t$other\t$anno{$id}\n";
	}else{
		print OUT "$id\t$other\t$str\n";
	}
}
close IN;
close IN1;
close OUT;

#######################################################################################
print STDOUT "\nDone. Total elapsed time : ",time()-$BEGIN_TIME,"s\n";
#######################################################################################

# ------------------------------------------------------------------
# sub function
# ------------------------------------------------------------------
################################################################################################################

sub ABSOLUTE_DIR{ #$pavfile=&ABSOLUTE_DIR($pavfile);
	my $cur_dir=`pwd`;chomp($cur_dir);
	my ($in)=@_;
	my $return="";
	if(-f $in){
		my $dir=dirname($in);
		my $file=basename($in);
		chdir $dir;$dir=`pwd`;chomp $dir;
		$return="$dir/$file";
	}elsif(-d $in){
		chdir $in;$return=`pwd`;chomp $return;
	}else{
		warn "Warning just for file and dir\n";
		exit;
	}
	chdir $cur_dir;
	return $return;
}

################################################################################################################

sub GetTime {
	my ($sec, $min, $hour, $day, $mon, $year, $wday, $yday, $isdst)=localtime(time());
	return sprintf("%4d-%02d-%02d %02d:%02d:%02d", $year+1900, $mon+1, $day, $hour, $min, $sec);
}

################################################################################################################
sub USAGE {
	my $usage=<<"USAGE";
 ProgramName:
     Version:	$version
     Contact:	nie <niepy\@biomarker.com.cn> 
Program Date:	2019
      Modify:	
 Description:	This program is used to ......
       Usage:
		Options:
		-anno		<file>	Integrated_Function.annotation.xls
		-cluster	<file>	kmeans_cluster.txt
		-out		<file>	kmeans_cluster_anno.xls

		-h		help

USAGE
	print $usage;
	exit;
}
