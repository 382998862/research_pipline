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
my ($in,$out);

GetOptions(
				"help|?" =>\&USAGE,
				"out|o:s"=>\$out,
				"in|i:s"=>\$in,
				) or &USAGE;
&USAGE unless ($in and $out);
# ------------------------------------------------------------------
# 
# ------------------------------------------------------------------

$in = abs_path($in);
$out = abs_path($out);

open (IN,$in) or die $!;
open (OUT,">$out") or die $!;

my %hash;
while (<IN>) {
	chomp;
	next if (/\#|^id/) ;
	my @line = split /\s+/,$_;
	if(defined $hash{$line[1]}){
		$hash{$line[1]} .= ",$line[0]";
	}else{
		$hash{$line[1]}=$line[0];
	}
}

print OUT "#cluster_id\tgene_id\n";
foreach my $cluster (sort {$a<=>$b} keys %hash){
	print OUT "$cluster\t$hash{$cluster}\n";
}

close IN;
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
		-in <file>	input file,xxx format,forced

		-out <file>	output file,optional

		-h		help

USAGE
	print $usage;
	exit;
}
