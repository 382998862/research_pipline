#!/usr/bin/perl -w
use strict;
use warnings;
use Getopt::Long;
use Data::Dumper;
use FindBin qw($Bin $Script);
use File::Basename qw(basename dirname);
my $BEGIN_TIME=time();
my $version="1.0.0";
#######################################################################################

# ------------------------------------------------------------------
# GetOptions
# ------------------------------------------------------------------
my ($fIn,$fIn1,$fOut);
GetOptions(
				"help|?" =>\&USAGE,
				"o:s"=>\$fOut,
				"i:s"=>\$fIn,
				"gtf:s"=>\$fIn1,
				) or &USAGE;
&USAGE unless ($fIn and $fIn1 and $fOut );
open (IN,$fIn) or die $!;
my %hgeneid = () ;
while (<IN>){
	chomp ;
	next if(/^#/);
	# get id
	my $id = (split /\s+/,$_)[0] ;
	$hgeneid{$id} = 1 ;
}
close(IN) ;

open (GTF, $fIn1) || die $! ;
open (OUT, ">$fOut") || die $! ;
$/ = "\n" ;
while (<GTF>){
	chomp ;
	next if (/\#/||/^$/) ;
	$_=~/transcript_id\s\"([^\"]+)\";/;
	my $id =$1;
	$id =~s/"//g;
	$id =~s/\;//g;
	if (defined $hgeneid{$id}){
		print OUT "$_ \n" ;
	}
}
$/ = "\n" ;
close(GTF) ;
close(OUT) ;


#######################################################################################
print STDOUT "\nDone. Total elapsed time : ",time()-$BEGIN_TIME,"s\n";
#######################################################################################

# ------------------------------------------------------------------
# sub function
# ------------------------------------------------------------------
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

sub GetTime {
	my ($sec, $min, $hour, $day, $mon, $year, $wday, $yday, $isdst)=localtime(time());
	return sprintf("%4d-%02d-%02d %02d:%02d:%02d", $year+1900, $mon+1, $day, $hour, $min, $sec);
}


sub USAGE {#
	my $usage=<<"USAGE";
Program:
Version: $version
Contact:Deng Dejing <zhoum\@biomarker.com.cn> 
Description:Abstract_fasta_seq_by_id
Usage:
  Options:
  -i  <file>      id-list file          must be given
  -gtf<file>      gtf file            must be given
  -o  <file>      outfile               must be given
 
  -h         Help

USAGE
	print $usage;
	exit;
}
