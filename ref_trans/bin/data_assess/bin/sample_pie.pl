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
#
#######################################################################################
# ------------------------------------------------------------------
# GetOptions
# ------------------------------------------------------------------
my ($fIn,$fOut);

GetOptions(
				"help|?" =>\&USAGE,
				"od:s"=>\$fOut,
				"i:s"=>\$fIn,
				) or &USAGE;
&USAGE unless ($fIn and $fOut);

# ------------------------------------------------------------------
open IN,$fIn or die "cannot open file $fIn, $!\n";
while (<IN>) {
	chomp;
	next if (/^$/||/\#/);
	my @lines = split /\t/,$_;
	open OUT, ">$fOut/$lines[1].data.stat" or die $! ;
	print OUT "#Type\tRatio\nAdapter Related\t$lines[8]\nLow Quality\t$lines[10]\nClean Data\t",100-$lines[8]-$lines[10],"\n";
	close OUT;
	`perl $Bin/Just_Pie.pl -i $fOut/$lines[1].data.stat -o $fOut/Raw_data_ratio_$lines[1].svg -w 500 -css $Bin/pie12.css `;
	`rm $fOut/$lines[1].data.stat `;
}  
close IN;


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

################################################################################################################
sub cut_str {
	my $string = shift;
	my @str = split /\s+/,$string;
	if (@str > 2) {
		return "$str[0] $str[1]"
	}else{
		return $string;
	}
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
     Contact:	Wang Yajing <wangyj\@biomarker.com.cn> 
Program Date:	2015.07.22
      Modify:	
 Description:	This program is used to ......
       Usage:
		Options:
		-i <file>	AllSample.data.stat

		-od <file>	output

		-h		help

USAGE
	print $usage;
	exit;
}
