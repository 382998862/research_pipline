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
my ($fIn,$fOut,$key);

GetOptions(
				"help|?" =>\&USAGE,
				"od:s"=>\$fOut,
				"i:s"=>\$fIn,
				"key:s"=>\$key,
				) or &USAGE;
&USAGE unless ($fIn and $key and  $fOut);
# ------------------------------------------------------------------
# 
# ------------------------------------------------------------------
my @array;
my @input=split(/,/,$fIn);
my %dict;
foreach my $run (@input) {
	open (IN,$run) or die $!;
	while (<IN>) {
		next if ($_=~/^\#|^\s*$/);
		next if ($_!~/$key\-/) ;
		my @lines=split /\t/,$_;
		my $sample=(split(/-|_/,$lines[0]))[2];
		if (exists $dict{$sample}) {
			my @result=@{$dict{$sample}};
			my $total_read=$lines[1]+$result[1];
			my $adapter=sprintf("%.2f",($lines[1]*$lines[2]+$result[1]*$result[2])/($lines[1]+$result[1]));
			my $rRNA=sprintf("%.2f",($lines[1]*$lines[3]+$result[1]*$result[3])/($lines[1]+$result[1]));
			my $inferior=sprintf("%.2f",($lines[1]*$lines[4]+$result[1]*$result[4])/($lines[1]+$result[1]));
			@array=($sample,$total_read,$adapter,$rRNA,$inferior);
			$dict{$sample}=[@array];

		} else {
			@array=($sample,$lines[1],$lines[2],$lines[3],$lines[4]);
			$dict{$sample}=[@array];
		}

	}
}


open (OUT,">$fOut/AllSample.data.stat") or die $!;
close IN;
print OUT "#Sample_ID\tRaw_reads\tAdapter%\trRNA%\tInferior\n";
foreach my $sample (sort keys %dict) {
	print OUT join "\t",@{$dict{$sample}};
	print OUT "\n";

}
#print OUT Dumper %dict;
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

sub max{#&max(lists or arry);
	#求列表中的最大值
	my $max=shift;
	my $temp;
	while (@_) {
		$temp=shift;
		$max=$max>$temp?$max:$temp;
	}
	return $max;
}

################################################################################################################

sub min{#&min(lists or arry);
	#求列表中的最小值
	my $min=shift;
	my $temp;
	while (@_) {
		$temp=shift;
		$min=$min<$temp?$min:$temp;
	}
	return $min;
}

################################################################################################################

sub revcom(){#&revcom($ref_seq);
	#获取字符串序列的反向互补序列，以字符串形式返回。ATTCCC->GGGAAT
	my $seq=shift;
	$seq=~tr/ATCGatcg/TAGCtagc/;
	$seq=reverse $seq;
	return uc $seq;			  
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
     Contact:	Simon Young <yangxh\@biomarker.com.cn> 
Program Date:	2012.07.02
      Modify:	Wang Yajing <wangyj\@biomarker.com.cn>
 Description:	This program is used to statics
       Usage:
		Options:this scripts only can be used when all samples of the project sequeced in the same RUN;
		-i <file>	input file,
					only a RUN:  All.sample.Total.xls
					multi  RUN:  All.sample.Total.xls,All.sample.Total.xls ,forced
		-od <dir>	output dirname,'./',optional
		
		-key <str>  ID of BMK Contract,      forced,(e.g. R96,T60)
		-h		help

USAGE
	print $usage;
	exit;
}
