#!/usr/bin/perl -w
use strict;
use warnings;
use Getopt::Long;
use Data::Dumper;
use FindBin qw($Bin $Script);
use File::Basename qw(basename dirname);
use newPerlBase;
my $BEGIN_TIME=time();
my $version="1.0.0";
#######################################################################################
#
#######################################################################################
# ------------------------------------------------------------------
# GetOptions
# ------------------------------------------------------------------
my ($fIn,$fId,$fOut);

GetOptions(
				"help|?" =>\&USAGE,
				"o:s"=>\$fOut,
				"l:s"=>\$fId,
				"f:s"=>\$fIn,
				) or &USAGE;
&USAGE unless ($fIn and $fId and $fOut);

my @fIn = split /,/,$fIn;
for (my $i=0; $i<@fIn; $i++) {
    if (-f "$fIn[$i]") {
        $fIn[$i] = &ABSOLUTE_DIR($fIn[$i]);
    } else {
        print "ERROR: file $fIn[$i] not exists.\n";
        &USAGE;
    }
}

# ------------------------------------------------------------------
# 
# ------------------------------------------------------------------
my %ids;
open (ID,$fId) or die $!;

while (<ID>) {
	chomp;
	next if (/^#/);
	next if (/^\s+$/);
	my @col=split /\t/;
#	$ids{$col[1]}=1;#默认基因ID在第二列
	$ids{$col[0]}=1;#默认基因ID在第1列
}

close ID;

# ------------------------------------------------------------------
# 
# ------------------------------------------------------------------
open (IN,"cat @fIn|") or die $!;
open (OUT,">$fOut") or die $!;
$/=">";
<IN>;

while (<IN>) {
	chomp;
	next if (/^\s+$/);
	my @lines=split /\n/;
	my ($id,$seq);

	for (my $l=0;$l<@lines;$l++) {
		if ($l==0) {
			$id=(split /\s+/,$lines[$l])[0];
		} else {
			$seq.=$lines[$l]."\n";
		}
	}

	print OUT ">$lines[0]\n$seq" if (exists $ids{$id});
}

close IN;
close OUT;

#######################################################################################
&timeLog("$Script Done. Total elapsed time : ".time()-$BEGIN_TIME."s");
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
      Modify:	
 Description:	This program is used to ......
       Usage:
		Options:
		-l <file>	id list file,tab format,forced

		-f <file>	input file,fasta format,forced

		-o <file>	output file,optional

		-h		help
		Example:
			perl abstractFabyId.pl -l geneId.list -f Vitis_vinifera.cds.fa -o Grape_miRNA_Target.cds.fa

USAGE
	print $usage;
	exit;
}
