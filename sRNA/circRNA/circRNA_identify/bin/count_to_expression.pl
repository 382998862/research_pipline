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
my ($fIn,$o);
GetOptions(
				"help|?" =>\&USAGE,
				"i:s"=>\$fIn,
				"o:s"=>\$o,
				) or &USAGE;
&USAGE unless ($fIn and $o);


my %H;
open (IN,$fIn) or die $!;
my $head_line=<IN>;
my @Samples=split/\t/,$head_line;
shift @Samples;
pop @Samples;
while (<IN>) {
	my @A=split/\t/,$_;
	shift @A;
	pop @A;
	for (my $i=0;$i<@A ;$i++) {
		$H{$Samples[$i]}+=$A[$i];
	}
}
close (IN) ;

open (IN,$fIn) or die $!;
open (OUT,">$o/All_gene_fpkm.list") or die $!;
open EXP,">$o/All_gene_expression.list" or die $!;
my $line=join "\t",@Samples;
my @sam;
foreach my $s (@Samples) {
	push @sam,$s."_count";
	push @sam,$s."_FPKM";
}
my $exp_line=join "\t",@sam;
print OUT "#ID\t$line\n";
print EXP "\#ID\t$exp_line\n";
<IN>;
while (<IN>) {
	my @A=split/\t/,$_;
	my $name=shift @A;
	my $len=pop @A;
	my @B;
	for (my $i=0;$i<@A ;$i++) {
		$B[$i]="$A[$i]\t".$A[$i]/($len/1000)/($H{$Samples[$i]}/1000000);
		$A[$i]=$A[$i]/($len/1000)/($H{$Samples[$i]}/1000000);
	}
	$line=join "\t",@A;
	$exp_line=join "\t",@B;
	print OUT "$name\t$line\n";
	print EXP "$name\t$exp_line\n";
}
close IN;
close OUT;
close EXP;


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


sub USAGE {#
	my $usage=<<"USAGE";
ProgramName:
Version:	$version
Contact:	Zhang XueChuan <zhangxc\@biomarker.com.cn> 
Program Date:   2013.10.11
Usage:
  Options:
  -i    <file>  input file,gene count file,forced 
  
  -o    <file>  out direction for gene expression file,forced 
  
  -h         Help

USAGE
	print $usage;
	exit;
}
