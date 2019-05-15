#!/usr/bin/perl -w
use strict;
use warnings;
use Getopt::Long;
use Data::Dumper;
use FindBin qw($Bin $Script);
use File::Basename qw(basename dirname);
use newPerlBase;
my $BEGIN_TIME=time();
my $version="1.0.1";
#######################################################################################

# ------------------------------------------------------------------
# GetOptions
# ------------------------------------------------------------------
my ($fIn,$collapsed,$od,$prefix);
GetOptions(
				"help|?" =>\&USAGE,
				"i:s"=>\$fIn,
				"c:s"=>\$collapsed,
				"od:s"=>\$od,
				"prefix:s"=>\$prefix,
				) or &USAGE;
&USAGE unless ($fIn);

my %CFG=%{readconf("$Bin/../CFG")};

if (!defined $collapsed){
	$collapsed=0;
} 
mkdir "$od" unless -d "$od";
open (OUT,">$od/$prefix.stat") or die $!;
print OUT "#length\treads\n";
my %H;
$/=">";
open (IN,$fIn) or die $!;
<IN>;
while (<IN>) {
	chomp;
	my ($head,$info)=split/\n/,$_;
	my $len=length $info;
	if ($collapsed==1) {
		my $deep=(split/_x/,$head)[1];
		$H{$len}+=$deep;
	}
	else{
		$H{$len}+=1;}
	
}
close (IN) ;

foreach my $key (sort {$a<=>$b} keys %H) {
	print  OUT "$key\t$H{$key}\n";
}

close OUT;


runOrDie("$CFG{Rscript} $Bin/draw_length_distribution.r $od/$prefix.stat  $od/$prefix.length ");


#######################################################################################
#print STDOUT "\nDone. Total elapsed time : ",time()-$BEGIN_TIME,"s\n";
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
Program Date:   2012.8.23
Modify:         liuhui <liuhui860307\@126.com>  2010-7-27  ##move adapter recongnize programs to my bin dir
Description:	this program is used to ......
Usage:
  Options:
  -i    <file>  input file,fasta format,'seq_0000009_x43832',forced 
  -c    <num> 	file format ,id  coolapsed  ,default = 1
  -od   <dir>	file dir,
  -prefix <str> prefix of out file,
  -h         Help

USAGE
	print $usage;
	exit;
}
