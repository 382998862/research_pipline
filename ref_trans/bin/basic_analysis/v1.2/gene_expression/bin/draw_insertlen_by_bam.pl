#!/usr/bin/perl -w
use strict;
use warnings;
use Getopt::Long;
use Data::Dumper;
use FindBin qw($Bin $Script);
use File::Basename qw(basename dirname);
use newPerlBase;
my %config=%{readconf("$Bin/../../../../../config/db_file.cfg")}; 

my $BEGIN_TIME=time();
my $version="1.0.0";
#######################################################################################

# ------------------------------------------------------------------
# GetOptions
# ------------------------------------------------------------------
my ($fIn,$index,$od);
GetOptions(
				"help|?" =>\&USAGE,
				"i:s"=>\$fIn,
				"index:s"=>\$index,
				"od:s"=>\$od,
				) or &USAGE;
&USAGE unless ($fIn and $index and $od);

mkdir $od unless -d $od;
$od=&ABSOLUTE_DIR($od);

my %hash;
open (IN,"samtools view $fIn|") or die $!;
while (<IN>) {
	my $insert=(split/\t/,$_)[8];
	next if $insert<=0;
	$hash{$insert}++;
}
close (IN) ;

	open (OUT,">","$od/$index.insertSize") or die $!;
	foreach my $insertSize (sort {$a <=> $b} keys %hash) {
		print OUT $insertSize,"\t",$hash{$insertSize},"\n";
	}
	close (OUT) ;

	my $max=(sort {$b <=> $a}values %hash)[0];
	open (OUT ,">","$od/$index.insertSize.psvg") or die $!;
	print OUT "Type:Line","\n";
	print OUT "Width:600","\n";
	print OUT "Height:400","\n";
	print OUT "WholeScale:0.9","\n";
	print OUT "FontSize:25","\n";
	print OUT "X:InsertSize","\n";
	print OUT "Y:Number of Reads","\n";
	print OUT "XStep:100","\n";
	print OUT "YStep:",int($max/10),"\n";
	print OUT "XStart:0","\n";
	print OUT "YStart:0","\n";
	print OUT "XEnd:800","\n";
	print OUT "YEnd:",$max+200,"\n\n";
	print OUT "Color:red","\n";

	foreach my $key (sort {$a <=> $b} keys %hash) {
		print OUT $key,":",$hash{$key},"\n";
	}
	close (OUT);

	chdir $od;
	system ("$config{distributing_svg1} $index.insertSize.psvg $index.insertSize.svg");
	system ("$config{svg2xxx1} $index.insertSize.svg");


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
Usage:
  Options:
  -i        <file>         input file,bam format,forced 
  
  -index    <str>          index of outfile,forced 
  
  -od       <dir>          output dir,forced 
  
  -h                       Help

USAGE
	print $usage;
	exit;
}
