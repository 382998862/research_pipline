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
my ($o);
my @sample;
GetOptions(
				"help|?" =>\&USAGE,
				"i:s{,}"=>\@sample,
				"o:s"=>\$o,
				) or &USAGE;
&USAGE unless ($o and @sample);


my %H;
foreach my $sam (@sample) {
	my $name=basename $sam;
	$name=~/^([^\.]+)\./;
	$name=$1;
	open (IN,$sam) or die $!;
	while (<IN>) {
		chomp;
		next if /^\s*$/;
		next if /^\#/;
		my @A=split/\t/,$_;
		$H{$name}{$A[0]}=$A[1];
	}
	close IN;
}

open (OUT,">All_Type.r") or die $!;
print OUT "#!/share/nas2/genome/biosoft/R/3.1.1/lib64/R/bin/Rscript\n";
print OUT "library(ggplot2)\n";
my @S;
my @T;
my @P;
foreach my $key (sort keys %H) {
	foreach my $t (sort keys %{$H{$key}}) {
		push @S,$key;
		push @T,$t;
		push @P,$H{$key}{$t};
	}
}
my $S=join '","',@S;
my $T=join '","',@T;
my $P=join ',',@P;
print OUT 'Sample <- c("'.$S.'")'."\n";
print OUT 'Type <- c("'.$T.'")'."\n";
print OUT 'Percent <- c('.$P.")\n";
print OUT "type <- data.frame(Sample,Type,Percent)\n";
print OUT "p<-ggplot(data=type,aes(x=Sample,y=Percent,fill=Type))+ geom_bar(position='fill',stat=".'"identity"'.",alpha=0.6)\n";
print OUT 'png(filename="'.$o.'", height = 900, width = 1200, res = 150, units = "px")'."\n";
print OUT "print(p)\n";
print OUT "dev.off()\n";

`$config{Rscript} All_Type.r`;
`rm All_Type.r`;




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
  -i <str>   input files,T1.type.stat T3.type.stat T2.type.stat, forced
  
  -o <file>  output file,png format,forced 
  
  -h         Help

USAGE
	print $usage;
	exit;
}
