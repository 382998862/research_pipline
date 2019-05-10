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
my ($id);
GetOptions(
				"help|?" =>\&USAGE,
				"id:s"=>\$id,
				) or &USAGE;
&USAGE unless ($id);


if (defined $id) {
	my @DEG_dir=glob "$id/*";
	my %Stat;
	foreach my $dir (@DEG_dir) {
		if (-d $dir) {
			$dir=~m/.*\/(\S+)/;
			my $nam=$1;
			if ($dir=~/_vs_/){
				my $file=(glob "$dir/*.DEG.final.xls")[0];
				unless(defined $file){
					$Stat{$nam}{up}=0;
					$Stat{$nam}{down}=0;
					$Stat{$nam}{total}=0;
					next;
				}
				open (IN,$file) or die $!;
				while (<IN>) {
					chomp;
					next if /^\#/;
					my $type=(split/\s+/,$_)[-1];
					$Stat{$nam}{up}++ if $type eq 'up';
					$Stat{$nam}{down}++ if $type eq 'down';
					$Stat{$nam}{total}++;
				}
				close IN;
			}
		}
	}
	print "\n#DEG_Set\tAll_DEG\tup-regulated\tdown-regulated\n";
	foreach my $key (sort keys %Stat) {
		$Stat{$key}{up}||=0;
		$Stat{$key}{down}||=0;
		print "$key\t$Stat{$key}{total}\t$Stat{$key}{up}\t$Stat{$key}{down}\n";
	}
	print "\n";
}


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
Usage:
  Options:
  -id <dir>  input dir,forced 
  
 
  -h         Help

USAGE
	print $usage;
	exit;
}
