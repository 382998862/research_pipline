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
my ($i,$o);
GetOptions(
				"help|?" =>\&USAGE,
				"i:s"=>\$i,
				"o:s"=>\$o,
				) or &USAGE;
&USAGE unless ($i and $o);


my %Info;
open (IN,$i) or die $!;
while (<IN>) {
	chomp;
	next if (/\#/||/^$/) ;
	my ($chr,$type,$start,$end,$strand,$ID)=(split /\t/,$_)[0,2,3,4,6,8];
	next unless $type=~/exon/;
	$ID=~/transcript_id\s\"([^\"]+)\";\sgene_id\s\"([^\"]+)\";/;
	my $gene=$2;
	my $trans=$1;
	$Info{$chr}{$gene}{INFO}{$trans}{INFO}{$start}="$chr\tunknown\texon\t$start\t$end\t.\t$strand\t.\t";
#	$Info{$chr}{$gene}{INFO}{$trans}{INFO}="$chr\tunknown\texon\t$start\t$end\t.\t$strand\t.\t";
	$Info{$chr}{$gene}{STRAND}=$strand;
	$Info{$chr}{$gene}{START}||=1000000000000000;
	$Info{$chr}{$gene}{START}=($Info{$chr}{$gene}{START}>$start)?$start:$Info{$chr}{$gene}{START};
	$Info{$chr}{$gene}{END}||=0;
	$Info{$chr}{$gene}{END}=($Info{$chr}{$gene}{END}>$end)?$Info{$chr}{$gene}{END}:$end;
	$Info{$chr}{$gene}{INFO}{$trans}{START}||=1000000000000000;
	$Info{$chr}{$gene}{INFO}{$trans}{START}=($Info{$chr}{$gene}{INFO}{$trans}{START}>$start)?$start:$Info{$chr}{$gene}{INFO}{$trans}{START};
	$Info{$chr}{$gene}{INFO}{$trans}{END}||=0;
	$Info{$chr}{$gene}{INFO}{$trans}{END}=($Info{$chr}{$gene}{INFO}{$trans}{END}>$end)?$Info{$chr}{$gene}{INFO}{$trans}{END}:$end;
}
close (IN) ;

open (OUT,">$o") or die $!;
foreach my $chr (sort keys %Info) {
	foreach my $gene (sort {$Info{$chr}{$a}{START}<=>$Info{$chr}{$b}{START}} keys %{$Info{$chr}}) {
		print OUT "$chr\tunknown\tgene\t$Info{$chr}{$gene}{START}\t$Info{$chr}{$gene}{END}\t.\t$Info{$chr}{$gene}{STRAND}\t.\tID=$gene\n";
		foreach my $trans (sort {$Info{$chr}{$gene}{INFO}{$a}{START}<=>$Info{$chr}{$gene}{INFO}{$b}{START}} keys %{$Info{$chr}{$gene}{INFO}}) {
			print OUT "$chr\tunknown\tmRNA\t$Info{$chr}{$gene}{INFO}{$trans}{START}\t$Info{$chr}{$gene}{INFO}{$trans}{END}\t.\t$Info{$chr}{$gene}{STRAND}\t.\tID=$trans;Parent=$gene\n";
			my $limit=0;
			foreach my $key (sort {$a<=>$b} keys %{$Info{$chr}{$gene}{INFO}{$trans}{INFO}}) {
				$limit++;
				print OUT $Info{$chr}{$gene}{INFO}{$trans}{INFO}{$key}."ID=exon:$trans:$limit;Parent=$trans\n";
			}
		}
	}
}
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


sub USAGE {#
	my $usage=<<"USAGE";
ProgramName:
Version:	$version
Contact:	Zhang XueChuan <zhangxc\@biomarker.com.cn> 
Usage:
  Options:
  -i <file>  input file,gtf format,forced 
  
  -o <file>  output file,gff format,forced 
  
  -h         Help

USAGE
	print $usage;
	exit;
}
