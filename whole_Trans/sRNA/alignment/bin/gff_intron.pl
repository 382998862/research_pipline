#!/usr/bin/perl -w
use strict;
use warnings;
use Getopt::Long;
use File::Basename qw(basename dirname);
my ($i,$o);
GetOptions(
				"help|?" =>\&USAGE,
				"in:s"=>\$i,
				"out:s"=>\$o,
				) or &USAGE;
&USAGE unless ($i and $o);

$i=ABSOLUTE_DIR($i);
#$o=ABSOLUTE_DIR($o);
open(IN,"$i");
open(OU,">$o");
#chomp(my @inp=<IN>);
my @inp;
while (<IN>) {
	chomp;
	next if /^\s*$/;
	next if /^#/;
	my @line = split "\t";
	next unless ($line[2] eq 'exon' || $line[2] eq 'gene' || $line[2] eq 'mRNA' || $line[2] eq 'transcript' );
	push @inp,$_;
}

for (my $i=0;$i<=$#inp-1;$i++){
	next if($inp[$i]=~/^\s+/);
	next if($inp[$i]=~/^#/);
	my $lar=$inp[$i];
	my $sma=$inp[$i+1];
	#print "$lar\n$sma\n";
	my @lar=split(/\s+/,$lar);
	my @sma=split(/\s+/,$sma);
	
	if($lar[2]=~/exon/ and $sma[2]=~/exon/){
		my $start;my $end;
		if(abs($lar[4]-$sma[3])<=2 || abs($lar[3]-$sma[4])<=2){
			next;
		}
		elsif($lar[4]<$sma[3]){
			$start=$lar[4]+1;
			$end=$sma[3]-1;
		}
		elsif($lar[3]>$sma[4]){
			$start=$sma[4]+1;
			$end=$lar[3]-1;
		}
		else{
			print "wrong gff3 file：$lar\n";
			die $!;
		}
		my $line;
		$lar[8]=~s/exon/intron/ig;
		$line="$lar[0]\t$lar[1]\tintron\t$start\t$end\t$lar[5]\t$lar[6]\t$lar[7]\t$lar[8]\n";
		print OU "$line";
	}	
	
}
close IN;close OU;


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

sub USAGE {#
	my $usage=<<"USAGE";
ProgramName:
produce intron gff file based on  gff3 file
Contact:	Liuxiaoshuang <liuxs\@biomarker.com.cn> 
Program Date:   2015/8/3
Usage:
  Options:
  -in <file>   input file
  
  
  -out  <file>   output file
  
  -h         Help

USAGE
	print $usage;
	exit;
}
