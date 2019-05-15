#!/usr/bin/perl
use strict;
use warnings;
use Getopt::Long;
use Data::Dumper;
use FindBin qw($Bin $Script);
use File::Basename qw(basename dirname);
use Cwd qw(abs_path getcwd);
my $BEGIN_TIME=time();
my $version="1.0.0";
my ($o,$order,$input);
GetOptions(
			"help|?"=>\&USAGE,
			"i:s"=>\$input,
			"o:s"=>\$o,
			"order:s"=>\$order,

			)or &USAGE;
&USAGE unless ($order and $o and $input);

my %sam_combine;
my $num=0;
my %group;
my @order=();
if ($order=~/;/){
	my @temp=split(/;/,$order);
	foreach my $group (@temp) {
		if($group=~/,/){
			my @temp2=split(/,/,$group);
			push @order,join("_",@temp2);
			foreach my $sam (@temp2) {
				$sam_combine{$sam}=join("_",@temp2);
				$group{$sam_combine{$sam}}=1;
			}
		}
		else {
				push @order,$group;
				$sam_combine{$group}=$group;
				$group{$sam_combine{$group}}=1;
		}
	}
}else{
	my @temp2=split(/,/,$order);
		foreach my $sam (@temp2) {
			push @order,$sam;
			$sam_combine{$sam}=$sam;
			$group{$sam_combine{$sam}}=1;
	}
}


my %change;

open (IN,$input)||die $!;
open (OUT,">$o") ||die $!;
while (<IN>) {
	chomp;
	next if ($_=~/^\s*$/);
	if ($_=~/\#/) {
		my @head=split(/\t/,$_);
		$num=0;
		print OUT "#ID\t",join("\t",@order),"\n";
		shift@head;
		foreach my $sam (@head) {
			if (exists $sam_combine{$sam}) {
				$change{$num}=$sam_combine{$sam};
			}
			$num++;
		}
		next;
	}
	my @temp3=split(/\t/,$_);
	my $gene=shift@temp3;
    my %fpkm=();
	for (my $i=0;$i<=$#temp3 ;$i++) {
		if (exists $change{$i}) {
			my $values=$change{$i};
			if (exists $fpkm{$values}) {
				$fpkm{$values}=[$fpkm{$values}->[0]+1,$fpkm{$values}->[1]+$temp3[$i]];
			}
			else {
				$fpkm{$values}=[1,$temp3[$i]];
			}
		}
	}
	print OUT $gene;
	foreach my $group (@order) {
		#print "$group\n";
		print OUT "\t",$fpkm{$group}->[1]/$fpkm{$group}->[0]; ##mean value for each sample 
	}
	print OUT "\n";
	
}
$/="\n";
close(IN);
close(OUT);
##################################################


################################################################################
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
ProgramName:  blast_supply
Version:	$version
Contact:	Wang Yajing <wangyj\@biomarker.com.cn> 
Program Date:   2015/08/17
Usage:
Options:
-i      <file>  All_gene_fpkm.list ,required 
-o      <file>  output file,required 
-order  <str>   sample name  used  and order 
           example: T01,T02;T03,T04;T05,T06 or T01,T02,T03
-h      Help

USAGE
	print $usage;
	exit;
}
