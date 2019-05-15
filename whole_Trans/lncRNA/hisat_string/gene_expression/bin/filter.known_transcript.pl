#!/usr/bin/perl -w
use strict;

use Getopt::Long;
use Data::Dumper;
use FindBin        qw($Bin $Script);
use File::Basename qw(basename dirname);
use Cwd            'abs_path';

my $BEGIN_TIME=time();

# ------------------------------------------------------------------
my ($gff, $known_trans_fa, $new_trans_fa , $out);
GetOptions(
				"help|?"      =>\&USAGE,
				"gff:s"	      =>\$gff,
				"kfa:s"       =>\$known_trans_fa,
				"nfa:s"       =>\$new_trans_fa,
				"out:s"     =>\$out,
				) or &USAGE;
&USAGE unless ($gff and $known_trans_fa and $new_trans_fa and $out) ;
################################

$gff = abs_path($gff);
my %hash;
open IN,"$gff" || die $!;
while (<IN>) {
	chomp;
	next if (/^$/ || /^\#/);
	my @tmp = split/\t+/,$_;
	my @info = split/;/,$tmp[8];
	if($tmp[2] eq "gene"){
		if($info[0] =~ /ID=(.*)/){
			$hash{$1}=1;
			print  "$1\n";
		}
	}
	
}

my %kh;
$/=">";
open(IN1,"$known_trans_fa") or die $!;
<IN1>;
open(OUT,">$out") or die $!;
while(<IN1>){
	chomp;
	my ($name,$fa)=split /\n/,$_,2;
	my ($id,$trans) = split /\s+/,$name,2;
	$kh{$id}=$fa;
	if(exists $hash{$id}){
		print OUT ">$id $trans\n$kh{$id}";
	}
}
$/="\n";

open(IN2,"$new_trans_fa") or die $!;
while(<IN2>){
	chomp;
	print OUT "$_\n";
}

close OUT;
close IN1;
close IN;



sub USAGE {#
	my $usage=<<"USAGE";

Usage:
	-gff
	-kfa
	-nfa
	-out               Human.newGene.track.list.info     must be given  ;
	
USAGE
	print $usage;
	exit;
}
