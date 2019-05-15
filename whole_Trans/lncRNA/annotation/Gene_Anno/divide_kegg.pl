#!/bin/env perl
#use strict;
use autodie;
if(@ARGV!=3){
	&usage();
}
my %class;
my %fileH;
my ($kegg_fa,$org,$outpath)=@ARGV;
open I,"$org";
while(<I>){
	chomp;
	my @tmp=split/\t/;
	my @type=split/;/,$tmp[3];
	$class{$tmp[1]}=$type[1];
	$fileH{$type[1]}=1;
}
close I;
unless(-d $outpath){
	mkdir $outpath;
}
for my $type (keys %fileH){
	open $type,">$outpath/$type.fa";	
}
open I,"$kegg_fa";
my $fh;
while(<I>){
	chomp;
	if(m/^>([^:]+)/){
		if(!exists $class{$1}){
			die "Err:can not find the class about species $1\n";
		}
		$fh=$class{$1};
		print $fh "$_\n";
	}
	else{
		print $fh "$_\n";
	}
}
close I;

for my $type (keys %fileH){
	close $type;
}
