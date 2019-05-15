#!/usr/bin/perl -w

use strict;
use warnings;
use Data::Dumper;

#######################################################################################

if (@ARGV != 3) {
	print "Usage: perl $0 <xxx.miRanda.txt> <output_dir> <output_prefix>\n";
	exit;
}

my ($fIn, $dOut, $prefix) = @ARGV;
die "wrong format of prefix, should not contain '/'\n" if ($prefix =~/\//) ;

open (OUT,">$dOut/$prefix.txt") or die $!;
open (IN,"$fIn") or die $!;

print OUT "#Seq1,Seq2,Tot Score,Tot Energy,Max Score,Max Energy,Strand,Len1,Len2,Positions\n";
while (<IN>) {
	chomp;
	next if (/^$/) ;
	next unless (/^\s+Forward/) ;

	my $align_str = "";

	while (<IN>) {
		chomp;
		last if (/^Complete/) ;
		next if (/^$/ || /^\s+Energy/ || /^Score/ || /^Seq1,Seq2/ || /^\s+Forward/) ;
		$align_str .= $_."\n";;
	}

	my @lines = split /\n/, $align_str;
	foreach my $line (reverse @lines) {
		print OUT $line, "\n";
	}
	print OUT "\n";
}

close (IN) ;
close (OUT) ;

