#!/usr/bin/perl -w
use strict;

if (@ARGV != 3)
{
	print "\n";
	print "    Function: Change splice gene name for ASprofile\n\n";
	print "    Usage: perl changename.pl <cuffcmp.transcripts.gtf.tmap> <transcripts.gtf> <outdir>\n\n";
	exit;
}

my %Alt_gene;

open (IN, $ARGV[0]) || die "Open $ARGV[0] failed!\n";
while (<IN>) {
	chomp;
	next if (/^$/ || /^\#/);
	
	my ($alt_gene,$code,undef,undef,$alt_trans) = split /\s+/, $_;
	next if ($alt_gene eq '-' or $code eq 'u') ;
	$Alt_gene{$alt_trans}=$alt_gene;
}
close IN;

open (IN, $ARGV[1]) || die "Open $ARGV[1] failed!\n";
open (OUT,">$ARGV[2]/transcripts.gtf") || die "Open $ARGV[2]/transcripts.gtf failed!\n";
while (<IN>) {
	chomp;
	next if (/^$/ || /^\#/);
	
	my @inform = split /\"/, $_;
	
	if (exists $Alt_gene{$inform[3]}) {
		$inform[1]=$Alt_gene{$inform[3]};
		print OUT join('"',@inform),"\n";
	}
	else {next;}
}
close IN;
close OUT;

open (IN, $ARGV[0]) || die "Open $ARGV[0] failed!\n";
open (OUT,">$ARGV[2]/cuffcmp.transcripts.gtf.tmap") || die "Open $ARGV[2]/cuffcmp.transcripts.gtf.tmap failed!\n";

while (<IN>) {
	chomp;
	next if (/^$/ || /^\#/);
	
	my @inform = split /\t/, $_;
	
	if (exists $Alt_gene{$inform[4]}) {
		$inform[3]=$Alt_gene{$inform[4]};
		print OUT join("\t",@inform),"\n";
	}
	else{next;}
}
close IN;
close OUT;


