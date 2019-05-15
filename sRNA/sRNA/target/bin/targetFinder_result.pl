#/usr/bin/perl -w

use strict;

if (@ARGV!=3) {
	print "Useage: <Targetfinder_out.txt> <outdir> <out_prefix> \n";
	exit;
}

my ($Targetfinder_out,$od,$out_prefix)=@ARGV;

open IN,"$Targetfinder_out" || die $!;
open OUT,">$od/$out_prefix.TargetFinder.aln.txt" || die $!;
while (<IN>) {
	chomp;
	s/\s+$//;
	next if (/^$/) ;
	my $info=$_;
	my ($miRNA,$mRNA,$score,$range);
	my ($miRNA_match,$mRNA_match,$match);
	if ($info=~/^query=/) {
		<IN>;
		$mRNA_match=<IN>;chomp $mRNA_match;
		$match=<IN>;chomp $match;
		$match=~s/:/\|/g;
		$match=~s/\./:/g;
		$miRNA_match=<IN>;chomp $miRNA_match;
		$miRNA_match=~s/^query/miRNA/;
		$info=~/^query=(.*),\starget=(.*),\sscore=(.*),\srange=(\d+-\d+),/;
		print OUT ">$1\t$2\t$3\t$4\n";
		print OUT "$miRNA_match\n";
		print OUT "$match\n";
		print OUT "$mRNA_match\n\n";
	}
}
close OUT;