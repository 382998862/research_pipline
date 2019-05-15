#/usr/bin/perl -w

use strict;

if (@ARGV!=3) {
	print "Useage: <RNAhybrid_out.txt> <outdir> <out_prefix> \n";
	exit;
}

my ($RNAhybrid_out,$od,$out_prefix)=@ARGV;

open IN,"$RNAhybrid_out" || die $!;
open OUT,">$od/$out_prefix.RNAhybrid.aln.txt" || die $!;
while (<IN>) {
	chomp;
	s/\s+$//;
	next if (/^$/) ;
	my ($miRNA,$miRNA_len,$mRNA,$mRNA_len,$mfe,$position,$p);
	my ($miRNA_unmatch,$miRNA_match,$mRNA_unmatch,$mRNA_match);
	if (/^target:\s(.*)/) {
		$mRNA=$1;
	}
	my $tmp1=<IN>;chomp $tmp1;
	if ($tmp1=~/^length:\s(\d+)/) {
		$mRNA_len=$1;
	}
	my $tmp2=<IN>;chomp $tmp2;
	if ($tmp2=~/^miRNA\s:\s(.*)/) {
		$miRNA=$1;
	}
	my $tmp3=<IN>;chomp $tmp3;
	if ($tmp3=~/^length:\s(\d+)/) {
		$miRNA_len=$1;
	}
	my $tmp4=<IN>;
	my $tmp5=<IN>;chomp $tmp5;
	if ($tmp5=~/^mfe:\s(.*)/) {
		$mfe=$1;
	}
	my $tmp6=<IN>;chomp $tmp6;
	if ($tmp6=~/^p-value:\s(.*)/) {
		$p=$1;
	}
	my $tmp7=<IN>;
	my $tmp8=<IN>;chomp $tmp8;
	if ($tmp8=~/^position\s+(\d+)/) {
		$position=$1;
	}

	my $mRNA_unmatch=<IN>;chomp $mRNA_unmatch;
	$mRNA_unmatch=~s/\s/-/g;
	my @m_unmatch=split//,$mRNA_unmatch;
	$mRNA_match=<IN>;chomp $mRNA_match;
	$mRNA_match=~s/\s/-/g;
	my @m_match=split//,$mRNA_match;
	$miRNA_match=<IN>;chomp $miRNA_match;
	$miRNA_match=~s/\s/-/g;
	my @mi_match=split//,$miRNA_match;
	$miRNA_unmatch=<IN>;chomp $miRNA_unmatch;
	$miRNA_unmatch=~s/\s/-/g;
	my @mi_unmatch=split//,$miRNA_unmatch;
	
	my ($mi_str,$match_str,$m_str);
	$mi_str =   "miRNA  3' ";
	$m_str =    "target 5' ";
	$match_str ="          ";
	foreach (10..$#mi_unmatch-3) {
		if ($m_unmatch[$_] ne "-" || $mi_unmatch[$_] ne "-") {
			$mi_str .= "$mi_unmatch[$_]";
			$m_str .= "$m_unmatch[$_]";
			$match_str .= " ";
		}
		else {
			$mi_str .= "$mi_match[$_]";
			$m_str .= "$m_match[$_]";
			if (($mi_match[$_] eq "G" && $m_match[$_] eq "U") || ($mi_match[$_] eq "U" && $m_match[$_] eq "G")) {
				$match_str .= ":";
			}
			else {
				$match_str .= "|";
			}
		}
	}
	$mi_str .= " 5'";
	$m_str .= " 3'";
	print OUT ">$miRNA\t$miRNA_len\t$mRNA\t$mRNA_len\t$position\t$mfe\t$p\n";
	print OUT "$mi_str\n$match_str\n$m_str\n\n";
}
close OUT;