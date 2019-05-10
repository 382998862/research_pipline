#!/usr/bin/perl -w
use strict;
use Data::Dumper;
use Getopt::Long;
my (%gene,%sum);
my (@gene_id,@group);
my $i = 0;
my ($In,$Outfile,$dis_num);
GetOptions(
			"i:s"=>\$In,
			"o:s"=>\$Outfile,
			"n:i"=>\$dis_num);
if (!defined $In || !defined $Outfile || !defined $dis_num){
	print "\n\tUsage:\tperl $0 -i fa_length_file -o out_prefix -n number\n\n" ;
	exit(1);
}
open IN,"$In" or die $!;
while (<IN>) {
	next if (/^\#/); 
	my @string = split/\t/;
	$gene{$string[0]}=$string[1];
}
close IN;
@gene_id = sort {$gene{$b} <=> $gene{$a}} keys %gene;
map {$group[$_] = [];} 0..($dis_num-1);
map {$sum{$_} = 0;} 0..($dis_num-1);
foreach  (@gene_id) {
	my ($the_min_num) = sort {$sum{$a} <=> $sum{$b}} 0..($dis_num-1);
	push @{$group[$the_min_num]},$_;
	$sum{$the_min_num} += $gene{$_};
}
foreach  (@group) {
	$i++;
	open OUT,">$Outfile.$i.list" or die $!;
	my @data = @{$_};
	print OUT join "\n",@data ;
	print OUT "\n";
	close OUT;
}

