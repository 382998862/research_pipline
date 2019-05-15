#!/usr/bin/perl
use strict;
use warnings;

open IN,"$ARGV[0]" or die $!;
open OUT,">$ARGV[1]"or die $!;
my %hash;
while(<IN>)
{
	chomp;
	next if(/#/);
	my @line = split /\t/,$_;
	my $gene_id = $line[-1];
	next if($gene_id =~"n/a");
	$gene_id =~s/,$//;
	my @gene = split/,/,$gene_id;
	foreach my $gene (@gene)
	{
		if(!exists $hash{$gene})
		{
			$hash{$gene}=1;
			print OUT "$gene\n";
		}
	}

}
close IN;
close OUT;

