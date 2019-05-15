#!/usr/bin/perl
use strict;
use warnings;

print "Hello, World...\n";
open IN,"$ARGV[0]";
open OUT ,">$ARGV[1]";
while(<IN>)
{
	chmod;
	my @line = split/\t/,$_;
	my @targets = split /;/,$line[1];
	next if($#targets>10);
	foreach my $target(@targets)
	{
		print OUT "$line[0]\t$target\n";
	}
}

close IN;
close OUT;