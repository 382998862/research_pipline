#!/usr/bin/perl -w
use strict;
use Data::Dumper;

if (@ARGV != 2) {
	print "\n\tFunction: Tackle Cuffmerge track2name.list , Create new formated list\n\n";
	print "\tInfile: AAA.newGene.track.list & AAA.track2name.list\n\n";
	print "\tUsage: perl Track2name_format.pl <in_index> <newGene_prefix>\n\n";
	exit;
}



my %name_track;
my %track_name;
my %info;


open IN,"$ARGV[0].track2name.list" || die $!;
while (<IN>) {
	chomp;
	next if (/^$/ || /^\#/);
	my @tmp = split/\t+/,$_;
	next if ($tmp[4]!~/=/ && $tmp[4]!~/j/);
	$name_track{$tmp[1]}{$tmp[0]}=1;
	$track_name{$tmp[0]}{$tmp[1]}=1;
	$info{$tmp[0]}{$tmp[2]}=$_;
}
close IN;

my %multitrack2name;
my %multiname2track;

my %multi_all;

open OUT1,">$ARGV[0].uniq_track2name.list.info" || die $!;
foreach my $gene_track (sort keys %track_name) {
	my @gene_names=(sort keys %{$track_name{$gene_track}});
	if ($#gene_names==0) {
		my @tracks=(sort keys %{$name_track{$gene_names[0]}});
		if ($#tracks==0) {
			foreach (sort keys %{$info{$gene_track}}) {
				print OUT1 "$info{$gene_track}{$_}\n";
			}
		}
		else {
			my $j = 2;
			foreach my $T (@tracks) {
				my $main_index = 0;
				foreach (sort keys %{$info{$T}}) {
					if ($info{$T}{$_}=~/Glyma(.*)\.1\t=/) {
						$main_index=1;
					}
				}
				if ($main_index==1) {
					$multitrack2name{$gene_names[0]}{$T}=1;
				}
				else {
					my $new_gene_name="$gene_names[0]"."-$j";
					$multitrack2name{$new_gene_name}{$T}=1;
					$j++;
				}
			}
		}
	}
	else {
		my ($name_str,$main_name);
		my $main_index = 0;
		foreach (sort @gene_names) {
			my @tracks=(sort keys %{$name_track{$_}});
			if ($#tracks==0) {
				$main_index=1;
				$name_str.="$_".",";
				$main_name = (sort @gene_names)[0];
			}
			else {
				$multi_all{$gene_track}{$_}=1;
			}
		}
		if ($main_index==1) {
			$name_str=~s/,$//;
			$multiname2track{$gene_track}{$main_name}=$name_str;
		}
	}
}
close OUT1;


open OUT1,">$ARGV[0].multi_track2name.list.info" || die $!;
foreach my $gene_name (sort keys %multitrack2name) {
	foreach my $gene_track (sort keys %{$multitrack2name{$gene_name}}) {
		foreach (sort keys %{$info{$gene_track}}) {
			my @tmp=split /\t+/,$info{$gene_track}{$_};
			$tmp[1]=$gene_name;
			print OUT1 join "\t",@tmp;
			print OUT1 "\n";
		}
	}
}
close OUT1;

open OUT1,">$ARGV[0].multi_name2track.list.info" || die $!;
foreach my $gene_track (sort keys %multiname2track) {
	foreach my $gene_name (sort keys %{$multiname2track{$gene_track}}) {
		foreach (sort keys %{$info{$gene_track}}) {
			my @tmp=split /\t+/,$info{$gene_track}{$_};
			$tmp[1]=$gene_name;
			print OUT1 join "\t",@tmp;
			print OUT1 "\t$multiname2track{$gene_track}{$gene_name}";
			print OUT1 "\n";
		}
	}
}
close OUT1;