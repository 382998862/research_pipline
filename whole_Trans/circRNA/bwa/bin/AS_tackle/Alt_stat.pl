#!/usr/bin/perl -w
use strict;
my $version= "1.0";

if (@ARGV!=3) {
	print "\n\tFunction: Extract AS stat Result from SpliceGrapher gff\n\n";
	print "\n\tVersion: $version\n\n";
	print "\tUsage: perl Alt_graphpdf.pl <SpliceGrapher.gff> <out_dir> <Sample> \n\n";
	exit;
}

open IN,"$ARGV[0]" || die $!;
my %graph_info;
my %non_AS_array;
my %AS_array;
my $graph_name;
while (<IN>) {
	chomp;
	next if (/^$/ || /^\#/) ;
	my @tmp=split /\t+/,$_;
	if ($tmp[2]=~/graph/) {
		$tmp[-1]=~/ID=(.*)/;
		$graph_name=$1;
		$graph_info{$tmp[0]}{$1}{Start}=$tmp[3];
		$graph_info{$tmp[0]}{$1}{END}=$tmp[4];
		$graph_info{$tmp[0]}{$1}{Strand}=$tmp[6];
	}
	else {
		if ($tmp[-1]!~/AltForm/) {
			push @{$non_AS_array{$tmp[0]}{$graph_name}},($tmp[3],$tmp[4]);
		}
		else {
			$tmp[-1]=~/ID=(\S+);*.AltForm=(.*);Isoforms/;
			my $alt_form=$2;
			if ($alt_form =~ /SE/) {
				$AS_array{$tmp[0]}{$graph_name}{SE}{$tmp[3]}=$tmp[4];
			}
			if ($alt_form =~ /IR/) {
				$AS_array{$tmp[0]}{$graph_name}{IR}{$tmp[3]}=$tmp[4];
			}
			if ($alt_form =~ /A3/) {
				$AS_array{$tmp[0]}{$graph_name}{A3}{$tmp[3]}=$tmp[4];
			}
			if ($alt_form =~ /A5/) {
				$AS_array{$tmp[0]}{$graph_name}{A5}{$tmp[3]}=$tmp[4];
			}
			if ($alt_form =~ /AI/) {
				$AS_array{$tmp[0]}{$graph_name}{AI}{$tmp[3]}=$tmp[4];
			}
			if ($alt_form =~ /AT/) {
				$AS_array{$tmp[0]}{$graph_name}{AT}{$tmp[3]}=$tmp[4];
			}
		}
	}
}
close IN;

open OUT,">$ARGV[1]/$ARGV[2].AltSplice_info.xls" || die;
open OUT1,">$ARGV[1]/$ARGV[2].IntronRetention.xls" || die;
open OUT2,">$ARGV[1]/$ARGV[2].skippedExon.xls" || die;
open OUT3,">$ARGV[1]/$ARGV[2].Alt3splice.xls" || die;
open OUT4,">$ARGV[1]/$ARGV[2].Alt5splice.xls" || die;
open OUT5,">$ARGV[1]/$ARGV[2].AltFirstExon.xls" || die;
open OUT6,">$ARGV[1]/$ARGV[2].AltLastExon.xls" || die;
#open OUT7,">$ARGV[1]/$ARGV[2].MutuallyExon.xls" || die;

open TE,">$ARGV[1]/$ARGV[2].Specialgene.list" || die;

my %AS_stat;
my $AS_total=0;
foreach my $chr (sort keys %AS_array) {
	foreach my $gene (sort keys %{$AS_array{$chr}}) {
		my $Afe_Ale=0;
		my %AS_info;
		my @exon=sort {$a<=>$b} @{$non_AS_array{$chr}{$gene}} if (defined @{$non_AS_array{$chr}{$gene}}) ;
		foreach my $AS_type (sort keys %{$AS_array{$chr}{$gene}}) {
			foreach my $start (sort{$a<=$b} keys %{$AS_array{$chr}{$gene}{$AS_type}}) {
				if ($AS_type eq "SE") {
					$AS_info{SE}{$start}=$AS_array{$chr}{$gene}{$AS_type}{$start};
				}
				if ($AS_type eq "IR") {
					for (my $i=1; $i<=($#exon-2); $i+=2) {
						if ($exon[$i] >= $start && $exon[$i+1] <= $AS_array{$chr}{$gene}{$AS_type}{$start}) {
							$AS_info{IR}{$exon[$i]} = $exon[$i+1];
						}
					}
				}
				if ($AS_type eq "A3") {
					if ($graph_info{$chr}{$gene}{Strand} eq "+") {
						push @{$AS_info{A3}{$AS_array{$chr}{$gene}{$AS_type}{$start}}},$start;
					}
					else {
						push @{$AS_info{A3}{$start}},$AS_array{$chr}{$gene}{$AS_type}{$start};
					}
				}
				if ($AS_type eq "A5") {
					if ($graph_info{$chr}{$gene}{Strand} eq "+") {
						push @{$AS_info{A5}{$start}},$AS_array{$chr}{$gene}{$AS_type}{$start};
					}
					else {
						push @{$AS_info{A5}{$AS_array{$chr}{$gene}{$AS_type}{$start}}},$start;
					}
				}
				if ($AS_type eq "AI") {
					if ($graph_info{$chr}{$gene}{Strand} eq "+") {
						push @{$AS_info{AFE}{$Afe_Ale}{New}},($start,$AS_array{$chr}{$gene}{$AS_type}{$start});
						push @{$AS_info{AFE}{$Afe_Ale}{Old}},($exon[0],$exon[1]);
						$Afe_Ale++;
					}
					else {
						push @{$AS_info{AFE}{$Afe_Ale}{New}},($start,$AS_array{$chr}{$gene}{$AS_type}{$start});
						push @{$AS_info{AFE}{$Afe_Ale}{Old}},($exon[-2],$exon[-1]);
						$Afe_Ale++;
					}
				}
				if ($AS_type eq "AT") {
					if ($graph_info{$chr}{$gene}{Strand} eq "+") {
						push @{$AS_info{ALE}{$Afe_Ale}{New}},($start,$AS_array{$chr}{$gene}{$AS_type}{$start});
						push @{$AS_info{ALE}{$Afe_Ale}{Old}},($exon[-2],$exon[-1]);
						$Afe_Ale++;
					}
					else {
						push @{$AS_info{ALE}{$Afe_Ale}{New}},($start,$AS_array{$chr}{$gene}{$AS_type}{$start});
						push @{$AS_info{ALE}{$Afe_Ale}{Old}},($exon[0],$exon[1]);
						$Afe_Ale++;
					}
				}
			}
		}
		foreach my $AS_type (sort keys %AS_info) {
			if ($AS_type eq "SE") {
				foreach my $start (sort keys %{$AS_info{SE}}) {
					print OUT "$gene\t$chr\t$graph_info{$chr}{$gene}{Strand}\tskippedExon\t$start\t$AS_info{SE}{$start}\n";
					print OUT2 "$gene\t$chr\t$graph_info{$chr}{$gene}{Strand}\tskippedExon\t$start\t$AS_info{SE}{$start}\n";
					$AS_stat{skippedExon}++;
					$AS_total++;
				}
			}
			if ($AS_type eq "IR") {
				foreach my $start (sort keys %{$AS_info{IR}}) {
					print OUT "$gene\t$chr\t$graph_info{$chr}{$gene}{Strand}\tIntronRetention\t$start\t$AS_info{IR}{$start}\n";
					print OUT1 "$gene\t$chr\t$graph_info{$chr}{$gene}{Strand}\tIntronRetention\t$start\t$AS_info{IR}{$start}\n";
					$AS_stat{IntronRetention}++;
					$AS_total++;
				}
			}
			if ($AS_type eq "A3") {
				if ($graph_info{$chr}{$gene}{Strand} eq "+") {
					foreach my $end (sort{$a<=>$b} keys %{$AS_info{A3}}) {
						my @starts =sort {$a<=>$b} @{$AS_info{A3}{$end}};
						if ($#starts<1) {
							print TE "$gene\tAlt3' site\t$starts[0]\t$end\n";
						}
						elsif ($#starts==1) {
							my $alt_site = $starts[1]-1;
							print OUT "$gene\t$chr\t$graph_info{$chr}{$gene}{Strand}\tAlt3' site\t$starts[0]\t$end\t$starts[1]\t$end\t$starts[0]:$alt_site\n";
							print OUT3 "$gene\t$chr\t$graph_info{$chr}{$gene}{Strand}\tAlt3' site\t$starts[0]\t$end\t$starts[1]\t$end\t$starts[0]:$alt_site\n";
							$AS_stat{"Alt3' site"}++;
							$AS_total++;
						}
						else {
							for (1..$#starts) {
								my $alt_site = $starts[$_]-1;
								print OUT "$gene\t$chr\t$graph_info{$chr}{$gene}{Strand}\tAlt3' site\t$starts[0]\t$end\t$starts[$_]\t$end\t$starts[0]:$alt_site\n";
								print OUT3 "$gene\t$chr\t$graph_info{$chr}{$gene}{Strand}\tAlt3' site\t$starts[0]\t$end\t$starts[$_]\t$end\t$starts[0]:$alt_site\n";
								$AS_stat{"Alt3' site"}++;
								$AS_total++;
							}
						}
					}
				}
				else {
					foreach my $start (sort{$a<=>$b} keys %{$AS_info{A3}}) {
						my @ends =sort {$a<=>$b} @{$AS_info{A3}{$start}};
						if ($#ends<1) {
							print TE "$gene\tAlt3' site\t$ends[0]\t$start\n";
						}
						elsif ($#ends==1) {
							my $alt_site = $ends[0]+1;
							print OUT "$gene\t$chr\t$graph_info{$chr}{$gene}{Strand}\tAlt3' site\t$start\t$ends[0]\t$start\t$ends[1]\t$alt_site:$ends[1]\n";
							print OUT3 "$gene\t$chr\t$graph_info{$chr}{$gene}{Strand}\tAlt3' site\t$start\t$ends[0]\t$start\t$ends[1]\t$alt_site:$ends[1]\n";
							$AS_stat{"Alt3' site"}++;
							$AS_total++;
						}
						else {
							for (1..$#ends) {
								my $alt_site = $ends[0]+1;
								print OUT "$gene\t$chr\t$graph_info{$chr}{$gene}{Strand}\tAlt3' site\t$start\t$ends[0]\t$start\t$ends[$_]\t$alt_site:$ends[$_]\n";
								print OUT3 "$gene\t$chr\t$graph_info{$chr}{$gene}{Strand}\tAlt3' site\t$start\t$ends[0]\t$start\t$ends[$_]\t$alt_site:$ends[$_]\n";
								$AS_stat{"Alt3' site"}++;
								$AS_total++;
							}
						}
					}
				}
			}
			if ($AS_type eq "A5") {
				if ($graph_info{$chr}{$gene}{Strand} eq "+") {
					foreach my $start (sort keys %{$AS_info{A5}}) {
						my @ends =sort {$a<=>$b} @{$AS_info{A5}{$start}};
						if ($#ends<1) {
							print TE "$gene\tAlt3' site\t$ends[0]\t$start\n";
						}
						elsif ($#ends==1) {
							my $alt_site = $ends[0]+1;
							print OUT "$gene\t$chr\t$graph_info{$chr}{$gene}{Strand}\tAlt5' site\t$start\t$ends[0]\t$start\t$ends[1]\t$alt_site:$ends[1]\n";
							print OUT3 "$gene\t$chr\t$graph_info{$chr}{$gene}{Strand}\tAlt5' site\t$start\t$ends[0]\t$start\t$ends[1]\t$alt_site:$ends[1]\n";
							$AS_stat{"Alt5' site"}++;
							$AS_total++;
						}
						else {
							for (1..$#ends) {
								my $alt_site = $ends[0]+1;
								print OUT "$gene\t$chr\t$graph_info{$chr}{$gene}{Strand}\tAlt5' site\t$start\t$ends[0]\t$start\t$ends[$_]\t$alt_site:$ends[$_]\n";
								print OUT3 "$gene\t$chr\t$graph_info{$chr}{$gene}{Strand}\tAlt5' site\t$start\t$ends[0]\t$start\t$ends[$_]\t$alt_site:$ends[$_]\n";
								$AS_stat{"Alt5' site"}++;
								$AS_total++;
							}
						}
					}
				}
				else {
					foreach my $end (sort{$a<=>$b} keys %{$AS_info{A5}}) {
						my @starts =sort {$a<=>$b} @{$AS_info{A5}{$end}};
						if ($#starts<1) {
							print TE "$gene\tAlt3' site\t$starts[0]\t$end\n";
						}
						elsif ($#starts==1) {
							my $alt_site = $starts[1]-1;
							print OUT "$gene\t$chr\t$graph_info{$chr}{$gene}{Strand}\tAlt5' site\t$starts[0]\t$end\t$starts[1]\t$end\t$starts[0]:$alt_site\n";
							print OUT4 "$gene\t$chr\t$graph_info{$chr}{$gene}{Strand}\tAlt5' site\t$starts[0]\t$end\t$starts[1]\t$end\t$starts[0]:$alt_site\n";
							$AS_stat{"Alt5' site"}++;
							$AS_total++;
						}
						else {
							for (1..$#starts) {
								my $alt_site = $starts[$_]-1;
								print OUT "$gene\t$chr\t$graph_info{$chr}{$gene}{Strand}\tAlt5' site\t$starts[0]\t$end\t$starts[$_]\t$end\t$starts[0]:$alt_site\n";
								print OUT4 "$gene\t$chr\t$graph_info{$chr}{$gene}{Strand}\tAlt5' site\t$starts[0]\t$end\t$starts[$_]\t$end\t$starts[0]:$alt_site\n";
								$AS_stat{"Alt5' site"}++;
								$AS_total++;
							}
						}
					}
				}
			}
			if ($AS_type eq "AFE") {
				foreach my $num (sort keys %{$AS_info{AFE}}) {
					print OUT "$gene\t$chr\t$graph_info{$chr}{$gene}{Strand}\tAltFirstExon\t$AS_info{AFE}{$num}{New}[0]\t$AS_info{AFE}{$num}{New}[1]\t$AS_info{AFE}{$num}{Old}[0]\t$AS_info{AFE}{$num}{Old}[1]\n";
					print OUT5 "$gene\t$chr\t$graph_info{$chr}{$gene}{Strand}\tAltFirstExon\t$AS_info{AFE}{$num}{New}[0]\t$AS_info{AFE}{$num}{New}[1]\t$AS_info{AFE}{$num}{Old}[0]\t$AS_info{AFE}{$num}{Old}[1]\n";
					$AS_stat{AltFirstExon}++;
					$AS_total++;
				}
			}
			if ($AS_type eq "ALE") {
				foreach my $num (sort keys %{$AS_info{ALE}}) {
					print OUT "$gene\t$chr\t$graph_info{$chr}{$gene}{Strand}\tAltLastExon\t$AS_info{ALE}{$num}{New}[0]\t$AS_info{ALE}{$num}{New}[1]\t$AS_info{ALE}{$num}{Old}[0]\t$AS_info{ALE}{$num}{Old}[1]\n";
					print OUT6 "$gene\t$chr\t$graph_info{$chr}{$gene}{Strand}\tAltLastExon\t$AS_info{ALE}{$num}{New}[0]\t$AS_info{ALE}{$num}{New}[1]\t$AS_info{ALE}{$num}{Old}[0]\t$AS_info{ALE}{$num}{Old}[1]\n";
					$AS_stat{AltLastExon}++;
					$AS_total++;
				}
			}
		}
	}
}
close OUT;
close OUT1;
close OUT2;
close OUT3;
close OUT4;
close OUT5;
close OUT6;
#close OUT7;

close TE;

open OUT,">$ARGV[1]/$ARGV[2].Altsplice.stat.xls" || die;
foreach my $key (sort keys %AS_stat) {
	print OUT "$key\t$AS_stat{$key}\n";
}
print OUT "total_AS\t$AS_total\n";
close OUT;
