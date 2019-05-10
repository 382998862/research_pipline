#!/usr/bin/perl -w
use	strict;
use	Getopt::Long;
use	Data::Dumper;
use FindBin qw($Bin $Script);
use File::Basename qw(basename dirname);
my $program_name=basename($0);

my $ver="1.0";
############################################
my %opts;
GetOptions(\%opts,"i=s","r=s","l=s","s=s","o=s","h");
if (!defined($opts{i})||!defined($opts{r})||!defined($opts{l})||!defined($opts{s})||!defined($opts{o})||defined($opts{h})) {
	
	print << "	Usage End.";
	Description:SNP quannity,
		version:$ver
	Usage:

		-i      infile or directory(snp result)        must be given;

		-r      infile(reference  chromosome length)   must be given;

		-l      length of partitioned                  must be given;

		-s      step length                            must be given;

		-o     outfile key                             must be given;

	Usage End.
		exit;
}
###############Time
my $Time_Start;
$Time_Start = sub_format_datetime(localtime(time()));
print "\nStart $program_name Time :[$Time_Start]\n\n";
############################################
my $in=$opts{i};
my $out=$opts{o};
my $ref=$opts{r};
my $step=$opts{s};
my $dis=$opts{l};
my %chr;
my %snpn;
my %dis;
my $snp_dis="$out.dis";
open(RL,"$ref")||die"can't open the file [$ref]";
while (<RL>) {
	chomp;
	next if(/^$/||/\#/);
	$chr{(split/\s+/,$_)[0]}=(split/\s+/,$_)[1];
}
close RL;

#C1695093      504    A     G     29     19  R
my @in_file;
if (-f "$in") {
	push @in_file,$in;
}
else{
	@in_file=glob("$in/*");
}

foreach  my $file(sort @in_file) {
	open(IN,"$file")||die"can't open the [$file]";
	my $chr;
	while (<IN>) {
		chomp;
		next if (/^$/ || /\#/);
		my @snp=split(/\s+/, $_);
		my $num=$dis/$step;
		my $i=int($snp[1]/$step)+1;
		for (my $k=$i;$k>$i-$num;$k--){
			last if ($k<=0);
			$dis{$snp[0]}{$k}++;
			if ($snp[6]!~/[ATCG]/) {
				$dis{$snp[0]}{"$k/Heterozygosity"}++;
			}
			my $base=$snp[2];
			$base=~tr/ATGC/GCAT/;
			if (($snp[2] ne $snp[3] && $snp[2]=~/[AG]/ && $snp[3]=~/[CT]/)||($snp[2] ne $snp[3] && $snp[2]=~/[CT]/ && $snp[3]=~/[AG]/)) {
				$dis{$snp[0]}{"$k/transversion"}++;
			}elsif ($snp[2] ne $snp[3] && $base eq $snp[3]) {
				$dis{$snp[0]}{"$k/transition"}++;
			}else {
				print $_;
			}
		}
		$snpn{$snp[0]}{"all"}++;
		if (($snp[2] ne $snp[3] && $snp[2]=~/[AG]/ && $snp[3]=~/[CT]/)||($snp[2] ne $snp[3] && $snp[2]=~/[CT]/ && $snp[3]=~/[AG]/)) {
			$snpn{$snp[0]}{"transversion"}++;
		}
		my $bise=$snp[2];
		$bise=~tr/ATGC/GCAT/;
		if ($snp[2] ne $snp[3] && $bise eq $snp[3]) {
			$snpn{$snp[0]}{"transition"}++;
		}
		if ($snp[6]!~/[ATCG]/) {
			$snpn{$snp[0]}{"Heterozygosity"}++;
		}
	}
	close IN;
}

open(OUTd,">$snp_dis");
my $all=0;
my $transition=0;
my $transversion=0;
my $Heterozygosity=0;
foreach  (sort {$a cmp $b}keys %snpn) {
	my $portion= $chr{$_}%$step>($dis-$step) ? int($chr{$_}/$step) : int($chr{$_}/$step)+1 ;
	for (my $j=1;$j<=$portion ;$j++) {
		$dis{$_}{$j}=0 if (!defined($dis{$_}{$j}));
		$dis{$_}{"$j/transversion"}=0 if (!defined($dis{$_}{"$j/transversion"}));
		$dis{$_}{"$j/transition"}=0 if (!defined($dis{$_}{"$j/transition"}));
		$dis{$_}{"$j/Heterozygosity"}=0 if (!defined($dis{$_}{"$j/Heterozygosity"}));
		printf OUTd ("$_\t".(($j-1)*$step+1).'-'.($j*$step+($dis-$step))."\t$j\t".$dis{$_}{$j}."\t%5.2f\t%5.2f\t%5.2f\t%5.2f\n",$dis{$_}{$j}*100/($snpn{$_}{'all'}+0.0000000000000000001),$dis{$_}{"$j/transition"}*100/($dis{$_}{$j}+0.0000000000000000001),$dis{$_}{"$j/transversion"}*100/($dis{$_}{$j}+0.000000000000000000001),$dis{$_}{"$j/Heterozygosity"}*100/($dis{$_}{$j}+0.000000000000000000001));
	}
	$snpn{$_}{'all'}=0 if (!defined($snpn{$_}{'all'}));
	$snpn{$_}{'transition'}=0 if (!defined($snpn{$_}{'transition'}));
	$snpn{$_}{'transversion'}=0 if (!defined($snpn{$_}{'transversion'}));
	$snpn{$_}{'Heterozygosity'}=0 if (!defined($snpn{$_}{'Heterozygosity'}));
	$all+=$snpn{$_}{'all'};
	$transition+=$snpn{$_}{'transition'};
	$transversion+=$snpn{$_}{'transversion'};
	$Heterozygosity+=$snpn{$_}{'Heterozygosity'};
}
close OUTd;
################Time
my $Time_End;
$Time_End = sub_format_datetime(localtime(time()));
print "\nEnd $program_name Time :[$Time_End]\n\n";

###############Subs
sub sub_format_datetime {#Time calculation subroutine
    my($sec, $min, $hour, $day, $mon, $year, $wday, $yday, $isdst) = @_;
	$wday = $yday = $isdst = 0;
    sprintf("%4d-%02d-%02d %02d:%02d:%02d", $year+1900, $mon+1, $day, $hour, $min, $sec);
}
