#!/usr/bin/perl -w

use strict;
use warnings;
use Getopt::Long;
use Data::Dumper;
use FindBin qw($Bin $Script);
use File::Basename qw(basename dirname);
use Cwd qw(abs_path);
use newPerlBase;
my %config=%{readconf("$Bin/../../../../config/db_file.cfg")};


my $BEGIN_TIME=time();
my $version="1.0.0";
my @Original_ARGV=@ARGV;
#######################################################################################

# ==============================================================
# Get Options
# ==============================================================
my ($fIn1,$fIn2,$cfg,$od);

GetOptions(
				"help|?" =>\&USAGE,
				"snp:s"=>\$fIn1,
				"gff:s"=>\$fIn2,
				"od:s"=>\$od,
				) or &USAGE;
&USAGE unless ($fIn1 and $fIn2 );
#===============================================================
# Default optional value 
#===============================================================
$od||="./stat";
&MKDIR($od);
&begining;

my %type = ('AG'=>'R', 'GA'=>'R', 'CT'=>'Y', 'TC'=>'Y', 'GT'=>'K', 'TG'=>'K', 'AC'=>'M', 'CA'=>'M', 
	'CG'=>'S', 'GC'=>'S', 'AT'=>'W', 'TA'=>'W', 'A'=>'A', 'C'=>'C', 'G'=>'G', 'T'=>'T', 
	'CGT'=>'B', 'AGT'=>'D', 'ACT'=>'H', 'ACG'=>'V', 'ACGT'=>'N',
);
my %Retype;
foreach my $key (keys(%type)) {
	$Retype{$type{$key}}=$key;
}
#===============================================================
# Process
#===============================================================
my @header;
my %base;
my %series;
my %count;
my %snp_site;
my %gene_snp;
my %snp_rato;
my %all_snp;
my $infer="";
my %SNP_type;
my %use_loc;
my (%more,%less);
#Chr    Pos     Gene_id Ref     Alt     T03     Depth   AlleDp  T04     Depth   AlleDp  T05     Depth   AlleDp  T06     Depth   AlleDp  Effect  Codon_change
open IN1,"$fIn1" or die $!;
while (<IN1>) {
	chomp;
	my @line=split /\t/,$_;
	if ($_=~/\#/){
		@header=@line;
		next;
	}
	my $str="$line[0]:$line[1]";
	my $line0=$line[0];
	my $line1=$line[1];
	   $snp_site{$line0}{$line1}{'All'}='Intergenic';
	my $variance_type = variance_type($line[3], $line[4]);
	my $snp_tp="$line[3]->$line[4]";
		$SNP_type{$variance_type}{'All'}++;
		$SNP_type{$snp_tp}{'All'}++;
		$series{'All'}=1;
	for (my $i=5;$i<@line ;$i=$i+3) {
		my $x=$line[$i];
		next if($header[$i]=~/effect/i);
		next if($x eq "N");
		next if ($x eq $line[3]);
		$series{$header[$i]}=1;
		$count{$header[$i]}{"homo"}=0 if(!exists $count{$header[$i]}{"homo"}); 
		$count{$header[$i]}{"hete"}=0 if(!exists $count{$header[$i]}{"hete"});
		$count{$header[$i]}{$variance_type}=0 if(!exists $count{$header[$i]}{$variance_type});
		$SNP_type{$snp_tp}{$header[$i]}++;
		$SNP_type{$variance_type}{$header[$i]}++;
		$count{$header[$i]}{$variance_type}++;
		$count{$header[$i]}{"homo"}++ if($x eq 'A'  or $x eq 'T'  or $x eq 'G'  or $x eq 'C');
		$count{$header[$i]}{"hete"}++ if($x ne 'A' and $x ne 'T' and $x ne 'G' and $x ne 'C');
		$snp_site{$line[0]}{$line[1]}{$header[$i]}='Intergenic';
	}

}
close IN1;


open (IN2, $fIn2) or die $!;
while (<IN2>) {
	chomp;
	next if (/^$/ || /^\#/) ;

	my ($chro, undef, $feture, $start, $end, undef, $strand, undef, $attribute) = split;
	next unless (($feture eq "gene") or ($feture eq "Gene")  ) ;
	my ($geneID) = $attribute =~/ID=([^;]+)/;
	($start, $end) = sort {$a <=> $b} ($start, $end);
    my $gene_len = $end-$start+1;
	&find_snps_in_gene($chro,$start,$end,$geneID,$gene_len);
}
close(IN2);

#my @types=('transition','A->G','G->A','C->T','T->C','transversion','A->C','C->A','A->T','T->A','C->G','G->C','G->T','T->G');
my @types=('A->G','G->A','C->T','T->C','A->C','C->A','A->T','T->A','C->G','G->C','G->T','T->G');

open STAT,">$od/All.snp_type.stat";
print STAT "#SNP type\t",join ("\t",sort keys %series),"\n";
foreach my $tp(@types) {
	print STAT $tp;
	foreach my $samples (sort keys %series) {
		my $value=(exists $SNP_type{$tp}{$samples})?$SNP_type{$tp}{$samples}:0;
		print STAT "\t$value";
	}
	print STAT "\n";
}
close STAT;


open STAT,">$od/AllSample.snp.stat";
print STAT "#BMK-ID\tSNP Number\tGenic SNP\tIntergenic SNP\tTransition\tTransversion\tHeterozygosity\n";
foreach my $sam (sort keys %count) {
	my $sum=$count{$sam}{homo}+$count{$sam}{hete};
	my $intergenic=$sum-$gene_snp{$sam};
	my $transition = sprintf("%.2f", $count{$sam}{'transition'} / $sum * 100)."%";
	my $transversion=sprintf("%.2f", $count{$sam}{'transversion'} / $sum * 100)."%";
	my $hete= sprintf("%.2f", $count{$sam}{hete} / $sum * 100)."%";
	print STAT "$sam\t$sum\t$gene_snp{$sam}\t$intergenic\t$transition\t$transversion\t$hete\n";
}
close STAT;


open (STAT,">$od/AllSample.SNP_density.stat") ||die $!;
print STAT "#Sample\tInterval\tGeneNum\n";
foreach my $sam ( sort keys %snp_rato) {
	foreach my $types (sort keys  %{$snp_rato{$sam}}) {
		print STAT "$sam\t$types\t$snp_rato{$sam}{$types}\n";
	}
}
close (STAT);
my $Rscript = $config{Rscript};
my $cmd = "$Rscript $Bin/dodgedBar.r --infile $od/AllSample.SNP_density.stat --outfile $od/AllSample.SNP_density.png --x.col 2 --group.col 1 --y.col 3 --x.lab \"SNP Number per Kb\" --group.lab \"Sample\" --y.lab \"Number of Gene\" --title.lab \"SNP Density\" --legend.col 1 >/dev/null 2>&1 ";
#print "$cmd\n";
system $cmd;


&end;
# ==============================================================
# sub function
# ==============================================================
sub MKDIR{
	my ($dirname)=@_;
	if (-e $dirname) {
		print "$dirname exists!";
	}else{
		mkdir $dirname;
	}
}
sub begining{
	my $cmd=$0." ";
	$cmd.=join(" ",@Original_ARGV);
	my $time=&GetTime;
	print STDOUT "$cmd\n START:$time\n";
}
sub end{
	print STDOUT "\nDone. Total elapsed time : ",time()-$BEGIN_TIME,"s\n";
}
sub GetTime {
	my ($sec, $min, $hour, $day, $mon, $year, $wday, $yday, $isdst)=localtime(time());
	return sprintf("%4d-%02d-%02d %02d:%02d:%02d", $year+1900, $mon+1, $day, $hour, $min, $sec);
}
sub ReadConfig{
	my($cfg)=@_;
	open (IN,"$cfg") || die "$!";
	my %para;
	my %sample;
	while (<IN>) {
		chomp;
		s/\s+//;
		s/\s+$//;
		next if (/\#/);
		next if (/$/);
		my @tmp=split /\s+/,$_;
		if ($tmp[0]=~m/Sample/) {
			my $fq1=<IN>;chomp $fq1;
			my $fq2=<IN>;chomp $fq2;
			my @fq_1=split /\s+/,$fq1;
			$sample{$tmp[1]}{FQ1}=$fq_1[1];
			my @fq_2=split /\s+/,$fq2;
			$sample{$tmp[1]}{FQ2}=$fq_2[1];
		}
		$para{$tmp[0]}=$tmp[1];
	}
	close IN;

	return (\%para,\%sample);
}
sub find_snps_in_gene {#
	my ($chro, $start, $end,$geneID,$gene_len) = @_;
	my %lishi;
#	my $start_int = pos_unit($start);
#	my $end_int   = pos_unit($end);
	my $start_int = $start;
	my $end_int   = $end;
	for(my $i=$start_int; $i<=$end_int; $i++)
	{
		next unless exists $snp_site{$chro}{$i};
        if (exists $use_loc{$chro}{$i}) {
            foreach my $samples (sort keys  %{$snp_site{$chro}{$i}}) {
                $snp_site{$chro}{$i}{$samples}=$geneID;
                $lishi{$geneID}{$samples}++;
            }
        }else {
            foreach my $samples (sort keys  %{$snp_site{$chro}{$i}}) {
                $snp_site{$chro}{$i}{$samples}=$geneID;
                $lishi{$geneID}{$samples}++;
                $gene_snp{$samples}++;
                $use_loc{$chro}{$i}=1;
            }        
        }

	}
	foreach my $samples (sort keys  %{$lishi{$geneID}}) {
		my $rate=$lishi{$geneID}{$samples}*1000/$gene_len;

		$snp_rato{$samples}{'0-1'}++ if ($rate<1);
		$snp_rato{$samples}{'1-2'}++ if ($rate>=1 and $rate<2);
		$snp_rato{$samples}{'2-3'}++ if ($rate>=2 and $rate<3);
		$snp_rato{$samples}{'3-4'}++ if ($rate>=3 and $rate<4);
		$snp_rato{$samples}{'4-5'}++ if ($rate>=4 and $rate<5);
		$snp_rato{$samples}{'5-6'}++ if ($rate>=5 and $rate<6);
		$snp_rato{$samples}{'6-7'}++ if ($rate>=6 and $rate<7);
		$snp_rato{$samples}{'7-8'}++ if ($rate>=7 and $rate<8);
		$snp_rato{$samples}{'8~'}++  if ($rate>8);
	}
}

sub variance_type
{
	my ($ref, $variance) = @_;
	my $type;
	my $genotype = $ref.$variance;
	if($genotype eq 'AG' || $genotype eq 'CT' || $genotype eq 'GA' || $genotype eq 'TC')
	{
		$type = 'transition';
	}else
	{
		$type = 'transversion';
	}
	return $type;
}


=c
sub pos_unit {
	my ($pos) = @_;
	my $unit  = 1e4;
	my $pos_int = int($pos/$unit);
	return $pos_int;
}
=cut

sub USAGE {#
	my $usage=<<"USAGE";
Program: $0
Version: $version
Contact: Shi Tongwei <shitw\@biomarker.com.cn> 
=====================================================================================================
Discription:
=====================================================================================================

Usage:
  Options:
  -snp	<file>	required	parent snplist file
  -gff  <file>  known_and_newgene.gff
  -od	<str>	optional	Directory where output file produced,optional,default [./]
  -h		Help

USAGE
	print $usage;
	exit;
}
