#!/usr/bin/perl -w
my $version="1.0.0";
use strict;
use Getopt::Long;
use Data::Dumper;
use FindBin qw($Bin $Script);
use File::Basename qw(basename dirname);
######################请在写程序之前，一定写明时间、程序用途、参数说明；每次修改程序时，也请做好注释工作
my ($In,$out,$gff);
GetOptions(
			"help|?"=>\&USAGE,
			"i:s"=>\$In,
			"gff:s"=>\$gff,
			"o:s"=>\$out,

			) or &USAGE;
&USAGE unless ($In and $out and $gff);
#############################################################################################
my $BEGIN_TIME=time();
print "\nStart Time :[$BEGIN_TIME]\n\n";
#############################################################################################
open (IN,"$In") || die "Can't open $In\n";
open (GFF,"$gff") || die "Can't open $gff\n";
open (OUT,">$out") || die "Can't creat $out\n";
my %gff = ();
while(<GFF>)
{
	s/[\r\n]$//;
	next unless(/gene/);
	my @line=split/\t/,$_;
	my $gene_id = pop @line;
	my $chr = shift @line;
	$gene_id =~s/(ID=)//;
	$gene_id =~s/(;)//;
	$gff{$gene_id}{$chr}= "$line[2]\t$line[3]";
}
close GFF;

my %circ = ();
<IN>;
while(<IN>)
{
	s/[\r\n]$//;
	next if(/#/);
	my @line=split/\t/,$_;
	my $circRNA_id = $line[0];
	my $chr =(split /:/,$circRNA_id)[0];
	my $start_end = (split /:/,$circRNA_id)[1];
	my $circRNA_start=(split /\|/,$start_end)[0];
	my $circRNA_end=(split /\|/,$start_end)[1];
	my $gene = $line[-1];
	$gene =~s/,$//;
	if($gene=~/,/)
	{
		my @gene_list = split/,/,$gene;
		foreach my $gene_id(@gene_list)
		{
			$circ{$gene_id}{$chr}{$circRNA_id} = "$circRNA_start\t$circRNA_end";
		}
	}
	else
	{
		$circ{$gene}{$chr}{$circRNA_id} = "$circRNA_start\t$circRNA_end";
	}
}

close IN;
print OUT "circRNA_ID\tchr\tcircRNA_start\tcircRNA_end\tgene_id\tstart\tend\n";
foreach my $gene_id ( keys %circ){
	foreach my $chr (keys %{$circ{$gene_id}}) {
		next unless(exists $gff{$gene_id}{$chr});
		my $mRNAPostion = $gff{$gene_id}{$chr};
		foreach my $circRNA_id (keys %{$circ{$gene_id}{$chr}})
		{
			my $circPostion = $circ{$gene_id}{$chr}{$circRNA_id};
			print OUT "$circRNA_id\t$chr\t$circPostion\t$gene_id\t$mRNAPostion\n";
		}

	}
}
close OUT;
#######################################################################################
print STDOUT "\nDone. Total elapsed time : ",time()-$BEGIN_TIME,"s\n";
##################################Subs######################################################
sub sub_format_datetime {#Time calculation subroutine
    my($sec, $min, $hour, $day, $mon, $year, $wday, $yday, $isdst) = @_;
	$wday = $yday = $isdst = 0;
    sprintf("%4d-%02d-%02d %02d:%02d:%02d", $year+1900, $mon+1, $day, $hour, $min, $sec);
}

######
sub USAGE {#
my $usage=<<"USAGE";
	ProgramName:
	Version:	$version
	Contact:	lv ran <lvr\@biomarker.com.cn>
	Program Date:   2014/9/3
	Usage:
	Options:
		-i			<file> input file,forced
		-gff		<file> gff file,forced
		-o			<file> output file,forced

USAGE
		print $usage;
		exit;
}
