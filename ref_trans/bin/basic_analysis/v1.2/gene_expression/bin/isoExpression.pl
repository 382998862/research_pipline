#!/usr/bin/perl -w
use strict;
use warnings;
use Getopt::Long;
use Data::Dumper;
use FindBin qw($Bin $Script);
use File::Basename qw(basename dirname);
my $BEGIN_TIME=time();
my $version="1.0.0";
#######################################################################################
#
#######################################################################################
# ------------------------------------------------------------------
# GetOptions
# ------------------------------------------------------------------
my ($cuffdiff_dir,$final_track_dir,$outdir);

GetOptions(
				"help|?" =>\&USAGE,
				"c:s"=>\$cuffdiff_dir,
				"f:s"=>\$final_track_dir,
				"o:s"=>\$outdir,
				) or &USAGE;
&USAGE unless ($cuffdiff_dir and $final_track_dir and $outdir);
$final_track_dir=&ABSOLUTE_DIR($final_track_dir);
$cuffdiff_dir=&ABSOLUTE_DIR($cuffdiff_dir);
`mkdir $outdir` unless (-d $outdir);
$outdir=&ABSOLUTE_DIR($outdir);
# ------------------------------------------------------------------
# get saved genes and transcripts track list
# ------------------------------------------------------------------
#KB742382.1      XLOC_000001     ENSAPLG00000016114      5801    32652   +       TCONS_00000003  ENSAPLT00000016844      17476   32485   2005
#KB742382.1      XLOC_000003     Duck_newGene_1  121709  128447  +       TCONS_00000010  Duck_newGene_1.2        121709  128447  1092
my $raw_track_file=`ls $final_track_dir/*.final_tracking.list`;
my $new_track_file=`ls $final_track_dir/*.newGene_filtered_final_tracking.list`;
chomp($raw_track_file);
chomp($new_track_file);
my %track;

open (TRACK,"$raw_track_file") or die $!;

while (<TRACK>) {
	chomp;
	my @col=split /\t/;
	next if ($col[7]=~/CUFF\./);
	$track{$col[6]}{'gene_id'}=$col[2];
	$track{$col[6]}{'stand'}=$col[5];
	$track{$col[6]}{'iso_id'}=$col[7];
	$track{$col[6]}{'length'}=$col[10];
	$track{$col[6]}{'locus'}="$col[0]:$col[3]-$col[4]";
}

close TRACK;
open (TRACK,"$new_track_file") or die $!;

while (<TRACK>) {
	chomp;
	my @col=split /\t/;
	next if ($col[7]=~/CUFF\./);
	$track{$col[6]}{'gene_id'}=$col[2];
	$track{$col[6]}{'stand'}=$col[5];
	$track{$col[6]}{'iso_id'}=$col[7];
	$track{$col[6]}{'length'}=$col[10];
	$track{$col[6]}{'locus'}="$col[0]:$col[3]-$col[4]";
}

close TRACK;

# ------------------------------------------------------------------
# abstract data from isoforms.fpkm_tracking and isoforms.count_tracking
# ------------------------------------------------------------------
#tracking_id     class_code      nearest_ref_id  gene_id gene_short_name tss_id  locus   length  coverage        T1_FPKM T1_conf_lo      T1_conf_hi      T1_status       T2_FPKM
#TCONS_00000003  =       ENSAPLT00000016844      XLOC_000001     ENSAPLG00000016114      TSS3    KB742382.1:5800-32652   2005    -       0.743331        0       3.85708 OK
my $fpkm_file="$cuffdiff_dir/isoforms.fpkm_tracking";
my $count_file="$cuffdiff_dir/isoforms.count_tracking";
my @sample;

open (FPKM,$fpkm_file) or die $!;
chomp(my $head=<FPKM>);
my @attributes_fpkm=(split /\t/,$head);
for (my $i=9;$i<@attributes_fpkm;$i+=4) {
	my ($sample)=$attributes_fpkm[$i]=~/(\S+)_FPKM$/;
	push @sample,$sample if ($sample ne 'NULL');
}

while (<FPKM>) {
	if (exists $track{(split /\t/)[0]}) {
		chomp;
		my @col=split /\t/;

		for (my $i=0;$i<(@col-9)/4;$i++) {
			$track{$col[0]}{$sample[$i]."_FPKM"}=$col[9+4*$i] if (defined $sample[$i]);
		}
	}
}

close FPKM;

open (COUNT,$count_file) or dir $!;
<COUNT>;

while (<COUNT>) {
	if (exists $track{(split /\t/)[0]}) {
		chomp;
		my @col=split /\t/;

		for (my $i=0;$i<(@col-1)/5;$i++) {
			$track{$col[0]}{$sample[$i]."_count"}=$col[1+5*$i] if (defined $sample[$i]) ;
		}
	}
}

close COUNT;

# ------------------------------------------------------------------
# output 
# ------------------------------------------------------------------
my $fpkm_counts;

if (@sample>0) {
	for (my $i=0;$i<@sample;$i++) {
		$fpkm_counts.="$sample[$i]\_FPKM\t$sample[$i]\_count\t";
	}

	$fpkm_counts=~s/\t$//;
}

open (LOG,">$outdir/err.log") or die $!;
open (OUT,">$outdir/AllSample.isoforms_expression.xls") or die $!;
print OUT "#transcript_id\tlength\tstrand\tgene_id\tlocus\t$fpkm_counts\n";

for my $i (sort {$track{$a}{'gene_id'} cmp $track{$b}{'gene_id'}} keys %track) {
	my $basic_inf=join "\t",($track{$i}{'iso_id'},$track{$i}{'length'},$track{$i}{'stand'},$track{$i}{'gene_id'},$track{$i}{'locus'});
	my $print=$basic_inf;

	for (my $s=0;$s<@sample;$s++) {
		if (!defined $track{$i}{$sample[$s]."_FPKM"} || !defined $track{$i}{$sample[$s]."_count"}) {
			delete $track{$i};
			print LOG "WARNING: $i in sample $sample[$s] can't get expression enrichment value.\n";
			next;
		}
		my $sample_inf=join "\t",($track{$i}{$sample[$s]."_FPKM"},$track{$i}{$sample[$s]."_count"});
		$print=$print."\t".$sample_inf;
	}

	print OUT "$print\n" if (exists $track{$i});
}

close OUT;
close LOG;

for (my $s=0;$s<@sample;$s++) {

	open (EXP,">$outdir/$sample[$s].isoExpression.xls") or die $!;
	print EXP "#transcript_id\tlength\tstrand\tgene_id\tlocus\tFPKM\tcount\n";

	for my $i (sort {$track{$a}{'gene_id'} cmp $track{$b}{'gene_id'}} keys %track) {
		my $basic_inf=join "\t",($track{$i}{'iso_id'},$track{$i}{'length'},$track{$i}{'stand'},$track{$i}{'gene_id'},$track{$i}{'locus'});
		my $print=$basic_inf;

		my $sample_inf=join "\t",($track{$i}{$sample[$s]."_FPKM"},$track{$i}{$sample[$s]."_count"});
		$print=$print."\t".$sample_inf;
		print EXP "$print\n";
	}

	close EXP;
}



#######################################################################################
print STDOUT "\nDone. Total elapsed time : ",time()-$BEGIN_TIME,"s\n";
#######################################################################################

# ------------------------------------------------------------------
# sub function
# ------------------------------------------------------------------
################################################################################################################

sub ABSOLUTE_DIR{ #$pavfile=&ABSOLUTE_DIR($pavfile);
	my $cur_dir=`pwd`;chomp($cur_dir);
	my ($in)=@_;
	my $return="";
	if(-f $in){
		my $dir=dirname($in);
		my $file=basename($in);
		chdir $dir;$dir=`pwd`;chomp $dir;
		$return="$dir/$file";
	}elsif(-d $in){
		chdir $in;$return=`pwd`;chomp $return;
	}else{
		warn "Warning just for file and dir\n";
		exit;
	}
	chdir $cur_dir;
	return $return;
}

################################################################################################################

sub GetTime {
	my ($sec, $min, $hour, $day, $mon, $year, $wday, $yday, $isdst)=localtime(time());
	return sprintf("%4d-%02d-%02d %02d:%02d:%02d", $year+1900, $mon+1, $day, $hour, $min, $sec);
}

################################################################################################################
sub USAGE {
	my $usage=<<"USAGE";
 ProgramName:	$0
     Version:	$version
     Contact:	Simon Young <yangxh\@biomarker.com.cn> 
Program Date:	2014.04.15
      Modify:	
 Description:	This program is used to abstact transcripts expression enrichment per sample.
       Usage:
		Options:
		-f <str>	input directory,final_track directory,forced

		-c <str>	input directory,Cuffdiff directory,forced

		-o <str>	output directory,forced

		-h		help

        Example:
		perl $0 -f geneExpression/final_track/ -c Tophat_Cufflinks/Cuffdiff/ -o isoExpression/

USAGE
	print $usage;
	exit;
}
