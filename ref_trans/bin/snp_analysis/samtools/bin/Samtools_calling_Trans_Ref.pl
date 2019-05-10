#!/usr/bin/perl -w
use strict;

use Getopt::Long;
use FindBin qw($Bin $Script);
use File::Basename;
use Cwd 'abs_path';
use newPerlBase;
my %config=%{readconf("$Bin/../../../../config/db_file.cfg")};
my ($fRef, $tophatDir, $gff, $dOut, $snp_C, $snp_D, $snp_Q, $snp_M, $step);
GetOptions(
				"help|?" =>\&USAGE,

				"ref:s"       =>\$fRef,
				"tophatDir:s" =>\$tophatDir,
				"gff:s"       =>\$gff,
				"d:s"         =>\$dOut,

				"snp_C:s"=>\$snp_C,
				"snp_D:s"=>\$snp_D,
				"snp_Q:s"=>\$snp_Q,
				"snp_M:s"=>\$snp_M,

				"step:s"=>\$step,
				) or &USAGE;
&USAGE unless ($fRef and $tophatDir and $gff);

$dOut ||= "./";               mkdir $dOut unless (-d $dOut);
$dOut = abs_path($dOut);
$fRef = abs_path($fRef);

unless (-e "$fRef.fai") {
    print "WARNING: please make index for reference sequence by following command first:\nsamtools faidx $fRef\n\n";
    exit;
}

$snp_C ||= 50;
$snp_D ||= "5-100";
$snp_Q ||= 20;
$snp_M ||= 5;
$step  ||= 1;


#
# step 1: calling snp using samtools and bcftools 
#
if ($step == 1) 
{
	my @sample = glob("$tophatDir/*");

	#
	# write shell 
	#
	mkdir "$dOut/work_sh" unless (-d "$dOut/work_sh") ;
	my $sh = "$dOut/work_sh/snp_calling.sh";
	open (SH, ">$sh") or die $!;
	foreach my $sam (@sample) 
	{
		my $samBase = basename($sam);
		my $outDir  = $dOut."/".$samBase;
		mkdir $outDir unless (-d $outDir);
		
		my $bamfile = "$sam/$sam.HISAT_aln.sorted.bam";
		next unless (-f $bamfile);
		$bamfile    = abs_path($bamfile);
		
		my $cmd = "samtools mpileup -f $fRef -D -C $snp_C -g -s -u $bamfile|$Bin/bin/bcftools/bcftools view -cegNv - >$outDir/$samBase.raw_snp.vcf";
		#my $cmd = "samtools mpileup  -f $fRef -D -C $snp_C  -u -b $dOut/filelist.txt|$Bin/bin/bcftools/bcftools view -cegNv - >$dOut/SNP/All_sample.raw_snp.vcf";
		print SH $cmd, " &&\n";
	}
	close (SH) ;
	
	#
	# qsub
	#
	#&qsubOrDie("$sh","$config{Queue_type2}",12,"15G");
	&qsubOrDie("$sh","middle.q",12,"15G");
	&qsubCheck("$sh");
	#Cut_shell_qsub  ($sh, 12, "15G", "middle.q");
	#Check_qsub_error($sh);

	$step++;
}

#
# step 2: filter false positive snps 
#
if ($step == 2) 
{
	my @vcfFiles = glob("$dOut/*/*.raw_snp.vcf");
	mkdir "$dOut/work_sh" unless (-d "$dOut/work_sh");
	my $sh = "$dOut/work_sh/snp_filtering.sh";
	open (SH, ">", $sh) || die "Open $sh failed!\n";
	foreach my $vcf (@vcfFiles)
	{
		my ($sample)   = basename($vcf) =~ /(\S+)\.raw_snp\.vcf$/;
		my $vcfDir     = dirname($vcf);
		my $filterDir = "$vcfDir/filter";
		mkdir $filterDir unless (-d $filterDir);

		my $cmd = "perl $Bin/bin/SNPfilter.pl -i $vcf -key $sample -q $snp_Q -d $snp_D -r $snp_M -o $filterDir";
		print "$cmd\n";
		`$cmd`;
		print SH $cmd, " &&\n";
	}
	close SH;

	#
	# merge sample snp result 
	#
	## get snp files 
	my @snpFiles = glob("$dOut/*/filter/*.snp");
	die "too much or less snp files\n" if (@snpFiles != @vcfFiles);
	
	my $cmd = "perl $Bin/bin/merge_sample_snp2.pl -i @snpFiles -gff $gff -d $dOut"; 
	`$cmd`;

	$step++;
}

#
# step 3: statistic 
#
if ($step == 3) 
{

	my @statFile = glob("$dOut/*/filter/*.stat.xls");
	my $cmd      = "perl $Bin/bin/SNPstat.pl -i @statFile -s $dOut/sam.merge.snp -d $dOut";
	`$cmd`;
}

sub Cut_shell_qsub
{#Cut shell for qsub 1000 line one file
	my ($shell, $cpu, $vf, $queue) = @_;

	`ssh cluster qsub-sge.pl --queue $queue --maxproc $cpu --resource vf=$vf --independent --reqsub $shell`;
}

sub Check_qsub_error 
{# Check The qsub process if error happend 
	my ($sh) = @_;

	my @Check_file = glob "$sh*.qsub/*.Check";
	my @sh_file    = glob "$sh*.qsub/*.sh";
	if ($#sh_file != $#Check_file) 
	{
		print "Their Some Error Happend in $sh qsub, Please Check..\n";
		die;
	}
	else 
	{
		print "$sh qsub is Done!\n";
	}
}

sub USAGE {#
	my $usage=<<"USAGE";
Discription:
		Call SNP in species with genome references by RNA-seq .
Usage:
  -ref			<file>	Reference genome file, fasta format, forced
  -tophatDir		<dir>	Tophat alignment result, forced
  -gff			<file>	GFF file, forced
  -d			<str>	Directory where output file produced,optional [default ./]

  -snp_C		<int>	Parameter for adjusting mapQ; 0 to disable;C50 to reduce the effect of reads with excessive mismatches for BWA-short, optional [default 50]
  -snp_D		<str>	Depth range for filtering SNP, optional [default 5-100]
  -snp_Q		<int>	Minimum mapping quality value, optional [default 20]
  -snp_M		<int>	Minimum multi_SNP range, optional [default 5]
  
  -step			<int>	step, optional [default 1]

  -h		Help

USAGE
	print $usage;
	exit;
}
