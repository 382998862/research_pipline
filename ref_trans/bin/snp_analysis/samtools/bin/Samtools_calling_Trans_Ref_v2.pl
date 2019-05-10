#!/usr/bin/perl -w
use strict;

use Getopt::Long;
use FindBin qw($Bin $Script);
use File::Basename;
use Cwd 'abs_path';

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
my ($min_depth, $max_depth) = $snp_D =~ /(\d+)-(\d+)/;

#
# step 1: calling snp using samtools and bcftools 
#
if ($step == 1) 
{
	my @sample = glob("$tophatDir/*/accepted_hits.bam");
	
	# write shell 
	#
	mkdir "$dOut/aligndir" unless (-d "$dOut/aligndir");
	mkdir "$dOut/SNP" unless (-d "$dOut/SNP");
	mkdir "$dOut/work_sh" unless (-d "$dOut/work_sh") ;
	my $sh = "$dOut/work_sh/snp_calling.sh";
	open (SH, ">$sh") or die $!;
	open(FB,">$dOut/filelist.txt")||die $!;
	foreach my $bam (@sample) 
	{
		my $sample=basename(dirname($bam));
		`samtools sort $bam $dOut/aligndir/$sample.sort`;  
	   #` cp  $bam $dOut/aligndir/$sample.bam` unless (-f "$dOut/aligndir/$sample.bam");
	   print FB "$dOut/aligndir/$sample.sort.bam\n";
	
	}
	close FB;
	#my $cmd = "samtools mpileup -f $fRef -D -C $snp_C -g -s -u -b $dOut/aligndir/filelist.txt|$Bin/bin/bcftools/bcftools view -cegINv - >$dOut/SNP/All_sample.raw_snp.vcf";
	my $cmd = "samtools mpileup -Q 13 -q $snp_Q -d $max_depth -B -f $fRef -D -C $snp_C  -u -b $dOut/filelist.txt|$Bin/bin/bcftools/bcftools view -cegNv - >$dOut/SNP/All_sample.raw_snp.vcf";
	close (SH) ; 
	print "$cmd\n";
	system("$cmd");
	$step++;
}

#
# step 2: filter false positive snps 
#
if ($step == 2) 
{
	my $vcf="$dOut/SNP/All_sample.raw_snp.vcf";
	my $cmd = "perl $Bin/bin/SNPfilter_v2.pl -i $vcf -q $snp_Q -d $snp_D -r $snp_M -o $dOut/SNP";
	print "$cmd\n";
	system($cmd);

	$step++;
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

sub step_cmd_process {
    my ($cmd, $step_hash, $step_n, $sh_dir,$queue,$cpu,$vf) = @_;
    my $sh_file = "$sh_dir/Step$step_n.$step_hash->{$step_n}.sh";
    my $log_file = "$sh_file.log";
    my $flag = 0;
    my $start_time = time();
    &log_current_time("step$step_n. $step_hash->{$step_n} start ...");
    &log_current_time("CMD: $cmd");

    if (-e $sh_file) {
        system "cat $sh_file >> $sh_file.bak";
        open (SH, ">$sh_file") or die "$!: $sh_file\n";
        print SH "$cmd\n";
        close SH;
    } else {
        open (SH, ">$sh_file") or die "$!: $sh_file\n";
        print SH "$cmd\n";
        close SH;
    }
    system("sh $sh_file >$log_file 2>&1");
   # qsubOrDie("$sh_file",$queue,$cpu,$vf);	
    #qsubCheck("$sh_file"); 
}


sub USAGE {#
	my $usage=<<"USAGE";
Modified by Liuxs(liuxs\@biomarker.com) 2016.09.05
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
