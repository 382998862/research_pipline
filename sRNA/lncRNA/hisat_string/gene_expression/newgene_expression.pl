#!/usr/bin/perl -w
use strict;

use Getopt::Long;
use Data::Dumper;
use FindBin        qw($Bin $Script);
use File::Basename qw(basename dirname);
use Cwd            'abs_path';

my $BEGIN_TIME=time();

# ------------------------------------------------------------------
my ($gtf, $genome, $cuffnorm, $cufflinks, $od, $index, $step,$oneStepOnly,$singleexon,$queue);
GetOptions(
				"help|?"      =>\&USAGE,
				"gtf:s"       =>\$gtf,
				"genome:s"    =>\$genome,
				"cuffnorm:s"  =>\$cuffnorm,
				"cufflinks:s" =>\$cufflinks,
				"od:s"        =>\$od,
				"index:s"     =>\$index,
				"singleexon"=>\$singleexon,
				"s:s"         =>\$step,
				"queue:s"     =>\$queue,
				"oneStepOnly:s"=>\$oneStepOnly,
				) or &USAGE;
&USAGE unless ($gtf and  $genome  and $index and $od) ;
################################
my $notename=`hostname`;chomp $notename;
$genome=abs_path($genome);
$gtf = abs_path($gtf);
$od ||= "./";            mkdir $od unless (-d $od); 
$od  = abs_path($od);

my $log_dir      = "$od/track_list_log";     mkdir $log_dir      unless (-d $log_dir);
my $track_dir    = "$od/final_track";        mkdir $track_dir    unless (-d $track_dir);
my $new_gene_dir = "$od/new_gene";           mkdir $new_gene_dir unless (-d $new_gene_dir);
my $old_gene_dir = "$od/old_gene";           mkdir $old_gene_dir unless (-d $old_gene_dir);
my $sh_dir       = "$od/work_sh";            mkdir $sh_dir       unless (-d $sh_dir);

$step = $step || 1;
####################################### Make info for the program
if ($step==1) {
	open OUT,">$od/work_sh/Step1_get_NewGene_track.sh" || die $!;
	print OUT "perl $Bin/bin/get_track.pl -gtf  $gtf -index $index  -out $log_dir/$index.newGene.track.list.info \n";
	close OUT;
	system "sh /share/nas2/genome/bmksoft/tool/qsub_sge_plus/v1.0/qsub_sge.plus.sh --reqsub --independent --queue $queue $od/work_sh/Step1_get_NewGene_track.sh";
	$step++ unless ($oneStepOnly) ;
}

if ($step == 2) {
	open OUT,">$od/work_sh/Step2_newGene_format.sh" || die $!;
	print OUT "cd $od/track_list_log && perl $Bin/bin/newGene_format.pl -gtf  $gtf -track $log_dir/$index.newGene.track.list.info -out1 $od/final_track/$index.newGene_final_tracking.list -out2 $od/final_track/$index.newGene_final.gff \n";
	close OUT;
	system "sh /share/nas2/genome/bmksoft/tool/qsub_sge_plus/v1.0/qsub_sge.plus.sh --reqsub --independent --queue $queue $od/work_sh/Step2_newGene_format.sh";
	$step++ unless ($oneStepOnly) ;
}

if ($step == 3) {
	open OUT,">$od/work_sh/Step3_filter.newGene.sh" || die $!;
	print OUT "perl $Bin/bin/filter.new_gene.pl -s $index -g $genome -newdir $od/final_track/ ";
	if (defined $singleexon) {print OUT "-singleexon";}
	print OUT "&& perl $Bin/bin/get_cds.pl $genome $od/final_track/$index.newGene_final.filtered.gff $od/final_track/$index.newGene\n";
	close OUT;
	system "sh /share/nas2/genome/bmksoft/tool/qsub_sge_plus/v1.0/qsub_sge.plus.sh --reqsub --independent --queue $queue $od/work_sh/Step3_filter.newGene.sh";
	#$step++ unless ($oneStepOnly) ;
}

######################################### subs

sub USAGE {#
	my $usage=<<"USAGE";

Usage:
	-gtf              Compare.gtf fil_trackle       must be given
	-genome           reference genome      must be given

	-index            Species name          must be given
	-singleexon			do not filter the single exon newgene,
	-od               output dir            option,default ./;
	-s                step of the program   option,default 1;

USAGE
	print $usage;
	exit;
}
