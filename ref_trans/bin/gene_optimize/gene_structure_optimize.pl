#!/usr/bin/perl -w
use strict;

use Getopt::Long;
use Data::Dumper;
use FindBin        qw($Bin $Script);
use File::Basename qw(basename dirname);
use Cwd            'abs_path';

my $BEGIN_TIME=time();

# ------------------------------------------------------------------
my ($gtf, $gff, $od, $index, $step,$oneStepOnly);
GetOptions(
				"help|?"      =>\&USAGE,
				"gtf:s"       =>\$gtf,
				"gff:s"    =>\$gff,
				"od:s"        =>\$od,
				"index:s"     =>\$index,
				"s:s"         =>\$step,
				"oneStepOnly:s"=>\$oneStepOnly,
				) or &USAGE;
&USAGE unless ($gtf and  $gff  and $index and $od) ;
################################
my $notename=`hostname`;chomp $notename;
$gff=abs_path($gff);
$gtf = abs_path($gtf);
$od ||= "./";            mkdir $od unless (-d $od); 
$od  = abs_path($od);

my $log_dir      = "$od/track_list_log";     mkdir $log_dir      unless (-d $log_dir);
my $track_dir    = "$od/final_track";        mkdir $track_dir    unless (-d $track_dir);
my $sh_dir       = "$od/work_sh";            mkdir $sh_dir       unless (-d $sh_dir);

$step = $step || 1;
####################################### Make info for the program
if ($step==1) {
	open OUT,">$od/work_sh/Step1.sh" || die $!;
	print OUT "cd $log_dir  &&   perl $Bin/track2name.pl $gtf  $index  ";
	close OUT;
	system "sh $od/work_sh/Step1.sh";
	$step++ unless ($oneStepOnly) ;
}

if ($step == 2) {
	open OUT,">$od/work_sh/Step2.sh" || die $!;
	print OUT "cd $log_dir && perl $Bin/Track2name_format.pl $index $index \n\n";
	close OUT;
	system "sh $od/work_sh/Step2.sh";
	$step++ unless ($oneStepOnly) ;
}

if ($step == 3) {
	open OUT,">$od/work_sh/Step3.sh" || die $!;
	print OUT "cd $od/track_list_log  &&  perl $Bin/final_track_info.pl $gtf $index $od/final_track    ";
	close OUT;
	system "sh $od/work_sh/Step3.sh";
	$step++ unless ($oneStepOnly) ;
}
if ($step == 4) {
	open OUT,">$od/work_sh/Step4.sh" || die $!;
	print OUT "perl $Bin/structure_optimize.pl  $od/final_track/$index.final.gff $gff  $od/$index    ";
	close OUT;
	system "sh $od/work_sh/Step4.sh";
	$step++ unless ($oneStepOnly) ;
}

######################################### subs

sub USAGE {#
	my $usage=<<"USAGE";

Usage:
	-gtf              Compare.gtf fil_trackle       must be given
	-gff                ref_ann 
	-index            Species name          must be given
	-od               output dir            option,default ./;
	-s                step of the program   option,default 1;

USAGE
	print $usage;
	exit;
}
