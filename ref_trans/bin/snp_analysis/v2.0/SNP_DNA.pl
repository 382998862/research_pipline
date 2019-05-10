#!/usr/bin/perl
use strict;
use warnings;
use Getopt::Long;
use Data::Dumper;
use FindBin qw($Bin $Script);
use File::Basename qw(basename dirname);
use Cwd qw(abs_path);
use newPerlBase;
my %config=%{readconf("$Bin/../../../config/db_file.cfg")};

my $BEGIN_TIME=time();
my $version="1.1.0";
my @Original_ARGV=@ARGV;

# ==============================================================
# Get Options
# ==============================================================
my ($cfg,$dOut,$tophat,$step,$gff,$qphred);

GetOptions(
				"help|?" =>\&USAGE,
				"cfg:s"=>\$cfg,
				"tophat:s"=>\$tophat,
				"gff:s"=>\$gff,
				"od:s"=>\$dOut,
				"step:s"=>\$step,
				"qphred:s"=>\$qphred,
				) or &USAGE;
&USAGE unless ($cfg and $tophat);
#===============================================================
# Default optional value 
#===============================================================
###### software dir ###########
#my $_mapping_	=	"$Bin/STAR_align/v1.1/STAR_align.pl";
my $_SNP_		=	"$Bin/bin/GATK_calling_Trans_Ref.pl";
$dOut||="./SNP_Analysis ";
$step||=1;
$qphred||=33;
###### global value ###########
my %para;
my %sample;
my $cmd;
my $notename = `hostname`; chomp $notename;
$dOut=abs_path($dOut);
$tophat=abs_path($tophat);
#===============================================================
# pipeline
#===============================================================
log_current_time("$Script start...");
&read_config($cfg);
&MKDIR("$dOut");
my $genome=basename($para{'Ref_seq'});

=c
&MKDIR("$dOut/STAR");
$cmd = "perl $_mapping_ -cfg $cfg -od $dOut/STAR >/dev/null ";
log_current_time("STAR alignment start...");
log_current_time("CMD: $cmd");
system $cmd;
log_current_time("STAR alignment done.");
=cut
my $Rscript = $config{Rscript};

if ($step==1) {
 
	&MKDIR("$dOut/aligndir");
	&MKDIR("$dOut/work_sh");
	my @bam=glob("$tophat/*/accepted_hits.bam");

	open (OUT ,">$dOut/work_sh/bam_process.sh") or die $!;
	foreach my $bam (@bam) {
		print OUT "perl $Bin/bin/bin/bam_process.pl -od $dOut/aligndir -bam  $bam \n";
	}
	close (OUT);

	&qsubOrDie("$dOut/work_sh/bam_process.sh",$para{'Queue_type1'},"30",$para{'gatk_memory'});
	&qsubCheck("$dOut/work_sh/bam_process.sh");
	system "perl $Bin/bin/bin/sort_fa.pl -fa $para{'Ref_seq'} -od $dOut/aligndir ";
	$step++;
}         
if ($step==2) {
	$cmd = "perl $_SNP_ -ref $dOut/aligndir/$genome.fa -aligndir $dOut/aligndir -ploidy $para{'ploidy'} -win $para{'window'} -clu $para{'cluster'} -QD $para{'QD'} -FS $para{'FS'} -od $dOut/SNP -doRecal $para{'Recal'} -doIndel $para{'ReAlignIndel'} -queue $para{'Queue_type1'} -vf $para{'gatk_memory'} ";
    $cmd.= " --qphred " if ($qphred==64);
	$cmd.= " >/dev/null ";
	log_current_time("GATK calling start...");
	log_current_time("CMD: $cmd");
	system $cmd;
	log_current_time("GATK calling done.");

	##SNP sites stat
	log_current_time("pairwised SNP abstrct and SNP density plot start...");
	if (defined $gff ){
	                 while(! (-e $gff)){
                                sleep(600);
                        }
		        $cmd=`cat $para{'Ref_ann'} $gff > $dOut/Integrated.gff` ;
		        $cmd="perl $Bin/util/get_gene_fa.pl $dOut/SNP/Genome/$genome.fa $gff $dOut/New_gene";
		        print "$cmd\n";
		        system $cmd;
	}
	else {
		$cmd=`cp $para{'Ref_ann'}  $dOut/Integrated.gff` ;
	}

	$cmd="perl $Bin/util/get_gene_fa.pl $dOut/SNP/Genome/$genome.fa $dOut/Integrated.gff $dOut/All_gene";
	print "$cmd\n";
	system $cmd;

	$cmd="perl $Bin/SNP_indel_anno/snp_indel_anno.pl -id $dOut -r $dOut/SNP/Genome/$genome.fa -s $para{'Project_key'} -queue $para{'Queue_type1'} >/dev/null ";
	print "$cmd\n";
	system $cmd;

	$cmd="perl $Bin/util/SNP_stat.pl -snp $dOut/stat/final.snp.anno.gatk.all.list  -gff $dOut/Integrated.gff -od $dOut/stat/ >/dev/null ";
	print "$cmd\n";
	system $cmd;

	$cmd="$Rscript $Bin/bin/bin/final_SNP_type.r  $dOut/stat/All.snp_type.stat  $dOut/stat/ >/dev/null ";
	print "$cmd\n";
	system $cmd;

	#$cmd="perl $Bin/util/pairwised_SNP_abstrct_density_plot.pl --ref $dOut/aligndir/$genome.fa --snp $dOut/stat/final.snp.anno.gatk.all.list --od $dOut/ >/dev/null ";
	#print "$cmd\n";
	#system $cmd;
	log_current_time("pairwised SNP abstrct and SNP density plot done.");
}

#######################################################################################
my $elapsed_time = (time()-$BEGIN_TIME).'s';
log_current_time("$Script done. elapsed time: $elapsed_time.");
####################################################################################################
sub read_config {#
	my($cfg)=@_;
	open (IN,"$cfg") || die "$!";
	while (<IN>) {
		chomp;
		s/\r$//;
		s/\s+$//;
		next if (/\#/ || /^$/);
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
}

####################################################################################################
sub log_current_time {
     # get parameter
     my ($info) = @_;

     # get current time with string
     my $curr_time = date_time_format(localtime(time()));

     # print info with time
     print "[$curr_time] $info\n";
}

sub check_file {#
	my ($file) = @_;
	if (-e $file) {
		return 1;
	}else{
		return 0;
	}
}
sub Check_qsub_error {#
	# Check The qsub process if error happend 
	my $sh=shift;
	my @Check_file=glob "$sh*.qsub/*.Check";
	my @sh_file=glob "$sh*.qsub/*.sh";

	if ($#sh_file!=$#Check_file) {
		print "Their Some Error Happend in $sh qsub, Please Check..\n";
		die;
	}
	else {
		print "$sh qsub is Done!\n";
	}
}

####################################################################################################
sub date_time_format {
    my ($sec, $min, $hour, $day, $mon, $year, $wday, $yday, $isdst)=localtime(time());
    return sprintf("%4d-%02d-%02d %02d:%02d:%02d", $year+1900, $mon+1, $day, $hour, $min, $sec);
}

sub MKDIR
{ # &MKDIR($out_dir);
	my ($dir)=@_;
	mkdir($dir) if(!-d $dir);
}
sub USAGE {#
	my $usage=<<"USAGE";
    Program: $0
    Version: $version
    Contact: Shi Tongwei <shitw\@biomarker.com.cn> 
Discription:
      Usage:
        Options:
        -cfg    <file>  required, ref_trans.detail.cfg
        -od     <path>  optional, directory where output file produced, default [./SNP_Trans]
        -tophat     <path>  required, directory where the tophat result						(template: ./Basic_Analysis/Tophat_Cufflinks/Tophat)
        -step		<num>	optional	run from this step,default [1]
						1 :  Preparation bam and genome fa;
						2 :  SNP with GATK
        -gff  <file>  new_gene.gff ,optinal
        -help           help

USAGE
	print $usage;
	exit;
}
