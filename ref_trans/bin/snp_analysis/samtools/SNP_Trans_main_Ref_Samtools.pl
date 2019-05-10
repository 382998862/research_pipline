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
my ($cfg,$dOut,$tophat,$step,$gff);

GetOptions(
				"help|?" =>\&USAGE,
				"cfg:s"=>\$cfg,
				"tophat:s"=>\$tophat,
				"gff:s"=>\$gff,
				"od:s"=>\$dOut,
				"step:s"=>\$step,
				) or &USAGE;
&USAGE unless ($cfg and $tophat);
#===============================================================
# Default optional value 
#===============================================================
###### software dir ###########
#my $_mapping_	=	"$Bin/STAR_align/v1.1/STAR_align.pl";
my $_SNP_		=	"$Bin/bin/Samtools_calling_Trans_Ref.pl";
$dOut||="./SNP_Analysis ";
$step||=1;

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
$para{SNP_C}||= 50;
$para{SNP_D}||= "5-100";
$para{SNP_Q}||= 20;
$para{SNP_M}||= 5;
if ($step==1) {
 
	$cmd = "perl $_SNP_ -ref $para{Ref_seq} -tophatDir $tophat -gff $para{Ref_ann} -d $dOut -snp_C $para{SNP_C} -snp_D $para{SNP_D} -snp_Q $para{SNP_Q} -snp_M $para{SNP_M}";
    log_current_time("Samtool calling start...");
	log_current_time("CMD: $cmd");
	system $cmd;
	log_current_time("Samtools calling done.");
	$step++;
}         
if ($step==2) {

	##SNP sites stat
	log_current_time("pairwised SNP abstrct and SNP density plot start...");
	if (defined $gff ){
	                 while(! (-e $gff)){
                                sleep(600);
                        }
		        $cmd=`cat $para{'Ref_ann'} $gff > $dOut/Integrated.gff` ;
                &MKDIR("$dOut/aligndir");
                
                system "perl $Bin/bin/bin/sort_fa.pl -fa $para{'Ref_seq'} -od $dOut/aligndir ";
		        $cmd="perl $Bin/util/get_gene_fa.pl $dOut/aligndir/$genome.fa $gff $dOut/New_gene";
		        print "$cmd\n";
		        system $cmd;
	}
	else {
		$cmd=`cp $para{'Ref_ann'}  $dOut/Integrated.gff` ;
	}
##
	$cmd="perl $Bin/util/get_gene_fa.pl $dOut/aligndir/$genome.fa $dOut/Integrated.gff $dOut/All_gene";
	print "$cmd\n";
	system $cmd;
    my @bamfile=glob("$tophat/*/*sorted.bam");
    die "no alignment bam file!!!" if(@bamfile<1);
    foreach my $bamfile(@bamfile){
        my $sample=basename(dirname($bamfile));
        `ln -s $bamfile $dOut/aligndir/$sample.bam`;
    }
    &MKDIR("$dOut/data");
    &MKDIR("$dOut/SNP");
	$cmd="perl $Bin/SNP_indel_anno/snp_indel_anno.pl -id $dOut -r $dOut/aligndir/$genome.fa -s $para{'Project_key'} -queue $para{'Queue_type2'}  ";
	print "$cmd\n";
	system $cmd;

	$cmd="perl $Bin/util/SNP_stat.pl -snp $dOut/stat/final.snp.anno.Samtools.all.list  -gff $dOut/Integrated.gff -od $dOut/stat/  ";
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