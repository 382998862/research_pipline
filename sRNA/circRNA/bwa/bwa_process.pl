#!/usr/bin/perl -w
#use strict;
use warnings;
use Getopt::Long;
use Data::Dumper;
use FindBin qw($Bin $Script);
use File::Basename qw(basename dirname);
use newPerlBase;
my $BEGIN_TIME=time();
my $path = $Bin;
my %config=%{readconf("$Bin/../project.cfg")};
my $Title=$config{Title};												#流程的名称，必填
my $version=$config{version};
my @Original_ARGV=@ARGV;
#######################################################################################
# ------------------------------------------------------------------
# GetOptions, Check format, Output info.
# ------------------------------------------------------------------
my ($cfg1, $cfg2, $od, $step, $log, $oneStepOnly);
GetOptions(
				"help|?" =>\&USAGE,
				"cfg1:s"  =>\$cfg1,
				"cfg2:s"  =>\$cfg2,
				"od:s"   =>\$od,
				"s:s"    =>\$step,
				"oneStepOnly:s"=>\$oneStepOnly,
				) or &USAGE;
&USAGE unless ($cfg1 and $cfg2 and $od) ;
################################

my $notename = `hostname`; chomp $notename;

$cfg1 = &ABSOLUTE_DIR($cfg1);
$cfg2 = &ABSOLUTE_DIR($cfg2);
&MKDIR($od);
$od  = &ABSOLUTE_DIR($od);     &MKDIR("$od/work_sh");

$step = $step || 1;
my $thread ||= 10;
#
# log file
#
my $startTime = GetTime();
my $user = `whoami`;  chomp $user;
my $workDir = `pwd`;  chomp $workDir;
my $task = "perl $Bin/$Script ".join("\t", @Original_ARGV);

open ($log,">", "$od/Bwa.".time().".log") or die $!;

print $log "######################################\n";
print $log "$startTime\n";
print $log "user:	$user\n";
print $log "work directory:	$workDir\n";
print $log "job:	$task\n";
print $log "######################################\n\n";

print $log "data config file:  $cfg1\n";
print $log "detail config file:  $cfg2\n";
print $log "output dir:  $od\n";
print $log "start from step: $step\n\n";

#==================================================================
# bins
#==================================================================


my $GFFREAD_BIN     = $config{gffread};     # 2014-12-17 ~
my $bwa_dir = $config{bwa};
my $CUFFQUANT_BIN   = $config{cuffquant};   # 2015-08-04
my $CUFFNORM_BIN    = $config{cuffnorm};    # 2015-08-04
my $samtools = $config{samtools};
#==================================================================
# load config file
#==================================================================

my %total_read;
my %para;
my %sample;

open (IN,"cat $cfg2 $cfg1|") || die "$!\n";
while (<IN>) {
	chomp;
	s/\r$//;s/^\s+//;s/\s+$//;
	next if (/^\#/ || /^$/);

	my @tmp=split /\s+/,$_;
	if ($tmp[0] eq "Sample") {
		my $fq1=<IN>;  chomp $fq1;
		my $fq2=<IN>;  chomp $fq2;
		my @fq_1=split /\s+/,$fq1;
		$sample{$tmp[1]}{FQ1}=$fq_1[1];
		if (!-f "$od/totalRead.stat.xls") {
			my $total_line=`less -S $sample{$tmp[1]}{FQ1} | wc -l `;chomp $total_line;
			$total_read{$tmp[1]}=$total_line/2;
		}
		my @fq_2=split /\s+/,$fq2;
		$sample{$tmp[1]}{FQ2}=$fq_2[1];
	}
	$para{$tmp[0]}=$tmp[1];
}
close IN;

if (!-f "$od/totalRead.stat.xls") {
	open OUT,">$od/totalRead.stat.xls" || die;
	foreach my $sam (sort keys %total_read) {
		print OUT "$sam\t$total_read{$sam}\n";
	}
	close OUT;
}else {
	open IN,"$od/totalRead.stat.xls" || die;
	while (<IN>) {
		chomp;
		next if (/^$/);
		my @tmp=split /\s+/,$_;
		$total_read{$tmp[0]}=$tmp[1];
	}
}

#==================================================================
# pipeline
#==================================================================

#######################################
#
# step 1: check and build bowtie index
#
######

my $genome;
my $idx_prefix;
my $gtf;
my $gff;

if ($step!=1) {
    $genome = basename($para{Ref_seq});
    $idx_prefix = basename($para{Ref_seq});
    $idx_prefix =~s/.fa$//;
    $gff = basename($para{Ref_ann});
    $gtf = basename($para{Ref_ann});
    $gtf =~s/\.gff3?$//i;
    $gtf.= ".gtf";
    print "$gtf\n";
}
if ($step==1) {
	print STDOUT "=== check and build bwa index ===\n";
	print $log  "=== check and build bwa index ===\n";

	mkdir("$od/Ref_Genome") if(!-d "$od/Ref_Genome");
	chdir "$od/Ref_Genome";
	$genome = basename($para{Ref_seq});
	$idx_prefix = $para{Ref_seq};

	if (-e "$idx_prefix.sa" and -e "$idx_prefix.pac" and -e "$idx_prefix.ann" and -e "$idx_prefix.amb" and -e "$idx_prefix.bwt" ){
		system "ln -s $idx_prefix* ./";
	}
	elsif(!-e "bwa_index.finish") {		
		system "ln -s $idx_prefix ./ ";
		&run_or_die("$bwa_dir index -a bwtsw $genome");
		system "touch bwa_index.finish";
    	}
	################################## gff2gtf
    $gff = basename($para{Ref_ann});
    $gtf = basename($para{Ref_ann});
    $gtf =~s/\.gff3?$//i;
    $gtf.= ".gtf";
	system "ln -s $para{Ref_ann} ./" unless (-f $gff);
	system "$GFFREAD_BIN $gff -T -o $gtf";
	chdir "../";

	$step++ unless ($oneStepOnly) ;
	print STDOUT "\n";
	print $log "\n";

}

####################################
#
# step 2: Align the RNA-seq read to genome using Bwa
#
#########

if ($step==2) {
	print STDOUT "=== Align the RNA-seq read to genome using Bwa  ===\n";
	print STDOUT "Bwa mapping shell file: $od/work_sh/Bwa.sh\n";
	print $log "=== Align the RNA-seq read to genome using Bwa  ===\n";
	print $log "Bwa mapping shell file: $od/work_sh/Bwa.sh\n";

	#
	# write shell
	#
	open OUT,">$od/work_sh/Bwa.sh" || die;
	open SAM,">$od/work_sh/samtools.sh" || die;
	&MKDIR("$od/Bwa");
	foreach my $sam (sort keys %sample) {
		&MKDIR("$od/Bwa/$sam");
		print OUT "cd $od/Bwa/$sam && ";
		print OUT "$bwa_dir mem -T 19 -t $thread $od/Ref_Genome/$genome $sample{$sam}{FQ1} $sample{$sam}{FQ2} > $sam.sam \n";
		print SAM " $samtools view -bS $od/Bwa/$sam/$sam.sam -o $od/Bwa/$sam/$sam.bam && $samtools  sort  $od/Bwa/$sam/$sam.bam  $od/Bwa/$sam/$sam.sort\n";
	}
	close OUT;
	close SAM;

	$para{Memory}||="15G";
	$para{Queue_type}||="medical.q";

	&Cut_shell_qsub("$od/work_sh/Bwa.sh",18,"$para{Memory}","$para{Queue_type}");
	&Check_qsub_error("$od/work_sh/Bwa.sh");

	&Cut_shell_qsub("$od/work_sh/samtools.sh",18,"$para{Memory}","$para{Queue_type}");
	&Check_qsub_error("$od/work_sh/samtools.sh");
	$step++ unless ($oneStepOnly) ;
	print STDOUT "\n";
	print $log "\n";
}

#####################################
#
# step 3: Statistic bam files
#
########

if ($step==3) {
	print STDOUT "=== Statistic bam files  ===\n";
	print STDOUT "shell file: $od/work_sh/Bwa_bam_stat.sh\n";
	print STDOUT "shell file: $od/work_sh/genome_bam2depth.sh\n";
	print STDOUT "shell file: $od/work_sh/genome_Checkgraph.sh\n";
	print $log "=== Statistic bam files  ===\n";
	print $log "shell file: $od/work_sh/Bwa_bam_stat.sh\n";
	print $log "shell file: $od/work_sh/genome_bam2depth.sh\n";
	print $log "shell file: $od/work_sh/genome_Checkgraph.sh\n";
	$para{Queue_type}||="medical.q";
	open OUT1, ">$od/work_sh/Bwa_bam_stat.sh"   || die;
	open OUT2, ">$od/work_sh/genome_bam2depth.sh"  || die;
	open OUT3, ">$od/work_sh/genome_Checkgraph.sh" || die;
	open OUT4, ">$od/work_sh/plot_ReadDensity.sh"  || die;
	&MKDIR("$od/Map_Stat");
	my @str_type_stat;
	foreach my $sam (sort keys %sample) {
		push @str_type_stat, "$od/Map_Stat/$sam.type.stat";
		print OUT1 "perl $Bin/bin/bam2map_stat.pl -i $sam -bam $od/Bwa/$sam/$sam.bam -totalRead $total_read{$sam} -od $od/Map_Stat\n";
		print OUT2 "$samtools depth $od/Bwa/$sam/$sam.sort.bam >$od/Bwa/$sam/$sam.sort.bam.depth && $samtools faidx $od/Ref_Genome/$genome\n";
		print OUT3 "perl $Bin/bin/get_percent_of_exon_intro_inter_by_gff_v1.4.pl -gff $od/Ref_Genome/$gff -i $od/Bwa/$sam/$sam.sort.bam.depth -od $od/Map_Stat -index $sam &&\n";
		print OUT4 "perl $Bin/bin/plotReadDensity2.pl -i $od/Bwa/$sam/$sam.sort.bam -a $od/Ref_Genome/$genome.fai -f bam -o $od/Map_Stat/ -k $sam  ";
		print OUT4 " -medical $para{medical} "	if(exists $para{medical});
		print OUT4 "\n";

	}
	close OUT1;
	close OUT2;
	close OUT3;
	close OUT4;

	`$config{qsub} --queue $para{Queue_type} --resource vf=80G --reqsub --independent $od/work_sh/Bwa_bam_stat.sh`;
	&Check_qsub_error("$od/work_sh/Bwa_bam_stat.sh");

	&Cut_shell_qsub("$od/work_sh/genome_bam2depth.sh",  18, "50G", $config{queue});
	&Check_qsub_error("$od/work_sh/genome_bam2depth.sh");

	&Cut_shell_qsub("$od/work_sh/genome_Checkgraph.sh", 18, "50G",  $config{queue});
	&Check_qsub_error("$od/work_sh/genome_Checkgraph.sh");

	&Cut_shell_qsub("$od/work_sh/plot_ReadDensity.sh",  18, "50G",  $config{queue});
	&Check_qsub_error("$od/work_sh/plot_ReadDensity.sh");
	
	my $str_type_stat = join " ", @str_type_stat;
	`perl $Bin/bin/map_total_type.pl -i $str_type_stat -o $od/Map_Stat/Total_Type.png`;
	print "perl $Bin/bin/map_total_type.pl -i $str_type_stat -o $od/Map_Stat/Total_Type.png\n";
	$step++ unless ($oneStepOnly) ;
	print STDOUT "\n";
	print $log "\n";
}
#######################################################################################
print STDOUT "\nDone. Total elapsed time : ",time()-$BEGIN_TIME,"s\n";
print $log  "\nDone. Total elapsed time : ",time()-$BEGIN_TIME,"s\n";
#######################################################################################
close($log);
#==================================================================
# subs
#==================================================================
sub LOAD_PARA {
	my $para_file= shift;
	my $para= shift;

	my $error_status = 0;
	open IN,$para_file || die "fail open: $para_file";
	while (<IN>) {
		chomp;
		s/^\s+//;s/\s+$//;s/\r$//;
		next if (/^$/ or /^\#/) ;
		my ($para_key,$para_value) = split(/\s+/,$_);
		$para->{$para_key} = $para_value;
		if (!defined $para_value) {
			warn "Non-exist: $para_key\n";
			$error_status = 1;
		}
	}
	close IN;
	die "\nExit due to error Undefined parameter\n" if ($error_status) ;
}

sub Cut_shell_qsub {#Cut shell for qsub 1000 line one file
	# &Cut_shell_qsub($shell,$cpu,$vf,$queue);
	my $shell = shift;
	my $cpu = shift;
	my $vf = shift;
	my $queue = shift;

	my $line = `less -S $shell |wc -l `;
	if ($line<=1000) {
		system "$config{qsub} --queue $queue --maxproc $cpu --resource vf=$vf --reqsub --independent $shell";
	}
	if ($line>=1000) {
		my @div=glob "$shell.div*";
		foreach (@div) {
			if (-e $_) {
				system "rm $_";
			}
		}
		@div=();
		my $div_index=1;
		my $line_num=0;
		open IN,"$shell" || die;
		while (<IN>) {
			chomp;
			open OUT,">>$shell.div.$div_index" || die;
			if ($line_num<1000) {
				print OUT "$_\n";
				$line_num++;
			}
			else {
				print OUT "$_\n";
				$div_index++;
				$line_num=0;
				close OUT;
			}
		}
		if ($line_num!=0) {
			close OUT;
		}
		@div=glob "$shell.div*";
		foreach my $div_file (@div) {
			system "$config{qsub} --queue $queue --maxproc $cpu --resource vf=$vf --reqsub --independent $div_file";
		}
	}
}

sub Check_qsub_error {#
	# Check The qsub process if error happend
	my $sh=shift;
	my @qsub_dir = glob "$sh*.qsub";
	print "@qsub_dir\n";
	if(@qsub_dir=0){print "Your qsub has not been run, Please Check..\n";}
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

sub GetTMR {#
	#Get Total Mapped Reads from file line which contain string "Mapped Reads\s"
	my $fStat=shift;
	open (IN,"<",$fStat) or die $!;
	while (<IN>) {
		if (/^Mapped Reads\s(\d+)/) {
			close (IN) ;
			return $1;
		}
	}
	close (IN) ;
	die "Error Reads Stat file.\n";
}

sub MKDIR
{ # &MKDIR($out_dir);
	my ($dir)=@_;
	rmdir($dir) if(-d $dir);
	mkdir($dir) if(!-d $dir);
}

sub ABSOLUTE_DIR
{ #$pavfile=&ABSOLUTE_DIR($pavfile);
	my $cur_dir=`pwd`;chomp($cur_dir);
	my ($in)=@_;
	my $return="";

	if(-f $in)
	{
		my $dir=dirname($in);
		my $file=basename($in);
		chdir $dir;$dir=`pwd`;chomp $dir;
		$return="$dir/$file";
	}
	elsif(-d $in)
	{
		chdir $in;$return=`pwd`;chomp $return;
	}
	else
	{
		warn "Warning just for file and dir in [sub ABSOLUTE_DIR]\n";
		exit;
	}

	chdir $cur_dir;
	return $return;
}

sub GetTime {
	my ($sec, $min, $hour, $day, $mon, $year, $wday, $yday, $isdst)=localtime(time());
	return sprintf("%4d-%02d-%02d %02d:%02d:%02d", $year+1900, $mon+1, $day, $hour, $min, $sec);
}
sub run_or_die()
{
        my ($cmd) = @_ ;
        &show_log($cmd);
        my $flag = system($cmd) ;
        if ($flag != 0){
                &show_log("Error: command fail: $cmd");
                exit(1);
        }
        &show_log("done.");
        return ;
}
sub show_log()
{
        my ($txt) = @_ ;
        my $time = time();
        my($sec, $min, $hour, $day, $mon, $year, $wday, $yday, $isdst) = localtime($time);
        $wday = $yday = $isdst = 0;
        my $Time=sprintf("%4d-%02d-%02d %02d:%02d:%02d", $year+1900, $mon+1, $day, $hour, $min, $sec);
        print "$Time:\t$txt\n" ;
}


sub USAGE {#
	my $usage=<<"USAGE";
Program: Bwa_Analysis Procedure
Version: $version
Contact: Mengmeng Song <songmm\@biomarker.com.cn>

Description:
	This program is a Procedure deal with Multiple Samples RNA_Analysis
	The program will Generate SAM file

Usage:
    -cfg1             data config, rawdata & refseq path
    -cfg2             detail config, analysis parameters
	-od               output dir            must be given
	-s                step of the program   option,default 1;
                            1  Check the index of Genome & build index for alignment
                            2  run Bwa analysis
                            3  stat Bwa sam out with Genome
                            4  Cat All Samples SAM file
							5  run Cuffquant
                            6  run Cuffnorm
	-oneStepOnly      perform one of the pipeline steps above mentioned
USAGE
	print $usage;
	exit;
}
