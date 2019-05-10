#!/usr/bin/perl -w
use strict;
use warnings;
use Getopt::Long;
use Data::Dumper;
use FindBin qw($Bin $Script);
use File::Basename qw(basename dirname);
use threads;
use newPerlBase;
my $BEGIN_TIME=time();
my $version="1.0.0";
my %config=%{readconf("$Bin/../config/db_file.cfg")}; 

my @Original_ARGV=@ARGV;
#######################################################################################
# ------------------------------------------------------------------
# GetOptions, Check format, Output info.
# ------------------------------------------------------------------
my ($bam, $cfg2, $od, $step, $oneStepOnly);
GetOptions(
				"help|?" =>\&USAGE,
				"bam:s"  =>\$bam,
				"cfg2:s"  =>\$cfg2,
				"od:s"   =>\$od,
				"s:s"    =>\$step,
				"oneStepOnly:s"=>\$oneStepOnly,
				) or &USAGE;
&USAGE unless ($bam and $cfg2 and $od) ;
################################

my $notename = `hostname`; chomp $notename;

$bam = &ABSOLUTE_DIR($bam);
$cfg2 = &ABSOLUTE_DIR($cfg2);
&MKDIR($od) unless (-d $od);
$od  = &ABSOLUTE_DIR($od);   
&MKDIR("$od/work_sh") unless (-d "$od/work_sh");

$step = $step || 1;


#
# log file 
#
my $startTime = GetTime();
my $user = `whoami`;  chomp $user;
my $workDir = `pwd`;  chomp $workDir;
my $task = "perl $Bin/$Script ".join("\t", @Original_ARGV);


#==================================================================
# bins 
#==================================================================
my $GFFREAD		="$config{cufflink}/gffread";     
my $CUFFCOMPARE		="$config{cufflink}/cuffcompare"; 
my $GFFCOMPARE		="$Bin/stringtie/gffcompare";
my $STRINGTIE		="$Bin/stringtie/stringtie";


my %para;
my %sample;
open (IN,"$cfg2 ") || die "$!\n";
while (<IN>) {
	chomp;
	s/\r$//;s/^\s+//;s/\s+$//;
	next if (/^\#/ || /^$/);
	my @tmp=split /\s+/,$_;
	$para{$tmp[0]}=$tmp[1];
}
close IN;

my @bams =glob("$bam/*/*HISAT_aln.sorted.bam");
my $gtf =(glob("$bam/../Ref_Genome/*.gtf"))[0];
my $idx_prefix    =(glob("$bam/../Ref_Genome/*.fa"))[0];

foreach my $path (sort @bams) {
	my $sam=basename(dirname $path);
	$sample{$sam}=$path;
}
#################################### 
# step 1: Assembly Analysis
########
if ($step==1) {
	open OUT,">$od/work_sh/StringTie.sh" || die;
	&MKDIR ("$od/StringTie");
	foreach my $sam (sort keys %sample) {
		&MKDIR("$od/StringTie/$sam");
		print OUT "$STRINGTIE $sample{$sam} -G $para{Ref_ann} -p 2 -l $sam -o $od/StringTie/$sam/StringTie_asm.gtf \n";
	}
	close OUT;
	$para{Memory}||="15G";
	$para{Queue_type1}||="general.q";
	&qsubOrDie("$od/work_sh/StringTie.sh","$para{Queue_type1}",40,"$para{Memory}");
	&qsubCheck("$od/work_sh/StringTie.sh");
	$step++ unless ($oneStepOnly) ;
	print STDOUT "\n";
}

#################################### 
# step 2: Cuffcompare
########
if ($step==2) {
	open OUT1,">$od/work_sh/Cuffcompare.sh" || die;
	&MKDIR("$od/Cuffcompare");
	foreach my $sam (sort keys %sample) {
		&MKDIR("$od/Cuffcompare/$sam");
		print OUT1 "cd $od/Cuffcompare/$sam && $CUFFCOMPARE -r  $gtf $od/StringTie/$sam/StringTie_asm.gtf \n";
	}
	close OUT1;
	&qsubOrDie("$od/work_sh/Cuffcompare.sh","$para{Queue_type1}",30,"5G");
	&qsubCheck("$od/work_sh/Cuffcompare.sh");
	$step++ unless ($oneStepOnly) ;
	print STDOUT "\n";
}
############cuffmerge
my $index = $para{'Project_key'};
my $gene_expression_dir  = "$od/geneExpression";
my $final_gtf="$od/geneExpression/final_track/final.gtf";
my $final_gff="$od/geneExpression/final_track/final.gff";

if ($step==3) {
	open OUT3,">$od/work_sh/StringTie_Merged.sh" || die;
	&MKDIR("$od/StringTie_Merged");
	print OUT3 "$STRINGTIE  --merge -G $para{Ref_ann} -F 0.1 -T 0.1 -i -o $od/StringTie_Merged/StringTie_merged.gtf ";
	foreach my $sam (sort keys %sample){
		print OUT3 " $od/StringTie/$sam/StringTie_asm.gtf ";
	}
	close OUT3;
        &qsubOrDie("$od/work_sh/StringTie_Merged.sh","$para{Queue_type1}",30,$para{Memory});
	&qsubCheck("$od/work_sh/StringTie_Merged.sh");
	$step++ unless ($oneStepOnly) ;
	print STDOUT "\n";

	`mkdir $od/Compare`	unless(-d "$od/Compare");
	&run_or_die("cd $od/Compare && $GFFCOMPARE -r $para{Ref_ann} $od/StringTie_Merged/StringTie_merged.gtf");

	####get newgene.gff
	mkdir $gene_expression_dir if(!-d $gene_expression_dir);
	my $cmd  ="perl $Bin/basic_analysis/new_gene/newgene_expression.pl -gtf $od/Compare/gffcmp.annotated.gtf -genome $idx_prefix -index $index -od $gene_expression_dir";
	$cmd .= " -singleexon "if(exists $para{SingleExon} && $para{SingleExon} eq "T");
	&run_or_die($cmd);
	&run_or_die("cat $para{'Ref_ann'} $gene_expression_dir/final_track/$index.newGene_final.filtered.gff >$gene_expression_dir/final_track/final.gff");
	&run_or_die("$GFFREAD $gene_expression_dir/final_track/final.gff -T -o $gene_expression_dir/final_track/final.gtf");
}


##################################### 
# step 4: quantify Analysis
#####################################
if ($step==4) {
	open OUT,">$od/work_sh/Ballgown.sh" || die;
	&MKDIR("$od/Ballgown");	
	foreach my $sam (sort keys %sample) {
		&MKDIR("$od/Ballgown/$sam");
		print OUT "cd $od/Ballgown/$sam && $STRINGTIE $sample{$sam} -eB -G $final_gtf -p 6  -o $od/Ballgown/$sam/tringTie_asm.gtf -A $od/Ballgown/$sam/gene_abundence.tab -C $od/Ballgown/$sam/known.cov_refs.gtf \n";
	}
	close OUT;
	&qsubOrDie("$od/work_sh/Ballgown.sh","$para{Queue_type1}",30 ,$para{Memory});
	&qsubCheck("$od/work_sh/Ballgown.sh");
	$step++ unless ($oneStepOnly) ;
	print STDOUT "\n";
}
my $prepDE_dir="$od/geneExpression/prepDE";
if($step ==5){
	&MKDIR($prepDE_dir);
	$para{'Read_length'}=150	if(!exists $para{'Read_length'});
	open OUT,">$od/work_sh/prepDE.sh" || die $!;
	print OUT "/share/nas2/genome/biosoft/Python/2.7.8/bin/python $Bin/stringtie/prepDE.py -i $od/Ballgown -g $prepDE_dir/Total_gene_count.csv -t $prepDE_dir/Total_transcript_count.csv -l $para{'Read_length'} && /usr/bin/dos2unix $prepDE_dir/Total_gene_count.csv $prepDE_dir/Total_transcript_count.csv\n";
	print OUT "/share/nas2/genome/biosoft/Python/2.7.8/bin/python $Bin/stringtie/prepDE_fpkm.py -i $od/Ballgown -g $prepDE_dir/Total_gene_fpkm.csv -t $prepDE_dir/Total_transcript_fpkm.csv -l $para{'Read_length'} && /usr/bin/dos2unix $prepDE_dir/Total_gene_fpkm.csv $prepDE_dir/Total_transcript_fpkm.csv\n";
	close(OUT);
	&qsubOrDie("$od/work_sh/prepDE.sh","$para{Queue_type1}",30 ,$para{Memory});
        &qsubCheck("$od/work_sh/prepDE.sh");
        $step++ unless ($oneStepOnly) ;
        print STDOUT "\n";
}

######################### Subsequence
if ($step==6){
	&run_or_die("$config{Rscript}  $Bin/basic_analysis/fpkm_saturation_Qv1.0.R --count $prepDE_dir/Total_gene_count.csv --fpkm $prepDE_dir/Total_gene_fpkm.csv --o $od/Map_Stat");
	my $new_gene_gff = "$od/geneExpression/final_track/$index.newGene_final.filtered.gff";
	my $new_gene = "$od/geneExpression/final_track/$index.newGene";
	&run_or_die("perl $Bin/basic_analysis/v1.2/bin/get_cds.pl $idx_prefix $new_gene_gff $new_gene");
}

print STDOUT "\nDone. Total elapsed time : ",time()-$BEGIN_TIME,"s\n";

#==================================================================
# subs 
#==================================================================
sub run_or_die{
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
sub show_log{
        my ($txt) = @_ ;
        my $time = time();
        my($sec, $min, $hour, $day, $mon, $year, $wday, $yday, $isdst) = localtime($time);
        $wday = $yday = $isdst = 0;
        my $Time=sprintf("%4d-%02d-%02d %02d:%02d:%02d", $year+1900, $mon+1, $day, $hour, $min, $sec);
        print "$Time:\t$txt\n" ;
}

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
sub step_cmd_process {
    my ($cmd, $sh_name, $sh_dir) = @_;
    my $sh_file = "$sh_dir/$sh_name";
    my $log_file = "$sh_file.log";
    my $flag = 0;
    my $start_time = time();
    &log_current_time("$sh_name start...");
#    &log_current_time("CMD: $cmd");

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

    $flag = system("sh $sh_file > $log_file");
    if ($flag != 0){
        log_current_time("Error: command failed: $cmd");
        exit(1);
    } else {
        my $escaped_time = (time()-$start_time)."s";
        &log_current_time("$sh_name done, escaped time: $escaped_time.");
    }
}
sub MKDIR{ # &MKDIR($out_dir);
	my ($dir)=@_;
	rmdir($dir) if(-d $dir);
	mkdir($dir) if(!-d $dir);
}
sub log_current_time {
     my ($info) = @_;
     my $curr_time = date_time_format(localtime(time()));
     print "[$curr_time] $info\n";
}
sub date_time_format {
    my ($sec, $min, $hour, $day, $mon, $year, $wday, $yday, $isdst)=localtime(time());
    return sprintf("%4d-%02d-%02d %02d:%02d:%02d", $year+1900, $mon+1, $day, $hour, $min, $sec);
}

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
		print "$in\n";
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

sub USAGE {#
	my $usage=<<"USAGE";
Program: Cufflinks &Quantification_Analysis Procedure
Contact: Meng Fei <mengf\@biomarker.com.cn>

Description:
	This program is a Procedure deal with Multiple Samples RNA_Analysis
	Cufflinks &Quantification Combination Designed for RNA Analysis with a Reference Genome

	The program will calculate the Junctions, Transcripts(RABT) Assembly, FPKM of genes & isoforms

Usage:
	-bam		bam dir(tophat or hisat)
	-cfg2		detail config, analysis parameters
	-od		output dir            must be given
	-s		step of the program   option,default 1;
                            1  run Stringtie assembly Transcripts & FPKM
                            2  run Cuffcompare analysis Alt splice
                            3  run Stringtie_merge to merge Transcripts
                            4  run ballgown
                            5  run preyDE
                            6  fpkm_saturation  
	-oneStepOnly      perform one of the pipeline steps above mentioned
USAGE
	print $usage;
	exit;
}
