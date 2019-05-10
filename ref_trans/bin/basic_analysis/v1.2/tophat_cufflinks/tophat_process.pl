#!/usr/bin/perl -w
use strict;
use warnings;
use Getopt::Long;
use Data::Dumper;
use FindBin qw($Bin $Script);
use File::Basename qw(basename dirname);
use newPerlBase;
my $BEGIN_TIME=time();
my $version="1.0.0";
my %config=%{readconf("$Bin/../../../../config/db_file.cfg")}; 

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

#
# log file 
#
my $startTime = GetTime();
my $user = `whoami`;  chomp $user;
my $workDir = `pwd`;  chomp $workDir;
my $task = "perl $Bin/$Script ".join("\t", @Original_ARGV);

open ($log,">", "$od/Tophat_cufflinks.".time().".log") or die $!;

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

#my $GFFREAD_BIN    = "/share/nas2/genome/biosoft/cufflinks/2.1.1/gffread";     # ~ 2014-12-17
my $GFFREAD_BIN     = "$config{cufflink}/gffread";     # 2014-12-17 ~ 
my $TOPHAT_BIN      = "$Bin/tophat/v2.0.13/tophat2";     # 2015-01-28 ~ 
## cufflinks suite
my $CUFFLINKS_BIN   = "$config{cufflink}/cufflinks";   # 2014-12-17 ~ 
my $CUFFCOMPARE_BIN = "$config{cufflink}/cuffcompare"; # 2014-12-17 ~ 
my $CUFFMERGE_BIN   = "$config{cufflink}/cuffmerge";   # 2014-12-17 ~ 
my $CUFFDIFF_BIN    = "$config{cufflink}/cuffdiff";    # 2014-12-17 ~ 
my $CUFFQUANT_BIN   = "$config{cufflink}/cuffquant";   # 2015-08-04 
my $CUFFNORM_BIN    = "$config{cufflink}/cuffnorm";    # 2015-08-04 
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
	if ($tmp[0]=~m/Sample/) {
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
}
else {
	open IN,"$od/totalRead.stat.xls" || die;
	while (<IN>) {
		chomp;
		next if (/^$/);
		my @tmp=split/\s+/,$_;
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
my $index = $para{'Project_key'};

if ($step!=1) {
    $genome = basename($para{Ref_seq});
    $idx_prefix = basename($para{Ref_seq});
    $idx_prefix =~s/.fa$//;
    $gff = basename($para{Ref_ann});
    $gtf = basename($para{Ref_ann});
    $gtf =~s/\.gff3?$//i;
    $gtf.= ".gtf";
}
if ($step==1) {
	print STDOUT "=== check and build bowtie index ===\n";
	print $log  "=== check and build bowtie index ===\n";

	&MKDIR("$od/Ref_Genome");
	chdir "$od/Ref_Genome";
	$genome = basename($para{Ref_seq});
    $idx_prefix = basename($para{Ref_seq});
    $idx_prefix =~s/.fa$//;

    if (!-f "$idx_prefix.1.bt2" and !-f "$idx_prefix.1.bt2l" ) {
        my $genome_path = dirname($para{Ref_seq});
    
        if (-e "$genome_path/$idx_prefix.1.bt2") {
            system "ln -s $genome_path/$genome $genome_path/$idx_prefix.*.bt2 ./";
        } elsif (-e "$genome_path/$idx_prefix.1.bt2l") {   # for large genomes more than about 4 billion nucleotides in length
            system "ln -s $genome_path/$genome $genome_path/$idx_prefix.*.bt2l ./";
        } else {
            system "ln -s $genome_path/$genome ./";
            system "$config{bowtie2_build} $genome $idx_prefix";
        }
	}

	################################## gff2gtf
    $gff = basename($para{Ref_ann});
    $gtf = basename($para{Ref_ann});
    $gtf =~s/\.gff3?$//i;
    $gtf.= ".gtf";
	system "ln -s $para{Ref_ann} ./" unless (-f $gff);
	system "$GFFREAD_BIN $gff -T -o $gtf";
	system "ln -s $para{Ref_seq} ./" unless (-f "$idx_prefix.fa");
	system "$config{samtools} faidx $idx_prefix.fa ";
	system "grep '>' $idx_prefix.fa >$index.hdrs";

	chdir "../";


	$step++ unless ($oneStepOnly) ;
	print STDOUT "\n";
	print $log "\n";

}

#################################### 
#
# step 2: Align the RNA-seq read to genome using Tophat2 + bowtie2
#
#########

if ($step==2) {
	print STDOUT "=== Align the RNA-seq read to genome using Tophat2 + bowtie2  ===\n";
	print STDOUT "Tophat mapping shell file: $od/work_sh/Tophat.sh\n";
	print $log "=== Align the RNA-seq read to genome using Tophat2 + bowtie2  ===\n";
	print $log "Tophat mapping shell file: $od/work_sh/Tophat.sh\n";
	
	#
	# write shell
	#
	open OUT,">$od/work_sh/Tophat.sh" || die;
	&MKDIR("$od/Tophat");
	foreach my $sam (sort keys %sample) {
		&MKDIR("$od/Tophat/$sam");
		print OUT "cd $od/Tophat/$sam && ";
#		print OUT "$TOPHAT_BIN -p 6 -G $od/Ref_Genome/$gtf --read-edit-dist $para{Mismatch} --bowtie1 -N $para{Mismatch} -r $para{Insert_size} --library-type $para{Lib_type} --fusion-search -o ./ $od/Ref_Genome/$genome $sample{$sam}{FQ1} $sample{$sam}{FQ2} && ";
        print OUT "$TOPHAT_BIN --output-dir ./ --read-mismatches $para{Mismatch} --read-edit-dist $para{Mismatch} --max-intron-length 5000000 --library-type $para{Lib_type} --num-threads 8 --GTF $od/Ref_Genome/$gtf --mate-inner-dist $para{Insert_size} $od/Ref_Genome/$idx_prefix $sample{$sam}{FQ1} $sample{$sam}{FQ2} && ";
		print OUT "$config{samtools} index accepted_hits.bam accepted_hits.bam.bai &&\n";
	}
	close OUT;
	#
	# qsub
	#
	$para{Memory}||="15G";
	$para{Queue_type1}||="general.q";
	&qsubOrDie("$od/work_sh/Tophat.sh","$para{Queue_type1}",18,"$para{Memory}");
	&qsubCheck("$od/work_sh/Tophat.sh");
	

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
	print STDOUT "shell file: $od/work_sh/Tophat_bam_stat.sh\n";
	print STDOUT "shell file: $od/work_sh/genome_bam2depth.sh\n";
	print STDOUT "shell file: $od/work_sh/genome_Checkgraph.sh\n";
	print $log "=== Statistic bam files  ===\n";
	print $log "shell file: $od/work_sh/Tophat_bam_stat.sh\n";
	print $log "shell file: $od/work_sh/genome_bam2depth.sh\n";
	print $log "shell file: $od/work_sh/genome_Checkgraph.sh\n";
	
	#
	# write shell 
	#
	open OUT1, ">$od/work_sh/Tophat_bam_stat.sh"   || die;
	open OUT2, ">$od/work_sh/genome_bam2depth.sh"  || die;
	open OUT3, ">$od/work_sh/genome_Checkgraph.sh" || die;
	open OUT4, ">$od/work_sh/plot_ReadDensity.sh"  || die;
	&MKDIR("$od/Map_Stat");
	my @str_type_stat;
	foreach my $sam (sort keys %sample) {
#		&MKDIR("$od/Map_Stat/$sam");
		push @str_type_stat, "$od/Map_Stat/$sam.type.stat";
		print OUT1 "perl $Bin/bin/bam2map_stat.pl -i $sam -bam $od/Tophat/$sam/accepted_hits.bam -totalRead $total_read{$sam} -od $od/Map_Stat\n";
		print OUT2 "$config{bedtools} genomecov -ibam $od/Tophat/$sam/accepted_hits.bam -g $idx_prefix.fa.fai -d -strand +  |awk '\$3!=0' >$od/Tophat/$sam/$sam.sort.bam.plus.depth && \n bedtools genomecov -ibam $od/Tophat/$sam/accepted_hits.bam -g $idx_prefix.fa.fai -d -strand -  |awk '\$3!=0' >$od/Tophat/$sam/$sam.sort.bam.minus.depth &&\n samtools depth $od/Tophat/$sam/accepted_hits.bam >$od/Tophat/$sam/$sam.sort.bam.depth &&\n";
		print OUT3 "perl $Bin/bin/get_percent_of_exon_intro_inter_by_gff_v1.4.pl -gff $od/Ref_Genome/$gff -i $od/Tophat/$sam/$sam.sort.bam.depth -od $od/Map_Stat -index $sam &&\n";
#		print OUT4 "perl $Bin/bin/plotReadDensity.pl -i $od/Tophat/$sam/$sam.sort.bam.depth -o $od/Map_Stat -k $sam -png -n 10 &&\n";
        print OUT4 "perl $Bin/bin/plotReadDensity2.pl -i $od/Tophat/$sam/$sam.sort.bam.plus.depth:$od/Tophat/$sam/$sam.sort.bam.minus.depth -o $od/Map_Stat/ -k $sam  &&\n";
	}
	close OUT1;
	close OUT2;
	close OUT3;
	close OUT4;


	# qsub
	#
	&qsubOrDie("$od/work_sh/Tophat_bam_stat.sh", "$para{Queue_type2}",   18, "20G");
	#`sh $od/work_sh/Tophat_bam_stat.sh`;
	&qsubOrDie("$od/work_sh/genome_bam2depth.sh", "$para{Queue_type2}",  18, "15G");
	&qsubOrDie("$od/work_sh/genome_Checkgraph.sh",  "$para{Queue_type2}", 18, "5G");
	&qsubOrDie("$od/work_sh/plot_ReadDensity.sh",  "$para{Queue_type2}",  18, "5G");
	#&qsubCheck("$od/work_sh/Tophat_bam_stat.sh");
	
	#my $str_type_stat = join " ", @str_type_stat;
	#`perl $Bin/bin/map_total_type.pl -i $str_type_stat -o $od/Map_Stat/Total_Type.png`;

	$step++ unless ($oneStepOnly) ;
	print STDOUT "\n";
	print $log "\n";
}

#################################### 
#
# step 4: Cufflinks Assembly Analysis
#
########

if ($step==4) {
	print STDOUT "=== Cufflinks Assembly Analysis ===\n";
	print STDOUT "cufflinks assembly shell file: $od/work_sh/Cufflinks.sh\n";
	print $log "=== Cufflinks Assembly Analysis ===\n";
	print $log "cufflinks assembly shell file: $od/work_sh/Cufflinks.sh\n"; 
	
	#
	# write shell 
	#
	open OUT,">$od/work_sh/Cufflinks.sh" || die;
	&MKDIR ("$od/Cufflinks");
	foreach my $sam (sort keys %sample) {
		&MKDIR("$od/Cufflinks/$sam");
		print OUT "cd $od/Cufflinks/$sam && ";
		print OUT "$CUFFLINKS_BIN -o ./ -p 6 -g $od/Ref_Genome/$gtf --library-type $para{Lib_type} -u -L $sam $od/Tophat/$sam/accepted_hits.bam && \n";
	}
	close OUT;
	#
	# qsub
	#
	$para{Memory}||="15G";
	$para{Queue_type1}||="general.q";
	&qsubOrDie("$od/work_sh/Cufflinks.sh","$para{Queue_type1}",18,"$para{Memory}");
	&qsubCheck("$od/work_sh/Cufflinks.sh");
	
	$step++ unless ($oneStepOnly) ;
	print STDOUT "\n";
	print $log "\n";
}

#################################### 
#
# step 5 & 6: Cuffcompare & Cuffmerge Analysis
#
########

if ($step==5) {
	print STDOUT "=== Cuffcompare: compare gtf of each sample with known gene model ===\n";
	print STDOUT "cuffcompare shell file: $od/work_sh/Cuffcompare.sh\n";
	print $log "=== Cuffcompare: compare gtf of each sample with known gene model ===\n";	
	print $log "cuffcompare shell file: $od/work_sh/Cuffcompare.sh\n";
	
	#
	# write shell 
	#
	open OUT1,">$od/work_sh/Cuffcompare.sh" || die;
	&MKDIR("$od/Cuffcompare");
	open OUT2,">$od/assembly_GTF_list.txt" || die;
	foreach my $sam (sort keys %sample) {
		&MKDIR("$od/Cuffcompare/$sam");
		print OUT1 "cd $od/Cuffcompare/$sam && ";
		print OUT1 "$CUFFCOMPARE_BIN -r $od/Ref_Genome/$gtf $od/Cufflinks/$sam/transcripts.gtf \n";
		print OUT2 "$od/Cufflinks/$sam/transcripts.gtf\n";
	}
	close OUT1;
	close OUT2;
	
	#
	# qsub
	#
	&qsubOrDie("$od/work_sh/Cuffcompare.sh","$para{Queue_type1}",18,"5G");
	&qsubCheck("$od/work_sh/Cuffcompare.sh");
	
	$step++ unless ($oneStepOnly) ;
	print STDOUT "\n";
	print $log "\n";
}

if ($step==6) {
	print STDOUT "=== Cuffmerge: merge gtf files of all samples ===\n";
	print STDOUT "cuffmerge shell file: $od/work_sh/Cuffmerge.sh\n";
	print $log "=== Cuffmerge: merge gtf files of all samples ===\n";
	print $log "cuffmerge shell file: $od/work_sh/Cuffmerge.sh\n";

	open OUT3,">$od/work_sh/Cuffmerge.sh" || die;
	&MKDIR("$od/Cuffmerge");
	#print OUT3 "$CUFFMERGE_BIN -o $od/Cuffmerge -g $od/Ref_Genome/$gtf -s $od/Ref_Genome/$genome -p 6 $od/assembly_GTF_list.txt\n";
	print OUT3 "$CUFFMERGE_BIN -o $od/Cuffmerge -g $od/Ref_Genome/$gtf  -p 6 $od/assembly_GTF_list.txt\n";
	close OUT3;

	#system "$CUFFMERGE_BIN -o $od/Cuffmerge -g $od/Ref_Genome/$gtf -s $od/Ref_Genome/$genome -p 6 $od/assembly_GTF_list.txt\n";
	system "$CUFFMERGE_BIN -o $od/Cuffmerge -g $od/Ref_Genome/$gtf  -p 6 $od/assembly_GTF_list.txt\n";

	$step++ unless ($oneStepOnly) ;
	print STDOUT "\n";
	print $log "\n";
}

##################################### 
#
# step 7: Cuffquant Analysis
#
#########

if ($step==7) {
	print STDOUT "=== Cuffquant: estimate isoform expression ===\n";
	print STDOUT "Cuffquant shell file: $od/work_sh/Cuffquant.sh\n";
	print $log "=== Cuffquant: estimate isoform expression  ===\n";
	print $log "Cuffquant shell file: $od/work_sh/Cuffquant.sh\n";
	
	#
	# write shell
	#
	open OUT,">$od/work_sh/Cuffquant.sh" || die;
	&MKDIR("$od/Cuffquant");
	foreach my $sam (sort keys %sample) {
		&MKDIR("$od/Cuffquant/$sam");
		print OUT "cd $od/Cuffquant/$sam && ";
		print OUT "$CUFFQUANT_BIN -o ./ -p 6 -b $od/Ref_Genome/$idx_prefix.fa -u  --max-bundle-frags  50000000  --library-type $para{Lib_type} $od/Cuffmerge/merged.gtf  $od/Tophat/$sam/accepted_hits.bam && \n";
	}
	close OUT;
	

	&qsubOrDie("$od/work_sh/Cuffquant.sh","$para{Queue_type1}",18,$para{Memory});
	&qsubCheck("$od/work_sh/Cuffquant.sh");
	
	$step++ unless ($oneStepOnly) ;
	print STDOUT "\n";
	print $log "\n";
}

if ($step==8) {
	print STDOUT "=== Cuffnorm:  generate tables of expression values  ===\n";
	print STDOUT "Cuffnorm shell file: $od/work_sh/Cuffnorm.sh\n";
	print $log "=== Cuffnorm:  generate tables of expression values  ===\n";
	print $log "Cuffnorm shell file: $od/work_sh/Cuffnorm.sh\n";
	

	open OUT1,">$od/work_sh/Cuffnorm.sh" || die;
	&MKDIR("$od/Cuffnorm");
	if ((keys %sample) != 1) {
		print OUT1 "$CUFFNORM_BIN --library-type $para{Lib_type} ";
		print OUT1 " --output-format cuffdiff " ;
		print OUT1 " -o $od/Cuffnorm -q -p 6   ";
		print OUT1 "-L ";
		my $lable;
		my $sam_str;
		foreach my $sam (sort keys %sample) {
			$lable.="$sam".",";
			$sam_str.="$od/Cuffquant/$sam/abundances.cxb ";
		}
		$lable=~s/,$//;
		$sam_str=~s/\s+$//;
		print OUT1 "$lable $od/Cuffmerge/merged.gtf  $sam_str\n";
		close OUT1;
	}else {
		print OUT1 "$CUFFNORM_BIN --library-type $para{Lib_type} ";
		print OUT1 " --output-format cuffdiff " ;
		print OUT1 " -o $od/Cuffnorm -q -p 6   ";
		print OUT1 "-L ";
		my $lable;
		my $sam_str;
		foreach my $sam (sort keys %sample) {
			$lable="$sam".","."NULL";
			$sam_str="$od/Cuffquant/$sam/abundances.cxb "."$od/Cuffquant/$sam/abundances.cxb";
		}
		$lable=~s/,$//;
		$sam_str=~s/\s+$//;
		print OUT1 "$lable $od/Cuffmerge/merged.gtf  $sam_str\n";
		close OUT1;		
	
	
	}

	#
	# qsub 
	#
	&qsubOrDie("$od/work_sh/Cuffnorm.sh","$para{Queue_type1}",18,$para{Memory});
	&qsubCheck("$od/work_sh/Cuffnorm.sh");
	
	$step++ unless ($oneStepOnly) ;
	print STDOUT "\n";
	print $log "\n";
}

###################################### Extract Gene Expression For DEG Analysis
#if ($step==8) {
#	open OUT,">$od/work_sh/GeneExpression.sh" || die;
#	&MKDIR("$od/geneExpression");
#	print OUT "perl $Bin/bin/FPKM_extract.pl -list $od/NewGene/$para{SPE}.Cuffmerge.ID.list -id $od/Cuffdiff -od $od/geneExpression -i $para{SPE}\n";
#	close OUT;
#
#	system "perl $Bin/bin/FPKM_extract.pl -list $od/NewGene/$para{SPE}.Cuffmerge.ID.list -id $od/Cuffdiff -od $od/geneExpression -i $para{SPE}";
#	$step=9;
#}

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

sub USAGE {#
	my $usage=<<"USAGE";
Program: Tophat&Cufflinks_Analysis Procedure
Version: $version
Contact: Meng Fei <mengf\@biomarker.com.cn>

Description:
	This program is a Procedure deal with Multiple Samples RNA_Analysis
	Tophat+Cufflinks Combination£ºDesigned for RNA Analysis with a Reference Genome

	The program will calculate the Junctions, Transcripts(RABT) Assembly, FPKM of genes & isoforms

Usage:
    -cfg1             data config, rawdata & refseq path
    -cfg2             detail config, analysis parameters
	-od               output dir            must be given
	-s                step of the program   option,default 1;
                            1  Check the index of Genome & build index for alignment
                            2  run Tophat analysis
                            3  stat Tophat bam out with Genome
                            4  run Cufflinks assembly Transcripts & FPKM
                            5  run Cuffcompare analysis Alt splice
                            6  run Cuffmerge to merge Transcripts
                            7  run Cuffquant
                            8  run Cuffnorm
	-oneStepOnly      perform one of the pipeline steps above mentioned
USAGE
	print $usage;
	exit;
}
