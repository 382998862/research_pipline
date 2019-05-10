#!/usr/bin/perl -w
use strict;
use Cwd qw(abs_path);
use Getopt::Long;
use Data::Dumper;
use FindBin qw($Bin $Script);
use File::Basename qw(basename dirname);
use threads;
use newPerlBase;
my %config=%{readconf("$Bin/../config/db_file.cfg")}; 
my $Title="Ref_Trans";  
my $version="v2.9"; 
my $BEGIN_TIME=time();
# ------------------------------------------------------------------
# GetOptions
# ------------------------------------------------------------------
my ($data_cfg, $detail_cfg, $bamdir, $odir, $step, $fusion, $test);

GetOptions(
    "cfg1:s" =>\$data_cfg,
    "cfg2:s" =>\$detail_cfg,
    "bamdir:s" =>\$bamdir,
    "od:s"   =>\$odir,
    "step:s" =>\$step,
    "help|h" =>\&USAGE,
    ) or &USAGE;
&USAGE unless ($detail_cfg and $odir and $bamdir);
    

######################################################
mkdirOrDie("$odir") unless (-d $odir);	

$odir		=abs_path($odir);
$data_cfg	=abs_path($data_cfg)	if (defined $data_cfg);
$detail_cfg	=abs_path($detail_cfg);
$bamdir		=abs_path($bamdir);

# ------------------------------------------------------------------
# read configs and steps
# ------------------------------------------------------------------
my (%data_cfg, %detail_cfg);

# read data config
&data_cfg_read($data_cfg,\%data_cfg) if (defined $data_cfg );

# read detail config
&detail_cfg_read($detail_cfg,\%detail_cfg);

# read steps
my %step;
$step ||=  join ',',(1..3);
&steps_process($step,\%step);
die "Error :You need provide  data_cfg for Gene Fusion" if (defined $step{4}  &&  !defined $data_cfg );
# ------------------------------------------------------------------
# main pipline
# ------------------------------------------------------------------
# work shell backup dir
my $sh_dir = "$odir/work_sh";

mkdirOrDie("$sh_dir") unless (-d $sh_dir);
my $cmd;

my $index = $detail_cfg{Project_key};
my $project_id = $detail_cfg{Project_id};

$detail_cfg{Call_SNP}||="GATK";
open (SH,">$sh_dir/structure_and_expression.sh") || die $! ;

######################### SNP Analysis
if ($step{1}) {
        if ($detail_cfg{Call_SNP} eq "GATK") {
                $cmd = "perl $Bin/snp_analysis/v2.0/SNP_Trans_main_Ref.pl -cfg $detail_cfg  -tophat $bamdir  -gff $odir/geneExpression/final_track/$index.newGene_final.filtered.gff -od $odir/SNP_Analysis   -qphred $data_cfg{Qphred} ";
                print SH $cmd,"\n";
        }elsif($detail_cfg{Call_SNP} eq "Samtools"){
                $cmd = "perl $Bin/snp_analysis/samtools/SNP_Trans_main_Ref_Samtools.pl -cfg $detail_cfg  -tophat $bamdir  -gff $odir/geneExpression/final_track/$index.newGene_final.filtered.gff -od $odir/SNP_Analysis ";
                print SH $cmd,"\n";       
        }else{
                print "Wrong Parameter Call_SNP: $detail_cfg{Call_SNP}\n";
                die;
        }
        
}

######################## Assembly_Quantification_Annoation_Difference 
if ($step{2}) {
     $cmd = "perl $Bin/Step3.2_basic_analysis_main.pl  --bamdir $bamdir  --cfg2 $detail_cfg --od $odir";
     print "$cmd\n";
     print SH $cmd,"\n";
}

#################################### Map stat 
if ($step{3}) {
    $cmd = "perl $Bin/Step3.3_map_stat.pl -bamdir $bamdir -queue  $detail_cfg{Queue_type2}   -od  $odir/Map_Stat  " ;
    $cmd.="-medical $detail_cfg{medical} "	if(defined $detail_cfg{medical});
     print SH $cmd,"\n";
     
}

#################################### Gene Fusion
if (exists $detail_cfg{medical}) {  ###融合基因分析 fusionmap
	my %genotype=(
		"GRCh37"        =>"Human.hg19",
                "GRCh38"        =>"Human.B38",
                "GRCm38"        =>"Mouse.B38",
                "Rnor_6.0"      =>"Rat.B6.0",
                "Rnor_5.0"      =>"Rat.B5.0",
                );
	$cmd="perl $Bin/gene_fusion/fusionmap3.pl --indir $bamdir --type $genotype{$detail_cfg{medical}}  --od  $odir/Gene_Fusion --script_cfg $Bin/gene_fusion/script.cfg   --cfg2 $detail_cfg ";
	print SH $cmd,"\n";
}


close(SH);
&qsubOrDie("$odir/work_sh/structure_and_expression.sh",  $detail_cfg{Queue_type2} ,  5, "20G");

#######################################################################################
print STDOUT "\nDone. Total elapsed time : ",time()-$BEGIN_TIME,"s\n";
#######################################################################################
# ------------------------------------------------------------------
# sub function
# ------------------------------------------------------------------
#############################################################################################################
sub data_cfg_read {
    &log_current_time("data config check:");
    my ($cfg_file, $data_cfg) = @_;
    my $sample_id;

    open (CFG, $cfg_file) or die "$!: $cfg_file\n";
    while (<CFG>) {
        chomp;
        next if (/^\s+/ or /^#/);
        if (/^Qphred/) {
            $data_cfg->{Qphred} = (split /\s+/,$_)[1];
        }
        if (/^Sample/) {
            $sample_id=(split/\s+/,$_)[1];
        }
        if ($_=~/^fq1/ || $_=~/^fq2/) {
            my $file=(split/\s+/,$_)[1];
            die "$file is not exist!\n" unless (-e $file);

            $data_cfg->{rawdata}{$sample_id}{fq1}=$file if $_=~/^fq1/;
            $data_cfg->{rawdata}{$sample_id}{fq2}=$file if $_=~/^fq2/;
        }
		if ($_=~/^Basecall_stat/) {
            my $file=(split/\s+/,$_)[1];
            die "$file is not exist!\n" unless (-e $file);
            $data_cfg->{Basecall_stat}= $file;
        }
    }
    close CFG;

    if (defined $data_cfg->{Qphred}) {
        print "Qphred: $data_cfg->{Qphred}\n";
    } else {
        $data_cfg->{Qphred} = 33;
        print "Qphred: $data_cfg->{Qphred} [ASCII encoding type of quality score of rawdata is unknown, and default is 33.]\n";
    }

    $data_cfg->{sample_num} = scalar keys %{$data_cfg->{rawdata}};
    print "sample_number: $data_cfg->{sample_num}\n";

    for my $s (sort keys %{$data_cfg->{rawdata}}) {
        print "${s}_fq1: $data_cfg->{rawdata}{$s}{fq1}\n${s}_fq2: $data_cfg->{rawdata}{$s}{fq2}\n";
    }
    &log_current_time("data config check done.\n");
}

#############################################################################################################
sub detail_cfg_read {
    &log_current_time("detail config check:");
    my ($cfg_file, $detail_cfg) = @_;

    open (CFG,$cfg_file ) or die "$!: $cfg_file\n";
    while (<CFG>) {
        chomp;
        next if (/^\s+/ or /^#/);
        my ($key, $value) = (split /\s+/)[0,1];
	next unless(defined $key and $value);
        $detail_cfg->{$key} = $value;
        if ($key eq 'Project_name' or $key eq 'Customer_info' or $key eq 'Project_id' or $key eq 'Project_key') {
            $detail_cfg->{$key} = $value;
        }
        if ($key eq 'Known_unigene' or $key eq 'Known_pep') {
            die "$key: $value is not exist!\n" unless (-e $value);
            $detail_cfg->{$key} = $value;
        }

        if ($key eq 'Known_anno') {
            #die "$key: $value is not illegal!\n" unless (-e "$value/02.gene-annotation" and -e "$value/Result");
            $detail_cfg->{$key} = $value;
        }
        if ($key eq 'Ref_seq' or $key eq 'Ref_ann') {
            die "$key: $value is not exist!\n" unless (-e $value);
            $detail_cfg->{$key} = $value;
        }
        if ($key =~/^SNP_/) {
            $detail_cfg->{$key} = $value;
        }
        if ($key eq 'nr' or $key eq 'Swissprot' or $key eq 'Kegg' or $key eq 'Pfam' or $key eq 'Cog' or $key eq 'Kog') {
            die "$key: $value is not exist!\n" unless (-e $value);
            $detail_cfg->{anno}{$key} = $value;
        }
        if ($key eq 'Queue_type1') {
            $detail_cfg->{$key} = $value;
        }
        if ($key eq 'Queue_type2') {
            $detail_cfg->{$key} = $value;
        }
        print "$key: $value\n" if (exists $detail_cfg->{$key});
    }
    close CFG;
    die "Must choose Queue_type1 and Queue_type2 !\n" unless (exists $detail_cfg->{Queue_type1} or exists $detail_cfg->{Queue_type2});
    &log_current_time("detail config check done.");
}

#############################################################################################################
sub steps_process {
    print "\n";
    &log_current_time("step check:");
    my ($step_str,  $step_hash) = @_;
    my %_step_ = (
        '1' => 'SNP_Analysis',
        '2' => 'Assembly_Quantification_Annoation_Difference ',
        '3' => 'Map_Stat',
        '4' => 'Gene_Fusion',
    );
        for my $s (split /,/,$step_str) {
            if ($s =~/^[1-4]$/) {
                $step_hash->{$s} = $_step_{$s};
            } else {
                print "ERROR: illegal steps specified by --step.\n";
                die;
            }
    }

    print "steps_to_run: ".(join ", ",(map {sprintf "$_.$step_hash->{$_}"} sort keys %{$step_hash})).".\n";
    &log_current_time("step check done.\n");
}

#############################################################################################################
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
    qsubOrDie("$sh_file",$queue,$cpu,$vf);	
    qsubCheck("$sh_file"); 
}

#############################################################################################################
#############################################################################################################
sub log_current_time {
     # get parameter
     my ($info) = @_;

     # get current time with string
     my $curr_time = date_time_format(localtime(time()));

     # print info with time
     print "[$curr_time] $info\n";
}

#############################################################################################################
sub date_time_format {
    my ($sec, $min, $hour, $day, $mon, $year, $wday, $yday, $isdst)=localtime(time());
    return sprintf("%4d-%02d-%02d %02d:%02d:%02d", $year+1900, $mon+1, $day, $hour, $min, $sec);
}


#############################################################################################################
sub USAGE {
	my $usage=<<"USAGE";
----------------------------------------------------------------------------------------------------------------------------
   Program: $Script
   Version: $version
   Contact: Simon Young <yangxh\@biomarker.com.cn> 
      Date: 2014-11-13
   Modifier: Wang Yajing <wangyj\@biomarker.com.cn> 
      Date: 2015-07-16
   Modifier: Wang Yajing <wangyj\@biomarker.com.cn> 
      Date: 2016-06-20
     Usage:
            --cfg1      <FILE>  data config, rawdata path for gene fusion analysis 
            --cfg2      <FILE>  detail config, analysis parameters & refseq path
            --od        <DIR>   analysis output directory
            --bamdir	<DIR>	bam dir
            --step      <INT>   steps to run (split by comma) [1] default : 1,2,3
                          1     SNP_Analysis
                          2     Assembly_Quantification_Annoation_Difference
                          3     Map_Stat
                          4     Gene_Fusion
            --cloud             analysis at biocloud
            --PL        <STR>   abbr. of Project Leader\'s name
            --CSE     <STR>   abbr. of Customer Service Executive\'s name
            --h                 help documents

   Example:
            perl $Script --cfg1 ref_trans.data.cfg --cfg2 ref_trans.detail.cfg  --bamdir Tophat --od Structure_and_Expression  -PL wangyj 

----------------------------------------------------------------------------------------------------------------------------
USAGE
	print $usage;
	exit;
}
