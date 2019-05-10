#!/usr/bin/perl -w

use strict;
use Cwd qw(abs_path);
use Getopt::Long;
use Data::Dumper;
use FindBin qw($Bin $Script);
use File::Basename qw(basename dirname);
use threads;
use newPerlBase;
my %config=%{readconf("$Bin/config/db_file.cfg")}; 
my $Title="Ref_Trans";  
my $version="v2.9.7"; 
# ------------------------------------------------------------------
# GetOptions
# ------------------------------------------------------------------
my ($data_cfg,$detail_cfg, $odir, $step, $step_by_step, $fusion);

GetOptions(
    "cfg1:s" =>\$data_cfg,
    "cfg2:s" =>\$detail_cfg,
    "od:s"   =>\$odir,
    "step:s" =>\$step,
    "sbs"    =>\$step_by_step,
    "help|h" =>\&USAGE,
    ) or &USAGE;
&USAGE unless ($data_cfg and $detail_cfg and $odir);
mkdirOrDie("$odir") unless (-d $odir);
# read detail config
my (%data_cfg, %detail_cfg);
&detail_cfg_read($detail_cfg,\%detail_cfg);

######################################################
$odir=abs_path($odir);
$data_cfg=abs_path($data_cfg);
$detail_cfg=abs_path($detail_cfg);
my $cmd;

# read data config
=head
print "------Check Phred value-----\n";
mkdirOrDie("$odir/Data_Assess") unless (-d "$odir/Data_Assess");
$cmd="perl $Bin/bin/data_assess/Phred_Change/Phred_Change.pl -cfg $data_cfg -od $odir/Data_Assess";
print "$cmd\n";
system "$cmd";
=cut

&data_cfg_read($data_cfg,\%data_cfg);

# read steps
my %step;
$step ||= ($step_by_step) ? '1' : join ',',(1..4);
&steps_process($step,$step_by_step,\%step);


# ------------------------------------------------------------------
# main pipline
# ------------------------------------------------------------------
my $cfg_dir = "$odir/Config";
mkdirOrDie("$cfg_dir") unless (-d $cfg_dir);
system "cp $data_cfg  $cfg_dir";
system "cp $detail_cfg  $cfg_dir";

# work shell backup dir
my $sh_dir = "$odir/work_sh";
mkdirOrDie("$sh_dir") unless (-d $sh_dir);

my $index = $detail_cfg{Project_key};
my $project_id = $detail_cfg{Project_id};

open LOG, ">$odir/bioinfor_pipeline.log";
&show_log("This project will be analysed in 4 steps.");
######################### Data Assesss
if ($step{1}) {
	$cmd = "perl $Bin/bin/data_assess/rna_seq_data_assess.pl -config $data_cfg -outdir $odir/Data_Assess  -queue $detail_cfg{Queue_type2} ";
	stepStart(1," Data assess");
	&show_log("step_1: Data assess start.");
	&step_cmd_process($cmd,\%step,1,$sh_dir,$detail_cfg{Queue_type2},3,"2G");
	stepTime(1);	
	&show_log("step_1: Data assess finished.");
}

#########################Tophat
if ($step{2}) {
	$cmd = "perl $Bin/bin/Step2_map.pl   -cfg1  $data_cfg -cfg2  $detail_cfg  -od  $odir/Mapping ";
	stepStart(2," Mapping ");
	&show_log("step_2: mapping start.");
	&step_cmd_process($cmd,\%step,2,$sh_dir,$detail_cfg{Queue_type2},1,"2G");
	stepTime(2);	
	&show_log("step_2: mapping finished.");
}

#########################  Structure_and_Expression_Analysis
if ($step{3}) {            
	$cmd = "perl $Bin/bin/Step3_structure_and_expression.pl  --cfg1  $data_cfg --cfg2  $detail_cfg --od $odir/Structure_and_Expression  --bamdir  $odir/Mapping/Hisat ";                         
	stepStart(3,"Gene Structure & Expression Analysis ");
	&show_log("step_3: Gene Structure & Expression Analysis   start.");
	&step_cmd_process($cmd,\%step,3,$sh_dir,$detail_cfg{Queue_type2},5,"20G");
	&show_log("step_3: Gene Structure & Expression Analysis   finished.");
	stepTime(3);
}


######################### Analysis Reports and Keep Genome information
if ($step{4}) {
        $cmd = "perl $Bin/bin/analysis_report/analysis_report.pl --cfg1 $data_cfg --cfg2 $detail_cfg --dir $odir ";
	stepStart(4," Analysis Reports and Keep Genome information"); 
	&show_log("step_4: Analysis Reports and Keep Genome information start.");
	&step_cmd_process($cmd,\%step,4,$sh_dir,$detail_cfg{Queue_type2},3,"2G");
	
	###get allGene protein seq and position by xugl 2016-10-13
	my $new_cds_predict=(glob("$odir/Structure_and_Expression/Anno_Integrate/New_Anno/Unigene_CDS_Predict/*.newGene.longest_transcript.pep.fa"))[0];
	my $known_gene_prt="$detail_cfg{Known_anno}/Known.longest_transcript.pep.fa";
	my $newGene_fa=(glob("$odir/Web_Report/BMK_2_NewGene/*longest_transcript.fa"))[0];

	$detail_cfg{Pfam}="/share/nas2/database/pfam/201703/Pfam-A.hmm"	if(exists $detail_cfg{Pfam});
	$cmd ="perl $Bin/bin/ORF_predict_integrate.pl -fa $newGene_fa -od $odir/Needed_Data/Allgene_annoNseq -known_pep $known_gene_prt -pfam $detail_cfg{Pfam} -prefix $detail_cfg{Project_key} -new_cds_predict $new_cds_predict";
	$cmd.=" -ss "if($detail_cfg{Lib_type} ne "fr-unstranded");
	&runOrDie($cmd);
	stepTime(4);
	&show_log("step_4: Analysis Reports and Keep Genome information finished.");
}

&show_log("Analysis of this project is completed.");
close LOG;
totalTime();	
########################################################################
# sub function
########################################################################
sub data_cfg_read {
    &log_current_time("data config check:");
    my ($cfg_file, $data_cfg) = @_;
    my $sample_id;

    open (CFG, $cfg_file) or die "$!: $cfg_file\n";
    open(O_CFG,">$cfg_file.temp")||die $!;
    while (<CFG>) {
        chomp;
        next if (/^\s+$/ or /^#/);        
        if (/^Qphred/) {
            $data_cfg->{Qphred} = (split /\s+/,$_)[1];
        }
        if (/^Sample/) {
            $sample_id=(split/\s+/,$_)[1];
        }
        if ($_=~/^fq1/ || $_=~/^fq2/) {
            my $file=(split/\s+/,$_)[1];
            $data_cfg->{rawdata}{$sample_id}{fq1}=$file if $_=~/^fq1/;
            $data_cfg->{rawdata}{$sample_id}{fq2}=$file if $_=~/^fq2/;
        }
	if ($_=~/^Basecall_stat/) {
            my $file=(split/\s+/,$_)[1];
            $data_cfg->{Basecall_stat}= $file;
        }
        print O_CFG "$_\n";
    }
    close CFG;
    system "mv $cfg_file.temp $cfg_file";
    if (defined $data_cfg->{Qphred}) {
        print "Qphred: $data_cfg->{Qphred}\n";
    } else {
        $data_cfg->{Qphred} = 33;
        print "Qphred: $data_cfg->{Qphred} [ASCII encoding type of quality score of rawdata is unknown, and default is 33.]\n";
    }

    $data_cfg->{sample_num} = scalar keys %{$data_cfg->{rawdata}};

    &log_current_time("data config check done.\n");
}

#############################################################################################################
sub detail_cfg_read {
    &log_current_time("detail config check:");
    my ($cfg_file, $detail_cfg) = @_;

    open (CFG,$cfg_file ) or die "$!: $cfg_file\n";
    while (<CFG>) {
        chomp;
        next if (/^\s+$/ or /^#/);
        my ($key, $value) = (split /\s+/)[0,1];
        next unless(defined $key);
        $detail_cfg->{$key} = $value;
        if ($key eq 'Project_name' or $key eq 'Customer_info' or $key eq 'Project_id' or $key eq 'Project_key') {
            $detail_cfg->{$key} = $value;
        }
        if ($key eq 'Known_unigene' or $key eq 'Known_pep') {
            die "$key: $value is not exist!\n" unless (-e $value);
            $detail_cfg->{$key} = $value;
        }

        if ($key eq 'Known_anno') {
            $detail_cfg->{$key} = $value;
        }
        if ($key eq 'Ref_seq' ) {
            die "$key: $value is not exist!\n" unless (-e $value);
            $detail_cfg->{$key} = $value;
        }
	 if ($key eq 'Ref_ann') {
            die "$key: $value is not exist!\n" unless (-e $value);
	    &Check_Ref_ann($value);
            $detail_cfg->{$key} = $value;
        }

        if ($key =~/^SNP_/) {
            $detail_cfg->{$key} = $value;
        }
        if ($key eq 'nr' or $key eq 'Swissprot' or $key eq 'Kegg' or $key eq 'Pfam' or $key eq 'Cog' or $key eq 'Kog'  or $key eq 'eggNOG'  ) {
            die "$key: $value is not exist!\n" unless (-e $value);
            $detail_cfg->{anno}{$key} = $value;
        }
        if ($key eq 'Queue_type1') {
            $detail_cfg->{$key} = $value;
        }
        if ($key eq 'Queue_type2') {
            $detail_cfg->{$key} = $value;
        }
    }
    close CFG;
    die "Must choose Queue_type1 and Queue_type2 !\n" unless (exists $detail_cfg->{Queue_type1} or exists $detail_cfg->{Queue_type2});
    &log_current_time("detail config check done.");
}

#############################################################################################################
sub Check_Ref_ann {
	my ($gff)=@_;
	print "$gff\n";
	my $gtf_name=basename($gff);
	$gtf_name =~ s/\.gff$/\.gtf/ig if ($gff =~/\.gff$/);
	$gtf_name =~ s/\.gff3$/\.gtf/ig if ($gff =~/\.gff3$/);
	$gtf_name=basename($gff) if ($gff =~/\.gtf$/);
	mkdir "$odir/log" unless -d "$odir/log";
	system "$config{cufflink}/gffread $gff -T -o $odir/log/$gtf_name " unless ($gff =~/\.gtf$/);
	system "cp $gff $odir/log/ " if ($gff =~/\.gtf$/);
	system "perl $Bin/bin/util/validate_gtf.pl $odir/log/$gtf_name > $odir/log/gtf.log";
	open IN,"cat $odir/log/gtf.log |" or die $!;
	my $n=0;
	while (<IN>){
		if (/<start> field is greater than <stop> field/ || /Missing 'gene_id' attribute in <attributes> field/ || /Inconsistent <strand> value across gene_id/){$n++;warn "the $gff have wrong format !!!!\n $_ \n";}
	}
	close IN;
	
	system "rm $odir/log/$gtf_name" if ($n==0);
	die "please check the $gff format !!!\n get information form the $odir/log/gtf.log !!" unless ($n==0);	
	
}
sub steps_process {
    print "\n";
    &log_current_time("step check:");
    my ($step_str, $step_by_step, $step_hash) = @_;
    my %_step_ = (
        '1' => 'Data_Assess',
        '2' => 'Map',
        '3' => 'Structure_and_Expression_Analysis',
        '4' => 'Analysis_Report',
        '5' => 'Backup',
    );

    if ($step_by_step) {
        print "step-by-step: ON\n";
        if ($step_str =~/^[1-5]$/) {
            for my $s ($step_str..5) {
                $step_hash->{$s} = $_step_{$s};
            }
        } else {
            print "ERROR: illegal steps specified by --step, or options specified by --step and --sbs clashed. \n";
            die;
        }
    } else {
        if ($step eq join ',',(1..5)) {
            print "step-by-step: ON\n";
        } else {
            print "step-by-step: OFF\n";
        }

        for my $s (split /,/,$step_str) {
            if ($s =~/^[1-5]$/) {
                $step_hash->{$s} = $_step_{$s};
            } else {
                print "ERROR: illegal steps specified by --step.\n";
                die;
            }
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
    $flag=system("sh $sh_file >$log_file 2>&1");
    if ($flag != 0){
	print "Error: command fail: $sh_file\n";
	exit(1);
    }
}

#############################################################################################################
# &show_log("cmd")
sub show_log()
{
    my ($txt) = @_ ;
    my $time = time();
    my $Time = date_time_format(localtime(time()));
    print LOG "$Time:\t$txt\n" ;
    return ($time) ;
}

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
   Modifier: Liu xiaoshuang <Liuxs\@biomarker.com.cn>
      Date: 2016-08-10
     Usage:
            --cfg1      <FILE>  data config, rawdata path
            --cfg2      <FILE>  detail config, analysis parameters & refseq path
            --od        <DIR>   analysis output directory

            --step      <INT>   step to start from or steps to run (split by comma) [1]
                          1     Data_Assess
                          2     Mapping and Assembly
                          3     Structure_and_Expression_Analysis
                          4     Analysis_Report
            --sbs               step-by-step analysis, start from the step specified by -step
                                or not, only run the steps specified by -step
            --h                 help documents

   Example:
            perl $Script --cfg1 ref_trans.data.cfg --cfg2 ref_trans.detail.cfg --od Analysis/

----------------------------------------------------------------------------------------------------------------------------
USAGE
	print $usage;
	exit;
}
