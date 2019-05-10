#!/usr/bin/perl -w
use strict;
use Cwd qw(abs_path);
use Getopt::Long;
use Data::Dumper;
use FindBin qw($Bin $Script);
use File::Basename qw(basename dirname);
my $BEGIN_TIME=time();
my $version="1.0.0";

# ------------------------------------------------------------------
# GetOptions
# ------------------------------------------------------------------
my ($data_cfg, $detail_cfg, $analysis_dir, $step);

GetOptions(
    "cfg1:s" =>\$data_cfg,
    "cfg2:s" =>\$detail_cfg,
    "dir:s" =>\$analysis_dir,
    "step:i" =>\$step,
    "help|h" =>\&USAGE,
    ) or &USAGE;
&USAGE unless ($detail_cfg and $analysis_dir);

$detail_cfg=abs_path($detail_cfg);
$data_cfg=abs_path($data_cfg);
$analysis_dir=abs_path($analysis_dir);

&log_current_time("$Script start...");
$step ||= 1;

my %detail_cfg;
my %data_cfg;
&detail_cfg_read($data_cfg,\%data_cfg);
&detail_cfg_read($detail_cfg,\%detail_cfg);

# make primary dir
my $sh_dir = "$analysis_dir/work_sh";
my $qc_report_dir = "$analysis_dir/QC_Report";
my $web_dir  = "$analysis_dir/Web_Report";
mkdir $sh_dir unless (-d $sh_dir);
mkdir $qc_report_dir unless (-d $qc_report_dir);
mkdir $web_dir unless (-d $web_dir);

my $project_id = $detail_cfg{Project_id};
my $cmd;
my $qc_report_script = "/share/nas1/niepy/Tools/module_reftrans_qc_report/make_qc_report.pl";
# ------------------------------------------------------------------
# main pipline
# ------------------------------------------------------------------
######################### quality control report
if ($step==1) {
	$cmd = "perl $qc_report_script -cfg1 $data_cfg -cfg2 $detail_cfg -idir $analysis_dir -odir $qc_report_dir ";
	&step_cmd_process($cmd,"8.1.qc_report.sh",$sh_dir);
	$step++;
}

######################### analysis result abstract
if ($step==2) {
	$cmd = "perl $Bin/web_report/Result_extract_local.pl --id $analysis_dir --od $web_dir --species $detail_cfg{Project_key} --genename $detail_cfg{Ref_seq}" ;
	&step_cmd_process($cmd,"8.2.web_report.sh",$sh_dir);
	if(defined $detail_cfg{medical}){
		my @files=();
		push @files,glob("$web_dir/BMK_4_geneExpression/BMK_3_Expression_Statistics/All_gene_*.list");
		push @files,glob("$web_dir/BMK_5_DEG_Analysis/BMK_1_All_DEG/*.xls");
		push @files,glob("$web_dir/BMK_5_DEG_Analysis/BMK_*_vs_*/BMK_1_*/*.DEG_final.xls");
		push @files,glob("$web_dir/BMK_5_DEG_Analysis/BMK_*_vs_*/BMK_1_*/*_vs_*.all.xls");
		push @files,glob("$web_dir/BMK_5_DEG_Analysis/BMK_*_vs_*/BMK_2_*/BMK_1_Anno/*.GO.list.txt");
		push @files,glob("$web_dir/BMK_5_DEG_Analysis/BMK_*_vs_*/BMK_2_*/BMK_1_Anno/*annotation.xls");
		push @files,glob("$web_dir/BMK_5_DEG_Analysis/BMK_*_vs_*/BMK_2_*/BMK_2_Enrichment/*.list");
		push @files,glob("$web_dir/BMK_5_DEG_Analysis/BMK_*_vs_*/BMK_3_GSEA/*GSEA.xls");
		push @files,glob("$web_dir/BMK_5_DEG_Analysis/BMK_*_vs_*/BMK_5_diff_AS_analysis/*.JC.xls");
		push @files,glob("$web_dir/BMK_5_DEG_Analysis/BMK_4_TF_Analysis/BMK_1_TFBS_Analysis/each_DEgeneRes/*txt");
		push @files,glob("$web_dir/BMK_5_DEG_Analysis/BMK_4_TF_Analysis/BMK_1_TFBS_Analysis/*xls");
		push @files,glob("$web_dir/BMK_5_DEG_Analysis/BMK_4_TF_Analysis/BMK_2_TF_activity/TFs_activity_grn.xls");
		push @files,glob("$web_dir/BMK_6_Alt_splice/*.AS.list.xls");
		push @files,glob("$web_dir/BMK_6_Alt_splice/*.fpkm");
		push @files,glob("$web_dir/BMK_8_Gene_Structure_Optimize/*.geneStructure.optimize.xls");
		push @files,glob("$web_dir/BMK_7_SNP_Analysis/final.*");
		open(LIST,">$analysis_dir/gene_symbol.list")||die $!;
		foreach my $f(@files){
			print LIST "$f\n";
		}
		close(LIST);
		my $cmd="perl $Bin/convert_symbol.pl -i $analysis_dir/gene_symbol.list ";
		my $path=dirname $detail_cfg{Ref_seq};		
		if(-e "$path/id_name.list"){$cmd .=" -id $path/id_name.list ";}else{ $cmd .=" -db $detail_cfg{medical} ";}
		print "$cmd\n";	`$cmd`;
		
		@files=glob("$web_dir/BMK_5_DEG_Analysis/BMK_2_DEG_PPI/*.sif");		
		push @files,glob("$web_dir/BMK_5_DEG_Analysis/BMK_2_DEG_PPI/*.txt");
		if(@files>0){
			my $cmd="perl $Bin/convert_symbol_ppi.pl -i ".join(",",@files);
			if(-e "$path/id_name.list"){$cmd .=" -id $path/id_name.list ";}else{ $cmd .=" -db $detail_cfg{medical} ";}
			print "$cmd\n";	`$cmd`;
		}
	}
	$step++;
}


######################### xml report
if ($step==3) {
=head
	$cmd ="perl $Bin/xml_report/Ref_xml_report_v1.6.2.pl --id $web_dir  --cfg $detail_cfg  --pp ";
	$cmd.=" --only_analysis "if(!exists $data_cfg{Basecall_stat});
=cut
	print "------local xml Start-----\n";
	$cmd ="perl $Bin/xml_report/build_xml.pl --id $web_dir --cfg1 $data_cfg  --cfg2 $detail_cfg ";
	$cmd.="&& perl $Bin/package.pl -in $analysis_dir";
	&step_cmd_process($cmd,"8.3.local_xml_report.sh",$sh_dir);
	print "------local xml End-----\n\n";
	$step++;
}

if ($step==4) {
=head
        $cmd ="perl $Bin/xml_report/Ref_xml_report_v1.6.2.pl --id $web_dir  --cfg $detail_cfg  --pp ";
        $cmd.=" --only_analysis "if(!exists $data_cfg{Basecall_stat});
=cut
        print "------biocloud xml Start-----\n";
        $cmd ="perl $Bin/xml_report/build_xml.pl --id $web_dir --cfg1 $data_cfg  --cfg2 $detail_cfg --cloud ";
        $cmd.="&& perl $Bin/package.pl -in $analysis_dir";
        &step_cmd_process($cmd,"8.4.biocloud_xml_report.sh",$sh_dir);
        print "------biocloud xml End-----\n\n";
        $step++;
}
#else{
	#print "Step 4 is used to generate xml for cloud,if you want to do this ,please set you options if '-step 4 -cloud'\n";
        #die;
#}
#######################################################################################
my $elapse_time = (time()-$BEGIN_TIME)."s";
&log_current_time("$Script done. Total elapsed time: $elapse_time.");
#######################################################################################
# ------------------------------------------------------------------
# sub function
# ------------------------------------------------------------------
#############################################################################################################
sub data_cfg_read {
    my ($cfg_file, $data_cfg) = @_;
    my $sample_id;

    open (CFG, $cfg_file) or die "$!: $cfg_file\n";
    while (<CFG>) {
        chomp;
        next if (/^\s*$/ or /^#/);
        $_=~s/^\s+//g;

        if (/^Sample/) {
            $sample_id=(split/\s+/,$_)[1];
        }
        if ($_=~/^fq1/ || $_=~/^fq2/) {
            my $file=(split/\s+/,$_)[1];
            die "$file is not exist!\n" unless (-e $file);

            $data_cfg->{rawdata}{$sample_id}{fq1}=$file if $_=~/^fq1/;
            $data_cfg->{rawdata}{$sample_id}{fq2}=$file if $_=~/^fq2/;
        }
    }
    close CFG;
    &log_current_time("data config done.");
}

#############################################################################################################
sub detail_cfg_read {
    my ($cfg_file, $detail_cfg) = @_;

    open (CFG,$cfg_file ) or die "$!: $cfg_file\n";
    while (<CFG>) {
        chomp;
        next if (/^\s*$/ or /^#/);
        $_=~s/^\s+//g;
        my ($key, $value) = (split /\s+/)[0,1];
	next unless(defined $key);
	$detail_cfg->{$key} = $value;

    }
    close CFG;
    &log_current_time("detail config done.");
}

#############################################################################################################
sub step_cmd_process {
    my ($cmd, $sh_name, $sh_dir) = @_;
    my $sh_file = "$sh_dir/$sh_name";
    my $log_file = "$sh_file.log";
    my $flag = 0;
    my $start_time = time();
    &log_current_time("$sh_name start...");
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

    $flag = system("sh $sh_file > $log_file");
    if ($flag != 0){
        log_current_time("Error: command failed: $cmd");
        exit(1);
    } else {
        my $escaped_time = (time()-$start_time)."s";
        &log_current_time("$sh_name done, escaped time: $escaped_time.");
    }
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
----------------------------------------------------------------------------------------------------
   Program: $Script
   Version: $version
   Contact: Simon Young <yangxh\@biomarker.com.cn> 
      Date: 2014-11-13

     Usage:
            --cfg1       <FILE>  data config, basecall stat
            --cfg2       <FILE>  detail config, analysis parameters
            --dir       <DIR>   analysis output directory

            --step      <INT>   step to start from                  [1]
                          1     QC report
                          2     web report
                          3     local xml report
                          4     xml report for biocloud, step4 is used to generate xml for biocloud

	    --PL        <STR>   abbr. of Project Leader\'s name
            --CSE       <STR>   abbr. of Customer Service Executive\'s name
            --h                 help documents

   Example:
            perl $Script --cfg detail.cfg --dir Analysis/

----------------------------------------------------------------------------------------------------
USAGE
	print $usage;
	exit;
}
