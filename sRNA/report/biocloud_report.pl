#!/usr/bin/perl -w
use strict;
use warnings;
use Getopt::Long;
use Data::Dumper;
use FindBin qw($Bin $Script);
use File::Basename qw(basename dirname);
my $BEGIN_TIME=time();
my $version="1.0.0";

#######################################################################################
# ------------------------------------------------------------------
# GetOptions
# ------------------------------------------------------------------
my ($indir,$data_cfg,$detail_cfg);

GetOptions(
                "help|?" =>\&USAGE,
                "in:s"=>\$indir,
		"data_cfg:s"=>\$data_cfg,
		"detail_cfg:s"=>\$detail_cfg,
                ) or &USAGE;
&USAGE unless ($indir and $data_cfg and $detail_cfg);

$indir = &ABSOLUTE_DIR("$indir");
$data_cfg = &ABSOLUTE_DIR("$data_cfg");
$detail_cfg = &ABSOLUTE_DIR("$detail_cfg");

my $web_report = "$indir/Web_Report";
my $index_html = "$indir/Web_Report/index.html";
my $need_dir = "$indir/Needed_Data";
my $html_report_zip = "$indir/Needed_Data/biomarker_htmlReport.zip";
my $data_zip = "$indir/Needed_Data/biomarker_Web_Report.zip";
`mkdir $need_dir` unless (-d "$need_dir");
`mkdir $need_dir/Config` unless(-d "$need_dir/Config");
&cmd("cp $data_cfg $need_dir/Config/data.cfg");
&cmd("cp $detail_cfg $need_dir/Config/detail.cfg");

if(-d "$indir/Basic_Analysis/sRNA_Analysis"){
	&cmd("cp $detail_cfg $need_dir/Config/free_pipeline.cfg");
}

if(-f $html_report_zip || -f $data_zip){
	&cmd("rm -r $indir/Needed_Data/*zip");
	print "######If You had have *zip file，the program will delete the Needed Dir because the program thinks you want to generate a new report!!!\n";
}

if(!-d $web_report || !-e $index_html){
	print "######Please Check whether Your Report had generated or not!!!\n";
	die;
}else{
	&cmd("cd $web_report && zip -q -r $html_report_zip HTML/gene HTML/lncRNA HTML/dataassess HTML/circRNA HTML/miRNA HTML/combine src index.html");
	&cmd("cd $web_report && zip -q -r $data_zip BMK_* HTML/gene HTML/lncRNA HTML/dataassess HTML/circRNA HTML/miRNA HTML/combine src index.html Readme.pdf *trans*");
}

if(-d "$indir/Basic_Analysis/Hisat_Stringtie"){
	my $anno_dir = "$indir/Anno_Integrate";
	my $All_longest_transcript = "$indir/Basic_Analysis/Hisat_Stringtie/genePredict/All.longest_transcript.fa";
	my $gff_compare_result = "$indir/Basic_Analysis/Hisat_Stringtie/Compare/gffcmp.annotated.gtf";
	my $all_lncRNA_gtf = "$indir/Basic_Analysis/Hisat_Stringtie/LncPredict/All_LncRNA.gtf";
	my $all_RNA_gtf = "$indir/Basic_Analysis/Hisat_Stringtie/Ballgown/All_RNA.gtf";
	&cmd("cp $indir/PPI/PPI.txt $need_dir");
	&cmd("cp $All_longest_transcript $gff_compare_result $all_lncRNA_gtf $all_RNA_gtf $need_dir");
	&cmd("cp -r $anno_dir $need_dir/");
	&cmd("rm -r $need_dir/Anno_Integrate/*/work_sh");
	&cmd("rm -r $need_dir/Anno_Integrate/New_Anno/*Dir");
	&cmd("rm -r $need_dir/Anno_Integrate/New_Anno/mid");
}


###generate html.cloud for biocloud
&cmd("perl $Bin/html_biocloud.pl -in $indir/Web_Report");

#######################################################################################
print STDOUT "\nDone. Total elapsed time : ",time()-$BEGIN_TIME,"s\n";
#######################################################################################

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
        warn "Warning just for file and dir\n";
        exit;
    }
    chdir $cur_dir;
    return $return;
}

sub cmd {
	my $Cmd = shift;
	system("$Cmd");
	print "$Cmd\n";
}
################################################################################################################

sub GetTime {
    my ($sec, $min, $hour, $day, $mon, $year, $wday, $yday, $isdst)=localtime(time());
    return sprintf("%4d-%02d-%02d %02d:%02d:%02d", $year+1900, $mon+1, $day, $hour, $min, $sec);
}

################################################################################################################
sub USAGE {
    my $usage=<<"USAGE";
ProgramName:
     Version:   $version
     Contact:   niepy <niepy\@biomarker.com.cn> 
Program Date:   2018
     Modify:   
Description:   This program is used to generate Needed_Data dir and zip, generate html.cloud for biocloud
     Usage:
       Options:
	-in <dir>   input directory,xxx format,forced
	-data_cfg	data.cfg,forced
	-detail_cfg	detail.cfg,forced
        -h      help
 Example: perl $0 -in Analysis -data_cfg data.cfg -detail_cfg detail.cfg

USAGE
    print $usage;
    exit;
}
