#!/usr/bin/perl -w
use strict;
use warnings;
use Getopt::Long;
use Data::Dumper;
use FindBin qw($Bin $Script);
use File::Basename qw(basename dirname);
use newPerlBase;
my %config=%{readconf("$Bin/../../config/db_file.cfg")};
my $BEGIN_TIME=time();
my $version="1.0.0";
my $indir;

# ------------------------------------------------------------------
# GetOptions
# ------------------------------------------------------------------
my ($fIn,$fOut);
GetOptions(
                "help|?" =>\&USAGE,
                "in:s"=>\$indir,
                ) or &USAGE;
&USAGE unless ($indir);
$indir = ABSOLUTE_DIR("$indir");

my $data_zip = "$indir/Needed_Data/biomarker_Web_Report.zip";
my $html_zip = "$indir/Needed_Data/biomarker_htmlReport.zip";

if(!-e $data_zip){
	&run_or_die("cd $indir/Web_Report && zip -q -r $indir/Needed_Data/biomarker_Web_Report.zip BMK* index.html Readme.pdf ref_trans_full_table.xls");}
if(!-e $html_zip){
	&run_or_die("cd $indir/Web_Report/ && zip -q -r $indir/Needed_Data/biomarker_htmlReport.zip BMK_*_html index.html");
}
system "$config{'python1'} $config{'generate_small_png'} -web_report $indir/Web_Report/";


my @deu=glob("$indir/Web_Report/BMK_5_DEG_Analysis/*vs*/BMK_*_DEXSeqReport/files");
if ($#deu !=0){
	foreach my $dir (@deu){
		system "rm -rf $dir " if (-d $dir);
	}
}
&run_or_die("perl $Bin/xml_report/xml_biocloud.pl -i $indir/Web_Report -o $indir/Web_Report");


#######################################################################################
print STDOUT "\nDone. Total elapsed time : ",time()-$BEGIN_TIME,"s\n";
#######################################################################################
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
     Contact:   Yang nan <yangn\@biomarker.com.cn> 
Program Date:   2016.03.08
      Modify:   
 Description:   This program is used to package the zip file to Need_Data......
       Usage:
        Options:
        -in <dir>   input directory,xxx format,forced
        -h      help

USAGE
    print $usage;
    exit;
}
