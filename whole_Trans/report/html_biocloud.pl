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
my ($indir);

GetOptions(
                "help|?" =>\&USAGE,
                "in:s"=>\$indir,
                ) or &USAGE;
&USAGE unless ($indir);

$indir = &ABSOLUTE_DIR("$indir");

my @html = glob "$indir/HTML/*/*html";
foreach my $html (@html){
	my $cloud_html = "$html.cloud";
	open(IN,"$html") or die $!;
	open(OUT,">$cloud_html") or die $!;
	while(<IN>){
		chomp;
		s/src\/images\/logo\.jpg/logo_image_path/ if (/src\/images\/logo\.jpg/);
		if($_=~/src\/images/){
			$_=~s/src\/images/project_Template_Path\/src\/images/g;
		}
		if($_=~/DEG_Analysis/){
			$_=~s/\.\.\/\.\.\/BMK_/project_Template_Path\/\.\.\/BMK_/;
		}
		print OUT "$_\n";
	}
	close OUT;
	close IN;
}




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
Description:	generate html.cloud for biocloud
     Usage:
       Options:
	-in <dir>   input directory,xxx format,forced
        -h      help
 Example: perl $0 -in Web_Report

USAGE
    print $usage;
    exit;
}
