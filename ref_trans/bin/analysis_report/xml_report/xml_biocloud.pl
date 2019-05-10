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
#
#######################################################################################
# ------------------------------------------------------------------
# GetOptions
# ------------------------------------------------------------------
my ($fIn,$fOut);

GetOptions(
                "help|?" =>\&USAGE,
                "o:s"=>\$fOut,
                "i:s"=>\$fIn,
                ) or &USAGE;
&USAGE unless ($fIn and $fOut);

$fIn=&ABSOLUTE_DIR($fIn);
system "mkdir -p $fOut" unless (-d $fOut) ;
$fOut=&ABSOLUTE_DIR($fOut);

#修改*.KEGG.list.html.cloud的logo和二次链接路径
my @kegg_list=glob("$fIn/BMK_9_html/*/*.html");
foreach my $list_file(@kegg_list) {
       	my $list_file_cloud="$list_file.cloud";
	my $name=basename (dirname $list_file);
        open IN,"<",$list_file;
        open OUT,">",$list_file_cloud;
        while(<IN>){
        		s/src\/images\/logo\.jpg/logo_image_path/ if (/src\/images\/logo\.jpg/);
			if($_=~/src\/images/){
				$_=~s/src\/images/project_Template_Path\/src\/images/g;
			}
			if($_=~/BMK_5_DEG_Analysis/){
				$_=~s/\.\.\/\.\.\/BMK_5_DEG_Analysis/project_Template_Path\/\.\.\/BMK_5_DEG_Analysis/;
			}
        		print OUT $_;
        	}
        	close OUT;
        	close IN;
	}
my @deus=glob("$fIn/BMK_5_DEG_Analysis/BMK_*_vs_*/BMK_4_DEXSeqReport/*testForDEU.html");
foreach my $deu(@deus){
	`cp $deu $deu.cloud`;
	my $cmd="sed -i 's/files/project_DEXSeqReport_Path\\/files/' $deu.cloud";
	`$cmd`;
}
`cp $fIn/configtest_raw.xml $fOut/configtest.xml`;

#######################################################################################
print STDOUT "\nDone. Total elapsed time : ",time()-$BEGIN_TIME,"s\n";
#######################################################################################

# ------------------------------------------------------------------
# sub function
# ------------------------------------------------------------------
################################################################################################################

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
     Contact:   Liu Tao <liut\@biomarker.com.cn> 
Program Date:   2016.03.11
      Modify:   
 Description:   This program is used to convert Web report from bio to that of cloud
       Usage:
        Options:
        -i <DIR>    input directory of Web report directory of bio, forced

        -o <file>   output directory, forced

        -h      help
   Example:
            perl $Script -i ./Web_Report -o ./Web_Report

----------------------------------------------------------------------------------------------------------------------------

USAGE
    print $usage;
    exit;
}

