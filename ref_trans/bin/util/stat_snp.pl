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
my ($id);
GetOptions(
				"help|?" =>\&USAGE,
				"id:s"=>\$id,
				) or &USAGE;
&USAGE unless ($id);

my @snp_files=sort glob "$id/*.snp.stat.xls";

if (@snp_files) {
    print "\n",join("\t","#Sample","Total","Genic","Intergenic","Transition","Transversion","Heterozygosity"),"\n";

    foreach my $snp_file (@snp_files) {
        my $file=basename($snp_file);
        my ($sample)=(split /\./,$file)[0];
        my $total_line=`tail -1 $snp_file`;
        chomp $total_line;
        my @data=split /\t+/,$total_line;
        shift @data;
        unshift @data,$sample;

        $data[$_]= sprintf("%.2f",$data[$_])."%" foreach (4..6);

        print join("\t",@data),"\n";
    }

    print "\n";
} else {
    die "ERROR: illigal input directory, without files *.snp.stat.xls. \n";
}




#######################################################################################
#print STDOUT "\nDone. Total elapsed time : ",time()-$BEGIN_TIME,"s\n";
#######################################################################################

# ------------------------------------------------------------------
# sub function
# ------------------------------------------------------------------
################################################################################################################
sub USAGE {#
    my $usage=<<"_USAGE_";
ProgramName: $Script
    Version: $version
    Contact: Simon Young <yangxh\@biomarker.com.cn> 
      Usage: 
             $Script -id your_analysis_od/SNP_Analysis

_USAGE_
    print $usage;
    exit;
}
