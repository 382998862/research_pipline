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

my @ftable=sort glob "$id/geneExpression/T*.mappedStat.xls";

if (@ftable) {
    print "\nmapping results statistic:\n";
    print join("\t","#Sample_ID","Total_Reads","Mapped_Reads","Mapped_Ratio","Uniq_Mapped_Reads","Uniq_Mapped_Ratio"),"\n";

    for my $ftable (@ftable) {
        my ($sample_id) = $ftable =~/\/([^\/]+).mappedStat.xls$/;
        my ($clean_num, $mapped_num, $mapped_ratio, $uniq_mapped_num, $uniq_mapped_ratio) = (0,0,'',0,0);

        open (TAB,$ftable) or die $!;
        while (<TAB>) {
            chomp;
            my @col = split /\t/;
            next unless (@col);
            if ($col[0] =~/Total Reads/) {
                $clean_num = format_figure($col[1]);
            }
            if ($col[0] =~/mapped Reads/) {
                $mapped_num = format_figure($col[1]);
                ($mapped_ratio) = $col[2] =~/([0-9\.]+)%$/;
                $mapped_ratio = format_figure($mapped_ratio)."%";
            }
            if ($col[0] =~/Uniq Map/) {
                my $clean_num_tmp = $clean_num;
                $clean_num_tmp =~s/,//g;
                $uniq_mapped_num = format_figure($col[1]);
                $uniq_mapped_ratio = format_figure($col[1]*100/$clean_num_tmp)."%";
            }
        }
        close TAB;

        print join("\t", $sample_id, $clean_num, $mapped_num, $mapped_ratio, $uniq_mapped_num, $uniq_mapped_ratio),"\n";
    }

    print "\n";
} else {
    die "ERROR: illigal input directory, without files *.snp.stat.xls. \n";
}

if (-e (glob "$id/geneExpression/final_track/*.newGene.longest_transcript.fa")[0] ) {
    my $new_gene_num = `grep ">" $id/geneExpression/final_track/*.newGene.longest_transcript.fa |wc -l `;
    chomp $new_gene_num;
    print "number of novel genes assembled by Cufflinks: ".format_figure($new_gene_num).".\n\n";
}

#######################################################################################
#print STDOUT "\nDone. Total elapsed time : ",time()-$BEGIN_TIME,"s\n";
#######################################################################################

# ------------------------------------------------------------------
# sub function
# ------------------------------------------------------------------
################################################################################################################
#整数格式：3位一个逗号
sub Integer_Three_Digit{#
	my $interger = shift;
	$interger=~s/(?<=\d)(?=(\d\d\d)+$)/,/g;
	return $interger;
}

################################################################################################################
#整数格式：3位一个逗号
#小数格式：小数点后两位
sub format_figure{#
	my $figure = shift;
	if (!defined $figure) {
		die;
	}
	if ($figure=~/\./) {
        if ($figure == 100) {
            $figure = 100;
        } else {
            $figure = sprintf("%.2f",$figure);
        }
	}else{
		$figure = Integer_Three_Digit($figure);
	}
	return $figure;
}

################################################################################################################
sub USAGE {#
    my $usage=<<"_USAGE_";
ProgramName: $Script
    Version: $version
    Contact: Simon Young <yangxh\@biomarker.com.cn> 
      Usage: 
             $Script -id your_analysis_od/Basic_Analysis

_USAGE_
    print $usage;
    exit;
}
