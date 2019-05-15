#!/usr/bin/perl -w
use strict;
no strict 'refs';
use warnings;
use Getopt::Long;
use Data::Dumper;
use autodie;
use FindBin qw($Bin $Script);
use File::Basename qw(basename dirname);
my $BEGIN_TIME=time();
my $version="1.0.0";

#######################################################################################

# ------------------------------------------------------------------
# GetOptions
# ------------------------------------------------------------------
my ($idir,$o);

GetOptions(
				"help|?" =>\&USAGE,
				"id:s"=>\$idir,
				"o:s"=>\$o,
				) or &USAGE;
&USAGE unless ($idir and $o);
`rm $idir/All.mappedStat.xls` if(-e "$idir/All.mappedStat.xls");
my %table_info;
push @{$table_info{"table"}}, [("BMK-ID","Total Reads","Mapped Reads","Uniq Map Reads","Multiple Map Reads","Reads Map to '+'","Reads Map to '-'")];
my @ftable = glob("$idir/*.mappedStat.xls");
    for my $ftable (@ftable) {
        my ($sample_id) = $ftable =~/\/([^\/]+).mappedStat.xls$/;
        my ($clean_num, $mapped_num, $mapped_ratio, $uniq_mapped_num, $uniq_mapped_ratio,$mui_mapped_num, $mui_mapped_ratio,$plus,$plus_ratio,$minus,$minus_ratio) = (0,0,'',0,'',0,'',0,'',0,'');
        open (TAB,$ftable) or die $!;
        $/="\n";
        while (<TAB>) {
	chomp;
	my @col = split /\t/;
	next unless (@col);
	if ($col[0] =~/Total Reads/) {
		$clean_num = format_figure($col[1]);
		next;
	}
	if ($col[0] =~/mapped Reads/) {
                $mapped_num = format_figure($col[1]);
                ($mapped_ratio) = $col[2] =~/([0-9\.]+)%$/;
                $mapped_ratio = $mapped_num.'('.format_figure($mapped_ratio).'%)';
		next;
	}
	if ($col[0] =~/Uniq Map/) {
                $uniq_mapped_num = format_figure($col[1]);
		($uniq_mapped_ratio) = $col[2] =~/([0-9\.]+)%$/;
                $uniq_mapped_ratio =$uniq_mapped_num.'('.format_figure($uniq_mapped_ratio).'%)';
		next;
	}
	if ($col[0] =~/Multiple Map/) {
                $mui_mapped_num = format_figure($col[1]);
		($mui_mapped_ratio) = $col[2] =~/([0-9\.]+)%$/;
                $mui_mapped_ratio = $mui_mapped_num.'('.format_figure($mui_mapped_ratio).'%)';
		next;
	}
	if ($col[0] =~/Map Plus Strand/) {
                $plus = format_figure($col[1]);
		($plus_ratio) = $col[2] =~/([0-9\.]+)%$/;
                $plus_ratio = $plus.'('.format_figure($plus_ratio).'%)';
		next;
	}
	if ($col[0] =~/Map Minus Strand/) {
                $minus = format_figure($col[1]);
		($minus_ratio) = $col[2] =~/([0-9\.]+)%$/;
                $minus_ratio = $minus.'('.format_figure($minus_ratio).'%)';
		next;
	}
	}
        push @{$table_info{"table"}}, [($sample_id, $clean_num, $mapped_ratio, $uniq_mapped_ratio, $mui_mapped_ratio, $plus_ratio,$minus_ratio)];
        close (TAB) ;
    }
open O,">$o";
for(my $i=0;$i<@{$table_info{"table"}};$i++){
	my $line="";
	for (my $j=0;$j<@{$table_info{"table"}[$i]};$j++){
		$line.=$table_info{"table"}[$i][$j]."\t";
	}
	$line=~s/\t$//;
	print O  "$line\n";
}
close O;



######################
#整数格式：3位一个逗号
sub Integer_Three_Digit{#
	my $interger = shift;
	$interger=~s/(?<=\d)(?=(\d\d\d)+$)/,/g;
	return $interger;
}

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


sub USAGE {
        my $usage=<<"USAGE";
#-----------------------------------------------------------------------------------------
  Program: $Script
  Version: $version
  Contact: songmm<songmm\@biomarker.com.cn>
     Date: 
 Modifier: songmm <songmm\@biomarker.com.cn>
     Date: 2016-07-07
    Usage:
      Options:
      --id  <dir>   input dir
      --o <dir>   output file

#-----------------------------------------------------------------------------------------
USAGE
        print $usage;
        exit;
}

