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
my (@fIn,$fOut);

GetOptions(
                "help|?" =>\&USAGE,
                "o:s"=>\$fOut,
                "i:s{1,}"=>\@fIn,
                ) or &USAGE;
&USAGE unless (@fIn and $fOut);
# ------------------------------------------------------------------
# 
# ------------------------------------------------------------------
my $cfg="ref_trans.detail.cfg";
map {&USAGE unless (-e "$_")} @fIn;
my $out="gene_pos.list";
#my $out=basename($gff);
#$out=~s/\.gff3?$/_pos\.list/g;
open (OUT,">$fOut/$out") or die $!;
print OUT "#Gene\tChr\tStart\tEnd\tStrand\n";
for my $gff (@fIn){
open (INN,$gff) or die $!;
while (<INN>) {
    chomp;
    next if (/^$/||/\#/);
    my ($chr,$type,$stat,$end,$strand,$info)=(split /\t/,$_)[0,2,3,4,6,-1];
    if ($type=~/gene/i) {
        my $id=(split /\;/,$info)[0];
        $id=~s/ID=//ig;
        print OUT "$id\t$chr\t$stat\t$end\t$strand\n";
    }
}
close INN;
}
close OUT;
if (-e "$fOut/$out") {
    system "touch $fOut/gff_format.Check";
}
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

sub max{#&max(lists or arry);
    #求列表中的最大值
    my $max=shift;
    my $temp;
    while (@_) {
        $temp=shift;
        $max=$max>$temp?$max:$temp;
    }
    return $max;
}

################################################################################################################

sub min{#&min(lists or arry);
    #求列表中的最小值
    my $min=shift;
    my $temp;
    while (@_) {
        $temp=shift;
        $min=$min<$temp?$min:$temp;
    }
    return $min;
}

################################################################################################################

sub revcom(){#&revcom($ref_seq);
    #获取字符串序列的反向互补序列，以字符串形式返回。ATTCCC->GGGAAT
    my $seq=shift;
    $seq=~tr/ATCGatcg/TAGCtagc/;
    $seq=reverse $seq;
    return uc $seq;           
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
     Contact:   ziy <ziy\@biomarker.com.cn> 
Program Date:   2016.01.20
      Modify:   
 Description:   This program is used to ......
       Usage:   perl $Script
        Options:
        -i <dir>   input gff file,multi gff sep by blank,forced

        -o <dir>   output dir,forced

        -h      help

USAGE
    print $usage;
    exit;
}
