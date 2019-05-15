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
my ($ref,$circ,$fOut);

GetOptions(
                "help|?" =>\&USAGE,
                "o:s"=>\$fOut,
                "ref:s"=>\$ref,
				"circ:s"=>\$circ,
                ) or &USAGE;
&USAGE unless ($ref and $circ and $fOut);
# ------------------------------------------------------------------
# 
# ------------------------------------------------------------------
my %circfa; 
## 参考文件
open (IN,$ref) or die $!;
open (FILE,$circ) or die $!;
<FILE>;
open (OUT,">$fOut") or die $!;
while(<FILE>)
{
	chomp;
	next if (/\#/);
	my @line=split/\t/,$_;
	my $circRNA_id=$line[0];
	my $chr = (split /:/,$circRNA_id)[0];
 	my $start_end = (split /:/,$circRNA_id)[1];
        my $start=(split /\|/,$start_end)[0];
        my $end=(split /\|/,$start_end)[1];
	$chr =~s/chr//;
	push @ {$circfa{$chr}{$circRNA_id}} , ($start,$end);
	#push @{ $circ{$gene_id}{$chr} }, "$circRNA_id\t$circRNA_start\t$circRNA_end";
}
close FILE;



my %fa;
my $chr;
while (<IN>) {
    chomp;
	if(/>/)
	{
		$_ =~/>(\S+)\s*/;
		$chr = $1;
		next;
	}
	$_=~s/[\r\n]$//;
	$fa{$chr}.=$_;
}
close IN;

foreach my $chr ( keys %fa){
	my $tmp = $chr;
	if($tmp=~/^chr/)
	{

		$tmp=~s/^chr//;
	}
	next unless(exists $circfa{$tmp});
	foreach my $circRNA_id (keys %{$circfa{$tmp}}) {
		my @postion = @{$circfa{$tmp}{$circRNA_id}};
		my $refseq = $fa{$chr};
		print OUT ">$circRNA_id\n";
		my $circseq = substr($refseq,$postion[0]-1,$postion[1]-$postion[0]+1);
		while(length($circseq)>500)
		{
			my $part = substr($circseq,0,500);
			$circseq = substr($circseq,500);
			print OUT "$part\n";
		}
		if(length($circseq)<=500)
		{
			print OUT "$circseq\n";
		}
	}
}
close OUT;
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
     Contact:   Simon Young <yangxh\@biomarker.com.cn> 
Program Date:   2012.07.02
      Modify:   
 Description:   This program is used to ......
       Usage:
        Options:
        -o <file>   circRNA fa file,xxx format,forced
        -ref <file>   reference fa file,forced
        -circ <file>   circRNA result file,forced
        -h      help
USAGE
    print $usage;
    exit;
}
