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
my (@fIn,$fOut,$know_miRNA_pos,$expression);

GetOptions(
                "help|?" =>\&USAGE,
                "o:s"=>\$fOut,
                "k:s"=>\$know_miRNA_pos,
                "e:s"=>\$expression,
                "i:s{1,}"=>\@fIn,
                ) or &USAGE;
&USAGE unless (@fIn and $fOut and $know_miRNA_pos and $expression);
# ------------------------------------------------------------------
# 
# ------------------------------------------------------------------

### creat know miRNA position hash
my %track;
open (IN,$know_miRNA_pos) or die $!;
while (<IN>) {
    chomp;
    next if (/^$/||/\#/);
    my ($miRNA_id,$chr,$start,$end,$strand,$precursor_ID,$precursor_Start,$precursor_End,$precursor_Strand)=split /\t/;
	#my $miRNA_id = $col[0];
	my $message="$chr\t$start\t$end\t$strand\t$precursor_ID\t$precursor_Start\t$precursor_End\t$precursor_Strand";
	$track{$miRNA_id} = $message;
}
close IN;

my %track_exp;
open (IN,$expression) or die $!;
while (<IN>) {
    chomp;
    next if (/^$/||/\#/||/^miRNA\t/);#unconservative_X_63824
	my $miRNA_all_id = (split /\t/)[0];
	if ($miRNA_all_id =~ /conservative_/) {
		#my ($miRNA_id,$chr,$start,$end,$strand,$precursor_ID,$precursor_Start,$precursor_End,$precursor_Strand)=split /\t/;
		#my $miRNA_id = $col[0];
		my @all = split /_/,$miRNA_all_id;
		my $head = $all[0];
		my $tail = join "_",@all[1..@all-1];
		#my $message="$chr\t$start\t$end\t$strand\t$precursor_ID\t$precursor_Start\t$precursor_End\t$precursor_Strand";
		$track_exp{$tail} = $miRNA_all_id;
	}else{
		$track_exp{$miRNA_all_id} = $miRNA_all_id;
	}
}
close IN;

map {&USAGE unless (-e "$_")} @fIn;
my $miRNA_pos="miRNA_pos.list";
my $miRNA_mature="miRNA_mature.fa";
my $miRNA_pre="miRNA_pre_seq.fa";
#my $out=basename($gff);
#$out=~s/\.gff3?$/_pos\.list/g;
open (OUT,">$fOut/$miRNA_pos") or die $!;
open (OUT1,">$fOut/$miRNA_mature") or die $!;
open (OUT2,">$fOut/$miRNA_pre") or die $!;
#print OUT "#miRNA\tChr\tStart\tEnd\tStrand\n";
print OUT "#miRNA_ID\tChr\tStart\tEnd\tStrand\tprecursor_ID\tprecursor_Start\tprecursor_End\tprecursor_Strand\n";


for my $miRNAsum (@fIn){
    open (INN,$miRNAsum) or die $!;
    while (<INN>) {
        chomp;
        next if (/^$/||/^\#/||/^miRNA\t/);
        my ($miRNA_id,$chr,$start,$end,$strand)=(split /\t/,$_)[0,4,6,7,5];
        if($track{$miRNA_id}){
			if ($track_exp{$miRNA_id}) {
				print OUT "$track_exp{$miRNA_id}\t$track{$miRNA_id}\n";
			}else{
				print OUT "$miRNA_id\t$track{$miRNA_id}\n";
			}
		}else{
			if ($track_exp{$miRNA_id}) {
				print OUT "$track_exp{$miRNA_id}\t$chr\t--\t--\t$strand\t--\t$start\t$end\t$strand\n";
			}else{
				print OUT "$miRNA_id\t$chr\t--\t--\t$strand\t--\t$start\t$end\t$strand\n";
			}
		}
        my ($mature_seq,$pre_seq)=(split /\t/,$_)[2,8];
		my $mature_seq_uc = uc $mature_seq;
		my $pre_seq_uc = uc $pre_seq;
		if ($track_exp{$miRNA_id}) {
			print OUT1 ">$track_exp{$miRNA_id}\n$mature_seq_uc\n";
			print OUT2 ">$track_exp{$miRNA_id}\n$pre_seq_uc\n";
		}else{
			print OUT1 ">$miRNA_id\n$mature_seq_uc\n";
			print OUT2 ">$miRNA_id\n$pre_seq_uc\n";
		}
    }
    close INN;
}
close OUT;
close OUT1;
close OUT2;

if (-e "$fOut/$miRNA_pos") {
    system "touch $fOut/miRNA_position.Check";
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
    #���б��е����ֵ
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
    #���б��е���Сֵ
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
    #��ȡ�ַ������еķ��򻥲����У����ַ�����ʽ���ء�ATTCCC->GGGAAT
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
     Contact:   linhj <linhj\@biomarker.com.cn> 
Program Date:   2017.03.15
      Modify:   
 Description:   This program is used to get miRNA position and miRNA mature seq and miRNA pre seq.
       Usage:   perl $Script
        Options:
        -i <dir>   input miRNA summary file,multi summary file sep by blank  : forced
        -k <file>  know miRNA position file,[eq: /share/nas1/linhj/bin/tools/know_miRNA_pos.list ]
        -o <dir>   output dir   : forced
        -e <file>  miRNA expression file : forced 

        -h      help

USAGE
    print $usage;
    exit;
}
