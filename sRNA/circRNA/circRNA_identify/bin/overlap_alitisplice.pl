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
# ------------------------------------------------------------------
#
# ------------------------------------------------------------------
open (IN,$fIn) or die $!;

my %hash;
#$/=">";
my $i=1;
my $head=<IN>;
while (<IN>) {
    chomp;
    my @lines=split /\t/,$_;
    my $chr = (split /:/,$lines[0])[0];
    my $start_end = (split /:/,$lines[0])[1];
    my $start=(split /\|/,$start_end)[0];
    my $end=(split /\|/,$start_end)[1];

    my $gene = $lines[-1];
	$gene =~s/,$//;
    $gene = $chr if(/intergenic/);
	if($gene=~/,/)
	{
		my @gene_list = split/,/,$gene;
		foreach my $gene_id(@gene_list)
		{
			push @{ $hash{$gene_id} }, [$start, $end, $_];
		}
	}
	else
	{
		push @{ $hash{$gene} }, [$start, $end, $_];
	}
}

close IN;

my %cluster = ();
my $index = 0 ;
foreach  my $gene (keys %hash) {
    my @array =sort {$a->[0] <=> $b->[0]} @{$hash{$gene}} ;
    push @{$cluster{"$index"}}, $array[0][2] ;
    for (my $i=1;$i<= $#array;$i++) {
        if ($array[$i][0]-$array[$i-1][1]<=0) {
            push @{$cluster{"$index"}}, $array[$i][2] ;
        }else{
            $index++ ;
            push @{$cluster{"$index"}}, $array[$i][2] ;
        }
    }
$index++ ;
}
open (OUT,">$fOut") or die $!;
print OUT "$head";
foreach my $key (sort {$a<=>$b} keys %cluster) {
    if (@{$cluster{$key}} >1 ) {
        foreach my $temp (@{$cluster{$key}}) {
            print OUT "$temp\n";
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
        -i <file>   input file,xxx format,forced

        -o <file>   output file,optional

        -h      help

USAGE
    print $usage;
    exit;
}
