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
my ($fIn,$o,$label,$mappedReads,$readLength);
GetOptions(
				"help|?" =>\&USAGE,
				"i:s"=>\$fIn,
				"o:s"=>\$o,
				"m:s"=>\$mappedReads,
				"len:s"=>\$readLength,
				"l:s"=>\$label
				) or &USAGE;
&USAGE unless ($fIn and $o);

$label||='SRPBM';

my %mapped;
my %length;
if (defined $mappedReads) {
    open(M, $mappedReads) or die $!;
	while (<M>) {
		next if ($.==1);
		my @line = split /\t/,$_;
		my $mapped = $line[2];
		$mapped=~s/,//g;
		$mapped=~s/\(\S+\)//g;
		$mapped{$line[0]}=$mapped;
	}
	close M;
}

if(defined $readLength)
{
	open(M, $readLength) or die $!;
	 while (<M>) {
		chomp;
                my @line = split /\t/,$_;
                $length{$line[0]}=$line[1];
        }
        close M;
}


my %H;
my %FPKM;
open (IN,$fIn) or die $!;
my $head_line=<IN>;
my @Samples=split/\t/,$head_line;
shift @Samples;
pop @Samples;
while (<IN>) {
	my @A=split/\t/,$_;
	shift @A;
	pop @A;
	for (my $i=0;$i<@A ;$i++) {
		$H{$Samples[$i]}+=$A[$i];
	}
}
close (IN) ;

my $line=join "\t",@Samples;
my @sam;
foreach my $s (@Samples) {
	push @sam,$s."_count";
	push @sam,$s."_$label";
}
my $exp_line=join "\t",@sam;

open (OUT,">$o/All_gene_expression.list") or die $!;
open EXP,">$o/All_gene_expression_detail.list" or die $!;

print OUT "#ID\t$line\n";
print EXP "\#ID\t$exp_line\n";
open (IN,$fIn) or die $!;
<IN>;
my @total;
while (<IN>) {
	chomp;
	my @A=split/\t/,$_;
	my $name=shift @A;
	my $len=pop @A;
	my @C;
	my @B;
	for (my $i=0;$i<@A;$i++) {
		if ($label eq 'TPM') {
			$C[$i]="$A[$i]\t".EXPRESSION($A[$i],$H{$Samples[$i]},$len,$label);
            		$FPKM{$Samples[$i]} += EXPRESSION($A[$i],$H{$Samples[$i]},$len,$label);
        	}
		else
		{
			my $value;
			if($label eq 'FPKM'){$value = EXPRESSION($A[$i],$H{$Samples[$i]},$len,$label);}
			elsif($label eq 'SRPBM'){$len=$length{$Samples[$i]};$value = EXPRESSION($A[$i],$mapped{$Samples[$i]},$len,$label);}
			$B[$i]="$A[$i]\t$value";
			$A[$i]=$value;
		}
	}
	push @C,$name;
	push @total,[@C];
	if ($label ne 'TPM') {  
		$line=join "\t",@A;
		$exp_line=join "\t",@B;
		print OUT "$name\t$line\n";
		print EXP "$name\t$exp_line\n";
    }
}
close IN;
if ($label eq 'TPM') {
	for(my $i=0;$i<@total;$i++)
	{
		print OUT "$total[$i][-1]";
		print EXP "$total[$i][-1]";
		for(my $j=0;$j<@{$total[$i]}-1;$j++)
		{
			my ($count,$count_nor)=split /\t/,$total[$i][$j],2;
			if ($FPKM{$Samples[$j]}==0) {
				print OUT "\t0";
				print EXP "\t$count\t0";
			}
			else
			{
				my $TPM = $count_nor*1000000/$FPKM{$Samples[$j]};
				print OUT "\t$TPM";
				print EXP "\t$count\t$TPM";
			}
		}
		print OUT "\n";
		print EXP "\n";
	}
}

close OUT;
close EXP;


#######################################################################################
print STDOUT "\nDone. Total elapsed time : ",time()-$BEGIN_TIME,"s\n";
#######################################################################################

# ------------------------------------------------------------------
# sub function
# ------------------------------------------------------------------
################################################################################################################


sub EXPRESSION
{
	my($value,$total,$length,$method)=@_;
	my $expression;
	if ($method eq 'TPM') {
        	$expression=$value*1000/$length;
   	}	
	elsif($method eq 'FPKM')
	{
		$expression=$value/($length/1000)/($total/1000000);
	}
	elsif($method eq 'SRPBM')
    	{
		$expression=$value*1000000000/$total/$length;
	}
	return $expression;
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


sub USAGE {#
	my $usage=<<"USAGE";
ProgramName:
Version:	$version
Contact:	Zhang XueChuan <zhangxc\@biomarker.com.cn> 
Program Date:   2013.10.11
Usage:
  Options:
  -i    <file>  input file,gene count file,forced 
  -l    normalization method(TPM,FPKM,SRPBM)
  -m    All.mappedStat.xls file,if use  SRPBM normalization,must be given
  -o    <file>  out direction for gene expression file,forced
  -len  <file>  sample read length stat 
  -h         Help

USAGE
	print $usage;
	exit;
}
