#!/usr/bin/env perl -w 
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
my ($deg,$ce,$o);

GetOptions(
				"help|?" =>\&USAGE,
				"o:s"=>\$o,
				"deg:s"=>\$deg,
				"ce:s"=>\$ce,
				) or &USAGE;
&USAGE unless ($deg and $ce and  $o);
# ------------------------------------------------------------------
# 
# ------------------------------------------------------------------
open (IN,$deg) or die $!;
open (INN,$ce) or die $!;
open (OUT,">$o") or die $!;

my $name = basename $o;
my $dir = dirname $o;
my ($type1,$group,$type,$other)=split /\./,$name,4;
open (Inter,">$dir/$type1.$group.$type.Interaction.list");
open (Att,">$dir/$type1.$group.$type.Attribution.list");
my %hash;
while (<IN>) {
	chomp;
	if (/^#/) {next;}
	my @tmp =split /\t/,$_;
	$hash{$tmp[0]}=$tmp[-1];
}
#ceRNA1 ceRNA2  ceRNA1_num      ceRNA2_num      overlap_num     ceRNA1_miRNA    ceRNA2_miRNA    overlap_miRNA
my %h;
while (<INN>) {
	chomp;
	
	if(/^#/){
		my @B = split /\t/,$_;
		my $title = join("\t",@B[0..7]);
		print OUT "$title\n";
		next;
	}
	my $co_exp;
	my @A = split /\t/,$_;
	my ($t1,$id1)=split /:|;/,$A[0],2;
	my ($t2,$id2)=split /:|;/,$A[1],2;
	my @a = split /\;|,/,$A[5];
	my @b = split /\;|,/,$A[6];
	my @c = split /\;|,/,$A[7];
	my $list = join("\t",$A[0],$A[1]);
	if(exists $hash{$id1} && exists $hash{$id2}){
		my $c1=0;my $c2=0;my $c3=0;my(%h1,%h2,%h3);
		my $flag =0;
		for (my $i=0;$i<@a;$i++){
			if(exists $hash{$a[$i]}){
				$flag++;
				$c1++;
				if($flag ==1){
					$h1{$list} = $a[$i];
				}else{
					$h1{$list} .= ";".$a[$i];
				}
			}
		}
		if($c1==0){$h1{$list}="--";}
		$flag =0;
        	for (my $j=0;$j<@b;$j++){
                	if(exists $hash{$b[$j]}){
                        	$flag++;
	                        $c2++;
        	                if($flag ==1){
                	                $h2{$list} = $b[$j];
                        	}else{
                                	$h2{$list} .= ";".$b[$j];
	                        }
        	        }
        	}
		if($c2==0){$h2{$list}="--";}
		$flag =0;
        	for (my $z=0;$z<@c;$z++){
                	if(exists $hash{$c[$z]}){
	                        $flag++;
        	                $c3++;
                	        if($flag ==1){
                        	        $h3{$list} = $c[$z];
	                        }else{
        	                        $h3{$list} .= ";".$c[$z];
                	        }
	                }
        	}
		
		if ($c3>0){
			my $tc1 = $c1+$c3;
			my $tc2 = $c2+$c3;
			print OUT "$list\t$tc1\t$tc2\t$c3\t$h1{$list}\t$h2{$list}\t$h3{$list}\n";
			print Inter "$id1\t$id2\n";
			if (!exists $h{$id1}){
				$h{$id1}=1;
				print Att "$id1\t$hash{$id1}\t$t1\n";
			}
			if (!exists $h{$id2}){
				$h{$id2}=1;
				print Att "$id2\t$hash{$id1}\t$t2\n";
			}
			my @mi = split /;|,/,$h3{$list};
			for(my $i=0;$i<@mi;$i++){
				if(exists $hash{$mi[$i]}){
					print Inter "$id1\t$mi[$i]\n";
					if(!exists $h{$mi[$i]}){
						$h{$mi[$i]}++;
						print Att "$mi[$i]\t$hash{$mi[$i]}\tmiRNA\n";
					}
					print Inter "$id2\t$mi[$i]\n";
				}
			}
		}
	}
}
close IN;
close INN;
close OUT;

#######################################################################################
print STDOUT "\nDone. Total elapsed time : ",time()-$BEGIN_TIME,"s\n";
#######################################################################################

# ------------------------------------------------------------------
# sub function
# ------------------------------------------------------------------
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
     Version:	$version
     Contact:	Simon Young <yangxh\@biomarker.com.cn> 
Program Date:	2012.07.02
      Modify:	
 Description:	This program is used to ......
       Usage:
		Options:
		-deg <file>	*_vs_*.all.DEG.xls
		-ce <file>	ceRNA_pair_adjust_p_Sig.txt or Coexpression_ceRNA_pair_adjust_p_Sig.txt
		-o <file>	*_vs_*.lncRNA-miRNA-mRNA.xls
		-h		help

USAGE
	print $usage;
	exit;
}
