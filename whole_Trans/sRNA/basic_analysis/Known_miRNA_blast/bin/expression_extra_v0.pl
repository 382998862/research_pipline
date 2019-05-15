#!usr/bin/perl -w
use strict;
use warnings;
use Getopt::Long;
use Data::Dumper;
use FindBin qw($Bin $Script);
use File::Basename qw(basename dirname);
my $BEGIN_TIME=time();
my $version="1.0.0";
##***********
#This script is used to parser .mrd file.
#4 files will be generated.
#1.mature.mrd
#2.novel_conservative.mrd
#3.novel_unconservative.mrd
#4.miRNA_expression.list
#************
my($exp_file,$fa_file,$out_exp,$out_fa,$pre,$out_pre);
GetOptions(
				"help|?"	=>\&USAGE,
				"exp:s"		=>\$exp_file,
				"fa:s"		=>\$fa_file,
				"pre:s"		=>\$pre,	
				"oexp:s"	=>\$out_exp,
				"ofa:s"		=>\$out_fa,
				"opre:s"	=>\$out_pre,
				) or &USAGE;
&USAGE unless ($exp_file and $fa_file and $pre and $out_exp and $out_fa and $out_pre);
open EXP, "<$exp_file";	    			######### csv file
open FA , "<$fa_file";      			######### All_miRNA.fa
open OUTEXP, ">$out_exp";				######### All_miRNA.expression.list
open (PRE,$pre) or die $!;				######### All_miRNA_Pre.fa

############################################## FA Files
my %gene;
$/=">";
<FA>;
while (<FA>) {
	chomp;
	my ($name,$seq) =(split/\n/,$_);
	$gene{$name} = $seq;
}
$/="\n";
close FA;

my %PRE;
$/=">";
<PRE>;
while(<PRE>){
	chomp;
	my ($pre_id,$pre_seq)=(split /\n/,$_)[0,1];
	$pre_seq=~s/\s+$//;
	$PRE{$pre_id}=$pre_seq;
}
$/="\n";
close PRE;


##############################################表达量输出

my $head = <EXP>;
chomp $head;
my @head=(split/\t/,$head);
my @index;
my @samples ;

#############################################去除样品名限制，修改从表头读取样品名方法
#表头格式：
##miRNA	read_count	precursor	total	S01	S02	S01(norm)	S02(norm)
##miRNA	read_count	precursor	total	S01	S02	S03	S04	S01(norm)	S02(norm)	S03(norm)	S04(norm)

for (my $i=4;$i<@head;$i++) {
	unless ($head[$i] =~ /\(norm\)$/) {
		push @index ,$i;
		push @samples ,$head[$i];
	}
}

push  @samples ,'geneLength';
unshift @samples,$head[0];
unshift @index ,0;

print OUTEXP  join "\t",@samples ;
print OUTEXP "\n";

my %miRNA;
my %PRE_exp;
while (<EXP>) {
	chomp;
	next if (/^#/);
	my @line = split/\t/, $_;
	my $miRNAName = $line[0];
	$miRNAName=~s/\s+$//;
	if($line[1]>0){
		my $string=shift;
		my $exp_total = 0;
		foreach my $index (@index) {
			if ($index == 0) {
				$string .= "$line[$index]\t";
			}
			else {
				my $exp_int = int($line[$index]);
				$exp_total+=$exp_int;
				$string .= "$exp_int\t";
			}
		}
		next if ($exp_total==0);

		my $genelen = length $gene{$miRNAName};
		$string .=$genelen;
		$miRNA{$miRNAName}=$string;
		$PRE_exp{$line[2]}++;
	}
}
close EXP;


open OUTFA, ">$out_fa";					######### All_miRNA.expressed.fa
foreach my $mat(sort keys %miRNA){
	print OUTFA ">$mat\n$gene{$mat}\n";
	print OUTEXP "$miRNA{$mat}\n";
}
close OUTFA;
close OUTEXP;

open (OUTPRE,">$out_pre") or die $!;	######### All_miRNA_Pre.expressed.fa
foreach my $pre(sort keys %PRE_exp){
	print OUTPRE ">$pre\n$PRE{$pre}\n";
}
close OUTPRE;

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

sub USAGE {#
	my $usage=<<"USAGE";
Program:	$0
Version:	$version
Usage:
  Options:

	-exp  <str> input exp file,miRNAs_expressed_all_samples_*.csv
	-fa   <str> All_miRNA.fa ,forced 
	-pre  <str> All_miRNA_Pre.fa,forced
	-oexp <str> expression filter file ,All_miRNA.count.list
	-ofa  <str> expressed miRNA fa file,All_miRNA.expressed.fa
	-opre <str> All_miRNA_Pre.expressed.fa
	-h          Help

USAGE
	print $usage;
	exit;
}
