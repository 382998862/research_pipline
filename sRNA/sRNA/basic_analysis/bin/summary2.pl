#!/usr/bin/perl
use strict;
use warnings;
use Getopt::Long;
use FindBin qw($Bin $Script);
use File::Basename qw(basename dirname);

#################################################
##summary the results produced by Free_Pineline##
#################################################

my  ($re_dir,$od,$sample,$ref);
GetOptions(
	"help|?"   =>\&USAGE,
	"re_dir:s" =>\$re_dir,
	"od:s"     =>\$od,
	"ref:s"    =>\$ref,
	) or &USAGE;

&USAGE unless ($re_dir and $od);
$ref||="no";
mkdir $od unless -d $od;
$re_dir=&ABSOLUTE_DIR($re_dir);
$od=&ABSOLUTE_DIR($od);


###############结果文件##############
my $work_dir=&check_time("$re_dir/miRDeep2/mirdeep_runs/");
my $output_mrd="$re_dir/miRDeep2/mirdeep_runs/$work_dir/output.mrd";
my $exp_list="$re_dir/miRDeep2/Result/miRNA_expression.list";         #######come frome output.mrd filtered
my $name_list="$re_dir/miRDeep2/Result/chname.list"; # yaling at 17-12-14
my $pre_coords="$re_dir/miRDeep2/mirdeep_runs/$work_dir/tmp/precursors.coords";
my $expression="$re_dir/miRDeep2/miRNA_Quantify/All_miRNA_expression.list";   #######quantified expression
#my $target=(glob("$re_dir/Target_Prediction/*.mir2target.list"))[0];
my $mature=(glob("$re_dir/../sRNA_Alignment/Ref_Database/*.mature.fa"))[0];
my $energy="$re_dir/miRDeep2/mirdeep_runs/$work_dir/tmp/precursors.str";

#################################################
my ($mi_name,%score,%seq,%star_seq,%genomeID,%strand,%start,%end,%pre_seq,%pre_stru,%target);
my %tpm;
my %energy;
open(IN,$exp_list)|| die "check whether $exp_list exist!\n";
<IN>;
while(<IN>){
	chomp;
	my @temp=split/\s+/,$_;
	$mi_name=$temp[1];
	my @temp1=split/,/,$temp[2];
	if($#temp1>0){
		foreach(@temp1){
			my $mi_name=$_;
			$score{$mi_name}=$temp[7];
			$pre_seq{$mi_name}=$temp[6];
			$pre_stru{$mi_name}=$temp[17];
			$seq{$mi_name}=$temp[5];
			$star_seq{$mi_name}=$temp[13];
			$strand{$mi_name}=$temp[14];
			$start{$mi_name}=$temp[15];
			$end{$mi_name}=$temp[16];
			$energy{$mi_name}=$temp[18];
		}
	}
	$score{$mi_name}=$temp[7];
	$pre_seq{$mi_name}=$temp[6];
	$pre_stru{$mi_name}=$temp[17];
	$seq{$mi_name}=$temp[5];
	$star_seq{$mi_name}=$temp[13];
	$strand{$mi_name}=$temp[14];
	$start{$mi_name}=$temp[15];
	$end{$mi_name}=$temp[16];
	$energy{$mi_name}=$temp[18];
}

my %name;
open(LIST, $name_list) || die $!;
while (<LIST>) {
	next if /^\s+$/;
	my @a = split;
	$name{$a[1]} = $a[0];
}
close LIST;

##All_miRNA_expression.list
open(IN,"$expression")|| die "check whether $expression exist!\n";
my %tpm_known;
$sample=<IN>;
my @sample=split/\s+/,$sample;
shift @sample;
$sample=join("\t",@sample);
while(<IN>){
	next if /^#/;
	my @line=split(/\s+/,$_);
	my $n=$line[0];
	if ($n =~ /novel/) {
	# if($n=~/conservative_/){
		# $n=~s/^.*conservative_//g;
		shift @line;
		@{$tpm{$n}}=@line;
	}
	else{
		shift @line;
		@{$tpm_known{$n}}=@line;
	}

}
close IN;

=c
open(IN,"$target")||die "check whether $target exist!\n";
<IN>;
while(<IN>){
	my @line=split(/\s+/,$_);
	$line[0]=~s/unconservative_//g;
	$line[0]=~s/conservative_//g;
	$target{$line[0]}=$line[1];
}
close IN;
=cut

=cut
open(IN,"$energy")||die "check whether $energy exist!\n";
$/=">";<IN>; my %energy;
while(<IN>){
	chomp;
	my @line=split(/\n/,$_);
	#my $len=length($line[2]);
	if($line[2]=~/.+\((.+)\)/){
		$energy{$line[0]}=$1;
	}
}
close IN;
=cut


open(OUT,">$od/sum_predicted_miRNA.txt")||die $!;
print OUT "miRNA\tscore_total\tmature_sequence\tstar_sequence\tGenomeID\tstrand\tstart\tend\tpre_seq\thairpin_stru\thairpin_energy\t$sample\n";
foreach (keys %tpm){
	print OUT "$_\t";
	my $old_name = $name{$_};
	if (exists $score{$old_name}) {print OUT "$score{$old_name}\t";} else {print OUT "\t";}
	if (exists $seq{$old_name}) {print OUT "$seq{$old_name}\t";} else {print OUT "\t";}
	if (exists $star_seq{$old_name}) {print OUT "$star_seq{$old_name}\t";} else {print OUT "\t";}
	if ($old_name =~ /^(.+)_/) {print OUT "$1\t";} else {print OUT "\t";}
	if (exists $strand{$old_name}) {print OUT "$strand{$old_name}\t";} else {print OUT "\t";}
	if (exists $start{$old_name}) {print OUT "$start{$old_name}\t";} else {print OUT "\t";}
	if (exists $end{$old_name}) {print OUT "$end{$old_name}\t";} else {print OUT "\t";}
	if (exists $pre_seq{$old_name}) {print OUT "$pre_seq{$old_name}\t";} else {print OUT "\t";}
	if (exists $pre_stru{$old_name}) {print OUT "$pre_stru{$old_name}\t";} else {print OUT "\t";}
	if (exists $energy{$old_name}) {print OUT "$energy{$old_name}\t";} else {print OUT "\t";}
	if (exists $tpm{$_}){
		my $tpm=join("\t",@{$tpm{$_}});
		print OUT "$tpm\t";
	}
	else{
		my @t=split(/\t/,$sample);
		#my @t=@sample;
		for(my $i=0;$i<=$#t;$i++){
			print OUT "\t";
		}
	}
	#if(exists $target{$_}){print OUT "$target{$_}";} else {print OUT "none";}
	print OUT "\n";
}
close OUT;
###################################################################

open(IN,"$mature")||die "check whether $mature exist!\n";
$/=">";<IN>;my %seq_known;
while(<IN>){
	chomp;
	my @line=split(/\n/,$_);
	$seq_known{$line[0]}=$line[1];
}
close IN;
$/="\n";


if($ref eq "yes"){
open(OUT,">$od/sum_known_miRNA.txt")||die $!;
print OUT "miRNA\tscore_total\tmature_sequence\tstar_sequence\tGenomeID\tstrand\tstart\tend\tpre_seq\thairpin_stru\thairpin_energy\t$sample\n";
foreach(keys %tpm_known){
	print OUT "$_\t";
	if (exists $score{$_}){print OUT "$score{$_}\t";} else{print OUT "--\t";}

	if (exists $seq_known{$_}){print OUT "$seq_known{$_}\t";} else{print OUT "--\t";}

	if(exists $star_seq{$_}){print OUT "$star_seq{$_}\t";} else {print OUT "--\t";}
	if(/^(.+)_/){print OUT "$1\t";} else {print OUT "--\t";}
	if(exists $strand{$_}){print OUT "$strand{$_}\t";} else {print OUT "--\t";}
	if(exists $start{$_}){print OUT "$start{$_}\t";} else {print OUT "--\t";}
	if(exists $end{$_}){print OUT "$end{$_}\t";} else {print OUT "--\t";}
	if(exists $pre_seq{$_}){print OUT "$pre_seq{$_}\t";} else {print OUT "--\t";}
	if(exists $pre_stru{$_}){print OUT "$pre_stru{$_}\t";} else {print OUT "--\t";}
	if(exists $energy{$_}){print OUT "$energy{$_}\t";} else {print OUT "--\t";}
	if(exists $tpm_known{$_}){
		my $tpm=join("\t",@{$tpm_known{$_}});
		print OUT "$tpm\t";
	}
	else{
		my @t=split("\t",$sample);
		for(my $i=0;$i<=$#t;$i++){
			print OUT "\t";
		}
	}
#	if(exists $target{$_}){print OUT "$target{$_}";} else {print OUT "none";}
	print OUT "\n";
}
close OUT;
}


##########################sub####################
sub check_time{
	my $dir=shift;
	my $check;
	my @Old_time=(0,0,0,0,0,0);
	foreach my $d (glob "$dir/run*") {
		next unless -d $d;
		my $run_name=basename $d;
		my @Time=(split/_/,$run_name)[3,2,1,5,6,7];
		foreach my $i (0..5) {
			if ($Time[$i]>$Old_time[$i]) {
				@Old_time=@Time;last;
			}
		}
	}
	return "run_$Old_time[2]_$Old_time[1]_$Old_time[0]_t_$Old_time[3]_$Old_time[4]_$Old_time[5]";
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


sub USAGE{
	my $usage=<<"USAGE";
ProgramName:
Usage:
  Options:
  -re_dir		<dir>   result dir produced by Free_Pineline ,forced

  -od			<dir>   output dir ,forced

  -h			Help

  -ref			<str>	"yes" if known mature exist; default "no";

USAGE
	print $usage;
	exit;
}
