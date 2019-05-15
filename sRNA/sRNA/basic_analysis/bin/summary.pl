use strict;
use warnings;
use Getopt::Long;
use FindBin qw($Bin $Script);
use File::Basename qw(basename dirname);

#################################################
##summary the results produced by Free_Pineline##
#################################################
my ($imrd,$ilist,$icoords,$imature,$ipre,$iexpression);
my  ($od,$sample,$ref);
GetOptions(
	"help|?"  		=>\&USAGE,
	#"imrd:s"		=>\$imrd,			### output.mrd
	"ilist:s"		=>\$ilist,			### miRNA_expression.list,come frome output.mrd filtered
	#"icoords:s"		=>\$icoords,		### precursors.coords
	"iexpression:s"	=>\$iexpression,	## All_miRNA_expression.list
	"imature:s"		=>\$imature,		## *.mature.fa
	"ipre:s"		=>\$ipre,			## *.pre.fa
	#"istr:s"		=>\$istr,			## precursors.str
	"od:s"			=>\$od,
	"ref:s"			=>\$ref,
	) or &USAGE;
	
&USAGE unless ($od);
$ref||="no";
mkdir $od unless -d $od;
$od=&ABSOLUTE_DIR($od);


###############结果文件##############
#my $work_dir=&check_time("$re_dir/miRDeep2/mirdeep_runs/");
#my $output_mrd="$re_dir/miRDeep2/mirdeep_runs/$work_dir/output.mrd";

#my $icoords="$re_dir/miRDeep2/mirdeep_runs/$work_dir/tmp/precursors.coords";
#my $iexpression="$re_dir/miRDeep2/miRNA_Quantify/All_miRNA_expression.list";   #######quantified expression  
#my $target=(glob("$re_dir/Target_Prediction/*.mir2target.list"))[0];
#my $imature=(glob("$re_dir/Alignment/Ref_Database/*.mature.fa"))[0];
#my $energy="$re_dir/miRDeep2/mirdeep_runs/$work_dir/tmp/precursors.str";

#################################################
my ($mi_name,%score,%seq,%star_seq,%genomeID,%strand,%start,%end,%pre_seq,%pre_stru,%tpm,%target);
my %energy;
open(IN,$ilist)|| die "check whether $ilist exist!\n";
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

###mrd###
=cut
open(IN,$output_mrd)|| die "check whether $output_mrd exist!\n";
$/=">";<IN>;
while(<IN>){
	chomp;
	next if /^\s*$/;
	($mi_name)=$_=~/^(\S+)/;
	($score{$mi_name})=$_=~/score\s+total\s+(\S+)/;
	my ($pre_seq)=$_=~/\npri_seq\s+(\S+)/;
	$pre_seq{$mi_name}=$pre_seq;
	($pre_stru{$mi_name})=$_=~/\npri_struct\s+(\S+)\s*#MM/;
	my ($obs_seq)=$_=~/\nobs\s+(\S+)/;
	my ($exp_seq)=$_=~/\nexp\s+(\S+)/;
	my $map_seq=($obs_seq?$obs_seq:$exp_seq);
	my $M_num  =($map_seq=~s/M/M/g);
	#print "$map_seq\n$M_num\n";
	my $MMMMM  ="M"x$M_num;
	my $index  =index($map_seq,$MMMMM);
	#print "$index\n$pre_seq\n";
	$seq{$mi_name}=substr($pre_seq,$index,$M_num);
	my $S_num=($map_seq=~s/S/S/g);
	my $SSSSS="S"x$S_num;
	$index=index($map_seq,$SSSSS);
	$star_seq{$mi_name}=substr($pre_seq,$index,$S_num);
}
close IN;
$/="\n";
open(IN,"$icoords")||die "check whether $icoords exist!\n";
#my (%pre_str,%pre_start,%pre_end);
while(<IN>){
	my @line=split(/\s+/,$_);
	$line[0]=~s/>//g;
	$strand{$line[0]}=$line[1];
	$start{$line[0]}=$line[2];
	$end{$line[0]}=$line[3];
}
close IN;
=cut


open(IN,"$iexpression")|| die "check whether $iexpression exist!\n"; my %tpm_known;
$sample=<IN>;
my @sample=split/\s+/,$sample;
shift @sample;
$sample=join("\t",@sample);
while(<IN>){
	next if /^#/;
	my @line=split(/\s+/,$_);
	my $n=$line[0];
	if($n=~/conservative_/){
		$n=~s/^.*conservative_//g;
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

open(OUT,">$od/sum_predicted_miRNA.txt")||die $!;
print OUT "miRNA\tscore_total\tmature_sequence\tstar_sequence\tGenomeID\tstrand\tstart\tend\tpre_seq\thairpin_stru\thairpin_energy\t$sample\n";
foreach(keys %tpm){
	print OUT "$_\t";
	if (exists $score{$_}){print OUT "$score{$_}\t";} else{print OUT "\t";}
	if (exists $seq{$_}){print OUT "$seq{$_}\t";} else{print OUT "\t";}
	if(exists $star_seq{$_}){print OUT "$star_seq{$_}\t";} else {print OUT "\t";}
	if(/^(.+)_/){print OUT "$1\t";} else {print OUT "\t";}
	if(exists $strand{$_}){print OUT "$strand{$_}\t";} else {print OUT "\t";}
	if(exists $start{$_}){print OUT "$start{$_}\t";} else {print OUT "\t";}
	if(exists $end{$_}){print OUT "$end{$_}\t";} else {print OUT "\t";}
	if(exists $pre_seq{$_}){print OUT "$pre_seq{$_}\t";} else {print OUT "\t";}
	if(exists $pre_stru{$_}){print OUT "$pre_stru{$_}\t";} else {print OUT "\t";}
	if(exists $energy{$_}){print OUT "$energy{$_}\t";} else {print OUT "\t";}
	if(exists $tpm{$_}){
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

open(IN,"$imature")||die "check whether $imature exist!\n";
$/=">";<IN>;my %seq_known;
while(<IN>){
	chomp;
	my @line=split(/\n/,$_);
	$seq_known{$line[0]}=$line[1];
}
close IN;
$/="\n";

open(IN,$ipre) or die $!;
$/=">";
<IN>;
my %PRE_KNOWN;
while(<IN>){
	chomp;
	my ($id,$seq)=(split /\n+/,$_,2)[0,1];
	$id=~s/MIR/miR/;
	$seq=~s/\s+//g;
	$PRE_KNOWN{$id}=$seq;
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
	#if(exists $pre_seq{$_}){print OUT "$pre_seq{$_}\t";} else {print OUT "--\t";}
	if(exists $PRE_KNOWN{$_}){print OUT "$PRE_KNOWN{$_}\t";} else {print OUT "--\t";}
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
  -ilist       miRNA_expression.list , com from output.mrd filtered
  
  -imature     *.mature.fa
  
  -iexpression All_miRNA.expression.list   TPM file
  
  -ipre        *.pre.fa 
  
  -od          <dir>   output dir ,forced
  
  -h           Help
  
  -ref         <str>	"yes" if known mature exist; default "no";

USAGE
	print $usage;
	exit;
}
