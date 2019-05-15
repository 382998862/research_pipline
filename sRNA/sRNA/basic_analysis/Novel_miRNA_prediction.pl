#!/usr/bin/perl -w
use strict;
use Getopt::Long;
use Data::Dumper;
use FindBin qw($Bin $Script);
use File::Basename qw(basename dirname);
use Cwd qw(abs_path);
use newPerlBase;
my $BEGIN_TIME=time();
my $version="1.0";

my($unquantified_reads,$idir,$main_cfg,$min,$max,$step,$sbs,$genome,$plant,$max_map);
GetOptions(
	"help|?"	=>\&USAGE,
	"uq:s"		=>\$unquantified_reads,
	"idir:s"	=>\$idir,
	"main_cfg:s"	=>\$main_cfg,
	"min:s"		=>\$min,
	"max:s"		=>\$max,
	"step:s"	=>\$step,
	"sbs"		=>\$sbs,
) or &USAGE;
&USAGE unless ( (defined $unquantified_reads and -f $unquantified_reads) and $idir and $main_cfg);
###############################################
#&USAGE unless (defined $fyes or defined $fno); #This script, with no -yes and -no options, only predict novel miRNAs. (The Known miRNAs should be identified before novel miRNAs prediction, by blastn against miRBase.)
###############################################

################default values
$step ||=0;
$min ||=15;
$max ||=35;
$plant||=0;
############### read CFG
my %CFG=%{readconf("$Bin/../CFG")};
my $miRBase = $CFG{miRBase};

$main_cfg=abs_path($main_cfg);
$idir=abs_path($idir);
my $miRDeep2_path="$idir/miRDeep2";
mkdir $miRDeep2_path unless (-d $miRDeep2_path);

#################read main_cfg
my %Raw_FQ;

open(IN,$main_cfg) or die $!;
while(<IN>){
	chomp;
	next if(/^#/||/^$/);
	if(/^NAME/){
		my $name=(split/\s+/,$_)[1];
		my $fq=<IN>;
		$Raw_FQ{$name}=(split/\s+/,$fq)[1];
		unless (-e $Raw_FQ{$name}) {
			print "Check Your config file:$main_cfg!Expecially $name\n\n";die;
		}
	}
	if(/^SPECIES_TYPE/){
		my @tmp=split(/\s+/,$_);
		$plant=$tmp[1];		
	}
	if(/^max_map/){
		my @tmp=split(/\s+/,$_);
		$max_map=$tmp[1];
	}
	if(/^GENOME/){
		my @tmp=split(/\s+/,$_);
		$genome=$tmp[1];
	}
}
close IN;

$max_map||=5;
my $genome_name=basename($genome);
$genome="$idir/../sRNA_Alignment/Ref_Database/$genome_name";

################## divide unquantified reads
if($step==0){
	stepStart(0,"Divide unquantified reads");
	runOrDie("perl $Bin/bin/divide_unquantified_reads.pl -uq $unquantified_reads -od $idir/miRDeep2 ");
	stepTime(0);
	$step++ unless ($sbs);
}

################# miRDeep2 mapper.pl
if($step==1){
	stepStart(1,"miRDeep2 mapper.pl");
	unless((-f "$idir/miRDeep2/Total_reads_collapsed.fa") || (-f "$idir/miRDeep2/Total_reads_collapsed_vs_genome.arf")){
		my %FQ;
		open (OUT,">$idir/miRDeep2/mapper.cfg") or die $!;
		foreach my $name (keys %Raw_FQ){
			die "[Error]: There is no file: $idir/miRDeep2/$name.mapper.fa.\n" unless (-f "$idir/miRDeep2/$name.mapper.fa");
			print OUT "$name.mapper.fa\t$name\n";
		}
		close OUT;
		runOrDie("cd $idir/miRDeep2 && perl $Bin/mirdeep_animal/mapper.pl $idir/miRDeep2/mapper.cfg -r $max_map -d -c -j -l $min -m -p $genome -s Total_reads_collapsed.fa -t Total_reads_collapsed_vs_genome.arf -o 6 -v -n 2>mapper.out");
	}

	stepTime(1);
	$step++ unless($sbs);
}

#################  miRDeep2 miRDeep2.pl
if($step==2){
	stepStart(2,"miRDeep2 miRDeep2.pl");

	if($plant==0){
		runOrDie("cd $idir/miRDeep2 && perl $Bin/mirdeep_animal/miRDeep2.pl Total_reads_collapsed.fa $genome Total_reads_collapsed_vs_genome.arf none none none -d -g -1 2> report.log");
	}
	elsif($plant==1){
		runOrDie("cd $idir/miRDeep2 && perl $Bin/mirdeep_plant/miRDeep2.pl Total_reads_collapsed.fa $genome Total_reads_collapsed_vs_genome.arf none none none -g 50000 -l 250 -d -m 10 -v -P -n d 2> report.log");
	}
	else{
		print "wrong parameter! plant=$plant \n";
		die $!;
	}

	stepTime(2);
	$step++ unless ($sbs);
}

#################### miRNA Filter
if($step==3){
	stepStart(3,"miRNA Filter");
	my $work_dir=&check_time("$idir/miRDeep2/mirdeep_runs");
	my $seq=join ",",keys %Raw_FQ;
	$/=">";
	open(IN,"$idir/miRDeep2/mirdeep_runs/$work_dir/output.mrd") or die $!;
	<IN>;
	my $count;
	while(<IN>){
		$count++;
	}
	$/="\n";
	my $parse_version;
	if($count<=1000){$parse_version="parser_mrd_v4.pl";}
	else{$parse_version="parser_mrd_v5.pl";}

	runOrDie("perl $Bin/bin/$parse_version -str $idir/miRDeep2/mirdeep_runs/$work_dir/tmp/precursors.str -coord $idir/miRDeep2/mirdeep_runs/$work_dir/tmp/precursors.coords -mrd $idir/miRDeep2/mirdeep_runs/$work_dir/output.mrd -samples $seq -od $idir/miRDeep2/Result");

	stepTime(3);
	$step++ unless ($sbs);
}

##############extract miRNA
if($step==4){
	stepStart(4,"extract miRNA");
	my %MIR;
	open(IN,"$idir/miRDeep2/Result/miRNA_expression.list ") or die $!;	#######score total >=0 (output.mrd)
	my $head=<IN>;
	chomp $head;
	$head=(split /\t/,$head,14)[-1];
	$head=~s/<mature read count:>//g;
	$head=~s/>//g;
	while(<IN>){
		chomp;
		my @Info=split /\t/,$_,14;
		if($Info[2] ne 'novel'){
			next if exists $MIR{Known}{$Info[1]};
			$MIR{Known}{$Info[1]}{len}=$Info[4];
			$MIR{Known}{$Info[1]}{mature}=$Info[5];
			$MIR{Known}{$Info[1]}{precursor}=$Info[6];
		}
		else{
			$MIR{novel}{$Info[1]}{len}=$Info[4];
			$MIR{novel}{$Info[1]}{mature}=$Info[5];
			$MIR{novel}{$Info[1]}{precursor}=$Info[6];
		}
	}
	close IN;
	open (Nmature,">$idir/miRDeep2/Result/Novel_mature.fa") or die $!;
	open (Npre,">$idir/miRDeep2/Result/Novel_Pre.fa") or die $!;
	open (INFO,">$idir/miRDeep2/Result/chname.list") or die $!;
	my $count = 0;
	foreach my $key (sort keys %{$MIR{novel}}){
		$count ++;
		print INFO "$key\tnovel_miR_$count\n";
		print Nmature ">novel_miR_$count\n$MIR{novel}{$key}{mature}\n";
		print Npre ">novel_miR_$count\n$MIR{novel}{$key}{precursor}\n";
	}
	close Nmature;
	close Npre;
	close INFO;

	runOrDie("cat $idir/miRDeep2/Result/Novel_mature.fa >$idir/miRDeep2/Result/All_miRNA.fa");
	runOrDie("cat $idir/miRDeep2/Result/Novel_Pre.fa >$idir/miRDeep2/Result/All_miRNA_Pre.fa");
	stepTime(4);
	$step++ unless ($sbs);
}

######################## Novel miRNA Expression Quantify
if($step==5){
	stepStart(5,"Novel miRNA Expression Quantify");
	my $Expression="$idir/miRDeep2/Novel_miRNA_Quantify";
	mkdir $Expression unless (-d $Expression);
	if(-d "$idir/miRDeep2/Novel_miRNA_Quantify/expression_analyses"){
		&cmd_call("rm -r $idir/miRDeep2/Novel_miRNA_Quantify/expression_analyses");
		&cmd_call("rm $idir/miRDeep2/Novel_miRNA_Quantify/miRNAs_expressed_all_samples_*.csv");
	}
	runOrDie("cd $Expression && perl $Bin/mirdeep_animal/quantifier.pl -p $idir/miRDeep2/Result/All_miRNA_Pre.fa -m $idir/miRDeep2/Result/All_miRNA.fa -r $idir/miRDeep2/Total_reads_collapsed.fa -d -g 1 -e 2 -f 5 -W ");
	my $expression = (glob "$idir/miRDeep2/Novel_miRNA_Quantify/miRNAs_expressed_all_samples_*.csv")[0];
	runOrDie("perl $Bin/bin/expression_extra_v0.pl -exp $expression -fa $idir/miRDeep2/Result/All_miRNA.fa -pre $idir/miRDeep2/Result/All_miRNA_Pre.fa -oexp $idir/miRDeep2/Novel_miRNA_Quantify/All_Novel_miRNA.count.list -ofa $idir/miRDeep2/Novel_miRNA_Quantify/All_Novel_miRNA.expressed.fa -opre $idir/miRDeep2/Novel_miRNA_Quantify/All_Novel_miRNA_Pre.expressed.fa");
	stepTime(5);
	$step++ unless ($sbs);
}



######################## All miRNA Expression Combination
if($step==6){
	stepStart(6,"All miRNA Expression Combination");

	# keep final miRNA Quantify directory and necessary files same as original directory structure
	# v
	my $Expression="$idir/miRDeep2/miRNA_Quantify";
	mkdir $Expression unless (-d $Expression);
	mkdir "$Expression/expression_analyses" unless (-d "$Expression/expression_analyses");
	my $novel_csv = (glob "$idir/miRDeep2/Novel_miRNA_Quantify/expression_analyses/expression_analyses_*/miRNA_expressed.csv")[0];
	my $known_csv = (glob "$idir/Known_miRNA_identify/expression_analyses/expression_analyses_*/miRNA_expressed.csv")[0];
	my $novel_mature_converted = (glob "$idir/miRDeep2/Novel_miRNA_Quantify/expression_analyses/expression_analyses_*/mature.converted")[0];
	my $known_mature_converted = (glob "$idir/Known_miRNA_identify/expression_analyses/expression_analyses_*/mature.converted")[0];
	my $novel_precursor_converted = (glob "$idir/miRDeep2/Novel_miRNA_Quantify/expression_analyses/expression_analyses_*/precursor.converted")[0];
	my $known_precursor_converted = (glob "$idir/Known_miRNA_identify/expression_analyses/expression_analyses_*/precursor.converted")[0];
	my $novel_mature2hairpin = (glob "$idir/miRDeep2/Novel_miRNA_Quantify/expression_analyses/expression_analyses_*/mature2hairpin")[0];
	my $known_mature2hairpin = (glob "$idir/Known_miRNA_identify/expression_analyses/expression_analyses_*/mature2hairpin")[0];
	my $novel_expression = (glob "$idir/miRDeep2/Novel_miRNA_Quantify/miRNAs_expressed_all_samples_*.csv")[0];
	my $known_expression = (glob "$idir/Known_miRNA_identify/miRNAs_expressed_all_samples_*.csv")[0];
	my $novel_mrd = (glob "$idir/miRDeep2/Novel_miRNA_Quantify/expression_analyses/*/miRBase.mrd")[0];
	my $known_mrd = (glob "$idir/Known_miRNA_identify/expression_analyses/*/miRBase.mrd")[0];
	my $novel_arf = (glob "$idir/miRDeep2/Novel_miRNA_Quantify/expression_analyses/expression_analyses_*/All_miRNA_mapped.arf")[0];
	my $known_arf = (glob "$idir/Known_miRNA_identify/expression_analyses/expression_analyses_*/Known_Mature_miRNA_mapped.arf")[0];
	die "At least one of these files do not exists:\n$novel_csv\n$known_csv\n$novel_mature_converted\n$known_mature_converted\n$novel_precursor_converted\n$known_precursor_converted\n$novel_mature2hairpin\n$known_mature2hairpin\n$novel_expression\n$known_expression\n$novel_mrd\n$known_mrd\n$novel_arf\n$known_arf\n" unless (-f $novel_csv and -f $known_csv and -f $novel_mature_converted and -f $known_mature_converted and -f $novel_precursor_converted and -f $known_precursor_converted and -f $novel_mature2hairpin and -f $known_mature2hairpin and -f $novel_expression and -f $known_expression and -f $novel_mrd and -f $known_mrd and -f $novel_arf and -f $known_arf);
	my $time = time();
	mkdir "$Expression/expression_analyses/expression_analyses_$time" unless (-d "$Expression/expression_analyses/expression_analyses_$time");
	`grep -v "^#" $novel_csv | cat $known_csv - > $Expression/expression_analyses/expression_analyses_$time/miRNA_expressed.csv`;
	`cat $known_mature_converted $novel_mature_converted > $Expression/expression_analyses/expression_analyses_$time/mature.converted`;
	`cat $known_precursor_converted $novel_precursor_converted > $Expression/expression_analyses/expression_analyses_$time/precursor.converted`;
	`cat $known_mature2hairpin $novel_mature2hairpin > $Expression/expression_analyses/expression_analyses_$time/mature2hairpin`;
	`grep -v "^#" $novel_expression | cat $known_expression - > $Expression/miRNAs_expressed_all_samples_$time.csv`;
	`cat $known_mrd $novel_mrd > $Expression/expression_analyses/expression_analyses_$time/miRBase.mrd`;
	`cat $known_arf $novel_arf > $Expression/expression_analyses/expression_analyses_$time/All_miRNA_mapped.arf`;
	`grep -v "^#" $idir/miRDeep2/Novel_miRNA_Quantify/All_Novel_miRNA.count.list | cat $idir/Known_miRNA_identify/Known_miRNA_result/All_Known_miRNA.count.list - > $Expression/All_miRNA.count.list`;
	`cat $idir/Known_miRNA_identify/Known_miRNA_result/All_Known_miRNA.expressed.fa  $idir/miRDeep2/Novel_miRNA_Quantify/All_Novel_miRNA.expressed.fa > $Expression/All_miRNA.expressed.fa`;
	`cat $idir/Known_miRNA_identify/Known_miRNA_result/All_Known_miRNA_Pre.expressed.fa  $idir/miRDeep2/Novel_miRNA_Quantify/All_Novel_miRNA_Pre.expressed.fa > $Expression/All_miRNA_Pre.expressed.fa`;

	stepTime(6);
	$step++ unless ($sbs);
}


######################### extract the quantified miRNA
if($step==7){
	stepStart(7,"Extract the quantified miRNA");
	open(Nmature,">$idir/miRDeep2/miRNA_Quantify/Novel_mature_expressed.fa") or die $!;
	open(Kmature,">$idir/miRDeep2/miRNA_Quantify/Known_mature_expressed.fa") or die $!;
	open(IN,"$idir/miRDeep2/miRNA_Quantify/All_miRNA.expressed.fa") or die $!;  #########come frome the quantified csv file
	$/=">";
	<IN>;
	while(<IN>){
		chomp;
		my ($name,$mature_seq)=(split /\n/,$_,2)[0,1];
		$name=~s/\s+$//;
		$mature_seq=~s/\s+//g;
		if($name=~/novel/){
			print Nmature ">$name\n$mature_seq\n";
			#print Npre ">$name\n$preMI_seq{$name}\n";
		}
		else{
			print Kmature ">$name\n$mature_seq\n";
			#print Kpre ">$name\n$preMI_seq{$name}\n";
		}
	}
	$/="\n";
	close IN;
	close Nmature;
	close Kmature;

	open(Kpre,">$idir/miRDeep2/miRNA_Quantify/Known_Pre_expressed.fa") or die $!;
	open(Npre,">$idir/miRDeep2/miRNA_Quantify/Novel_Pre_expressed.fa") or die $!;
	open(NPS,"$idir/miRDeep2/miRNA_Quantify/All_miRNA_Pre.expressed.fa") or die $!;
	$/=">";
	<NPS>;
	while(<NPS>){
		chomp;
		my ($name,$seq)=(split /\n/,$_,2)[0,1];
		$name=~s/\s+$//;
		$seq=~s/\s+//g;
		if($name=~/novel/){
			print Npre ">$name\n$seq\n";
		}else{
			print Kpre ">$name\n$seq\n";
		}
	}
	$/="\n";
	close NPS;
	close Npre;
	close Kpre;
	stepTime(7);
	$step++ unless ($sbs);
}

######################### miRNA base distribution
if($step==8){
	stepStart(8,"base distribution");
	runOrDie("perl $Bin/bin/plot_first_base_dist.pl -i $idir/miRDeep2/miRNA_Quantify/Novel_mature_expressed.fa -o $idir/miRDeep2/miRNA_Quantify -p Novel");
	runOrDie("perl $Bin/bin/plot_base_dist_by_pos.pl -i $idir/miRDeep2/miRNA_Quantify/Novel_mature_expressed.fa -o $idir/miRDeep2/miRNA_Quantify -p Novel");
		runOrDie("perl $Bin/bin/length_stat.pl -i $idir/miRDeep2/miRNA_Quantify/Novel_mature_expressed.fa -od $idir/miRDeep2/miRNA_Quantify -prefix Novel_miRNA");
	my $num_known=`less "$idir/miRDeep2/miRNA_Quantify/Known_mature_expressed.fa"|grep '>' |wc -l`; chomp $num_known;
	if($num_known>1){
		runOrDie("perl $Bin/bin/plot_first_base_dist.pl -i $idir/miRDeep2/miRNA_Quantify/Known_mature_expressed.fa -o $idir/miRDeep2/miRNA_Quantify -p Known");
		runOrDie("perl $Bin/bin/plot_first_base_dist.pl -i $idir/miRDeep2/miRNA_Quantify/All_miRNA.expressed.fa -o $idir/miRDeep2/miRNA_Quantify -p All");
		runOrDie("perl $Bin/bin/plot_base_dist_by_pos.pl -i $idir/miRDeep2/miRNA_Quantify/Known_mature_expressed.fa -o $idir/miRDeep2/miRNA_Quantify -p Known");
		runOrDie("perl $Bin/bin/plot_base_dist_by_pos.pl -i $idir/miRDeep2/miRNA_Quantify/All_miRNA.expressed.fa -o $idir/miRDeep2/miRNA_Quantify -p All");
		runOrDie("perl $Bin/bin/length_stat.pl -i $idir/miRDeep2/miRNA_Quantify/Known_mature_expressed.fa -od $idir/miRDeep2/miRNA_Quantify -prefix Known_miRNA");
		runOrDie("perl $Bin/bin/length_stat.pl -i $idir/miRDeep2/miRNA_Quantify/All_miRNA.expressed.fa -od $idir/miRDeep2/miRNA_Quantify -prefix All_miRNA");
	}
	stepTime(8);
	$step++ unless($sbs);
}


###################### Normalize, TPM
if($step==9){
	stepStart(9,"Normalize the miRNA expression");
	runOrDie("perl $Bin/bin/count_to_expression.pl -i $idir/miRDeep2/miRNA_Quantify/All_miRNA.count.list -o $idir/miRDeep2/miRNA_Quantify/All_miRNA_expression.list");   ##### counts to expression
	stepTime(9);
	$step++ unless($sbs);
}


########################  summary
if($step==10){
	stepStart(10,"miRNA summary");
	runOrDie("perl $Bin/bin/summary2.pl -re_dir $idir -od $idir/miRDeep2/miRNA_Quantify -ref yes");
	my $seq=join ",",keys %Raw_FQ;
	my $fyes_p=1;
	runOrDie("perl $Bin/bin/summary_len.pl -dir $idir -fyes $fyes_p -min 15 -max 35 -seq $seq");
	stepTime(10);
	$step++ unless($sbs);
}


###################### creat pdf
if($step==11){
	stepStart(11,"creat pdf");
	runOrDie("perl $Bin/mirdeep_animal/pdf_create.pl -indir $idir/miRDeep2");
	stepTime(11);
}

totalTime();


####################
sub cmd_call {
	print "@_\n";
	system(@_) == 0 or die "system @_ failed: $?";
}


#########################sub
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


sub USAGE {#
	my $usage=<<"USAGE";
Program: $Script
Version: $version
Usage:
  Options:
	-uq		unquantified_reads.fa, fasta file. The ID should be this format: >SampleName_\\d+_x\\d+
	-idir		input directory
	-main_cfg	main_cfg
	-min		the min length of miRNA, default 18
	-max		the max length of miRNA, default 30
	-genome		reference genome fa file
	-plant		animal 0; plant 1
	-max_map	default 5
	-step		<num> default 0
		0		Divide unquantified reads
		1		miRDeep2 mapper.pl
		2		miRDeep2 miRDeep2.pl
		3		miRNA Filter
		4		extract miRNA
		5		miRNA Expression Quantify
		6		All miRNA Expression Combination
		7		extract the quantified miRNA
		8		miRNA base distribution
		9		Normalize the miRNA expression
		10		summary
		11		creat pdf
	-sbs		only one step
	-h		Help

USAGE
	print $usage;
	exit;
}

=c
	-yes	<str>	if the species have miRBase,this should be Latin abbreviations,for example vvi
	-no	<str>	if the species do not have miRBase,this should be All,Viridiplantae,Metazoa,Viruses,Chromalveolata,Mycetozoa
=cut

