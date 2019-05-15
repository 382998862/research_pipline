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

my($idir,$main_cfg,$min,$max,$fyes,$fno,$step,$sbs);
GetOptions(
	"help|?"	=>\&USAGE,
	"idir:s"	=>\$idir,
	"main_cfg:s"=>\$main_cfg,
	"min:s"		=>\$min,
	"max:s"		=>\$max,
#	"yes:s"		=>\$fyes, #This script, with no -yes and -no options, only predict novel miRNAs.
#	"no:s"		=>\$fno, #This script, with no -yes and -no options, only predict novel miRNAs.
	"step:s"	=>\$step,
	"sbs"		=>\$sbs,
) or &USAGE;
&USAGE unless ($idir and $main_cfg);
#&USAGE unless (defined $fyes or defined $fno) ; #This script, with no -yes and -no options, only predict novel miRNAs.

################default values
$step ||=1;
$min ||=15;
$max ||=35;

############### read CFG
my %CFG=%{readconf("$Bin/../CFG")};
my $miRBase = $CFG{miRBase};

$main_cfg=abs_path($main_cfg);
$idir=abs_path($idir);
my $miRDeep2_path="$idir/miRDeep2";
mkdir $miRDeep2_path unless (-d $miRDeep2_path);

#################read main_cfg
my $genome;
my $plant=1;
my $max_map=5;
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
	if (/^GENOME/) {
		$genome=(split/\s+/,$_)[1];
		unless (-e $genome) {
			print "Check Your FASTQ File $genome!";die;
		}
	}
	if (/^SPECIES_TYPE\s+(\S+)/) {
		$plant=$1;
	}
	if (/^max_map\s+(\S+)/) {
		$max_map=$1;
	}
}
close IN;

my $genome_name=basename($genome);
$genome="$idir/Alignment/Ref_Database/$genome_name";


=c
################## Creat miRBase file
if($step==0){
	stepStart(0,"Creat miRBase file");
	if(defined $fyes){
		#runOrDie("perl $Bin/bin/select_fa.pl -i $fyes -fa $miRBase/mature.fa -o $idir/Alignment/Ref_Database/$fyes.mature.fa >>$idir/Alignment/Ref_Database/select_mature.out");
		#runOrDie("perl $Bin/bin/select_fa.pl -i $fyes -fa $miRBase/hairpin.fa -o $idir/Alignment/Ref_Database/hairpin_$fyes.fa >>$idir/Alignment/Ref_Database/select_hairpin.out");
		&cmd_call("perl $Bin/bin/get_mature_pre.pl -abbre $fyes -ilist $Bin/bin/All_miRNA.list -omature $idir/Alignment/Ref_Database/$fyes.mature.fa -opre $idir/Alignment/Ref_Database/hairpin_$fyes.fa -olist $idir/Alignment/Ref_Database/$fyes.miRNA.list");
	}
	if(defined $fno){
		runOrDie("perl $Bin/bin/select_fa.pl -i $fno -fa $miRBase/mature.fa -o $idir/Alignment/Ref_Database/$fno.mature.fa >>$idir/Alignment/Ref_Database/select_mature.out");
	}
	stepTime(0);
	$step++ unless ($sbs);
}
=cut


################# miRDeep2 mapper.pl
if($step==1){
	stepStart(1,"miRDeep2 mapper.pl");
	unless((-f "$idir/miRDeep2/Total_reads_collapsed.fa") || (-f "$idir/miRDeep2/Total_reads_collapsed_vs_genome.arf")){
		my %FQ;
		open (OUT,">$idir/miRDeep2/mapper.cfg") or die $!;
		foreach my $name (keys %Raw_FQ){
			$FQ{$name}="$idir/Alignment/$name/genome_map/$name.genome_map.fa";
			#$FQ{$name}="$idir/Alignment/$name/$name.no_rRNA_rep.fa";
			unless(-e "$idir/miRDeep2/$name.mapper.fa"){
				`ln -s $FQ{$name} $idir/miRDeep2/$name.mapper.fa`;
			}
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
		if(defined $fyes){
			runOrDie("cd $idir/miRDeep2 && perl $Bin/mirdeep_animal/miRDeep2.pl Total_reads_collapsed.fa $genome Total_reads_collapsed_vs_genome.arf $idir/Alignment/Ref_Database/$fyes.mature.fa none $idir/Alignment/Ref_Database/hairpin_$fyes.fa -d -g -1 2>report.log");
		}
		elsif(defined $fno){
			runOrDie("cd $idir/miRDeep2 && perl $Bin/mirdeep_animal/miRDeep2.pl Total_reads_collapsed.fa $genome Total_reads_collapsed_vs_genome.arf none $idir/Alignment/Ref_Database/$fno.mature.fa none -d -g -1 2> report.log");
		}
		else {
            runOrDie("cd $idir/miRDeep2 && perl $Bin/mirdeep_animal/miRDeep2.pl Total_reads_collapsed.fa $genome Total_reads_collapsed_vs_genome.arf none none none -d -g -1 2> report.log");
        }
	}
	elsif($plant==1){
		if(defined $fyes){
			runOrDie("cd $idir/miRDeep2 && perl $Bin/mirdeep_plant/miRDeep2.pl Total_reads_collapsed.fa $genome Total_reads_collapsed_vs_genome.arf $idir/Alignment/Ref_Database/$fyes.mature.fa none $idir/Alignment/Ref_Database/hairpin_$fyes.fa -d -g 50000 -l 250  -m 10 -v -P -n d  2> report.log");
		}
		elsif(defined $fno){
			runOrDie("cd $idir/miRDeep2 && perl $Bin/mirdeep_plant/miRDeep2.pl Total_reads_collapsed.fa $genome Total_reads_collapsed_vs_genome.arf none $idir/Alignment/Ref_Database/$fno.mature.fa none -g 50000 -l 250 -d -m 10 -v -P -n d 2> report.log")
		}
		else {
            runOrDie("cd $idir/miRDeep2 && perl $Bin/mirdeep_plant/miRDeep2.pl Total_reads_collapsed.fa $genome Total_reads_collapsed_vs_genome.arf none none none -g 50000 -l 250 -d -m 10 -v -P -n d 2> report.log");
        }
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

	if(defined $fyes){
		runOrDie("perl $Bin/bin/$parse_version -str $idir/miRDeep2/mirdeep_runs/$work_dir/tmp/precursors.str -coord $idir/miRDeep2/mirdeep_runs/$work_dir/tmp/precursors.coords -mrd $idir/miRDeep2/mirdeep_runs/$work_dir/output.mrd -samples $seq -od $idir/miRDeep2/Result -known $idir/Alignment/Ref_Database/$fyes.mature.fa");
	}
	else{
		runOrDie("perl $Bin/bin/$parse_version -str $idir/miRDeep2/mirdeep_runs/$work_dir/tmp/precursors.str -coord $idir/miRDeep2/mirdeep_runs/$work_dir/tmp/precursors.coords -mrd $idir/miRDeep2/mirdeep_runs/$work_dir/output.mrd -samples $seq -od $idir/miRDeep2/Result");
	}

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
	}
	close IN;
	if (defined $fyes){
		open(Kmature,">$idir/miRDeep2/Result/Known_mature.fa") or die $!;
		open(Kpre,">$idir/miRDeep2/Result/Known_pre.fa") or die $!;
		foreach my $key (keys %{$MIR{Known}}){
			print Kmature ">$key\n$MIR{Known}{$key}{mature}\n";
			print Kpre ">$key\n$MIR{Known}{$key}{precursor}\n";
		}
		close Kmature;
		close Kpre;
	}
	open (Nmature,">$idir/miRDeep2/Result/Novel_mature.fa") or die $!;
	open (Npre,">$idir/miRDeep2/Result/Novel_Pre.fa") or die $!;
	open (INFO,">$idir/miRDeep2/Result/chname.list") or die $!;
	my $count = 0;
	foreach my $key (keys %{$MIR{novel}}){
		$count ++;
		print INFO "$key\tnovel_miR_$count\n";
		print Nmature ">Novel_miR_$key\n$MIR{novel}{$key}{mature}\n";
		print Npre ">Novel_miR_$key\n$MIR{novel}{$key}{precursor}\n";
	}
	close Nmature;
	close Npre;
	close INFO;
	if(defined $fyes){
		runOrDie("cat $idir/Alignment/Ref_Database/$fyes.mature.fa $idir/miRDeep2/Result/Novel_mature.fa >$idir/miRDeep2/Result/All_miRNA.fa");
		runOrDie("cat $idir/Alignment/Ref_Database/hairpin_$fyes.fa $idir/miRDeep2/Result/Novel_Pre.fa >$idir/miRDeep2/Result/All_miRNA_Pre.fa");
	}else{
		runOrDie("cat $idir/miRDeep2/Result/Novel_mature.fa >$idir/miRDeep2/Result/All_miRNA.fa");
		runOrDie("cat $idir/miRDeep2/Result/Novel_Pre.fa >$idir/miRDeep2/Result/All_miRNA_Pre.fa");
	}
	stepTime(4);
	$step++ unless ($sbs);
}

######################## miRNA Expression Quantify
if($step==5){
	stepStart(5,"miRNA Expression Quantify");
	my $Expression="$idir/miRDeep2/miRNA_Quantify";
	mkdir $Expression unless (-d $Expression);
	if(-d "$idir/miRDeep2/miRNA_Quantify/expression_analyses"){
		&cmd_call("rm -r $idir/miRDeep2/miRNA_Quantify/expression_analyses");
		&cmd_call("rm $idir/miRDeep2/miRNA_Quantify/miRNAs_expressed_all_samples_*.csv");
	}
	runOrDie("cd $Expression && perl $Bin/mirdeep_animal/quantifier.pl -p $idir/miRDeep2/Result/All_miRNA_Pre.fa -m $idir/miRDeep2/Result/All_miRNA.fa -r $idir/miRDeep2/Total_reads_collapsed.fa -d -g 1 -e 2 -f 5 -W ");
	my $expression = (glob "$idir/miRDeep2/miRNA_Quantify/miRNAs_expressed_all_samples_*.csv")[0];
	runOrDie("perl $Bin/bin/expression_extra_v0.pl -exp $expression -fa $idir/miRDeep2/Result/All_miRNA.fa -pre $idir/miRDeep2/Result/All_miRNA_Pre.fa -oexp $idir/miRDeep2/miRNA_Quantify/All_miRNA.count.list -ofa $idir/miRDeep2/miRNA_Quantify/All_miRNA.expressed.fa -opre $idir/miRDeep2/miRNA_Quantify/All_miRNA_Pre.expressed.fa");
	stepTime(5);
	$step++ unless ($sbs);
}

######################### extract the quantified miRNA
if($step==6){
	stepStart(6,"Extract the quantified miRNA");
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
	stepTime(6);
	$step++ unless ($sbs);
}

######################### miRNA base distribution
if($step==7){
	stepStart(7,"base distribution");
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
	stepTime(7);
	$step++ unless($sbs);
}


###################### Normalize, TPM
if($step==8){
	stepStart(8,"Normalize the miRNA expression");
	runOrDie("perl $Bin/bin/count_to_expression.pl -i $idir/miRDeep2/miRNA_Quantify/All_miRNA.count.list -o $idir/miRDeep2/miRNA_Quantify/All_miRNA_expression.list");   ##### counts to expression
	stepTime(8);
	$step++ unless($sbs);
}


########################  summary
if($step==9){
	stepStart(9,"miRNA summary");
	if(defined $fyes){
		runOrDie("perl $Bin/bin/summary2.pl -re_dir $idir -od $idir/miRDeep2/miRNA_Quantify -ref yes");
	}else{
		runOrDie("perl $Bin/bin/summary2.pl -re_dir $idir -od $idir/miRDeep2/miRNA_Quantify -ref no");
	}
	my $seq=join ",",keys %Raw_FQ;
	my $fyes_p=0;
	if (defined $fyes){$fyes_p=1};
	runOrDie("perl $Bin/bin/summary_len.pl -dir $idir -fyes $fyes_p -min 15 -max 35 -seq $seq");
	stepTime(9);
	$step++ unless($sbs);
}


###################### creat pdf
if($step==10){
	stepStart(10,"creat pdf");
	runOrDie("perl $Bin/mirdeep_animal/pdf_create.pl -indir $idir/miRDeep2");
	stepTime(10);
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
	-idir		input directory
	-main_cfg	main_cfg
	-min		the min length of miRNA, default 18
	-max		the max length of miRNA, default 30
	-step		<num> default 0
		1		miRDeep2 mapper.pl
		2		miRDeep2 miRDeep2.pl
		3		miRNA Filter
		4		extract miRNA
		5		miRNA Expression Quantify
		6		extract the quantified miRNA
		7		miRNA base distribution
		8		Normalize the miRNA expression
		9		summary
		10		creat pdf
	-sbs		only one step
	-h		Help

USAGE
	print $usage;
	exit;
}
