#!/usr/bin/perl -w
use strict;
use Getopt::Long;
use Data::Dumper;
use FindBin qw($Bin $Script);
use File::Basename qw(basename dirname);
use Cwd qw(abs_path);
use newPerlBase;
my $Title="data_assess.pl";
my $version="1.0";

#-------------------------------------------------
# GetOptions
#-------------------------------------------------
my($od,$main_cfg,$word_cfg,$Q,$min,$max,$step,$sbs);
GetOptions(
	"help|?"		=>\&USAGE,
	"od:s"			=>\$od,
	"main_cfg:s"	=>\$main_cfg,
	"word_cfg:s"	=>\$word_cfg,
	"Q:s"			=>\$Q,
	"min:s"			=>\$min,
	"max:s"			=>\$max,
	"sbs:s"			=>\$sbs,
	"step:s"		=>\$step,
) or &USAGE;
&USAGE unless ($od and $main_cfg and $word_cfg);

$step ||=0;
$Q   ||=33;
$min ||=15;
$max ||=35;

$main_cfg=abs_path($main_cfg);
$word_cfg=abs_path($word_cfg);
$od=abs_path($od);
mkdir $od unless (-d $od);
my $work_sh="$od/work_sh";
mkdir $work_sh unless(-d $work_sh);
my $PNG="$od/PNG";
mkdir $PNG unless (-d $PNG);

#createLog($Title,$version,$$,"$od/log/","test");

########################## bin
my %config=%{readconf("$Bin/../CFG")};
my $fastq_filter_by_Qxx = $config{fastq_filter_by_Qxx};
my $fastq_filter_by_Ns  = $config{fastq_filter_by_Ns};


########## read the config file
my %Raw_FQ;
my $adapter_limit=0;
my $adapter;
my %Read;
my %sample_info;


open(IN,$main_cfg) or die $!;
while(<IN>){
	chomp;
	next if(/^#/||/^$/);
	if(/^NAME/){
		my $name=(split /\s+/,$_)[1];
		my $fq=<IN>;
		$Raw_FQ{$name}=(split /\s+/,$fq)[1];
		unless (-e $Raw_FQ{$name}){
			print "Check your fastq file: $name\n";
			die;
		}
		$Read{$name}{Total}=0;
	}
	if(/^SampleID/){
		my ($sample,$bmk)=(split/\s+/,$_)[1,2];
		$sample_info{$bmk}=$sample;
	}
}
close IN;

open (IN,$word_cfg) or die $!;
while(<IN>){
	chomp;
	next if (/^#/||/^$/);
	if(/^adapter/){
		$adapter_limit=1;
		$adapter=(split /\s+/,$_)[1];
	}
}
close IN;

################ filter the low quality sequence
if($step==0){
	#&log_current_time("Low quality filter begin:");
	stepStart(0,"Low quality filter");
	open(SH,">$work_sh/step0_quality.sh") or die $!;
	foreach my $name(sort keys %Raw_FQ){
		print SH "$fastq_filter_by_Qxx -a $Raw_FQ{$name} -f $name -o $od -w 30 -q 0.85 -Q $Q\n";
	}
	close SH;
	qsubOrDie("$od/work_sh/step0_quality.sh",$config{queue},$config{cpu},$config{vf});
	stepTime(0);
	#log_current_time("Low quality filter done!");
	$step++ unless ($sbs);
}

##################filter N
if($step==1){
	#&log_current_time("N filter begin:");
	stepStart(1,"N filter");
	open(SH,">$work_sh/step1_filter_N.sh") or die $!;
	foreach my $name(sort keys %Raw_FQ){
		print SH "$fastq_filter_by_Ns -a $od/${name}_good_1.fq -u 0.1 -c $od/${name}_N.fq\n";
	}
	close SH;
	qsubOrDie("$od/work_sh/step1_filter_N.sh",$config{queue},$config{cpu},$config{vf});
	#&log_current_time("N filter dong!");
	stepTime(1);
	$step++ unless ($sbs);
}

################# fq2fa
if($step==2){
	#&log_current_time("FQ2FA begin:");
	stepStart(2,"FQ2FA");
	open(SH,">$work_sh/step2_fq2fa.sh") or die $!;
	foreach my $name (keys %Raw_FQ){
		print SH "perl $Bin/bin/fastq2fasta.pl $od/${name}_N.fq > $od/${name}_N.fa \n";
	}
	close SH;
	qsubOrDie("$od/work_sh/step2_fq2fa.sh",$config{queue},$config{cpu},$config{vf});
	#&log_current_time("FQ2FA done!");
	stepTime(2);
	$step++ unless ($sbs);
}

################ filter 3' adapter
if($step==3){
	#&log_current_time("Filter 3' adapter begin:");
	stepStart(3,"Filter 3' adapter");
	if ($adapter_limit==1){
		open(SH,">$work_sh/step3_adapter3.sh") or die $!;
		foreach my $name (keys %Raw_FQ){
			print SH "perl $Bin/bin/clip_adapters.pl $od/${name}_N.fa $adapter > $od/${name}_reads_clip.fa\n";
		}
		close SH;
		qsubOrDie("$od/work_sh/step3_adapter3.sh",$config{queue},$config{cpu},$config{vf});
	}
	else{
		foreach my $name (keys %Raw_FQ){
			open(SH,">$work_sh/step3_adapter3.sh") or die $!;
			print SH "cp -r $od/${name}_N.fa $od/${name}_reads_clip.fa\n";
			runOrDie("$od/work_sh/step3_adapter3.sh");
		}
	}
	#&log_current_time("Filter 3' adapter done!");
	stepTime(3);
	$step++ unless($sbs);
}


################# length filter
if($step==4){
	#&log_current_time("Length Filter begin:");
	stepStart(4,"Length Filter");
	foreach my $name(keys %Raw_FQ){
		$/=">";
		open(IN,"$od/${name}_reads_clip.fa") or die $!;
		open(OUT,">$od/$name.clean.fa") or die $!;
		<IN>;
		while(<IN>){
			chomp;
			my $len=(split /\n/,$_)[1];
			$len=length $len;
			if($len>=$min && $len<=$max){
				print OUT ">$_";
			}
		}
		close IN;
		close OUT;
		$/="\n";
	}
	#&log_current_time("Length Filter done!");
	stepTime(4);
	$step++ unless($sbs);
}


########### Clean collapse
if($step==5){
	#&log_current_time("Clean collapse begin:");
	stepStart(5,"Clean collapse");
	open (SH,">$work_sh/step4_clean_collapse.sh") or die $!;
	my @samples = (sort keys %Raw_FQ);
	foreach my $name(@samples){
		print SH "perl $Bin/bin/collapse_reads_md.pl $od/$name.clean.fa $name >$od/$name.collapse.fa\n";
	}
	close SH;
	qsubOrDie("$od/work_sh/step4_clean_collapse.sh",$config{queue},$config{cpu},$config{vf});
	#&log_current_time("Clean collapse done!");
	stepTime(5);
	$step++ unless($sbs);
}


############Data assess
if($step==6){
	#&log_current_time("Data assess begin:");
	stepStart(6,"Data assess");
	mkdir "$od/QC" unless (-d "$od/QC");
	open(OUT,">$work_sh/assess.cfg") or die $!;
	foreach my $name(sort keys %Raw_FQ){
		print OUT "Sample\t$name\n";
		print OUT "fq1\t$od/${name}_N.fq\n";
	}
	close OUT;
	open(SH,">$od/work_sh/step5_data_assess.sh") or die $!;
	print SH "perl $Bin/rna_seq_data_assess.pl -Q $Q -config $work_sh/assess.cfg -outdir $od/QC";
	close SH;
	runOrDie("$od/work_sh/step5_data_assess.sh");
	runOrDie("cp $od/QC/PNG/* $od/PNG");
	#&log_current_time("Data assess done!");
	stepTime(6);
	$step++ unless($sbs);
}

############ Error assess
if($step==7){
	#&log_current_time("Error assess begin:");
	stepStart(7,"Error assess");
	mkdir "$od/PNG" unless (-d "$od/PNG");
	open(SH,">$work_sh/step6_error_draw.sh") or die $!;
	foreach my $name(keys %Raw_FQ){
		print SH "cd $od &&";
		print SH "\tperl $Bin/bin/Error_Ration.pl -quality $od/QC/$name.quality -prefix $name &&";
		print SH "\tcp $od/$name.error.png $od/PNG/$name.error.png\n";
	}
	close SH;
	qsubOrDie("$od/work_sh/step6_error_draw.sh",$config{queue},$config{cpu},$config{vf});
	#&log_current_time("Error assess done!");
	stepTime(7);
	$step++ unless($sbs);
}

############# Q stat
if($step==8){
	#&log_current_time("Q stat begin:");
	stepStart(8,"Q stat");
	foreach my $name(keys %Raw_FQ){
		$Read{$name}{L}=0;
		$Read{$name}{M}=0;
		$Read{$name}{C}=0;
		$/=">";
		open(IN,"$od/${name}_reads_clip.fa") or die $!;
		<IN>;
		while(<IN>){
			chomp;
			my $len=(split /\n/,$_)[1];
			$len=length $len;
			if($len<$min){
				$Read{$name}{L}++;
			}
			elsif($len>$max){
				$Read{$name}{M}++;
			}
			else{
				$Read{$name}{C}++;
			}
		}
		close IN;
		$/="\n";
	}

	my %Q30;
	open(GC,"$od/QC/AllSample_GC_Q.stat") or die $!;
	<GC>;
	while(<GC>){
		chomp;
		my ($name,$q)=(split/\s+/,$_)[0,-1];
		$Q30{$name}=$q;
	}
	close GC;
	open(FILTER,">$od/All_sample_filter.stat") or die $!;
	print FILTER "#BMK-ID\tRaw_reads\tLow_quality\tContaining'N'reads\tLength<$min\tLength>$max\tClean_reads\tQ30(%)\n";
	foreach my $name (sort keys %Raw_FQ) {
		$Read{$name}{Q}=0;
		$Read{$name}{N}=0;
		$Read{$name}{Q}=(split/\s+/,`wc -l $od/${name}_good_1.fq`)[0]/4;
		$Read{$name}{N}=(split/\s+/,`wc -l $od/${name}_N.fq`)[0]/4;
		$Read{$name}{Total}=(split/\s+/,`wc -l $Raw_FQ{$name}`)[0]/4;

		open (OUT,">$od/$name.filter.stat") or die $!;
		print OUT "#Type\tNumber\tPercentage\n";
		print OUT "Total_reads_number\t$Read{$name}{Total}\t100.00%\n";

		my ($q_n,$q_p,$n_n,$n_p,$l_p,$m_p,$c_p)=(0,0,0,0,0,0,0);

		$q_n=$Read{$name}{Total}-$Read{$name}{Q};
		$q_p=sprintf "%.2f",$q_n/$Read{$name}{Total}*100;
		print OUT "Filter_low_quality_reads\t$q_n\t$q_p%\n";

		$n_n=$Read{$name}{Q}-$Read{$name}{N};
		$n_p=sprintf "%.2f",$n_n/$Read{$name}{Total}*100;
		print OUT "Filter_having_'N'_reads\t$n_n\t$n_p%\n";

		$l_p=sprintf "%.2f",$Read{$name}{L}/$Read{$name}{Total}*100;
		print OUT "Length<$min\t$Read{$name}{L}\t$l_p%\n";

		$m_p=sprintf "%.2f",$Read{$name}{M}/$Read{$name}{Total}*100;
		print OUT "Length>$max\t$Read{$name}{M}\t$m_p%\n";

		$c_p=sprintf "%.2f",$Read{$name}{C}/$Read{$name}{Total}*100;
		print OUT "Clean_Reads\t$Read{$name}{C}\t$c_p%\n";

		close OUT;

		print FILTER "$name\t$Read{$name}{Total}\t$q_n\t$n_n\t$Read{$name}{L}\t$Read{$name}{M}\t$Read{$name}{C}\t$Q30{$name}\n";
	}
	close FILTER;
	#&log_current_time("Q stat done!");
	stepTime(8);
}

totalTime();

################################ sub
sub log_current_time{
        # get parameter
     my ($info) = @_;

     # get current time with string
     my $curr_time = date_time_format(localtime(time()));

     # print info with time
     print "[$curr_time] $info\n";
}


sub date_time_format{
        my ($sec, $min, $hour, $day, $mon, $year, $wday, $yday, $isdst)=localtime(time());
    return sprintf("%4d-%02d-%02d %02d:%02d:%02d", $year+1900, $mon+1, $day, $hour, $min, $sec);
}



sub USAGE {#
	my $usage=<<"USAGE";
Program: $Script
Version: $version
Usage:
  Options:
	-od		output directory
	-main_cfg	detail.cfg
	-word_cfg	word.cfg
	-Q		forced,default 33
	-min		the min length of miRNA ,forced, default 18
	-max		the max length of miRNA ,forced, default 30
	-step		<num> default 0
		0		filter the low quality sequence
		1		filter N
		2		fq2fa
		3		filter 3' adapter
		4		length filter
		5		clean read collapsed
		6		Data assess
		7		Error assess
		8		Q stat
	-sbs		only one step
	-h		Help

USAGE
	print $usage;
	exit;
}
