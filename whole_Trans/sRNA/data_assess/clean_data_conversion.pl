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
my($od,$main_cfg,$word_cfg,$Q,$min,$max,$step,$sbs,$stat_file);
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
	"stat:s"		=>\$stat_file,
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
my %config=%{readconf("$Bin/../../CFG")};


########## read the config file
my %Raw_FQ;
my $adapter_limit=0;
my $adapter;
my %Read;
my %sample_info;
my %fqName_2_sampleName;

open(IN,$main_cfg) or die $!;
while(<IN>){
	chomp;
	next if(/^#/||/^$/);
	if(/^NAME/){
		my $name=(split /\s+/,$_)[1];
		my $fq=<IN>;
		my $fq_file = (split /\s+/,$fq)[1];
		$Raw_FQ{$name}=$fq_file;
		unless (-e $Raw_FQ{$name}){
			print "Check your fastq file: $name\n";
			die;
		}
		my ($fq_index) = $fq_file =~ /([^\/]+)_good_\d\.fq/;
		$fqName_2_sampleName{$fq_index} = $name;
		$Read{$name}{Total}=0;
	}
	if(/^adapter/){
		$adapter_limit=1;
		$adapter=(split /\s+/,$_)[1];
	}
}
close IN;

open (IN,$word_cfg) or die $!;
while(<IN>){
	chomp;
	next if (/^#/||/^$/);
	if(/^SampleID/){
		my ($sample,$bmk)=(split/\s+/,$_)[1,2];
		$sample_info{$bmk}=$sample;
	}
}
close IN;


################# fq2fa
if($step==0){
	#&log_current_time("FQ2FA begin:");
	stepStart(0,"FQ2FA");
	open(SH,">$work_sh/step0_fq2fa.sh") or die $!;
	foreach my $name (keys %Raw_FQ){
		print SH "perl $Bin/bin/fastq2fasta.pl $Raw_FQ{$name} > $od/${name}.clean.fa \n";
	}
	close SH;
	qsubOrDie("$od/work_sh/step0_fq2fa.sh",$config{queue},$config{cpu},$config{vf});
	#&log_current_time("FQ2FA done!");
	stepTime(0);
	$step++ unless ($sbs);
}

########### Clean collapse
if($step==1){
	#&log_current_time("Clean collapse begin:");
	stepStart(1,"Clean collapse");
	open (SH,">$work_sh/step1_clean_collapse.sh") or die $!;
	my @samples = (sort keys %Raw_FQ);
	foreach my $name(@samples){
		print SH "perl $Bin/bin/collapse_reads_md.pl $od/$name.clean.fa $name >$od/$name.collapse.fa\n";
	}
	close SH;
	qsubOrDie("$od/work_sh/step1_clean_collapse.sh",$config{queue},$config{cpu},$config{vf});
	#&log_current_time("Clean collapse done!");
	stepTime(1);
	$step++ unless($sbs);
}


############ All sample stat
if($step==2){
	#&log_current_time("All sample stat begin:");
	stepStart(2,"All sample stat");
	if ($stat_file and -f $stat_file) {
		my $SRNA_Rec_dir = dirname($stat_file);
		open IN,$stat_file;
		open STAT,">$od/All_sample_filter.stat";
		my $title = <IN>;
		my ($len_l,$len_g) = $title =~ /Length\<(\d+).*Length\>(\d+)/;
		print STAT "Samples\tBMK-ID\tRaw_reads\tLow_quality\tContaining\'N\'reads\tLength<$len_l\tLength>$len_g\tClean_reads\tQ30(\%)\n";
		while (<IN>) {
			chomp;
			next if (/^#/ or /^\s*$/);
			my @line = split "\t";
			next unless (defined $fqName_2_sampleName{$line[0]});
			print STAT "$sample_info{$fqName_2_sampleName{$line[0]}}\t$fqName_2_sampleName{$line[0]}\t$line[1]\t$line[6]\t$line[5]\t$line[3]\t$line[4]\t$line[7]\t$line[8]\n";
		}
		close IN;
		close STAT;
		if (-d "$SRNA_Rec_dir/Data_Assess/") {
			mkdir "$od/QC" unless (-d "$od/QC");
			my $sh_line = "";
			my $QC_stat_file = glob "$SRNA_Rec_dir/Data_Assess/AllSample_GC_Q.stat";
			open IN,$QC_stat_file;
			open OUT,">$od/QC/AllSample_GC_Q.stat";
			print OUT "SampleID\tReadSum\tBaseSum\tGC(\%)\tN(\%)\tQ20(\%)\tCycleQ20(\%)\tQ30(\%)\n";
			while (<IN>) {
				chomp;
				next if (/^#/ or /^\s*$/);
				my @line = split "\t";
				next unless (defined $fqName_2_sampleName{$line[0]});

				my @files = glob "$SRNA_Rec_dir/Data_Assess/$line[0]*";
				foreach my $file (@files) {
					my $cp_file = basename($file);
					$cp_file =~ s/$line[0]/$fqName_2_sampleName{$line[0]}/;
					runOrDie("cp $file $od/QC/$cp_file");
				}

				my @pngs = glob "$SRNA_Rec_dir/Data_Assess/PNG/$line[0]*";
				foreach my $png (@pngs) {
					my $cp_png = basename($png);
					$cp_png =~ s/$line[0]/$fqName_2_sampleName{$line[0]}/;
					runOrDie("cp $png $PNG/$cp_png");
				}

				$line[0] = $fqName_2_sampleName{$line[0]};
				my $new_line = join "\t",@line;
				print OUT "$new_line\n";

				if (-f "$od/QC/$line[0].quality") {
					$sh_line .= "cd $od &&";
					$sh_line .= "\tperl $Bin/bin/Error_Ration.pl -quality $od/QC/$line[0].quality -prefix $line[0] &&";
					$sh_line .= "\tcp $od/$line[0].error.png $od/PNG/$line[0].error.png\n";
				}
			}
			close IN;
			close OUT;

			if ($sh_line) {
				open(SH,">$work_sh/step2_error_draw.sh") or die $!;
				print SH "$sh_line";
				close SH;
				qsubOrDie("$od/work_sh/step2_error_draw.sh",$config{queue},$config{cpu},$config{vf});
			}

		}
	}
	
	#&log_current_time("All sample stat done!");
	stepTime(2);
	$step++ unless($sbs);
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
	-stat		/*/SRNA_Rec/All_Sample_stat.xls  file path, to get clean data qc stat info. unforced
			Example: /share/nas1/liux/testing/clean_fq_qc/SRNA_Rec/All_Sample_stat.xls

	-step		<num> default 0
		0		fq2fa
		1		clean read collapsed
		2		All sample stat
	-sbs		only one step
	-h		Help

USAGE
	print $usage;
	exit;
}
