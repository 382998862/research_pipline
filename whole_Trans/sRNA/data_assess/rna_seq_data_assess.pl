#!/usr/local/bin/perl -w
#
# Copyright (c) BMK 2009
# Writer:         Yangsh <yangsh@biomarker.com.cn>
# Program Date:   2010.
# Modifier:       lium <yangsh@biomarker.com.cn>
# Last Modified:  2013/9/16.
my $version="1.5";
my $BEGIN=time();

use strict;
use Cwd;
use Getopt::Long;
use Data::Dumper;
use FindBin qw($Bin $Script);
use File::Basename qw(basename dirname);
use newPerlBase;


# ------------------------------------------------------------------
# GetOptions
# ------------------------------------------------------------------
my %opts;
GetOptions(\%opts,"Q=i","config=s","outdir=s","byraw=i","h" ) or &USAGE;
# check
&USAGE if(!defined($opts{Q}) || !defined($opts{"outdir"}) || !defined($opts{"config"}) || defined($opts{h}));
if( !defined($opts{"byraw"}) ) { $opts{"byraw"} = 0; }



###############Time
my $Time_Start;
$Time_Start = sub_format_datetime(localtime(time()));
print "\nStart Time :[$Time_Start]\n\n";



################ get argv and outdir
my $pwd=`pwd`;chomp($pwd);
my $programe_dir=basename($0);
my $path=dirname($0);
my ($Q,$cfg,$od) = ($opts{Q}, $opts{"config"}, $opts{"outdir"});
`mkdir $od/work_sh -p` if(!-d "$od/work_sh");
$od=ABSOLUTE_DIR($od,$pwd);

#createLog($Title,$version,$$,"$od/log/","test");
################ read cfg
my %CFG=%{readconf("$Bin/../CFG")};


################ program path ################
our $GC_Q_SVG_pl 	= "$Bin/bin/GCqual_Solexa_check_svg.pl";
our $cycle_Q_svg_pl	= "$Bin/bin/Cycle_Q_SVG.pl";
our $svg2png     	= $CFG{svg2png};
our $fastq_qc_stat 	= $CFG{fastq_qc_stat};




################ pipeline ####################
# get sample information from config file
my %fqpair;
if( $opts{"byraw"} == 0 ) {get_sample_info_from_config_by_fq($cfg,\%fqpair);} 	# by fq*
else {get_sample_info_from_config_by_raw_fq($cfg,\%fqpair);} 			# by raw_fq*

# check sample information
check_sample_information(\%fqpair);

# create work shell
create_work_shell($od, \%fqpair, $Q);

# calling qsub
calling_qsub($od);

# summary_qc_stat_info
summary_qc_stat_info($od);

# convert svg into png
`cd $od && perl $svg2png ./ PNG && cd ../ `;

###############Time
my $Time_End;
$Time_End = sub_format_datetime(localtime(time()));
print "\nEnd Time :[$Time_End]\n\n";
&Runtime($BEGIN);





#+------------------------------------
# create work shell for data quality stat,
# and drawing quality graph
#+------------------------------------
sub create_work_shell{
	# get parameters
	my ($od, $SamData, $Q) = @_;

	# create shell script file
	open (OUT,">$od/work_sh/GC_Q_svg.sh");

	# iter every sample
	foreach my $sample_name (sort keys %$SamData) {
		# for pair-end sequencing
		if( exists $SamData->{"$sample_name"}->{"2"} ){
			my $fq1 = $SamData->{"$sample_name"}->{"1"};
			my $fq2 = $SamData->{"$sample_name"}->{"2"};
			print "fq1: $fq1\n";
			print "fq2: $fq2\n";
			print OUT "cd  $od && $fastq_qc_stat -Q $Q -a $fq1 -b $fq2 -f $sample_name -q 45&& ";
		# for single sequencing
		} else {
			my $fq1 = $SamData->{"$sample_name"}->{"1"};
			print OUT "cd  $od && $fastq_qc_stat -Q $Q -a $fq1 -f $sample_name  -q 45&& ";
		}
		# drawing cmd
		print OUT "perl $GC_Q_SVG_pl  $od/$sample_name.quality $od/$sample_name.acgtn $od/ &&";
		print OUT "perl $cycle_Q_svg_pl -i $od/$sample_name.quality -o $od/$sample_name &&\n";
	}

	# close file
	close(OUT);
}




#+------------------------------------
# calling qsub
#+------------------------------------
sub calling_qsub{
	# get parameters
	my ($od) = @_;

	# get current number of qsub
	my $qstatN = `qstat|wc`;
	$qstatN=~s/^\s+//;
	$qstatN=(split/\s+/,$qstatN)[0];

	# calculate new max process number
	my $qsubNum=$qstatN+10;

	# calling qsub
	#qsubOrDie("$od/work_sh/GC_Q_svg.sh",$CFG{queue},$qsubNum,$CFG{vf});
	qsubOrDie("$od/work_sh/GC_Q_svg.sh",$CFG{queue},$CFG{cpu},$CFG{vf});
	#system "sh /share/nas2/genome/bmksoft/tool/qsub_sge_plus/v1.0/qsub_sge.plus.sh $od/work_sh/GC_Q_svg.sh --reqsub --maxproc $qsubNum" ;
}


#+------------------------------------
# summary all quality stat information
#+------------------------------------
sub summary_qc_stat_info{
	# get parameters
	my ($od) = @_;

	# get all sample stat file
	my @STAT=glob("$od/*.stat");
	my @Paraf;
	foreach my $stat (@STAT) {
		push @Paraf,$stat unless $stat=~/.+cycleQ\.stat/;
	}

	# open summary stat file
	open (OUT,">$od/AllSample_GC_Q.stat");
	print OUT "SampleID\tReadSum\tBaseSum\tGC(%)\tN(%)\tQ20(%)\tCycleQ20(%)\tQ30(%)\n";
	foreach my $paraf (@Paraf) {
		# get file name and basename
		my $file=basename($paraf);
		print "file: $file\n";
		next if($file=~/AllSample_GC_Q.stat$/); 	# ignore self
		$file=~s/\.stat$//;

		# get last line information
		open(IN,"$paraf")||die"can't open $paraf\n";
		<IN>;
		my $now;
		while (<IN>)
		{
			chomp;
			my @A=split/\s+/,$_,2;
			$now="$file\t$A[1]\n";
		}
		close(IN);
		print OUT $now;
	}
	close(OUT);
}





#+------------------------------------
# check sample information
# if a sample doesnot have raw_fq1, then error
# NOTE: raw_fq2 not exists for single read
#+------------------------------------
sub check_sample_information {
	# get parameters
	my ($SamData) = @_;

	# foreach
	foreach my $sample_name (keys %$SamData) {
		# check raw_fq1
		if( exists $SamData->{$sample_name}->{"1"} ) {
			next;
		}
		# error
		print "sample error: ($sample_name) deesnot have raw_fq1 or fq1\n";
		exit;
	}
}



#+------------------------------------
# get sample information by raw_fq*
#+------------------------------------
sub get_sample_info_from_config_by_raw_fq {
	# get parameters
	my ($config_file, $SamData) = @_;

	# open file
	my $error_status = 0;
	open IN,$config_file || die "fail open: $config_file";
	# init variables
	my $line = "";
	my $sample_name = "";
	my $control = 0; 	# for record reading
	while( $line = <IN> ) {
		# process empty and '\n'
		chomp($line);
		$line =~ s/\s+$//;
		$line =~ s/\r$//;

		# ignore empty line and comment
		next if( ($line=~/^$/) or ($line=~/^\#/) );

		# reading sample name
		if ( $line =~ /^Sample/ ) {
			$sample_name = (split(/\s+/, $line))[1];
			if( $sample_name eq "" ) {
				print "config error: extract sample name failed!\n"; exit;
			}
		}
		# reading raw_fq1
		if( $line =~ /^raw_fq1/ ){
			# check sample name
			if( $sample_name eq "" ) {
				print "config error: sample name is empty!\n"; exit;
			}
			# check repeat sample name
			if( exists $SamData->{$sample_name}->{"1"} ){
				print "config error: repeat sample name $sample_name\n"; exit;
			}
			# store file path
			$SamData->{$sample_name}->{"1"} = (split(/\s+/, $line))[1];
		}
		# reading raw_fq2
		if( $line =~ /^raw_fq2/ ){
			# check sample name
			if( $sample_name eq "" ) {
				print "config error: sample name is empty!\n"; exit;
			}
			# check repeat sample name
			if( exists $SamData->{$sample_name}->{"2"} ){
				print "config error: repeat sample name $sample_name\n"; exit;
			}
			# store file path
			$SamData->{$sample_name}->{"2"} = (split(/\s+/, $line))[1];
		}
	}
}


#+------------------------------------
# get sample information by fq*
#+------------------------------------
sub get_sample_info_from_config_by_fq {
	# get parameters
	my ($config_file, $SamData) = @_;

	# open file
	my $error_status = 0;
	open IN,$config_file || die "fail open: $config_file";
	# init variables
	my $line = "";
	my $sample_name = "";
	my $control = 0; 	# for record reading
	while( $line = <IN> ) {
		# process empty and '\n'
		chomp($line);
		$line =~ s/\s+$//;
		$line =~ s/\r$//;

		# ignore empty line and comment
		next if( ($line=~/^$/) or ($line=~/^\#/) );

		# reading sample name
		if ( $line =~ /^Sample/ ) {
			$sample_name = (split(/\s+/, $line))[1];
			if( $sample_name eq "" ) {
				print "config error: extract sample name failed!\n"; exit;
			}
		}
		# reading fq1
		if( $line =~ /^fq1/ ){
			# check sample name
			if( $sample_name eq "" ) {
				print "config error: sample name is empty!\n"; exit;
			}
			# check repeat sample name
			if( exists $SamData->{$sample_name}->{"1"} ){
				print "config error: repeat sample name $sample_name\n"; exit;
			}
			# store file path
			$SamData->{$sample_name}->{"1"} = (split(/\s+/, $line))[1];
		}
		# reading fq2
		if( $line =~ /^fq2/ ){
			# check sample name
			if( $sample_name eq "" ) {
				print "config error: sample name is empty!\n"; exit;
			}
			# check repeat sample name
			if( exists $SamData->{$sample_name}->{"2"} ){
				print "config error: repeat sample name $sample_name\n"; exit;
			}
			# store file path
			$SamData->{$sample_name}->{"2"} = (split(/\s+/, $line))[1];
		}
	}
}





sub MKDIR { # &MKDIR($out_dir);
	my ($dir)=@_;
	rmdir($dir) if(-d $dir);
	mkdir($dir) if(!-d $dir);
}

sub Runtime
{ # &Runtime($BEGIN);
	my ($t1)=@_;
	my $t=time()-$t1;
	print "Total elapsed time: ${t}s\n";
}
sub Done
{
	my ($out)=@_;
	print "$out Done ~ ......\n";
}
sub sub_format_datetime {#Time calculation subroutine
    my($sec, $min, $hour, $day, $mon, $year, $wday, $yday, $isdst) = @_;
	$wday = $yday = $isdst = 0;
    sprintf("%4d-%02d-%02d %02d:%02d:%02d", $year+1900, $mon+1, $day, $hour, $min, $sec);
}
sub ABSOLUTE_DIR
{
        my ($in,$cur_dir)=@_;
        my $return="";

        if(-f $in)
        {
                my $dir=dirname($in);
                my $file=basename($in);
                chdir $dir;$dir=`pwd`;chomp $dir;
                $return="$dir/$file";
        }
        elsif(-d $in)
        {
                chdir $in;$return=`pwd`;chomp $return;
        }
        else
        {
                warn "Warning just for file and dir\n";
                exit;
        }
        chdir $cur_dir;
        return $return;
}



################################################################################################################
sub USAGE {#
	my $usage=<<"USAGE";
ProgramName:
Version: Version$version
Contact: Liu Min <lium\@biomarker.com.cn>
Program Date:   2013.9.16
Modify: none
Description: this program is used to assess the quality of RNA-SEQ raw data
Usage:
  Example: nohup perl rna_seq_data_assess.pl -Q 33 -config config.txt -outdir /path/to/outdir
  Options:
  -Q 		<int>   the cutoff of quality encode: 33 or 64 [forced]
  -config 	<file> 	the filename of config file [forced]
  -outdir 	<file> 	the directory of output [forced]
  -byraw 	<int> 	0: using fq? get seq file; 1: using raw_fq? get seq file [optional][default 0]
  -h 		<none> 	Help

USAGE
	print $usage;
	exit;
}










