#!/usr/bin/perl -w
# 
# Copyright (c) BIO_MARK 2014
# Writer:         Dengdj <dengdj@biomarker.com.cn>
# Program Date:   2014
# Modifier:       Dengdj <dengdj@biomarker.com.cn>
# Last Modified:  2014.
# Modifier:       Wangyj <wangyj@biomarker.com.cn>
# Last Modified:  2014.
my $ver="1.2";

use strict;
use Getopt::Long;
use Data::Dumper;
use FindBin qw($Bin $Script);
use File::Basename qw(basename dirname);
use newPerlBase;
my %config=%{readconf("$Bin/../../../../config/db_file.cfg")};

######################请在写程序之前，一定写明时间、程序用途、参数说明；每次修改程序时，也请做好注释工作

my %opts;
GetOptions(\%opts,"s=s","r=s","id=s","queue=s","h" );

#&help()if(defined $opts{h});
if(!defined($opts{id}) || !defined($opts{s}) || !defined($opts{r}) || !defined($opts{queue})|| defined($opts{h}))
{
	print <<"	Usage End.";
	Description:
		
		Version: $ver
		v1.2:    annotation, then select each sample.

	Usage:

		-id          inputdir                               <infile>                   must be given
		-r           ref file                               <infile>                   must be given
        -queue       the queue is used for qsub jobs        <string>                   must be given
        -s           species                                <string>                   must be given
		-h           Help document

	Usage End.

	exit;
}

###############Time
my $Time_Start;
$Time_Start = sub_format_datetime(localtime(time()));
print "\nStart Time :[$Time_Start]\n\n";
################
# get parameter
my $snpdir = $opts{id} ;
$snpdir = &ABSOLUTE_DIR($snpdir) ;
my $species = $opts{s} ;
my $reffile = $opts{r} ;
my $queue=$opts{queue} ;
$reffile = &ABSOLUTE_DIR($reffile);

my $snpEff = "$Bin/v3_6/snpEff/snpEff.jar" ;
my $Rscript = $config{Rscript};
###############Start
`mkdir -p $snpdir/data/$species `  unless (-d "$snpdir/data/$species" ) ;
`ln -s $reffile $snpdir/data/$species/sequences.fa ` unless (-e "$snpdir/data/$species/sequences.fa" ) ;
`ln -s $snpdir/Integrated.gff $snpdir/data/$species/genes.gff ` unless (-e "$snpdir/data/$species/genes.gff" ) ;

open CONFIG, ">$snpdir/snpEff.config" or die $!;
	print CONFIG "#$species\n$species",'.genome:',"$snpdir/data/$species/sequences.fa";
close CONFIG;


open LIST,">$snpdir/SNP/sample.list" or die $!;
	my @samples=glob("$snpdir/aligndir/*.bam");
	foreach my $sam (@samples) {
		chomp$sam;
		my $temp=basename($sam);
		$temp=~s/\.bam//;
		print LIST "$temp\n";
	}
close LIST;

my $cmd=" java -XX:ParallelGCThreads=5 -jar  $snpEff  build  -v  -c  $snpdir/snpEff.config  -gff3 $species >$snpdir/data/$species/build.log 2>&1";
print $cmd,"\n";
system $cmd;

my $fa=(glob("$snpdir/SNP/Genome/*.fa"))[0];

$cmd="perl $Bin/SNPAnnotation_v1.2.pl -i $snpdir/SNP/final.snp.vcf -od $snpdir/SNP/snp_anno -o $species -r  $fa  -s $species -mode SNP -db $snpdir -n $snpdir/SNP/sample.list ";
print $cmd,"\n";
system $cmd;



$cmd="perl $Bin/SNPAnnotation_v1.2.pl -i $snpdir/SNP/final.indel.vcf -od $snpdir/SNP/indel_anno -o $species -r  $fa  -s $species -mode Indel -db $snpdir -n $snpdir/SNP/sample.list -queue $queue ";
print $cmd,"\n";
system $cmd;



`mkdir -p $snpdir/stat/snp_anno `  unless (-d "$snpdir/stat/snp_anno");
`mkdir -p $snpdir/stat/indel_anno `  unless (-d "$snpdir/stat/indel_anno");
system "cp $snpdir/SNP/snp_anno/result/anno_stat.xls  $snpdir/stat/snp_anno/final_SNP.anno.stat ";
system "cp $snpdir/SNP/indel_anno/result/anno_stat.xls $snpdir/stat/indel_anno/final_Indel.anno.stat ";
	
$cmd="$Rscript $Bin/snp_indel_anno_plot.r $snpdir/stat/snp_anno/final_SNP.anno.stat $snpdir/stat/snp_anno  snp";
print $cmd,"\n";
system $cmd;

$cmd="$Rscript $Bin/snp_indel_anno_plot.r $snpdir/stat/indel_anno/final_Indel.anno.stat  $snpdir/stat/indel_anno indel";
print $cmd,"\n";
system $cmd;

my @vcffile=glob("$snpdir/SNP/snp_anno/anno_stat/*.vcf");
foreach my $vcfs (@vcffile) {
	my $outfile = basename$vcfs;
	$outfile=~s/\.vcf$//;
	if ($outfile=~/gatk$/){
		`cp  $vcfs $snpdir/stat ` unless (-e "$snpdir/stat/$outfile.vcf");
		$outfile.='.all';
		$cmd="perl $Bin/vcf_to_snplist_v1.5.1.pl -i $vcfs -o $snpdir/stat/$outfile.list -ref 1 ";
		print $cmd,"\n";
		system $cmd;
	}
	else {
	$cmd="perl $Bin/vcf_to_snplist_v1.5.1.pl -i $vcfs -o $snpdir/stat/snp_anno/$outfile.list -ref 1 ";
	print $cmd,"\n";
	system $cmd;
	}
}


@vcffile=glob("$snpdir/SNP/indel_anno/anno_stat/*.vcf");
foreach my $vcfs (@vcffile) {
	my $outfile = basename$vcfs;
	$outfile=~s/\.vcf$//;
	if ($outfile=~/gatk$/){
		`cp $vcfs $snpdir/stat ` unless (-e "$snpdir/stat/$outfile.vcf");
		$outfile.='.all';
		$cmd="perl $Bin/vcf_to_indellist_v1.4.pl -i $vcfs -o $snpdir/stat/$outfile.list -ref 1 ";
		print $cmd,"\n";
		system $cmd;
	}
	else{
	$cmd="perl $Bin/vcf_to_indellist_v1.4.pl -i $vcfs -o $snpdir/stat/indel_anno/$outfile.list -ref 1 ";
	print $cmd,"\n";
	system $cmd;
	}
}


###############Time
my $Time_End;
$Time_End = sub_format_datetime(localtime(time()));
print "\nEnd Time :[$Time_End]\n\n";

###############Subs
sub sub_format_datetime {#Time calculation subroutine
    my($sec, $min, $hour, $day, $mon, $year, $wday, $yday, $isdst) = @_;
	$wday = $yday = $isdst = 0;
    sprintf("%4d-%02d-%02d %02d:%02d:%02d", $year+1900, $mon+1, $day, $hour, $min, $sec);
}

sub ABSOLUTE_DIR
{ #$pavfile=&ABSOLUTE_DIR($pavfile);
	my $cur_dir=`pwd`;
	$cur_dir =~ s/\n$//;
	my ($in)=@_;
	my $return="";
	
	if(-f $in)
	{
		my $dir=dirname($in);
		my $file=basename($in);
		chdir $dir;$dir=`pwd`;
		$dir =~ s/\n$// ;
		$return="$dir/$file";
	}
	elsif(-d $in)
	{
		chdir $in;$return=`pwd`;
		$return =~ s/\n$// ;
	}
	else
	{
		warn "Warning just for file and dir\n";
		exit;
	}
	
	chdir $cur_dir;
	return $return;
}

# &show_log("txt")
sub show_log()
{
	my ($txt) = @_ ;
	my $time = time();
	my $Time = &sub_format_datetime(localtime($time));
	print "$Time:\t$txt\n" ;
	return ($time) ;
}

#&run_or_die($cmd);
sub run_or_die()
{
	my ($cmd) = @_ ;
	my $start_time = &show_log($cmd);
	my $flag = system($cmd) ;
	if ($flag != 0){
		my $end_time = &show_log("Error: command fail: $cmd");
		&show_log("Elaseped time: ".($end_time - $start_time)."s\n");
		exit(1);
	}
	my $end_time = &show_log("done.");
	&show_log("Elaseped time: ".($end_time - $start_time)."s\n");
	return ;
}




