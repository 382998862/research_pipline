#!/usr/bin/perl -w
# 
# Copyright (c) BIO_MARK 2014
# Writer:         Dengdj <dengdj@biomarker.com.cn>
# Program Date:   2014
# Modifier:       Dengdj <dengdj@biomarker.com.cn>
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
GetOptions(\%opts,"i=s","o=s","s=s","r=s","od=s","n=s","m=s","mode=s","min=s","max=s","queue=s","db=s","h" );

#&help()if(defined $opts{h});
if(!defined($opts{i}) || !defined($opts{s}) || !defined($opts{od}) || !defined($opts{mode}) || !defined($opts{r}) || !defined($opts{db}) || !defined($opts{o}) || defined($opts{h}))
{
	print <<"	Usage End.";
	Description:
		
		Version: $ver
		v1.2:    annotation, then select each sample.

	Usage:

		-i           vcf file                               <infile>                   must be given
		-od          outdir                                 <outdir>                   must be given
		-o           out prefix                             <outfile>                  must be given
		-r           ref file                               <infile>                   must be given
		-s           species                                <string>                   must be given
		-mode        SNP/Indel                              <string>                   must be given
		-db          database dir                           <dir>                      must be given
		-n           do annotation for given sample         [sample list file]         optional
		-m           max process for qsub                   [int]                      optional [25]
		-min         minimum indel length for draw          [int]                      optional [-10]
		-max         maximum indel length for draw          [int]                      optional [10]
        -queue       the queue is used for qsub jobs                                   optional [general.q]
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
my $vcffile = $opts{i} ;
$vcffile = &ABSOLUTE_DIR($vcffile) ;
my $outprefix = $opts{o} ;
my $outdir = $opts{od} ;
mkdir $outdir ;
$outdir = &ABSOLUTE_DIR($outdir) ;
my $database_dir = $opts{db} ;
$database_dir = &ABSOLUTE_DIR($database_dir);
my $species = $opts{s} ;
my $reffile = $opts{r} ;
$reffile = &ABSOLUTE_DIR($reffile);
my $mode = $opts{mode} ;
my $maxproc = defined $opts{m} ? $opts{m} : 25 ;
my $min_length = defined $opts{min} ? $opts{min} : -10 ;
my $max_length = defined $opts{max} ? $opts{max} : 10 ;
my $sample_list_file = defined $opts{n} ? $opts{n} : "" ;
my $queue = defined $opts{queue} ? $opts{queue} : "general.q" ;
my $shdir = "$outdir/sh_dir" ;
mkdir $shdir ;
#my $java_dir = "/usr/java/latest/bin/" ;
my $gatk = $config{gatk} ;
my $snpEff = "$Bin/v3_6/snpEff/snpEff.jar" ;
my $sEffConfig="$database_dir/snpEff.config";
# annotation for variants
my $anno_file = &annotation_for_variants($vcffile);

# get sample list
my @sample_files = ();
my @samples = ();

#if (!defined $opts{n}){
#	push @sample_files, $anno_file ;
#	@samples = ("all");
#}
#else{
	$sample_list_file = &ABSOLUTE_DIR($sample_list_file);
	my $asamples = &get_sample_files(\@sample_files, $anno_file, $sample_list_file,$queue);
	@samples = @{$asamples} ;
#}

# result statistic
&statistic_annotation_result(\@sample_files, $anno_file, $mode);


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



#my $asamples = &get_sample_files(\@sample_files, $vcffile, $sample_list_file);
sub get_sample_files()
{
	my ($asample_files, $vcffile, $sample_list_file,$queue) = @_ ;
	my $shfile = "$shdir/get_sample_vcf.sh" ;
	my @samples = ();
	open (IN, $sample_list_file) || die "Can't open $sample_list_file, $!\n" ;
	open (SH, ">$shfile") || die "Can't creat $shfile, $!\n" ;
	while(<IN>){
		chomp ;
		next if (m/# \#/ || m/# \s*$/);
		my $sample = (split)[0] ;
		my $dir = "$outdir/$sample" ;
		mkdir "$dir" ;
		my $basename = basename($vcffile);
		(my $outfile = "$dir/$basename") =~ s/.vcf$/.$sample.vcf/ ;
		print SH "java -XX:ParallelGCThreads=5 -jar $gatk -T SelectVariants -R $reffile --variant $vcffile -o $outfile -sn $sample -env \n" ;
		push @{$asample_files}, $outfile ;
		push @samples, $sample ;
	}

	push @{$asample_files},$vcffile ;
	push @samples, 'all' ;

	close(IN);
	close(SH);
	&qsubOrDie($shfile, $queue,10,"25G");

	return(\@samples) ;
}

#my $anno_file = &annotation_for_variants($vcffile);
sub annotation_for_variants()
{
	my ($vcffile) = @_ ;
	my $basename = basename($vcffile) ;
	(my $outfile = "$outdir/$basename") =~ s/.vcf$/.anno.vcf/ ;
	my $cmd = "" ;
	if (!-f $outfile){
		$cmd = "ln -s $vcffile $outdir" ;
		&run_or_die($cmd);
	}

	$cmd = "cd $outdir && java -XX:ParallelGCThreads=5 -Xmx10240m -jar $snpEff $species -v $vcffile -c $sEffConfig -o gatk > $outfile" ;
	&run_or_die($cmd);

	(my $gatk_anno = $outfile) =~ s/.vcf$/.gatk.vcf/ ;
	#$cmd = "java -XX:ParallelGCThreads=5 -jar $gatk -T VariantAnnotator -R $reffile -A SnpEff --variant $vcffile --snpEffFile $outfile -L $vcffile -o $gatk_anno" ;
	$cmd = "perl $Bin/extract_oneEff_anno.pl -i $outfile -o $gatk_anno" ;
	&run_or_die($cmd);

	return ($gatk_anno);
}

#&annotation_for_variants_old(\@sample_files, \@samples);

#&statistic_annotation_result($aanno_files, $anno_file, $mode);
sub statistic_annotation_result()
{
	my ($aanno_files, $anno_file, $mode) = @_ ;
	my $dir = "$outdir/anno_stat" ;
	mkdir $dir ;
	&link_to_dir($aanno_files, $dir);
	mkdir "$outdir/result" ;
	my $outfile = "$outdir/result/anno_stat.xls" ;
	my $cmd = "perl $Bin/anno_stat_v1.1.pl -id $dir -o $outfile" ;
	&run_or_die($cmd);

	# merge files
#	$outfile = "$outdir/result/$outprefix.combine.$mode.vcf" ;
#	$cmd = "java -XX:ParallelGCThreads=5 -jar $gatk -T CombineVariants -R $reffile -nt 8 -o $outfile";
#	for (my $i=0; $i<@{$aanno_files}; $i++){
#		my $vcffile = $aanno_files->[$i] ;
#		$cmd .= " -V $vcffile" ;
#	}
#	&run_or_die($cmd);
	# statistic indel
	if ($mode =~ /Indel/i){
		my $out_stat_file = "$outdir/result/$outprefix.combine.$mode" ;
		$cmd = "perl $Bin/indel_stat_v1.1.pl -i $anno_file -o $out_stat_file -m $min_length -x $max_length" ;
		&run_or_die($cmd);
	}

	return ;
}

#&link_to_dir($aanno_files, $dir);
sub link_to_dir()
{
	my ($aanno_files, $dir) = @_ ;
	for (my $i=0; $i<@{$aanno_files}; $i++){
		my $cmd = "ln -s ".$aanno_files->[$i]." $dir/" ;
		&run_or_die($cmd);
	}

	return ;
}

