#!/usr/bin/perl -w
use strict;
use warnings;
use Getopt::Long;
use Data::Dumper;
use FindBin qw($Bin $Script);
use File::Basename qw(basename dirname);
use newPerlBase;
my $BEGIN_TIME=time();
my $version="1.0.0";
my %config=%{readconf("$Bin/../config/db_file.cfg")}; 

my @Original_ARGV=@ARGV;
#######################################################################################
# ------------------------------------------------------------------
# GetOptions, Check format, Output info.
# ------------------------------------------------------------------
my ($bamdir, $od,$queue,$medical);
GetOptions(
				"help|?" =>\&USAGE,
				"bamdir:s"  =>\$bamdir,
				"medical:s" =>\$medical,
				"queue:s"  =>\$queue,				
				"od:s"   =>\$od,
				) or &USAGE;
&USAGE unless ($bamdir  and $od) ;
################################

my $notename = `hostname`; chomp $notename;

$queue||="medical.q";
&MKDIR($od) unless (-d $od);
$od  = &ABSOLUTE_DIR($od);    
&MKDIR("$od/work_sh")  unless (-d "$od/work_sh");

my %total_read;
my %sample;

my @bam		=glob("$bamdir/*/*sorted.bam");
my $gff		=(glob("$bamdir/../Ref_Genome/*.gff*"))[0];
my $idx_prefix	=(glob("$bamdir/../Ref_Genome/*.fa"))[0];


open IN,"$bamdir/../totalRead.stat.xls" || die;
while (<IN>) {
	chomp;
	next if (/^$/);
	my @tmp=split/\s+/,$_;
	$total_read{$tmp[0]}=$tmp[1];
}
close(IN);
	
foreach my $path (sort @bam) {
	my $sam=basename(dirname($path));
	$sample{$sam}=$path;
}

#==================================================================
	open OUT1, ">$od/work_sh/Tophat_bam_stat.sh"   || die;
	open OUT2, ">$od/work_sh/genome_bam2depth.sh"  || die;
	open OUT3, ">$od/work_sh/genome_Checkgraph.sh" || die;
	open OUT4, ">$od/work_sh/plot_ReadDensity.sh"  || die;
	foreach my $sam (sort keys %sample) {
		print OUT1 "perl $Bin/basic_analysis/v1.2/tophat_cufflinks/bin/bam2map_stat.pl -i $sam -bam $sample{$sam} -totalRead $total_read{$sam} -od $od\n";
		print OUT2 "$config{bedtools} genomecov -ibam $sample{$sam} -g $idx_prefix.fai -d -strand +  |awk '\$3!=0' >$sample{$sam}.plus.depth && \n bedtools genomecov -ibam $sample{$sam} -g $idx_prefix.fai -d -strand -  |awk '\$3!=0' >$sample{$sam}.minus.depth &&\n samtools depth $sample{$sam} >$sample{$sam}.depth \n";
		print OUT3 "perl $Bin/basic_analysis/v1.2/tophat_cufflinks/bin/get_percent_of_exon_intro_inter_by_gff_v1.4.pl -gff $gff -i  $sample{$sam}.depth -od $od -index $sam \n";

		my $cmd="perl $Bin/basic_analysis/v1.2/tophat_cufflinks/bin/plotReadDensity2.pl -i $sample{$sam}.plus.depth:$sample{$sam}.minus.depth -o $od/ -k $sam ";
		$cmd.="-medical $medical "      if(defined $medical);
                print OUT4 "$cmd\n";
	}
	close OUT1;	close OUT2;	close OUT3;	close OUT4;
	&qsubOrDie("$od/work_sh/genome_bam2depth.sh", $queue ,  15, "15G");
	&qsubOrDie("$od/work_sh/genome_Checkgraph.sh",  $queue , 15, "15G");
	&qsubOrDie("$od/work_sh/plot_ReadDensity.sh",  $queue ,  15, "20G");
	&qsubOrDie("$od/work_sh/Tophat_bam_stat.sh", $queue,   15, "20G");

my $draw = "perl $Bin/basic_analysis/v1.2/gene_expression/detail_visualization/draw_total_random.pl -id $od -od $od";print "$draw\n";
`$draw`;


#######################################################################################
print STDOUT "\nDone. Total elapsed time : ",time()-$BEGIN_TIME,"s\n";
#######################################################################################
#==================================================================
# subs 
#==================================================================
sub LOAD_PARA {
	my $para_file= shift;
	my $para= shift;

	my $error_status = 0;
	open IN,$para_file || die "fail open: $para_file";
	while (<IN>) {
		chomp;
		s/^\s+//;s/\s+$//;s/\r$//;
		next if (/^$/ or /^\#/) ;
		my ($para_key,$para_value) = split(/\s+/,$_);
		$para->{$para_key} = $para_value;
		if (!defined $para_value) {
			warn "Non-exist: $para_key\n";
			$error_status = 1;
		}
	}
	close IN;
	die "\nExit due to error Undefined parameter\n" if ($error_status) ;
}


sub GetTMR {#
	#Get Total Mapped Reads from file line which contain string "Mapped Reads\s"
	my $fStat=shift;
	open (IN,"<",$fStat) or die $!;
	while (<IN>) {
		if (/^Mapped Reads\s(\d+)/) {
			close (IN) ;
			return $1;
		}
	}
	close (IN) ;
	die "Error Reads Stat file.\n";
}

sub MKDIR
{ # &MKDIR($out_dir);
	my ($dir)=@_;
	rmdir($dir) if(-d $dir);
	mkdir($dir) if(!-d $dir);
}

sub ABSOLUTE_DIR
{ #$pavfile=&ABSOLUTE_DIR($pavfile);
	my $cur_dir=`pwd`;chomp($cur_dir);
	my ($in)=@_;
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
		warn "Warning just for file and dir in [sub ABSOLUTE_DIR]\n";
		exit;
	}
	
	chdir $cur_dir;
	return $return;
}

sub GetTime {
	my ($sec, $min, $hour, $day, $mon, $year, $wday, $yday, $isdst)=localtime(time());
	return sprintf("%4d-%02d-%02d %02d:%02d:%02d", $year+1900, $mon+1, $day, $hour, $min, $sec);
}

sub USAGE {#
	my $usage=<<"USAGE";
Program: Tophat&Cufflinks_Analysis Procedure
Version: $version
Contact: Meng Fei <mengf\@biomarker.com.cn>

Description:
	Tophat Map stat  for RNA Analysis with a Reference Genome

Usage:
    -tophat         tophat dirname    must be given
    -queue        queue for qsub jobs. default :  general.q 
    -od                output dir            must be given
    -medical		medical project ,given GRCh37/GRCm38
    -help             help documents
USAGE
	print $usage;
	exit;
}
