#!/usr/bin/perl -w
use strict;
use warnings;
use Getopt::Long;
use Data::Dumper;
use newPerlBase;
use FindBin qw($Bin $Script);
use File::Basename qw(basename dirname);
my $BEGIN_TIME    = time();
my @Original_ARGV = @ARGV;
my $version="2.4.0";
#######################################################################################
# ------------------------------------------------------------------
# GetOptions, Check format, Output info.
# ------------------------------------------------------------------
my ( $cfg1, $od ,$queue,$cpu,$vf);
GetOptions(
    "help|?"        => \&USAGE,
    "cfg1:s"        => \$cfg1,
    "od:s"          => \$od,
    "queue:s"=>\$queue,
    "cpu:s"=>\$cpu,
    "vf:s"=>\$vf,
    
) or &USAGE;
&USAGE unless ( $cfg1 and $od );
$cfg1=&ABSOLUTE_DIR($cfg1);
&MKDIR($od);
$od = &ABSOLUTE_DIR($od);
&MKDIR("$od/Mapped");
&MKDIR("$od/work_sh");
$queue= $queue|| "general.q";
$cpu =$cpu || 50;
$vf =$vf|| "15G";
open(OUT,">$od/work_sh/Total_Reads.sh") or die $!;
open( IN, "$cfg1" ) || die "$!\n";
while (<IN>) {
    chomp;
    s/\r$//;
    s/^\s+//;
    s/\s+$//;
    next if ( /^\#/ || /^$/ );

    my @tmp = split /\s+/, $_;
    if ( $tmp[0] eq "Sample" ) {
        my $sam=$tmp[1];
        my $fq1 = <IN>;
        chomp $fq1;
        my $fq2 = <IN>;
        chomp $fq2;
        my $fq_1 = (split /\s+/, $fq1)[1];
        my $fq_2 = (split /\s+/, $fq2)[1];
        print OUT "wc -l $fq_1 |awk  '{print \"$sam\\t\"\$1/2}' > $od/Mapped/$sam.reads.txt \n";
    }

}
close IN;
close OUT;
qsubOrDie("$od/work_sh/Total_Reads.sh",$queue,$cpu,$vf);
qsubCheck("$od/work_sh/Total_Reads.sh");
`cat $od/Mapped/*.reads.txt > $od/totalRead.stat.xls`;
#my @tmp = glob ("$od/Mapped/*.reads.txt");
#open (OUT,">$od/totalRead.stat.xls");
#foreach my $k (@tmp){
#	open (IN,"$k");
#	while(<IN>){
#		chomp;
#		print OUT "$_\n";
#	}
#	close(IN);
#}
#close (OUT);



sub MKDIR {    # &MKDIR($out_dir);
    my ($dir) = @_;

    #	rmdir($dir) if(-d $dir);
    mkdir($dir) if ( !-d $dir );
}

sub ABSOLUTE_DIR {    #$pavfile=&ABSOLUTE_DIR($pavfile);
    my $cur_dir = `pwd`;
    chomp($cur_dir);
    my ($in) = @_;
    my $return = "";

    if ( -f $in ) {
        my $dir  = dirname($in);
        my $file = basename($in);
        chdir $dir;
        $dir = `pwd`;
        chomp $dir;
        $return = "$dir/$file";
    }
    elsif ( -d $in ) {
        chdir $in;
        $return = `pwd`;
        chomp $return;
    }
    else {
        warn "Warning just for file and dir in [sub ABSOLUTE_DIR]\n";
        exit;
    }

    chdir $cur_dir;
    return $return;
}

sub GetTime {
    my ( $sec, $min, $hour, $day, $mon, $year, $wday, $yday, $isdst ) =
      localtime( time() );
    return sprintf(
        "%4d-%02d-%02d %02d:%02d:%02d",
        $year + 1900,
        $mon + 1, $day, $hour, $min, $sec
    );
}

sub USAGE {    #
    my $usage = <<"USAGE";
Program: Tophat&Cufflinks_Analysis Procedure
Version: $version
Contact: Meng Fei <mengf\@biomarker.com.cn>

Description:
	This program is a Procedure deal with Multiple Samples RNA_Analysis
	Tophat+Cufflinks Combination£ºDesigned for RNA Analysis with a Reference Genome

	The program will calculate the Junctions, Transcripts(RABT) Assembly, FPKM of genes & isoforms

Usage:
	-cfg1	data config, rawdata & refseq path
	-queue	default general.q
	-cpu	default 50
	-vf	default 15G
	-od	output dir ,must be given
	
USAGE
    print $usage;
    exit;
}
