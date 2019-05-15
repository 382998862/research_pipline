#!/usr/bin/perl -w
use strict;
use warnings;
use Getopt::Long;
use Data::Dumper;
use FindBin qw($Bin $Script);
use File::Basename qw(basename dirname);
use newPerlBase;
my %config=%{readconf("$Bin/../../lncRNA_pip.cfg")};
#######################################################################################
# ------------------------------------------------------------------
# GetOptions, Check format, Output info.
# ------------------------------------------------------------------
my (  $cfg2, $od,$in );
GetOptions(
    "help|?"        => \&USAGE,
    "cfg2:s"        => \$cfg2,
    "in:s"	    => \$in,
    "od:s"          => \$od,
) or &USAGE;
&USAGE unless ( $cfg2 and $in and $od );
################################

$cfg2 = &ABSOLUTE_DIR($cfg2);
&MKDIR($od);
$od = &ABSOLUTE_DIR($od);
$in = &ABSOLUTE_DIR($in);
&MKDIR("$od/work_sh");

#==================================================================
# load config file
#==================================================================

my %para;
open(CFG,$cfg2)||die $!;
while(<CFG>){
	chomp;next if($_=~/^#|^$|^\s+/);
	my @tmp=split(/\s+/,$_);
	$para{$tmp[0]}=$tmp[1];
}
close(CFG);

    $para{Memory}     ||= "15G";
    $para{Queue_type} ||= "general.q";
    $para{CPU} ||= "30";

    open OUT2, ">$od/work_sh/genome_bam2depth.sh"  || die;
    open OUT3, ">$od/work_sh/genome_Checkgraph.sh" || die;
    open OUT4, ">$od/work_sh/plot_ReadDensity.sh"  || die;
    my @str_type_stat;
    
    my @bams=glob("$in/Hisat/*/*HISAT_aln.sorted.bam");
    my @samples=map{(split(/\./,basename $_))[0]} @bams;
    my $fai=(glob("$in/Ref_Genome/*.fai"))[0]; 
    foreach my $sam ( sort{$a cmp $b} @samples) {
        push @str_type_stat, "$od/$sam.type.stat";

        print OUT2 "$config{samtools} depth $in/Hisat/$sam/$sam.HISAT_aln.sorted.bam >$in/Hisat/$sam/$sam.HISAT_aln.sorted.bam.depth &&\n";
        print OUT3 "perl $Bin/bin/get_percent_of_exon_intro_inter_by_gff_v1.4.pl -gff $para{Ref_ann} -i $in/Hisat/$sam/$sam.HISAT_aln.sorted.bam.depth -od $od -index $sam  &&  ";    #randcheck
        print OUT3 "perl $Bin/bin/draw_total_random.pl -id $od -od $od  \n";
        print OUT4 "$config{Rscript} $Bin/bin/pie.R infile=$od/$sam.type.stat outfile=$od/$sam.type.png legend.col=1  value.col=2 skip=1 sep=t\n";    #2015/08/10,modify bu niulg,type
	print OUT4 "perl $Bin/bin/plotReadDensity.2.pl -q $para{Queue_type} -a $fai -f bam -i $in/Hisat/$sam/$sam.HISAT_aln.sorted.bam -o $od -k $sam";
	print OUT4 " -medical $para{medical} "	if(exists $para{medical});
	print OUT4 "\n";
    }
    close OUT2;
    close OUT3;
    close OUT4;

&qsub("$od/work_sh/genome_bam2depth.sh");
&qsub("$od/work_sh/genome_Checkgraph.sh");
&qsub("$od/work_sh/plot_ReadDensity.sh");

&run_or_die("$config{Rscript} $Bin/bin/fpkm_saturation_Qv1.0.R --count $in/prepDE/All_gene_counts.list --fpkm $in/prepDE/All_gene_fpkm.list --o $in/Mapped/");
#"$config{Rscript} $Bin/bin/fpkm_saturation_v3.R --i /share/nas1/wanghao/project/BMK160104-B225-0201-20-Human-fats-Lnc/Analysis-lnc/Basic_Analysis/geneExpression/L01.geneExpression.xls --o $od/$sam.Saturation.png";

my $str_type_stat = join " ", @str_type_stat;
&run_or_die("perl $Bin/bin/map_total_type.pl -i $str_type_stat -o $od/Total_Type.png");
 

#==================================================================
# subs
#==================================================================
sub qsub()
{
        my $shfile= shift;
	my $queue="medical.q";
	$queue=$para{Queue_type}	if(exists $para{Queue_type});
        my $cmd = "sh /share/nas2/genome/bmksoft/tool/qsub_sge_plus/v1.0/qsub_sge.plus.sh --reqsub --independent $shfile --queue $queue ";
        &run_or_die($cmd);
        return ;
}
sub run_or_die()
{
        my ($cmd) = @_ ;
        &show_log($cmd);
        my $flag = system($cmd) ;
        if ($flag != 0){
                &show_log("Error: command fail: $cmd");
                exit(1);
        }
        &show_log("done.");
        return ;
}
sub show_log()
{
        my ($txt) = @_ ;
        my $time = time();
        my($sec, $min, $hour, $day, $mon, $year, $wday, $yday, $isdst) = localtime($time);
        $wday = $yday = $isdst = 0;
        my $Time=sprintf("%4d-%02d-%02d %02d:%02d:%02d", $year+1900, $mon+1, $day, $hour, $min, $sec);
        print "$Time:\t$txt\n" ;
        return ($time) ;
}


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
Contact: Meng Fei <mengf\@biomarker.com.cn>

Description:
	The program will statistic the bam file and draw the picture 

Usage:				bam file  Map Stat
    -cfg2             detail config, analysis parameters
    -od               output dir Map_stat            	    must be given
    -in               input dir Hisat_Stringtie             must be given
	
                            
USAGE
    print $usage;
    exit;
}
