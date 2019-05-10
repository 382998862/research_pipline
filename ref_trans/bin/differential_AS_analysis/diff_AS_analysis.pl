#!/usr/bin/perl -w
use Getopt::Long;
use Getopt::Std;
use Config::General;
use FindBin qw($Bin $Script);
use newPerlBase;
use Cwd qw(abs_path getcwd);
use threads::shared qw($threads);
use File::Basename qw(basename dirname);
#######################################################################################

# ------------------------------------------------------------------
# GetOptions
# ------------------------------------------------------------------
my ($version,$GTF,$od,$t,$readLength,$cfg,$anchorLength,$libType,$cstat,$tstat,$nthread,$statoff,$bamdir);
$version="4.0.2";

GetOptions(
        "h|?"           =>\&USAGE,
	"version:s"	=>\$version,
	"cfg2:s"		=>\$cfg,
	"od:s"     	=>\$od,
	"bamdir:s"      =>\$bamdir,
	
)or &USAGE;
&USAGE unless ($cfg and $od and $bamdir);

#########rMATS parameter
$cfg=abs_path($cfg);
$anchorLength||=1;
$cstat||=0.0001;
if ($cstat<0 || $cstat>1) {
	print "ERROR: --cstat parameter should be: 0 ≤ cstat < 1";
}

$nthread||=1;
$tstat = $nthread;
$libType||="fr-unstranded";
$t||="paired";
$readLength||=150;

`mkdir $od` unless -d $od;
$od=abs_path($od);

##############rMATS standard parameter

my $gffread="/share/nas2/genome/bin/gffread";
my @coms=();my @seps=();
open(CFG,$cfg)||die $!;
while(<CFG>){
	chomp;next if(/^#/);
	if(/^Com/){
		my @tmp=split(/\s+/,$_);
		push @coms,$tmp[1];
	}elsif(/^Sep/){
		my @tmp=split(/\s+/,$_);
		push @seps,$tmp[1];
	}
	if(/^Ref_ann/){
		my $GFF=(split/\s+/,$_)[1];$GFF=abs_path($GFF);
		open(GTF,">$od/genome.gtf")||die $!;
		`$gffread $GFF -T -o $od/genome.gtf`;
		close(GTF);
	}
	if(/^Lib_type/){
		$libType=(split/\s+/,$_)[1];
	}
	if (/^ReadLength/) {
                $readLength=(split/\s+/,$_)[1];
        }
        if (/^Readtype/) {
                $t=(split/\s+/,$_)[1];
        }

}
close(CFG);


$GTF="$od/genome.gtf";$GTF=abs_path($GTF);
my $select="--gtf $GTF -t $t --readLength $readLength";
################rMATS analysis

my $Python="/share/nas2/genome/biosoft/Python/2.7.8/bin/python";
my $rMATS="/share/nas2/genome/biosoft/Python/rMATS/4.0.2/rMATS-turbo-Linux-UCS2";

mkdir "$od/work_sh" unless -d "$od/work_sh";

open(MATS,">$od/work_sh/rmats.sh")||die $!;
if(@coms>0){
	#my $odir="$od/diff_AS";
	foreach my $com(@coms){
		my $bam_file=&getdecfg($com);
		if(!defined $statoff) {
		print MATS "$Python $rMATS/rmats.py $select --b1 $od/$bam_file/b1.txt --b2 $od/$bam_file/b2.txt --od $od/$bam_file --anchorLength $anchorLength --libType $libType --cstat $cstat --tstat $tstat --nthread $nthread \n";
		}else{
		print MATS "$Python $rMATS/rmats.py $select --b1 $od/$bam_file/b1.txt --b2 $od/$bam_file/b2.txt --od $od/$bam_file --anchorLength $anchorLength --libType $libType --cstat $cstat --tstat $tstat --nthread $nthread --statoff \n";
		}
	}
}
if(@seps>0){
	#my $odir="$od/diff_AS";
	foreach my $sep(@seps){
		my $bam_file=&getdecfg($sep);
		if(!defined $statoff) {
                print MATS "$Python $rMATS/rmats.py $select --b1 $od/$bam_file/b1.txt --b2 $od/$bam_file/b2.txt --od $od/$bam_file --anchorLength $anchorLength --libType $libType --cstat $cstat --tstat $tstat --nthread $nthread \n";
		}else{
		print MATS "$Python $rMATS/rmats.py $select --b1 $od/$bam_file/b1.txt --b2 $od/$bam_file/b2.txt --od $od/$bam_file --anchorLength $anchorLength --libType $libType --cstat $cstat --tstat $tstat --nthread $nthread --statoff \n";
		}
	}	
}
close(MATS);
&qsub("$od/work_sh/rmats.sh");

#######################rMATS result processing
#my @JCEC_file=glob("$od/");
##########
sub getdecfg{
        my $vs=shift;
        my (@g1,@g2,$b1,$b2);
	my %sample;
	#my $bamdir="$od/Mapping/Hisat";
	my @bam=glob("$bamdir/*/*sorted.bam");
	foreach my $path (sort @bam) {
        	my $sam=basename(dirname($path));
        	$sample{$sam}=$path;
	}	
        if($vs=~/\;/){
		my @path_b1;my @path_b2;
                my @tmp=split(/\;/,$vs,2);
                @g1=split(/,/,$tmp[0]);
		foreach my $sam(sort @g1) {
			my $each_path=$sample{$sam};
			my $each_bam=abs_path($each_path);
			push @path_b1,$each_bam;
		}
		$b1=join(",",@path_b1);
                @g2=split(/,/,$tmp[1]);
		foreach my $sam(sort @g2) {
                        my $each_path=$sample{$sam};
			my $each_bam=abs_path($each_path);
                        push @path_b2,$each_bam;
                }
		$b2=join(",",@path_b2);
        }else{
                my @tmp=split(/,/,$vs,2);
		$b1=abs_path($sample{$tmp[0]});
		$b2=abs_path($sample{$tmp[1]});	
                push @g1,$tmp[0];
                push @g2,$tmp[1];
        }
        my $base=join("_",@g1)."_vs_".join("_",@g2);

        `mkdir -p $od/$base`       unless(-d "$od/$base");
        open(B1,">$od/$base/b1.txt")||die $!;
        print B1 $b1,"\n";
        close (B1);
        open(B2,">$od/$base/b2.txt")||die $!;
        print B2 $b2,"\n";
        close (B2);
        return $base;
}

sub qsub{
        my $shfile= shift;
	my $queue="medical.q";
        my $cmd = "sh /share/nas2/genome/bmksoft/tool/qsub_sge_plus/v1.0/qsub_sge.plus.sh --reqsub --independent $shfile --queue $queue ";
        &run_or_die($cmd);              
        return ;
}
sub run_or_die{
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
sub show_log{
        my ($txt) = @_ ;
        my $time = time();
        my($sec, $min, $hour, $day, $mon, $year, $wday, $yday, $isdst) = localtime($time);
        $wday = $yday = $isdst = 0;
        my $Time=sprintf("%4d-%02d-%02d %02d:%02d:%02d", $year+1900, $mon+1, $day, $hour, $min, $sec);
        print "$Time:\t$txt\n" ;
        return ($time) ;
}
sub USAGE{
	my $usage=<<"USAGE";

Description: Multivariate Analysis of Transcript Splicing (MATS)
Version: 1.0.0
Contact: lijj\@biomarker.com.cn
Usage:
	

	-h		help;
	--version	Version 4.0.2;
	--cfg2	file	Ref_trans detail config file, forced;
	--od	path	forced; output folder of post step;
	--bamdir    <DIR>   bam dir;
	
Example:

	perl diff_AS_analysis.pl --cfg ref_trans.detail.cfg --bamdir Tophat/Hisat_bamdir_path --od outDir 
USAGE
	print $usage;
	exit;
}
