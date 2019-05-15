#!/usr/bin/perl -w
use strict;
use Cwd qw(abs_path);
use Getopt::Long;
use Data::Dumper;
use FindBin qw($Bin $Script);
use File::Basename qw(basename dirname);
use newPerlBase;
my $BEGIN_TIME=time();
my $path = $Bin;
$path = substr($path,0,index($path,"bin"));
my %config=%{readconf("$path/project.cfg")};
my $Title=$config{Title};												#流程的名称，必填
my $version=$config{version};
# ------------------------------------------------------------------
# GetOptions
# ------------------------------------------------------------------
my ($bwa_dir, $ref,$gff,$odir);

GetOptions(
	"samFile:s" =>\$bwa_dir,
	"ref:s"=>\$ref,
	"gff:s"=>\$gff,
	"od:s"   =>\$odir,
	"help|h" =>\&USAGE,
	) or &USAGE;
&USAGE unless ($bwa_dir and $ref and $gff and $odir);

system "mkdir -p $odir" unless (-d $odir);
$odir=abs_path($odir);
$bwa_dir = abs_path $bwa_dir;
$ref = abs_path $ref;
$gff = abs_path $gff;
my $threads ||= 5;
&log_current_time("$Script start...");

# ------------------------------------------------------------------
# main pipline
# ------------------------------------------------------------------
# work shell backup dir
my $sh_dir = "$odir/work_sh";
mkdir $sh_dir unless (-d $sh_dir);
my $cmd;
my @bwa_files = glob("$bwa_dir/*/*.sam");
my $prefix = basename($gff);
$prefix=~s/gtf/txt/;
my $head = '#chrom\\tstart\\tend\\tname\\tscore\\tstrand\\tthickStart\\tthickEnd\\titemRgb\\texonCount\\texonSizes\\texonOffsets\\treadNumber\\tcircType\\tgeneName\\tisoformName\\tindex\\tflankIntron';
system "$config{gtfToGenePred} -genePredExt $gff $odir/refFlat.tmp.txt && awk 'BEGIN{FS=\"\t\"};{print \$12\"\t\"\$1\"\t\"\$2\"\t\"\$3\"\t\"\$4\"\t\"\$5\"\t\"\$6\"\t\"\$7\"\t\"\$8\"\t\"\$9\"\t\"\$10}' $odir/refFlat.tmp.txt > $odir/$prefix";
open SH,">$sh_dir/CIRCexplorer2.sh";
foreach my $file(@bwa_files){
	my $outcirc = basename($file);
	$outcirc =~s/.sam//;
	`rm -r $odir/$outcirc` if(-d "$odir/$outcirc");
	print SH "export PYTHONPATH=/share/nas1/songmm/software/CIRCexplorer2-2.1.0/CIRCexplorer2/lib/python2.7/site-packages/ && $config{CIRCexplorer} parse -t BWA $file -o $odir/$outcirc && $config{CIRCexplorer} annotate -r $odir/$prefix -g $ref $odir/$outcirc > $odir/CIRCexplorer2_annotate.log && sed -i 'N;2i\\$head' $odir/$outcirc/annotate/circ_fusion.txt && perl $Bin/filter_CIRC.pl -c $odir/$outcirc/annotate/circ_fusion.txt -o $odir/$outcirc/annotate/circ_fusion.txt.filter && mv $odir/$outcirc/annotate/circ_fusion.txt.filter $odir/$outcirc/annotate/circ_fusion.txt\n";
}
close SH;
`$config{qsub} --resource vf=50G --queue \"general.q\" --independent $sh_dir/CIRCexplorer2.sh --reqsub`;
#######################################################################################
my $elapse_time = (time()-$BEGIN_TIME)."s";
&log_current_time("$Script done. Total elapsed time: $elapse_time.");
#######################################################################################
# ------------------------------------------------------------------
# sub function
# ------------------------------------------------------------------
#############################################################################################################
sub data_cfg_read {
    &log_current_time("data config check:");
    my ($cfg_file, $data_cfg) = @_;
    my $sample_id;

    open (CFG, $cfg_file) or die "$!: $cfg_file\n";
    while (<CFG>) {
        chomp;
        next if (/^\s+/ or /^#/);
        if (/^Qphred/) {
            $data_cfg->{Qphred} = (split /\s+/,$_)[1];
        }
        if (/^Sample/) {
            $sample_id=(split/\s+/,$_)[1];
        }
        if ($_=~/^fq1/ || $_=~/^fq2/) {
            my $file=(split/\s+/,$_)[1];
            die "$file is not exist!\n" unless (-e $file);

            $data_cfg->{rawdata}{$sample_id}{fq1}=$file if $_=~/^fq1/;
            $data_cfg->{rawdata}{$sample_id}{fq2}=$file if $_=~/^fq2/;
        }
    }
    close CFG;

    if (defined $data_cfg->{Qphred}) {
        print "Qphred: $data_cfg->{Qphred}\n";
    } else {
        $data_cfg->{Qphred} = 33;
        print "Qphred: $data_cfg->{Qphred} [ASCII encoding type of quality score of rawdata is unknown, and default is 33.]\n";
    }

    $data_cfg->{sample_num} = scalar keys %{$data_cfg->{rawdata}};
    print "sample_number: $data_cfg->{sample_num}\n";

    for my $s (sort keys %{$data_cfg->{rawdata}}) {
        print "${s}_fq1: $data_cfg->{rawdata}{$s}{fq1}\n${s}_fq2: $data_cfg->{rawdata}{$s}{fq2}\n";
    }
    &log_current_time("data config check done.\n");
}

#############################################################################################################

sub creat_new_config {
	my $config = shift;
	my $key = shift;
	my $fq1 = shift;
	my $fq2 = shift;
	my $para = shift;
	my $conf = shift;

	$conf .= "Index\t$key\n";
	$conf .= "FQ1\t$fq1\n";
	$conf .= "FQ2\t$fq2\n\n";

	foreach my $para_meter (sort keys %{$para}) {
		$para_meter =~ /^para_(.*)/;
		$conf .= "$1\t$para->{$para_meter}\n";
	}
	open OUT,">$config" || die $!;
	print OUT "$conf";
	close OUT;
}

#############################################################################################################
sub detail_cfg_read {
    &log_current_time("detail config check:");
    my ($cfg_file, $detail_cfg) = @_;

    open (CFG,$cfg_file ) or die "$!: $cfg_file\n";
    while (<CFG>) {
        chomp;
        next if (/^\s+/ or /^#/);
        my ($key, $value) = (split /\s+/)[0,1];

        if ($key eq 'Project_name' or $key eq 'Customer_info' or $key eq 'Project_id' or $key eq 'Project_key') {
            $detail_cfg->{$key} = $value;
        }
        if ($key eq 'Known_unigene' or $key eq 'Known_pep') {
            die "$key: $value is not exist!\n" unless (-e $value);
            $detail_cfg->{$key} = $value;
        }
        if ($key eq 'Known_anno') {
            die "$key: $value is not illegal!\n" unless (-e "$value/02.gene-annotation" and -e "$value/Result");
            $detail_cfg->{$key} = $value;
        }
        if ($key eq 'Ref_seq' or $key eq 'Ref_ann') {
            die "$key: $value is not exist!\n" unless (-e $value);
            $detail_cfg->{$key} = $value;
        }
        if ($key =~/^SNP_/) {
            $detail_cfg->{$key} = $value;
        }
        if ($key eq 'nr' or $key eq 'Swissprot' or $key eq 'Kegg' or $key eq 'Pfam' or $key eq 'Cog' or $key eq 'Kog') {
            die "$key: $value is not exist!\n" unless (-e $value);
            $detail_cfg->{anno}{$key} = $value;
        }

        print "$key: $value\n" if (exists $detail_cfg->{$key});
    }
    close CFG;

    &log_current_time("detail config check done.");
}

#############################################################################################################
sub step_cmd_process {
    my ($cmd, $step,$sh_dir) = @_;
    my $sh_file = "$sh_dir/$step.sh";
    my $log_file = "$sh_file.log";
    my $flag = 0;
    my $start_time = time();
    &log_current_time("$step step start ...");
    &log_current_time("CMD: $cmd");

    if (-e $sh_file) {
        system "cat $sh_file >> $sh_file.bak";
        open (SH, ">$sh_file") or die "$!: $sh_file\n";
        print SH "$cmd\n";
        close SH;
    } else {
        open (SH, ">$sh_file") or die "$!: $sh_file\n";
        print SH "$cmd\n";
        close SH;
    }

    $flag = system("$cmd > $log_file");
    if ($flag != 0){
        log_current_time("Error: command failed: $cmd");
        exit(1);
    } else {
        my $escaped_time = (time()-$start_time)."s";
        &log_current_time("$step step done, escaped time: $escaped_time.\n");
    }
}

#############################################################################################################
sub log_current_time {
     # get parameter
     my ($info) = @_;

     # get current time with string
     my $curr_time = date_time_format(localtime(time()));

     # print info with time
     print "[$curr_time] $info\n";
}

#############################################################################################################
sub date_time_format {
    my ($sec, $min, $hour, $day, $mon, $year, $wday, $yday, $isdst)=localtime(time());
    return sprintf("%4d-%02d-%02d %02d:%02d:%02d", $year+1900, $mon+1, $day, $hour, $min, $sec);
}



#############################################################################################################
sub USAGE {
	my $usage=<<"USAGE";
----------------------------------------------------------------------------------------------------------------------------
   Program: $Script
   Version: $version
   Contact: Simon Young <yangxh\@biomarker.com.cn>
      Date: 2014-11-13

     Usage:
            --samFile      <FILE>  sam file directory
            --ref      <FILE>  reference geneome fasta file
            --gff      <FILE>    annotation file
            --re       <PAR>
			--od       <DIR>     result dir
            --h                 help documents

   Example:
            perl $Script --samFile bwa_dir --ref *.fa --gff *.gtf --od CircRNA_identify/

----------------------------------------------------------------------------------------------------------------------------
USAGE
	print $usage;
	exit;
}
