#!/usr/bin/perl -w
use strict;
use Cwd qw(abs_path);
use Getopt::Long;
use Data::Dumper;
use FindBin qw($Bin $Script);
use File::Basename qw(basename dirname);
use newPerlBase;
my $BEGIN_TIME=time();
my %config=%{readconf("$Bin/../../project.cfg")};
my $Title=$config{Title};												#流程的名称，必填
my $version=$config{version};
# ------------------------------------------------------------------
# GetOptions
# ------------------------------------------------------------------
my ($datacfg,$detailcfg,$odir);

GetOptions(
	"datacfg:s" =>\$datacfg,
	"detailcfg:s"=>\$detailcfg,
    	"od:s"   =>\$odir,
    	"help|h" =>\&USAGE,
    ) or &USAGE;
&USAGE unless ($datacfg and $detailcfg and $odir);

system "mkdir -p $odir" unless (-d $odir);
$odir=abs_path($odir);
my $thread ||= 10;
&log_current_time("$Script start...");
my (%data_cfg, %detail_cfg);
&data_cfg_read($datacfg,\%data_cfg);
&detail_cfg_read($detailcfg,\%detail_cfg);

my $sh_dir = "$odir/work_sh";
mkdir $sh_dir unless (-d $sh_dir);
my $bowtie2_index = "$odir/bowtie2_index";
mkdir $bowtie2_index unless (-d $bowtie2_index);
`cd $bowtie2_index && $config{bowtie2_build} $detail_cfg{Ref_seq} REF &>bt2_build.log && touch bowtie2_build.finish` if(!-e "$bowtie2_index/bowtie2_build.finish");
my $ref = basename($detail_cfg{Ref_seq});
mkdir "$bowtie2_index/chr" unless (-d "$bowtie2_index/chr");
&cut_fa_file_to_dir($detail_cfg{Ref_seq},"$bowtie2_index/chr");

open (SH1,">$sh_dir/S1.Bowtie_Mapping.sh");
open (SH2,">$sh_dir/S2.Unmapped.sh");
open (SH3,">$sh_dir/S3.Candidate_CircRNA.sh");
my $head = '#chrom\\tstart\\tend\\tname\\tn_reads\\tstrand\\tn_uniq\\tuniq_bridges\\tbest_qual_left\\tbest_qual_right\\ttissues\\ttiss_counts\\tedits\\tanchor_overlap\\tbreakpoints\\tsignal\\tstrandmatch\\tcategory';

foreach my $sample(sort keys %{$data_cfg{rawdata}}){
	my $result = "$odir/${sample}";
	mkdir $result unless (-d $result);
	my $fq1 = $data_cfg{rawdata}{$sample}{fq1};
	my $fq2 = $data_cfg{rawdata}{$sample}{fq2};
	print SH1 "cd $result && ln -snf $detail_cfg{Ref_seq} ./ && $config{bowtie2} -p3 --very-sensitive --mm --score-min=C,-15,0 -x $bowtie2_index/REF --rf -q -U $fq1,$fq2 2> bt2_firstpass.log | $config{samtools} view -hbuS - | $config{samtools} sort - $sample \n";
	print SH2 "cd $result && $config{samtools} view -hf 4 $sample.bam | $config{samtools} view -Sb - > $sample.unmapped.bam && ";
	print SH2 "$config{python} $Bin/bin/find_circ/unmapped2anchors.py $sample.unmapped.bam > $sample.qfa \n";
	print SH3 "cd $result && $config{bowtie2} -p3 --reorder --mm --score-min=C,-15,0 -q -x $bowtie2_index/REF --rf -U $sample.qfa 2> bt2_secondpass.log | $config{python} $Bin/bin/find_circ/find_circ.py -G $ref -p 5 -n $sample -s sites.log -R $sample.spliced_reads.fa > $sample.splice_sites.bed && ";
	print SH3 "grep CIRCULAR $sample.splice_sites.bed |grep -v chrM |awk '\$5>=2' |grep UNAMBIGUOUS_BP | grep ANCHOR_UNIQUE |$config{python} $Bin/bin/find_circ/maxlength.py 100000 > circ_candidates.bed && sed -i 'N;2i\\$head' circ_candidates.bed \n";
}
close(SH1);
close(SH2);
close(SH3);
`$config{qsub} --resource vf=50G --queue $detail_cfg{Queue_type} --independent $sh_dir/S1.Bowtie_Mapping.sh --reqsub`;
&qsubCheck("$sh_dir/S1.Bowtie_Mapping.sh");
`$config{qsub} --resource vf=50G --queue $detail_cfg{Queue_type} --independent $sh_dir/S2.Unmapped.sh --reqsub`;
&qsubCheck("$sh_dir/S2.Unmapped.sh");
`$config{qsub} --resource vf=50G --queue $detail_cfg{Queue_type} --independent $sh_dir/S3.Candidate_CircRNA.sh --reqsub`;
&qsubCheck("$sh_dir/S3.Candidate_CircRNA.sh");


#######################################################################################
my $elapse_time = (time()-$BEGIN_TIME)."s";
&log_current_time("$Script done. Total elapsed time: $elapse_time.");
#######################################################################################
# ------------------------------------------------------------------
# sub function
# ------------------------------------------------------------------
#############################################################################################################
sub cut_fa_file_to_dir
{
	my($mRNA,$odir)=@_;
	$/=">";
	open IN,$mRNA;
	while(<IN>)
	{
		chomp;
		next if /^\s*$/;
		my @line = split/\n/,$_,2;
		my $name = $line[0];
		$name=~/(\w+)\s*/;
		open OUT,">$odir/$1.fa";
		print OUT "$_\n";
		close OUT;
	}
}
sub data_cfg_read {
    &log_current_time("data config check:");
    my ($cfg_file, $data_cfg) = @_;
    my $sample_id;

    open (CFG, $cfg_file) or die "$!: $cfg_file\n";
    while (<CFG>) {
        chomp;
        next if (/^\s+/ or /^#/ or /^$/);
    	my @tmp=split(/\s+/,$_);
	if ($tmp[0] eq "Qphred") {
            $data_cfg->{Qphred} = $tmp[1];
        }
        if ($tmp[0] eq "Sample") {
            $sample_id=$tmp[1];
        }
        if ($tmp[0] eq "fq1" || $tmp[0] eq "fq2") {
            my $file=$tmp[1];
            die "$file is not exist!\n" unless (-e $file);

            $data_cfg->{rawdata}{$sample_id}{fq1}=$file if $tmp[0] eq "fq1";
            $data_cfg->{rawdata}{$sample_id}{fq2}=$file if $tmp[0] eq "fq2";
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
        next if (/^\s+/ or /^#/ or /^$/);
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
	if ($key eq 'Queue_type' or $key eq 'Queue_type1'){
	    die "$key: $value is not exist!\n" unless (defined $value);
	    $detail_cfg->{$key} = $value;	
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
sub qsubCheck{
        my $sh = shift;
        my @Check_file = glob "$sh*.qsub/*.Check";
        my @sh_file    = glob "$sh*.qsub/*.sh";
	if(@sh_file==0){
		print "qsub is abnormal,Please Check whether system is abnormal or not..\n";
		die;
	}
        if ( $#sh_file != $#Check_file ) {
                print "Their Some Error Happend in $sh qsub, Please Check..\n";
                die;
        }else {
                print "$sh qsub is Done!\n";
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
   Contact: songmm <songmm\@biomarker.com.cn>
      Date: 2016-08-11

     Usage:
            --datacfg      <FILE>  bwa mem generate sam file
            --detailcfg      <FILE>  reference geneome fasta file
	    --od       <DIR>     result dir
            --h                 help documents

   Example:
            perl $Script --datacfg data.cfg --detailcfg detail.cfg --od CircRNA_identify/find_circ

----------------------------------------------------------------------------------------------------------------------------
USAGE
	print $usage;
	exit;
}
