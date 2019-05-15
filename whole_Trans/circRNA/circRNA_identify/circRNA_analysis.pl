#!/usr/bin/perl -w
use strict;
use Cwd qw(abs_path);
use Getopt::Long;
use Data::Dumper;
use FindBin qw($Bin $Script);
use File::Basename qw(basename dirname);
use newPerlBase;
my $BEGIN_TIME=time();
my %config=%{readconf("$Bin/../project.cfg")};
my $Title=$config{Title};											
my $version=$config{version};
# ------------------------------------------------------------------
# GetOptions
# ------------------------------------------------------------------
my ($data_cfg, $detail_cfg, $odir, $idir,$gtf,$gff,$step);

GetOptions(
    "data_cfg:s"	=>\$data_cfg,
    "detail_cfg:s"	=>\$detail_cfg,
    "gtf:s"	=>\$gtf,
    "gff:s"	=>\$gff,
    "od:s"	=>\$odir,
    "id:s"	=>\$idir,
    "step:s"	=>\$step,	
    "help|h"	=>\&USAGE,
    ) or &USAGE;
&USAGE unless ($data_cfg and $detail_cfg and $odir);

system "mkdir -p $odir" unless (-d $odir);
$odir=abs_path($odir);
$idir="$odir/../circRNA_Bwa/Bwa";
if(!-d $idir){
	print "$idir do not exists! Please give the correct dir!\n";
	die;
}
$data_cfg=abs_path($data_cfg);
$detail_cfg=abs_path($detail_cfg);

&log_current_time("$Script start...");

# ------------------------------------------------------------------
# read configs and steps
# ------------------------------------------------------------------

my (%data_cfg, %detail_cfg);
&data_cfg_read($data_cfg,\%data_cfg);
&detail_cfg_read($detail_cfg,\%detail_cfg);

###known circRNA annotation database
my $Circ2diasease = "/share/nas1/niepy/Database/Circ2Disease/Circ2Disease_Association.txt";
my $circBank = "/share/nas1/niepy/Database/circBank/circBank_annotation.txt";

# ------------------------------------------------------------------
# main pipline
# ------------------------------------------------------------------

my $sh_dir = "$odir/work_sh";
mkdir $sh_dir unless (-d $sh_dir);
my $cmd;
my $Rscript = $config{Rscript};
my $GFFREAD_BIN     = $config{gffread};
$step ||= 1;
my $circ;
my @samples = (sort keys %{$data_cfg{rawdata}});
my $sample = join(",",@samples);
######################### circRNA identify
if ($step==1) {
	system "mkdir -p $odir/circRNA_identify" unless (-d "$odir/circRNA_identify");
	if(!defined $gtf){
		$gff||=$detail_cfg{Ref_ann};
		$gtf = "$odir/circRNA_identify/".basename $gff;
		$gtf=~s/gff\d*/gtf/;
		`$GFFREAD_BIN $detail_cfg{Ref_ann} -T -o $gtf`;
	}
	open SH,">$sh_dir/circRNA_identify.sh";
	$cmd = "perl $Bin/bin/circRNA_identify.v3.pl -samFile  $idir -ref $detail_cfg{Ref_seq} -gff $gtf -od $odir/circRNA_identify";
	$cmd.=" -re $detail_cfg{RE} " if (exists $detail_cfg{RE});
	$cmd.="\n";
	$cmd .="perl $Bin/bin/find_circ.pl -datacfg $data_cfg -detailcfg $detail_cfg -od $odir/circRNA_identify/find_circ\n" if(exists $detail_cfg{find_circ} && $detail_cfg{find_circ}==1);
	$cmd .="perl $Bin/bin/CIRCexplorer2.pl -samFile $idir -ref $detail_cfg{Ref_seq} -gff $gtf -od $odir/circRNA_identify/CIRCexplorer\n" if(exists $detail_cfg{CIRCexplorer} && $detail_cfg{CIRCexplorer}==1);
	print SH "$cmd";
	close SH;
	&run_or_die("$config{qsub} --resource vf=50G --queue $detail_cfg{Queue_type} --independent $sh_dir/circRNA_identify.sh");
	&run_or_die("perl $Bin/bin/intersect.pl --indir $odir/circRNA_identify/ --odir $odir/circRNA_identify");
	$step = 2;
}

my @circs=glob "$odir/circRNA_identify/*.intersect";

if(scalar(@circs)!=scalar(@samples)){
	print "Something must be wrong!\n";
	print scalar(@circs)," samples have intersect file, while we have ",scalar(@samples)," samples in total!\n";
	die;
}
$circ = join ",",(glob "$odir/circRNA_identify/*.intersect");

######################### circRNA expression analysis
if ($step==2) {
 	mkdir("$odir/expression/") if(!-d "$odir/expression/");
	$cmd = "$config{Rscript} $Bin/bin/R/combination.r --infile $circ --sample $sample --outfile $odir/expression/ ";
	&step_cmd_process($cmd,"expression",$sh_dir);
	$step = 3;
}
######################### circRNA sequence
if ($step ==3) {
	$gff||=$detail_cfg{Ref_ann};
	$cmd ="perl $Bin/bin/exactcircRNA_CIRI.pl -gff $gff -c $odir/expression/All_gene_counts_detail.xls_tmp -ref $detail_cfg{Ref_seq} -od $odir/circRNA_identify/ && mv $odir/circRNA_identify/All_gene_counts_detail.xls $odir/circRNA_identify/All_gene_counts.list $odir/expression/ && mv $odir/circRNA_identify/*.fa $odir/circRNA.fa && perl $Bin/bin/exact_circRNA_gene.pl $odir/expression/All_gene_counts_detail.xls $odir/circRNA_identify/circRNA_gene.list && perl $Bin/bin/circRNA_vennFile.pl -i $odir/expression/All_gene_counts.list -o $odir/expression ";
	&step_cmd_process($cmd,"circRNA_sequence",$sh_dir);
	my @venn = glob("$odir/expression/*.venn.txt");
        my $venn="";
        my $name;
        foreach(@venn)
        {
                my $tmp = basename($_);
                $venn.="--lst $_ ";
                $tmp =~s/.venn.txt//;
                $name .="--lab $tmp ";
        }
        $cmd="perl $Bin/bin/draw_venn_diagram.pl $venn $name --pre sample --od $odir/expression";
        `$cmd && touch $odir/expression/venn.ok` if(@venn>=2 && @venn<=5 && !-e "$odir/expression/venn.ok");
	$step = 4;
}
######################### circRNA statistics
if ($step ==4) {
	mkdir("$odir/statistics") if(!-d "$odir/statistics");
	$cmd = "perl $Bin/bin/sep_circRNA_result.pl --i $circ --s $sample --o $odir/statistics/circRNA.type.stat  && $Rscript $Bin/bin/R/pie.r $odir/statistics/circRNA.type.stat $odir/statistics/ $sample && $Rscript $Bin/bin/R/simpleBar.r --infile $odir/expression/All_gene_counts.list --outfile $odir/statistics/ --x.col 1 --x.lab \"chromosome\" --y.lab \"the number of junction reads\" --sample $sample --no.grid && $Rscript $Bin/bin/R/statistics.r $odir/expression/All_gene_counts.list $odir/statistics/statistic.xls ";
	my $length_files = (glob("$odir/circRNA_identify/*_length.stat"))[0];
	my $prefix = "All_circRNA";
	$cmd .=" && perl $Bin/bin/lnc_stat.pl -i $length_files -o $odir/statistics/ -k $prefix";
    	&step_cmd_process($cmd,"circRNA_statistics",$sh_dir);
	if(exists $detail_cfg{circBase})
	{	
		`mkdir $odir/statistics/Known`	unless(-d "$odir/statistics/Known");
		&run_or_die("perl $Bin/bin/newCircRNA.pl -ref $detail_cfg{circBase} -c $odir/circRNA.fa -od $odir/statistics/Known");
		&run_or_die("$Rscript $Bin/bin/R/pie.r $odir/statistics/Known/forpie.list $odir/statistics/Known/circ_pred");
		if(exists $detail_cfg{medical} && $detail_cfg{medical}=~/GRCh37|GRCh38/i){ 
			&run_or_die("perl $Bin/bin/known_circRNA_annotation.pl -b $odir/statistics/Known/ratio.out -cd $Circ2diasease -cb $circBank -od $odir/statistics/Known/");
		}
	}
	$step = 5;
}
######################### circRNA alternative splicing analysis
if ($step ==5) {
	mkdir("$odir/overlap_alitisplice/") if(!-d "$odir/overlap_alitisplice/");
	$cmd = "perl $Bin/bin/overlap_alitisplice.pl -i $odir/expression/All_gene_counts_detail.xls -o $odir/overlap_alitisplice/overlap_alitisplice.xls";
	&step_cmd_process($cmd,"alternative_splicing",$sh_dir);
	$step = 6;
}
######################### circRNA Rename with gene_id
if ($step ==6) {
	mkdir("$odir/new_name/") if(!-d "$odir/new_name/");
	my $anno = "$detail_cfg{Known_anno}/Result/Integrated_Function.annotation.xls";
	if(!exists $detail_cfg{Symbol}){
		my $dir=dirname $detail_cfg{Ref_seq};
		$detail_cfg{Symbol}="$dir/id_name.list"	if(-e "$dir/id_name.list");
	}
	$cmd = "perl $Bin/bin/circRNARename.pl -i $odir/expression/All_gene_counts_detail.xls -o $odir/new_name/circRNA_newname.xls -gff $gff -anno $anno ";
	$cmd .=" -label $detail_cfg{miRBase} "	if(exists $detail_cfg{miRBase});	
	$cmd .=" -symbol $detail_cfg{Symbol} "	if(exists $detail_cfg{Symbol});
	print "$cmd\n";
	&step_cmd_process($cmd,"circRNA_Rename",$sh_dir);
	$step++;
}

if($step ==7){
	###count to exp , circos for count
	my $normalization;
	if(exists $detail_cfg{normalization} && $detail_cfg{normalization} eq "SRPBM"){
		$normalization = 'SRPBM';
		my $mapstat = "$odir/../circRNA_Bwa/Map_Stat/All.mappedStat.xls";
		my $countlen = "$odir/../circRNA_Bwa/Map_Stat/All.readLength.xls";
		$cmd = "perl $Bin/exp/count_to_expression.pl -i $odir/expression/All_gene_counts.list -o $odir/expression -l $normalization -m $mapstat -len $countlen ";
	}elsif(exists $detail_cfg{normalization} && $detail_cfg{normalization} eq "TPM") {
		$normalization = 'TPM';
		$cmd ="perl $Bin/exp/count_to_expression.pl -i $odir/expression/All_gene_counts.list -o $odir/expression -l $normalization ";
	}elsif(exists $detail_cfg{normalization} && $detail_cfg{normalization} eq "FPKM"){
		$normalization = 'FPKM';
		$cmd ="perl $Bin/exp/count_to_expression.pl -i $odir/expression/All_gene_counts.list -o $odir/expression -l $normalization ";
	}
	$cmd.="&& perl $Bin/circos/draw_circos_circRNA.pl -c $odir/expression/All_gene_counts.list -e $odir/expression/All_gene_expression.list -od $odir/circos";
	if(exists $detail_cfg{medical}){
		$cmd.=" -genome $Bin/../../tools/circos/$detail_cfg{medical}.txt -s chr";
	}else{
		$cmd.=" -genome $detail_cfg{Ref_seq} ";
		$cmd.=" -s $detail_cfg{label} "		if(exists $detail_cfg{label});
	}
	&run_or_die($cmd);
}
#######################################################################################
my $elapse_time = (time()-$BEGIN_TIME)."s";
&log_current_time("$Script done. Total elapsed time: $elapse_time.");
#######################################################################################
# ------------------------------------------------------------------
# sub function
# ------------------------------------------------------------------
#############################################################################################################
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
}
############################
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
        next if (/^\s*$/ or /^#/);
        my ($key, $value) = (split /\s+/)[0,1];
		$detail_cfg->{$key} = $value;
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
        #print "$key: $value\n" if (exists $detail_cfg->{$key});
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
	if(-e "$sh_file.finish")
	{
		return;
	}
	else
	{
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
			`touch $sh_file.finish`;
		}
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
            --data_cfg      <FILE>  data config, rawdata path
            --detail_cfg      <FILE>  detail config, analysis parameters & refseq path
	    --gff	<FILE>	If given will not use Ref_ann in detail_cfg
	    --gtf	<FILE>	If given will not use Ref_ann in detail_cfg
            --od        <DIR>   analysis output directory
            --h                 help documents
	    --step	<int>	begin from which step , default 1

   Example:
            perl $Script --data_cfg data.cfg --detail_cfg detail.cfg --od CircRNA_identify/

----------------------------------------------------------------------------------------------------------------------------
USAGE
	print $usage;
	exit;
}
