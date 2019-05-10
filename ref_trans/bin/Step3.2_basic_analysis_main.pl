#!/usr/bin/perl -w
use strict;
use Cwd qw(abs_path);
use Getopt::Long;
use Data::Dumper;
use FindBin qw($Bin $Script);
use File::Basename qw(basename dirname);
use newPerlBase;
my %config=%{readconf("$Bin/../config/db_file.cfg")}; 

my $BEGIN_TIME=time();
my $version="1.0.0";

# ------------------------------------------------------------------
# GetOptions
# ------------------------------------------------------------------

my ($bamdir, $detail_cfg, $odir, $step, $log);
GetOptions(
				"help|?"	=>\&USAGE,
				"bamdir:s"	=>\$bamdir,
				"cfg2:s"	=>\$detail_cfg,
				"od:s"		=>\$odir,
	                        "step:i"	=>\$step,
				) or &USAGE;
&USAGE unless ($bamdir and $detail_cfg and $odir) ;


system "mkdir -p $odir" unless (-d $odir);
$odir=abs_path($odir);
$bamdir=abs_path($bamdir);
$detail_cfg=abs_path($detail_cfg);

&log_current_time("$Script start...");
# ------------------------------------------------------------------
# preparation
# ------------------------------------------------------------------
my  %detail_cfg;
# read detail config
&detail_cfg_read($detail_cfg,\%detail_cfg);
my $index = $detail_cfg{Project_key};
my $project_id = $detail_cfg{Project_id};
my $cmd;
my $sh_dir = "$odir/work_sh";
my %step;
$step ||=  join ',',(1..6);
&steps_process($step,\%step);

# ------------------------------------------------------------------
# main pipline
# ------------------------------------------------------------------
#########################  Cufflinks
if ($step{1}) {
    $cmd = "perl $Bin/assembly_and_quantification.pl  -bam $bamdir  -cfg2 $detail_cfg -od $odir  ";
    print $cmd,"\n";
    system ($cmd);
    delete $step{1};    
}


open (SH,">$sh_dir/After_Cufflinks_analysis.sh") if ((keys %step)>0);
######################### Function Annotation & Integrate
if ($step{2}) {
    my $novel_unigene = "$odir/geneExpression/final_track/$index.newGene.longest_transcript.fa";
        mkdirOrDie ("$odir/PPI") unless (-d "$odir/PPI");
        $cmd="$config{python1} $Bin/annotation/PPI/PPI.py -fa $novel_unigene -cfg $detail_cfg -of $odir/PPI/PPI.txt -q $detail_cfg{Queue_type2} -p $detail_cfg{Project_key} \n ";
	mkdirOrDie("$odir/Anno_Integrate/New_Anno") unless (-d "$odir/Anno_Integrate/New_Anno");
	my $known_unigene = $detail_cfg{Known_unigene};
	my $all_unigene="$odir/geneExpression/final_track/All.longest_transcript.fa";
	`cat $known_unigene $novel_unigene >$all_unigene `;
	`cp $detail_cfg  $odir/Anno_Integrate/New_Anno/Gene_Func_Anno_Pipline.cfg`;
	my $anno_cfg="$odir/Anno_Integrate/New_Anno/Gene_Func_Anno_Pipline.cfg";
	`echo mRNA  $novel_unigene  >>$anno_cfg `;
	$cmd.= "perl $Bin/annotation/Gene_Anno/Gene_Func_Anno_Pipline.pl  --cfg  $anno_cfg --od $odir/Anno_Integrate/New_Anno  -queue  $detail_cfg{Queue_type2} ";
    
    # load specified database
    $cmd.= "--nr " if (exists $detail_cfg{anno}{nr});
    $cmd.= "--swissprot " if (exists $detail_cfg{anno}{Swissprot});
    $cmd.= "--kog " if (exists $detail_cfg{anno}{Kog});
    $cmd.= "--cog " if (exists $detail_cfg{anno}{Cog});
    $cmd.= "--kegg " if (exists $detail_cfg{anno}{Kegg});
    $cmd.= "--GO " if (exists $detail_cfg{anno}{nr});
    $cmd.= "--pfam " if (exists $detail_cfg{anno}{Pfam});
    $cmd.= "--eggNOG " if (exists $detail_cfg{anno}{eggNOG});
    my $judge=0;
    if(exists $detail_cfg{Known_anno}) {
        unless(-d  "$detail_cfg{Known_anno}/02.gene-annotation" || -d "$detail_cfg{Known_anno}/Result"){
            $judge=1;
        }
    }
    if (!exists $detail_cfg{Known_anno} || $judge==1){
        mkdirOrDie("$odir/Anno_Integrate/Known_anno") unless (-d "$odir/Anno_Integrate/Known_anno");
        `cp $detail_cfg  $odir/Anno_Integrate/Known_anno/Gene_Func_Anno_Pipline.cfg`;
        my $anno_cfg="$odir/Anno_Integrate/Known_anno/Gene_Func_Anno_Pipline.cfg";
        `echo mRNA  $known_unigene  >>$anno_cfg `;
        $cmd.="&& perl $Bin/annotation/Gene_Anno/Gene_Func_Anno_Pipline.pl  --cfg  $anno_cfg --od $odir/Anno_Integrate/Known_anno  -queue  $detail_cfg{Queue_type2}  ";
        $cmd.= "--nr " if (exists $detail_cfg{anno}{nr});
        $cmd.= "--swissprot " if (exists $detail_cfg{anno}{Swissprot});
        $cmd.= "--kog " if (exists $detail_cfg{anno}{Kog});
        $cmd.= "--cog " if (exists $detail_cfg{anno}{Cog});
        $cmd.= "--kegg " if (exists $detail_cfg{anno}{Kegg});
        $cmd.= "--GO " if (exists $detail_cfg{anno}{nr});
        $cmd.= "--pfam " if (exists $detail_cfg{anno}{Pfam});
        $cmd.= "--eggNOG " if (exists $detail_cfg{anno}{eggNOG});
        $cmd.= "--trembl " if (exists $detail_cfg{anno}{TrEMBL});
        $cmd.= "--nt " if (exists $detail_cfg{anno}{nt});
        $detail_cfg{Known_anno}="$odir/Anno_Integrate/Known_anno";
        if (exists $detail_cfg{Annotation_file}) {
            `echo "$detail_cfg{Ref_seq}\t$detail_cfg{Known_anno}\n" >>$detail_cfg{Annotation_file}`;
        }
        else{
             print "Error, Annotation_file should be offered\n Write info to Annotation_file\n$detail_cfg{Ref_seq}\t$detail_cfg{Known_anno}\n";
        }
    }
    
    my $known_anno=$detail_cfg{Known_anno};
    $cmd.= "&& perl $Bin/annotation/annotation_integrate/anno_integrated2.pl -i ${index}_Unigene -gene $known_unigene,$novel_unigene -anno $known_anno,$odir/Anno_Integrate/New_Anno -od $odir/Anno_Integrate/Allgene_Anno -cfg $detail_cfg ";
    if ($step{3}){
        print SH $cmd;
    }else{
        print SH $cmd,"\n";            
    }
}
######################### DEG Analysis
if ($step{3}) {
	my($new_gene_gff);
	$new_gene_gff = "$odir/geneExpression/final_track/$index.newGene_final.filtered.gff";
	`mkdir $odir/DEG_Analysis`	unless(-d "$odir/DEG_Analysis");
	`cp $odir/geneExpression/prepDE/Total_gene_fpkm.csv $odir/DEG_Analysis/All_gene_fpkm.list`;
	`cp $odir/geneExpression/prepDE/Total_gene_count.csv $odir/DEG_Analysis/All_gene_counts.list`;
	$cmd = "perl $Bin/DEG_Analysis/DEG_Analysis.pl -conf $detail_cfg -count $odir/DEG_Analysis/All_gene_counts.list -fpkm $odir/DEG_Analysis/All_gene_fpkm.list  -type FPKM -anno $odir/Anno_Integrate/Allgene_Anno/Result -ppi $odir/PPI/PPI.txt -gff $new_gene_gff -od $odir/DEG_Analysis ";
	$cmd.= "\nperl $Bin/DEG_Analysis/DEU/hisat_DEU_analysis.pl  -i $odir/DEG_Analysis/All_gene_counts.list -cfg $detail_cfg -od $odir/DEU_Analysis -hisat $bamdir -gtfdir $bamdir/../Ref_Genome";

        print SH "&&" if ($step{2});
	print "$cmd\n";
        print SH $cmd,"\n";
}

######################### Alternative Splicing Analysis
if ($step{4}) {
	$cmd = "perl $Bin/altsplice_analysis/v1.2/altsplice_analysis.pl -cufflinks  $odir/StringTie  -Ref_Genome $bamdir/../Ref_Genome -od $odir/Alitsplice_Analysis  -queue $detail_cfg{Queue_type2} ";
	print SH $cmd,"\n";
}

########################Diff_AS_analysis
if ($step{5}) {
	$cmd = "perl $Bin/differential_AS_analysis/diff_AS_analysis.pl --cfg $detail_cfg --bamdir $bamdir --od $odir/differential_AS_analysis ";
	#$cmd.= "&& perl $Bin/differential_AS_analysis/diff_AS_res_deal.pl --cfg $detail_cfg --od $odir/differential_AS_analysis";
	print SH $cmd,"\n";
}

######################### Genic Structure Optimize
if ($step{6}) {
	$cmd = "perl $Bin/gene_optimize/gene_structure_optimize.pl -gtf $odir/Compare/gffcmp.annotated.gtf -gff $detail_cfg{Ref_ann} -od $odir/Gene_Structure_Optimize -index $index ";
	print SH $cmd,"\n";   
}
close (SH);
&qsubOrDie("$sh_dir/After_Cufflinks_analysis.sh",  $detail_cfg{Queue_type2} ,  6, "20G");


#######################################################################################
my $elapse_time = (time()-$BEGIN_TIME)."s";
&log_current_time("$Script done. Total elapsed time: $elapse_time.");
#######################################################################################
# ------------------------------------------------------------------
# sub function
# ------------------------------------------------------------------
#############################################################################################################
sub data_cfg_read {
    my ($cfg_file, $data_cfg) = @_;
    my $sample_id;

    open (CFG, $cfg_file) or die "$!: $cfg_file\n";
    while (<CFG>) {
        chomp;
        next if (/^\s+$/ or /^#/);

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
    &log_current_time("data config done.");
}
#################################################################
sub steps_process {
    print "\n";
    &log_current_time("step check:");
    my ($step_str,  $step_hash) = @_;
    my %_step_ = (
        '1' => 'Cufflinks_Expression',
        '2' => 'Annoation',
        '3' => 'Difference_Enrichment',
        '4' => 'Alternative_splicing',
	'5' => 'differential_AS_analysis',
	'6' => 'Gene_Optimization',
    );
        for my $s (split /,/,$step_str) {
            if ($s =~/^[1-6]$/) {
                $step_hash->{$s} = $_step_{$s};
            } else {
                print "ERROR: illegal steps specified by --step.\n";
                die;
            }
    }

    print "steps_to_run: ".(join ", ",(map {sprintf "$_.$step_hash->{$_}"} sort keys %{$step_hash})).".\n";
    &log_current_time("step check done.\n");
}

#############################################################################################################
sub detail_cfg_read {
    &log_current_time("detail config check:");
    my ($cfg_file, $detail_cfg) = @_;

    open (CFG,$cfg_file ) or die "$!: $cfg_file\n";
    while (<CFG>) {
        chomp;
        next if (/^\s+$/ or /^#/);
        my ($key, $value) = (split /\s+/)[0,1];
	next unless(defined $key);
        $detail_cfg->{$key} = $value;
        if ($key eq 'Project_name' or $key eq 'Customer_info' or $key eq 'Project_id' or $key eq 'Project_key') {
            $detail_cfg->{$key} = $value;
        }
        if ($key eq 'Known_unigene' or $key eq 'Known_pep') {
            die "$key: $value is not exist!\n" unless (-e $value);
            $detail_cfg->{$key} = $value;
        }

        if ($key eq 'Known_anno') {
            $detail_cfg->{$key} = $value;
        }
        if ($key eq 'Ref_seq' or $key eq 'Ref_ann') {
            die "$key: $value is not exist!\n" unless (-e $value);
            $detail_cfg->{$key} = $value;
        }
        if ($key =~/^SNP_/) {
            $detail_cfg->{$key} = $value;
        }
        if ($key eq 'nr' or $key eq 'Swissprot' or $key eq 'Kegg' or $key eq 'Pfam' or $key eq 'Cog' or $key eq 'Kog'  or $key eq 'eggNOG'  ) {
            die "$key: $value is not exist!\n" unless (-e $value);
            $detail_cfg->{anno}{$key} = $value;
        }
        if ($key eq 'Queue_type1') {
            $detail_cfg->{$key} = $value;
        }
        if ($key eq 'Queue_type2') {
            $detail_cfg->{$key} = $value;
        }
        print "$key: $value\n" if (exists $detail_cfg->{$key});
    }
    close CFG;
    die "Must choose Queue_type1 and Queue_type2 !\n" unless (exists $detail_cfg->{Queue_type1} or exists $detail_cfg->{Queue_type2});
    &log_current_time("detail config check done.");
}

#############################################################################################################
sub step_cmd_process {
    my ($cmd, $sh_name, $sh_dir) = @_;
    my $sh_file = "$sh_dir/$sh_name";
    my $log_file = "$sh_file.log";
    my $flag = 0;
    my $start_time = time();
    &log_current_time("$sh_name start...");
#    &log_current_time("CMD: $cmd");

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

    $flag = system("sh $sh_file > $log_file");
    if ($flag != 0){
        log_current_time("Error: command failed: $cmd");
        exit(1);
    } else {
        my $escaped_time = (time()-$start_time)."s";
        &log_current_time("$sh_name done, escaped time: $escaped_time.");
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
----------------------------------------------------------------------------------------------------
   Program: $Script
   Version: $version
   Contact: Simon Young <yangxh\@biomarker.com.cn> 
      Date: 2014-11-13

     Usage:
            --bamdir	<DIR>   bam dir
            --cfg2      <FILE>  detail config, analysis parameters
            --od         <DIR>  analysis output directory

            --step      <INT>   steps to run (split by comma) [1] default : 1,2,3,4,5,6
                          1     Cufflinks & Gene Expression & Subsequence
                          2     Gene Annoation
                          3     Difference & Enrichment
                          4     Alternative splicing
			  5	Differential_AS_analysis
                          6     Gene Optimization
            --h                 help documents

   Example:
            perl $Script  --bamdir  Tophat/Hisat --cfg2 detail.cfg --od Structure_and_Expression 

----------------------------------------------------------------------------------------------------
USAGE
	print $usage;
	exit;
}
