#!/usr/bin/perl -w
use strict;
use warnings;
use Getopt::Long;
use Data::Dumper;
use FindBin qw($Bin $Script);
use File::Basename qw(basename dirname);
my $BEGIN_TIME=time();
use newPerlBase;
my $Title="LncRNA";  
my $version="v2.0.3"; 
my @Original_ARGV=@ARGV;
#######################################################################################
# ------------------------------------------------------------------
# GetOptions, Check format, Output info.
# ------------------------------------------------------------------
my ( $gtf,$gff , $od, $step, $key, $type, $db, $log,$oneStepOnly,$cfg,$prep);
GetOptions(
				"help|?" =>\&USAGE,
				"gtf:s"  =>\$gtf,
				"gff:s"  =>\$gff,
				"prep:s" =>\$prep,
				"od:s"   =>\$od,
				"cfg:s"    =>\$cfg,
				"type:s"=>\$type,
				"db:s"=>\$db,
				"step:s"=>\$step,
				"oneStepOnly"=>\$oneStepOnly,
				) or &USAGE;
&USAGE unless ($gtf and $gff and $cfg and $od and $prep) ;
################################

#########################
my $notename = `hostname`; chomp $notename;
$cfg = &ABSOLUTE_DIR($cfg);
$gtf = &ABSOLUTE_DIR($gtf);
$gff = &ABSOLUTE_DIR($gff);

&MKDIR($od);
$od  = &ABSOLUTE_DIR($od);     &MKDIR("$od/work_sh");
&MKDIR("$od/lnc_predict");
&show_log2("step_3: lncRNA analysis start.");

$step = $step || 1;
my %para;
open (IN,"$cfg") || die "$!\n";
while (<IN>) {
        chomp;
        s/\r$//;s/^\s+//;s/\s+$//;
        next if (/^\#/ || /^$/);

        my @tmp=split /\s+/,$_;
        $para{$tmp[0]}=$tmp[1];
}
if (exists $para{'db'}){
	$db=$para{'db'};
}elsif (exists $para{'nr'}){
	my $nr_name=basename($para{'nr'});
	$nr_name=~s/nr_//;
	$db=$nr_name;
}else{$db="MAM";}

if (exists $para{'type'}){
	$type=$para{'type'};
}elsif (defined $db and $db =~/INV|MAM|PRI|ROD|SYN|VRT|ENV|UNA/){
	$type="ve";
}elsif(defined $db and $db =~/BCT|PLN/){
	$type="pl";
}else{$type="ve";}
$para{'Queue_type'}="general.q" unless (exists $para{'Queue_type'});
$para{'Memory'} = "15G" unless (exists $para{'Memory'});
$para{'CPU'} = "50" unless (exists $para{'CPU'});

my $startTime = GetTime();
my $user = `whoami`;  chomp $user;
my $workDir = `pwd`;  chomp $workDir;
my $task = "perl $Bin/$Script ".join("\t", @Original_ARGV);

open ($log,">", "$od/lncRNA_identfy.".time().".log") or die $!;

print $log "######################################\n";
print $log "$startTime\n";
print $log "user:	$user\n";
print $log "work directory:	$workDir\n";
print $log "job:	$task\n";
print $log "######################################\n\n";

print $log "start from step: $step\n\n";


#==================================================================
&MKDIR("$od/lnc_predict/Lnc_filter");
&MKDIR("$od/lnc_predict/code_filter");
my $cmd;
#my $lnc_predict="/share/nas2/genome/bmksoft/pipeline/LncRNA_pipeline/v3.1.6/bin/basic_analysis/lncRNA_Analysis/lnc_predict/Lnc_predict.v6.pl";
my $lnc_predict = "$Bin/lnc_predict/Lnc_predict.v6.pl";

my $lnc_fa="$od/lnc_predict/Candicate_lncRNA.fa";
if ($step==1) {
	print STDOUT "===  lncRNA prediction  ===\n";
	print STDOUT "lncRNA prediction shell file: $od/work_sh/lncRNA prediction1.sh\n";
	print $log "=== lncRNA prediction ===\n";
	print $log "lncRNA prediction shell file: $od/work_sh/lncRNA prediction1.sh\n";
	my $fpkm="$prep/Candicate_lncRNA_fpkm.list";

	open OUT,">$od/work_sh/lncRNA_prediction.sh" || die;
	my $cmd= "perl $Bin/lnc_predict/get_Candicate_lncRNA_fa.pl -gtf $gtf -fa $para{'Ref_seq'} -fpkm $fpkm -o $lnc_fa";
	$cmd  .= "&& perl $lnc_predict -cfg $cfg -fa $lnc_fa -type $type -db $db -od $od/lnc_predict/code_filter -key lnc_code_filter";
	$cmd  .= " -pfam $para{'DATA'}	"if(exists $para{'DATA'});
	print OUT "$cmd \n";
        close OUT;
        &Cut_shell_qsub("$od/work_sh/lncRNA_prediction.sh",$para{'CPU'},$para{'Memory'},$para{'Queue_type'});
        &Check_qsub_error("$od/work_sh/lncRNA_prediction.sh");

        $step++ unless ($oneStepOnly) ;
        print STDOUT "\n";
        print $log "\n";
}
if ($step==2) {
        print STDOUT "===  lncRNA prediction and statistics  ===\n";
        print $log "=== lncRNA statistics ===\n";
        open OUT,">$od/work_sh/lncRNA_statistics.sh" || die;
#	my $cmd = "perl $Bin/lnc_predict/lnc_predict_veen.pl -cpc $od/lnc_predict/code_filter/CPC/lnc_code_filter.result.txt -cnci $od/lnc_predict/code_filter/CNCI/CNCI.index -pfam $od/lnc_predict/code_filter/Pfam/Pfam_result.txt -cpat $od/lnc_predict/code_filter/CPAT/cpat.txt -od $od/lnc_predict/code_filter/ -fa $od/lnc_predict/Candicate_lncRNA.fa && "; 
	my $cmd .= "grep -v 'NA' $od/lnc_predict/code_filter/list.txt |cut -f 1 >$od/lnc_predict/lnc_filter_id.list  && ";		###get overlap lncRNA
	$cmd.= "perl $Bin/lnc_predict/abstract_gtf_seq_by_transid.pl -i $od/lnc_predict/lnc_filter_id.list -gtf $gtf -o $od/lnc_predict/filter_final.gtf && perl $Bin/lnc_predict/gtf_to_gff.pl -i $od/lnc_predict/filter_final.gtf -o $od/lnc_predict/filter_final.gff && perl $Bin/lnc_predict/get_final_lnc_fa.pl $lnc_fa $od/lnc_predict/filter_final.gff $od/lnc_predict/lnc_filter_final.fa &&";	###get novel lncRNA gtf/gff/fa

	if(exists $para{Lnc_ann}){
		$cmd.=" cat $od/lnc_predict/filter_final.gtf $od/../Hisat_Stringtie/LncPredict/Known_lncRNA.gtf >$od/lnc_predict/All_filter_final.gtf && cat $od/../Hisat_Stringtie/LncPredict/Known_lncRNA.gff $od/lnc_predict/filter_final.gff >$od/lnc_predict/All_filter_final.gff && ";
	}else{
		$cmd.=" cp $od/lnc_predict/filter_final.gtf $od/lnc_predict/All_filter_final.gtf && cp $od/lnc_predict/filter_final.gff $od/lnc_predict/All_filter_final.gff &&";
	}
	$cmd .=" perl $Bin/lnc_predict/get_lncRNA_fa.pl $od/lnc_predict/All_filter_final.gtf $para{'Ref_seq'}  $od/lnc_predict/All_filter_final.fa && ";

	$cmd.="perl $Bin/lnc_predict/lnc.stat.class_code.pl -i $od/lnc_predict/filter_final.gtf -o $od/lnc_predict/Lnc_filter  &&";
	$cmd.="perl $Bin/lnc_predict/lncRNA_class.pl -gtf $od/lnc_predict/filter_final.gtf -out $od/lnc_predict/lncRNA_class_id.list &&";
	$cmd.=" perl $Bin/lnc_predict/lncRNA_draw_circos.pl -genome $od/../Hisat_Stringtie/Ref_Genome/genome_size.txt -gff $od/lnc_predict/filter_final.gff -gtf  $od/lnc_predict/filter_final.gtf -od $od/lnc_predict/circos  &&  ";
	$cmd.= " perl $Bin/lnc_predict/get_final_lncRNA_exp.pl $prep $od/lnc_predict/lnc_filter_id.list $od/../Hisat_Stringtie/prepDE";
	print OUT "$cmd \n";
	close OUT;
	&Cut_shell_qsub("$od/work_sh/lncRNA_statistics.sh",$para{'CPU'},$para{'Memory'},$para{'Queue_type'});
	&Check_qsub_error("$od/work_sh/lncRNA_statistics.sh");
	
	$step++ unless ($oneStepOnly) ;
	print STDOUT "\n";
	print $log "\n";
}

if ($step==3) {
	print STDOUT "=== lncRNA targetgene prediction   ===\n";
	print STDOUT "lncRNA targetgene prediction shell file: $od/work_sh/Lnc_target_predict.sh\n";
	print $log "=== lncRNA_targetgene_predict ===\n";
	print $log "lncRNA_targetgene_predict shell file: $od/work_sh/Lnc_target_predict.sh\n";

	my $sample_num=&get_Sample("$prep/All_lncRNA_fpkm.list") -1 ;
	open OUT,">$od/work_sh/Lnc_target_predict.sh" || die;
	&MKDIR("$od/Lnc_target_predict");

	#########################Cis target######################
	my $cmd;
	if(exists $para{'Lnc_ann'}){
		$cmd .="cat $para{'Lnc_ann'} $od/lnc_predict/filter_final.gff >$od/Lnc_target_predict/All_lncRNA.gff && ";
	}else{	$cmd .="cp $od/lnc_predict/filter_final.gff >$od/Lnc_target_predict/All_lncRNA.gff && ";}
	
	$cmd .="perl $Bin/predict_target/Cis_target.pl -lncgff $od/Lnc_target_predict/All_lncRNA.gff -genegff $gff -output $od/Lnc_target_predict/Cis_target_gene.xls -cis $para{'Cis_dist'}\n";	

	########################Trans target#####################
	if($sample_num>=5){
		$para{'Trans_p'}	||="0.01";
		$para{'Trans_cor'}	||="0.9";
		$para{'Trans_method'}	||="pearson";
		$cmd .="perl $Bin/predict_target/co_expression.pl -lncRNA $prep/All_lncRNA_fpkm.list -mRNA $prep/All_gene_fpkm.list -od $od/Lnc_target_predict/Trans -cor $para{Trans_cor} -pvalue $para{Trans_p} -method $para{Trans_method} -queue $para{Queue_type}\n";
	}
	#########################WGCNA###########################
	if($sample_num>=15){
		$cmd .="Rscript $Bin/predict_target/WGCNA_v1.2.R --indir $prep/ --outdir $od/Lnc_target_predict/Trans_target/ --meanFPKM 5 -f 0.2 --samplenum $sample_num\n";
	}
	print OUT "$cmd";
	close OUT;
	&Cut_shell_qsub("$od/work_sh/Lnc_target_predict.sh",$para{'CPU'},$para{'Memory'},$para{'Queue_type'});
	&Check_qsub_error("$od/work_sh/Lnc_target_predict.sh");
	$step++ unless ($oneStepOnly) ;
	print STDOUT "\n";
	print $log "\n";
}
if($step==4){
	if(exists $para{medical}){
		open OUT,">$od/work_sh/Lnc_disease.sh" || die;
		my $cmd="perl $Bin/Lnc_Disease/lnc_disease_anno.pl $od/lnc_predict/All_filter_final.fa $od/Lnc_disease $para{medical} ";
		print OUT "$cmd\n";
		close(OUT);
		&Cut_shell_qsub("$od/work_sh/Lnc_disease.sh",$para{'CPU'},$para{'Memory'},$para{'Queue_type'});
		&Check_qsub_error("$od/work_sh/Lnc_disease.sh");
	}

}
	

&show_log2("step_3: lncRNA analysis finished.");
#        $cmd .=" && perl $Bin/combine/random/random_galk_anno.pl -i $od/Combine/ceRNA/ceRNA_pair_adjust_p_Sig_diff.txt -type $od/Combine/ceRNA/ceRNA_pair_adjust_p_Sig_diff.node -ratio 0.05 -num 5 -anno $od/Anno_Integrate/Allgene_Anno/Result/ -od $od/Combine/random\n";
#######################################################################################
print STDOUT "\nDone. Total elapsed time : ",time()-$BEGIN_TIME,"s\n";
print $log  "\nDone. Total elapsed time : ",time()-$BEGIN_TIME,"s\n";
#######################################################################################
close($log);
#==================================================================
# subs 
#==================================================================
sub get_Sample{
	my $file=shift;
	my $head=`head -n 1 $file`;
	chomp($head);
	my @tmp=split(/\s+/,$head);
	return scalar(@tmp);
}
sub Cut_shell_qsub {#Cut shell for qsub 1000 line one file
	# &Cut_shell_qsub($shell,$cpu,$vf,$queue);
	my $shell = shift;
	my $cpu = shift;
	my $vf = shift;
	my $queue = shift;

	my $status;
	my $line = `less -S $shell |wc -l `;
	chomp $line;
	if ($line<=1000) {
		if ($notename=~/login\-0\-4/) {
			$status =system "sh /share/nas2/genome/bmksoft/tool/qsub_sge_plus/v1.0/qsub_sge.plus.sh $shell --queue $queue --maxproc $cpu --resource vf=$vf --reqsub --independent";
		}else{
			$status =system "sh /share/nas2/genome/bmksoft/tool/qsub_sge_plus/v1.0/qsub_sge.plus.sh $shell --queue $queue --maxproc $cpu --resource vf=$vf --reqsub --independent";
		}
	}
	if ($line>1000) {
		my @div=glob "$shell.div*";
		foreach (@div) {
			if (-e $_) {
				system "rm $_";
			}
		}
		@div=();
		my $div_index=1;
		my $line_num=1;
		open IN,"$shell" || die;
		while (<IN>) {
			chomp;
			open OUT,">>$shell.div.$div_index.sh" || die;
			if ($line_num<1000) {
				print OUT "$_\n";
				$line_num++;
			}
			else {
				print OUT "$_\n";
				$div_index++;
				$line_num=1;
				close OUT;
			}
		}
		if ($line_num!=1) {
			close OUT;
		}
		@div=glob "$shell.div*";
		foreach my $div_file (@div) {
			if ($notename=~/login\-0\-4/) {
				$status = system "sh /share/nas2/genome/bmksoft/tool/qsub_sge_plus/v1.0/qsub_sge.plus.sh $div_file --queue $queue --maxproc $cpu --resource vf=$vf --reqsub --independent";
			}else{
				$status = system "sh /share/nas2/genome/bmksoft/tool/qsub_sge_plus/v1.0/qsub_sge.plus.sh $div_file --queue $queue --maxproc $cpu --resource vf=$vf --reqsub --independent";
			}
		}
	}
	$status ==0 or die "qsub die ,please check\n";
}

sub Check_qsub_error {#
	# Check The qsub process if error happend 
	my $sh=shift;
	my @Check_file=glob "$sh*.qsub/*.Check";
	my @sh_file=glob "$sh*.qsub/*.sh";
	if (@sh_file==0){
		print "qsub launch Error,Please Check..\n";die;
	}
	if ($#sh_file!=$#Check_file) {
		print "Their Some Error Happend in $sh qsub, Please Check..\n";
		die;
	}
	else {
		print "$sh qsub is Done!\n";
	}
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
#	rmdir($dir) if(-d $dir);
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
sub step_cmd_process {
    my ($cmd, $sh_name, $sh_dir) = @_;
    my $sh_file = "$sh_dir/$sh_name";
    my $log_file = "$sh_file.log";
    my $flag = 0;
    my $start_time = time();
    &log_current_time("$sh_name start...");

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
#################
sub sub_format_datetime {#Time calculation subroutine
    my($sec, $min, $hour, $day, $mon, $year, $wday, $yday, $isdst) = @_;
    $wday = $yday = $isdst = 0;
    sprintf("%4d-%02d-%02d %02d:%02d:%02d", $year+1900, $mon+1, $day, $hour, $min, $sec);
}

# &show_log2("cmd")
sub show_log2()
{
    open LOG, ">>$od/../../../bioinfor_pipeline.log";
	my ($txt) = @_ ;
    my $time = time();
    my $Time = &sub_format_datetime(localtime($time));
    print LOG "$Time:\t$txt\n" ;
    return ($time) ;
	close LOG;
}
####################


sub USAGE {#
	my $usage=<<"USAGE";

Description:
	This program is a Procedure deal with lncRNA prediction 


Usage:
	-cfg	detail.cfg ,	must be given;
	-gtf	Merged_filter.gtf (the compare gtf file with Lnc_ann), must be given;
	-gff	All_gene.gff, must be given
	-od 	output dir , ./Lnc_filter must be given;
	-prep	prep RNA quantify dir, must be given;
	-step		
		1:lncRNA prediction ;
		2:lncRNA statistics;
		3:lncRNA targetgene prediction;
		4:lncRNA disease
	-oneStepOnly
	-type   DataBase For CNCI predict:"ve" for vertebrate species, "pl" for plat species;
        -db     DataBase For CPC predict:"PLN,BCT,INV,MAM,PHG,PRI,ROD,SYN,UNA,VRL,VRT,ENV";
	
USAGE
	print $usage;
	exit;
}
