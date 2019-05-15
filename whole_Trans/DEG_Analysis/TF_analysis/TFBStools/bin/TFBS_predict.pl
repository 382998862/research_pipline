#!/usr/bin/perl -w

use strict;
use Cwd qw(abs_path);
use Getopt::Long;
use Getopt::Std;
use Data::Dumper;
use FindBin qw($Bin $Script);
use File::Basename qw(basename dirname);
use threads;
use newPerlBase;
# ------------------------------------------------------------------
# # GetOptions
# # ------------------------------------------------------------------
my ($gff, $genome, $odir, $detail_cfg, $newgene_gff, $score);
GetOptions(
	"gff2:s"        =>\$newgene_gff,
	"cfg:s"		=> \$detail_cfg,
	"od:s"		=>\$odir,
	"score:s"	=>\$score,
	"help|h"	=>\&USAGE,
) or &USAGE;
&USAGE unless ($odir and $detail_cfg);

########read files
#
`mkdir -p $odir`  unless (-d $odir);
$odir = abs_path($odir);
$detail_cfg = abs_path($detail_cfg);
$score||="90%";

open (CFG, "$detail_cfg")||die $!;
my ($spe_id, $medical, $key);
while(<CFG>){
        chomp;next if(/^#/);
        if(/^spe_id/) {
                $spe_id=(split/\s+/,$_)[1];
        }
        if(/^Ref_seq/) {
                $genome=(split/\s+/,$_)[1];$genome=abs_path($genome);
        }
        if(/^Ref_ann/) {
                $gff=(split/\s+/,$_)[1];$gff=abs_path($gff);
        }
	if(/^medical/) {
		$medical=(split/\s+/,$_)[1];
	}
	if(/^Project_key/) {
		$key=(split/\s+/,$_)[1];
	}
}


$newgene_gff = abs_path($newgene_gff) if(defined $newgene_gff);

my $work_sh="$odir/work_sh";
`mkdir $work_sh`                unless( -d $work_sh);
my $res_path="$odir/each_geneRes";
`mkdir $res_path`	unless( -d $res_path);

###################extract info
#
my $TFBSdb="/share/nas1/lijj/develop/TFBStools/$key"."_TFBS_"."$medical";

open (FA,"$genome") or die "Open $genome failed:$!\n";
open (GFF,"$gff") or die "Open $gff failed:$!\n";

my @gene_id=();
my %gene_pos=();
my %strand=();
if ((-e $TFBSdb) && ($score eq "90%")) {
	if(defined $newgene_gff){
        	open(GFF2,"$newgene_gff") or die $!;
	        while(<GFF2>){
        	        chomp;next if(/^#/);
                	my @lines = split(/\s+/,$_);
	                next if ($lines[2] ne "gene");
        	        my $name=(split/\;/,(split/\=/,$lines[-1])[-1])[0];
	                push @gene_id,$name;
        	        $gene_pos{$name}{chr}=$lines[0];
			if ($lines[6] eq "+"){
		                $gene_pos{$name}{start}=$lines[3]-1000;
        		}else{
                		$gene_pos{$name}{start}=$lines[4];
	        	}
        	        $strand{$name}=$lines[6];
        	}
	        close(GFF2);
	}else{
		print "need $newgene_gff, Please Check..\n";
		die;
	}
}else{
	while(<GFF>){
        	chomp;next if(/^#/);
	        my @lines = split(/\s+/,$_);
        	next if ($lines[2] ne "gene");
	       	my $name=(split/\;/,(split/\=/,$lines[-1])[-1])[0];
		push @gene_id,$name;
		$gene_pos{$name}{chr}=$lines[0];
		if ($lines[6] eq "+"){
	               $gene_pos{$name}{start}=$lines[3]-1000;
	        }else{
        	       $gene_pos{$name}{start}=$lines[4];
        	}
        	$strand{$name}=$lines[6];
	}

	if(defined $newgene_gff){
        	open(GFF2,"$newgene_gff") or die $!;
	        while(<GFF2>){
        	        chomp;next if(/^#/);
                	my @lines = split(/\s+/,$_);
	                next if ($lines[2] ne "gene");
        	        my $name=(split/\;/,(split/\=/,$lines[-1])[-1])[0];
                	push @gene_id,$name;
	                $gene_pos{$name}{chr}=$lines[0];
			if ($lines[6] eq "+"){
                                $gene_pos{$name}{start}=$lines[3]-1000;
                        }else{
                                $gene_pos{$name}{start}=$lines[4];
                        }
        	        $strand{$name}=$lines[6];
	        }
        	close(GFF2);
	}
}
my %seq=();
my $name;
while(<FA>){
	chomp;
	if (/>(.+)/) {
		$name=(split /\s+/,$1)[0];
	} else {
		$seq{$name}.=$_;
	}
	
}

close (FA);
close (GFF);
close (CFG);

###############processing data
#
my $tfbs_predict="$Bin/TFBSTools_individual.R";
my $Rscript="/share/nas2/genome/biosoft/R/3.5.2/bin/Rscript";
#
#
open (OUT,">$odir/genes_upstream.fa") or die "$!\n";
open(RUN_SH,">$work_sh/TFBS_predict.sh")||die $!;

my $num=1;
my $cmd;
`mkdir -p $odir/tmp`  unless(-d "$odir/tmp");

foreach my $id (@gene_id) {
	my $sequence="";
	my $s=$strand{$id};
	my $chomo=$gene_pos{$id}{chr};
	my $up_s=$gene_pos{$id}{start};
	my $gene_fa;
	if($s eq "+"){
		my $fa_need=$seq{$chomo};
		my $fa = substr($fa_need, $up_s, 1000);
		$sequence.=$fa;
		$gene_fa=$id.".fa";
		if(!(-e $gene_fa)){
			open(GENE_FA,">$odir/tmp/$gene_fa")||die $!;
			print GENE_FA ">",$id,"\n",$sequence;
			close(GENE_FA);
		}
		if($num%20 == 1) {
			$cmd="$Rscript $tfbs_predict $spe_id $odir/tmp/$gene_fa $s $score $res_path ";
			$num++;
			next;
		}elsif($num%20 != 0) {
			$cmd.="&& $Rscript $tfbs_predict $spe_id $odir/tmp/$gene_fa $s $score $res_path ";
			if($num == scalar(@gene_id)) {
				print RUN_SH $cmd,"\n";
			}
		}else{
			$cmd.="&& $Rscript $tfbs_predict $spe_id $odir/tmp/$gene_fa $s $score $res_path";
			print RUN_SH $cmd,"\n";
		}
	}else{
		my $fa_need=reverse($seq{$chomo});
		my $fa = substr($fa_need, $up_s, 1000);
		$sequence.=$fa;
		$gene_fa=$id.".fa";
		if(!(-e $gene_fa)){
			open(GENE_FA,">$odir/tmp/$gene_fa")||die $!;
        	        print GENE_FA ">",$id,"\n",$sequence;
	                close(GENE_FA);
		}
		if($num%20 == 1) {
                        $cmd="$Rscript $tfbs_predict $spe_id $odir/tmp/$gene_fa $s $score $res_path ";
                        $num++;
                        next;
                }elsif($num%20 != 0) {
                        $cmd.="&& $Rscript $tfbs_predict $spe_id $odir/tmp/$gene_fa $s $score $res_path ";
                        if($num == scalar(@gene_id)) {
                                print RUN_SH $cmd,"\n";
                        }
                }else{
                        $cmd.="&& $Rscript $tfbs_predict $spe_id $odir/tmp/$gene_fa $s $score $res_path";
                        print RUN_SH $cmd,"\n";
                }
	}
	print OUT ">",$id,"\n",$sequence,"\n";
	$num++;
}

close (OUT);
close(RUN_SH);

&qsub("$work_sh/TFBS_predict.sh");
&qsubcheck("$work_sh/TFBS_predict.sh");
#`rm -rf $odir/tmp`;


#########################merge results and plot seqlogo of each gene
#
my $seqlogo="$Bin/TFBS_seqlogo.R";
`mkdir -p $odir/seqLogo`  unless (-d "$odir/seqLogo");

if (-e $TFBSdb){
	open(ALLRES,">$odir/newGenes_TFBS_predictRes.txt") ||die $!;
	print ALLRES "Model_id\tseqname\tstart\tend\tscore\tstrand\tframe\tTF\tclass\tsequence\tPvalue\n";
}else{ 
	open(ALLRES,">$odir/allGenes_TFBS_predictRes.txt") ||die $!;
	print ALLRES "Model_id\tseqname\tstart\tend\tscore\tstrand\tframe\tTF\tclass\tsequence\tPvalue\n";
}
my @files=();
push @files,glob("$res_path/*predictRes.txt");

foreach(my $i=0;$i<@files;$i++) {
	my $f=$files[$i];$f=abs_path($f);
	##########merge results
	open (G,"$f")||die $!;
	my $gene=basename($f);
        my $geneID=(split /\./,$gene)[0];
        my $file=$geneID."_plot.txt";
        open (P,">$odir/seqLogo/$file") ||die $!;

	while(<G>){
		chomp;next if(/^Model_id/);
		my @array=();
		my @tmp=split(/\t/,$_);
		
		my $TF=(split/\=/,(split/\;/,$tmp[-2])[0])[1];
		my $class=(split/\=/,(split/\;/,$tmp[-2])[1])[1];
		my $sequ=(split/\=/,(split/\;/,$tmp[-2])[-1])[1];
		print P $TF,"\t",$sequ,"\n";

		push @array,$tmp[0],$tmp[1],$tmp[2],$tmp[3],$tmp[4],$tmp[5],$tmp[6],$TF,$class,$sequ,$tmp[-1];
		print ALLRES join("\t",@array),"\n";

	}

	close(P);
	close(G);
}

close(ALLRES);

#########plot
my @pfiles=();
push @pfiles,glob("$odir/seqLogo/*plot.txt");
open(PLOT_SH,">$work_sh/TFBS_plot.sh")||die $!;

foreach(my $i=0;$i<@pfiles;$i++) {
	my $plot=$pfiles[$i];$plot=abs_path($plot);
	my $lines=`wc -l $plot`;
	my $num=(split(/\s+/,$lines))[0];
	next if($num == 0);
	my $gene=basename($plot);
	my $id=(split /\./,$gene)[0];
	$cmd="$Rscript $seqlogo $plot $id $odir/seqLogo ";
#	system($cmd);

	my $pdf=$id."_seqlogo.pdf";
	my $png=$id."_seqlogo.png";
	$cmd.="&& convert $odir/seqLogo/$pdf $odir/seqLogo/$png ";
	print PLOT_SH "$cmd\n";
	&run_or_die($cmd);
}
close(PLOT_SH);
#`rm $odir/seqLogo/*plot.txt`;

######################
###########
sub qsub{
        my $shfile= shift;
        my $queue="medical.q";
        my $cmd = "sh /share/nas2/genome/bmksoft/tool/qsub_sge_plus/v1.0/qsub_sge.plus.sh --reqsub --independent $shfile --queue $queue ";
        &run_or_die($cmd);
        return ;
}
sub qsubcheck{
        my $sh = shift;
        my @Check_file = glob "$sh*.qsub/*.Check";
        my @sh_file    = glob "$sh*.qsub/*.sh";
        if ( $#sh_file != $#Check_file ) {
                print "Their Some Error Happend in $sh qsub, Please Check..\n";
                die;
        }else {
                print "$sh qsub is Done!\n";
        }
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
sub USAGE       #use function
{
        my $usage=<<"USAGE";
    Description:
	Contact: lijj\@biomarker.com.cn
	version:   v1.0.0 at 2018.10.10 
	Description:  The prediction of transcriptor factor bingding sites of genes you providing;
    Usage:   
	Options:
	-cfg            Configure file, <required>;
                        For example:$Bin/detail.cfg;
	-gff2		the gff3 files of new genes;
	-score		the score for predicting TFBS of genes;
	-od		the output path
	-h		Help, display this help info.
	
USAGE
        print $usage;
        exit;  #die program
}

