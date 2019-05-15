#!/usr/bin/perl -w
use strict;
use warnings;
use Getopt::Long;
use Data::Dumper;
use FindBin qw($Bin $Script);
use File::Basename qw(basename dirname);
use newPerlBase;
use Cwd qw(abs_path);
my $BEGIN_TIME=time();
my $version="1.0.0";
#######################################################################################
#
#######################################################################################
# ------------------------------------------------------------------
# GetOptions
# ------------------------------------------------------------------
my ($idir,$odir,$cfg,$ce);

GetOptions(
				"help|?" =>\&USAGE,
				"od:s"=>\$odir,
				"in:s"=>\$idir,
				"ce:s"=>\$ce,
				"cfg:s"=>\$cfg,
				) or &USAGE;
&USAGE unless ($idir and $cfg and $ce);
# ------------------------------------------------------------------
$cfg = abs_path($cfg);
$idir = abs_path($idir);
`mkdir $odir` unless (-d "$odir");
$odir = abs_path($odir);

`mkdir $odir/work_sh` unless (-d "$odir/work_sh");

my $work_sh = "$odir/work_sh";
$ce = abs_path($ce);
print"$ce\n";

my %sample = &relate($cfg);
my %config = ();
&readcfg($cfg);
my @diff;my %vs;
open (CFG,"$cfg") or die $!;
print "read cfg and diff group is:\n";
while (<CFG>){
	chomp;
	next if(/^#|^$|^\s$/);
	my @tmp = split /\s+/,$_;
	if($tmp[0] eq "Diff"){
		push @diff,$tmp[1];
		print "$tmp[1]\n";
	}
}
close(CFG);

my @RNA = keys %{@sample{keys %sample}};
print "@RNA\n";
my (%h,%LSM,%CSM);
#$vs{$group}{$rna}=$Group	$vs{W01_vs_W02}{$gene}=L01_vs_L02
if(scalar(@RNA)==4){
foreach my $k1 (@diff){
	&switch($k1);
	foreach my $k2 (@RNA){
		$h{$k1}{$k2}= "$idir/$k2.$vs{$k1}{$k2}.DEG_final.xls";
		print "$idir/$k2.$vs{$k1}{$k2}.DEG_final.xls\n";
	}
}

open(SH,">$work_sh/ce_cytoscape.sh");
foreach my $vs (@diff){
	system ("cat $h{$vs}{gene} $h{$vs}{lncRNA} $h{$vs}{sRNA} > $odir/$vs.lncRNA-miRNA-mRNA.all.DEG.xls");
#	print "cat $h{$vs}{gene} $h{$vs}{lncRNA} $h{$vs}{sRNA} > $odir/$vs.lncRNA-miRNA-mRNA.all.DEG.xls\n";
	system ("cat $h{$vs}{gene} $h{$vs}{circRNA} $h{$vs}{sRNA} > $odir/$vs.circRNA-miRNA-mRNA.all.DEG.xls");
#	print "cat $h{$vs}{gene} $h{$vs}{circRNA} $h{$vs}{sRNA} > $odir/$vs.circRNA-miRNA-mRNA.all.DEG.xls\n";
	`mkdir $odir/lncRNA-miRNA-mRNA` unless (-d "$odir/lncRNA-miRNA-mRNA");
	`mkdir $odir/circRNA-miRNA-mRNA` unless (-d "$odir/circRNA-miRNA-mRNA");
	print SH "perl $Bin/abstract.pl -deg $odir/$vs.lncRNA-miRNA-mRNA.all.DEG.xls -ce $ce -o $odir/lncRNA-miRNA-mRNA/$vs.lncRNA-miRNA-mRNA.ceRNA_pair_adjust_p_Sig_diff.xls\n";
	print SH "perl $Bin/abstract.pl -deg $odir/$vs.circRNA-miRNA-mRNA.all.DEG.xls -ce $ce -o $odir/circRNA-miRNA-mRNA/$vs.circRNA-miRNA-mRNA.ceRNA_pair_adjust_p_Sig_diff.xls\n";
}
close(SH);

`sh /share/nas2/genome/bmksoft/tool/qsub_sge_plus/v1.0/qsub_sge.plus.sh --reqsub --independent --queue medical.q $work_sh/ce_cytoscape.sh`;
#system("sh $work_sh/ce_cytoscape.sh > $work_sh/ce_cytoscape.sh.log");
}
######################################################################################l
print STDOUT "\nDone. Total elapsed time : ",time()-$BEGIN_TIME,"s\n";
#######################################################################################

# ------------------------------------------------------------------
# sub function
# ------------------------------------------------------------------
################################################################################################################
sub run_or_die(){
        my ($cmd) = @_ ;
        &show_log($cmd);
        my $flag = system($cmd) ;
        if ($flag != 0){
                &show_log("Error: command fail: $cmd");
                exit(1);
        }
        &show_log("done.");
}
sub show_log(){
        my ($txt) = @_ ;
        my $time = time();
        my($sec, $min, $hour, $day, $mon, $year, $wday, $yday, $isdst) = localtime($time);
        $wday = $yday = $isdst = 0;
        my $Time=sprintf("%4d-%02d-%02d %02d:%02d:%02d", $year+1900, $mon+1, $day, $hour, $min, $sec);
        print "$Time:\t$txt\n" ;
}
########################
sub readcfg{
        my $cfg=shift;
        open(CFG,$cfg)||die $!;
        while(<CFG>){
                next if($_=~/^#|^\s$/);
                chomp;my @tmp=split(/\s+/,$_);
                if($tmp[0] eq "Com"||$tmp[0] eq "Sep"){
                        push @{$config{$tmp[0]}},$tmp[1];
                }else{
                        $config{$tmp[0]}=$tmp[1];
                }
        }
        close(CFG);
}
#######################
sub relate{
        my $cfg=shift;
        my %sample=();
	my @rnas=();
        open(CFG,$cfg)||die $!;
        while(<CFG>){
                chomp;next if($_!~/^Sample/);
                my @tmp=split(/\s+/,$_);shift @tmp;
                if($tmp[0] eq "ID" && scalar(@rnas)==0){
                        shift @tmp;@rnas=@tmp;
                }else{
                        my $id=shift @tmp;
                        for(my $i=0;$i<@tmp;$i++){
                                $sample{$id}{$rnas[$i]}=$tmp[$i];
                        }
			if(exists $sample{$id}{lncRNA} && !exists $sample{$id}{gene}){
	                        $sample{$id}{gene}=  $sample{$id}{lncRNA};
			}
                }
        }
        close(CFG);
        return %sample; #relation{W01}{lncRNA}="L01";
}


sub switch {
	my $group = shift;
	my ($v1,$v2,@V1,@V2);
	
		($v1,$v2) = split /_vs_/,$group,2;
		@V1 = split /_/,$v1;
		@V2 = split /_/,$v2;
		foreach my $key ("sRNA","lncRNA","circRNA","gene"){
			my (@s1,@s2);
			for(my $i=0;$i<@V1;$i++){
				push @s1,$sample{$V1[$i]}{$key};
			}
			for(my $j=0;$j<@V2;$j++){
				push @s2,$sample{$V2[$j]}{$key};
			}
			my $S1 = join ("_",@s1);
			my $S2 = join ("_",@s2);
			my $Group = join("_vs_",$S1,$S2);
			$vs{$group}{$key}=$Group;
		}
}

sub getCyto {
	my $inter = shift;
	my $att = shift;
	my $line = shift;
	my $Cyto = shift;
	open (O1,">$inter");
	open (O2,">$att");
	my %ind;
	my ($center,$type)=split //,$line;
	foreach my $rna(sort keys %$Cyto){
		foreach my $stat (sort keys %{$$Cyto{$rna}}){
			if (!exists $ind{$rna}){
				print O2 "$rna\t$center\t$stat\n";
				$ind{$rna}=1;
			}
			my @p = split /\;/,$$Cyto{$rna}{$stat};
			for (my $i=0;$i<@p;$i++){
				my ($id,$ud,$k) = split /\(|\)/,$p[$i];
				print O1 "$rna\t$id\t$line\n";
				if (!exists $ind{$id}){
					print O2 "$id\t$type\t$ud\n";
					$ind{$id}=1;
				}
			}
		}
	}
	close(O1);close(O2);
	return ($inter,$att);
}

################################################################################################################

sub GetTime {
	my ($sec, $min, $hour, $day, $mon, $year, $wday, $yday, $isdst)=localtime(time());
	return sprintf("%4d-%02d-%02d %02d:%02d:%02d", $year+1900, $mon+1, $day, $hour, $min, $sec);
}

################################################################################################################
sub USAGE {
	my $usage=<<"USAGE";
 ProgramName:
     Version:	$version
     Contact:	Niepy <niepy\@biomarker.com.cn> 
Program Date:	2018.05.23
      Modify:	
 Description:	This program is used to abstract diff.lncRNA-miRNA-mRNA and diff.circRNA-miRNA-mRNA form ceRNA result
       Usage:
		Options:
		-in <dir>	input dir ,DEG_Analysis,force
		-od <dir>	output dir , Combine/Cytoscape,force
		-ce <file>	ceRNA result,Coexpression_ceRNA_pair_adjust_p_Sig.txt(sample>=5) or ceRNA_pair_adjust_p_Sig.txt
		-cfg <file>	cfg,force
		-h		help

USAGE
	print $usage;
	exit;
}
