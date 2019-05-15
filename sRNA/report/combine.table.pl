#!/usr/bin/perl -w
use strict;
use warnings;
use Getopt::Long;
use Data::Dumper;
use FindBin qw($Bin $Script);
use File::Basename qw(basename dirname);
use newPerlBase;
my $BEGIN_TIME=time();
my $version="1.0.0";
#######################################################################################
#
#######################################################################################
# ------------------------------------------------------------------
# GetOptions
# ------------------------------------------------------------------
my ($idir,$odir,$data_cfg,$detail_cfg,);

GetOptions(
				"help|?" =>\&USAGE,
				"od:s"=>\$odir,
				"in:s"=>\$idir,
				"data_cfg:s"=>\$data_cfg,
				"detail_cfg:s"=>\$detail_cfg,
				) or &USAGE;
&USAGE unless ($idir and $odir and $data_cfg and $detail_cfg);
# ------------------------------------------------------------------

$data_cfg = &ABSOLUTE_DIR($data_cfg);
$detail_cfg = &ABSOLUTE_DIR($detail_cfg);
$idir = &ABSOLUTE_DIR($idir);
`mkdir $odir` unless (-d "$odir");
$odir = &ABSOLUTE_DIR($odir);

my %sample = &relate($data_cfg);
my %config = ();
&readcfg($detail_cfg);
my @diff;my %vs;
open (CFG,"$detail_cfg") or die $!;
print ("read cfgi and diff group is:\n");
while (<CFG>){
	chomp;
	next if(/^#|^$|^\s$/);
	my @tmp = split /\s+/,$_;
	if($tmp[0] eq "Sep"){
		$tmp[1]=~s/\;/_vs_/;
		$tmp[1]=~s/,/_/g;
		push @diff,$tmp[1];
		print "$tmp[1]\n";
	}
	if($tmp[0] eq "Com"){
		$tmp[1]=~s/,/_vs_/;
		push @diff,$tmp[1];
		print "$tmp[1]\n";
	}
}
close(CFG);

my @RNA = ("gene","lncRNA","circRNA","sRNA");
my (%h,$gnum,$lnum,$cnum,$snum);
foreach my $k1 (@diff){
	foreach my $k2 (@RNA){
		$h{$k1}{$k2}= "$idir/$k2/$k1.DEG.xls";
	}
}

my (%File);my $flag4=0;my $flag3=0; my $flag2=0;
foreach my $group(@diff){
        if(exists $h{$group}{gene}){
                $gnum = `less $h{$group}{gene} |wc -l`;
                push @{$File{$group}{4}},$h{$group}{gene} if($gnum>3);
                push @{$File{$group}{3}},$h{$group}{gene} if($gnum==3);
                push @{$File{$group}{2}},$h{$group}{gene} if($gnum==2);
        }
        if(exists $h{$group}{lncRNA}){
                $lnum = `less $h{$group}{lncRNA} |wc -l`;
                push @{$File{$group}{4}},$h{$group}{lncRNA} if($lnum>3);
                push @{$File{$group}{3}},$h{$group}{lncRNA} if($lnum==3);
                push @{$File{$group}{2}},$h{$group}{lncRNA} if($lnum==2);
        }
        if(exists $h{$group}{circRNA}){
                $cnum = `less $h{$group}{circRNA} |wc -l`;
                push @{$File{$group}{4}},$h{$group}{circRNA} if($cnum>3);
                push @{$File{$group}{3}},$h{$group}{circRNA} if($cnum==3);
                push @{$File{$group}{2}},$h{$group}{circRNA} if($cnum==2);
        }
        if(exists $h{$group}{sRNA}){
                $snum = `less $h{$group}{sRNA} |wc -l`;
                push @{$File{$group}{4}},$h{$group}{sRNA} if($snum>3);
                push @{$File{$group}{3}},$h{$group}{sRNA} if($snum==3);
                push @{$File{$group}{2}},$h{$group}{sRNA} if($snum==2);
        }
}
foreach my $g(keys %File){
        foreach my $num (keys %{$File{$g}}){
                if($num==4){
                        if(scalar(@{$File{$g}{4}})==scalar(@RNA)){$flag4++;}
                }
                if($num==3){
                        if(scalar(@{$File{$g}{3}})==scalar(@RNA)){$flag3++;}
                }
                if($num==2){
                        if(scalar(@{$File{$g}{2}})==scalar(@RNA)){$flag2++;}
                }
        }
}

if($flag4>=1){
        foreach my $vs(keys %File){
                print "$vs\n";
                foreach my $type(@RNA){
                        `head -n 4 $h{$vs}{$type} > $odir/diff_$type.xls`;
                }
                last;
        }
}else{
        if($flag3>=1){
                foreach my $vs(keys %File){
                        print "$vs\n";
                        foreach my $type(@RNA){
                                `head -n 4 $h{$vs}{$type} > $odir/diff_$type.xls`;
                        }
                        last;
                }
        }else{
                if($flag2>=1){
                        foreach my $vs(keys %File){
                                print "$vs\n";
                                foreach my $type(@RNA){
                                        `head -n 4 $h{$vs}{$type} > $odir/diff_$type.xls`;
                                }
                        }
                }else{
                        foreach my $type(@RNA){
                                my $f=0;
                                foreach my $vs(keys %File){
                                        my $line = `less $h{$vs}{$type}|wc -l`;chomp($line);
                                        if($line>=2){`head -4 $h{$vs}{$type} > $odir/diff_$type.xls`;$f++;}
                                        last if($f>=1);
                                }
                        }
                }
        }
}

if(scalar(@RNA)==4){
my @ce = ("lncRNA-miRNA-mRNA","circRNA-miRNA-mRNA");
my (%h1,$lgnum,$cgnum);
foreach my $n (@diff){
        foreach my $m(@ce){
                $h1{$n}{$m}= "$idir/$m/$n.$m.ceRNA_pair_adjust_p_Sig_diff.xls";
        }
}

foreach my $ce(@ce){
        my $ff=0;
        foreach my $group1(@diff){
                my $l = `less $h1{$group1}{$ce}|wc -l`;chomp($l);
                if($l>1){`head -4 $h1{$group1}{$ce} > $odir/$ce.xls`;$ff++;}
                last if($ff>0);
        }
}
}

=cut
my $flag=0;
foreach my $group(@diff){
	$gnum = `less $h{$group}{gene} |wc -l `;
	$lnum = `less $h{$group}{lncRNA} |wc -l`;
	$cnum = `less $h{$group}{circRNA} |wc -l`;
	$snum = `less $h{$group}{sRNA} |wc -l`;

	if(($gnum >1)&&($lnum>1)&&($cnum>1)&&($snum>1)){
		$flag ++;
	}
	if($flag ==1){
		`head -n 4 $h{$group}{gene} > $odir/diff_gene.xls`;
		`head -n 4 $h{$group}{lncRNA} > $odir/diff_lncRNA.xls`;
		`head -n 4 $h{$group}{circRNA} > $odir/diff_circRNA.xls`;
		`head -n 4 $h{$group}{sRNA} > $odir/diff_miRNA.xls`;
	}
}

my @ce = ("lncRNA-miRNA-mRNA","circRNA-miRNA-mRNA");
my (%h1,$lgnum,$cgnum);
foreach my $n (@diff){
	foreach my $m(@ce){
		$h1{$n}{$m}= "$idir/$m/$n.$m.ceRNA_pair_adjust_p_Sig_diff.xls";
	}
} 

my $flag1=0;my $flag2=0;
foreach my $group1(@diff){
	$lgnum = `less $h1{$group1}{"lncRNA-miRNA-mRNA"} |wc -l`;
	$cgnum = `less $h1{$group1}{"circRNA-miRNA-mRNA"} |wc -l`;
	if($lgnum >1){
		$flag1++;
	}
	if($flag1 ==1){
		`head -n 4 $h1{$group1}{"lncRNA-miRNA-mRNA"} > $odir/lncRNA-miRNA-mRNA.xls`;
	}
	if($cgnum >1){
		$flag2++;
	}
	if ($flag2==1){
		`head -n 4 $h1{$group1}{"circRNA-miRNA-mRNA"} > $odir/circRNA-miRNA-mRNA.xls`;
	}
}
=cut
######################################################################################l
print STDOUT "\nDone. Total elapsed time : ",time()-$BEGIN_TIME,"s\n";
#######################################################################################

# ------------------------------------------------------------------
# sub function
# ------------------------------------------------------------------
################################################################################################################
sub qsub(){
        my $shfile= shift;
        my $cmd = "sh /share/nas2/genome/bmksoft/tool/qsub_sge_plus/v1.0/qsub_sge.plus.sh --reqsub --independent $shfile --queue $config{Queue_type}";
        &run_or_die($cmd);
}
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
        my $file=shift;
        open(REL,$file)||die $!;
        my %sample=();
        while(<REL>){
                chomp;
                next if($_ !~/^SampleID/);
                my @tmp=split(/\s+/,$_);
                $sample{$tmp[-1]}{sRNA}=$tmp[2];
                $sample{$tmp[-1]}{lncRNA}=$tmp[3];
                $sample{$tmp[-1]}{circRNA}=$tmp[3];
                $sample{$tmp[-1]}{gene}=$tmp[3];
        }
        close(REL);
        return %sample;
}

sub switch {
	my $group = shift;
	my ($v1,$v2,@V1,@V2);
	
		($v1,$v2) = split /_vs_/,$group,2;
		@V1 = split /_/,$v1;
		@V2 = split /_/,$v2;
		foreach my $key ("sRNA","lncRNA","circRNA","gene","lncRNA-miRNA-mRNA","circRNA-miRNA-mRNA"){
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

sub ABSOLUTE_DIR{ #$pavfile=&ABSOLUTE_DIR($pavfile);
	my $cur_dir=`pwd`;chomp($cur_dir);
	my ($in)=@_;
	my $return="";
	if(-f $in){
		my $dir=dirname($in);
		my $file=basename($in);
		chdir $dir;$dir=`pwd`;chomp $dir;
		$return="$dir/$file";
	}elsif(-d $in){
		chdir $in;$return=`pwd`;chomp $return;
	}else{
		warn "Warning just for file and dir\n";
		exit;
	}
	chdir $cur_dir;
	return $return;
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
		-data_cfg <file>	data_cfg,force
		-detail_cfg <file>	detail_cfg,force
		-h		help

USAGE
	print $usage;
	exit;
}
