use Getopt::Long;
use Getopt::Std;
use Config::General;
use File::Basename qw(basename dirname);
use FindBin qw($Bin $Script);
use Cwd qw(abs_path getcwd);
use List::Util qw/sum/;
my ($mirna,$target,$cfg,$odir);
GetOptions(
        "h|?"           =>\&USAGE,
	"mirna:s"	=>\$mirna,
	"target:s"	=>\$target,
	"cfg:s"     	=>\$cfg,
	"od:s"		=>\$odir,
)or &USAGE;
&USAGE unless ($odir and $mirna and $target and $cfg);

my %CFG=&readcfg("$Bin/CFG");
my $RNAHybrid_BIN	= $CFG{RNAHybrid};
my $testsvm		= $CFG{testsvm};
my $miranda_BIN		= $CFG{miranda};
my $targetscan		="$Bin/bin/targetscan_50.pl";
my $targetscan_context  ="$Bin/bin/targetscan_50_context_scores.pl";


`mkdir -p $odir`	unless(-d $odir);
$odir	=abs_path($odir);
$mirna	=abs_path($mirna);
$target	=abs_path($target);

$/=">";my %seqs=();
open(MIR,$mirna)||die $!;
open(NEW_MIR,">$odir/New.All_miRNA.expressed.fa");
while(<MIR>){
        chomp;next if($_=~/^$/);
        my ($id,$seq)=split(/\n/,$_,2);
	$seq=~tr/atcguU/ATCGTT/;
        $seqs{$id}=$seq;
	print NEW_MIR ">$id\n$seq";
}
close(MIR);
close(NEW_MIR);

$mirna = "$odir/New.All_miRNA.expressed.fa";
#########
&transSeq($target,"$odir/target.fa");
$target="$odir/target.fa";

#########read Config
$cfg  =abs_path($cfg);
my %para=&readcfg($cfg);

#########################
my $key=9606;
$key=$para{spe_id}	if(exists $para{spe_id});

$para{SPECIES_TYPE}	||= 0	if(!exists $para{SPECIES_TYPE});
## if the user forgot setting software options, choose targetscan and miRanda as default analysis tool
$para{RNAhybrid}	||= 0	if(!exists $para{RNAhybrid});
$para{miRanda}	  	||= 1	if(!exists $para{miRanda});
$para{targetscan} 	||= 1	if(!exists $para{targetscan});
$para{RNAhybrid_EN} 	||= -20	if(!exists $para{RNAhybrid_EN});
$para{miRanda_SCORE}	||= 50	if(!exists $para{miRanda_SCORE});
$para{miRanda_EN}	||= -20	if(!exists $para{miRanda_EN});
$para{miRanda_scale}	||= 4.0	if(!exists $para{miRanda_scale});
$para{miRanda_go}	||= -2	if(!exists $para{miRanda_go});
$para{miRanda_ge}	||= -8	if(!exists $para{miRanda_ge});

##
=cut
$/=">";my %seqs=();
open(MIR,$mirna)||die $!;
while(<MIR>){
	chop;next if($_=~/^$/);
	my ($id,$seq)=split(/\n/,$_,2);
	$seqs{$id}=$seq;
}
close(MIR);
##
=cut
my $work_sh="$odir/work_sh";
`mkdir $work_sh`	unless(-d $work_sh);
open(SH,">$work_sh/s1.miRNA_target.sh")||die $!;
if ($para{SPECIES_TYPE}==0){	
	my $miR_div="$odir/miR_div";
	`mkdir $miR_div`	unless(-d $miR_div);
	&run_or_die("split -l 1000 $mirna -a 4 -d $miR_div/miRNA_");
	my @miRNAs=glob("$miR_div/miRNA_*");	

	####################	RNAhybrid
	if($para{RNAhybrid} ==1) {
		`mkdir -p $odir/RNAhybrid/split`         unless(-d "$odir/RNAhybrid/split");
		my $gene_utr="$odir/RNAhybrid/target_utr.fa";
		&Gene_utr($target,$gene_utr);
		&run_or_die("split -l 1000 $gene_utr -a 4 -d $odir/RNAhybrid/split/target_");	
		my @tars=glob("$odir/RNAhybrid/split/target_*");
		for(my $i=0;$i<@tars;$i++){
			for(my $j=0;$j<@miRNAs;$j++){
				print SH "$RNAHybrid_BIN -m 50000 -d 1.9,0.28 -b 1 -e $para{RNAhybrid_EN} -t $tars[$i] -q $miRNAs[$j] >$odir/RNAhybrid/split/out.$i.$j && perl $Bin/bin/RNAhybrid_result.pl $odir/RNAhybrid/split/out.$i.$j $odir/RNAhybrid/split out.$i.$j && perl $testsvm $odir/RNAhybrid/split/out.$i.$j.RNAhybrid.aln.txt $odir/RNAhybrid/split/out.$i.$j.RNAhybrid.aln_svm.txt \n"
			}
		}
	}
	####################	targetscan
	if($para{targetscan} ==1){
		`mkdir -p $odir/targetscan/split`	unless(-d "$odir/targetscan/split");
		&targetFileFormat($target,"$odir/targetscan/target.file","UTR");
		&targetFileFormat($mirna,"$odir/targetscan/mirna.file","miRNA");
		&miR_context($mirna,"$odir/targetscan/miR_for_context_scores.txt");

		&run_or_die("split -l 500 $odir/targetscan/target.file -a 4 -d $odir/targetscan/split/target_");
		&run_or_die("split -l 500 $odir/targetscan/miR_for_context_scores.txt -a 4 -d $odir/targetscan/split/miRcontext_");
		&run_or_die("split -l 500 $odir/targetscan/mirna.file -a 4 -d $odir/targetscan/split/mirna_");
		my @tars=glob("$odir/targetscan/split/target_*");
		my @mirs=glob("$odir/targetscan/split/mirna_*");
		my @context=glob("$odir/targetscan/split/miRcontext_*");
		&run_or_die("sed -i 1\"imiRNA_family_ID\\tSpecies_ID\\tMiRBase_ID\\tMature_sequence\" $odir/targetscan/split/miRcontext_*");
		&run_or_die("sed -i 1\"imiRNA_family_ID\\tSpecies_ID\\tMiRBase_ID\\tMature_sequence\" $odir/targetscan/miR_for_context_scores.txt");

		for(my $i=0;$i<@tars;$i++){
			for(my $j=0;$j<@mirs;$j++){
				print SH "perl $targetscan $mirs[$j] $tars[$i] $odir/targetscan/split/out.$i.$j && perl $targetscan_context $context[$j] $tars[$i]  $odir/targetscan/split/out.$i.$j  $odir/targetscan/split/out.$i.$j.context\n";
			}
		}
	}
	###################	miranda
	if($para{miRanda} ==1){
        	`mkdir -p $odir/miranda/split` unless(-d "$odir/miranda/split");
		&run_or_die("split -l 1000 $target -a 4 -d $odir/miranda/split/target_ ");
		my $cmd="-sc $para{miRanda_SCORE} -en $para{miRanda_EN} -scale $para{miRanda_scale} -go $para{miRanda_go} -ge $para{miRanda_ge} -quiet ";
		$cmd .= " -strict " if (defined $para{miRanda_strict} && $para{miRanda_strict});
		my @tars=glob("$odir/miranda/split/target_*");
		for(my $i=0;$i<@tars;$i++){
			for(my $j=0;$j<@miRNAs;$j++){
				print SH "$miranda_BIN $miRNAs[$j] $tars[$i] -out $odir/miranda/split/out.$i.$j $cmd && perl $Bin/bin/miRanda_result.pl $odir/miranda/split/out.$i.$j $odir/miranda/split out.$i.$j.miRanda.aln \n"; 
			}
		}
	}
}else{
	`mkdir -p $odir/TargetFinder/split` unless(-d "$odir/TargetFinder/split");
	my $cmd=" -c $para{TargetFinder_score} ";
	$cmd .="-r "if($para{TargetFinder_rev}!=0);
	foreach my $id(keys %seqs){
		print SH "perl $Bin/bin/targetfinder.pl -s $seqs{$id} -d $target -q $id $cmd >$odir/TargetFinder/split/$id.TargetFinder.txt && perl $Bin/bin/targetFinder_result.pl $odir/TargetFinder/split/$id.TargetFinder.txt $odir/TargetFinder/split $id \n";
	}
}
close(SH);
&qsub("$work_sh/s1.miRNA_target.sh");

############ cat target prediction result
$/="\n";

my %mir=();
my %genes=();
if ($para{SPECIES_TYPE}==0) {
	my (%rnahybrid,%miranda,%targetscan);
	my $cmd="awk '{print \$1\"\\t\"\$2";
	if ($para{RNAhybrid} == 1){
		&run_or_die("cat $odir/RNAhybrid/split/*aln_svm.txt >$odir/RNAhybrid.aln_svm.txt");
		%rnahybrid=&RNAhybrid("$odir/RNAhybrid.aln_svm.txt");
		$cmd.="\"\\t\"\$3";
	}
	if($para{miRanda} == 1){
		&run_or_die("cat $odir/miranda/split/*.miRanda.aln.txt|grep -v '^#' > $odir/miRanda.aln.txt");
		&run_or_die("sed -i 1\"i#Seq1,Seq2,Tot Score,Tot Energy,Max Score,Max Energy,Strand,Len1,Len2,Positions\" $odir/miRanda.aln.txt");	
		%miranda=&Miranda("$odir/miRanda.aln.txt");
		$cmd.="\"\\t\"\$4";
	}
	if($para{targetscan}==1){
		&run_or_die("head -n 1 $odir/targetscan/split/out.0.0.context >$odir/targetscan.context.txt");
		&run_or_die("cat $odir/targetscan/split/out*context| grep -v '^Gene ID' >$odir/targetscan.context.txt.tmp && awk -F \$\"\\t\" '{if(\$10<0) {print \$0}}' $odir/targetscan.context.txt.tmp >> $odir/targetscan.context.txt && rm $odir/targetscan.context.txt.tmp");
		%targetscan=&scan("$odir/targetscan.context.txt");
		$cmd.="\"\\t\"\$5";
	}
	open(OUT,">$odir/miRNA_target_info.xls.tmp")||die $!;
	print OUT "#miRNA\ttarget\tRNAhybrid\tmiRanda\ttargetscan\n";
	foreach my $mi(keys %seqs){
		next if(!exists $miranda{$mi} && !exists $targetscan{$mi} && !exists $rnahybrid{$mi});
		foreach my $g(keys %genes){
			my($r1,$r2,$r3)=(0,0,0);
			$r1=1   if(exists $rnahybrid{$mi} && exists $rnahybrid{$mi}{$g});
			$r2=1	if(exists $miranda{$mi} && exists $miranda{$mi}{$g});
			$r3=1   if(exists $targetscan{$mi} && exists $targetscan{$mi}{$g});
			my $sum=$r1+$r2+$r3;
			print OUT "$mi\t$g\t$r1\t$r2\t$r3\n"	if($sum>0);
		}	
	}	
	close(OUT);
	&run_or_die("$cmd}' $odir/miRNA_target_info.xls.tmp >$odir/miRNA_target_info.xls ");
#	`rm $odir/miRNA_target_info.xls.tmp`;
	open(INFO,"$odir/miRNA_target_info.xls")||die $!;
	while(<INFO>){
		next if($_=~/^#/);
		my @tmp=split;
		my $mi=shift @tmp;my $g=shift @tmp;
		if((sum @tmp)==scalar(@tmp)){
			$mir{$mi}{$g}++;
		}	
	}
	close(INFO);
	
}else{
	&run_or_die("cat $miR_div/*.TargetFinder.aln.txt >$odir/TargetFinder.aln.txt");
	%mir=&RNAhybrid("$odir/TargetFinder.aln.txt");
}

%genes=();
## mir target list
open(OUT,">$odir/mir2target.list")||die $!;
print OUT "#ID\tTarget\n";
foreach my $mi(keys %mir){
	my @targets=keys %{$mir{$mi}};
	print OUT "$mi\t",join("\;",@targets),"\n";
	foreach my $g(@targets){$genes{$g}++;}
}
close(OUT);

## mir num and target num
open OUT,">$odir/mir2target.stat\n" || die $!;
print OUT "miR_num\t",scalar(keys %mir),"\n";
print OUT "target_gene_num\t",scalar(keys %genes),"\n";
close OUT;

###########subs
sub readcfg{
	my $cfg=shift;
	my %config=();
	$/="\n";
	open(IN,$cfg)||die $!;
	while(<IN>){
		chomp;next if($_=~/^#/);
		my @tmp=split(/\s+/,$_);
		$config{$tmp[0]}=$tmp[1];
	}
	close(IN);
	return %config;
}
sub transSeq{
        my($seq,$o)=@_;
        $/=">";
        open(SEQ,$seq)||die $!;
        open(OUT,">$o")||die $!;
        while(<SEQ>){
                chop;
                next if($_=~/^$/);
                my($id,$fa)=split(/\n/,$_,2);
                $id=(split(/\s+/,$id))[0];
                $fa=~s/\n//g;
                print OUT ">$id\n$fa\n";
        }
        close(SEQ);
        close(OUT);
}

sub Miranda{
        my $file=shift;
        my %mir=();
        open(ALN,$file)||die $!;
        while(<ALN>){
                chomp;next if($_!~/^>>/);
                my @tmp=split(/\s+/,$_);#>>miRNA	target
                $tmp[0]=~s/>>//g;
                $mir{$tmp[0]}{$tmp[1]}++;
		$genes{$tmp[1]}++;
        }
        return %mir;
}
sub RNAhybrid{
	my $file=shift;
	my %mir=();
	open(ALN,$file)||die $!;
	while(<ALN>){
		chomp;next if($_!~/^>/);
		my @tmp=split(/\s+/,$_);##>miRNA length target
		$tmp[0]=~s/>//g;
		$mir{$tmp[0]}{$tmp[2]}++;
		$genes{$tmp[2]}++;
	}
	return %mir;
}
sub scan{
	my $file=shift;
	my %mir=();
	open(FILE,$file)||die $!;
	while(<FILE>){
		chomp;
		next if($_=~/^Gene_ID/);
		my @tmp=split(/\s+/,$_);
		$mir{$tmp[2]}{$tmp[0]}++;
		$genes{$tmp[0]}++;
	}
	close(FILE);
	return %mir;
}

###########################################
sub targetFileFormat{
        my ($fa,$out,$type)=@_;
        $/=">";
        open(FA,$fa)||die $!;
        open(OUT,">$out")||die $!;
        while(<FA>){
                chop;
                next if($_=~/^$/);
                my ($id,$seq)=split(/\n/,$_,2);
		$seq=~s/\n//g;
                if($type eq "UTR"){
                        print OUT "$id\t$key\t$seq\n";
                }elsif($type eq "miRNA"){
                        print OUT "$id\t",substr($seq,1,7),"\t$key\n";
                }
        }
        close(FA);
        close(OUT);
}
sub miR_context{
        my ($i,$o)=@_;
        $/=">";
        open(IN,$i)||die $!;
        open(OUT,">$o")||die $!;
#       print OUT "miRNA_family_ID\tSpecies_ID\tMiRBase_ID\tMature_sequence\n";
        while(<IN>){
                chop;next if($_=~/^$/);
                my ($id,$seq)=split(/\n/,$_,2);
		$seq=~s/\n//g;
                print OUT "$id\t$key\t$id\t$seq\n";
                        
        }
        close(IN);
        close(OUT);
}

############################################
sub Gene_utr
{ 
  
  my ($gene_fa,$utr_fa)=@_;

  open IN,"$gene_fa" || die $!;
  open OUT,">$utr_fa" || die $!;
  $/='>';
  while (<IN>) {
      chop;
      next if (/^$/);
      my ($id,$seq)=split /\n+/,$_,2;
      $id=(split(/\s+/,$id))[0];
      $seq=~s/\n//g;
      my $length=length $seq;
      my $new_seq;
      if ($length>50000) {
        $new_seq=substr ($seq,0,50000);
        print OUT ">$id\n$new_seq\n";
      }
      else {
        print OUT ">$id\n$seq\n";
      }
  }
  close IN;
  close OUT;
  $/="\n";
}

sub intersection {#
   my ($A,$B)=@_;
   my %uniqA=map {$_,1} @{$A};
   my %uniqB=map {$_,1} @{$B};
   my %merge=();
   my %overlap=();
   foreach  (keys %uniqA,keys %uniqB) {
           $merge{$_}++ && $overlap{$_}++;
   }
   my @result = keys %overlap;
   return \@result;
}


sub qsub{
        my $shfile= shift;
	my $queue="medical.q";
	$queue=$para{Queue_type}	if(exists $para{Queue_type});
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
}

sub readConfig{
        my $configFile=shift;
        my $d=Config::General->new(-ConfigFile => "$configFile");
        my %config=$d->getall;
        return %config;
}
sub USAGE{
        my $usage=<<"USAGE";
Program: $0
Version: 1.0.0
Contact: wenyh\@biomarker.com.cn
Usage:
        Options:
        -miRNA          miRNA.fa
        -target         lncRNA/mRNA/cirRNA target.fa
        -od             output dir
        -cfg            config

        -h      Help

Example:

USAGE
        print $usage;
        exit;
}

