#!/usr/bin/perl -w
#
my $ver="1.0";
my $title = "fusionmap";
use strict;
use Getopt::Long;
use Data::Dumper;
use Config::General;
use FindBin qw($Bin $Script);
use File::Basename qw(basename dirname);
use newPerlBase;
# ------------------------------------------------------------------
# GetOptions
# ------------------------------------------------------------------
my ($indir,$type,$odir,$script_cfg,$test,$detail_cfg);

GetOptions(
	"indir:s"	=>\$indir,
	"cfg2:s"	=>\$detail_cfg,
	"script_cfg:s"	=>\$script_cfg,
	"od:s"		=>\$odir,
	"type:s"	=>\$type,
	"test"		=>\$test,
	"help|h"	=>\&USAGE,
    ) or &USAGE;
&USAGE unless ($indir and $type and $odir and $script_cfg and $detail_cfg );
#$script_cfg=&ABSOLUTE_DIR($script_cfg);
my %script_config=%{readconf("$script_cfg")};

#read detail config

my %detail_cfg=&readConfig($detail_cfg);
system "mkdir -p $odir" unless (-d $odir);

$odir=&ABSOLUTE_DIR($odir);
$indir=&ABSOLUTE_DIR($indir);

createLog($title,$ver,$$,"$odir/log/",1);
timeLog ();
my %para;
my %sample;
my @Bam  = glob ("$indir/*");
my @Bam_file;
my @ID;
for (my $i=0;$i<@Bam ;$i++) {
	my $id = (split /\//,$Bam[$i])[-1];
	my $bam = "$Bam[$i]/$id.HISAT_aln.sorted.bam";
	push @Bam_file ,$bam ;
	push @ID, $id;
	print "$bam\n";
}
my $bowtie_path =$script_config{bowtie};
my $fusionmap_path =$script_config{fusionmap};
my $mono_path =$script_config{mono};
my $samtools_path =$script_config{samtools};
my $queue = $script_config{que};
my $vf = $script_config{vf};
my $cpu = $script_config{cpu};
my $Ref	= $script_config{fusionmapRef};

my $step;
my $work_sh = "$odir/work_sh";
mkdir $work_sh;
stepStart(1,"bowtie");
# ------------------------------------------------------------------
# read configs and steps
# ------------------------------------------------------------------

stepStart(1,"fuisonmap_Analysis");
my $shfile = "$work_sh/fusionmap.sh" ;
open (SH, ">$shfile") || die "Can't creat $shfile, $!\n" ;
mkdirOrDie("$odir/fusionmap_result");
 `ln -sf $Bin/Fusion $odir/fusionmap_result/ `;
mkdirOrDie("$odir/fusionmap_result/ReferenceLibrary");
if ($type eq "Human.B38") {
	`ln -sf $Ref/Human.B38* $odir/fusionmap_result/ReferenceLibrary`;
}elsif ($type eq "Human.hg19") {
	`ln -sf $Ref/Human.hg19* $odir/fusionmap_result/ReferenceLibrary`;
}elsif ($type eq "Mouse.B38") {
	`ln -sf $Ref/Mouse.B38* $odir/fusionmap_result/ReferenceLibrary`;
}elsif ($type eq "Rat.B5.0") {
	`ln -sf $Ref/Rat.B5.0* $odir/fusionmap_result/ReferenceLibrary`;
}elsif ($type eq "Rat.B6.0") {
        `ln -sf $Ref/Rat.B6.0* $odir/fusionmap_result/ReferenceLibrary`;
}


#`ln -sf $Bin/ReferenceLibrary $odir/fusionmap_result/`;
for(my $i=0;$i<@Bam_file;$i++){
  mkdirOrDie("$odir/fusionmap_result/$ID[$i]");
  print SH "perl $Bin/fusion_config/FusionMapconfigCreater.pl --bamfile $Bam_file[$i] --outputpath $odir/fusionmap_result/$ID[$i] --outputname $ID[$i] --datatype PE && export PATH=/share/nas1/liaoqp/project/xudd/softwares/mono/2.10.9/bin:\$PATH &&";
  if ($type eq "Human.B38") {
      print SH "$mono_path $fusionmap_path --pereport $odir/fusionmap_result Human.B38 Ensembl.R83 $odir/fusionmap_result/$ID[$i]/FusionPE.config > $odir/fusionmap_result/$ID[$i]/log.txt";
  }elsif ($type eq "Human.hg19" ) {
	  print SH "$mono_path $fusionmap_path --pereport $odir/fusionmap_result Human.hg19 RefGene $odir/fusionmap_result/$ID[$i]/FusionPE.config > $odir/fusionmap_result/$ID[$i]/log.txt";
  }elsif ($type eq "Mouse.B38") {
	  print SH "$mono_path $fusionmap_path --pereport $odir/fusionmap_result Mouse.B38 Ensembl.R78 $odir/fusionmap_result/$ID[$i]/FusionPE.config > $odir/fusionmap_result/$ID[$i]/log.txt";
  }elsif ($type eq "Rat.B5.0") {
	  print SH "$mono_path $fusionmap_path --pereport $odir/fusionmap_result Rat.B5.0 Ensembl.R78 $odir/fusionmap_result/$ID[$i]/FusionPE.config > $odir/fusionmap_result/$ID[$i]/log.txt";
  }elsif ($type eq "Rat.B6.0") {
          print SH "$mono_path $fusionmap_path --pereport $odir/fusionmap_result Rat.B6.0 Ensembl.R84 $odir/fusionmap_result/$ID[$i]/FusionPE.config > $odir/fusionmap_result/$ID[$i]/log.txt";
  }
  print SH "&&cat $odir/fusionmap_result/$ID[$i]/$ID[$i]\.PairedEndFusionReport.txt | awk '\$4>20 {print \$0}'|egrep -v 'InBlackList|InFamilyList|InParalogueList|SameEnsemblGene' >$odir/fusionmap_result/$ID[$i]/$ID[$i]\.FusionReport_filter.txt&&sed  's/HISAT_aln.sorted.bam.//g' -i $odir/fusionmap_result/$ID[$i]/$ID[$i]\.FusionReport_filter.txt\n";
}
close SH;
#die;
qsubOrDie("$odir/work_sh/fusionmap.sh",$queue,$cpu,$vf);
stepTime(1,"fusionmap_Analysis");

stepStart(2,"circos_plot and result");
$shfile = "$work_sh/circos.sh";
open (SH,">$shfile") || die "Can't creat $shfile, $!\n" ;
mkdirOrDie("$odir/result/png");
mkdirOrDie("$odir/result/report");
#foreach my $id (sort keys %sample){
for (my $i=0;$i<@ID ;$i++) {
 mkdirOrDie("$odir/result/$ID[$i]");
 print SH "perl $Bin/circos/link.pl --cfg $Bin/circos/path.cfg -od $odir/result/$ID[$i] --fr $odir/fusionmap_result/$ID[$i]/$ID[$i]\.FusionReport_filter.txt --type $type ";
 print SH "&& cp $odir/result/$ID[$i]/circos\.png $odir/result/png/$ID[$i]\_circos\.png ";
 print SH "&& cp $odir/result/$ID[$i]/FusionReport.txt $odir/result/report/$ID[$i]\_FusionReport.txt\n";
}
close SH;
qsubOrDie("$odir/work_sh/circos.sh",$queue,$cpu,$vf);
stepTime(2,"circos_plot and result");
############################################3

=head
stepStart(3,"pfam annotation");
my $genome= $detail_cfg{Ref_seq};
$shfile = "$work_sh/pfam.sh";
mkdirOrDie("$odir/pfam");
mkdirOrDie("$odir/pfam/mid");
&run_or_die("cat $odir/result/report/*_FusionReport.txt | less | grep -v FusionID | sort | uniq > $odir/pfam/fusion.gtf");
open(TMPGTF,"$odir/pfam/fusion.gtf")||die $!;
my @gtffile=<TMPGTF>;
close(TMPGTF);
my $gene_num=@gtffile;
if ($gene_num<=30000){
&run_or_die("perl $Bin/pfam/get_fusion_fa.pl -fa $genome -gtf  $odir/pfam/fusion.gtf -o $odir/pfam/fusion.fa");
&run_or_die("perl $Bin/pfam/cut_fa_file_to_dir.pl -mRNA $odir/pfam/fusion.fa -od $odir/pfam/mid -cut 200");
open (SH,">$shfile") || die "Can't creat $shfile, $!\n" ;
print "$detail_cfg{Pfam}";
print SH "perl $Bin/pfam/Gene_Pfam_Anno.pl -idir  $odir/pfam/mid -pfam $detail_cfg{Pfam} --odir $odir/pfam/ --cpu 4 -queue  general.q";
close SH;
qsubOrDie("$odir/work_sh/pfam.sh",$queue,$cpu,$vf);
print "perl $Bin/pfam/extrac_pfam.pl -i $odir/pfam/Result/fusion.fa.Pfam.anno.details -o $odir/result/fusion.fa.Pfam.anno\n";
`perl $Bin/pfam/extrac_pfam.pl -i $odir/pfam/Result/fusion.fa.Pfam.anno.details -o $odir/result/Fusion_gene_Pfam.anno`;
}
else {print "the number of fusiongenes > 30000 and the pfam annotation is skipped";
}
stepTime(3,"pfam annotation");
=cut
#######################################################################################
# ------------------------------------------------------------------
# sub function
# ------------------------------------------------------------------
#############################################################################################################
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
##########################read config

sub readConfig{
	my $configFile=shift;
	my $d=Config::General->new(-ConfigFile => "$configFile");
	my %config=$d->getall;	
	return %config;
}
#############################################################################################################
sub steps_process {
    print "\n";
    &log_current_time("step check:");
    my ($step_str, $step_by_step, $step_hash) = @_;
    my %_step_ = (
        '1' => 'Data_Assess',
        '2' => 'Basic_Analysis',
        '3' => 'Anno_Integrate',
        '4' => 'DEG_Analysis',
        '5' => 'Alitsplice_Analysis',
        '6' => 'SNP_Analysis',
        '7' => 'Gene_Structure_Optimize',
        '8' => 'Analysis_Report',
        '9' => 'More',
    );

    if ($step_by_step) {
        print "step-by-step: ON\n";
        if ($step_str =~/^[1-9]$/) {
            for my $s ($step_str..9) {
                $step_hash->{$s} = $_step_{$s};
            }
        } else {
            print "ERROR: illegal steps specified by --step, or options specified by --step and --sbs clashed. \n";
            die;
        }
    } else {
        if ($step eq join ',',(1..9)) {
            print "step-by-step: ON\n";
        } else {
            print "step-by-step: OFF\n";
        }

        for my $s (split /,/,$step_str) {
            if ($s =~/^[1-9]$/) {
                $step_hash->{$s} = $_step_{$s};
            } else {
                print "ERROR: illegal steps specified by --step.\n";
                die;
            }
        }
    }

    print "steps_to_run: ".(join ", ",(map {sprintf "$_.$step_hash->{$_}"} sort keys %{$step_hash})).".\n";
    &log_current_time("step check done.\n");
}

#############################################################################################################
sub step_cmd_process {
    my ($cmd, $step_hash, $step_n, $sh_dir) = @_;
    my $sh_file = "$sh_dir/$step_n.$step_hash->{$step_n}.sh";
    my $log_file = "$sh_file.log";
    my $flag = 0;
    my $start_time = time();
    &log_current_time("step$step_n. $step_hash->{$step_n} start ...");
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
        &log_current_time("step$step_n. $step_hash->{$step_n} done, escaped time: $escaped_time.\n");
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

#############################################
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
        return ($time) ;
}


#############################################################################################################
sub USAGE {
	my $usage=<<"USAGE";
----------------------------------------------------------------------------------------------------------------------------
   Program: $title
   Version: $ver

     Usage:
            --indir      <DIR>  bam file path
            --type      <FILE>  the type Chr,such as : Human.B38  Human.hg19 Mouse.B38  Rat.B5.0
            --od        <DIR>   analysis output directory
            --script_cfg <FILE> script.cfg,contain software path and queue ;

            --h                 help documents

   Example:
  perl $Script --indir ../result/Basic_Analysis/Tophat_Cufflinks/Tophat/ --type B38 --od Analysis/ --script_cfg script.config
/**********/**********/**********/**********/  

                       __,,,,_ 
       _ ___.--'''`--''// ,-_ `-. 
   \`"' ' |  \  \ \\/ / // / ,-  `,_ 
  /'`  \   |  Y  | \|/ / // / -.,__ `-. 
 /<"\    \ \  |  | ||/ // | \/    |`-._`-._ 
/  _.-.  .-\,___|  _-| / \ \/|_/  |    `-._) 
`-'  f/ |       / __/ \__  / |__/ | 
     `-'       |  -| -|\__ \  |-' | 
            __/   /__,-'    ),'  |' 
           ((__.-'((____..-' \__,' 
/**********/**********/**********/**********/  

USAGE
	print $usage;
	exit;
}
