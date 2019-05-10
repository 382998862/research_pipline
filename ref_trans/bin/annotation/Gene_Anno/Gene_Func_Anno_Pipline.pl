#!/usr/bin/perl -w
use strict;
use warnings;
use Getopt::Long;
use Data::Dumper;
use FindBin qw($Bin $Script);
use threads;
use Thread::Semaphore;
use File::Basename qw(basename dirname);
my $programe_dir=basename($0);
my $path=dirname($0);
use Term::ANSIColor qw(:constants);
$Term::ANSIColor::AUTORESET=1;
my $version="1.4.0";
my $BEGIN_TIME=time();

use newPerlBase;
my %config=%{readconf("$Bin/db_file.cfg")}; 

#######################################################################################
# ------------------------------------------------------------------
# GetOptions, Check format, Output info
# ------------------------------------------------------------------
my ($cfg,$odir);
my ($nr,$nt,$Swissprot,$TrEMBL,$GO,$eggNOG,$Pfam,$Kegg,$all,$Cog,$Kog);
my ($Query,$Blastp,$Blast_cpu,$Blast_e,$blast_cut,$hmmscan_cpu);
my ($sh_dir,$Result_dir,$Tab_dir,$div_dir,$queue);
my $step;

GetOptions(
		"cfg:s"=>\$cfg,
		"all"=>\$all,
		"nr"=>\$nr,
		"nt"=>\$nt,
                "kog"=>\$Kog,
                "cog"=>\$Cog,
		"swissprot"=>\$Swissprot,
		"trembl"=>\$TrEMBL,
		"eggNOG"=>\$eggNOG,
		"pfam"=>\$Pfam,
		"kegg"=>\$Kegg,
		"GO"=>\$GO,
		"queue:s"=>\$queue,
		"od:s"=>\$odir,
                "step=i" =>\$step,
		"help"=>\&USAGE
	) or &USAGE;
&USAGE unless ($cfg and $odir) ;

# ------------------------------------------------------------------
# load arguments form config file & init parameters.
# ------------------------------------------------------------------
&MKDIR($odir);
$odir=&ABSOLUTE_DIR($odir);
$cfg=&ABSOLUTE_DIR($cfg);

# load arguments form config file
my %CFG;
&LOAD_PARA($cfg,\%CFG);
$Query=$CFG{mRNA};
$Blast_cpu=$CFG{blast_cpu};
$hmmscan_cpu=$CFG{hmmscan_cpu};
$Blast_e=$CFG{blast_e};
$blast_cut=$CFG{blast_cut};
# init parameters
my $Q_name=basename $Query;
$step ||= 1;
$queue||="general.q";
# deal with illegal parameters
if (( !defined $nr)&& ( !defined $nt) && ( !defined $Cog) && ( !defined $Kog) && ( !defined $eggNOG) && ( !defined $Pfam) && ( !defined $Swissprot) && ( !defined $TrEMBL) && ( !defined $Kegg) && (! defined $GO) && (!defined $all) ) {
	print BOLD RED "Error: You must choose database to align.\n";
	exit;
}
if (-f $Query) {
    print "[".&date_time_format(time())."] make sure your query FASTA file: $Query\n";
} else {
    print BOLD RED "Error: Can't find your query file to annotaion.\n";
    exit;
}
unless ($step>=1 and $step<=4) {
	print BOLD RED "Error: start step must belong to {1,2,3,4}.\n";
	exit;
}

#######################################################################################
print "\n[".&date_time_format($BEGIN_TIME)."] $Script start ...\n";
#######################################################################################
# ------------------------------------------------------------------
# preparation
# ------------------------------------------------------------------
# make basic directories
&MKDIR("$odir/mid");
my $mid_dir="$odir/mid";
&MKDIR("$odir/Result");
$Result_dir="$odir/Result";
&MKDIR("$odir/work_sh");
$sh_dir="$odir/work_sh";
chomp(my $who=`whoami`);
#system("cp -f /home/xugl/.kobasrc ~/") unless -e "/home/$who/.kobasrc";
#system("export LD_LIBRARY_PATH=/share/nas2/genome/biosoft/R/3.1.1/lib64/R/lib/:\$LD_LIBRARY_PATH");
# ------------------------------------------------------------------
# cut query FASTA file
# ------------------------------------------------------------------
if ($step == 1) {
	unless(-e "$mid_dir/cut_fa.Check"){
		print "[".&date_time_format(time())."] cut query FASTA file:\n";
		print "[".&date_time_format(time())."] perl $Bin/bin/cut_fa_file_to_dir.pl -mRNA $Query -od $mid_dir -cut $blast_cut\n\n";
		my $cmd_res=system "perl $Bin/bin/cut_fa_file_to_dir.pl -mRNA $Query -od $mid_dir -cut $blast_cut";
		if($cmd_res){
			print STDERR "Err:Cut fa failed!\n";
		}
		else{
			system "touch $mid_dir/cut_fa.Check";
		}
	}
	else{
		print "[".&date_time_format(time())."] cut query FASTA file has been finished:\n\n";
	}
    $step = 2;
}

# ------------------------------------------------------------------
# align with databases and alignment result format.
# ------------------------------------------------------------------
if ($step == 2) {
    print "[".&date_time_format(time())."] now align with databases:\n";
    my $SH="";
	my @work_sh;
    if (defined $nr||defined $all) {
		$SH="$odir/work_sh/NR_Anno.sh";
		open (SH,">$SH") or die $!;
        print "        nr: $CFG{nr} \n";
        print SH "perl $Bin/bin/Gene_Nr_Anno.pl -id $mid_dir -Database $CFG{nr} -od $odir -Blast_cpu $Blast_cpu -Blast_e $Blast_e -queue  $queue\n";
		close SH;	
		push @work_sh,$SH;
    }
	if (defined $Kegg||defined $all) {
        print "      KEGG: $CFG{Kegg} \n";
		$SH="$odir/work_sh/KEGG_Anno.sh";
		open (SH,">$SH") or die $!;
        print SH "perl $Bin/bin/Gene_KEGG_Anno.pl -id $mid_dir -Database $CFG{Kegg} -od $odir -Blast_cpu $Blast_cpu -Blast_e $Blast_e  -queue  $queue \n";
		close SH;
		push @work_sh,$SH;
    }
	if (defined $nt||defined $all) {
        print "        nt: $CFG{nt} \n";
		$SH="$odir/work_sh/NT_Anno.sh";
		open (SH,">$SH") or die $!;
        print SH "perl $Bin/bin/Gene_Nt_Anno.pl -id $mid_dir -Database $CFG{nt} -od $odir -Blast_cpu $Blast_cpu -Blast_e $Blast_e -queue  $queue \n";
		close SH;
		push @work_sh,$SH;
    }
    if (defined $eggNOG ||defined $all) {
        print "       eggNOG: $CFG{eggNOG} \n";
		$SH="$odir/work_sh/eggNOG_Anno.sh";
		open (SH,">$SH") or die $!;
        print SH "perl $Bin/bin/Gene_eggNOG_Anno.pl  -id $mid_dir -Database $CFG{eggNOG} -od $odir -cpu $Blast_cpu -queue  $queue -Blast_e $Blast_e\n";
		close SH;
		push @work_sh,$SH;
    }
    
        if (defined $Cog||defined $all) {

                 print "       COG: $CFG{Cog} \n";
                $SH="$odir/work_sh/COG_Anno.sh";
 	        open (SH,">$SH") or die $!;
         print SH "perl $Bin/bin/Gene_Cog_Anno.pl -id $mid_dir -Database $CFG{Cog} -od $odir -Blast_cpu $Blast_cpu -Blast_e $Blast_e -queue   $queue\n";
		close SH;
		push @work_sh,$SH;
    }
    if (defined $Kog||defined $all) {
        print "       KOG: $CFG{Kog} \n";
                $SH="$odir/work_sh/KOG_Anno.sh";
                open (SH,">$SH") or die $!;
                print SH "perl $Bin/bin/Gene_Kog_Anno.pl -id $mid_dir -Database $CFG{Kog} -od $odir -Blast_cpu $Blast_cpu -Blast_e $Blast_e -queue   $queue \n";
                close SH;
	push @work_sh,$SH;
    }
   
  
    if (defined $Pfam||defined $all) {
        print "      Pfam: $CFG{Pfam} \n";
		$SH="$odir/work_sh/Pfam_Anno.sh";
		open (SH,">$SH") or die $!;
        print SH "perl $Bin/bin/Gene_Pfam_Anno.pl --idir $mid_dir --pfam $CFG{Pfam} --odir $odir --cpu $hmmscan_cpu -queue  $queue -lib_type $CFG{Lib_type}\n";
		close SH;
		push @work_sh,$SH;
    }
    if (defined $Swissprot||defined $all) {
        print "Swiss-Prot: $CFG{Swissprot} \n";
		$SH="$odir/work_sh/Swissprot_Anno.sh";
		open (SH,">$SH") or die $!;
        print SH "perl $Bin/bin/Gene_Swissprot_Anno.pl -id $mid_dir -Database $CFG{Swissprot} -od $odir -Blast_cpu $Blast_cpu -Blast_e $Blast_e -queue  $queue\n";
		close SH;
		push @work_sh,$SH;
    }
    if (defined $TrEMBL||defined $all) {
        print "    TrEMBL: $CFG{TrEMBL} \n";
		$SH="$odir/work_sh/TrEMBL_Anno.sh";
		open (SH,">$SH") or die $!;
        print SH "perl $Bin/bin/Gene_TrEMBL_Anno.pl -id $mid_dir -Database $CFG{TrEMBL} -od $odir -Blast_cpu $Blast_cpu -Blast_e $Blast_e -queue  $queue\n";
		close SH;
		push @work_sh,$SH;
    }
    &Shell_run (\@work_sh);
    $step = 3;
}

# ------------------------------------------------------------------
# Blast2GO based on nr database blast result
# ------------------------------------------------------------------
if ($step == 3 && ($GO||$all)) {
    print "[".&date_time_format(time())."] start Blast2GO based on nr blast result:\n";
    print "[".&date_time_format(time())."] perl $Bin/bin/Blast2GO.pl -id $odir/Nr_Dir/ -od $Result_dir -sh_dir $sh_dir -k $Q_name -queue  $queue\n\n";
    system  "perl $Bin/bin/Blast2GO.pl -id $odir/Nr_Dir/ -od $Result_dir -sh_dir $sh_dir -k $Q_name -queue  $queue";
    $step = 4;
}

# ------------------------------------------------------------------
# annotation result tackle, plot and stat
# ------------------------------------------------------------------
if ($step == 4) {
    ## result integrate
    if((defined $eggNOG || defined $Cog||defined $Kog || defined $nr || defined $nt || defined $Swissprot || defined $TrEMBL || defined $Kegg) || defined $all){
        print "[".&date_time_format(time())."] Anno Integrate:\n";
        print "[".&date_time_format(time())."] perl $Bin/bin/anno_integrate.pl -gene $Query -id $Result_dir -od $Result_dir -key $Q_name \n\n";
        system "perl $Bin/bin/anno_integrate.pl -gene $Query -id $Result_dir -od $Result_dir -key $Q_name ";
    }

    ## GO histogram plot
    if ($GO||$all) {
        print "[".&date_time_format(time())."] Draw GO Graph:\n";
        print "[".&date_time_format(time())."] perl $Bin/bin/draw_GO_graph.pl -i $Result_dir/All_Database_annotation.xls -k $Q_name -od $Result_dir \n\n";
        system "perl $Bin/bin/draw_GO_graph.pl -i $Result_dir/All_Database_annotation.xls -k $Q_name -od $Result_dir ";
    }

    ## nr pie plot
    if ($nr||$all) {
		print "[".&date_time_format(time())."] perl $Bin/bin/util/nr_pie_stat.pl -i $Result_dir/$Q_name.nr.anno.txt -o $Result_dir/$Q_name.nr.lib.stat -m 10\n\n";
        system "perl $Bin/bin/util/nr_pie_stat.pl -i $Result_dir/$Q_name.nr.anno.txt -o $Result_dir/$Q_name.nr.lib.stat -m 10 ";
    }
}

#######################################################################################
print "\n[".&date_time_format(time())."] $Script done. Total elapsed time: ".(time()-$BEGIN_TIME)."s\n";
#######################################################################################
# ------------------------------------------------------------------
# sub function
# ------------------------------------------------------------------
################################################################################################################
sub LOAD_PARA {
	my $para_file= shift;
	my $para= shift;

    # switch node in which store databases to search.
    # add node to configs file when needed.
    my %db_node;
    my $db_idx = 0;
    open(NODE,"$Bin/db_file.cfg") or die;
    while (<NODE>) {
        next unless (/^db_node/);
        my $node = (split /\s+/,$_)[1];
        $db_node{$db_idx++} = $node;
    }
    close NODE;
    my $random = int(rand(scalar keys %db_node));

	my $error_status = 0;
	open IN,$para_file || die "fail open: $para_file";
	while (<IN>) {
		chomp;
		s/^\s+//;s/\s+$//;s/\r$//;
		next if (/^$/ or /^\#/) ;
		my ($para_key,$para_value) = split(/\s+/,$_);
        $para_value =~s/\/lustre|\/share/$db_node{$random}/ unless ($para_key eq 'mRNA');
        $para->{$para_key} = $para_value;

		if (!defined $para_value) {
			warn "Non-exist: $para_key\n";
			$error_status = 1;
		}
        if ($para_key eq "nr" or $para_key eq "Cog" or $para_key eq "Kog" or $para_key eq "Pfam" or $para_key eq "Swissprot" or $para_key eq "Kegg" or $para_key eq "nt" or $para_key eq "TrEMBL") {
            unless (-f $para_value) {
                warn "Error: Can't find $para_key database: $para_value.\n";
                $error_status = 1;
            }
        }
	}
	close IN;
	die "\nExit due to error Undefined parameter\n" if ($error_status) ;
}

sub LOAD_SEQ {
	my ($fa,$info) = @_;

	open IN,"$fa" || die $!;
	$/='>';
	while (<IN>) {
		chomp;
		s/^\s+//;s/\s+$//;s/\r+$//;
		next if (/^$/ || /^\#/);
		my ($head,$seq)=split/\n+/,$_,2;
		my $id=(split/\s+/,$head)[0];
		$info->{$id}=$seq;
	}
	$/="\n";
	close IN;
}

sub Shell_run {
	my $sh = shift;
	my $cpu=2;
	my $semaphore=Thread::Semaphore->new($cpu);
	foreach my $work_sh(@{$sh}){
		$semaphore->down();
		my $thread=threads->new(sub{
			my $time=time();
			unless(-e "$work_sh.finish"){
				my $cmd_res=system "sh $work_sh>$work_sh.log 2>&1";
				if($cmd_res){
					die "$work_sh doesnot run successfully,$!\n";
				}
				else{
					system"touch $work_sh.finish";
				}
			}
			my $runTime=time()-$time;
			print "$work_sh run time:$runTime"."s\n";
			$semaphore->up();
		},$work_sh);
		my $tid=$thread->tid();
		print "$work_sh thread id is $tid\n";
		foreach my $thread (threads->list(threads::joinable)){
			$thread->join();
			$tid=$thread->tid();
			print "The thread id $tid has joined,";
			if(my $error=$thread->error()){
				print "but $error";
			}
			print "\n";
		}
	}
	for my $thread (threads->list(threads::all)){
		$thread->join();
		my $tid=$thread->tid();
		print "The thread id $tid has joined";
		if(my $error=$thread->error()){
			print "but $error";
		}
		print "\n";
	}
}

sub cut_str {
	my $string = shift;
	my @str = split /\s+/,$string;
	if (@str > 2) {
		return "$str[0] $str[1]"
	}else{
		return $string;
	}
}

sub parse_config { # load config file
	my $config_file= shift;
	my $DataBase= shift;
	
	my $error_status = 0;
	open IN,$config_file || die "fail open: $config_file";
	while (<IN>) {
		chomp;
		s/\s+//;s/\s+$//;s/\r$//;
		next if(/$/ or /\#/);
		my ($software_name,$software_address) = split(/\s+/,$_);
		$DataBase->{$software_name} = $software_address;
		if (! -e $software_address){
			warn "Non-exist:  $software_name  $software_address\n"; 
			$error_status = 1;
		}
	}
	close IN;
	die "\nExit due to error of software configuration\n" if($error_status);
}

sub Runtime { # &Runtime($BEGIN);
	my ($t1)=@_;
	my $t=time()-$t1;
	print "Total elapsed time: ${t}s\n";
}

sub date_time_format {#Time calculation subroutine
    my ($sec, $min, $hour, $day, $mon, $year, $wday, $yday, $isdst)=localtime(time());
    return sprintf("%4d-%02d-%02d %02d:%02d:%02d", $year+1900, $mon+1, $day, $hour, $min, $sec);
}

sub MKDIR { # &MKDIR($out_dir);
	my ($dir)=@_;
	rmdir($dir) if(-d $dir);
	mkdir($dir) if(!-d $dir);
}

sub ABSOLUTE_DIR { #$pavfile=&ABSOLUTE_DIR($pavfile);
	my $cur_dir=`pwd`;chomp($cur_dir);
	my ($in)=@_;
	my $return="";

	if (-f $in) {
		my $dir=dirname($in);
		my $file=basename($in);
		chdir $dir;$dir=`pwd`;chomp $dir;
		$return="$dir/$file";
	} elsif(-d $in) {
		chdir $in;$return=`pwd`;chomp $return;
	} else {
		warn "Warning just for file or dir in $in\n";
		exit;
	}

	chdir $cur_dir;
	return $return;
}

sub USAGE {
    my $usage =<<"_USAGE_";
#-------------------------------------------------------------------------------------------------
    Program: $Script
    Version: $version
     Writer: Meng Fei <mengf\@biomarker.com.cn>
       Data: 2012-00-00
   Modifier: Simon Young <simonyoung8824\@gmail.com>
       Data: 2014-09-28
Description: alignment with fuctional database and annotate genes. v1.1.0-based,
             1) add KOG and Pfam database annotation, additionally predict CDS;
             2) deal with disturbed work_sh output and excessive qsub;
             3) fix several Blast2GO empty annotation bug;
             4) fix excel overflow bug by spliting files with more than 65,000 lines.

      Usage:
             --all             search against all database
             
               --nr            search against nr database
               --swissprot     search against Swiss-Prot database
               --eggNOG        search against eggNOG database
               --kegg          search against KEGG database
               --pfam          search against Pfam database
               --GO            run Blast2GO start GO Annotation (--nr is required)
               --nt            search against nt database
               --trembl        search against TrEMBL database
	       --cog           search against COG database
	       --kog	       search against KOG database
             --step            program start step, (1|2|3|4)
               1               start from the beginning, default
               2               start from alignment with databases and alignment result format
               3               start from Blast2GO
               4               start from result tackle, plot and stat

             --cfg             query and database config
             --queue        queue for qsub jobs
             --od              output dir
             --help            show help document

	Example:
	  perl $Script --all --cfg Gene_Func_Anno_Pipline.cfg --od ./

#-------------------------------------------------------------------------------------------------
_USAGE_
    print $usage;
    exit;
}
