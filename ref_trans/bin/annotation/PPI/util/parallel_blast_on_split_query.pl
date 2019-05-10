#!/usr/bin/perl -w
use strict;
use warnings;
use Getopt::Long;
use Data::Dumper;
use FindBin qw($Bin $Script);
use File::Basename qw(basename dirname);
use Term::ANSIColor qw(:constants);
use newPerlBase;
$Term::ANSIColor::AUTORESET=1;
my $BEGIN_TIME=time();
my $version="1.0.0";
#######################################################################################
#
#######################################################################################
# ------------------------------------------------------------------
# GetOptions
# ------------------------------------------------------------------
my ($query_file, $database_file, $index, $odir, $cutnum, $evalue, $xml, $notab, $verbose,$queue);

GetOptions(
    "help|?"    => \&USAGE,
    "query=s"   => \$query_file,
    "database=s"=> \$database_file,
    "index=s"   => \$index,
    "odir=s"    => \$odir,
    "cutnum=i"  => \$cutnum,
    "evalue=f"  => \$evalue,
    "xml"       => \$xml,
    "notab"     => \$notab,
    "verbose"   => \$verbose,
    "queue=s"   => \$queue,
    ) or &USAGE;
&USAGE unless ($query_file and $database_file);

# ------------------------------------------------------------------
# init 
# ------------------------------------------------------------------
$index  ||= basename($query_file);
$odir   ||= './';
$cutnum ||= 200;
$evalue ||= 1e-5;

system "mkdir -p $odir" unless (-d $odir);
$query_file    = &ABSOLUTE_DIR($query_file);
$database_file = &ABSOLUTE_DIR($database_file);
$odir          = &ABSOLUTE_DIR($odir);
#$index =~s/\.fa(sta)?$//i;
my $work_sh = "$odir/work_sh";
system "mkdir -p $work_sh" unless (-d $work_sh);

# utils 
my $blast = ${&readconf("$Bin/../../../../config/db_file.cfg")}{blast};  # 2015-05-26
my $formatdb = "$blast/makeblastdb";  # 2015-05-26
my $merge_blast_xml = "$Bin/merge_blast_xml.py";    # 2015-05-26
my $blast_parser = "$Bin/blast_parser.pl";          # 2015-05-26
my $cmd;

&timeLog("$Script start...");
# ------------------------------------------------------------------
# split query file and formatdb if needed
# ------------------------------------------------------------------
my $program;
my $blastdb;
my $query_biotype = fasta_format_check($query_file);
my $sbjct_biotype = fasta_format_check($database_file);

# make sure the BLAST program to be used 
if ($query_biotype eq 'nucleotide' && $sbjct_biotype eq 'aminoacide') {
    $program = 'blastx';
    #$blastdb = "$database_file.psq";
	$blastdb = (glob "$database_file*.psq")[0];
} elsif ($query_biotype eq 'aminoacide' && $sbjct_biotype eq 'aminoacide') {
    $program = 'blastp';
    $blastdb = (glob "$database_file*.psq")[0];
} elsif ($query_biotype eq 'nucleotide' && $sbjct_biotype eq 'nucleotide') {
    $program = 'blastn';
    $blastdb = (glob "$database_file*.nsq")[0];
} else {
    die "ERROR: only suit for BLAST program blastx, blastp or blastn!\n";
}

# formatdb if needed 
if (not defined $blastdb or (defined $blastdb && !-f $blastdb)) {
    mkdir "$odir/blastdb" unless (-d "$odir/blastdb");
    $cmd = "cd $odir/blastdb/ && ln -s $database_file* ./ && " unless (-s "$odir/blastdb/".basename($database_file));
    $database_file = "$odir/blastdb/".basename($database_file);

    # nt or pep 
   #if ($sbjct_biotype eq 'aminoacide') {
   #    $cmd.= "$formatdb -i $database_file -p T ";
   #    $blastdb = (glob $database_file.'*.psi')[0];
   #} elsif ($sbjct_biotype eq 'nucleotide') {
   #    $cmd.= "$formatdb -i $database_file -p F ";
   #    $blastdb = (glob $database_file.'*.nsi')[0];
   #}

    &run_or_die($cmd) if (not defined $blastdb or (defined $blastdb && !-f $blastdb));
}

# split query file
my %query;
my $query_dir = "$odir/query_dir";
mkdir $query_dir unless (-d $query_dir);

&load_fasta($query_file,\%query);
&cut_fasta(\%query, $query_dir, $cutnum, $index);

# ------------------------------------------------------------------
# BLAST 
# ------------------------------------------------------------------
# blast 
my $blast_shell_file = "$work_sh/$index.blast.sh";
my @subfiles = glob "$query_dir/$index*.fa";
mkdir "$odir/aln_subdir/" unless (-d "$odir/aln_subdir/");

# creat shell file
open OUT,">$blast_shell_file" or die $!;
foreach my $subfile (@subfiles) {
    my $name = basename($subfile);
    if ($xml) {
		print OUT "$blast/$program -task blastx-fast -max_target_seqs 100 -outfmt 5 -evalue $evalue -db $database_file -query $subfile -num_threads 2 -out $odir/aln_subdir/$name.blast.xml &&\n";
        #print OUT "$blast/$program -b 50 -v 50 -p $program -e $evalue -F F -d $database_file -i $subfile -m 7 -a 4 -o $odir/aln_subdir/$name.blast.xml && \n";
    } else {
		print OUT "$blast/$program -task blastx-fast -max_target_seqs 100  -evalue $evalue -db $database_file -query $subfile -num_threads 2 -out $odir/aln_subdir/$name.blast &&\n";
        #print OUT "$blastall -b 50 -v 50 -p $program -e $evalue -F F -d $database_file -i $subfile -m 0 -a 4 -o $odir/aln_subdir/$name.blast && \n";
    }
}
close OUT;

#run the shell file
&qsubOrDie("$blast_shell_file",$queue,80,"6G");

# ------------------------------------------------------------------
# merge BLAST result 
# parse the BLAST result and convert to tabular format 
# ------------------------------------------------------------------

if ($xml) {
    $cmd = "python $merge_blast_xml $odir/$index.blast.xml $odir/aln_subdir/*.blast.xml ";
} else {
    $cmd = "cat $odir/aln_subdir/*.blast > $odir/$index.blast ";
}

&runOrDie($cmd);

# ------------------------------------------------------------------
# convert alignment result to tabular format
# ------------------------------------------------------------------
# convert to tabular format
if ($xml) {
    $cmd = "perl $blast_parser -eval $evalue -tophit 1  -m 7 -topmatch 1 $odir/$index.blast.xml > $odir/$index.blast.tab.best && ";
    $cmd.= "perl $blast_parser -eval $evalue -tophit 50 -m 7 -topmatch 1 $odir/$index.blast.xml > $odir/$index.blast.tab.best ";
} else {
    $cmd = "perl $blast_parser -eval $evalue -tophit 1  -m 0 -topmatch 1 $odir/$index.blast > $odir/$index.blast.tab.best && ";
    $cmd.= "perl $blast_parser -eval $evalue -tophit 50 -m 0 -topmatch 1 $odir/$index.blast > $odir/$index.blast.tab ";
}

&runOrDie($cmd) unless ($notab);

# remove intermediate file 
unless ($verbose) {
    system "rm -r $query_dir $odir/aln_subdir ";
    system "rm -r $odir/blastdb " if (-d "$odir/blastdb");
}

#######################################################################################
my $elapse_time = (time()-$BEGIN_TIME)."s";
&timeLog("$Script done. Total elapsed time: $elapse_time.");
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

#############################################################################################################
sub load_fasta {
    &log_current_time("load FASTA file.");
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

#############################################################################################################
sub cut_fasta {
    &log_current_time("cut query sequence.");
    my ($fa,$od,$cut,$name) = @_;
    my %seq=%$fa;
    my @aa=sort(keys %seq);
    my $index=0;
    my $filenum = int(@aa/$cut) + 1;
    my $decimals=int(log($filenum)/log(10))+1;

    LAB: for (my $i=1;;) {
        my $num=0;
        $i=sprintf "%0${decimals}d",$i;
        open OUT,">$od/$name.$i.fa" || die $!;
        for ($index..$#aa) {
            $index++;
            if ($num<$cut) {
                print OUT ">$aa[$_]\n$seq{$aa[$_]}\n";
                $num++;
            }
            if ($num>=$cut) {
                $num=0;
                $i++;
                close OUT;
                if ($index==$#aa+1) {
                    last;
                } else {
                    next LAB;
                }
            }
        }

        close OUT if ($num);
        last;
    }
}

#############################################################################################################
sub Cut_shell_qsub {#Cut shell for qsub 1000 line one file
    # &Cut_shell_qsub($shell,$cpu,$vf,$queue);
    my $shell = shift;
    my $cpu = shift;
    my $vf = shift;
    my $queue = shift;
    my $line = system "wc -l $shell";
    my $notename=`hostname`;chomp $notename;

    if ($line<=1000) {
		die "no sub function\n";
        #system "sh /share/nas2/genome/bmksoft/tool/qsub_sge_plus/v1.0/qsub_sge.plus.sh --queue $queue --maxproc $cpu --resource vf=$vf --independent --reqsub $shell ";
    } else {
        my @div = glob "$shell.div*.sh";
        foreach (@div) {
            system "rm $_" if (-e $_);
        }
        my $div_index=1;
        my $line_num=1;

        open IN,"$shell" or die $!;
        while (<IN>) {
            chomp;
            open OUT,">>$shell.div.$div_index.sh" or die $!;

            if ($line_num<1000) {
                print OUT "$_\n";
                $line_num++;
            } else {
                print OUT "$_\n";
                $div_index++;
                $line_num=0;
                close OUT;
            }
        }
        close OUT if ($line_num!=0);

        @div=glob "$shell.div*.sh";
        foreach my $div_file (@div) {
			die"no sub function\n";
           # system "sh /share/nas2/genome/bmksoft/tool/qsub_sge_plus/v1.0/qsub_sge.plus.sh --queue $queue --maxproc $cpu --resource vf=$vf --reqsub $div_file ";
        }
    }
}

#############################################################################################################
#&run_or_die($cmd);
sub run_or_die() {
    my ($cmd) = @_ ;
    &log_current_time($cmd);
    my $flag = system($cmd);
#    my $flag = 0;

    if ($flag){
        &log_current_time("Error: command fail: $cmd");
        exit(1);
    }
}

#############################################################################################################
sub fasta_format_check {
    my $file = shift;
    my $context = `grep -v "\^#" $file|head -n 2 `; chomp $context;
    my ($id, $seq) = (split /\n/,$context);
    my $biotype;

    if ($id=~/^>/) {
        $biotype = ($seq =~/[^ATCGUN]/i) ? 'aminoacide' : 'nucleotide'; # error when nucleotide with degeneracy 
    } else {
        die "ERROR: this may be not a FASTA format file, please check $file !\n";
    }

    return $biotype;
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
	my $usage = <<"USAGE";
 ProgramName: $Script
     Version: $version
     Contact: Simon Young <yangxh\@biomarker.com.cn> 
Program Date: 2015-05-26
      Modify: 
 Description: This script is used to BLASTA on split query parallelly. Only suit for BLAST program blastx, blastp or blastn.
              
       Usage: 
        Options:
        --query     <FILE>  query file, FASTA format, nucleotide or amino acid sequence, required
        --database  <FILE>  database file, FASTA format, nucleotide or amino acid sequence, required

        --index     <STR>   prefix of output files, optional, default basename of query file
        --odir      <DIR>   output directory, optional                                                  [./]
        --cutnum    <INT>   queried sequence number per splited single file                             [200]
        --evalue    <FLOAT> expectation value, optional                                                 [1e-5]
        --xml               report BLAST alignment result as XML format, default BLAST standard format
        --notab             not convert alignment result to tabular format
        --verb              verbose, save intermediate result

        Examples:
            perl $Script --query Maize.Known.longest_transcript.fa --database protein.sequences.v10.fa 
            perl $Script --query Bee.Unigene.fa --database protein.sequences.v10.fa --odir PPI/align_STRING --cutnum 100 --evalue 1e-10

USAGE
	print $usage;
	exit;
}


