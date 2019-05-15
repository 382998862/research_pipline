#!/usr/bin/perl

use warnings;
use strict;
use Getopt::Long;
use File::Path;
use File::Basename;
use POSIX qw(floor);
use threads;
use Thread::Semaphore;

####################################### USAGE ####################################################
my $usage =
"
This script is a wrapper of TAPIR for prediction of miRNA targets.
Options:

Preprocessing/mapping:
-i <string>     miRNA file in fasta format
-t <string>     target transcriptome file in fasta format
-o <string>     output file
-s <float>      score cutoff in TAPIR. [default: 4]        
-r <float>      mfe ratio cutoff. [Default: 0.65]
-b              tabular report the result.
-m <int>        number of threads to lanuch. [Default: 1]

Other:
-q              Quiet mode .. no log/progress information to STDERR
-h              Print help message and quit 

Example of use:

$0 -i ath_mature.fa -t TAIR10.fa -m 4 -q -b -o target_prediction.csv
";

# if there are no arguments, return the help message and quit
unless($ARGV[0]) {
    die "No arguments specified!\n$usage";
}
###################################### INPUT #######################################################
my ($opt_i,$opt_h,$opt_q,$opt_t,$opt_s,$opt_r,$opt_m,$opt_b,$opt_o);

GetOptions(
    'input|i=s'         => \$opt_i,
    'target|t=s'        => \$opt_t,
    'output|o=s'        => \$opt_o,
    'score|s=f'         => \$opt_s,
    'threads|m=i'       => \$opt_m,
    'ratio|r=f'         => \$opt_r,
    'tabular|b!'        => \$opt_b,
    'help|h!'           => \$opt_h,
    'quiet|q!'          => \$opt_q,
);

#################################### GLOBAL VARIABLES ################################################
my $threads=1;
$threads=$opt_m if $opt_m;

my $score=4;
$score=$opt_s if $opt_s;

my $ratio=0.65;
$ratio=$opt_r if $opt_r;

my $mir_file=$opt_i if $opt_i;
my $target_file=$opt_t if $opt_t;

my $max=0; #maximum length of target, used in tapir_hybrid
####################################### MAIN ########################################################
# Check options
check_options();
unless($opt_q){
    print STDERR "########Predition target of miRNA using TAPIR########\n";
}

# Check dependencies
my $tapir_check = check_tapir();
unless($tapir_check) {
    die "FAIL: please install TAPIR software.\n\n$usage";
}

# make tmp directory
my $dir="tapir_tmp";
unless(-d $dir){
    mkdir($dir);
}

# running TAPIR...
unless($opt_q) {
    print STDERR "Running TAPIR for predicting...\n"; 
}
my $semaphore=new Thread::Semaphore($threads);
my $ref_all=patition_for_prediction($target_file);
foreach(@$ref_all){
    $semaphore->down();
    my $thr = threads->new(\&tapir, \$mir_file,\$_);
    $thr->detach();
}
waitquit(\$threads);

foreach(@$ref_all){
    system("cat $_.result >>$dir/tapir_prediction.result");
}
if($opt_b){
    convert_tab("$dir/tapir_prediction.result",$opt_o);
}
else{
    system("mv $dir/tapir_prediction.result $opt_o");
}
# Parse the result of tapir_hybrid
# foreach(@$ref_all){
    # my $file=$_.".result";
    # hybrid_parser($file);
# }
rmtree($dir);
unless($opt_q) {
    print STDERR "Prediction down\n"; 
}
exit;
######################################### SUBS #####################################################
sub check_options{
    if($opt_h) {
        die "$usage";
    }
    unless($opt_i){
        die "FATAL: please give a miRNA file in fasta format\n";
    }
    unless($opt_t){
        die "FATAL: please give a miRNA file in fasta format\n";
    }
    unless($opt_o){
        die "FATAL: please specify the output file\n";
    }
    unless(-f $opt_i){
        die "FATAL: can not find file $opt_i\n";
    }
    unless(-f $opt_t){
        die "FATAL: can not find file $opt_t\n";
    }
    unless(($ratio > 0) and ($ratio <= 1)) {
        die "FATAL: option r must be greater than 0 and less than or equal to 1\n\n$usage";
    }
    unless($threads=~/^\d+$/){
        die "FATAL: threads number must be an integer\n\n$usage";
    }
    
    ## check number of cores on the system and threads to be used
    my $cores=`grep -ic ^processor /proc/cpuinfo`;
    if($cores !~ /^\d+$/){
        $cores=`sysctl -n hw.physicalcpu`;
        if($cores !~ /^\d+$/){
            $cores=`sysctl -n hw.logicalcpu`;
        }
    }
    if($cores !~ /^\d+$/){
        $cores=1;
    }
    
    if($threads > $cores){ 
        $threads=$cores;
    }
}

sub check_tapir {
    #check tapir_hybrid
    (open(TRH, "trh -h |")) || return 0;
    <TRH>;
    my $tapir = <TRH>;
    close TRH;
    (open(PAR, "hybrid_parser -h |")) || return 0;
    my $parser = <PAR>;
    close PAR;
    if($tapir =~ /^Usage/ and $parser =~ /^Usage/) {
        return 1;
    }
    else {
        return 0;
    }
    # should never be here
    return 0;
}

sub tapir{
    my($mir,$target)=@_;
    my $cmd="trh -u 5 -b 10 -d 2.38,0.19 -c -m $max";
    $cmd.=" -q ".$$mir;
    $cmd.=" -t $$target ";
    #$cmd.=" >$$target.result";
    $cmd.=" |hybrid_parser --score $score --mfe_ratio $ratio >$$target.result";
    system($cmd);
    $semaphore->up();
}
sub hybrid_parser{
    my $file=shift;
    my $cmd="cat $file | ";
    $cmd.="hybrid_parser --score $score --mfe_ratio $ratio";
    $cmd.=" >>$dir/tapir_prediction.result";
    system($cmd);
}

sub patition_for_prediction{
    my $file=shift;
    open IN,$file || die "Can not open $file";
    my $seq;
    my $tmp="";
    my @all;
    my @pat_files;
    while(<IN>){
        if(/^>/){
            if($tmp){
                push(@all,$tmp);
                $seq=~s/\s+//g;
                my $len=length($seq);
                $max=$len if $len>$max;
            }
            $tmp="";
            $seq="";
        }
        else{
            $seq.=$_;
        }
        $tmp.=$_;
    }
    push(@all,$tmp);
    $seq=~s/\s+//g;
    my $len=length($seq);
    $max=$len if $len>$max;
    close IN;
    my $num=floor(@all/$threads);
    $file=basename($file);
    if($num==0){
        foreach(1..@all){
            open TMP,">$dir/${file}"."_tmp".$_ || die "Cannot create ${file}_tmp$_";
            push @pat_files,"$dir/${file}"."_tmp".$_ ;
            print TMP $all[$_-1];
            close TMP;
        }
    }
    else{
        my $begin=0;
        foreach(1..$threads-1){
            open TMP,">$dir/${file}"."_tmp".$_ || die "Cannot create ${file}_tmp$_";
            push @pat_files,"$dir/${file}"."_tmp".$_ ;
            my $end=$_*$num-1;
            print TMP (@all[$begin..$end]);
            $begin=$end+1;
            close TMP;
        }
        open TMP,">$dir/${file}"."_tmp".$threads  || die "Cannot create ${file}_tmp$threads";
        push @pat_files,"$dir/${file}"."_tmp".$threads;
        print TMP (@all[$begin..$#all]);
        close TMP;
    }
    return \@pat_files;
}
sub waitquit{
    my $r_threads=shift;
    my $num=0;
    while($num<$$r_threads){
        $semaphore->down();
        $num++;
    }
}
sub convert_tab{
    my($file,$outfile)=@_;
    open IN, $file or die "can not open file $file\n";
    open OUT, ">$outfile" or die "can not create file $outfile\n";
    print OUT join("\t","#miRNA", "target", "score", "mfe", "mfe_ratio", "start")."\n";
    while(<IN>){
        next if /^#/;
        chomp;
        if(/^miRNA\s+(\S+)/){
            print OUT $1."\t";
        }
        elsif(/^target\s+(\S+)/){
            print OUT $1."\t";
        }
        elsif(/^score\s+(\S+)/){
            print OUT $1."\t";
        }
        elsif(/^mfe\s+(\S+)/){
            print OUT $1."\t";
        }
        elsif(/^mfe_ratio\s+([\d\.]+)/){
            print OUT $1."\t";
        }
        elsif(/^start\s+(\d+)/){
            print OUT $1."\n";
        }
    }
    close IN;
    close OUT;
}