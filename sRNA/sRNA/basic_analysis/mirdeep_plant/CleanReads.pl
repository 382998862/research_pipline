#!/usr/bin/perl -w

use strict;
use Getopt::Long;
use File::Basename;

####################################### USAGE ####################################################
my $help =
"
Clean raw sRNA-seq data and output the uniq reads in fasta format

Usage: CleanReads.pl -q -f fasta -i reads1.fa=rd1 -i reads2.fa=rd2

Options:
-h/--help                   Print help message and quit
-v/--version                Print version and quit
-q/--quiet                  Quiet mode .. no log/progress information to STDERR
-i/--intput [input1.fa=rd1] input files, should be in format like reads1.fa=rd1, 'rd1' denotes 
                            the three-letter prefix of the input file, 
                            and this option can exist many times
-a/--adap3 [string] :       3' adapter sequence (defaul 'TCGTATGCCGTCTTCTGCTTG')
-g/--adap5 [string] :       5' adapter sequence (defaul '^GTTCAGAGTTCTACAGTCCGACGATC')
-e/--error [float >0..1] :  Maximum allowed error rate (no. of errors divided by
                            the length of the matching region) (default: 0.1)
-m/--min [integer] :        Discard trimmed reads that are shorter than m. (default: 18)
-x/--max [integer] :        Discard trimmed reads that are longer than x. (default: 30)
-f/--format [string] :      Input file format; can be either 'fastq', 'fasta'. (default: 'fastq')

Dependencies (must be in PATH):
  cutadapt 1.3

Documentation: perldoc CleanReads.pl
";

# if there are no arguments, return the help message and quit
unless($ARGV[0]) {
    die "No arguments specified!\n$help";
}

###################################### INPUT #######################################################
# declare global variables
my ($opt_h,$opt_v,$opt_q,$opt_a,$opt_e,$opt_g,$opt_m,$opt_x,$opt_f,%opt_i);

# Initialize opt_p and opt_a and opt_m defaults
$opt_a = "TCGTATGCCGTCTTCTGCTTG";                   #3' adaptor sequence
$opt_g = "^GTTCAGAGTTCTACAGTCCGACGATC";             #5' adaptor sequence
$opt_m = 18;                                        #minimal length of cleaned reads
$opt_x = 30;                                        #maxmal length of cleaned reads
$opt_e = 0.1;                                       #error rate
$opt_f = "fastq";                                   #format must be "fastq" "fq" or "fasta" or "fa"
# get options
#getopts('a:g:m:x:e:f:p:qhv');
GetOptions(
    'input|i=s'     =>  \%opt_i,
    'adap3|a=s'     =>  \$opt_a,
    'adap5|g=s'     =>  \$opt_g,
    'min|m=i'       =>  \$opt_m,
    'max|x=i'       =>  \$opt_x,
    'error|e=f'     =>  \$opt_e,
    'format|f=s'    =>  \$opt_f,
    'help|h!'       =>  \$opt_h,
    'version|v!'    =>  \$opt_v,
    'quiet|q!'      =>  \$opt_q,
);
####################################### MAIN ########################################################
#check parameters
check_param();
unless($opt_q){
    print STDERR "############Clean reads and collapse to unique tag############\n";
}
# Start speaking to user, unless quiet mode is on
unless($opt_q) {
    print STDERR "\nCleanReads.pl version $version\n";
    print STDERR "Checking Dependencies\n";
}
# Check dependencies
unless(check_cutadapt()) {
        die "FAIL\n\n$help";
}

# Running CleanReads.pl for each of the input files
#my @all_cleaned;
foreach(keys %opt_i){
    my $prefix=$opt_i{$_};
    #run cutadapt
    unless($opt_q) {
        print STDERR "Clean and collapse reads for file $_\n";
    }
    my $tmp_file=make_clean_file($_);
    #create the uniq tag file
    my $ret_file=make_uniq_file($tmp_file,$prefix);
    #push(@all_cleaned,$ret_file);
    system "rm -f $tmp_file";
}
# Merge all the collapse files into one
#system("cat @all_cleaned");
unless($opt_q){
        print STDERR "Down\n";
        print STDERR "#" x 60;
        print STDERR "\n";
}
exit;

######################################### SUBS #####################################################
# Here be sub-routines

sub check_param{
    
    # If user wants help or version, give it, and die
    if($opt_h) {
        die "$help";
    }
    if($opt_v) {
        die "CleanReads.pl version $version\n";
    }
    # check input files
    unless(%opt_i){
        die "FATAL: Input files should not be empty!\n";
    }
    foreach(keys %opt_i){
        unless(-f $_){
            die "FATAL: Can not find file $_\n";
        }
    }
    # prefix should be three letter, and cannot be the same
    my %pref;
    foreach(values %opt_i){
        unless(/^\w\w\w$/ and not(/_/)){
#            die "prefix $_ does not contain exactly three alphabet letters\n";
        }
        $pref{$_}++;
    }
    foreach(values %pref){
        die "FATAL: The prefix must be different\n" unless $_==1;
    }
    # Option e must be a number between 0 and 1
    unless($opt_e =~ /^[\.\d]+$/) {
        die "FATAL: Invalid entry for option e. Must be a number >0 and <= 1\n\n$help";
    }
    unless(($opt_e > 0) and ($opt_e <= 1)) {
        die "FATAL: Invalid entry for option e. Must be a number >0 and <= 1\n\n$help";
    }
    # Option m and x must be an integer
    unless(($opt_m =~ /^\d+$/) and ($opt_x =~ /^\d+$/)) {
        die "FATAL: Option m and x must be an integer\n\n$help";
    }
    if($opt_m > $opt_x){
        die "FATAL: Option x must be larger than option m\n\n$help";
    }
    unless(($opt_f eq "fastq") or ($opt_f eq "fasta")){
        die "FATAL: Option f must be 'fastq' or 'fasta' or 'fa'\n\n$help";
    }
}
sub check_cutadapt {
    unless($opt_q) {
        print STDERR "\tcutadapt: ";
    }
    (open(RD, "cutadapt --version |")) || return 0;
    my $v = <RD>;
    close RD;
    unless($opt_q) {
        print STDERR "PASS $v\n";
    }
    return 1;
}
sub make_clean_file {
    my $file=shift;
    my @prefix=qw(.fa .fasta .fastq .fq);
    my $cleaned=basename($file,@prefix).".tmp";
    
    my $command = "cutadapt";
    if($opt_a) {
        $command .= " -a '$opt_a'";
    }
    if($opt_g) {
        $command .= " -g '$opt_g'";
    }
    if($opt_m){
        $command .= " -m $opt_m";
    }
    if($opt_x){
        $command .= " -M $opt_x";
    }
    if($opt_e){
        $command .= " -e $opt_e";
    }
    if($opt_f){
        $command .= " -f $opt_f";
    }
    $command .= " -O 5 -n 2";
    $command .= " $file >$cleaned";
    $command .= " 2>/dev/null";
    
    system "$command";
    return $cleaned;
}

sub make_uniq_file{
    my ($file,$prefix)=@_;
    my $outfile="collapsed_".$file.".fa";
    $outfile=~s/.tmp//;
    my %hash;
    (open IN,$file) || die "Cannot open $file!";
    (open OUT,">$outfile") || die "Cannot open $file!";
    my $i=0;
    my($mod,$left);
    if($opt_f eq "fastq"){
        $mod=4;
        $left=2;
    }
    elsif($opt_f eq "fasta"){
        $mod=2;
        $left=0;
    }
    while(<IN>){
        chomp;
        $i++;
        next unless $i%$mod==$left;
        next if /N/;
        $hash{$_}++;
    }
    close IN;
    my ($uniq);
    foreach(sort {$hash{$b}<=>$hash{$a}} keys %hash){
        $uniq++;
        print OUT ">".$prefix."_".$uniq."_x".$hash{$_}."\n";
        print OUT $_."\n";
    }
    close OUT;
    return $outfile;
}
