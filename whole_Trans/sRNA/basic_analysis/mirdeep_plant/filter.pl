#!/usr/bin/perl

use warnings;
use strict;
use Getopt::Long;
use File::Path;
use File::Basename;

####################################### USAGE ####################################################
my $usage =
"
This script takes as input a file with deep sequencing reads (these can be in
different formats, see the options below). The script then processes the reads
and/or maps them to the reference genome, as designated by the options given.
Options:

Preprocessing/mapping:
-i <string>     collapsed reads file in fasta format
-l <string>     library for reads filtering, users can appoint one or more libraries,
                and the library must be in fasta format
-v <int>        Report alignments with at most <int> mismatches. [default: 2]        
-k              keep the matched alignments
-m <int>        number of threads to lanuch. [Default: 1]

Other:
-q              Quiet mode .. no log/progress information to STDERR
-h              Print help message and quit 

Example of use:

$0 -l Rfam.fasta -l Repbase.fasta -i collapsed_reads.fasta -v 2 -k -m 4
";

# if there are no arguments, return the help message and quit
unless($ARGV[0]) {
    die "No arguments specified!\n$usage";
}
###################################### INPUT #######################################################
my @opt_libs=();
my ($opt_i,$opt_h,$opt_q,$opt_v,$opt_k,$opt_m);

GetOptions(
    'input|i=s'         => \$opt_i,
    'help|h!'           => \$opt_h,
    'quiet|q!'          => \$opt_q,
    'mismatches|v=i'    => \$opt_v,
    'threads|m=i'       => \$opt_m,
    'keep|k!'           => \$opt_k,
    'lib|l=s'           => \@opt_libs,
);

#################################### GLOBAL VARIABLES ################################################
my $threads=1;
$threads=$opt_m if $opt_m;

my $mismatches=2;
$mismatches=$opt_v if $opt_v;

my @prefix=qw(.fa .fasta);
my $base_file=basename($opt_i,@prefix);

####################################### MAIN ########################################################
unless($opt_q){
    print STDERR "############Filtering reads############\n";
}
# mkdir bowtie_index
my $index="bowtie_index";
if(not -d $index){
        mkdir($index);
}
# make tmp directory for filtering
my $dir="filter_reads_tmp";
mkdir $dir if not -d $dir;

# Check options
check_options();
# Check dependencies
unless($opt_q){
    print STDERR "Checking dependencies\n";
}
my $bowtie_check = check_bowtie_version();
unless($bowtie_check) {
    die "FAIL\n\n$usage";
}
my $bowtie_build_check = check_bowtie_build_version();
unless($bowtie_build_check) {
    die "FAIL\n\n$usage";
}



# transform RNA to DNA
my @dna_libs=rna2dna(@opt_libs);

my $latest_read=$opt_i;
my $latest_unmapped=$dir."/".$base_file."_filter";
my $latest_mapped=$base_file;
foreach(@dna_libs){
    # build bowtie index
    my $lib_index=basename($_);
    my $ebwt_check = check_ebwt("$index/$lib_index");
    # If absent, call bowtie-build
    unless($ebwt_check) {
        bowtie_build($_, "$index/$lib_index");
    }

    # call bowtie
    # Keep it quiet if needed
    unless($opt_q) {
        print STDERR "Running bowtie for mapping to $lib_index...\n"; 
    }
    $latest_unmapped.="_".$lib_index;
    $latest_mapped.="_vs_".$lib_index;
    system "bowtie -f -v $mismatches --best -a -p $threads --norc --un $latest_unmapped --al $latest_mapped $index/$lib_index $latest_read >/dev/null";
    $latest_read=$latest_unmapped;
}
system("mv $latest_unmapped filtered_${base_file}.fa");
rmtree($dir);

unless($opt_q){
        print STDERR "Down\n";
        print STDERR "#" x 60;
        print STDERR "\n";
}
exit;
######################################### SUBS #####################################################
sub check_options{
    if($opt_h) {
        die "$usage";
    }
    foreach(@opt_libs,$opt_i){
        die "Can not find file $_\n$usage" unless "-f $_";
    }
    unless($mismatches=~/^\d+$/){
        die "FATAL: allowed mismatches must be an integer\n\n$usage";
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
sub rna2dna{
    my @files = @_;
    my @ret;
    foreach(@files){
        my $dna_out=$dir."/".basename($_);
        push(@ret,$dna_out);
        open DNA,">$dna_out" or die "can not create $dna_out\n";
        open (FASTA, "<$_") or die "can not open $_\n";
        while (<FASTA>){
            if (/^>/){
                print DNA $_;
            }
            else{
                s/U/T/ig;
                print DNA uc($_);
            }
        }
        close DNA;
        close FASTA;
    }
    
    return @ret;
}
sub check_bowtie_version {
    unless($opt_q) {
        print STDERR "\tbowtie: ";
    }
    (open(BTV, "bowtie --version |")) || return 0;
    my $vline = <BTV>;
    close BTV;
    my $version;
    if($vline =~ /^bowtie version (\S+)/) {
        $version = $1;
        if(($version =~ /^0\.12/) or ($version =~ /^1\./)) {
            unless($opt_q) {
            print STDERR "PASS version $version\n";
            }
            return 1;
        } 
        else {
            return 0;
        }
    } 
    else {
        return 0;
    }
    # should never be here
    return 0;
}

sub check_bowtie_build_version {
    unless($opt_q) {
        print STDERR "\tbowtie-build: ";
    }
    (open(BBV, "bowtie-build --version |")) || return 0;
    my $vline = <BBV>;
    close BBV;
    if($vline =~ /^bowtie-build version (\S+)/) {
        unless($opt_q) {
            print STDERR "PASS version $1\n";
        }
        return 1;
    } 
    else { 
        return 0;
    }
}

sub bowtie_build {
    my($genome,$index) = @_;
    unless($opt_q) {
        print STDERR "Building bowtie index for file $index ...\n";
    }
    system("bowtie-build $genome $index > /dev/null");
    return 1;
}
sub merge_libs{
    my @files=@_;
    my $ret=$dir."/filter_library.fasta";
    system("cat @files > $ret");
    return $ret;
}
sub check_ebwt {
    my($base) = @_;
    unless($opt_q) {
        print STDERR "\tBowtie index files for $base ";
    }
    my @ebwt_suffixes = qw(.1.ebwt .2.ebwt .3.ebwt .4.ebwt .rev.1.ebwt .rev.2.ebwt);
    my $ebwt_count = 0;
    my $ebwt_file;
    foreach my $ebwt_suffix (@ebwt_suffixes) {
        $ebwt_file = "$base" . "$ebwt_suffix";
        if(-r $ebwt_file) {
            ++$ebwt_count;
        }
    }
    if($ebwt_count == 6) {
        unless($opt_q) {
            print STDERR "PRESENT\n";
        }
        return 1;
    } 
    elsif ($ebwt_count > 0) {
        unless($opt_q) {
            print STDERR "INCOMPLETE only $ebwt_count of the 6 expeted index files were found.\n";
        }
        die "\tABORTING .. please clean up old incompleted ebwt files for transcript $base and try again\n";
    } 
    else {
        unless($opt_q) {
            print STDERR "ABSENT\n";
        }
        return 0;
    }
}