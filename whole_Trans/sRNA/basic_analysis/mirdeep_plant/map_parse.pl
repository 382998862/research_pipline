#!/usr/bin/perl -w
use strict;
use Getopt::Std;
use File::Basename;

####################################### USAGE ####################################################

my $usage =
"
The script processes the reads and/or maps them to the reference genome, as designated by the options given.
The mapped file is then converted to an 'arf' file, which used in miRDeep2. 
Options:

-g genome         The genome file where input reads file will be mapped to
-l [int]          The seed length [default: 18]
-n [int]          Mismatches allowed in the seed [default: 2]
-m [int]          Suppress all alignments for a particular read or pair if more than <int> 
                  reportable alignments exist for it. [default: 15]
-q                Quiet mode .. no log/progress information to STDERR
-p [int]          Number of threads to use for bowtie 
-h                Print help message and quit
Example of use:

$0 collapse_reads.fa -g genome.fa -l 18 -n 2 -q -p 8 > reads_vs_genome.arf
";

###################################### INPUT #######################################################
unless($ARGV[0]) {
    die "$usage";
}
my $file_reads=shift or die $usage;
my %options=();
getopts("g:l:n:m:p:qh",\%options);

#################################### GLOBAL VARIABLES ################################################

my $mismatches_seed=2;
my $seed_len=18;
my $max_map=15;
my $threads=1;
my $genome;
$mismatches_seed=$options{'n'} if exists $options{'n'};
$seed_len=$options{'l'} if exists $options{'l'};
$max_map=$options{'m'} if exists $options{'m'};
$threads=$options{'p'} if exists $options{'p'};
$genome=$options{'g'} if exists $options{'g'};

# Check parameters
check_param();

####################################### MAIN ########################################################

unless($options{q}){
    print STDERR "############Map reads and parse to arf format############\n";
}
# mkdir bowtie_index
my $index="bowtie_index";
if(not -d $index){
        mkdir($index);
}
# Check dependencies
unless($options{q}){
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
# Main part
my @suffix=qw(.fa .fasta);
my $base_genome=basename($genome,@suffix);
my $base_file=basename($file_reads,@suffix);

my $ebwt_check = check_ebwt("bowtie_index/$base_genome");
# If absent, call bowtie-build
unless($ebwt_check) {
    bowtie_build($genome,"bowtie_index/$base_genome");
}
map_reads($file_reads,"$base_file.bwt");
convert_bowtie("$base_file.bwt",$base_file."_vs_genome.arf");

unless($options{q}){
        print STDERR "Down\n";
        print STDERR "#" x 60;
        print STDERR "\n";
}
exit;
######################################### SUBS #####################################################
sub check_param{
    
    # If user wants help or version, give it, and die
    if($options{'h'}) {
        die "$usage";
    }
    # check read file
    unless(-f "$file_reads"){
        die "No reads file $file_reads could be found\n$usage";
    }
     unless(-f $genome){
        die "No genome file $genome could be found\n$usage";
    }
    # check integer
    unless(($mismatches_seed =~ /^\d+$/) and ($seed_len =~ /^\d+$/) 
        and ($max_map=~/^\d+$/) and ($threads=~/^\d+$/)) {
        die "FATAL: Option m, l, n, p must be an integer\n\n$usage";
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
        print STDERR "More threads specified than cores on the system. Reducing the number of threads to $cores\n"; 
        $threads=$cores;
    }
}

sub check_bowtie_version {
    unless($options{'q'}) {
        print STDERR "\tbowtie: ";
    }
    (open(BTV, "bowtie --version |")) || return 0;
    my $vline = <BTV>;
    close BTV;
    my $version;
    if($vline =~ /^bowtie version (\S+)/) {
        $version = $1;
        if(($version =~ /^0\.12/) or ($version =~ /^1\./)) {
            unless($options{'q'}) {
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
    unless($options{'q'}) {
        print STDERR "\tbowtie-build: ";
    }
    (open(BBV, "bowtie-build --version |")) || return 0;
    my $vline = <BBV>;
    close BBV;
    if($vline =~ /^bowtie-build version (\S+)/) {
        unless($options{'q'}) {
            print STDERR "PASS version $1\n";
        }
        return 1;
    } 
    else { 
        return 0;
    }
}
sub check_ebwt {
    my($base) = @_;
    unless($options{'q'}) {
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
        unless($options{'q'}) {
            print STDERR "PRESENT\n";
        }
        return 1;
    } 
    elsif ($ebwt_count > 0) {
        unless($options{'q'}) {
            print STDERR "INCOMPLETE only $ebwt_count of the 6 expeted index files were found.\n";
        }
        die "\tABORTING .. please clean up old incompleted ebwt files for transcript $base and try again\n";
    } 
    else {
        unless($options{'q'}) {
            print STDERR "ABSENT\n";
        }
        return 0;
    }
}

sub bowtie_build {
    my($genome,$index) = @_;
    unless($options{'q'}) {
        print STDERR "Building bowtie index for file $genome ...\n";
    }
    system("bowtie-build $genome $index > /dev/null");
    return 1;
}

sub map_reads{
    #map reads to genome
    my($file_reads,$out_mapped)=@_;
    unless($options{q}){
        print STDERR "mapping reads to genome index...\n";
    }
    my $ret_mapping=`bowtie -p $threads -f -n $mismatches_seed -e 80 -l $seed_len -a -m $max_map --best --strata $index/$base_genome $file_reads $out_mapped`;
}

sub convert_bowtie{
    my ($mapped,$converted)=@_;
    my @line;
    my @gseq;
    my @changes;
    my $pos;
    my $ont;
    my $mm=0;
    my @edit="";
    my $reverse = 0;
    unless($options{q}){
        print STDERR "Converting mapped bwt file to arf format...\n";
    }
    open IN,"<$mapped" or die "Cannot open file $mapped\n";
    open OUT,">$converted" or die "Cannot create file $converted\n";
    while(<IN>){
        @line = split(/\t/);
        if($line[1] eq "-" ){
            $line[4] = reverse($line[4]);
            $line[4] =~ tr/ACGTN/TGCAN/;
        }
        @gseq = split(//,lc $line[4]);
        @edit= split(//,("m" x length($line[4])));
        $mm=0;
        
        if($line[7]){
            
            @changes = split(/,/,$line[7]);
            #$mm=scalar @changes;
            foreach(@changes){
                if(/(\d+):(\w+)\>\w+/){
                    $mm++;
                    $gseq[$1] = lc $2;
                    ## this was not outcommented before
                    #$gseq[$1] =~ tr/acgt/tgca/;
                    $edit[$1] = "M";   
                }
            }
        }
        
        my @id =split(/\s/,$line[0]);
        my @db = split(/\s/,$line[2]);
        my @removed=remove_trailing_nts(length($line[4]),length($line[4]),lc $line[4],length($line[4]),$line[3]+length($line[4]),join("",@gseq),$mm,join("",@edit));
        print OUT join("\t",$id[0],$removed[0],1,$removed[1],$removed[2],$db[0],$removed[3],($line[3]+1),$removed[4],$removed[5],$line[1],$removed[6],$removed[7])."\n";
    }
    close IN;
    close OUT;
}

sub remove_trailing_nts{

    my($query_map_lng,$query_end,$query_seq,$db_map_lng,$db_end,$db_seq,$edits,$edit_string)=@_;

    while($edit_string=~/M$/){

	$query_map_lng--;
	$query_end--;
	chop $query_seq;
	$db_map_lng--;
	$db_end--;
	chop $db_seq;
	$edits--;
	chop $edit_string;
    }
    return ($query_map_lng,$query_end,$query_seq,$db_map_lng,$db_end,$db_seq,$edits,$edit_string);
}

sub cat_to{
    
    my($file_1,$file_2)=@_;

    open OUT, ">>$file_2" or die "cannot print to $file_2\n";

    open IN, "<$file_1" or die "cannot read from $file_1\n";
    
    while(my $line = <IN>){

	print OUT "$line";
    }

    close IN;

    close OUT;

    return;
}

