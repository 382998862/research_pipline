#!/usr/bin/perl

use warnings;
use strict;
use Getopt::Std;
use File::Copy;
use File::Path;
use threads;
use Thread::Semaphore;


my $usage=
"$0 file_command_line file_structure rounds_controls
-m  threads number
-o  output permuted file
-a   Output progress to screen
";

my $file_command_line=shift or die $usage;
my $file_structure=shift or die $usage;
my $rounds=shift or die $usage;


#options
my %options=();
getopts("m:o:a",\%options);
die $usage unless $options{'o'};

my $threads=1;
$threads=$options{'m'} if(exists $options{'m'});


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

if($threads > $cores){ print STDERR "More threads specified than cores on the system. Reducing the number of threads to $cores\n"; $threads=$cores;}

my $semaphore=new Thread::Semaphore($threads);

my %hash_desc;
my %hash_seq;
my %hash_struct;
my %hash_mfe;

my $ltime=time();
my $dir="dir_perform_controls$ltime";
mkdir $dir;

my $command_line=parse_file_command_line($file_command_line);

if($options{a}){print STDERR "total number of rounds controls=$rounds\n";}

parse_file_struct($file_structure);
perform_controls();
rmtree($dir);

sub perform_controls{
    my $round=1;
    my @all_tmp;
    while($round<=$rounds){
        if($options{a}){print STDERR "$round\r";}
        my $command_line_tmp=$command_line;
        permute_order("$dir/precursors_permuted.str$round");
        my $tmp_out="$dir/output_permuted.mrd$round";
        push @all_tmp,$tmp_out;
        $command_line_tmp=~s/precursors_permuted.str/precursors_permuted.str$round/;
        $command_line_tmp.=" >> ".$tmp_out;
        $command_line_tmp.=" 2> /dev/null";
        
        `echo 'permutation $round\n\n' >$tmp_out`;
        $semaphore->down();
        my $thr = threads->new(\&run_control, \$command_line_tmp);
        $thr->detach();
        
        $round++;
    }
    waitquit(\$threads);
    cat2(@all_tmp,$options{'o'});
    
    if($options{a}){print STDERR "controls performed\n\n";}
}
sub run_control{
    my($command)=@_;
    #open OUT,">$$out";
    #my $ret=`$$command`;
    #print OUT "permutation $$round\n\n";
    #print OUT "$ret";
    #close OUT;
    system($$command);
    $semaphore->up();
}
sub cat2{
    my $out=pop;
    open OUT,">$out";
    foreach(@_){
        open IN,$_;
        while(<IN>){
            print OUT $_;
        }
        close IN;
    }
    close OUT;
}
sub waitquit{
    my $r_threads=shift;
    my $num=0;
    while($num<$$r_threads){
        $semaphore->down();
        $num++;
    }
}

sub parse_file_command_line{

    my ($file) = @_;

    open (FILE, "<$file") or die "can not open $file\n";
    while (my $line=<FILE>){

	if($line=~/(\S+)/){
	    
	    chomp $line;

	    $line=~s/$file_structure/$dir\/precursors_permuted.str/;
	    
	    $line=~s/>.+//;

	    return $line;

	}
    }
    die "$file is empty\n";
}
sub permute_order{
    my $permute_file=shift;
    open OUT,">$permute_file";
    my @ids=sort keys %hash_seq;
    my $number_ids=scalar @ids;

    my @ids_cur=@ids;
    my $number_ids_cur=$number_ids;

    for(my $i=0; $i<$number_ids; $i++){

        my $id=$ids[$i];

        my $rand=int(rand($number_ids_cur));
        my $id_cur=$ids_cur[$rand];
        
        my $seq=$hash_seq{$id_cur};
        my $struct=$hash_struct{$id_cur};
        my $mfe=$hash_mfe{$id_cur};

        splice(@ids_cur,$rand,1);
    #	delete($hash_seq{$id_cur});

        $number_ids_cur--;

        print OUT ">$id\n$seq\n$struct ($mfe)\n";
    }
    close OUT;
}

sub parse_file_struct{
    #parses the output from RNAfoldand reads it into hashes

    my($file) = @_;
    my($id,$desc,$seq,$struct,$mfe) = ();

    open (FILE_STRUCT, "<$file") or die "can not open $file\n";
    while (<FILE_STRUCT>)
    {
        chomp;
        if (/^>(\S+)\s*(.*)/)
	{
	    $id          = $1;
	    $desc        = $2;
	    $seq         = "";
	    $struct      = "";
	    $mfe         = "";
	    while (<FILE_STRUCT>){
                chomp;
                if (/^>(\S+)\s*(.*)/){
		    $hash_desc{$id}   = $desc;
		    $hash_seq{$id}    = $seq;
		    $hash_struct{$id} = $struct;
		    $hash_mfe{$id}    = $mfe;

		    $id          = $1;
		    $desc        = $2;
		    $seq         = "";
		    $struct      = "";
		    $mfe         = "";

		    next;
                }
		if(/^\w/){
#		    tr/uU/tT/;
		    $seq .= $_;
		}if(/((\.|\(|\))+)/){
		    $struct .=$1;
		}
		if(/\((\s*-\d+\.\d+)\)/){
		    $mfe = $1;
		}
	    
	    }
        }
    }

    $hash_desc{$id}        = $desc;
    $hash_seq{$id}         = $seq;
    $hash_struct{$id}      = $struct;
    $hash_mfe{$id}         = $mfe;

    close FILE_STRUCT;
    return;
}

