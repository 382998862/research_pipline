#!/usr/bin/perl -w
use strict;
use warnings;
use Getopt::Long;
use Data::Dumper;
use FindBin qw($Bin $Script);
use File::Basename qw(basename dirname);

use lib ("$Bin/PerlLib");
use Gene_obj;
use Nuc_translator;
use Fasta_reader;
use Longest_orf;
use List::Util qw (min max);

my $BEGIN_TIME=time();
my $version="1.0.0";

# ------------------------------------------------------------------
# GetOptions
# ------------------------------------------------------------------
my ($transcripts_file, $prefix, $odir, $pfam_db, $genetic_code, $min_pep_len, $retain_long_ORFs, $top_ORFs_train, $cup_threads_num, $top_strand_only, $verbose);

GetOptions(
    "trans:s"=>\$transcripts_file,
    "pref:s" =>\$prefix,
    "odir:s" =>\$odir,
    "pfam:s" =>\$pfam_db,
    "code:s" =>\$genetic_code,
    "minL:i" =>\$min_pep_len,
    "saveL:i"=>\$retain_long_ORFs,
    "hmmn:i" =>\$top_ORFs_train,
    "cup:i"  =>\$cup_threads_num,
    "ss"     =>\$top_strand_only,
    "verb"   =>\$verbose,
    "help|h" =>\&HELP,
    ) or &USAGE;
&USAGE unless ($transcripts_file and $prefix);

chomp(my $wd=`pwd`);
$odir             ||= $wd;
$pfam_db          ||= "/share/nas19/yangxh/Database/pfam/27.0/Pfam-AB.hmm";
$genetic_code     ||= 'universal';
$min_pep_len      ||= 50;
$retain_long_ORFs ||= 900;
$top_ORFs_train   ||= 500;
$cup_threads_num  ||= 4;
$top_strand_only = (defined $top_strand_only)? 1 : 0;

system "mkdir -p $odir" unless (-d $odir);
$odir = &ABSOLUTE_DIR($odir);
die "ERROR: Don't understand options: @ARGV.\n" if (@ARGV);

my $util_dir = "$Bin/util";
$|++; #Output autoflush,Êä³ö»º³å¿ØÖÆ,Á¢¼´Ë¢ÐÂ»º³åÇø
our $SEE = $verbose;
my $NO_REPLACE = 0; #default replace

#my $hmmscan_prog = `which hmmscan`; # check for hmmscan utility
#$hmmscan_prog =~ s/\s//g;
#die "ERROR: Cannot locate pfam database at: $pfam_db.\n" unless (-s $pfam_db);
#die "ERROR: Cannot locate 'hmmscan' program in your PATH setting, so can't search pfam.\n" unless ($hmmscan_prog);
my $hmmscan_prog ="/share/nas2/genome/biosoft/hmmer/3.1b2/bin/hmmscan";

#######################################################################################
print "\n[".&GetTime($BEGIN_TIME)."] $Script start ...\n";
#######################################################################################
# ------------------------------------------------------------------
# get longest candidate CDS sequences.
# ------------------------------------------------------------------
&Nuc_translator::use_specified_genetic_code($genetic_code) unless ($genetic_code eq 'universal');

my $cds_file = "$odir/$prefix.longest_orfs.cds";
my $gff3_file = "$odir/$prefix.longest_orfs.gff3";
my $pep_file = "$odir/$prefix.longest_orfs.pep";
my %orf_lengths;

unless ($NO_REPLACE && -s $pep_file && -s $cds_file && -s $gff3_file) {
    open (PEP, ">$pep_file") or die $!;
    open (CDS, ">$cds_file") or die $!; 
    open (GFF, ">$gff3_file") or die $!;

    my $counter = 0;
    my $fasta_reader = new Fasta_reader($transcripts_file);

    while (my $seq_obj = $fasta_reader->next()) {
        my $acc = $seq_obj->get_accession();
        my $sequence = $seq_obj->get_sequence();
        my $orf_counter = 0;

        my $longest_orf_finder = new Longest_orf();
        $longest_orf_finder->allow_5prime_partials();
        $longest_orf_finder->allow_3prime_partials();
        $longest_orf_finder->forward_strand_only() if ($top_strand_only);

        my @orf_structs = $longest_orf_finder->capture_all_ORFs($sequence);
        @orf_structs = reverse sort {$a->{length}<=>$b->{length}} @orf_structs;

        while (@orf_structs) {
            my $orf = shift @orf_structs;
            my $start = $orf->{start};
            my $stop = $orf->{stop};

            $stop += 3 if ($stop <= 0); # edge issue

            my $nt_length = $orf->{length};
            my $aa_length = int($orf->{length}/3);
            my $orient = $orf->{orient};
            my $protein = $orf->{protein};

            next if ($aa_length < $min_pep_len);

            my $coords_href = { $start => $stop };
            my $gene_obj = new Gene_obj();

            $counter++;
            $orf_counter++;
            $gene_obj->populate_gene_object($coords_href, $coords_href);
            $gene_obj->{asmbl_id} = $acc;

            my $model_id = "$acc|orf$orf_counter|m.$counter";
            my $gene_id = "$acc|orf$orf_counter|g.$counter";

            $gene_obj->{TU_feat_name} = $gene_id;
            $gene_obj->{Model_feat_name} = $model_id;

            my $cds = $gene_obj->create_CDS_sequence(\$sequence);

            my $got_start = 0;
            my $got_stop = 0;

            $got_start = 1 if ($protein =~ /^M/);
            $got_stop = 1 if ($protein =~ /\*$/);

            my $prot_type = "";

            if ($got_start && $got_stop) {
                $prot_type = "complete";
            } elsif ($got_start) {
                $prot_type = "5prime_partial";
            } elsif ($got_stop) {
                $prot_type = "3prime_partial";
            } else {
                $prot_type = "internal";
            }

            ($aa_length, $nt_length) =($aa_length."aa", $nt_length."nt");
            ($start,$stop) = ($stop,$start) if ($orient =~/-/);

            $gene_obj->{com_name} = "type=$prot_type len=$nt_length";

            print PEP ">$model_id $gene_id type=$prot_type len=$aa_length loc=$acc:$start-$stop:$orient\n$protein\n";
            print CDS ">$model_id $gene_id type=$prot_type len=$nt_length loc=$acc:$start-$stop:$orient\n$cds\n";
            print GFF $gene_obj->to_GFF3_format() . "\n";

            $orf_lengths{$model_id} = length($cds);
        }
    }

    close PEP;
    close CDS;
    close GFF;
}

# ------------------------------------------------------------------
# Train a Markov model based on longest candidate CDS sequences, score all candidates.
# ------------------------------------------------------------------
# get longest entries
my $top_cds_file = "$odir/$prefix.longest_orfs.cds.top_${top_ORFs_train}_longest";
my $cmd = "$util_dir/get_top_longest_fasta_entries.pl $cds_file $top_ORFs_train > $top_cds_file";
&process_cmd($cmd) unless ($NO_REPLACE && -s $top_cds_file);

#get base frequencies
$cmd = "$util_dir/compute_base_probs.pl $transcripts_file $top_strand_only > $odir/$prefix.base_freqs.dat";
&process_cmd($cmd) unless ($NO_REPLACE && -s "$odir/$prefix.base_freqs.dat");

# get hexamer scores
$cmd = "$util_dir/seq_n_baseprobs_to_logliklihood_vals.pl $top_cds_file $odir/$prefix.base_freqs.dat > $odir/$prefix.hexamer.scores";
&process_cmd($cmd) unless ($NO_REPLACE && -s "$odir/$prefix.hexamer.scores");

# score all cds entries
$cmd = "$util_dir/score_CDS_liklihood_all_6_frames.pl $cds_file $odir/$prefix.hexamer.scores > $odir/$prefix.longest_orfs.cds.scores";
&process_cmd($cmd) unless ($NO_REPLACE && -s "$odir/$prefix.longest_orfs.cds.scores");

# ------------------------------------------------------------------
# run pfam search.
# ------------------------------------------------------------------
my %has_pfam_hit;
$cmd = "$hmmscan_prog --cpu $cup_threads_num --noali --cut_nc --acc --notextw --tblout $odir/$prefix.longest_orfs.pep.pfam.tab $pfam_db $odir/$prefix.longest_orfs.pep > $odir/$prefix.longest_orfs.pep.pfam.out";
&process_cmd($cmd) unless ($NO_REPLACE && -s "$odir/$prefix.longest_orfs.pep.pfam.tab");

open (PFAM, "$odir/$prefix.longest_orfs.pep.pfam.tab") or die $!;

while (<PFAM>) {
    chomp;
    my @col = split(/\s+/);
    my $orf_acc = $col[2];
    $has_pfam_hit{$orf_acc} = 1;
}
close PFAM;

# ------------------------------------------------------------------
# select the final set.
# ------------------------------------------------------------------
# get accs for best entries
open (OUT, ">$odir/$prefix.longest_orfs.cds.scores.selected") or die $!;
open (IN, "$cds_file.scores") or die $!;

while (<IN>) {
    chomp;
    my ($acc, @scores) = split /\t/;

    my $score_1 = shift @scores;
    my $max_score_other_frame = max(@scores);

    print OUT "$acc\n" if ($has_pfam_hit{$acc} || $orf_lengths{$acc} >= $retain_long_ORFs || ($score_1 > 0 && $score_1 > $max_score_other_frame));
}

close IN;
close OUT;

# index the current gff file:
$cmd = "$util_dir/index_gff3_files_by_isoform.pl $gff3_file";
&process_cmd($cmd);

# retrieve the best entries:
$cmd = "$util_dir/gene_list_to_gff.pl $odir/$prefix.longest_orfs.cds.scores.selected $gff3_file.inx > $odir/$prefix.best_candidates.gff3";
&process_cmd($cmd);

# exclude shadow orfs (smaller orfs in different reading frame that are eclipsed by longer orfs)
$cmd = "$util_dir/remove_eclipsed_ORFs.pl $odir/$prefix.best_candidates.gff3 > $odir/$prefix.best_candidates.final.gff3";
&process_cmd($cmd);

# ------------------------------------------------------------------
# output final set.
# ------------------------------------------------------------------
# make a BED file for viewing in IGV
$cmd = "$util_dir/gff3_file_to_bed.pl $odir/$prefix.best_candidates.final.gff3 > $odir/$prefix.best_candidates.final.bed";
&process_cmd($cmd);

# make a peptide file
$cmd = "$util_dir/gff3_file_to_proteins.pl $odir/$prefix.best_candidates.final.gff3 $transcripts_file > $odir/$prefix.best_candidates.final.pep.fa";
&process_cmd($cmd);

# make a CDS file
$cmd = "$util_dir/gff3_file_to_proteins.pl $odir/$prefix.best_candidates.final.gff3 $transcripts_file CDS > $odir/$prefix.best_candidates.final.cds.fa";
&process_cmd($cmd);

# ------------------------------------------------------------------
# abstract pfam annotation of final best candidates cds.
# ------------------------------------------------------------------
# abstract pfam annotation 
$cmd = "$util_dir/abstract_pfam_annotation.pl $odir/$prefix.best_candidates.final.gff3 $odir/$prefix.longest_orfs.pep.pfam.tab $odir/$prefix.best_candidates.final.pfam.tab $odir/$prefix.best_candidates.final.pfam.anno";
&process_cmd($cmd);

unless (defined $verbose) {
    system "rm $gff3_file $cds_file $pep_file $top_cds_file";
    system "rm $odir/$prefix.longest_orfs.pep.pfam.out $odir/$prefix.longest_orfs.pep.pfam.tab $odir/$prefix.longest_orfs.cds.scores $odir/$prefix.longest_orfs.cds.scores.selected";
    system "rm $odir/$prefix.base_freqs.dat $odir/$prefix.longest_orfs.gff3.inx $odir/$prefix.hexamer.scores $odir/$prefix.best_candidates.gff3";
}

#######################################################################################
print "\n[".&GetTime(time())."] $Script done. Total elapsed time: ".(time()-$BEGIN_TIME)."s\n";
#######################################################################################
# ------------------------------------------------------------------
# sub function
# ------------------------------------------------------------------
################################################################################################################
sub ABSOLUTE_DIR{ #$pavfile=&ABSOLUTE_DIR($pavfile);
    my $cur_dir=`pwd`; chomp($cur_dir);
    my ($in)=@_;
    my $return="";
    if(-f $in){
        my $dir=dirname($in);
        my $file=basename($in);
        chdir $dir; $dir=`pwd`; chomp $dir;
        $return="$dir/$file";
    }elsif(-d $in){
        chdir $in; $return=`pwd`; chomp $return;
    }else{
        warn "Warning just for file and dir\n";
        exit;
    }
    chdir $cur_dir;
    return $return;
}

################################################################################################################
sub process_cmd {
    my $cmd = shift;
    print "[".&GetTime(time())."] $cmd\n";
    my $ret = system $cmd;
    die "ERROR: CMD: $cmd died with ret $ret" if ($ret);
    return;
}

################################################################################################################
sub GetTime {
    my ($sec, $min, $hour, $day, $mon, $year, $wday, $yday, $isdst)=localtime(time());
    return sprintf("%4d-%02d-%02d %02d:%02d:%02d", $year+1900, $mon+1, $day, $hour, $min, $sec);
}

################################################################################################################
sub HELP {
    my $help=<<"_HELP_";
#-------------------------------------------------------------------------------------------------
    Program: $Script
    Version: $version
 Programmer: Brian Haas <bhaas\@broadinstitute.org> && Alexie Papanicolaou 
       Data: 2012-08-15
    Contact: Simon Young <simonyoung8824\@gmail.com>
       Data: 2014-09-27
    Fuction: the script is used to identify candidate coding regions within transcript sequences, such as those 
             generated by de novo RNA-Seq transcript assembly using Trinity, or constructed based on RNA-Seq 
             alignments to the genome using Tophat and Cufflinks.

Description: TransDecoder identifies likely coding sequences based on the following criteria:
             1) a minimum length open reading frame (ORF) is found in a transcript sequence.
             2) a log-likelihood score similar to what is computed by the GeneID software is > 0.
             3) the above coding score is greatest when the ORF is scored in the 1st reading frame as 
                compared to scores in the other 5 reading frames.
             4) if a candidate ORF is found fully encapsulated by the coordinates of another candidate ORF,
                the longer one is reported. However, a single transcript can report multiple ORFs (allowing
                for operons, chimeras, etc).
             5) the putative peptide has a match to a Pfam domain above the noise cutoff score.

      Usage: perl $Script --trans transcript.fa --pref prefix [options]
             --trans <STR>   input transcript file, FASTA format.
             --pref  <STR>   output files\' prefix.

    Options:
             --odir  <STR>   output directory.                                                   [./]
             --pfam  <STR>   path to Pfam database, *.hmm.                                       [/share/nas19/yangxh/Database/pfam/27.0/Pfam-AB.hmm]
                             using hmmscan to search Pfam domain as ORF retention criteria.
             --code  <STR>   genetic code, universal|Euplotes|Tetrahymena|Candida|Acetabularia.  ['universal']
             --minL  <INT>   minimun peptide length.                                             [50]
             --saveL <INT>   retain all ORFs found that are of minimum length in nucleotides.    [900]
             --hmmn  <INT>   top longest ORFs to train Markov Model.                             [500]
             --cup   <INT>   number of threads to use by hmmscan.                                [4]
             --ss            strand-specific, only analyzes top strand.
             --verb          verbose, save intermediate result.
             --help          show this help document.

    Example:
             perl $Script --trans Litchi.Transcripts.fa --pref Litchi --odir Orf_predict 
             perl $Script --trans Litchi.Transcripts.fa --pref Litchi --odir Orf_predict --minL 100 --cup 6 

     Result:
            Intermediate:
          x prefix.longest_orfs.gff3                : positions of all ORFs as found in the target transcripts.
          x prefix.longest_orfs.pep                 : all ORFs meeting the minimum length criteria, regardless of coding potential.
          x prefix.longest_orfs.pep.pfam.out        : Pfam domain search standard output of all detected ORFs.
          x prefix.longest_orfs.pep.pfam.tab        : Pfam domain search result table of all detected ORFs.
          x prefix.longest_orfs.cds                 : the nucleotide coding sequence for all detected ORFs.
          x prefix.longest_orfs.cds.top_500_longest : the top 500 longest ORFs, used for training a Markov model for CDS.
          x prefix.hexamer.scores                   : log likelihood score for each k-mer (coding/random).
          x prefix.longest_orfs.cds.scores          : the log likelihood sum scores for each ORF across each of the 6 reading frames.
          x prefix.longest_orfs.cds.scores.selected : the ids of the ORFs that were selected based on the scoring criteria.
          x prefix.best_candidates.gff3             : the positions of the selected ORFs in transcripts.

            Final:
            prefix.best_candidates.final.gff3       : positions within the target transcripts of the final selected ORFs.
            prefix.best_candidates.final.bed        : bed-formatted file describing ORF positions, best for viewing using IGV.
          * prefix.best_candidates.final.cds.fa     : nucleotide sequences for coding regions of the final candidate ORFs.
          * prefix.best_candidates.final.pep.fa     : peptide sequences for the final candidate ORFs.
            prefix.best_candidates.final.pfam.tab   : Pfam domain search result table of the final candidate ORFs.
          * prefix.best_candidates.final.pfam.anno  : Pfam annotation of the final candidate ORFs.

#-------------------------------------------------------------------------------------------------
_HELP_
    print $help;
    exit;
}

################################################################################################################
sub USAGE {
    my $usage=<<"USAGE";
#-------------------------------------------------------------------------------------------------
Version: $version
Fuction: the script is used to identify candidate coding regions within transcript sequences and search Pfam domain.

  Usage: perl $Script --trans transcript.fa --pref prefix [options]

Options:
        --odir  <STR>   output directory.                                                   [./]
        --pfam  <STR>   path to Pfam database, *.hmm.                                       [\$Bin/pfam/Pfam-AB.hmm]
                        using hmmscan to search Pfam domain as ORF retention criteria.
        --code  <STR>   genetic code, universal|Euplotes|Tetrahymena|Candida|Acetabularia.  ['universal']
        --minL  <INT>   minimun peptide length.                                             [50]
        --saveL <INT>   retain all ORFs found that are of minimum length in nucleotides.    [900]
        --hmmn  <INT>   number of longest ORFs to train Markov Model.                       [500]
        --cup   <INT>   number of threads to use to searche Pfam domain.                    [4]
        --ss            strand-specific, only analyzes top strand.
        --help          help document.

#-------------------------------------------------------------------------------------------------
USAGE
    print $usage;
    exit;
}
