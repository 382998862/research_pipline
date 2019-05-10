#!/usr/bin/perl -w
use strict;
use Getopt::Long;
use Cwd 'abs_path';
use FindBin qw($Bin $Script);
use File::Basename qw(basename dirname);
use Data::Dumper;
use Algorithm::Combinatorics qw(combinations permutations);

#===============================================================
#-------------------------------------------------
# get options
#-------------------------------------------------
my ($snp_file, $refSeq_file, $odir);

GetOptions(
        "help|?" =>\&HELP,
        "snp:s"  =>\$snp_file,
        "ref:s"  =>\$refSeq_file,
        "od:s"   =>\$odir,
        ) or &USAGE;
&USAGE unless ($snp_file and $refSeq_file);

#-------------------------------------------------
# preparation
#-------------------------------------------------
$odir ||= "./";
system "mkdir -p $odir/Pairwised_SNP/" unless (-d "$odir/Pairwised_SNP/");
$odir = abs_path($odir);

my %gene_len;
my @sample;
my %snp;
my %density_stat; #	$density_stat{sample_id}{density}=gene_num;
my %pairwise_snp_stat;
#-------------------------------------------------
# get ref seq len
#-------------------------------------------------
open (REF,$refSeq_file) or die $! ;
$/=">";<REF>;

while (<REF>) {
    my ($id,$seq) = (split /\n/,$_,2);
    $seq =~s/\s+//g;
    $gene_len{$id} = length $seq;
}

$/="\n";
close REF;

#-------------------------------------------------
# load merged snp file
#-------------------------------------------------
my @field;
my %snp_per_gene;

open (SNP, $snp_file) or die $!;

while (<SNP>) {
    chomp;

    if ($.==1) {
        @field = (split /\s+/);
        for (my $i=4; $i<@field-3; $i=$i+3) {
            push @sample, $field[$i];
        }
    } else {
        die "ERROR: the input SNP file may be illegal, please check.\n" if (@field==0);
        my @col = (split /\s+/);

        $snp{$col[0]}{$col[1]}->{'refAlle'} = $col[2];
        $snp{$col[0]}{$col[1]}->{'altAlle'} = $col[3];
        $snp_per_gene{'All'}{$col[0]}++;

        for (my $i=4; $i<@field-3; $i=$i+3) {
            $snp{$col[0]}{$col[1]}->{$field[$i]}->{'genotype'} = $col[$i];
            $snp{$col[0]}{$col[1]}->{$field[$i]}->{'totalDep'} = $col[$i+1];
            $snp{$col[0]}{$col[1]}->{$field[$i]}->{'alleDeps'} = $col[$i+2];

            $snp_per_gene{$field[$i]}{$col[0]}++ unless ($col[$i+1] eq '.');
        }
    }
}

close SNP;

if (@sample<2) {
    print "ERROR: less than 2 samples.\n";
    exit;
}

#-------------------------------------------------
# output pair-wised sample snp
#-------------------------------------------------
$,="\t";$\="\n";
my $iter = combinations(\@sample,2);

while (my $c = $iter->next) {
    my ($boy, $girl) = @$c;
    my $header = join "\t",('#GeneID', 'Postion', 'RefAlle', 'AltAlles', "$boy.genotype", "$boy.totalDep", "$boy.alleDeps", "$girl.genotype", "$girl.totalDep", "$girl.alleDeps");

    open (OUT, ">$odir/Pairwised_SNP/$boy.$girl.parwised_snp.list") or die $! ;
    print OUT $header;
    open (OUT1, ">$odir/Pairwised_SNP/$boy.$girl.parwised_snp.hete_hete.list") or die $! ;
    open (OUT2, ">$odir/Pairwised_SNP/$boy.$girl.parwised_snp.hete_homo.list") or die $! ;
    open (OUT3, ">$odir/Pairwised_SNP/$boy.$girl.parwised_snp.homo_hete.list") or die $! ;
    open (OUT4, ">$odir/Pairwised_SNP/$boy.$girl.parwised_snp.homo_homo.list") or die $! ;
    print OUT1 $header;
    print OUT2 $header;
    print OUT3 $header;
    print OUT4 $header;

    for my $g (sort {$a cmp $b} keys %snp) {
        for my $p (sort {$a <=> $b} keys %{$snp{$g}}) {
            if ($snp{$g}{$p}->{$boy}->{'totalDep'} ne '.' and
                $snp{$g}{$p}->{$girl}->{'totalDep'} ne '.' and
                $snp{$g}{$p}->{$boy}->{'genotype'} ne $snp{$g}{$p}->{$girl}->{'genotype'}) {
                my $boy_inf = join "\t",($snp{$g}{$p}->{$boy}->{'genotype'}, $snp{$g}{$p}->{$boy}->{'totalDep'}, $snp{$g}{$p}->{$boy}->{'alleDeps'});
                my $girl_inf = join "\t",($snp{$g}{$p}->{$girl}->{'genotype'}, $snp{$g}{$p}->{$girl}->{'totalDep'}, $snp{$g}{$p}->{$girl}->{'alleDeps'});
                print OUT $g, $p, $snp{$g}{$p}->{'refAlle'}, $snp{$g}{$p}->{'altAlle'}, $boy_inf, $girl_inf;
                $pairwise_snp_stat{"${boy}_vs_$girl"}{'Total'}++;

                if ($snp{$g}{$p}->{$boy}->{'genotype'} =~/[^ATCG]/i) {  #hete boy
                    if ($snp{$g}{$p}->{$girl}->{'genotype'} =~/[^ATCG]/i) { #hete girl
                        print OUT1 $g, $p, $snp{$g}{$p}->{'refAlle'}, $snp{$g}{$p}->{'altAlle'}, $boy_inf, $girl_inf;
                        $pairwise_snp_stat{"${boy}_vs_$girl"}{'S1.hete.S2.hete'}++;
                    } else {    #homo boy
                        print OUT2 $g, $p, $snp{$g}{$p}->{'refAlle'}, $snp{$g}{$p}->{'altAlle'}, $boy_inf, $girl_inf;
                        $pairwise_snp_stat{"${boy}_vs_$girl"}{'S1.hete.S2.homo'}++;
                    }
                } else {    # homo boy
                    if ($snp{$g}{$p}->{$girl}->{'genotype'} =~/[^ATCG]/i) { #hete girl
                        print OUT3 $g, $p, $snp{$g}{$p}->{'refAlle'}, $snp{$g}{$p}->{'altAlle'}, $boy_inf, $girl_inf;
                        $pairwise_snp_stat{"${boy}_vs_$girl"}{'S1.homo.S2.hete'}++;
                    } else {    #homo boy
                        print OUT4 $g, $p, $snp{$g}{$p}->{'refAlle'}, $snp{$g}{$p}->{'altAlle'}, $boy_inf, $girl_inf;
                        $pairwise_snp_stat{"${boy}_vs_$girl"}{'S1.homo.S2.homo'}++;
                    }
                }

            }
        }
    }

    close OUT;
    close OUT1;
    close OUT2;
    close OUT3;
    close OUT4;
}

open (STAT, ">$odir/Pairwised_SNP/parwised_snp.stat") or die $!;
print STAT '#Type', 'S1.homo.S2.homo', 'S1.hete.S2.homo', 'S1.homo.S2.hete', 'S1.hete.S2.hete', 'Total';

for my $couple (sort keys %pairwise_snp_stat) {
    for my $type ('S1.homo.S2.homo', 'S1.hete.S2.homo', 'S1.homo.S2.hete', 'S1.hete.S2.hete', 'Total') {
        $pairwise_snp_stat{$couple}{$type} = 0 unless (exists $pairwise_snp_stat{$couple}{$type} );
    }
    print STAT $couple, 
        $pairwise_snp_stat{$couple}{'S1.homo.S2.homo'}, 
        $pairwise_snp_stat{$couple}{'S1.hete.S2.homo'}, 
        $pairwise_snp_stat{$couple}{'S1.homo.S2.hete'}, 
        $pairwise_snp_stat{$couple}{'S1.hete.S2.hete'}, 
        $pairwise_snp_stat{$couple}{'Total'};
}

close STAT;

#-------------------------------------------------
# snp density stat & plot 
#-------------------------------------------------
# snp density stat
for my $s (@sample,'All') {
    for my $g (keys %gene_len) {
#        next unless (exists $snp_per_gene{$s}{$g});    #take gene with snp into account only
        $snp_per_gene{$s}{$g} = 0 unless (exists $snp_per_gene{$s}{$g});    #take gene without snp into account also
        my $snp_per_gene_per_kb = $snp_per_gene{$s}{$g}*1000/$gene_len{$g};

        if ($snp_per_gene_per_kb < 1) {
            $density_stat{$s}{'0-1'}++;
        } elsif ($snp_per_gene_per_kb < 2) {
            $density_stat{$s}{'1-2'}++;
        } elsif ($snp_per_gene_per_kb < 3) {
            $density_stat{$s}{'2-3'}++;
        } elsif ($snp_per_gene_per_kb < 4) {
            $density_stat{$s}{'3-4'}++;
        } elsif ($snp_per_gene_per_kb < 5) {
            $density_stat{$s}{'4-5'}++;
        } elsif ($snp_per_gene_per_kb < 6) {
            $density_stat{$s}{'5-6'}++;
        } elsif ($snp_per_gene_per_kb < 7) {
            $density_stat{$s}{'6-7'}++;
        } elsif ($snp_per_gene_per_kb < 8) {
            $density_stat{$s}{'7-8'}++;
        } else {
            $density_stat{$s}{'8~'}++;
        }
    }
}
=c
# output snp density
open (STAT,">$odir/AllSample.SNP_density.stat") or die $!;
print STAT '#Sample', 'Interval', 'GeneNum';

for my $s (sort keys %density_stat) {
    for my $c ('0-1','1-2','2-3','3-4','4-5','5-6','6-7','7-8','8~') {
        if (exists $density_stat{$s}{$c}) {
            print STAT $s, $c, $density_stat{$s}{$c};
        } else {
            print STAT $s, $c, '0';
        }
    }
}

close STAT;

# density plot
my $Rscript = "/share/nas2/genome/biosoft/R/3.1.1/bin/Rscript";
my $cmd = "$Rscript $Bin/dodgedBar.r --infile $odir/AllSample.SNP_density.stat --outfile $odir/AllSample.SNP_density.png --x.col 2 --group.col 1 --y.col 3 --x.lab \"SNP Number per Kb\" --group.lab \"Sample\" --y.lab \"Number of Gene\" --title.lab \"SNP Density\" --legend.col 1 >/dev/null 2>&1 ";
#print "$cmd\n";
system $cmd;
=cut

# ==============================================================
# sub function
# ==============================================================
sub USAGE {#
	my $usage=<<"_USAGE_";
Options:
    --ref   <file>  reference sequence file, FASTA format, forced
    --snp   <file>  marged SNP table files, tab-delimited , forced
    --od    <path>  output directory, optional, default CWD
    --h             help documents
Examples:
    perl $Script --ref Maize.Unigene.fa --snp SNP_Analysis/SNP/final.snp.list --od SNP_Analysis/

_USAGE_
	print $usage;
	exit;
}

sub HELP {#
	my $usage=<<"_HELP_";
 Program: $Script
 Version: v1.0
 Contact: Simon Young <yangxh\@biomarker.com.cn> 
Function: this script is used to abstract pair-wised SNP from marged SNP table,
          and survey and plot the density of SNP. suit for no_ref trans pipline.
 Options:
    --ref   <file>  reference sequence file, FASTA format, forced
    --snp   <file>  marged SNP table file, tab-delimited , forced
    --od    <path>  output directory, optional, default CWD
    --h             help documents
Examples:
    perl $Script --ref Maize.Unigene.fa --snp SNP_Analysis/SNP/final.snp.list --od SNP_Analysis/

 Details:
1. example of the input marged SNP table file:
--------------------------------------------------------------------------------
#Chr            Pos     Ref     Alt     T1      Depth   AlleDp  T2      Depth   AlleDp  T3      Depth   AlleDp  T4      Depth   AlleDp  Effect  Codon_change    Gene_id
comp1001_c1     665     T       C       N       .       .       N       .       .       Y       3       1,2     T       2       2,0     UNKNOWN
comp1004_c0     257     G       C       S       195     148,47  S       153     106,47  S       52      36,16   S       36      24,12   UNKNOWN
--------------------------------------------------------------------------------
This file must be a tab-delimited table and contain a header line start with "#" and some SNP lines of at lest 2 samples' info.
The SNP info lines must contain the following fields in turn:
gene id, snp potions, refernece allele, alter allele, and at lest 2 sample' genotype, total depth, and depth of every allele.

2. illustration of the output:
--------------------------------------------------------------------------------
YOUR_OUTPUT_DIR/
    |-- AllSample.SNP_density.stat      # snp density statistic result
    |-- AllSample.SNP_density.png       # snp density plot
    `-- Pairwised_SNP                   # pair-wised SNP table directory
        |-- parwised_snp.stat                   # pair-wised SNP statistic table
        |-- S1.S2.parwised_snp.hete_hete.list   # SNP table, both sample is heterozygous
        |-- S1.S2.parwised_snp.hete_homo.list   # SNP table, smaple S1 is heterozygous, while S2 is homozygous
        |-- S1.S2.parwised_snp.homo_hete.list   # SNP table, smaple S1 is homozygous, while S2 is heterozygous
        |-- S1.S2.parwised_snp.homo_homo.list   # SNP table, both sample is homozygous
        |-- S1.S2.parwised_snp.list             # SNP table, genotype is different between sample S1 and S2.
        `-- ... ...                             # etc.
--------------------------------------------------------------------------------
_HELP_
	print $usage;
	exit;
}