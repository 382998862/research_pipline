#!/usr/bin/env perl

use strict;
use warnings;
use FindBin;
use lib ("$FindBin::Bin/../PerlLib");
use Gene_obj;
use GFF3_utils;
use Carp;
use Nuc_translator;

my $usage = "\nusage: $0 gff3_file pfam_tab_in pfam_tab_out pfam_anno_out\n\n";

my $gff3_file     = $ARGV[0] or die $usage;
my $pfam_tab_in   = $ARGV[1] or die $usage;
my $pfam_tab_out  = $ARGV[2] or die $usage;
my $pfam_anno_out = $ARGV[3] or die $usage;

## get id list
my %id_list;

my $gene_obj_indexer_href = {};
## associate gene identifiers with contig id's.
my $contig_to_gene_list_href = &GFF3_utils::index_GFF3_gene_objs($gff3_file, $gene_obj_indexer_href);

foreach my $asmbl_id (sort keys %$contig_to_gene_list_href) {
    my @gene_ids = @{$contig_to_gene_list_href->{$asmbl_id}};
    foreach my $gene_id (@gene_ids) {
        my $gene_obj_ref = $gene_obj_indexer_href->{$gene_id};
        foreach my $isoform ($gene_obj_ref, $gene_obj_ref->get_additional_isoforms()) {
            my $isoform_id = $isoform->{Model_feat_name};
            $id_list{$isoform_id} = 1;
        }
    }
}


# output selected lines.
my %raw_id_list;

open (IN, "$pfam_tab_in") or die $!;
open (OUT, ">$pfam_tab_out") or die $!;

while (<IN>) {
    if (/^#/) {
        print OUT $_;
    } else {
        chomp;
        my @col = (split /\s+/,$_,19);
        my ($id, $pfam_id, $pfam_descript) = (split /\s+/,$_,19)[2,1,18];
        my $raw_id =$id;
        $raw_id =~ s/\|.+$//;
        $raw_id_list{$raw_id}{'pfam_ids'} .= "$pfam_id;";
        $raw_id_list{$raw_id}{'pfam_description'} .= "$pfam_descript;; ";

        die "ERROR: Illegal pfam_tab_in.\n" unless (defined $col[18]);
        print OUT "$_\n" if (exists $id_list{$id});
    }
}

close OUT;
close IN;

# output selected annotation.

open (OUT, ">$pfam_anno_out") or die $!;
print OUT "#Gene_ID\tPfam_IDs\tPfam_Description\n";

for my $id (sort keys %raw_id_list) {
    my $pfam_ids = $raw_id_list{$id}{'pfam_ids'};
    my $pfam_descript = $raw_id_list{$id}{'pfam_description'};

    $pfam_ids =~ s/;$//;
    $pfam_descript =~ s/;; $//;
    print OUT "$id\t$pfam_ids\t$pfam_descript\n";
}

close OUT;

exit(0);

