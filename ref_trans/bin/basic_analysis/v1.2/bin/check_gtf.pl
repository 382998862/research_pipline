#!/usr/bin/perl -w
use strict;
use Data::Dumper;

if (@ARGV != 3) {
	print "\n\tFunction: check and correct gtf file\n\n";
	print "\tUsage: <gtf> <gff> <out_gtf>\n\n";
	exit;
}

my $gtf=$ARGV[0];
my $gff=$ARGV[1];
my $out=$ARGV[2];

open(IN,"$gff")||die $!;

my %hash;
while (<IN>) {
    if (/mRNA/) {
        my $mRNA;
        my $gene;
        #my @line=split/\t+/,$_;
        if (/ID=(.+?)\;/) {
            #print "mRNA:$1\n";
            $mRNA=$1;
        }
        if (/Parent=(.+)/i) {
            my $gene=$1;
		$gene=~s/\;$//;
#            $gene=$1;
            #print "gene:$1\n";
        }
        $hash{$mRNA}=$gene;
        
    }
    
}
close IN;

open(OU,">$out")||die $!;
open(IN,"$gtf")||die $!;
my $count=0;
while (<IN>) {
    chomp;
    if ($_=~/transcript_id/ and $_=~/gene_id/) {
        print OU "$_\n";
    }
    elsif($_=~/transcript_id \"(.+?)\"/ and $_!~/gene_id/){
        my $line=$_;
        my $id=$1;
        unless (defined $hash{$id}){
            print "warning:no gene id found for $id\n";
            next;
        }
        $line.=" gene_id \"".$hash{$id}."\";";
        #print OU "$line\n";
        $count++;
        next;
    }
    else{
        print "unrecogize line for $_\n";
    }
    
}
close OU;
close IN;
print "warnings:there are $count transcripts removed\n";



