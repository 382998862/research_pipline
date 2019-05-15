#!/usr/bin/perl
use strict;
use warnings;


open( IN, "$ARGV[0]" )     || die "open file failed\n";
open( FA, "$ARGV[1] ")	|| die "open file failed\n";
open( FILE, "$ARGV[2]" )   || die "open file failed\n";
open( OUT,  ">$ARGV[3]" ) || die "open or create file failed\n";
my %hash;
while (<IN>) {
	chomp;
	my @line = split/\t/,$_;
	$hash{$line[0]}=$line[1];
}
close IN;

my %hc;
$/=">";
<FA>;
while(<FA>){
	chomp;
	my ($name,$seq)=(split /\n/,$_,2);
	my $id = (split /\s+/,$name)[0];
	s/\n//;s/\r//;
	my $len = length $seq;
	$len = $len-1;
	$hc{$id}=$len;
}
$/="\n";
my %ex;
print OUT "#circRNA_id\tcircBase_id\tcircRNA_length\tcircBase_length\tidentity\tblast_length\tblast_mismatch\tblast_gap\tblast_start1\tblast_end1\tblast_start2\tblast_end2\tblast_Evalue\tblast_bitscore\tratio1\tratio2\n";
while(<FILE>)
{
	chomp;
	next if(/#/);
	my @line = split/\t/,$_;
#	if(!exists $ex{$line[0]})
#	{
		if(exists $hash{$line[1]})
		{
			if(exists $hc{$line[0]}){
				my $length = $hash{$line[1]};
				my $length1 = $hc{$line[0]};
				my $ratio = (($line[3]-$line[4]-$line[5])/$length)*100;
				my $ratio1 = (($line[3]-$line[4]-$line[5])/$length1)*100;
				next if($ratio<90 or $line[2]<90 or $ratio1<90);
				$ex{$line[0]}=1;
				printf OUT ("%s\t%s\t%s\t%s\t%.2f\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%.2f\t%.2f\n", $line[0],$line[1],$length1,$length,$line[2],$line[3],$line[4],$line[5],$line[6],$line[7],$line[8],$line[9],$line[10],$line[11],$ratio1,$ratio);
			}
		}
#	}
	
}
