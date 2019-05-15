#!/usr/bin/perl

=head1 Name

  cog_parser.pl  --  extract the cog number.

=head1 Description

  This program is designed for extracting the cog number after blast.

=head1 Version

  Author: sunjuan, sunjuan@genomics.org.cn
  Version: 1.3,  Date: 2008-1-31

=head1 Usage
	
  perl cog_parser.pl <blast_tab>
  --nohead      do not show the first instruction line.
  --verbose     output verbose information to screen  
  --help        output help information to screen  

=head1 Exmple

  perl cog_parser.pl PLASMID.fasta.ori.glimmer3.pep.blast.tab

=cut

use strict;
use Getopt::Long;
use FindBin qw($Bin $Script);
use File::Basename qw(basename dirname); 
use Data::Dumper;
use newPerlBase;
my ($Nohead,$Verbose,$Help);

GetOptions(
	"nohead"=>\$Nohead,
	"verbose"=>\$Verbose,
	"help"=>\$Help
);

die `pod2text $0` if (@ARGV == 0 || $Help);

my $blast_tab = shift;

my $seq_name = $1 if ($blast_tab=~/([^\/]+)\.Cog\.blast/);
my $outdir = dirname($blast_tab);
#print "$seq_name\n";

my %file_db=%{readconf("$Bin/../../db_file.cfg")}; ;

my $whog = $file_db{cog_whog_file};
my $fun = $file_db{cog_fun_file};

my (%COG,%fun);
my $category;

##read fun.txt
open FUN,$fun || die "fail $fun";

while (<FUN>) {
	if (/(\[\w\])\s(.+?)\n$/) {
		my ($class,$function) = ($1,$2);
		$fun{$class}{defination} = $function;
		$fun{$class}{category} = $category;
	}else {
		chomp $_;
		$category = $_;
	}
}
close FUN;

##read whog
open WHOG,$whog || die "fail $whog";
while (<WHOG>) {
	chomp;
	my @line = split /,/,$_;
	my ($protein_id,$class,$cog_num,$cog_anno) = ($line[0],$line[8],$line[6],$line[9]);
	$cog_anno =~ s/\s+/ /g;
	$fun{$class}{category} = "--" unless (exists $fun{$class}{category});
	$fun{$class}{defination} = "--" unless (exists $fun{$class}{defination});
	$COG{$protein_id} = "$cog_num\t$cog_anno\t$class\t$fun{$class}{category}\t$fun{$class}{defination}";
}
close WHOG;

##read the tab file and create a file including cog info
open IN,$blast_tab || die "fail $blast_tab";
open OUT,">$outdir/$seq_name.Cog.class" || die "fail $seq_name.Cog.class";
#print OUT "#Gene name\tPortein ID in COG\tE_value\tIdentity\tScore\tCOG id\tCOG class defination\tFunction code\tFunctional categories\tFunction class defination\n" unless (defined $Nohead);
while (<IN>) {
	next if /^\#/;
	my @t = split /\t/;
        my @ids=split/\|/,$t[4];
	my $protein_ID=$ids[1];
	print OUT "$t[0]\t$t[4]\t$t[13]\t$t[8]\t$t[12]\t$COG{$protein_ID}\n" if (exists $COG{$protein_ID});
}
close OUT;
close IN;

sub file_cfg_read {
    my ($cfg_file, $file_db) = @_;

    open (CFG,$cfg_file ) or die "$!: $cfg_file\n";
    while (<CFG>) {
        chomp;
        next if (/^\s*$/ or /^#/);
        my ($key, $value) = (split /\s+/)[0,1];
            $file_db->{$key} = $value;
    }
    close CFG;
}
