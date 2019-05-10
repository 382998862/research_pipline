#!/usr/bin/perl -w
use strict;
use warnings;
use Getopt::Long;
use Data::Dumper;
use FindBin qw($Bin $Script);
use File::Basename qw(basename dirname);
use Term::ANSIColor qw(:constants);
use newPerlBase;
use newPerlBase;
my %config=%{readconf("$Bin/../../../../config/db_file.cfg")}; 

my ($table, $component_file, $function_file, $process_file,$outfile);

GetOptions(
    "help|?"    => \&USAGE,
    "table:s"   => \$table,
    "o:s"=>\$outfile,
    ) or &USAGE;
&USAGE unless ($table);

$component_file=$config{component_file};
$function_file=$config{function_file};
$process_file=$config{process_file};

$component_file	||="/share/nas2/database/GO/releases201703/cellular_component.txt";
$function_file	||="/share/nas2/database/GO/releases201703/molecular_function.txt";
$process_file	||="/share/nas2/database/GO/releases201703/biological_process.txt";
$outfile||="gene_classify.txt";


open(IN,"$table")||die $!;
my $head=<IN>;
chomp($head);
my @head=split/\t/,$head;
my $gene;my $go;
for(my $i=0;$i<=$#head;$i++){
    if ($head[$i]=~/ID/) {
        $gene=$i;
    }
    if ($head[$i]=~/GO_annotation/) {
        $go=$i;
    }
}
my %hash_go;
my %go_bef;
while (<IN>) {
    chomp;
    my @line=split/\t/,$_;
    $go_bef{$line[$gene]}=$line[$go];
    next if $line[$go] eq '--';
    my (@go_id)=$line[$go]=~/(GO\:\d+)/g;
    foreach my $go(@go_id){
        push @{$hash_go{$line[$gene]}},$go;
    }
}
close IN;

my %three_go_tree_hash=();
my %go_second;
my %go_hash;
my @three_go_tree_file = ($component_file, $function_file, $process_file);
foreach my $file (@three_go_tree_file) {
	open IN,$file or die "cannot open file $file, $!\n";
	my @array = ();
	my $name;
	my $key;
        my $go_sec;
	while (my $line = <IN>) {
		chomp($line);
		next if ($line =~ /^$|^\#/);
		next if ($line !~ /(^\s+>)|(^\s+%)/);
		if ($line =~ /^\s>(.*)\s;/) {
			$name = $1;
			$name =~ s/[_\-]/ /;
			next;
		}
		if ($line=~/^\s\s%/) {
			my @key_in=split(/\s+;\s+/,$line);
			$key_in[0]=~s/^\s\s%//;
			$key = $key_in[0];
			push @array, $key;
			my ($go_id) = $key_in[1] =~ /(GO:\d+)/;
                	$go_sec=$go_id;
                	push @{$three_go_tree_hash{$name}{$go_id}},$key;
                	if (exists $go_second{$name}{$key} and $go_second{$name}{$key} ne $go_id) {
                		print "two go id for $name:$key\n";
                	}
                	$go_second{$name}{$key}=$go_id;
                	push @{$go_hash{$go_id}},$go_sec;
		}else{
			my ($go_id) = $line =~ /(GO:\d+)/;
			next if (!defined $go_id);
                	push @{$three_go_tree_hash{$name}{$go_id}},$key;
                	push @{$go_hash{$go_id}},$go_sec;
		}
	}
	close IN;
}

open(OU,">$outfile.temp")||die $!;
foreach my $gene(keys %hash_go){
    print OU "$gene\t";
    my @go=@{$hash_go{$gene}};
    my %result;
    my $mark=0;
    foreach my $go(@go){
        my $mark=0;
        foreach my $first(keys %three_go_tree_hash){
            if(exists $three_go_tree_hash{$first}{$go}) {
                $mark=1;
                my @sec=@{$three_go_tree_hash{$first}{$go}};
                foreach my $cat(@sec){
                    $result{$first}{$cat}=$go_second{$first}{$cat};
                }
            }
        }
        if ($mark==0) {
#            print "$go has no secondary categary\n";
        }
        
    }
    foreach my $first(keys %result){
        foreach my $second(keys %{$result{$first}}){
            print OU "$first: $second ($result{$first}{$second});; ";
        }
    }
    print OU "\n";
}
close OU;
######################add###################
open(IN,"$table")||die $!;
open(IN1,"$outfile.temp")||die $!;
open(OU,">$outfile")||die $!;
my %go_sec;
<IN1>;
while (<IN1>) {
    chomp;
    my ($gene,$go)=split/\t/,$_,2;
    $go_sec{$gene}=$go;
}
close IN1;
$head=<IN>;
chomp($head);
print OU "$head\tGO_second_level_annotation\n";
while (<IN>) {
    chomp;
    my $gene=(split/\t/,$_)[0];
    $go_sec{$gene}||="--";
    print OU "$_\t$go_sec{$gene}\n";
}
close IN;
close OU;








sub USAGE {
	my $usage = <<"USAGE";
 ProgramName: $Script
     Version: v1.0
     Contact: Liuxiaoshuang <liuxs\@biomarker.com.cn> 
Program Date: 2017-01-04
      Modify: 
 Description: This script is used to get secondary GO classification.
              
       Usage: 
        Options:
        --table     <FILE>  ref_trans_full_table.xls, required
        --c  <FILE>  cellular_component.txt

        --f     <FILE>   molecular_function.txt
        --p      <FILE>   biological_process.txt
        --o    <FILE>   outfile
        Examples:
            perl $Script --table Ref_test/ceshi25/Web_Report/ref_trans_full_table.xls
USAGE
	print $usage;
	exit;
}
