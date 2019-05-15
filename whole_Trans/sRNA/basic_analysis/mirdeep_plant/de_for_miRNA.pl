#!/usr/bin/perl
use warnings;
use strict;
use Getopt::Long;
use Cwd 'getcwd';

####################################### USAGE ####################################################
my $usage =
"$0 -s result11.csv=sample1 -s result12.csv=sample1 -s result21.csv=sample2 -s result22.csv=sample2

Options:
-s/--sample     input files, should be in format like <file_path>=<sample_name>, and it can exist
                many times. The first two columns of the miRNA expression file must be miRNA name and its count number.
";

unless($ARGV[0]) {
    die "No arguments specified!\n";
}

my %opt_s;
my %all_sample;
my @samples;
my $factor;
my @names;
my %hash_de;

GetOptions('sample|s=s' => \%opt_s);
my $cur_path=getcwd(__FILE__);
####################################### MAIN ########################################################

foreach(keys %opt_s){
    push(@{$all_sample{$opt_s{$_}}},$_);
}
@samples=keys %all_sample;

die "There must be two sample categories!\n" unless @samples == 2;

foreach my $sample(@samples){
    for(my $i=0;$i<@{$all_sample{$sample}};$i++){
        $factor.="'".$sample."',";
        store_info($all_sample{$sample}->[$i],$sample.($i+1));
        push(@names,$sample.($i+1));
    }
}
chop($factor);

open OUT,">for_de_file.tmp" or die "can not create file\n";
print OUT join("\t","miRNA",@names)."\n";
foreach my $mirna(sort keys %hash_de){
    print OUT $mirna;
    foreach(@names){
        print OUT "\t".$hash_de{$mirna}{$_};
    }
    print OUT "\n";
}
close OUT;

create_R();
system("R --slave --vanilla <for_de.R 1>rlog.txt 2>&1");
clean();

sub store_info{
	my($file,$name)=@_;
	open IN,$file or die "can not open file $file";
	while(<IN>){
		next if /^#/;
		my($mirna,$read)=(split /\s+/)[0,1];
		$hash_de{$mirna}{$name}=$read;
	}
}

sub create_R{
    open OUTR,">for_de.R" or die "can not create file\n";
    my $out_file=<<TMP;
library(DESeq);
setwd("$cur_path");
data<-read.table("for_de_file.tmp",header=T,row.names=1);
condition<-factor(c($factor));
cds<-newCountDataSet(data,condition);
cds=estimateSizeFactors(cds);
if(length(condition)==2){
    cds=estimateDispersions(cds,method="blind", sharingMode="fit-only");
}else{
    cds=estimateDispersions(cds);
}
nbinomTest(cds,levels(condition))->result;
write.table(result,"dif_exp_miRNA.result",quote=F,sep="\\t",row.names=F);
TMP
    print OUTR $out_file;
    close OUTR
}
sub clean{
    unlink "for_de_file.tmp","for_de.R";
}