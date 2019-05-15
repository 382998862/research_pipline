#!/usr/bin/perl -w
my ($lnc_fa,$gff,$o)=@ARGV;
if(@ARGV!=3){
	print "perl $0\n1.Candicate lncRNA fa \n2.Final lncRNA gff \n3.Final lncRNA fa(output)\n";
	exit;
}
my %rnas=();
open(GFF,$gff)||die $!;
while(<GFF>){
	chomp;
	my @tmp=split(/\t/,$_);
	if($tmp[2] eq "mRNA")	{
		$tmp[8]=~/ID=(.*?)\;/;
		my $id=$1;
		$id=(split(/\;/,$id))[0];
		$rnas{$id}="$tmp[0] $tmp[3]-$tmp[4]";
	}
}
close(GFF);
$/=">";
open(IN,$lnc_fa)||die $!;
open(OUT,">$o")||die $!;
while(<IN>){
	chop;
	next if($_=~/^$/);
	my($id,$seq)=split(/\n/,$_,2);
	$seq=~s/\n//g;
	$id=(split(/\s+/,$id))[0];
	if(exists $rnas{$id}){
		print OUT ">$id $rnas{$id}\n$seq\n";
	}
}
close(OUT);
close(IN);

