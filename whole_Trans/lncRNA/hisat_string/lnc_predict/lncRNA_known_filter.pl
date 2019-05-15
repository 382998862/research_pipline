use FindBin qw($Bin $Script);
use File::Basename qw(basename dirname);

my ($known,$gtf,$od)=@ARGV;
if(@ARGV!=3){
	print "perl $0\n1. Known LncRNA gtf\n2.Candicate lncRNA gtf.\n3. output dir\n";
	exit;
}
chdir $od;
`$Bin/../bin/gffcompare -r $known $gtf`;
my $base=basename $gtf;
open(REF,"$od/gffcmp.$base.refmap")||die $!;
my %filter=();
while(<REF>){
	chomp;
	next if($_=~/^ref_gene_id/);
	my @tmp=split(/\t/,$_,4);
	my @a = split /,/,$tmp[3];
	foreach my $k(@a){
		my $id=(split(/\|/,$k))[1];
		print "$id\n";
		$filter{$id}++;
	}
}
close(REF);

open(GTF,$gtf)||die $!;
open(OUT,">$od/unKnown_lncRNA.gtf")||die $!;
while(my $line=<GTF>){
	chomp($line);
	my @tmp=split(/\t/,$line);	
	$tmp[8]=~/transcript_id (.*?)\;/;
	my $id=$1;
	$id=~s/\"//g;
	if(!exists $filter{$id}){
		print OUT "$line\n";
	}
}
close(OUT);
close(GTF);

