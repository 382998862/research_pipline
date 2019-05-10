my ($in,$gff,$od)=@ARGV;
if(@ARGV!=3){
	print "perl $0 1.in  2.gff  3.od\n";
	print "in: input Cuffnorm dir\n";
	print "od: ouput path will produce All_gene_counts.xls and All_gene_fpkm.xls\n";
	die;
}
`mkdir -p $od`	if(!-d $od);
my %gene=();
open(TRACK,"$in/genes.read_group_tracking")||die $!;
while(<TRACK>){
	chomp;
	next if($_=~/^tracking_id/);
	my ($g,$sample,$rep,$raw,$internal,$external,$FPKM,$len,$sta)=split(/\t/,$_);
	$gene{$sample}{$g}{raw}=$raw;
	$gene{$sample}{$g}{external}=$external;
	$gene{$sample}{$g}{FPKM}=$FPKM;
}

close(TRACK);
open(GFF,$gff)||die $!;
while(<GFF>){
	chomp;
	my @tmp=split(/\t/,$_);
	next if($tmp[2] ne "gene");
	$tmp[8]=~/ID=(.*?)$/;
	my $g=(split(/\;/,$1))[0];
	${$g}{locs}="$tmp[0]:$tmp[3]-$tmp[4]";
	${$g}{length}=$tmp[4]-$tmp[3]+1;
	${$g}{strand}=$tmp[6];
}
close(GFF);

my @samples=keys %gene;
foreach my $s(@samples){
	open(FPKM,">$od/$s.geneExpression.xls")||die $!;
	print FPKM "#GeneID\tLength\tFPKM\tLocus\tStrand\tCount\tNormalized\n";
	foreach my $g(sort{$a cmp $b} keys %{$gene{$s}}){
		print FPKM "$g\t${$g}{length}\t$gene{$s}{$g}{FPKM}\t${$g}{locs}\t${$g}{strand}\t$gene{$s}{$g}{raw}\t$gene{$s}{$g}{external}\n";
	}
	close(FPKM);
}

