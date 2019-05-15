my($i,$g2go)=@ARGV;
if(@ARGV!=2){
	print "perl $0 \n1.input Unigene_Annotation/Result/Known.longest_transcript.fa.GO.anno.txt \n2.outfile1 contained gene and GO id info\n";
	exit;
}
open(IN,$i)||die $!;
open(OUT,">$g2go")||die $!;

while(<IN>){
	chomp;
	next if($_=~/^#/);
	my @tmp=split(/\t/,$_); 
	my $gene=shift @tmp;
	my $num=shift @tmp;
	foreach my $t(@tmp){
		$t=~/(.*?): (.*?) \((GO:.*?)\)/;
		print OUT "$3\t$2\t$gene\t$1\n";
	}
	print OUT "";

}
close(OUT);
close(IN);

