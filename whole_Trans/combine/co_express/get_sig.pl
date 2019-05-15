my($i,$o,$cor,$p)=@ARGV;
if(@ARGV!=4){
	print "perl $0\n1. Total coexpression file\n2. Significant coexpression file\n3.Correlation coefficient \n4.pvalue\n";
	exit;
}
open(IN,$i)||die $!;
open(OUT,">$o")||die $!;
my $head=<IN>;chomp($head);
print OUT "$head\n";
while(<IN>){
	chomp;
	my ($r1,$r2,$co,$pvalue)=split(/\t/,$_);
	next if( $pvalue > $p);
	next if( abs($co) < $cor);
	print OUT "$_\n";
}
close(OUT);
close(IN);

