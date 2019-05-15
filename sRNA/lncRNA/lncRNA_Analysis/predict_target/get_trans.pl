my ($input,$output,$cor,$pvalue)=@ARGV;
if(@ARGV!=4){
	print "perl $0\n1.input file(Contained 4 column)#RNA1    RNA2    coefficient     pvalue\n2.output file \n3.coefficient \n4.pvalue\n";
	exit;
}
open(IN,$input)||die $!;
my %lncs=();
while(<IN>){
	chomp;
	next if($_=~/^#|RNA1/);
	my ($lnc,$gene,$c,$p)=split(/\t/,$_);
	next if($p>$pvalue);
	next if(abs($c)<$cor);
	$lncs{$lnc}++;
	push @{$lnc},$gene;
	
}
close(IN);
open(OUT,">$output")||die $!;
print OUT "#ID\tTrans_target_gene\n";
foreach my $lnc(keys %lncs){
	print OUT "$lnc\t",join("\;",@{$lnc}),"\n";
}
close(OUT);

