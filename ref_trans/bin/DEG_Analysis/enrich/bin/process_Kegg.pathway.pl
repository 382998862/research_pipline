my($i,$o)=@ARGV;
if(@ARGV!=2){
	print "perl $0\n1.input Mouse_Unigene.Kegg.pathway \n2.output KEGG.file\n";
	exit;
}

open(IN,$i)||die $!;
open(OUT,">$o")||die $!;
while(<IN>){
	chomp;
	next if($_=~/^#/);
	my($path,$ko,$num,$gene,$K)=split(/\t/,$_); 
	my @genes=split(/\;/,$gene);
	foreach my $g(@genes){
		next if($g=~/^$/);
		print OUT "$ko\t$path\t$g\n";
	}
}
close(OUT);
close(IN);
