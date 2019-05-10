my($i,$backgroud,$o,$index)=@ARGV;
if(@ARGV<3){
	print "perl $0 inputfile backgroupfile outputfile  index(not must,defalut 1)\n";
	exit;
}
$index||=1;
open(BACK,$backgroud)||die $!;
my %hash=();
my $head=<BACK>;chomp($head);my @header=split(/\t/,$head);shift @header;
my $num=scalar(@header);
my $temp="--";
for(my $j=1;$j<$num;$j++){
	$temp.="\t--";
}
while(<BACK>){
	chomp;
	my @tmp=split(/\t/,$_);
	my $id=shift @tmp;
	$hash{$id}=join("\t",@tmp);
}
close(BACK);

open(IN,$i)||die $!;
open(OUT,">$o")||die $!;
while(<IN>){
	chomp;next if($_=~/^$/);my @tmp=split(/\t/,$_);
	if($_=~/^#/){
		my $head="$_\t".join("\t",@header);
		print OUT "$head\n";
	}else{
		if(exists $hash{$tmp[$index-1]}){
			print OUT "$_\t$hash{$tmp[$index-1]}\n";
		}else{
			print OUT "$_\t$temp\n";
		}
	}
}
close(IN);
close(OUT);
