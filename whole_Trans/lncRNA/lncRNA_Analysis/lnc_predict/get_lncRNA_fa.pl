my ($gtf,$genome,$output) =@ARGV;

###########################
if(scalar(@ARGV)!=3){
	print "perl $0\n1.GTF\n2.genome.fa\n3.lncRNA.fa\n";
	exit;
}

###########################
$/=">";
open(GENOME,$genome)||die $!;
my %fa=();
while(<GENOME>){
	chop;
	my ($id,$seq)=split(/\n/,$_,2);
	$id=(split(/\s+/,$id))[0];
	$seq=~s/\n//g;
	$fa{$id}=$seq;
}
close(GENOME);

###########################
$/="\n";
my %lncs=();
open(GTF,$gtf)||die $!;
while(<GTF>){
	chomp;
	my @tmp=split(/\t/,$_);
	next if($tmp[2] ne "exon");
	$tmp[8]=~/transcript_id \"(.*?)\"\;/;
	my $id=$1;
	$lncs{$id}{chr}=$tmp[0];
	$lncs{$id}{strand}=$tmp[6];
	push @{$id},[$tmp[3],$tmp[4]];
}
close(GTF);
open(OUT ,">$output")||die $!;
foreach my $lnc(keys %lncs){
	my @exons=sort{$a[0] <=> $b[0]} @{$lnc};
	my $seq;
	for(my $i=0;$i<@exons;$i++){
		$seq .=substr($fa{$lncs{$lnc}{chr}},$exons[$i][0]-1,$exons[$i][1]-$exons[$i][0]+1);
	}
	$seq=&trans($seq)	if($lncs{$lnc}{strand} eq "-");
	print OUT ">$lnc\n$seq\n";
}
close(OUT);


sub trans{
	my $seq=shift;
	$seq=reverse $seq;
	$seq=~tr/ATCGatcg/TAGCtagc/;
	return $seq;
}
