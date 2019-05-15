my($gff,$gtf,$id_name)=@ARGV;
if(@ARGV<2||@ARGV>3){
	print "perl $0\n1.input gff file\n2.output gtf file\n3.ouput id_symbol relation(not must)\n";
	exit;
}
my $exon="CDS";
my $flag=`awk  -F \$\"\t\" '{print \$3}' $gff|grep exon|wc -l`;
chomp($flag);
$flag=(split(/\s+/,$flag))[0];
$exon="exon"	if($flag>0);

open(GFF,$gff)||die $!;
my %genes=();
my %info=();
my %sym=();
while(<GFF>){
	chomp;
	next if($_=~/^#/);
	my ($chr,$source,$type,$start,$end,$f1,$strand,$f2,$info)=split(/\t/,$_);
	my ($id,$parent);
	if($type =~/gene|processed_transcript|lncRNA_gene/){
		$id=(split(/\;/,$info))[0];
		$id=~s/ID=|gene://g;
		$genes{$id}++;

		$info=~/Name=(.*?)\;/;
		$sym{$id}=$1;
	}elsif($type =~/transcript|lincRNA|lncRNA|RNA/){
		$info=~/ID=(.*?)\;Parent=(.*?)$/;
		($id,$parent)=($1,$2);
		$parent=(split(/\;/,$parent))[0];
		$id=~s/transcript://g;
		$parent=~s/gene://g;
		push @{$parent},$id;

		$info=~/Name=(.*?)\;/;
		$sym{$id}=$1;
	}elsif($type =~/$exon/){
		$info=~/ID=(.*?)\;Parent=(.*?)$/;
		$id=$1;
		my $parent=(split(/\;/,$2))[0];
		$parent=~s/transcript://g;
		push @{$parent},[$start,$end];	
	}	
	$info{$id}=join("\t",($chr,$source,$start,$end,$f1,$strand,$f2));
}
close(GFF);

open(GTF,">$gtf")||die $!;
foreach my $g(keys %genes){
	my ($chr,$source,$start,$end,$f1,$strand,$f2)=split(/\t/,$info{$g});
	foreach my $t(@{$g}){
		my @exons=sort{$a[0] <=> $b[0]} @{$t};
		for(my $i=0;$i<@exons;$i++){
			print GTF "$chr\t$source\texon\t$exons[$i][0]\t$exons[$i][1]\t.\t$strand\t.\ttranscript_id \"$t\"\; gene_id \"$g\"\;\n";
		}

	}
}
close(GTF);
if(defined $id_name){
	open(ID,">$id_name")||die $!;
	print ID"#ID\tsymbol\tSource_Gene_ID\tSource_Gene_symbol\n";
	foreach my $g(keys %genes){
		foreach my $t(@{$g}){
			print ID "$t\t$sym{$t}\t$g\t$sym{$g}\n";
		}
	}
	close(ID);
}
