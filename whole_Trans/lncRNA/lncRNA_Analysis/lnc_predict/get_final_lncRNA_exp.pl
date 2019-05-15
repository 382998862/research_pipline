use Cwd qw(abs_path getcwd);

my ($prep,$lncRNA_id)=@ARGV;
if(@ARGV!=3 && @ARGV!=2){
	print "perl $0\n1. quantify prep dir\n2. Filter Final lncRNA id\n3. output dir(default prep dir). Two files(All_lncRNA_count.list, All_lncRNA_fpkm.list would be produced!)\n";
	exit;
}
$prep=abs_path($prep);
$lncRNA_id=abs_path($lncRNA_id);

my %ids=();
open(ID,$lncRNA_id)||die $!;
while(<ID>){
	chomp;
	next if($_=~/^#/);
	my @tmp=split(/\s+/,$_);
	$ids{$tmp[0]}++;
}
close(ID);


&writeOut("$prep/Candicate_lncRNA_counts.list","$prep/Candicate_lncRNA_counts.list.tmp");
&writeOut("$prep/Candicate_lncRNA_fpkm.list","$prep/Candicate_lncRNA_fpkm.list.tmp");


if(-e "$prep/Known_lncRNA_fpkm.list"){
	`grep -v '^#' $prep/Known_lncRNA_counts.list >$prep/Known_lncRNA_counts.list.tmp`;
	`grep -v '^#' $prep/Known_lncRNA_fpkm.list >$prep/Known_lncRNA_fpkm.list.tmp`;

	`cat $prep/Candicate_lncRNA_counts.list.tmp $prep/Known_lncRNA_counts.list.tmp|sed 's/ /\\t/g' > $prep/All_lncRNA_counts.list `;
	`cat $prep/Candicate_lncRNA_fpkm.list.tmp $prep/Known_lncRNA_fpkm.list.tmp|awk 'NF-=1'|sed 's/ /\\t/g' > $prep/All_lncRNA_fpkm.list `;
}else{
	`awk 'NF-=1' $prep/Candicate_lncRNA_fpkm.list.tmp  |sed 's/ /\\t/g'> $prep/All_lncRNA_fpkm.list`;
	`awk 'NF-=1' $prep/Candicate_lncRNA_counts.list.tmp|sed 's/ /\\t/g' > $prep/All_lncRNA_counts.list`;
}
`rm $prep/*tmp`;


sub writeOut{
	my ($i,$o)=@_;
	open(IN,$i)||die $!;
	open(OUT,">$o")||die $!;
	while(<IN>){
		chomp;
		print OUT "$_\n"	if($_=~/^#/);
		my @tmp=split(/\t/,$_);
		if(exists $ids{$tmp[0]}){
			print OUT "$_\n";
		}
	}
	close(IN);
	close(OUT);

}
