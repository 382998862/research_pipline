use File::Basename qw(basename dirname);
use FindBin qw($Bin $Script);
use Cwd qw(abs_path getcwd);

my($lnc_fa,$od,$db)=@ARGV;
if(@ARGV!=3){
	print "perl $0\n1. lncRNA.fa \n2.output dir\n3. db (GRCh38/GRCm38)\n";
	die;
}
`mkdir $od`	unless(-d $od);
$od=abs_path($od);
$lnc_fa=abs_path($lnc_fa);

my $spe;
if($db eq "GRCh38" ||$db eq "GRCh37" ||$db eq "hg19"||$db eq "hg38"){
	$spe="Human";
}elsif($db eq "GRCm38" || $db eq "mm9" || $db eq "mm10"){
	$spe="Mouse";
}else{
	print "The given $db is not correct for this analysis!\n";
	exit;
}
my $blastall="/share/nas2/genome/bin/blastall";
my $cmd ="$blastall -p blastn -d $Bin/lncRNA_Disease.fa -i $lnc_fa -S 2  -F F -m 8 -a 4 -e 0.00001 -o $od/blast.out\n";
print $cmd,"\n";
 `$cmd`;
my %length=();
$/=">";
open(FA,"$Bin/lncRNA_Disease.fa")||die $!;
while(<FA>){
        chop;
        next if($_=~/^$/);
        my ($id,$seq)=split(/\n/,$_,2);
        $seq=~s/\n//g;
        $length{$id}=length($seq);
}
close(FA);
$/="\n";
#Query id;Subject id,%identity,alignment length,mismatches,gap
my %info=();
open(BLAST,"$od/blast.out")||die $!;
while(<BLAST>){
	chomp;
	my @tmp=split(/\t/,$_);	
	my $ratio=$tmp[3]/$length{$tmp[1]};
	next if($ratio < 0.8);
	my ($transcript,$gene)=split(/\./,$tmp[1],2);
	$info{$transcript}=$gene;
	push @{$transcript},$tmp[0];
}
close(BLAST);

open(BACK,"$Bin/data_v2017.txt")||die $!;
open(OUT,">$od/LncRNA_disease.txt")||die $!;
print OUT "#LncRNA_id\tDB_lncRNA\tSource\tdisease\tregulate\tdescription\tchr\tstart\tend\tstrand\tpmid\n";
while(<BACK>){
	chomp;
	my @tmp=split(/\t/,$_);
	next if($spe ne $tmp[-4]);
	my @flag1=split(/\;\s/,$tmp[-3]);
	push @flag1,$tmp[-2];
	foreach my $f(@flag1){
		if(exists $info{$f}){
			my @lncs=@{$f};
			foreach my $l(@lncs){
				print OUT "$l\t$f\t$tmp[0]\t$tmp[1]\t$tmp[2]\t$tmp[3]\t$tmp[4]\t$tmp[5]\t$tmp[6]\t$tmp[7]\t$tmp[-1]\n";

			}
		}
	}
}
close(BACK);
close(OUT);
