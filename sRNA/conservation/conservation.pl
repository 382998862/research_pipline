use Getopt::Long;
use Data::Dumper;
use Config::General;
use FindBin qw($Bin $Script);
use File::Basename qw(basename dirname);
use Cwd qw(abs_path getcwd);
my $extract="/share/nas2/genome/biosoft/kentUtils-master/bin/bigWigAverageOverBed";
#######################################################################################
# ------------------------------------------------------------------
# GetOptions, Check format, Output info.
# ------------------------------------------------------------------
my ($bw, $gff,$type,$od,$gtf,$bed,$db);
GetOptions(
				"help|?"	=>\&USAGE,
				"bw:s"		=>\$bw,
				"db:s"		=>\$db,
				"bed:s"		=>\$bed,
				"gtf:s"		=>\$gtf,
				"gff:s"		=>\$gff,
				"od:s"		=>\$od,
				"key:s"		=>\$key,
				) or &USAGE;
&USAGE unless ($od) ;

$key||="Conservation";

mkdir $od unless(-d $od);
$od=abs_path($od);
if(!defined $bed && !defined $gff && !defined $gtf){
	&USAGE;
}
if(!defined $bw){
	if(defined $db){
		my %config=&readConfig("$Bin/conservation.cfg");
		$bw=$config{$db};
	}else{
		&USAGE;
	}
}



if(defined $gtf){
	$gtf=abs_path($gtf);
	open(GTF,$gtf)||die $!;
	open(BED,">$od/$key.bed")||die $!;

	while(<GTF>){
		chomp;next if($_=~/^#/);
		my @tmp=split(/\t/,$_);
		next if($tmp[2] ne "transcript");
		$tmp[8]=~/transcript_id "(.*?)"\;/;
		my $id=$1;
		print BED "chr$tmp[0]\t$tmp[3]\t$tmp[4]\t$id\n";
	}
	close(BED);
	close(GTF);
}

if(defined $gff){
	$gff=abs_path($gff);
	open(GFF,$gff)||die $!;
	open(BED,">$od/$key.bed")||die $!;
	while(<GFF>){
		chomp;next if($_=~/^#/);
		my @tmp=split(/\t/,$_);
		next if($tmp[8] =~ "Parent=");
		my $id=(split(/\;/,$tmp[8]))[0];
		$id=~s/ID=//;
		print BED "chr$tmp[0]\t$tmp[3]\t$tmp[4]\t$id\n";
	}
	close(BED);
	close(GFF);
}
if(defined $bed){
	$bed=abs_path($bed);
	my $flag=`grep ^chr $bed|wc -l`;
	chomp($flag);$flag=(split(/\s+/,$flag))[0];
	
	if($flag>0){
		`cp $bed $od/$key.bed`;
	}else{
		`sed 's/^/chr/' $bed >$od/$key.bed`;
	}
}
$bed="$od/$key.bed";
if(!-e $bed){
	print "Something must be wrong!\n";
	exit;
}

print "$extract $bw $bed $od/$key.PhastCons.score && sed -i \"1i#ID\\tsize\\tcovered\\tsum\\tmean0\\tmean\" $od/$key.PhastCons.score\n";
`$extract $bw $bed $od/$key.PhastCons.score && sed -i "1i#ID\\tsize\\tcovered\\tsum\\tmean0\\tmean" $od/$key.PhastCons.score`;
print "/share/nas2/genome/biosoft/R/3.1.1/bin/Rscript $Bin/phastcons_cdf.r $od/$key.PhastCons.score $od $key\n";
`/share/nas2/genome/biosoft/R/3.1.1/bin/Rscript $Bin/phastcons_cdf.r $od/$key.PhastCons.score $od $key`;
#name		- name field from bed, which should be unique
#size		- size of bed (sum of exon sizes
#covered	- # bases within exons covered by bigWig
#sum		- sum of values over all bases covered
#mean0		- average over bases with non-covered bases counting as zeroes
#mean		-average over just covered bases

#####sub
sub readConfig{
        my $configFile=shift;
        my $d=Config::General->new(-ConfigFile => "$configFile");
        my %config=$d->getall;  
        return %config;
}

sub USAGE {#
	my $usage=<<"USAGE";
Program: Get PhastCons from BigWig file downloaded from UCSC 
Version: $version
Contact: Wen yanhua

Description:
	This program is a Procedure deal with Get PhastCons from BigWig file downloaded from UCSC

Usage:
	-bw				BigWig file (some Species could be found in /share/nas1/wenyh/database/phastCons)
	-db				db(GRCh37,GRCh38,GRCm38,mm9,Rn_Celera,Rnor_5.0,Rnor_6.0)
					one of bw/db must be given.
	-gff                            gff file (would extract gene conservation) 
	-gtf				gtf file (would extract transcript conservation) 
	-bed				bed file (would extract region conservation given in bed)
					bed format(4 column): chr	start	end	name
					Only one of bed/gtf/gff file is needed!
	-key				output key word, default Region
	-od				output dir

 
USAGE
	print $usage;
	exit;
}
