use Getopt::Long;
use Data::Dumper;
use FindBin qw($Bin $Script);
use File::Basename qw(basename dirname);
my $BEGIN_TIME=time();
my $version="1.0.0";
my $bigWigSummary = "/share/nas2/genome/biosoft/kentUtils-master/bin/bigWigSummary";
#######################################################################################
# ------------------------------------------------------------------
# GetOptions, Check format, Output info.
# ------------------------------------------------------------------
my ($bw, $gff,$type,$od);
GetOptions(
				"help|?" =>\&USAGE,
				"bw:s"  =>\$bw,
				"circ:s"  =>\$circ,
				"type:s"   =>\$type,
				"gff:s"=>\$gff,
				"od:s"   =>\$od,
				) or &USAGE;
&USAGE unless ($bw and $circ and $od and $gff) ;
mkdir $od if(!-e $od);
$type ||="mean";
open(IN, "$circ") or die("Could not open input circular result file!\n"); 
open(OUTPUT, ">$od/conservation.result") or die("Could not open output file!\n");
my %genes;
while(<IN>)
{
	chomp;	
	next if(/circRNA_ID/);
	my @line = split/\t/,$_;
	my $trans = $line[0];
	my $chr = "chr".(split /:/,$line[0])[0];
	my ($start,$end) = split /\|/,(split /:/,$line[0])[1],2; 
	my $len= $end -$start;
	my $gene = $line[-1];
	$gene=~s/,//g;
	my @gene = split /;/,$gene;	
	foreach my $g(@gene){$genes{$g}=1};
	system("$bigWigSummary -type=$type $bw $chr $start $end 1 >$od/$$.$type.tmp 2>&1");
	open(TEMP, "$od/$$.$type.tmp") or die("Cannot open result temp file!\n");
	my $score = '';
	while(<TEMP>) {
		$score .= $_;
		chomp($score);
	}
	close(TEMP);
	if($score =~ /^no/) { 
		$score = 0;
	}
	print OUTPUT "$chr\tcircRNA\t$start\t$end\t$len\t$trans\t$score\n";
	system("rm $od/$$.$type.tmp");
}
close IN;


open GFF,"$gff";
while(<GFF>)
{
	chomp;
	my @line = split/\t/,$_;
	my $gene = $line[-1];
	$gene=~s/ID=//;
	$gene=~s/;$//;
	if($line[2] eq 'gene' && exists $genes{$gene})
	{
		my($chr,$start,$end)=("chr".$line[0],$line[3],$line[4]);
		my $len= $end -$start;
		system("$bigWigSummary -type=$type $bw $chr $start $end 1 >$od/$$.$type.tmp 2>&1");
        	open(TEMP, "$od/$$.$type.tmp") or die("Cannot open result temp file!\n");
        	my $score = '';
        	while(<TEMP>) {
                	$score .= $_;
                	chomp($score);
        	}
        	close(TEMP);
        	if($score =~ /^no/) {
               		$score = 0;
       		 }	
        	print OUTPUT "$chr\tgene\t$start\t$end\t$len\t$gene\t$score\n";
       		system("rm $od/$$.$type.tmp");
	}
}
close GFF;
close OUTPUT;
print "################################################\n";
print "Draw boxplot of circRNA scorephastcons\n";
if (-e "$od/scores.txt") {
	`rm $od/scores.txt`;
}
my $cmd1=qx(echo "Type    Score" >>$od/scores.txt);
my $cmd3 = qx(cut -f 2,7 $od/conservation.result >>$od/scores.txt);
#my $cmd4= qx(grep -v \'intergenic region\' $od/scores.txt > $od/scores.txt1);
my $cmd5 = "/share/nas2/genome/biosoft/R/3.1.1/lib64/R/bin/Rscript $Bin/R/phastcons_cdf.r $od/scores.txt $od ";
$cmd5 .= ">$od/plot_box.log 2>&1";
`$cmd5`;

#####################################################
#####sub

sub cmd_call {
	my ($shell,$cpu,$vf,$queue)=@_;
	system "sh /share/nas2/genome/bmksoft/tool/qsub_sge_plus/v1.0/qsub_sge.plus.sh --queue $queue --maxproc $cpu --resource vf=$vf --reqsub --independent $shell";
}


sub USAGE {#
	my $usage=<<"USAGE";
Program: Get PhastCons from BigWig file downloaded from UCSC 
Version: $version
Contact: Zhang Qiuxue <zhangqx\@biomarker.com.cn>

Description:
	This program is a Procedure deal with Get PhastCons from BigWig file downloaded from UCSC

Usage:
	-bw				BigWig file 
	-circ				circular result file
	-gff                            gff file 
	-od				output file 
	-type				phastcons type default mean;
 
USAGE
	print $usage;
	exit;
}
