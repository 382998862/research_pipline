use strict;
use Cwd;
use Getopt::Long;
use Data::Dumper;
use FindBin qw($Bin $Script);
use File::Basename qw(basename dirname);
my ($genome,$od,$win,$circRNA,$geneExpression,$species);
GetOptions(
	"genome:s"=>\$genome,
	"c:s"=>\$circRNA,
	"e:s"=>\$geneExpression,
	"od:s"=>\$od,
	"win:s"=>\$win,
	"s:s"=>\$species,
	"h|?"=>\&help,
)or &help;
&help unless ($genome and $od and $circRNA and $geneExpression);
$win =$win || 1000000;
$genome=&ABSOLUTE_DIR($genome);
&MKDIR($od);
$od=&ABSOLUTE_DIR($od);
my $dir = "$od/circosResult";
&MKDIR($dir);
my %position;
if ($genome !~/.txt$/ && !-e $genome) {
	print "karyotype file is not exists! please check the file or input a genome file to generate karyotype file\n";
	exit;
}
elsif($genome =~/.txt$/) {
	#### /share/nas2/genome/biosoft/circos/current/data/karyotype/
	open(IN,$genome) or die $!;
	while (<IN>)
	{
		chomp;
		last if/^band/;
		my @line=split /\s/,$_;
		my $chr=$line[3];
		$chr =~s/chr//i;
		$chr=uc($chr);
		my $length=$line[5];
		$position{$chr}=$length;
	}
	close IN;
	print "cp $genome $od/karyotype.txt\n";
	`cp $genome $od/karyotype.txt`;
}elsif ($genome =~/.fa$/) {
	open(IN,$genome) or die $!;
	open(OUT,">$od/karyotype.txt.tmp") or die $!;
	$/='>';
	<IN>;
	while (<IN>)
	{
		chomp;
		my ($head,$seq)=split /\n/,$_,2;
		my $chr=(split /\s+/,$head)[0];
		$chr =~s/chr//i;
		next if($chr!~/^\d+$/ and $chr!~/^X$/ and $chr!~/^Y$/);
		$seq=~s/\n//g;
		my @len=split //,$seq;
		my $length=@len;
		## file format:chr - hs1 1 0 249250621 chr1
		print OUT "chr - $species"."$chr $chr 0 $length chr".$chr."\n" if(defined $species);
		print OUT "chr - $chr $chr 0 $length chr".$chr."\n" unless(defined $species);
		$position{$chr}=$length;
	}
	$/="\n";
close IN;
close OUT;
`less $od/karyotype.txt.tmp|sort -n -k3,3 >$od/karyotype.txt && rm $od/karyotype.txt.tmp`;
}

##########################################
open IN,"$circRNA" or die $!;
my %cposition;
while (<IN>) {
	next if (/\#/);
	my @line = split /\t/, $_;
	my $circRNA_id = $line[0];
	my $chr = (split /:/,$circRNA_id)[0];
	$chr =~s/chr//i;
	my $start_end = (split /:/,$circRNA_id)[1];
	my $start =  (split /\|/,$start_end)[0];
	my $end =  (split /\|/,$start_end)[1];
	$cposition{$chr}{$start}{$end}=1;
}
close IN;
#####################################################
open OUT,">$od/circRNA.position.density.list" or die $!;
my %stat;

foreach my $chr(sort keys %position)
{
	if(exists $cposition{$chr})
	{
		my $index = 0;
		foreach my $start(sort {$a <=> $b} keys %{$cposition{$chr}})
		{
			for(my $i=$index;$i<=$position{$chr};$i+=$win)
			{
				$index=$i,last if($start>=$i && $start<=$i+$win);
			}
			foreach my $end(sort {$a <=> $b}keys %{$cposition{$chr}{$start}})
			{
				$stat{$chr}{$index}{$index+$win}++ if($end<$index+$win);
				$stat{$chr}{$index}{$index+$win}++ if($end>$index+$win);
			}
		}
	}
}


foreach my $chr(sort {$a <=> $b}keys%stat)
{
	foreach my $start(sort {$a <=> $b}keys %{$stat{$chr}})
	{
		foreach my $end(sort {$a <=> $b}keys %{$stat{$chr}{$start}})
		{
			print OUT "$species"."$chr $start $end $stat{$chr}{$start}{$end}\n" if(defined $species);
			print OUT "$chr $start $end $stat{$chr}{$start}{$end}\n"	unless(defined $species);
			
		}
	}
}
close OUT;
#########################################
###circ_DEG_Analysis/All_gene_fpkm.list
open IN,"$geneExpression" or die $!;
my $head = <IN>;
chomp $head;
my @line = split/\t/,$head;
shift @line;
close IN;
my $cmd ="";
for(my $i=0;$i<scalar(@line);$i++)
{
	my $col = $i+4;
	$cmd="touch $od/$line[$i].exp && less -S $geneExpression|grep -v \"#\"|sed 's/chr//ig'|awk -F \":\" '{print \$1\" \"\$2}'|awk -F \"|\" '{print \$1\" \"\$2\" \"\$3}'|awk -F \" \" '{print \"$species\"\$1\" \"\$2\" \"\$3\" \"\$$col\" fpkm=\"\$$col}' > $od/$line[$i].exp" if(defined $species);
	$cmd="touch $od/$line[$i].exp && less -S $geneExpression|grep -v \"#\"|sed 's/chr//ig'|awk -F \":\" '{print \$1\" \"\$2}'|awk -F \"|\" '{print \$1\" \"\$2\" \"\$3}'|awk -F \" \" '{print \$1\" \"\$2\" \"\$3\" \"\$$col\" fpkm=\"\$$col}' > $od/$line[$i].exp" unless(defined $species);
	`$cmd`;
}
#################################################################
my $chr_info = `head -n 3 $od/karyotype.txt|cut -f3`;
my @chr = split /\n/,$chr_info;
`cp $Bin/conf/ideogram.conf $Bin/conf/ticks.conf  $dir && sed -i 's/<pairwise 1 2>/<pairwise $chr[1] $chr[2]>/g' $dir/ideogram.conf`;
my $gap = 1/((scalar(@line)+2)*2);
print "$gap\n";
my $r1= 1.00-$gap;
my $r0= 1.00-$gap*2;
open(CONF,">$dir/circos.conf") or die $!;
print CONF <<END;
karyotype=$od/karyotype.txt
chromosomes_units=$win
chromosomes_display_default = yes
<<include ideogram.conf>>
<<include ticks.conf>>
<image>
#angle_offset* = -150
<<include etc/image.conf>>
</image>
<plots>
<plot>
type= histogram
r1= ${r1}r
r0= ${r0}r
file= $od/circRNA.position.density.list
color=vdgreen
thickness = 10p
<axes>
show= data
<axis>
color= vlgrey
spacing = 0.1r
thickness = 1
</axis>
<axis>
color     = grey
spacing   = 0.2r
thickness = 1
</axis>
</axes>
</plot>
END
my @expfiles = glob("$od/*.exp");
foreach my $exp (@expfiles)
{
	$r1= $r0-$gap;
	$r0= $r1-$gap;
	print CONF <<END;
<plot>
type       = heatmap
r1         = ${r1}r
r0         = ${r0}r
file       = $exp
<rules>
use       = yes
<rule>
condition  = (var(fpkm) >= 0 && var(fpkm) <= 1)
color = blue
</rule>
<rule>
condition  = (var(fpkm) >1 && var(fpkm) <= 10)
color = green
</rule>
<rule>
condition  = (var(fpkm) >10 && var(fpkm) <= 100)
color = orange
</rule>
<rule>
condition  = var(fpkm) >100
color = red
</rule>
</rules>
</plot>
END
}
print CONF <<END;
</plots>
<<include etc/housekeeping.conf>>
<<include etc/colors_fonts_patterns.conf>>
END

`/share/nas2/genome/biosoft/circos/current/bin/circos -conf $dir/circos.conf -outputdir $dir -outputfile circos`;

################################## sub funcions
sub MKDIR
{ # &MKDIR($out_dir);
	my ($dir)=@_;
	rmdir($dir) if(-d $dir);
	mkdir($dir) if(!-d $dir);
}

sub ABSOLUTE_DIR
{ #$pavfile=&ABSOLUTE_DIR($pavfile);
	my $cur_dir=`pwd`;chomp($cur_dir);
	my ($in)=@_;
	my $return="";

	if(-f $in)
	{
		my $dir=dirname($in);
		my $file=basename($in);
		chdir $dir;$dir=`pwd`;chomp $dir;
		$return="$dir/$file";
	}
	elsif(-d $in)
	{
		chdir $in;$return=`pwd`;chomp $return;
	}
	else
	{
		warn "Warning just for file and dir in [sub ABSOLUTE_DIR]\n";
		exit;
	}

	chdir $cur_dir;
	return $return;
}


sub max{#&max(lists or arry);
	#ÇóÁÐ±íÖÐµÄ×î´óÖµ
	my $max=shift;
	my $temp;
	while (@_) {
		$temp=shift;
		$max=$max>$temp?$max:$temp;
	}
	return $max;
}

################################################################################################################

sub min{#&min(lists or arry);
	#ÇóÁÐ±íÖÐµÄ×îÐ¡Öµ
	my $min=shift;
	my $temp;
	while (@_) {
		$temp=shift;
		$min=$min<$temp?$min:$temp;
	}
	return $min;
}

sub help
{
	my $usage =print <<"USAGE";
	Description:sefla
		Function : use draw the circos picture;
		Usage    :
		-genome
		genome fa or Ref_Genome/genome_size.txt ;must be given;
		-od
		out put dir ;must be given;
		-c
		circRNA position file;must be given;
		-e
		circRNA expression file;must be given;
		-win
		stastic unit of chrome ;default 1000000;choise
		-h
		Help document
	for example:
		perl draw_circos_circRNA.pl -genome genome.fa -c All_circRNA.list -e All.All_gene_fpkm.list -od ./circos -win 1000000
USAGE
	print $usage;
	exit;
}
