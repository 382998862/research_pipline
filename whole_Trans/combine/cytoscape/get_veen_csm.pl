#!/usr/bin/env perl
use strict;
use warnings;
use Getopt::Long;
use Data::Dumper;
use FindBin qw($Bin $Script);
use File::Basename qw(basename dirname);
use newPerlBase;
use Cwd qw(abs_path getcwd);
my $BEGIN_TIME=time();
my $version="1.0.0";
#######################################################################################
#
#######################################################################################
# ------------------------------------------------------------------
# GetOptions
# ------------------------------------------------------------------
my ($idir,$odir,$cfg,$mi2c,$mi2g,$chostgene);

GetOptions(
				"help|?" =>\&USAGE,
				"od:s"=>\$odir,
				"in:s"=>\$idir,
				"cfg:s"=>\$cfg,
				"mi2c:s"=>\$mi2c,
				"mi2g:s"=>\$mi2g,
				"chostgene:s"=>\$chostgene,
				) or &USAGE;
&USAGE unless ($idir and $cfg);
# ------------------------------------------------------------------
$odir||="./";
$cfg = abs_path($cfg);
$idir = abs_path($idir);
`mkdir $odir` unless (-d "$odir");
$odir = abs_path($odir);

$mi2c||="$idir/circRNA.mir2target.list";
$mi2g||="$idir/gene.mir2target.list";
$chostgene||="$idir/circRNA.source.xls";

my %sample = &relate($cfg);
my @diff;my %vs;
open (CFG,"$cfg") or die $!;
print ("read cfg and diff group is:");
while (<CFG>){
	chomp;
	next if(/^#|^$|^\s$/);
	my @tmp = split /\s+/,$_;
	if($tmp[0] eq "Diff"){
		push @diff,$tmp[1];
		print "$tmp[1]\n";
	}
}
close(CFG);

my (%mi2c,%c2mi,%mi2g,%g2mi,%c2g,%g2c);

###read miRNA2circRNA file

if(-e "$mi2c"){
	print ("3:MI2C...... %mi2c,%c2mi\n");
	open(MI2C,"$mi2c") or die $!;
	while(<MI2C>){
		chomp;
		next if(/^#|^$/);
		my @A = split /\t/,$_;
		my @B = split /\;|,/,$A[1];
		my $nm = join (";",@B);
		$mi2c{$A[0]}=$nm;
		for (my $i=0;$i<@B;$i++){
			if(exists $c2mi{$B[$i]}){
				$c2mi{$B[$i]} .= ";".$A[0];
			}else{
				$c2mi{$B[$i]} = $A[0];
			}
		}
	}
	close (MI2C);
}
###read miRNA2gene file
if(-e "$mi2g"){
	print ("4:MI2G...... %mi2g,%g2mi\n");
	open (MI2G,"$mi2g") or die $!;
	while(<MI2G>){
		chomp;
		next if(/^#|^$/);
		my @A = split /\t/,$_;
		my @B = split /\;|,/,$A[1];
		my $nm = join (";",@B);
		$mi2g{$A[0]}=$nm;
		for (my $i=0;$i<@B;$i++){
			if(exists $g2mi{$B[$i]}){
				$g2mi{$B[$i]} .= ";".$A[0];
			}else{
				$g2mi{$B[$i]} =$A[0];
			}
		}
	}
	close (MI2G);
}
###read circRNA hostgene file
if(-e "$chostgene"){
	print ("5:C2G...... %c2g,%g2c\n");
	open (C2G,"$chostgene") or die $!;
	while (<C2G>){
		chomp;
		next if(/^#|^$/);
		my @A = split /\t/,$_;
		my @B = split /\;|,/,$A[1];
		my $nm = join (";",@B);
		$c2g{$A[0]}=$nm;
		for (my $i=0;$i<@B;$i++){
			if(exists $g2c{$B[$i]}){
				$g2c{$B[$i]} .= ";".$A[0];
			}else{
				$g2c{$B[$i]} = $A[0];
			}
		}
	}
	close(C2G);
}
print ("read target file progress finish!!!\n\n");

print "get diff result file and creat hash\n\n";
print "#group\tRNA_type\tdiff_files\n";
my (%h,%hash);
foreach my $vs(@diff){
	&switch($vs);
	foreach my $k ("sRNA","circRNA","gene"){
		$h{$vs}{$k} = "$idir/$k.$vs{$vs}{$k}.DEG_final.xls";
		print "$vs\t$k\t$h{$vs}{$k}\n";
	}
}

print "\n%h===>\$hash{\$deg}{\$trans}{\$id}=\$regulate\n\n";
foreach my $deg (keys %h){
	foreach my $trans(keys %{$h{$deg}}){
		open (IN,"$h{$deg}{$trans}");
		while (<IN>){
			chomp;
			next if(/^#/);
			my @a = split /\t/,$_;
			my $info = $a[0]."(".$a[-1].")";
			$hash{$deg}{$trans}{$a[0]}=$info;
		}
		close(IN);
	}
}

print "core-circRNA\n";
`mkdir $odir/circRNA` unless (-d "$odir/circRNA");
my $circdir = "$odir/circRNA";
foreach my $circ (keys %h){
        open (IN,"$h{$circ}{circRNA}");
        my (@A,@B,@C,%tmp,%tmp1);
        while(<IN>){
                chomp;
                next if(/^#/);
                my @a = split;
                push @A,$a[0];
        }
        close (IN);
        open (IN,"$h{$circ}{gene}");
        while(<IN>){
                chomp;
                next if(/^#/);
                my @a = split;
                if(exists $g2c{$a[0]}){
                        my @t = split /\;/,$g2c{$a[0]};
                        for (my $i=0;$i<@t;$i++){
                                if(!exists $tmp{$t[$i]}){
                                        push @B,$t[$i];
                                        $tmp{$t[$i]}=1;
                                }
                        }
                }
        }
        close (IN);
        open(IN,"$h{$circ}{sRNA}");
        while(<IN>){
                chomp;
                next if(/^#/);
                my @a = split;
                if(exists $mi2c{$a[0]}){
                        my @t = split /\;/,$mi2c{$a[0]};
                        for (my $i=0;$i<@t;$i++){
                                if(!exists $tmp1{$t[$i]}){
                                        push @C,$t[$i];
                                        $tmp1{$t[$i]}=1;
                                }
                        }
                }
        }
        close(IN);
	my (@list1,@name1);
	if(scalar(@A)>=1){push @list1,[@A];push @name1,"DE_circRNA";}
	if(scalar(@B)>=1){push @list1,[@B];push @name1,"DE_Hostgene_circRNA";}
	if(scalar(@C)>=1){push @list1,[@C];push @name1,"DE_miRNA_TargetcircRNA";}
	if(@list1>=2){
		&Veen(\@list1,"$circdir/$circ",\@name1);
	}
}

print "core-mRNA\n";
`mkdir $odir/gene` unless (-d "$odir/gene");
my $gdir = "$odir/gene";
foreach my $g (keys %h){
        open (IN,"$h{$g}{gene}");
        my (@A,@B,@C,@D,@E,%tmp,%tmp1,%tmp2,%tmp3);
        while(<IN>){
                chomp;
                next if(/^#/);
                my @a = split;
                push @A,$a[0];
        }
        close (IN);
        open(IN,"$h{$g}{sRNA}");
        while(<IN>){
                chomp;
                next if(/^#/);
                my @a = split;
                if(exists $mi2g{$a[0]}){
                        my @t = split /\;/,$mi2g{$a[0]};
                        for (my $i=0;$i<@t;$i++){
                                if(!exists $tmp2{$t[$i]}){
                                        push @D,$t[$i];
                                        $tmp2{$t[$i]}=1;
                                }
                        }
                }
        }
        close(IN);
	open(IN,"$h{$g}{circRNA}");
        while(<IN>){
                chomp;
                next if(/^#/);
                my @a = split;
                if(exists $c2g{$a[0]}){
                        my @t = split /\;/,$c2g{$a[0]};
                        for (my $i=0;$i<@t;$i++){
                                if(!exists $tmp3{$t[$i]}){
                                        push @E,$t[$i];
                                        $tmp3{$t[$i]}=1;
                                }
                        }
                }
        }
        close(IN);
	my (@list1,@name1);
	if(scalar(@A)>=1){push @list1,[@A];push @name1,"DE_mRNA";}
	if(scalar(@D)>=1){push @list1,[@D];push @name1,"DE_miRNA_TargetmRNA";}
	if(scalar(@E)>=1){push @list1,[@E];push @name1,"DE_circRNA_Hostgene";}
	if(scalar(@list1)>=2){
	        &Veen(\@list1,"$gdir/$g",\@name1);
	}
}
print "core-miRNA\n";
`mkdir $odir/sRNA` unless (-d "$odir/sRNA");
my $sdir = "$odir/sRNA";
foreach my $mi (keys %h){
        open (IN,"$h{$mi}{sRNA}");
        my (@A,@B,@C,@D,%tmp,%tmp1,%tmp2);
        while(<IN>){
                chomp;
                next if(/^#/);
                my @a = split;
                push @A,$a[0];
        }
        close (IN);
        open (IN,"$h{$mi}{gene}");
        while(<IN>){
                chomp;
                next if(/^#/);
                my @a = split;
                if(exists $g2mi{$a[0]}){
                        my @t = split /\;/,$g2mi{$a[0]};
                        for (my $i=0;$i<@t;$i++){
                                if(!exists $tmp{$t[$i]}){
                                        push @B,$t[$i];
                                        $tmp{$t[$i]}=1;
                                }
                        }
                }
        }
        close (IN);
        open(IN,"$h{$mi}{circRNA}");
        while(<IN>){
                chomp;
                next if(/^#/);
                my @a = split;
                if(exists $c2mi{$a[0]}){
                        my @t = split /\;/,$c2mi{$a[0]};
                        for (my $i=0;$i<@t;$i++){
                                if(!exists $tmp2{$t[$i]}){
                                        push @D,$t[$i];
                                        $tmp2{$t[$i]}=1;
                                }
                        }
                }
        }
        close(IN);
#        my @list1 = ("@A","@B","@C","@D");
#        my @name1 = ("DE_miRNA","DE_mRNA_TargetmiRNA","DE_LncRNA_TargetmiRNA","DE_circRNA_TargetmiRNA");
	my (@list1,@name1);
	if(scalar(@A)>=1){push @list1,[@A];push @name1,"DE_miRNA";}
	if(scalar(@B)>=1){push @list1,[@B];push @name1,"DE_mRNA_TargetmiRNA";}
	if(scalar(@D)>=1){push @list1,[@D];push @name1,"DE_circRNA_TargetmiRNA";}
	if(scalar(@list1)>=2){
		&Veen(\@list1,"$sdir/$mi",\@name1);
	}
}

#######################################################################################
print STDOUT "\nDone. Total elapsed time : ",time()-$BEGIN_TIME,"s\n";
#######################################################################################

# ------------------------------------------------------------------
# sub function
# ------------------------------------------------------------------
################################################################################################################
sub relate{
        my $cfg=shift;
        my %sample=();
        my @rnas=();
        open(CFG,$cfg)||die $!;
        while(<CFG>){
                chomp;next if($_!~/^Sample/);
                my @tmp=split(/\s+/,$_);shift @tmp;
                if($tmp[0] eq "ID" && scalar(@rnas)==0){
                        shift @tmp;@rnas=@tmp;
                }else{
                        my $id=shift @tmp;
                        for(my $i=0;$i<@tmp;$i++){
                                $sample{$id}{$rnas[$i]}=$tmp[$i];
                        }
                        $sample{$id}{gene}=  $sample{$id}{lncRNA} if(exists $sample{$id}{lncRNA} && !exists $sample{$id}{gene});
                }
        }
        close(CFG);
        return %sample;
}

sub switch {
	my $group = shift;
	my ($v1,$v2,@V1,@V2);
	
		($v1,$v2) = split /_vs_/,$group,2;
		@V1 = split /_/,$v1;
		@V2 = split /_/,$v2;
		foreach my $key ("sRNA","lncRNA","circRNA","gene"){
			my (@s1,@s2);
			for(my $i=0;$i<@V1;$i++){
				push @s1,$sample{$V1[$i]}{$key};
			}
			for(my $j=0;$j<@V2;$j++){
				push @s2,$sample{$V2[$j]}{$key};
			}
			my $S1 = join ("_",@s1);
			my $S2 = join ("_",@s2);
			my $Group = join("_vs_",$S1,$S2);
			$vs{$group}{$key}=$Group;
		}
}

sub Veen{
	my ($List,$name,$id)=@_;
	my $dir = dirname $name;
	my $prefix = basename $name;
	my (%info,%venny, %com);
	open (SET, ">$name.degset.xls") or die $!;
	print SET "#DEG_Set\tDEG_Num\tDEG_IDs\n";
	my @label =@$id;
	my @list =@$List;
	for my $i (1..($#list+1)){
		my @ids=@{$list[$1-1]};
		for (my $j=0;$j<@ids;$j++){
			my $deg_id = $ids[$j];
			$info{$deg_id}{$label[$i-1]} = 1;
		}
		my ($id_num,$ids) = ($#ids+1, (join ";",@ids));
		print SET "$label[$i-1]\t$id_num\t$ids\n";
	}

	close SET;
	for my $e (sort keys %info) {
		my $com = join ",",(sort keys %{$info{$e}});
		$com{$com}++;
		$venny{$com}{$e} = 1;
	}
	open (VENN, ">$name.vennset.xls");
	print VENN "#Venn_Set\tElement_Num\tElement_IDs\n";

	for my $s (sort keys %venny) {
		my $elements = join ";",(sort keys %{$venny{$s}});
		print VENN "$s\t$com{$s}\t$elements\n";
	}
	my @color=("'cornflowerblue'","'green'","'yellow'","'darkorchid1'","'red'");
	my %DEG;
	my $list_content;
	my $label_content;
	my $color_content;
	for my $i (1..($#list+1)) {
		my @id = @{$list[$i-1]};
		for(my $j=0;$j<@id;$j++){
			my $deg_id = "'".$id[$j]."'";
			push @{$DEG{$label[$i-1]}}, $deg_id;
		}
	}
	my @name=qw(A B C D E F G M);
	my $name_content;
	for my $i (0..@label-1) {
		$list_content.= "$name[$i] <- c(".(join ", ",@{$DEG{$label[$i]}}).")\n";
		$label_content.= "$name[$i] = $name[$i], ";
		$name_content.="\"$label[$i]\",";
		$color_content.= "$color[$i], ";
	}
	$list_content =~ s/\n$//;
	$label_content =~ s/, $//;
	$color_content =~ s/, $//;
	$name_content=~s/,$//;
	my $more_opts = "";
	if (@label == 5) {
		$more_opts.= "    cex = 0.5,\n";
		$more_opts.= "    cat.cex = 0.6,\n";
		$more_opts.= "    margin = 0.1,\n";
		$more_opts.= "    cat.dist = c(0.20, 0.25, 0.20, 0.20, 0.25),\n";
		$more_opts.= "    scaled = FALSE,\n";
	}
	elsif (@label == 4) {
		$more_opts.= "    cex = 0.6,\n";
		$more_opts.= "    cat.cex = 0.7,\n";
		$more_opts.= "    margin = 0.08,\n";
		$more_opts.= "    scaled = FALSE,\n";
	}
	elsif (@label == 3) {
		$more_opts.= "    cex = 0.7,\n";
		$more_opts.= "    cat.cex = 0.7,\n";
		$more_opts.= "    margin = 0.06,\n";
		$more_opts.= "    euler.d = FALSE,\n";
		$more_opts.= "    cat.pos = c(0,0,180),\n";
		$more_opts.= "    scaled = FALSE,\n";
	}
	elsif (@label == 2) {
		$more_opts.= "    cex = 0.8,\n";
		$more_opts.= "    cat.cex = 0.5,\n";
		$more_opts.= "    margin = 0.05,\n";
		$more_opts.= "    cat.pos = 180,\n";
		$more_opts.= "    euler.d = TRUE,\n";
		$more_opts.= "    scaled = FALSE,\n";
	}
my $R_script = << "EOF";
#!/share/nas2/genome/biosoft/R/3.1.1/lib64/R/bin/Rscript
#import library
library(grid)
library(VennDiagram)
#init deg lists
$list_content
lst=list($label_content)
names(lst)=c($name_content)
#plot venn
venn.diagram(
	x = lst,
	filename = "$prefix.venn.png",
	imagetype = "png",
	fill = c($color_content),
	height=1300,
	width=1300,
	resolution=200,
	units='px',
	wd=1,
	$more_opts
);
	
EOF
	open (RS, ">$name.venn.r") or die;
		print RS $R_script;
	close RS;

	`cd $dir && /share/nas2/genome/biosoft/R/3.1.1/bin/Rscript $name.venn.r`;
}

################################################################################################################

sub GetTime {
	my ($sec, $min, $hour, $day, $mon, $year, $wday, $yday, $isdst)=localtime(time());
	return sprintf("%4d-%02d-%02d %02d:%02d:%02d", $year+1900, $mon+1, $day, $hour, $min, $sec);
}

################################################################################################################
sub USAGE {
	my $usage=<<"USAGE";
 ProgramName:
     Version:	$version
     Contact:	Niepy <niepy\@biomarker.com.cn> 
Program Date:	2018.05.23
      Modify:	
 Description:	This program is used to .... 
       Usage:
		Options:
		-in <dir>	input dir ,DEG_Analysis,force
		-od <dir>	output dir , Combine/Cytoscape,force
		-cfg <file>	cfg,force
		-mi2c <file>	miRNA2circRNA target file,
		-mi2g <file>	miRNA2gene target file,
		-chostgene <file>	circRNA hostgene file
		-h		help

USAGE
	print $usage;
	exit;
}
