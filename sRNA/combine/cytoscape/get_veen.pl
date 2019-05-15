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
my ($idir,$odir,$cfg,$l2mc,$l2mt,$mi2l,$mi2c,$mi2g,$chostgene,$l2m);

GetOptions(
				"help|?" =>\&USAGE,
				"od:s"=>\$odir,
				"in:s"=>\$idir,
				"cfg:s"=>\$cfg,
				"l2mc:s"=>\$l2mc,
				"l2mt:s"=>\$l2mt,
				"l2m:s"=>\$l2m,
				"mi2l:s"=>\$mi2l,
				"mi2c:s"=>\$mi2c,
				"mi2g:s"=>\$mi2g,
				"chostgene:s"=>\$chostgene,
				) or &USAGE;
&USAGE unless ($idir and $odir and $cfg);
# ------------------------------------------------------------------
$odir||="./";
$cfg = abs_path($cfg);
$idir = abs_path($idir);
$odir = abs_path($odir);
`mkdir $odir` unless (-d "$odir");
$l2mc||="$idir/Cis_target_gene.xls";
$l2mt||="$idir/Trans_target_gene.xls";
$l2m||="$idir/noval_target_gene.xls";
$mi2l||="$idir/lncRNA.mir2target.list";
$mi2c||="$idir/circRNA.mir2target.list";
$mi2g||="$idir/gene.mir2target.list";
$chostgene||="$idir/circRNA.source.xls";

my %sample = &relate($cfg);
my @diff;my %vs;
open (CFG,"$cfg") or die $!;
print ("read cfgi and diff group is:");
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

my (%l2m,%m2l,%mi2l,%l2mi,%mi2c,%c2mi,%mi2g,%g2mi,%c2g,%g2c);

###read lncRNA2mRNA file and then get %l2m and %m2l
print ("1:L2M...... %l2m,%m2l\n");
if ((defined $l2m) && (-e "$l2m")){
	open (L2M,"$l2m");
	while(<L2M>){
		chomp;
		next if(/^#|^$/);
		my @A = split /\t/,$_;
		my @B = split /\;|,/,$A[1];
		my $nm = join (";",@B);
		$l2m{$A[0]}=$nm;
		for (my $i=0;$i<@B;$i++){
			if(exists $m2l{$B[$i]}){
				$m2l{$B[$i]} .= ";".$A[0];
			}else{
				$m2l{$B[$i]} = $A[0];
			}
		}
	}
	close (L2M);
}else {
	open (L2MC,"$l2mc");
	while(<L2MC>){
		chomp;
		next if(/^#|^\s$/);
		my @A = split /\t/,$_;
		my @B = split /\;|,/,$A[1];
		my $nm = join (";",@B);
		$l2m{$A[0]}{Cis}=$nm;
		for (my $i=0;$i<@B;$i++){
			if(exists $m2l{$B[$i]}{Cis}){
				$m2l{$B[$i]}{Cis} .= ";".$A[0];
			}else{
				$m2l{$B[$i]}{Cis} = $A[0];
			}
		}
	}
	if(-e $l2mt){
		open (L2MT,"$l2mt");
		while (<L2MT>){
			next if(/^#|^\s$/);
			my @A = split /\t/,$_;
	                my @B = split /\;|,/,$A[1];
                	my $nm = join (";",@B);
        	        $l2m{$A[0]}{Trans}=$nm;
	                for (my $i=0;$i<@B;$i++){
                	        if(exists $m2l{$B[$i]}{Trans}){
        	                        $m2l{$B[$i]}{Trans} .= ";".$A[0];
	                        }else{
	                               	$m2l{$B[$i]}{Trans} = $A[0];
                        	}
                	}
		}
	}
}


###read miRNA2lncRNA file and then get %mi2l and %l2mi
if(-e "$mi2l"){
	print "2:MI2L...... %mi2l,%l2mi\n";
	open (MI2L,"$mi2l");
	while(<MI2L>){
		chomp;
		next if(/^#|^$/);
		my @A = split /\t/,$_;
		my @B = split /\;|,/,$A[1];
		my $nm = join (";",@B);
		$mi2l{$A[0]}=$nm;
		for (my $i=0;$i<@B;$i++){
			if(exists $l2mi{$B[$i]}){
				$l2mi{$B[$i]} .= ";".$A[0];
			}else{
				$l2mi{$B[$i]} = $A[0];
			}
		}
	}
	close (MI2L);
}
###read miRNA2circRNA file

if(-e "$mi2c"){
	print ("3:MI2C...... %mi2c,%c2mi\n");
	open(MI2C,"$mi2c");
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
	open (MI2G,"$mi2g");
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
	open (C2G,"$chostgene");
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
	foreach my $k ("sRNA","lncRNA","circRNA","gene"){
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
###################lncRNA core#############################
print "core-lncRNA\n";
`mkdir $odir/lncRNA` unless (-d "$odir/lncRNA");
my $lncdir = "$odir/lncRNA";
if((defined $l2m) && (-e "$l2m")){
	foreach my $lnc (keys %h){
        open (IN,"$h{$lnc}{lncRNA}");
        my (@A,@B,@C,%tmp,%tmp1);
        while(<IN>){
                chomp;
                next if(/^#/);
                my @a = split;
                push @A,$a[0];
        }
	close (IN);
        open (IN,"$h{$lnc}{gene}");
        while(<IN>){
                chomp;
                next if(/^#/);
                my @a = split;
                if(exists $m2l{$a[0]}){
                        my @t = split /\;/,$m2l{$a[0]};
                        for (my $i=0;$i<@t;$i++){
                                if(!exists $tmp{$t[$i]}){
                                        push @B,$t[$i];
                                        $tmp{$t[$i]}=1;
                                }
                        }
                }
        }
	close (IN);
        open(IN,"$h{$lnc}{sRNA}");
        while(<IN>){
                chomp;
                next if(/^#/);
                my @a = split;
                if(exists $mi2l{$a[0]}){
                        my @t = split /\;/,$mi2l{$a[0]};
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
	if($#A>=1){push @list1,@A;push @name1,"DE_LncRNA";}
	if($#B>=1){push @list1,@B;push @name1,"DE_mRNA_TargetLncRNA"};
	if($#C>=1){push @list1,@C;push @name1,"DE_miRNA_TargetLncRNA"};
	print "$#A\t$#B\t$#C\n";
        #my @list1 = ("@A","@B","@C");
        #my @name1 = ("DE_LncRNA","DE_mRNA_TargetLncRNA","DE_miRNA_TargetLncRNA");
        if(@list1>=2){
        	&Veen(\@list1,"$lncdir/$lnc",\@name1);
	}
	}
}else{
	foreach my $lnc (keys %h){
	open (IN,"$h{$lnc}{lncRNA}");
	my (@A,@B,@C,@D,%tmp,%tmp1,%tmp2);
	while(<IN>){
		chomp;
		next if(/^#/);
		my @a = split;
		push @A,$a[0];
	}
	close (IN);
	open (IN,"$h{$lnc}{gene}");
	while(<IN>){
		chomp;
		next if(/^#/);
		my @a = split;
		if(exists $m2l{$a[0]}{Cis}){
			my @t = split /\;/,$m2l{$a[0]}{Cis};
			for (my $i=0;$i<@t;$i++){
				if(!exists $tmp{$t[$i]}){
					push @B,$t[$i];
					$tmp{$t[$i]}=1;
				}
			}
		}
		if (exists $m2l{$a[0]}{Trans}){
			my @t = split /\;/,$m2l{$a[0]}{Trans};
			for (my $i=0;$i<@t;$i++){
				if(!exists $tmp1{$t[$i]}){
					push @C,$t[$i];
					$tmp1{$t[$i]}=1;
				}
			}
		}
	}
	close (IN);
	open(IN,"$h{$lnc}{sRNA}");
	while(<IN>){
		chomp;
		next if(/^#/);
		my @a = split;
		if(exists $mi2l{$a[0]}){
			my @t = split /\;/,$mi2l{$a[0]};
			for (my $i=0;$i<@t;$i++){
				if(!exists $tmp2{$t[$i]}){
					push @D,$t[$i];
					$tmp2{$t[$i]}=1;
				}
			}
		}
	}
	close(IN);
	my @list1 = ("@A","@B","@D");
	my @list2 = ("@A","@C","@D");
	my @name1 = ("DE_LncRNA","DE_Cis.mRNA_TargetLncRNA","DE_miRNA_TargetLncRNA");
	my @name2 = ("DE_LncRNA","DE_Trans.mRNA_TargetLncRNA","DE_miRNA_TargetLncRNA");
	&Veen(\@list1,"$lncdir/$lnc.Cis",\@name1);
	&Veen(\@list2,"$lncdir/$lnc.Trans",\@name2);
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
        my @list1 = ("@A","@B","@C");
        my @name1 = ("DE_circRNA","DE_Hostgene_circRNA","DE_miRNA_TargetcircRNA");
        &Veen(\@list1,"$circdir/$circ",\@name1);
}

print "core-mRNA\n";
`mkdir $odir/gene` unless (-d "$odir/gene");
my $gdir = "$odir/gene";
if((defined $l2m) && (-e "$l2m")){
	foreach my $g (keys %h){
        open (IN,"$h{$g}{gene}");
        my (@A,@B,@C,@D,%tmp,%tmp1,%tmp2);
        while(<IN>){
                chomp;
                next if(/^#/);
                my @a = split;
                push @A,$a[0];
        }
	close (IN);
        open (IN,"$h{$g}{lncRNA}");
        while(<IN>){
                chomp;
                next if(/^#/);
                my @a = split;
                if(exists $l2m{$a[0]}){
                        my @t = split /\;/,$l2m{$a[0]};
                        for (my $i=0;$i<@t;$i++){
                                if(!exists $tmp{$t[$i]}){
                                        push @B,$t[$i];
                                        $tmp{$t[$i]}=1;
                                }
                        }
                }
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
                                if(!exists $tmp1{$t[$i]}){
                                        push @C,$t[$i];
                                        $tmp1{$t[$i]}=1;
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
                                if(!exists $tmp2{$t[$i]}){
                                        push @D,$t[$i];
                                        $tmp2{$t[$i]}=1;
                                }
                        }
                }
        }
        close(IN);
        my @list1 = ("@A","@B","@C","@D");
        my @name1 = ("DE_mRNA","DE_LncRNA_TargetmRNA","DE_miRNA_TargetmRNA","DE_circRNA_Hostgene");
        &Veen(\@list1,"$gdir/$g",\@name1);
	}
}
else{
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
        open (IN,"$h{$g}{lncRNA}");
        while(<IN>){
                chomp;
                next if(/^#/);
                my @a = split;
                if(exists $l2m{$a[0]}{Cis}){
                        my @t = split /\;/,$l2m{$a[0]}{Cis};
                        for (my $i=0;$i<@t;$i++){
                                if(!exists $tmp{$t[$i]}){
                                        push @B,$t[$i];
                                        $tmp{$t[$i]}=1;
                                }
                        }
                }
		if(exists $l2m{$a[0]}{Trans}){
                        my @t = split /\;/,$l2m{$a[0]}{Trans};
                        for (my $i=0;$i<@t;$i++){
                                if(!exists $tmp1{$t[$i]}){
                                        push @C,$t[$i];
                                        $tmp1{$t[$i]}=1;
                                }
                        }
                }
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
        my @list1 = ("@A","@B","@D","@E");
	my @list2 = ("@A","@C","@D","@E");
        my @name1 = ("DE_mRNA","DE_LncRNA_Target.CismRNA","DE_miRNA_TargetmRNA","DE_circRNA_Hostgene");
	my @name2 = ("DE_mRNA","DE_LncRNA_Target.TransmRNA","DE_miRNA_TargetmRNA","DE_circRNA_Hostgene");
        &Veen(\@list1,"$gdir/$g.Cis",\@name1);
	&Veen(\@list2,"$gdir/$g.Trans",\@name2);
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
        open(IN,"$h{$mi}{lncRNA}");
        while(<IN>){
                chomp;
                next if(/^#/);
                my @a = split;
                if(exists $l2mi{$a[0]}){
                        my @t = split /\;/,$l2mi{$a[0]};
                        for (my $i=0;$i<@t;$i++){
                                if(!exists $tmp1{$t[$i]}){
                                        push @C,$t[$i];
                                        $tmp1{$t[$i]}=1;
                                }
                        }
                }
        }
        close(IN);
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
        my @list1 = ("@A","@B","@C","@D");
        my @name1 = ("DE_miRNA","DE_mRNA_TargetmiRNA","DE_LncRNA_TargetmiRNA","DE_circRNA_TargetmiRNA");
        &Veen(\@list1,"$sdir/$mi",\@name1);
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
        return %sample; #relation{W01}{lncRNA}="L01";
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

sub getCyto {
	my $inter = shift;
	my $att = shift;
	my $line = shift;
	my $Cyto = shift;
	open (O1,">$inter");
	open (O2,">$att");
	my %ind;
	my ($center,$type)=split //,$line;
	foreach my $rna(sort keys %$Cyto){
		foreach my $stat (sort keys %{$$Cyto{$rna}}){
			if (!exists $ind{$rna}){
				print O2 "$rna\t$center\t$stat\n";
				$ind{$rna}=1;
			}
			my @p = split /\;/,$$Cyto{$rna}{$stat};
			for (my $i=0;$i<@p;$i++){
				my ($id,$ud,$k) = split /\(|\)/,$p[$i];
				print O1 "$rna\t$id\t$line\n";
				if (!exists $ind{$id}){
					print O2 "$id\t$type\t$ud\n";
					$ind{$id}=1;
				}
			}
		}
	}
	close(O1);close(O2);
	return ($inter,$att);
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
	for my $i (1..@list){
		my @ids;
		
		for my $deg ($list[$i-1]){
			@ids = split /\s/,$deg;
			for (my $j=0;$j<@ids;$j++){
				my $deg_id = $ids[$j];
				$info{$deg_id}{$label[$i-1]} = 1;
			}
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
	for my $i (1..@list) {
		my @id;
		for my $deg ($list[$i-1]){
			@id = split /\s/,$deg;
			for(my $j=0;$j<@id;$j++){
				my $deg_id = "'".$id[$j]."'";
			push @{$DEG{$label[$i-1]}}, $deg_id;
			}
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

	"cd $dir && /share/nas2/genome/biosoft/R/3.1.1/bin/Rscript $name.venn.r";
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
		-l2m <file>	lncRNA2mRNA target file, l2m can exist with (l2mc and l2mt) at the same time.
		-l2mc <file>	Cis target file, (l2mc and l2mt) can exist with l2m at the same time.
		-l2mt <file>	Trans target file, l2mt can not exist 
		-mi2l <file>	miRNA2lncrNA target file,
		-mi2c <file>	miRNA2circRNA target file,
		-mi2g <file>	miRNA2gene target file,
		-chostgene <file>	circRNA hostgene file
		-h		help

USAGE
	print $usage;
	exit;
}
