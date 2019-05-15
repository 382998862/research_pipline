#!/usr/bin/env perl
use strict;
use warnings;
use Cwd qw(abs_path);
use Getopt::Long;
use Data::Dumper;
use FindBin qw($Bin $Script);
use File::Basename qw(basename dirname);
use newPerlBase;
my $BEGIN_TIME=time();
my $version="1.0.0";
#######################################################################################
my ($idir,$odir,$cfg,$mi2g,$nl);

GetOptions(
				"help|?" =>\&USAGE,
				"od:s"=>\$odir,
				"in:s"=>\$idir,
				"cfg:s"=>\$cfg,
				"mi2g:s"=>\$mi2g,
				"name2list:s"=>\$nl,
				) or &USAGE;
&USAGE unless ($idir and $odir and $cfg);
# ------------------------------------------------------------------
$odir||="./";
$cfg = abs_path($cfg);
$idir = abs_path($idir);
$odir = abs_path($odir);
`mkdir $odir` unless (-d "$odir");
$mi2g||="$idir/gene.mir2target.list";

my %name;
if(defined $nl){
	open (IN,$nl);
	while(<IN>){
		chomp;
		next if(/^#/);
		my @t = split /\t/,$_;
		$name{$t[0]}=$t[1];
		my $up = $t[1]."(up)";
		my $down = $t[1]."(down)";
		$name{"$t[0](up)"}=$up;
		$name{"$t[0](down)"}=$down;
	}
	close(IN);
}

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

my (%mi2g,%g2mi,);
###read miRNA2gene file
if(-e $mi2g){
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
print ("read target file progress finish!!!\n\n");

print "get diff result file and creat hash\n\n";
print "#group\tRNA_type\tdiff_files\n";
my (%h,%hash);
foreach my $vs(@diff){
	&switch($vs);
	foreach my $k ("sRNA","gene"){
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
print "core-miRNA\n";
`mkdir $odir/sRNA` unless (-d "$odir/sRNA");
my $sdir = "$odir/sRNA";
foreach my $mi (keys %h){
        open (IN,"$h{$mi}{sRNA}");
        open (OUT,">$sdir/$mi.DEG.xls");
        my $sginter = "$sdir/$mi.sm.Interaction.list";
        my $sgatt = "$sdir/$mi.sm.Attribution.list";
        my (%sgCyto);
        while(<IN>){
                chomp;
                if(/^#/){
			if(defined $nl){
				print OUT "$_\tdiff_target_gene\tSymbol\n";
			}else{
				print OUT "$_\tdiff_target_gene\n";
			}
			next;
		}
                my @A = split /\t/,$_;
                print OUT "$_\t";
                if(exists $mi2g{$A[0]}){
                        my @sg = split /\;/,$mi2g{$A[0]};
                        my $ginfo;my $flag =0;my $ninfo;
                        for (my $i=0;$i<@sg;$i++){
                                if(exists $hash{$mi}{gene}{$sg[$i]}){
                                        $flag++;
					if(defined $nl){
                                          if(exists $name{$hash{$mi}{gene}{$sg[$i]}}){
                                                if($flag==1){
                                                        $ginfo = $hash{$mi}{gene}{$sg[$i]};
                                                        $ninfo = $name{$hash{$mi}{gene}{$sg[$i]}};
                                                }else{
                                                        $ginfo .= ";".$hash{$mi}{gene}{$sg[$i]};
                                                        $ninfo .= ";".$name{$hash{$mi}{gene}{$sg[$i]}};
                                                }
                                          }else{
                                                if($flag==1){
                                                        $ginfo = $hash{$mi}{gene}{$sg[$i]};
                                                        $ninfo = $hash{$mi}{gene}{$sg[$i]};
                                                }else{
                                                        $ginfo .= ";".$hash{$mi}{gene}{$sg[$i]};
                                                        $ninfo .= ";".$hash{$mi}{gene}{$sg[$i]};
                                                }
                                          }
					}else{
						if($flag==1){
							$ginfo = $hash{$mi}{gene}{$sg[$i]};
						}else{
							$ginfo .= ";".$hash{$mi}{gene}{$sg[$i]};
						}
					}
                                }
                        }
                        if($flag >0){
				if(defined $nl){
	                                print OUT "$ginfo\t$ninfo\n";
				}else{
					print OUT "$ginfo\n";
				}
                                $sgCyto{$A[0]}{$A[-1]}=$ginfo;
                        }else{
				if(defined $nl){
	                                print OUT "--\t--\n";
				}else{
					print OUT "--\n";
				}
                        }
                }else{
			if(defined $nl){
	                        print OUT "--\t--\n";
			}else{
				print OUT "--\n";
			}
                }
        }
        close (IN);
        close (OUT);
        &getCyto($sginter,$sgatt,"sm",\%sgCyto,"1");
}

print "core-gene\n";
`mkdir $odir/gene` unless (-d "$odir/gene");
my $gdir = "$odir/gene";
foreach my $g (keys %h){
        open (IN,"$h{$g}{gene}");
        open (OUT,">$gdir/$g.DEG.xls");
	my $gsinter = "$gdir/$g.ms.Interaction.list";
        my $gsatt = "$gdir/$g.ms.Attribution.list";
        my (%gsCyto);
        while(<IN>){
                chomp;
                if(/^#/){
			my ($i,$d)=split /\t/,$_,2;
			if(defined $nl){
				print OUT "$i\tSymbol\t$d\tdiff_target_miRNA\n";
			}else{
				print OUT "$i\t$d\tdiff_target_miRNA\n";
			}
			next;
		}
                my @A = split /\t/,$_;
		my ($g,$fpkm)=split /\t/,$_,2;
		if(defined $nl){
		  if(exists $name{$g}){
			print OUT "$g\t$name{$g}\t$fpkm\t";
		  }else{
			print OUT "$g\t$g\t$fpkm\t";
		  }
		}else{
			print OUT "$g\t$fpkm\t";
		}
                if(exists $g2mi{$A[0]}){
                        my @gs = split /\;/,$g2mi{$A[0]};
                        my $sinfo;my $flag =0;
                        for (my $i=0;$i<@gs;$i++){
                                if(exists $hash{$g}{sRNA}{$gs[$i]}){
                                        $flag++;
                                        if($flag==1){
                                                $sinfo = $hash{$g}{sRNA}{$gs[$i]};
                                        }else{
                                                $sinfo .= ";".$hash{$g}{sRNA}{$gs[$i]};
                                        }
                                }
                        }
                        if($flag >0){
                                print OUT "$sinfo\n";
				$gsCyto{$A[0]}{$A[-1]}=$sinfo;
                        }else{
                                print OUT "--\n";
                        }
                }else{
                        print OUT "--\n";
                }
	}
        close (IN);
        close (OUT);
        &getCyto($gsinter,$gsatt,"ms",\%gsCyto,"1");
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
	my $flag = shift;
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
				if($flag==1){
					my ($id,$ud,$k) = split /\(|\)/,$p[$i];
					print O1 "$rna\t$id\t$line\n";
					if (!exists $ind{$id}){
						print O2 "$id\t$type\t$ud\n";
						$ind{$id}=1;
					}
				}else{
					print O1 "$rna\t$p[$i]\t$line\n";
					if(!exists $ind{$p[$i]}){
						print O2 "$p[$i]\t$type\ttarget\n";
					}
				}
			}
		}
	}
	close(O1);close(O2);
	return ($inter,$att);
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
		-mi2g <file>	miRNA2gene target file,
		-name2list <file>	name2list file
		-h		help

USAGE
	print $usage;
	exit;
}
