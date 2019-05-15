#!/usr/bin/perl -w
#Writer           songmm <songmm@biomarker.com.cn>
my $version="1.0.0";
my $BEGIN=time();

use strict;
use warnings;
use experimental 'smartmatch';
use Getopt::Long;
use Data::Dumper;
use FindBin qw($Bin $Script);
use File::Basename qw(basename dirname);
use Cwd 'abs_path';
my $programe_dir=basename($0);
my $path=dirname($0);
#######################################################################################
# ------------------------------------------------------------------
# GetOptions, Check format, Output info.
# ------------------------------------------------------------------
my ($gff,$circexplorer,$od,$ref);
GetOptions(
	"help|?" =>\&USAGE,
	"gff:s"=>\$gff,
	"c:s"=>\$circexplorer,
	"ref:s"=>\$ref,
	"od:s"=>\$od,
	) or &USAGE;
&USAGE unless ($gff and $circexplorer and $od and $ref) ;
mkdir $od if (!-d"$od");
$od=abs_path($od);
my $prefix = basename($circexplorer);
my %gffTable;
my %uniq;
open(GFF,"sort -k4n,4 -k5n,5 $gff|") or die $!;
while (<GFF>) {
    chomp;
    next if(/^#/);
    my @info = split /\t/,$_;
    my $chr = $info[0];
    #$chr=~s/chr//ig;
    my $type = $info[2];
    my $start = $info[3];
    my $end = $info[4];
    my @arr=($start,$end);
    if($type eq 'exon' && !exists $uniq{$chr}{"$start;$end"})
    {
        push @{$gffTable{$chr}},[@arr];
        $uniq{$chr}{"$start;$end"}=1;
    }
}
close GFF;
my %hash_ref;
$/=">";
open(REF, "$ref");
while (<REF>) {
    chomp;
    next if(/^\s*$/);
    my @line = split /\n/,$_,2;
    $line[0] =~/(\S+)\s*/;
    my $chr = $1;
    #$chr=~s/chr//ig;
    my $seq = (join "",(split/\n/,$line[1]));
    $hash_ref{$chr}=$seq;
}
close REF;
$/="\n";
open(IN, $circexplorer) or die $!;
my @head = split /\t/,<IN>;
open(STA, ">$od/$prefix.stat");
open COUNT,">$od/All_gene_counts.list";
open DETAIL,">$od/All_gene_counts_detail.xls";
open(LEN, ">$od/${prefix}_length.stat");
print LEN "circRNA_type\tgeneLength\n";
#circRNA_ID      T01     T02     circRNA_type    gene_id	strand
pop @head;
my $id = pop @head;
pop @head;
print DETAIL join("\t",@head)."\tgeneLength\t$id\n";
print COUNT join("\t",@head)."\tgeneLength\n";
open(SEQ, ">$od/$prefix.fa");
while (<IN>) {
    chomp;
    next if($.==1);
    my @info = split /\t/,$_;
    my @circRNA_id = split /:/,$info[0];
    my $chr = $circRNA_id[0];
    #$chr=~s/chr//ig;
    my($start,$end)=(split /\|/,$circRNA_id[1])[0..1];
    my $strand = pop @info;
    my $gene_id=pop @info;
    my $type =pop @info; 
    my @starts;
    my @ends;
    if(!exists $gffTable{$chr})
    {	
		 print LEN "$type\t".($end-$start+1)."\n";
	         if ($hash_ref{$chr}) {
        	 	my $seq = substr($hash_ref{$chr},$start-1,$end-$start+1);
			if($strand eq '-')
			{
				$seq =~ tr/ATCG/TAGC/;
                                $seq = reverse $seq;
			}
              		print SEQ ">$chr:$start|$end\n$seq\n";
	   	 	print COUNT join("\t",@info)."\t".($end-$start+1)."\n";
			print DETAIL join("\t",@info)."\t".($end-$start+1)."\t$gene_id\n";
  	        
		}
	next;
    }
    for(my $i=0;$i<@{$gffTable{$chr}};$i++){
        push @starts,$gffTable{$chr}[$i][0];
        push @ends,$gffTable{$chr}[$i][1];
    }
    my ($index1,$index2);
    if ($type eq 'intron') {
        print STA "$type-->$chr:$start|$end\t";
        print STA "\n";
        print LEN "$type\t".($end-$start+1)."\n";
        if ($hash_ref{$chr}) {
            my $seq = substr($hash_ref{$chr},$start-1,$end-$start+1);
	    if($strand eq '-')
            {
                $seq =~ tr/ATCG/TAGC/;
                $seq = reverse $seq;
            }

            print SEQ ">$chr:$start|$end\n$seq\n";
	    print COUNT join("\t",@info)."\t".($end-$start+1)."\n";
            print DETAIL join("\t",@info)."\t".($end-$start+1)."\t$gene_id\n";
        }
        
    }
    elsif($type eq 'exon')
    {
        print STA "$type--->$chr:$start|$end\n";
        my $length=0;
        my $seq="";
        TEST:for(my $i=0;$i<scalar(@starts);$i++)
        {
            my $current_start="";
            my $current_end="";
            if ($start>=$starts[$i] && $start<$ends[$i]) {
                my $s=$start;
                my $e=$ends[$i];
                LINE:while ($end>$e) {
                    print STA "$s-->$e;";
                    if($current_start eq '' && $current_end eq '')
                    {
                        $current_start=$s;
                        $current_end=$e;
                    }
                    else
                    {
                        if($s<=$current_end && $e>$current_end){$s = $current_end;}
                        elsif($s<=$current_end && $e<=$current_end)
                        {
                            $current_start=$s;
                            $i++;
                            $s=$starts[$i];
                            $e=$ends[$i];
                            print STA "(real:not count);";
                            next LINE;
                        }
                    }
                    $length+=$e-$s+1;
                    if ($hash_ref{$chr}) {
                        $seq .= substr($hash_ref{$chr},$s-1,$e-$s+1);
                    }
                    print STA "(real:$s-->$e);";
                    $current_start=$s;
                    $current_end=$e;
                    $i++;
                    $s=$starts[$i];
                    $e=$ends[$i];
                }
                if($end<=$e)
                {
                	if($end>$s)
                	{
                		print STA "$s-->$e;";
				print STA "(real:$s-->$end);";
				$length+=$end-$s+1;
				if ($hash_ref{$chr}) {
                           		$seq .= substr($hash_ref{$chr},$s-1,$end-$s+1);
                        	}

                	}
                	else
                	{
                		print STA "$s-->$e;";
                		print STA "(real:$ends[$i-1]-->$end);";
                		$length+=$end-$ends[$i-1]+1;
				if ($hash_ref{$chr}) {
                                        $seq .= substr($hash_ref{$chr},$ends[$i-1]-1,$end-$ends[$i-1]+1);
                                }
                	}
		}
                print STA "\n";
                print LEN "$type\t$length\n";
		if($strand eq '-')
                {
                     $seq =~ tr/ATCG/TAGC/;
                     $seq = reverse $seq;
                }
	        print SEQ ">$chr:$start|$end\n$seq\n";
		print COUNT join("\t",@info)."\t$length\n";
                print DETAIL join("\t",@info)."\t$length\t$gene_id\n";
                last TEST;
            }
            elsif($start<$starts[$i])
            {
                my $s=$start;
                my $e=$ends[$i];
                LINE1:while ($end>$e) {
                    print STA "$s-->$e;";
                    if($current_start eq '' && $current_end eq '')
                    {
                        $current_start=$s;
                        $current_end=$e;
                    }
                    else
                    {
                        if($s<=$current_end && $e>$current_end){$s = $current_end;}
                        elsif($s<=$current_end && $e<=$current_end)
                        {
                            $current_start=$s;
                            $i++;
                            $s=$starts[$i];
                            $e=$ends[$i];
                            print STA "(real:not count);";
                            next LINE1;
                        }
                    }
                    $length+=$e-$s+1;
                    if ($hash_ref{$chr}) {
                        $seq .= substr($hash_ref{$chr},$s-1,$e-$s+1);
                    }
                    print STA "(real:$s-->$e);";
                    $current_start=$s;
                    $current_end=$e;
                    $i++;
                    $s=$starts[$i];
                    $e=$ends[$i];
                }
                if($end<=$e)
                {
                	if($end>$s)
                	{
                		print STA "$s-->$e;";
				print STA "(real:$s-->$end);";
				$length+=$end-$s+1;
				if ($hash_ref{$chr}) {
                                        $seq .= substr($hash_ref{$chr},$s-1,$end-$s+1);
                                }
	
                	}
                	else
                	{
                		print STA "$s-->$e;";
                		print STA "(real:$ends[$i-1]-->$end);";
                		$length+=$end-$ends[$i-1]+1;
				if ($hash_ref{$chr}) {
                                        $seq .= substr($hash_ref{$chr},$ends[$i-1]-1,$end-$ends[$i-1]+1);
                                }

                	}

		}
                print STA "\n";
                print LEN "$type\t$length\n";
		if($strand eq '-')
                {
                    $seq =~ tr/ATCG/TAGC/;
                    $seq = reverse $seq;
                }

                print SEQ ">$chr:$start|$end\n$seq\n";
		print COUNT join("\t",@info)."\t$length\n";
                print DETAIL join("\t",@info)."\t$length\t$gene_id\n";
                last TEST;
            }
        }

    }
else
    {
        print LEN "$type\t".($end-$start+1)."\n";
        if ($hash_ref{$chr}) {
            my $seq = substr($hash_ref{$chr},$start-1,$end-$start+1);
	    if($strand eq '-')
            {
                $seq =~ tr/ATCG/TAGC/;
                $seq = reverse $seq;
            }

            print SEQ ">$chr:$start|$end\n$seq\n";
	    print COUNT join("\t",@info)."\t".($end-$start+1)."\n";
            print DETAIL join("\t",@info)."\t".($end-$start+1)."\t$gene_id\n";
        }
    }
}
close IN;
close STA;
close LEN;
close SEQ;
close COUNT;
close DETAIL;
=pod
sub find
{
    my($arr,$num)=@_;
    my ($ta,$tb) = (0,(scalar @$arr) - 1);
    while($ta<$tb)
    {
        my $mid = int( ($tb+$ta) / 2 );
        return $mid if ($$arr[$mid] == $num);
        if($$arr[$mid] > $num)
        {
            $tb = $mid;
        }else
        {
            $ta = $mid;
        }
    }
}
=cut 
sub find
{
    my($arr,$num)=@_;
    for(my $i=0;$i<scalar($arr);$i++)
    {
        return $i if($$arr[$i]==$num);
    }
}
sub find_scope
{
    my($start,$end,$num)=@_;
    for(my $i=0;$i<scalar(@$start)-1;$i++)
    {
        if(($num>$$start[$i] && $num<$$end[$i])||($num>$$end[$i] && $num<$$start[$i+1])){return $i;}
	#elsif($num==$$arr[$i+1]){return ($i+1);}
	#elsif($num>$$arr[$i] && $num<$$arr[$i+1]){return $i;}
    }
}

sub USAGE
{
    	my $usage=<<"USAGE";
Program: extract intron or exon circRNA
Version: $version
Contact: songmm <songmm\@biomarker.com.cn>

Description:

Usage:
  -gff  <file>  gff file must be given;
  -c    <file>  circRNA file must be given;
  -ref  <file>   genome file must be given;
  -od   <dir>   out directory must be given;
USAGE
	print $usage;
	exit;
}
