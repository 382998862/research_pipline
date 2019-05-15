my($len,$out)=@ARGV;
if(@ARGV!=2){
	print "perl $0\n1.genome len file(the first two column should be :chr	length)\n2.output ideogram file for circos plot\n";
	exit;
}

my $num=30;
my @color=1..$num;
        open(LEN,$len)||die $!;
        open(OUT,">$out")||die $!;
        my %length=();
        while(<LEN>){
                chomp;
                next if($_=~/^#/);
                my @tmp=split(/\s+/,$_); #chr ACGTN ACGT N
                $length{$tmp[0]}=$tmp[1];
        }
	close(LEN);
        if(exists $length{1}){
		my @chrs=1..22;
		push @chrs,("X","Y","MT");
		my $total=0;
                foreach my $c(@chrs){
                        if(exists $length{$c}){
                                print OUT "chr\t-\t$c\t$c\t0\t$length{$c}\tchr$color[$total]\n";
                                $total++;
                        }
                }
        }else{
                my @keys=sort{$length{$a} <=> $length{$b}} keys %length;
                my $max=$num;
                $max=@keys      if(@keys<$max);
                for(my $i=1;$i<=$max;$i++){
                        print OUT "chr\t-\t$keys[-$i]\t$keys[-$i]\t0\t$length{$keys[-$i]}\tchr$color[-$i]\n";
                }
        }
        close(OUT);

