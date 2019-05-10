#!/usr/bin/env perl
use strict;
my $usage="
	perl $0 1.xml 2.xml 3.xml... > merged_Trinity.xml
	or
	perl $0 *.xml > merged_Trinity.xml
	  #Note: please provide all the xml
";

# 查看是否给了所有可读文件
die $usage unless @ARGV;
for (@ARGV){
	die $usage."$_ is unreadable \n" unless (-s $_);
}

############
#	打印blast header section
##############
open HEADER,'<',$ARGV[0] or die "Can't read $ARGV[0]\n";
while(<HEADER>){
	last if /<Iteration>/;
	print;
}
close HEADER;

############
#	打印blast hits
#################

# 用于更新iteration 次数
my $count;
# 用于记录已有blast结果的记录
my %seen;
for (@ARGV){
	my $label_pattern='<Iteration>';	
	# 每次打开一个文件
	my $ref=&Single_Label_Reader($_,$label_pattern);
	#开始读取数据,每次处理一个记录,不耗费内存
	my $content;
	while($content=&$ref){
		my ($label,@content)=@$content;
		# 只处理完整的Blast记录
		next unless grep {/<\/Iteration>/} @content;
		# 只处理没处理过的记录
		my ($cur_id)=map {/Iteration_query-def>(.*)<\/Iteration_query-def/} @content;
		next if $seen{$cur_id}++;
		
		# 修饰一些数字记录，确保与一次比对产生的更为接近
		++$count;
		@content = map {s/(\s+<Iteration_iter-num>)\d+(<\/Iteration_iter-num>)/$1$count$2/gr;} @content;
		@content = map {s/(\s+<Iteration_query-ID>Query_)\d+(<\/Iteration_query-ID>)/$1$count$2/gr;} @content;
		# 去除可能存在的tail 因为可能每一个都完整
		my $num=0;
		# 查看末尾3行，是否需要去除
		for (-1,-2,-3){
			if ($content[$_]=~/<\/BlastOutput_iterations/){
				$num=0-$_;
				last;
			};
		}
		pop @content for (1..$num);
		
		# 打印完整的blast记录
		print $label,@content;
	}
}

############
#	打印blast tail
#################
print "</BlastOutput_iterations>
</BlastOutput>

";


######
#	子函数
###########
sub Single_Label_Reader {
	my $file=shift;
	my $label_pattern=shift;
	die "Error:$file unreachable.\n" unless (-s $file);
	open IN,'<',$file or die "Error:$file unreadable.\n";
	#开始处理各记录
	my $temp1;	
	# 去除标签出现前的内容 方便某类内容合并
	while($temp1=<IN>){
		last if $temp1=~/$label_pattern/;
	}
	my $temp2=$temp1;
	my @temp;
	return sub {
		return '' if (eof(IN));
		while(!eof(IN)){
			$temp1=$temp2;
			@temp=();
			while(<IN>){
				if (/$label_pattern/){
					$temp2=$_;		 
					last;
				}
				push @temp,$_;
			}
			return [$temp1,@temp];
		}
	}
}