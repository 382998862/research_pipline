#!/usr/bin/perl -w
# 
# Copyright (c) BIO_MARK 2012
# Writer:         Zhangyh <zhangyh@biomarker.com.cn>
# Program Date:   2012/8/28 16:20
# Modifier:       Zhangyh <zhangyh@biomarker.com.cn>
# Last Modified:  2012/8/30 13:44

#━━━━━━━━━━━━━━━━━━━━━━━━
#       说明文档
#	根据富集分析结果生成网页
#	
#	备注：输入目录结构
#	
#		|-- kegg_enrichment
#		|   |-- xxx.KEGG.stat
#		|   `-- xxx.KEGG.xls
#		`-- kegg_map
#			|-- ko00010.html
#			|-- ko00010.png
#				...
#			|-- koxxxxx.html
#			|-- koxxxxx.png
#
#━━━━━━━━━━━━━━━━━━━━━━━━

my $ver="1.0.0";


use Cwd;
use File::Basename qw(basename dirname);
use strict;
use Getopt::Long;
use Data::Dumper;

my %opts;
GetOptions(\%opts,
	"i=s",
	"o=s",
	"h" 
);
&help() if(!defined($opts{i}) || !defined($opts{o}) || defined($opts{h}));
###############Time
my $BEGIN=time();
my $Time_Start = sub_format_datetime(localtime(time()));
print "\nStart Time :[$Time_Start]\n\n";
################

my ($in,$out,$stat,$xls,@queue,%hash);
INIT: {
	&MKDIR($opts{o});
	$in   = ABSOLUTE_DIR($opts{i});
	$out  = ABSOLUTE_DIR($opts{o});
	$stat = glob("$in/kegg_enrichment/*stat");
	$xls  = glob("$in/kegg_enrichment/*xls");
}
open (OUT,">$out/pathway.html") || die "Can't creat $out,$!\n";

HTML: {
	&html_head();
	&html_head_to_table1();
	open (STAT,$stat) || die "Can't open $stat:$!\n";
	<STAT>;
	my $i = 1;
	while (<STAT>) {
		chomp;
		if ($. == 2) {
			my @num = map {/of\s(\d+?)\s/} (split /\t/,$_,-1);
			my $title = [
				'#',
				'Pathway',
				"DEGs with pathway annotation ($num[0])",
				"All genes with pathway annotation ($num[1])",
				'p_value',
				'corr_p_value',
				'Pathway ID',
			];
			&format_stat_title($title);
			my @arr = split/\t/;
			my $row = &split_row(\@arr,$i);
			&format_stat_row($row);
			$i += 1;
			next;
		}
		my @arr = split/\t/;
		my $row = &split_row(\@arr,$i);
		&format_stat_row($row);
		$i += 1;
	}
	close (STAT);
	&html_table1_to_table2();
	open (XLS,$xls) || die "Can't open $xls:$!\n";
	&format_xls_title();
	<XLS>;
	<XLS>;
	while (<XLS>){
		$hash{(split/\t/)[0]}{gene_id} = (split/\t/)[6];
		$hash{(split/\t/)[0]}{ko_id}   = (split/\t/)[1];
	}
	my $j = 1;
	for (@queue) {
		&format_xls_row(\$_,\$hash{$_}{gene_id},\$hash{$_}{ko_id},\$j);
		$j += 1;
	}
	close (XLS);
	&html_table2_to_end();
}
close (OUT);

###############Time
my $Time_End;
&Runtime($BEGIN);
$Time_End = sub_format_datetime(localtime(time()));
print "\nEnd Time :[$Time_End]\n\n";


#━━━━━━━━━━━━━━━━━━━━━━━┓
#                  Sub Routines                ┃
#━━━━━━━━━━━━━━━━━━━━━━━┛
sub sub_format_datetime {#Time calculation subroutine
    my($sec, $min, $hour, $day, $mon, $year, $wday, $yday, $isdst) = @_;
	$wday = $yday = $isdst = 0;
    sprintf("%4d-%02d-%02d %02d:%02d:%02d", $year+1900, $mon+1, $day, $hour, $min, $sec);
}
sub Runtime{ # &Runtime($BEGIN);
	my ($t1)=@_;
	my $t=time()-$t1;
	print "\nTotal elapsed time: ${t}s\n";
}
sub ABSOLUTE_DIR{ #$pavfile=&ABSOLUTE_DIR($pavfile);
	my $cur_dir=`pwd`;chomp($cur_dir);
	my ($in)=@_;
	my $return="";
	if(-f $in){
		my $dir=dirname($in);
		my $file=basename($in);
		chdir $dir;$dir=`pwd`;chomp $dir;
		$return="$dir/$file";
	}elsif(-d $in){
		chdir $in;$return=`pwd`;chomp $return;
	}else{
		warn "Warning just for file and dir\n";
		exit;
	}
	chdir $cur_dir;
	return $return;
}
sub MKDIR{ # &MKDIR($out_dir);
	my ($dir)=@_;
	rmdir($dir) if(-d $dir);
	mkdir($dir) if(!-d $dir);
}

####################
#	sub for webpage
####################
sub html_head{
print OUT <<'________HTML________';
<html>
	<head>
		<meta http-equiv="content-type" content="text/html; charset=utf-8">
		<title>Sample1vsSample2</title>
		<style type="text/css">
		body {background-color: #fff;}
		table {background-color: #000; border-collapse: collapse; border: solid #000 1px; margin: 0 0 50px 0;}
		tr {background-color: #fff;}
		th, td {border: solid #000 1px;}
		</style>
		<script type="text/javascript">
		<!--
		function reSize2() {
			parent.document.getElementsByTagName("iframe")[0].style.height = document.body.scrollHeight + 10;
			parent.parent.document.getElementsByTagName("iframe")[0].style.height = parent.document.body.scrollHeight;
		}

		preRow = null;
		preColor = null;
		function colorRow(trObj) {
			if (preRow != null) {
				preRow.style.backgroundColor = preColor;
			}
			preRow = trObj;
			preColor = trObj.style.backgroundColor;
			trObj.style.backgroundColor = "#ff0";
		}

		function diffColor(tables) {
			color = ["#fff", "#bbf"];
			for (i = 0; i < tables.length; i++) {
				trObj = tables[i].getElementsByTagName("tr");
				for (j = 1; j < trObj.length; j++) {
					trObj[j].style.backgroundColor = color[j % color.length];
				}
			}
		}

		function showPer(tableObj) {
			trObj = tableObj.getElementsByTagName("tr");
			if (trObj.length < 2) {
				return;
			}
			sum1 = trObj[0].cells[2].innerHTML.replace(/^.*\(([\d]+)\).*$/, "$1");
			sum2 = trObj[0].cells[3].innerHTML.replace(/^.*\(([\d]+)\).*$/, "$1");
			trObj[0].cells[2].innerHTML = "DEGs with pathway annotation (" + sum1 + ")";
			trObj[0].cells[3].innerHTML = "All genes with pathway annotation (" + sum2 + ")";
			for (i = 1; i < trObj.length; i++) {
				trObj[i].cells[2].innerHTML += " (" + (Math.round(trObj[i].cells[2].innerHTML * 10000/ sum1) / 100) + "%)";
				trObj[i].cells[3].innerHTML += " (" + (Math.round(trObj[i].cells[3].innerHTML * 10000/ sum2) / 100) + "%)";
			}
		}

		window.onload = function() {
			setTimeout("reSize2()", 1);
		}
		//-->
		</script>
	</head>
________HTML________
}
sub html_head_to_table1{
	print OUT <<________HTML________;
	<body>
		<table>
			<caption style='font-weight: 900;'>1. Sample1vsSample2</caption>
________HTML________
}
sub format_stat_title{
	my $title = shift @_;
	print OUT "			<tr>\n";
	print OUT "				<th>$_</th>\n" for (@{$title});
	print OUT "			</tr>\n";
}
sub split_row{
		my @arr = @{$_[0]};
		my $i =$_[1];
		push @queue,$arr[0];
		$arr[2] =~ s/^(\d+?)\s.*/$1/;
		$arr[3] =~ s/^(\d+?)\s.*/$1/;
		$arr[4] =sprintf ("%.4e",$arr[4]);
		$arr[5] = ($arr[5] == 1) ? 1 : sprintf("%.4e",$arr[5]);
		my $row = [$i,@arr[0,2,3,4,5,1]];
		return $row;
}
sub format_stat_row{
	my $row = $_[0];
print OUT <<________HTML________;
			<tr>
				<td>${$row}[0]</td>
				<td><a href='#gene${$row}[0]' title='click to view genes' onclick='javascript: colorRow(document.getElementsByTagName(\"table\")[1].rows[${$row}[0]]);'>${$row}[1]</a></td>
________HTML________
	print OUT "				<td>$_</td>\n" for @{$row}[2..6];
	print OUT "			</tr>\n";
}
sub html_table1_to_table2{
	print OUT <<________HTML________;
		</table>
		<table>
________HTML________
}
sub format_xls_title{
	print OUT <<________HTML________;
			<tr>
				<th>#</th>
				<th>Pathway</th>
				<th>Differentially expressed genes</th>
			</tr>
________HTML________
}
sub format_xls_row{
	my ($term,$gene_id,$ko_id,$i) = @_;
	$gene_id = join ("; ",split/;/,$$gene_id);
	if (-f "$in/kegg_map/$$ko_id.html") {
		print OUT <<________HTML________;
				<tr>
					<td>$$i</td>
					<td><a href='kegg_map/$$ko_id.html' title='click to view map' target='_blank'>$$term</a></td>
					<td><a name='gene$$i'></a>$gene_id</td>
				</tr>
________HTML________
	} else {
		print OUT <<________HTML________;
				<tr>
					<td>$$i</td>
					<td>$$term</td>
					<td><a name='gene$$i'></a>$gene_id</td>
				</tr>
________HTML________
	}
}
sub html_table2_to_end{
	print OUT <<________HTML________;
		</table>
		<script type='text/javascript'>showPer(document.getElementsByTagName('table')[0]);
			diffColor([document.getElementsByTagName('table')[0], document.getElementsByTagName('table')[1]]);
		</script>
	</body>
</html>
________HTML________
}
sub help {#usage
	print <<"	Usage End.";
	Description:
		Creat Kegg_pathway web page
		Version: $ver

	Usage:

		-i  <str>  input dir

		-o  <str>  output dir

	[notice]:input dir should contain:

		|-- kegg_enrichment
		|   |-- Tea.KEGG.stat
		|   `-- Tea.KEGG.xls
		`-- kegg_map
			|-- ko00010.html
			|-- ko00010.png
				...
			|-- koxxxxx.html
			|-- koxxxxx.png
			
	Usage End.

	exit;
}