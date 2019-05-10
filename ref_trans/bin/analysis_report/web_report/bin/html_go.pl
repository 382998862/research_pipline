#!/usr/bin/perl -w
# 
# Copyright (c) BIO_MARK 2012
# Writer:         Zhangyh <zhangyh@biomarker.com.cn>
# Program Date:   2012/9/3 9:25
# Modifier:       Zhangyh <zhangyh@biomarker.com.cn>
# Last Modified:  2012/9/3 13:56
my $ver="1.0.0";
sub help {#usage
	print <<"	Usage End.";
	
	Creat GO annotation webpage

	perl  html_go.pl  -i input  -o output  -k Gb1,Gb2
	
	-i           input dir
	-o           output dir
	-k           key words   example: -k Gb1,Gb2 
	-h           Help document

	[notice] Input dir should contain:
	
	xxx.GO.Biological.stat	xxx.GO.Biological.xls
	xxx.GO.Cellular.stat	xxx.GO.Cellular.xls
	xxx.GO.Molecular.stat	xxx.GO.Molecular.xls

	Usage End.
	exit;
}


use Cwd;
use File::Basename qw(basename dirname);
use strict;
use Getopt::Long;
use Data::Dumper;

my %opts;
GetOptions(\%opts,
	"i=s",
	"o=s",
	"k=s",
	"h" 
);
&help() if(!defined($opts{i}) || !defined($opts{o}) || !defined($opts{k}) || defined($opts{h}));

#©¥©¥©¥©¥©¥©¥©¥©¥©¥©¥©¥©¥  Time  ©¥©¥©¥©¥©¥©¥©¥©¥©¥©¥©¥©¥
my $BEGIN = time();
my $Time_Start = sub_format_datetime(localtime(time()));
print "\nStart Time :[$Time_Start]\n\n";
#©¥©¥©¥©¥©¥©¥©¥©¥©¥©¥©¥©¥©¥©¥©¥©¥©¥©¥©¥©¥©¥©¥©¥©¥©¥©¥©¥©¥

my ($in,$out,
    $c_stat,$c_xls,
    $f_stat,$f_xls,
    $p_stat,$p_xls,);

&MKDIR($opts{o});
$in  = ABSOLUTE_DIR($opts{i});
$out = ABSOLUTE_DIR($opts{o});

$c_stat = glob("$in/*.Cellular.stat");
$c_xls  = glob("$in/*.Cellular.xls");
$p_stat = glob("$in/*.Biological.stat");
$p_xls  = glob("$in/*.Biological.xls");
$f_stat = glob("$in/*.Molecular.stat");
$f_xls  = glob("$in/*.Molecular.xls");
my $key = join ("vs",split /,/,$opts{k},-1);
my $type = ['Component','Function','Process'];
my $file = ['C','F','P'];

&main_html(\$c_stat,\$c_xls,\$type->[0],\$file->[0]);
&main_html(\$p_stat,\$p_xls,\$type->[2],\$file->[2]);
&main_html(\$f_stat,\$f_xls,\$type->[1],\$file->[1]);

sub main_html{
	my ($stat,$xls,$type,$f) = @_;
	my(@column,%column);
	open OUT, ">", "$out/${key}_$${f}1.html" || die "Can't creat $out/{$key}_$${f}1.html:$!\n";
	open STAT, "<", "$$stat" || die "Can't open $$stat:$!\n";
	open XLS, "<", "$$xls" || die "Can't open $$xls:$!\n";
	
	&html_begin($f);
	while (<STAT>) {
		next if /\#/;
		chomp;
		my ($term,$c_frequency,$g_frequency,undef,$corr_value) = split /\t/;
		&go_stat(\$term,\$c_frequency,\$g_frequency,\$corr_value,$f);
		push @column,$term;
	}
	while (<XLS>) {
		next if /\#/;
		chomp;
		my ($go_term,$gene) = (split /\t/)[0,5];
		$column{$go_term} = $gene;
	}
	&stat2xls();
	&go_xls(\$_,\$column{$_}) for @column;
	&html_end();
	
	close STAT;
	close XLS;
	close OUT;
}


#©¥©¥©¥©¥©¥©¥©¥©¥©¥©¥©¥©¥  Time  ©¥©¥©¥©¥©¥©¥©¥©¥©¥©¥©¥©¥
&Runtime($BEGIN);
my $Time_End = sub_format_datetime(localtime(time()));
print "\nEnd Time :[$Time_End]\n\n";



#©¥©¥©¥©¥©¥©¥©¥©¥©¥©¥©¥©¥©¥©¥©¥©¥©¥©¥©¥©¥©¥©¥©¥©¥©¥©¥©¥©¥©·
#                     Sub Routines                       ©§
#©¥©¥©¥©¥©¥©¥©¥©¥©¥©¥©¥©¥©¥©¥©¥©¥©¥©¥©¥©¥©¥©¥©¥©¥©¥©¥©¥©¥©¿
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
sub MKDIR{ # &MKDIR($out_dir);
	my ($dir)=@_;
	rmdir($dir) if(-d $dir);
	mkdir($dir) if(!-d $dir);
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

#©¥©¥©¥©¥©¥©¥©¥©¥©¥©¥©¥©¥©¥©¥©¥©¥©¥©¥©¥©¥©¥©¥©¥©¥©¥©¥©¥©¥©·
#                        HTML Sub                        ©§
#©¥©¥©¥©¥©¥©¥©¥©¥©¥©¥©¥©¥©¥©¥©¥©¥©¥©¥©¥©¥©¥©¥©¥©¥©¥©¥©¥©¥©¿
sub html_begin{
	my $f = shift @_;
	print OUT <<'HTML';
<?xml version="1.0" encoding="iso-8859-1"?>
<!DOCTYPE html PUBLIC "-//W3C//DTD XHTML 1.0 Transitional//EN" "http://www.w3.org/TR/xhtml1/DTD/xhtml1-transitional.dtd">
<html xmlns="http://www.w3.org/1999/xhtml" lang="en-US" xml:lang="en-US">
	<head>
HTML
	print OUT <<"HTML";
		<title>Terms for ${key}_$$f</title>
HTML
	print OUT <<'HTML';
		<script type="text/javascript">//<![CDATA[
			<!--
			url = "";
			//-->

			//]]></script>

			<style type="text/css">
			<!--
			body {background-color: #fff;}
			table {background-color: #000; border-collapse: collapse; border: solid #000 1px;}
			tr {background-color: #fff;}
			th, td {border: solid #000 1px;}
			a.applet {cursor: pointer; line-height: 50px; text-decoration: underline;}
			-->
			</style>
		<script type="text/javascript">
			<!--
			function reSize2() {
				if (parent.parent.document.getElementsByTagName("iframe")[0])
					parent.parent.document.getElementsByTagName("iframe")[0].style.height = document.body.scrollHeight + 50 + "px";
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

			function embedApplet(file, obj) {
				path = (navigator.userAgent.indexOf("IE") >= 0 && top.location != self.location)? "./GO/need/" : "";
				code = '<applet code="edu/stanford/genetics/treeview/applet/ButtonApplet.class" archive="' + path + 'TreeViewApplet.jar,' + path + 'nanoxml-2.2.2.jar,' + path + 'plugins/Dendrogram.jar" width="150" height="50" alt="Your browser understands the &lt;applet&gt; tag but isn\'t running the applet, for some reason.">' + "\n" +
			'Your browser is completely ignoring the &lt;applet&gt; tag!<br />Please install Java Runtime Environment (JRE) first, you can get it from <a href="http://www.java.com/en/download/manual.jsp">http://www.java.com/en/download/manual.jsp</a>' + "\n" +
			'<param name="cdtFile" value="' + path + file + '">' + "\n" +
			'<param name="cdtName" value="Result">' + "\n" +
			'<param name="plugins" value="edu.stanford.genetics.treeview.plugin.dendroview.DendrogramFactory">' + "\n" +
			'</applet>';
				obj.innerHTML = code;
			}

			function addLink(tdObj) {
				tdObj.onmouseover = "";
				if (url.match(/^\s*$/)) {
					return;
				}
				links = tdObj.innerHTML.split(/,\s+/);
				if (temp = links[0].match(/^(.*)>([\d]+)$/)) {
					links[0] = temp[1] + "><a href='" + url.replace("<REPLACE_THIS>", temp[2]) + "' target='_blank'>" + temp[2] + "</a>";
				}
				for (i = 1; i < links.length; i++) {
					if (links[i].match(/^[\d]+$/)) {
						links[i] = "<a href='" + url.replace("<REPLACE_THIS>", links[i]) + "' target='_blank'>" + links[i] + "</a>";
					}
				}
				tdObj.innerHTML = links.join(", ");
			}

			window.onload = function() {
				tableObj = document.getElementsByTagName("table")[1];
				tableObj2 = document.getElementsByTagName("table")[2];
				tableObj2.style.marginTop = "50px";
				trObj = tableObj.getElementsByTagName("tr");
				trObj2 = tableObj2.getElementsByTagName("tr");
				for (i = 1; i < trObj.length; i++) {
					trObj[i].cells[0].innerHTML += "&nbsp;<a title='click to view genes' href='#gene" + i + "' onclick='javascript: colorRow(document.getElementsByTagName(\"table\")[2].rows[" + i + "]);'>view genes</a>";
					trObj2[i].cells[1].onmouseover = function() {addLink(this);};
					trObj2[i].cells[1].innerHTML = "<a name='gene" + i + "'></a>" + trObj2[i].cells[1].innerHTML;
				}
				diffColor([document.getElementsByTagName("table")[1], document.getElementsByTagName("table")[2]]);
				setTimeout("reSize2()", 1);
			}
			//-->
		</script>
	</head>
	<body>
HTML
	print OUT <<"HTML";
		<center><h2>Terms for ${key}_$$f</h2></center>
HTML
	print OUT <<'HTML';
		<hr /><a name="table" />
		<center><h3>Result Table</h3></center>
		<p />
		<table border="1" align="center" cellpadding="2" width="400">
			<tr>
				<td bgcolor="#FFCC99" align="center" nowrap width="100%"><b>Terms from the Component Ontology with p-value as good or better than 1</b></td>
			</tr>
		</table>
		<table border="2" align="center">
			<tr bgcolor="#CCCCFF">
				<th align="center">Gene Ontology term</th><th align="center">Cluster frequency</th>
				<th align="center">Genome frequency of use</th>
				<th align="center">Corrected P-value</th>
			</tr>
HTML
}
sub go_stat{
	my ($term,$c_fre,$g_fre,$value,$f) = @_;
	my ($term_str,$term_id) = $$term =~ /(.+)\s\((.+)\)/;
	my ($c_str,$c_num) = $$c_fre=~/^(.+?)\s(\d+\.?\d*)%$/;

	my ($g_str,$g_num) = $$g_fre=~/^(.+?)\s(\d+\.?\d*)%$/;
	$c_num = sprintf ("%.2f",$c_num);
	$g_num = sprintf ("%.2f",$g_num);
	$value = ($$value == 1) ? '1' : sprintf ("%.2e",$$value);
	print OUT <<HTML;
			<tr>
				<td><a target="infowin" href="http://amigo.geneontology.org/cgi-bin/amigo/go.cgi?action=query&amp;view=query&amp;query=${term_id}&amp;search_constraint=terms">$term_str</a></td>
				<td>$c_str $c_num%</td>
				<td>$g_str $g_num%</td>
				<td>$value</td>
			</tr>
HTML
# <td><a onclick="javascript: embedApplet('cluster/${key}_$${f}_${term_str}.cdt', this);" class="applet">show applet</a></td>
}
sub stat2xls{
	print OUT <<HTML;
		</table>
		<table border="2" align="center">
			<tr>
				<th align="center">Gene Ontology term</th>
				<th align="center">Genes annotated to the term</th>
			</tr>
HTML
}
sub go_xls{
	my ($term,$gene) = @_;
	my ($term_str,$term_id) = $$term =~ /(.+)\s\((.+)\)/;
	$gene = join (", ",split/;/,$$gene);
	print OUT <<HTML
			<tr>
				<td><a target="infowin" href="http://amigo.geneontology.org/cgi-bin/amigo/go.cgi?action=query&amp;view=query&amp;query=${term_id}&amp;search_constraint=terms">${term_str}</a></td>
				<td>$gene</td>
			</tr>
HTML
}
sub html_end{
	print OUT <<HTML;
		</table>
	</body>
</html>
HTML
}
