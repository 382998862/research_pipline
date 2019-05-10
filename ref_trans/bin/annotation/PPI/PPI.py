#!/usr/bin/env python
#coding:utf-8

import sys,os
import argparse
import re

def runOrDie(cmd):
    try:
        print cmd
        os.system(cmd)
    except Exception,err:
        print "Error:"+err

if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("-fa",dest="fasta",type=str,help = "the query fa file")
    parser.add_argument("-cfg","--config",dest="config",type=str,help="the ref pipeline detail config")
    parser.add_argument("-of","--out_file",dest="outFile",type=str,help="the out file to store the ppi")
    parser.add_argument("-q","--queue",dest="queue",type=str,default="general.q",help = "the queue for blast")
    parser.add_argument("-p","--prefix",dest="prefix",type=str)

    args = parser.parse_args()

    ##check the options
    if  not args.config or not args.outFile:
        print parser.print_help()
        sys.exit("please supply the right parameters")

    ##set variables
    clade = ""   ##the string database clade
    kegg_to_clades = {"Animals.fa":"eukaryota","Archaea.fa":"archaea","Bacteria.fa":"bacteria","Fungi.fa":"eukaryota","Plants.fa":"eukaryota","Protists.fa":"eukaryota"}
    ##get
    if args.fasta:
        args.fasta = os.path.abspath(args.fasta)
    args.config = os.path.abspath(args.config)
    args.outFile = os.path.abspath(args.outFile)
    outDir = os.path.dirname(args.outFile)+"/ppi_tmp"
    if not os.path.exists(outDir):
        os.makedirs(outDir, 0755)

    ##
    ##read the config
    CFG = open(args.config,"r")
    for line in CFG.readlines():
        if  re.match(r"Known_anno",line):
            known_gene_anno_dir = re.split("\s+",line)[1]
        if re.match(r"Kegg",line):
            kegg_clade = os.path.basename(re.split("\s+",line)[1])
    CFG.close()
    if not kegg_to_clades.has_key(kegg_clade):
        sys.exc_info("Err:the kegg class from detail.cfg is not in the one-to-one table")
        sys.exit(1)
    clade = kegg_to_clades[kegg_clade]
    ##annotation the fasta file for ppi
    #first get the fasta name list
    if args.fasta:
        fasta = open(args.fasta,"r")
        fastaList = open(outDir+"/gene.list","w")
        for line in fasta.readlines():
            match = re.match(r">(\S+)", line)
            if match:
                fastaList.write(match.group(1)+"\n")
        fasta.close()
        fastaList.close()

        cmd = os.path.dirname("perl "+sys.argv[0])+"/"+"ppi_network.pl --qseq " +args.fasta + " --queue " + args.queue + " --clade " +clade+ " --odir " +outDir +" -index newGene > "+outDir+"/newGene.ppi.log 2>&1 "
        runOrDie(cmd)
        runOrDie("tail -n +2 "+outDir+"/newGene.ppi.detail.txt >"+outDir+"/newGene.ppi.txt.noHead")

    ##known ppi txt
    known_ppi = known_gene_anno_dir + "/ppi.txt"
    if not os.path.isfile(known_ppi):
        known_geneid_fa = known_gene_anno_dir + "/Known.longest_transcript.fa"
        known_geneid_fa_fh = open(known_geneid_fa, "r")
        known_geneid_list_fh = open(outDir + "/known_geneid_list", "w")
        for line in known_geneid_fa_fh.readlines():
            match = re.match(r">(\S+)", line)
            if match:
                known_geneid_list_fh.write(match.group(1) + "\n")
        known_geneid_list_fh.close()
        known_geneid_fa_fh.close()

        cmd = "perl " + os.path.dirname(sys.argv[0]) + "/" + "ppi_network.pl --qseq " + known_geneid_fa + " --queue " + args.queue + " --clade " + clade + " --odir " + outDir + " -index "+args.prefix+".knownGene > " + outDir + "/knownGene.ppi.log 2>&1 "
        runOrDie(cmd)
        if args.fasta:
            os.system("cat "+outDir +"/"+args.prefix +".knownGene.ppi.detail.txt "+outDir+"/newGene.ppi.detail.txt.noHead >"+args.outFile)
        else:
            os.system("cp "+outDir+"/"+args.prefix+".knownGene.ppi.detail.txt "+args.outFile)
    else:
        if args.fasta:
            os.system("cat "+known_ppi +" "+outDir+"/newGene.ppi.txt.noHead >"+args.outFile)
        else:
            os.system("cp " + outDir + "/knownGene.ppi.detail.txt " + args.outFile)



