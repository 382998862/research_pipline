#!/usr/bin/env python

import sys
import os
import getopt
from Bio.Blast.Applications import NcbiblastxCommandline
from Bio.Blast import NCBIXML

# This script takes a file with fasta genes and uses NCBI BLAST+ to look for
# matches in the eggNOG BLAST database. It also uses eggNOG to assign functional
# categories to each fasta. Inputs are the directory where eggnog files are
# found and an input fasta.

def get_arguments(argv):
    if len(argv) == 0:
        usage()
        sys.exit(2)
    inFasta = None
    eggNOGdirectory = None
    blastcmd = None
    outdir = None
    try:
        opts, args = getopt.getopt(argv, "f:e:c:o:")
    except getopt.GetoptError:
        usage()
        sys.exit(2)
    for opt, arg in opts:
        if opt == '-f':
            inFasta = arg
        elif opt == '-e':
            eggNOGdirectory = arg
        elif opt == '-c':
            blastcmd = arg
        elif opt == '-o':
            outdir = arg
    return (inFasta, eggNOGdirectory,blastcmd,outdir)

def usage():
    print "eggNOGblast.py\n \
        -f <input fasta file>\n \
        -e <directory with eggNOG files>\n \
        -c <blastx path>\n \
        -o <directory with outfile >"
def blast_eggNOG(inFasta, eggNOGdirectory,blastcmd,outdir):
    blastx_cline = NcbiblastxCommandline(query = inFasta, cmd=blastcmd,task ='blastx-fast',
                    db = eggNOGdirectory + '/' + "eggnogv4.db",
                    evalue=1e-5,max_target_seqs = 5, outfmt = 5, num_threads = 2,
                    out = outdir + '/' + os.path.split(inFasta)[1] + ".xml")
    print(blastx_cline)
    stdout, stderr = blastx_cline()

def parse_blast(inFasta, eggNOGdirectory,outdir):
    result_handle = open( outdir + '/' + os.path.split(inFasta)[1] + '.xml')
    outfile = open( outdir + '/' + os.path.split(inFasta)[1] + '_blast.txt', 'w')
    outfile.write("Query\tMatch\tE-value\tStart\tEnd\tHSP_score\n")
    blast_records = NCBIXML.parse(result_handle)
    for rec in blast_records:
        for align in rec.alignments:
            for hsp in align.hsps:
                outfile.write("%s\t%s\t%f\t%i\t%i\t%i\n" \
                % (rec.query, align.hit_def, hsp.expect, 
                hsp.query_start, hsp.query_end,hsp.score))
    result_handle.close()
    outfile.close()

def assign_functions(inFasta, eggNOGdirectory,outdir):
    functionDict = {}
    proteinDict = {}
    membersFile = open(eggNOGdirectory + '/' + "AllNOG.members.txt", 'r')
    funcFile = open(eggNOGdirectory + '/' + "AllNOG.funccat.txt", 'r')
    descriptionFile = open(eggNOGdirectory + '/' + "AllNOG.description.txt", 'r')
    for i, line in enumerate(membersFile):
        if i > 0:
            line = line.strip()
            nog, protein, start, end = line.split('\t')
            protein = protein.split('.')[1]
            proteinDict[protein] = nog
    for line in funcFile:
        line = line.strip()
        nog, function = line.split('\t')
        functionDict[nog] = {}
        functionDict[nog]["function"] = function
    for line in descriptionFile:
        line = line.strip()
        if len(line.split('\t')) > 1:
            nog, description = line.split('\t')
            functionDict[nog]["description"] = description
        else:
            functionDict[line]["description"] = "no description"
    membersFile.close()
    funcFile.close()
    descriptionFile.close()

    blastFile = open( outdir + '/' + os.path.split(inFasta)[1] + "_blast.txt", 'r')
    outfile = open( outdir + '/' + os.path.split(inFasta)[1] + "_functions.txt", 'w')
    for i, line in enumerate(blastFile):
        if i > 0:
            line = line.strip()
            query, match, e, start, end,hsp_score = line.split('\t')
            try:
                nog = proteinDict[match.split('.')[1]]
                outfile.write("%s\t%s\t%s\t%s\t%s\t%s\n" % (query, match, nog, hsp_score, \
                functionDict[nog]["function"], 
                functionDict[nog]["description"]))
            except KeyError:
                print match.split('.')[1]  + " not in dictionary"
    blastFile.close()
    outfile.close()

if None in get_arguments(sys.argv[1:]):
    usage()
    sys.exit(2)
else:
    inFasta, eggNOGdirectory,blastcmd,outdir = get_arguments(sys.argv[1:])

blast_eggNOG(inFasta, eggNOGdirectory,blastcmd,outdir)
parse_blast(inFasta, eggNOGdirectory,outdir)
assign_functions(inFasta, eggNOGdirectory,outdir)
