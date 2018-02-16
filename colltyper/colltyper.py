#!/usr/bin/env python

''' Script for typing MTB using the Coll scheme. Takes in a VCF file and votes on lineage'''

# Take in FASTA or VCF

import sys
import argparse
import vcf
import csv
import os
from operator import itemgetter
from pkg_resources import resource_filename
from __init__ import __version__

COLL_VERSION = __version__

def CollArgumentParser():
    parser = argparse.ArgumentParser(
        description = 'Coll scheme typing of TB, version %s' % COLL_VERSION)
    parser.add_argument('vcf',
        help='VCF file to type. Reference must be H37Rv', default=None)
    parser.add_argument('--fasta',
        help='File is not VCF but FASTA. (H37Rv with variants mapped in.)',
        action='store_true',
        default=False)
    parser.add_argument('-v', '--version',
        help='Show version and exit.',
        action='version',
        version=COLL_VERSION)
    parser.add_argument('--scheme',
        help='Scheme classification dictionary. (If using non-default)',
        default=os.path.join(resource_filename(__name__,'data/Coll_scheme_classification.csv')))
    args = parser.parse_args()

    return args

def ReadClassification(scheme):
    try:
        myscheme = csv.reader(open(scheme,"rU"), skipinitialspace=False,
            delimiter=',')
    except Exception as e:
        sys.stdout.write(e)
        sys.stdout.write("Failed to open classification scheme\n")
        sys.exit(2)

    schemedic = {}
    header = myscheme.next()
    lineages = set()

    try:
        lineagecol = header.index("lineage")
        poscol = header.index("Position")
        allelecol = header.index("Allele change")
    except ValueError:
        sys.stdout.write("ERROR: The following columns must be in the scheme file: 'lineage', 'Position', 'Allele change'\n")
        sys.exit(3)

    for row in myscheme:
        # NOTE: 4.9 is reference, so REF/ALT is switched for that lineage and lineage 4
        if row[lineagecol].lstrip("lineage") == "4.9" or row[lineagecol].lstrip("lineage") == "4":
            schemedic[row[poscol]] = {"Lineage": row[lineagecol].lstrip("lineage"), "Allele": row[allelecol][0]}
        else:
            schemedic[row[poscol]] = {"Lineage": row[lineagecol].lstrip("lineage"), "Allele": row[allelecol][-1]}
        lineages.add(row[lineagecol].lstrip("lineage"))

    return schemedic, lineages
    
def Classify(vcf_reader, schemedic, lineages):
    lineagevote = {l: {"Value": "False"} for l in lineages}
    for record in vcf_reader:
        if record.INFO["TYPE"] != ["snp"]:
            continue
        pos = str(record.POS)
        if pos not in schemedic:
            continue

        ref = record.REF
        if len(record.ALT) > 1:
            sys.stdout.write("WARNING: Multiple variants at position %s, check results carefully\n" % pos)
        allele = record.ALT[0]
        sample = record.samples[0]
        samplename = sample.sample

        # Check if mutation votes for any particular lineage

        if (allele == schemedic[pos]["Allele"]):
            lineagevote[schemedic[pos]["Lineage"]]["Value"] = "True"
            lineagevote[schemedic[pos]["Lineage"]]["Read depth"] = sample["DP"]
            lineagevote[schemedic[pos]["Lineage"]]["Genotype likelihood"] = sample["GL"]

    return lineagevote

def sortresults(vote):
    # Only report true votes
    truevotes = {k:v for k,v in vote.items() if v["Value"] == "True"}
    #resultsdic = {}

    # Vote should be pruned if parent mutation is missing, e.g. 4.1.1 should not show up if 4.1 is missing
    # NOT YET IMPLEMENTED
    # Sort by genotype likelihood
    GL_list = []
    for res in truevotes:
        GL = truevotes[res]["Genotype likelihood"]
        GLpseudo = (GL[0] - 1.0) / (GL[1] - 1.0)
        GL_list.append( (res, truevotes[res]["Read depth"], GLpseudo) )
    sorted_GL_list = sorted(GL_list, key=itemgetter(2), reverse=True)
    return sorted_GL_list


def main():
    args = CollArgumentParser()

    if args.fasta:
        sys.stdout.write("FASTA not supported in this version\n")
        sys.exit(0)
    else:
        try:
            vcf_reader = vcf.Reader(open(args.vcf,'rU'))
        except:
            sys.stdout.write("ERROR: Failed to read VCF file. Exiting.\n")
            sys.exit(1)

    # Read Classification file
    try:
        schemedic, lineages = ReadClassification(args.scheme)
    except:
        sys.stdout.write(args.scheme)
        sys.stdout.write("ERROR: Failed to understand classification scheme\n")
        sys.exit(4)

    vote = Classify(vcf_reader, schemedic, lineages)
    sortedvote = sortresults(vote)

    print("Lineage\tRead depth\tRel. gen. lik.")
    for r in sortedvote:
        resultsstring = "\t".join([str(i) for i in r])
        print(resultsstring)

    sys.stdout.write("Thank you for using Colltyper\n")
    sys.exit(0)


if __name__ == '__main__':
    main()