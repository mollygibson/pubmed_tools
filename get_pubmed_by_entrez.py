#!/usr/bin/env python
import argparse
import Bio
import re
import math
from Bio import Entrez
from Bio import SeqIO
Entrez.email = "molly.gibson@wustl.edu"


parser = argparse.ArgumentParser(description='This program takes as input a list of NCT numbers and finds associated pubmed ids and returns them')
#Input
parser.add_argument('--search', type=str, dest='strTerm', help='Pubmed Search Terms')
parser.add_argument('--search_list', type=str, dest='strList', help='List of terms to search')
parser.add_argument('--nct', type=str, dest='strNCT',help='List of NCT numbers, one on each line.')
parser.add_argument('--out', type=str, dest='output',help='Output file.')
args = parser.parse_args()


if args.strNCT:
    nct_numbers = []
    for line in open(args.strNCT, 'r'):
        nct_numbers.append(line.rstrip())

    id_numbers = []
    for nct_number in nct_numbers:
        
        handle = Entrez.esearch(db="pubmed", term=nct_number, rettype="docsum")
        record = Entrez.read(handle)
        handle.close()
        
        if len(record["IdList"])>0:
            for idnum in record["IdList"]:
                id_numbers.append(idnum)

    handle = Entrez.efetch(db="pubmed", id=",".join(id_numbers), rettype="medline", retmode="text")

elif args.strList:
    term_list = []
    for line in open(args.strList, 'r'):
        term_list.append("\"" + line.rstrip() + "\"")    

    search_terms = " OR ".join(term_list)

    search_results = Entrez.esearch(db="pubmed", term=search_terms, retmax=200000)
    search_ids =  Entrez.read(search_results)["IdList"]

    num_records = len(search_ids)
    print("We found " + str(num_records) + " records matching your search")

    rounded = int(math.ceil(num_records / 10000.0)) * 10000

    for x in range(1,rounded,10000):
        handle = Entrez.efetch(db="pubmed", id=search_ids, rettype="abstract", retmode="xml", retstart=x)

        out_file = open(args.output + "_" + str(x) + ".txt", 'w')
        out_file.write(handle.read())
        out_file.close()

elif args.strTerm:

    search_terms = args.strTerm
    search_results = Entrez.esearch(db="pubmed", term=search_terms, retmax=200000)
    search_ids =  Entrez.read(search_results)["IdList"]

    num_records = len(search_ids)
    print("We found " + str(num_records) + " records matching your search")

    rounded = int(math.ceil(num_records / 10000.0)) * 10000

    for x in range(1,rounded,10000):
        handle = Entrez.efetch(db="pubmed", id=search_ids, rettype="abstract", retmode="xml", retstart=x)
        
        out_file = open(args.output + "_" + str(x) + ".txt", 'w')
        out_file.write(handle.read())
        out_file.close()

            
