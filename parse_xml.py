#!/usr/bin/env python
import argparse
import Bio
import re
import os
import importlib
from Bio import Entrez
from Bio import SeqIO

Entrez.email = "molly.gibson@wustl.edu"

import sys



parser = argparse.ArgumentParser(description='This program takes the xml output from pubmed search and outputs the title, abstract, and keywords')
#Input
parser.add_argument('--in', type=str, dest='input', help='Input file base name.')
parser.add_argument('--out', type=str, dest='output',help='Output file.')
parser.add_argument('--pmid_out', type=str, dest='pubs', help='Output of pmids.')
args = parser.parse_args()

# get all files with the base string
output = open(args.output, 'w')
pub_out = open(args.pubs, 'w')
pub_out.write("PMID\tJournal\tAuthor\tTitle\n")
for file_name in os.listdir('.'):
    if file_name.startswith(args.input):
        
        handle = open(file_name)
        records = Entrez.read(handle)


        for record in records['PubmedArticle']:
            try:
                #print(record['MedlineCitation']['Article']['AuthorList'][0]['LastName'])
                output.write("PMID:" + record['MedlineCitation']['PMID'] + "\n")
                pub_out.write(record['MedlineCitation']['PMID'] + "\t" + record['MedlineCitation']['Article']['Journal']['JournalIssue']['PubDate']['Year'] + "\t" + record['MedlineCitation']['Article']['Journal']['Title'] + "\t" + record['MedlineCitation']['Article']['AuthorList'][0]['LastName'] + "\t" + record['MedlineCitation']['Article']['ArticleTitle']  + "\n")
                output.write(record['MedlineCitation']['Article']['ArticleTitle'] + "\n")
                if 'Abstract' in record['MedlineCitation']['Article']:
                    for item in record['MedlineCitation']['Article']['Abstract']['AbstractText']:
                        output.write(item + "\n")
                if 'KeywordList' in record['MedlineCitation']:
                    for item in record['MedlineCitation']['KeywordList']:
                        for item2 in item:
                            output.write(item2 + "\n")
                output.write("RECORD\n")
            except:
                if 'BookDocument' in record:
                    print("Book Document")
                else:
                    print(record)
