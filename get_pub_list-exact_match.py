import re
import collections
import sys
import string
import argparse

parser = argparse.ArgumentParser(description='This program takes the xml output from pubmed search and outputs the title, abstract, and keywords')
#Input
parser.add_argument('--in', type=str, dest='input', help='Input file of parsed xml file.')
parser.add_argument('--out', type=str, dest='output',help='Output file.')
parser.add_argument('--diseases', type=str, dest='disease_in', help='List of diseases to look for')
args = parser.parse_args()

input = open(args.input, 'r')
disease_list = open(args.disease_in, 'r')

regex = re.compile('[%s]' % re.escape(string.punctuation))

# Pull in disease list and make a list of all diseases
diseases = []
for line in disease_list:
    new_line = line.rstrip()
    diseases.append(new_line.lower())


# Read in the entire input and split it into publications by splitting on RECORD
s = input.read()
output = open(args.output, 'w')
counts = collections.defaultdict(int)
publications = s.split("RECORD\n")

# For each publication, find which diseases are found 
for pub in publications:
    if pub:
        
        m = re.match("PMID:([0-9]*)",pub)
        pmid = m.group(1)
        
        pub_out = regex.sub(' ', pub)


        if not re.search(r'neuron', pub_out.lower()) and not re.search(r'olfactory bulb', pub_out.lower()) and not re.search(r'glomeruli', pub_out.lower()) and not re.search(r'olfactory epithel', pub_out.lower()):
            for disease in diseases:            
                if re.search(r'\b' + disease + r'\b', pub_out.lower()):
                    output.write(disease + "\t" + pmid + "\n")
                



