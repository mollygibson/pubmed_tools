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
counts = collections.defaultdict(int)
publications = s.split("RECORD\n")

# For each publication, find which diseases are found 
for pub in publications:
    if pub:
        pub_out = regex.sub(' ', pub)

        if not re.search(r'neuron', pub_out.lower()):
            phrase_list =  []
            
            
            for disease in diseases:            
                if re.search(r'\b' + disease + r'\b', pub_out.lower()):
                    if not disease in phrase_list:
                        phrase_list.append(disease)
                
            for phrase in phrase_list:
                counts[phrase] += 1


output = open(args.output, 'w')

for item in sorted(counts, key=counts.get, reverse=True):
    print_phrase = False
    for disease in diseases:
        if item == disease:
            print_phrase = True
    if print_phrase:
        output.write(item + "\t"+ str(counts[item]) + "\n")
