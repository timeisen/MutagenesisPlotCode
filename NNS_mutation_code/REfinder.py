from collections import defaultdict as dd
import sys
from Bio import SeqIO
from Bio.Seq import Seq
import regex
import math
from Bio.Restriction import *
import argparse
import re




def main():
	parser = argparse.ArgumentParser(description='find seqs that contain a restriction enzyme site')
	parser.add_argument('-i', '--input', help = 'Input file, .json format')
	parser.add_argument('-e', '--enzyme', help = 'Restriction enzyme to recode. Must be standard naming.')
	args = parser.parse_args()

	enzymeF = eval(args.enzyme)

	with open(args.input, 'r') as f1: #main code block that parses the fastq file (.extendedFrags )
		for record in SeqIO.parse(f1, "fasta"):
			if len(enzymeF.search(record.seq)) > 0: 
				pos_list = [m.start() for m in re.finditer(enzymeF.site, str(record.seq))]
				print(record.id, 'Positions:', pos_list)

if __name__ == "__main__":
    main()