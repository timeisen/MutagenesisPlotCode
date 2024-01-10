from Bio import SeqIO
import sys
from collections import defaultdict as dd

AA_letters = 'ACDEFGHIKLMNPQRSTVWY-X'
aa_dict = dd(lambda: {aa: 0 for aa in AA_letters})

with open(sys.argv[1], 'r') as f:
	for record in SeqIO.parse(f, 'fasta'):
		for pos, aa in enumerate(record.seq):
			aa_dict[pos][aa] += 1

with open(sys.argv[2], 'w+') as f:
	header = "pos\t" + "\t".join(AA_letters) + '\n'
	f.write(header)
	for pos, aa_count_dict in aa_dict.items():	
		newline = str(pos) + '\t' + "\t".join([str(aa_count_dict[aa]) for aa in AA_letters]) + '\n'
		f.write(newline)