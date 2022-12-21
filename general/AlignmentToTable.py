######### AlignmentToTable.py ##########
# Takes an alignment file (fasta format) and 
# returns universal position
# and ungapped positions as a table
# arg 1 is fasta, with gaps. arg 2 is table
# TJE 2022 11 23
########################################
from Bio import SeqIO
import sys
from collections import defaultdict as dd

all_seqs = {}
with open(sys.argv[1], 'r') as f1:
	for record in SeqIO.parse(f1, 'fasta'):
		all_seqs[record.id] = record.seq
		seq_length = len(record.seq)


pos_table = dd(lambda: dd(int))
ungapped_dict = {seqname:0 for seqname in all_seqs.keys()}
for pos in range(seq_length):
	for seqname, seq in all_seqs.items():
		if seq[pos] != '-': 
			ungapped_dict[seqname] += 1
			pos_table[pos][seqname] = ungapped_dict[seqname]
		else: pos_table[pos][seqname] = 'NA'
with open(sys.argv[2], 'w+') as f2:
	header = 'pos' + "\t" + "\t".join([i for i in ungapped_dict.keys()]) + '\n'
	f2.write(header)
	for pos, ungapped_pos_dict in pos_table.items():
		newline = str(pos) + '\t' + "\t".join([str(i) for i in ungapped_pos_dict.values()]) + '\n'
		f2.write(newline)
