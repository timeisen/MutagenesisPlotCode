##############################
# GenerateRecodedWT.py
# Timothy J. Eisen
# Written 2021 09 22
# Updated 2022 10 26
# Updated 2022 11 02, adding functionality for specific NNS positions in the .json list
# Generates Recoded wildtype seqs
##############################

#imports
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from collections import defaultdict as dd
import random 
import sys
import json
import argparse

class mutagenesis_region: #class of regions
		def __init__(self, wt_dna_seq, region_id, left_handle, right_handle, positions, begin_num):
			self.wt_dna_seq = SeqRecord(Seq(wt_dna_seq), id = region_id)
			self.wt_prot = self.wt_dna_seq.translate()
			self.left_handle = left_handle
			if positions:
				self.positions = [int(i) for i in positions] #convert to int
			else:
				self.positions = [int(i) for i in range(self.wt_dna_seq.translate())] #convert to int
			self.right_handle = right_handle
			self.begin_num = begin_num

def generate_codon_dicts(file):
	codon_list = [] #list of all codons, must be a file called "Codons.txt"
	codon_dict = dd(list)
	with open(file, 'r') as f:
		next(f) #skip the header
		for line in f:
			codon_list.append(line.split()[0].strip())
			codon_dict[line.split()[1].strip()].append(line.split()[0].strip())
	return(codon_list, codon_dict)

def recode(dna_seq, prot_seq, all_codon_dict, num_var):
	randposaa = random.sample(range(len(prot_seq)), k = num_var)
	newseq = dna_seq
	for idx in randposaa:
		wt_codon = dna_seq[(idx*3):(idx*3+3)]
		choice_list = [i for i in all_codon_dict[prot_seq[idx]] if i != wt_codon]
		if len(choice_list) < 1: alt_codon = wt_codon
		else: alt_codon = random.choice(choice_list)
		newseq = str(newseq[:(idx * 3)]) + alt_codon + str(newseq[((idx + 1) * 3):])
	return(str(newseq))

def populate_seq(region, codon_list, all_codon_dict, NUM_ALT_WT):
	aa_SEQ, names_SEQ, recode_SEQ, recoded_names = [], [], [], []
	idx = 0

	while(idx < NUM_ALT_WT):
		recoded_variant = recode(region.wt_dna_seq.seq, region.wt_prot, all_codon_dict, num_var = 5)
		if recoded_variant in recode_SEQ: continue #make sure I don't pick the same ones
		else:
			recode_SEQ.append(region.left_handle + recoded_variant + region.right_handle)
			recoded_names.append(region.wt_dna_seq.id + '_' + str(idx))
			idx += 1
	[aa_SEQ.append(i) for i in recode_SEQ]
	[names_SEQ.append(i) for i in recoded_names]

	#change asteriks to X for writing
	names_SEQ = [val.replace("*", "X") for val in names_SEQ]
	return(names_SEQ, aa_SEQ)

#write file.
def file_writer(file, aa_SEQ, names_SEQ):
	SeqCounter = 0

	#change * to X for writing, biopython defaults to asteriks. 
	names_SEQ = [val.replace("*", "X") for val in names_SEQ]

	with open(file, 'w+') as f:
		for idx, dna_seq in enumerate(aa_SEQ):
			f.write(">" + names_SEQ[idx] + "\n")
			f.write(dna_seq + "\n")
			SeqCounter += 1
	print("Sequences written: ", SeqCounter)

def create_seqs(regionA, filename, codon_file, NUM_ALT_WT):
	codon_list, codon_dict = generate_codon_dicts(codon_file)
	names_SEQ, aa_SEQ = populate_seq(regionA, codon_list, codon_dict, NUM_ALT_WT)	
	file_writer(filename, aa_SEQ, names_SEQ)

#parse the json file.
def read_json(json_file):
	with open(json_file, 'r') as f:
		all_regions = json.load(f)
		all_regions_class = {}
		for region_name, region in all_regions.items():
			all_regions_class[region_name] = mutagenesis_region(wt_dna_seq = region['wt_dna_seq'],\
				region_id = region['region_id'],\
				left_handle  = region['left_handle'],\
				right_handle  = region['right_handle'],\
				positions  = region['positions'],\
				begin_num = region['begin_num'])
	return(all_regions_class)


def main():
	parser = argparse.ArgumentParser(description='Generates recoded wildtype sequences.')
	parser.add_argument('-o', '--output', help = 'Output file, fasta')
	parser.add_argument('-i', '--input', help = 'Input file, .json format')
	parser.add_argument('-s', '--sequence_name', help = 'Name of sequence to query. Must correspond to json entry.')
	parser.add_argument('-n', '--number', type = int, help = 'Number of additional sequences to generate')
	parser.add_argument('-c', '--codons', help = 'tab delimited [codon \\t aa] codon file, with header')

	args = parser.parse_args()

	all_regions_class = read_json(args.input)
	create_seqs(all_regions_class[args.sequence_name], args.output, args.codons, args.number)


if __name__ == "__main__":
    main()

