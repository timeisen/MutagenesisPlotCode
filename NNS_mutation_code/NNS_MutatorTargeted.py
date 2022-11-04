##############################
# NNS_MutatorTargeted.py
# Timothy J. Eisen
# Written 2021 09 22
# Updated 2022 10 26
# Updated 2022 11 02, adding functionality for specific NNS positions in the .json list
# Updated 2022 11 03, adding functionality for argsparse
##############################

#imports
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from collections import defaultdict as dd
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

def populate_seq(region, codon_list):
	aa_SEQ, names_SEQ, = [], []
	#add in the WT sequence:
	aa_SEQ.append(str(region.left_handle + region.wt_dna_seq.seq + region.right_handle))
	names_SEQ.append(region.wt_dna_seq.id)

	#this is the main code block that performs the mutation
	for prot_idx in region.positions:
		idx = prot_idx - region.begin_num
		wt_codon = region.wt_dna_seq.seq[(idx*3):(idx*3+3)] #get the wt codon
		for alt_codon in codon_list: #for each codon in the codon list, generate an "alt codon"
			if wt_codon == alt_codon: #if the wt codon is the alt codon, continue through the loop
				continue
			else:
				new_variant = Seq(str(region.wt_dna_seq.seq[:(idx * 3)]) + alt_codon + str(region.wt_dna_seq.seq[((idx + 1) * 3):])) #generate the new variant
				variant_name = str(wt_codon + str(prot_idx) + alt_codon + '_' + wt_codon.translate() + str(prot_idx) + Seq(alt_codon).translate()) #generate the name
				aa_SEQ.append(region.left_handle + str(new_variant) + region.right_handle) #append the variant to a list, with handles
				names_SEQ.append(variant_name) #append the name to a list

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

def create_seqs(regionA, filename, codon_file):
	codon_list, codon_dict = generate_codon_dicts(codon_file)
	names_SEQ, aa_SEQ = populate_seq(regionA, codon_list)	
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
	parser = argparse.ArgumentParser(description = 'Generates codon variants of a sequence baesd on user-supplied positions. \n To mutate all position, supply an empty bracket [] in the .json file.') #argument parser
	parser.add_argument('-o', '--output', help = 'Output file, fasta')
	parser.add_argument('-i', '--input', help = 'Input file, .json format')
	parser.add_argument('-s', '--sequence_name', help = 'Name of sequence to query. Must correspond to json entry.')
	parser.add_argument('-c', '--codons', help = 'tab delimited [codon \\t aa] codon file, with header')

	args = parser.parse_args()

	#main run segment
	all_regions_class = read_json(args.input)
	create_seqs(all_regions_class[args.sequence_name], args.output, args.codons)


if __name__ == "__main__":
    main()
