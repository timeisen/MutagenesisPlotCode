from collections import defaultdict as dd
import sys
from Bio import SeqIO
from Bio.Seq import Seq
import regex
import math
from Bio.Restriction import *
import argparse



def recode(record, seq, enzymeF, all_codon_dict):
########## THIS FUNCTION IS NOT FINISHED #########
	newseq = dna_seq
	baseseq = seq.removesuffix(record.right_handle).removeprefix(record.left_handle) #just the variant
	pos = baseseq.find(enzymeF.site)
	prot_seq_var = baseseq.translate()
	codonnumber = math.floor(pos / 3)


	return(str(newseq))

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

def main():
	parser = argparse.ArgumentParser(description='Generates recoded wildtype sequences.')
	parser.add_argument('-o', '--output', help = 'Output file, fasta')
	parser.add_argument('-i', '--input', help = 'Input file, .json format')
	parser.add_argument('-s', '--sequence_name', help = 'Name of sequence to query. Must correspond to json entry.')
	parser.add_argument('-e', '--enzyme', help = 'Restriction enzyme to recode')
	args = parser.parse_args()

	enzymeF = eval(args.enzyme)
	region = read_json(args.input)[args.sequence_name]

	with open(args.input, 'r') as f1: #main code block that parses the fastq file (.extendedFrags )
		for record in SeqIO.parse(f1, "fasta"):
			if len(enzymeF.search(record.seq)) > 0: 
				print('recoding:', record.id)
				newseq = recode(record, seq, enzymeF, all_codon_dict)

if __name__ == "__main__":
    main()