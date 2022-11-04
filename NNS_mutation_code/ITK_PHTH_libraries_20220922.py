#BTK_regional_mutagenesis.py
#saraste dimer
#TJE 2021 09 22
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from collections import defaultdict as dd
import sys
import random 
import json
from types import SimpleNamespace
from Bio.Restriction import *

MAX_SEQ_SUM = 6000 #lib size
NUM_ALT_WT = 0

class mutagenesis_region:
		def __init__(self, wt_dna_seq, region_id, left_handle, right_handle, internal, begin_num):
			self.wt_dna_seq = SeqRecord(Seq(wt_dna_seq), id = region_id)
			self.wt_prot = self.wt_dna_seq.translate()
			self.left_handle = left_handle
			self.internal = internal
			self.right_handle = right_handle
			self.begin_num = begin_num

all_codon_dict = dd(list)
with open("hsap_codon_freq_T.txt", 'r') as f:
	next(f)
	for line in f:
		all_codon_dict[line.split()[1]].append(line.split()[0])

def populate_seq(region, codon_list, all_codon_dict, NUM_ALT_WT = NUM_ALT_WT):
	aa_SEQ, names_SEQ, recode_SEQ, recoded_names = [], [], [], []

	#add in the WT sequence:
	aa_SEQ.append(str(region.left_handle + region.wt_dna_seq.seq + region.right_handle))
	names_SEQ.append(region.wt_dna_seq.id)

	for idx in range(len(region.wt_prot)):
		wt_codon = region.wt_dna_seq.seq[(idx*3):(idx*3+3)]
		for alt_codon in codon_list:
			if wt_codon == alt_codon:
				continue
			else:
				new_variant = Seq(str(region.wt_dna_seq.seq[:(idx * 3)]) + alt_codon + str(region.wt_dna_seq.seq[((idx + 1) * 3):]))
				variant_name = str(wt_codon + str(idx + region.begin_num) + alt_codon + '_' + wt_codon.translate() + str(idx + region.begin_num) + Seq(alt_codon).translate())
				if len(BamHI.search(new_variant)) > 0: #note it's possible that this scheme creates e.g. XmaI sites when removing BsaI sites, but that's probably unlikely
					while len(BamHI.search(new_variant)) > 0: 
						new_variant = Seq(recode(new_variant, new_variant.translate(), all_codon_dict, num_var = 1))
					variant_name = variant_name + "_BamHI_Recode"
				if len(XmaI.search(new_variant)) > 0: 
					while len(XmaI.search(new_variant)) > 0: 
						new_variant = Seq(recode(new_variant, new_variant.translate(), all_codon_dict, num_var = 1))
					variant_name = variant_name + "_XmaI_Recode"
				if len(BsaI.search(new_variant)) > 0: 
					while len(BsaI.search(new_variant)) > 0: 
						new_variant = Seq(recode(new_variant, new_variant.translate(), all_codon_dict, num_var = 1))
					variant_name = variant_name + "_BsaI_Recode"
				aa_SEQ.append(region.left_handle + str(new_variant) + region.right_handle)
				names_SEQ.append(variant_name)
	if len(aa_SEQ) < len(region.wt_prot) * 20 * 2: additional_seqs = len(region.wt_prot) * 20 * 2 - len(aa_SEQ) + NUM_ALT_WT#Fill in any lost seqs with WT
	else: additional_seqs = NUM_ALT_WT
	idx = 0
	while(idx < additional_seqs):
		recoded_variant = recode(region.wt_dna_seq.seq, region.wt_prot, all_codon_dict, num_var = 5)
		if recoded_variant in recode_SEQ: continue #make sure I don't pick the same ones
		elif len(BamHI.search(Seq(recoded_variant))) > 0: continue
		elif len(XmaI.search(Seq(recoded_variant))) > 0: continue
		elif len(BsaI.search(Seq(recoded_variant))) > 0: continue
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
def file_writer(file, aa_SEQ, names_SEQ, MAX_SEQ_SUM):
	SeqCounter = 0
	with open(file, 'w+') as f:
		for idx, dna_seq in enumerate(aa_SEQ):
			if SeqCounter == MAX_SEQ_SUM: break
			f.write(">TJE_" + names_SEQ[idx] + "\n")
			f.write(dna_seq + "\n")
			SeqCounter += 1
	print("Sequences written: ", SeqCounter)

def recode(dna_seq, prot_seq, all_codon_dict, num_var, all = False):
	randposaa = random.sample(range(len(prot_seq)), k = num_var)
	if all: randposaa = [i for i in range(len(prot_seq))]
	newseq = dna_seq
	for idx in randposaa:
		wt_codon = dna_seq[(idx*3):(idx*3+3)]
		choice_list = [i for i in all_codon_dict[prot_seq[idx]] if i != wt_codon]
		if len(choice_list) < 1: alt_codon = wt_codon
		else: alt_codon = random.choice(choice_list)
		newseq = str(newseq[:(idx * 3)]) + alt_codon + str(newseq[((idx + 1) * 3):])
	return(str(newseq))



codon_list = []
with open("top_2_codons_hsap_by_usage.txt", 'r') as f:
	next(f)
	for line in f:
		codon_list.append(line.split()[0])

def create_seqs(regionA, regionB, regionC, filename):
	names_SEQ_A, aa_SEQ_A = populate_seq(regionA, codon_list, all_codon_dict)
	names_SEQ_B, aa_SEQ_B = populate_seq(regionB, codon_list, all_codon_dict)
	names_SEQ_C, aa_SEQ_C = populate_seq(regionC, codon_list, all_codon_dict)
	
	aa_SEQ = aa_SEQ_A + aa_SEQ_B + aa_SEQ_C
	names_SEQ = names_SEQ_A + names_SEQ_B + names_SEQ_C

	file_writer(filename, aa_SEQ, names_SEQ, MAX_SEQ_SUM)

#parse the json file.
def read_json(json_file):
	with open(json_file, 'r') as f:
		all_regions = json.load(f)
		all_regions_class = {}
		for region_name, region in all_regions.items():
			all_regions_class[region_name] = mutagenesis_region(wt_dna_seq = region['wt_dna_seq'],\
				region_id = region['region_id'],\
				internal = region['internal'],\
				left_handle  = region['left_handle'],\
				right_handle  = region['right_handle'],\
				begin_num = region['begin_num'])
	return(all_regions_class)

all_regions_class = read_json('ITK_PHTH.json')
create_seqs(all_regions_class['region_15'], all_regions_class['region_16'], all_regions_class['region_17'], 'ITK_PHTH_Seqs_Draft1_20220922.txt')


###Subu primers:
# 96	ATTCAAGGGT TGGACGACTC
# 154	TATCACGGAA GGACTCAACG
# 62	AGGACACCAG ACCAATGAAG
# 105	CCGGGAGGAA GATATAGCAC