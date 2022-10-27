##############################################
#FastaToMetadata.py
#Timothy J. Eisen
#2022 08 05
#Converts a fasta file to a metadata sheet for small libraries
##############################################

import sys
from Bio import SeqIO

wt_names = ['Qi_Contact', 'canon_site', 'periph_site', 'ITK_ActivationLoop', 'Region8']

def parser(seq_name):
	seq_name_mod = seq_name.replace("TJE_", "")
	if any([name in seq_name_mod for name in wt_names]):
		wt_codon = 'NA'
		pos = '0'
		var_codon = 'NA'
		wt_aa = 'NA'
		var_aa = 'WT'
		mut_short = 'WT'
		variant = 'FALSE'
		recode = 'NA'
	elif 'Recode' in seq_name_mod:
		wt_codon = seq_name_mod.split("_")[0][0:3]
		pos = seq_name_mod.split("_")[0][3:-3]
		var_codon = seq_name_mod.split("_")[0][-3:]
		wt_aa = seq_name_mod.split("_")[1][0]
		var_aa = seq_name_mod.split("_")[1][-1]
		mut_short = seq_name_mod.split("_")[1]
		variant = 'TRUE'
		recode = seq_name_mod.split("_")[2]	
	else:
		wt_codon = seq_name_mod.split("_")[0][0:3]
		pos = seq_name_mod.split("_")[0][3:-3]
		var_codon = seq_name_mod.split("_")[0][-3:]
		wt_aa = seq_name_mod.split("_")[1][0]
		var_aa = seq_name_mod.split("_")[1][-1]
		mut_short = seq_name_mod.split("_")[1]
		variant = 'TRUE'
		recode = 'NA'
	return([seq_name, wt_codon, pos, var_codon, wt_aa, var_aa, mut_short, variant, recode])

with open(sys.argv[1], 'r') as f1, open(sys.argv[2], 'w+') as f2:
	header = "seq_name\twt_codon\tpos\tvar_codon\twt_aa\tvar_aa\tmut_short\tvariant\trecode\n"
	f2.write(header)
	for record in SeqIO.parse(f1, 'fasta'):
		newline = "\t".join(parser(record.id)) + "\n"
		f2.write(newline)
