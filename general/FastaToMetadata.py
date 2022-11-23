##TJE 2022 11 21
#for generating tab-delimited metdata 
from Bio import SeqIO
import sys

with open(sys.argv[1], 'r') as f1, open(sys.argv[2], 'w+') as f2:
	write_line = "\t".join(['seq_name',	'wt_codon',	'pos',	'var_codon',	'wt_aa',	'var_aa',	'mut_short',	'variant recode', 'region'])
	f2.write(write_line + "\n")
	for record in SeqIO.parse(f1, 'fasta'):
		seq_name = record.id
		if 'TJE_ITK_PHTH' in seq_name:
			wt_codon = 'NA'
			pos = 'NA'
			var_codon = 'NA'
			wt_aa = 'NA'
			var_aa = 'NA'
			mut_short = 'NA'
			variant_recode = 'NA'
			region = seq_name.split("_")[4].strip()
		else:
			wt_codon = seq_name.split('_')[1][:3]
			pos = seq_name.split('_')[1][3:-3]
			var_codon = seq_name.split('_')[1][-3:]
			wt_aa = seq_name.split('_')[2][:1]
			var_aa = seq_name.split('_')[2][-1:]
			mut_short = wt_aa + pos + var_aa
			if int(pos) < 51: region = 'R15'
			elif int(pos) > 50 and int(pos) < 101: region = 'R16'
			else: region = 'R17'
			if 'recode' in seq_name: variant_recode = 'TRUE'
			else: variant_recode = 'FALSE'
		write_line = "\t".join([seq_name,	wt_codon,	pos,	var_codon,	wt_aa,	var_aa,	mut_short,	variant_recode, region])
		f2.write(write_line + "\n")
