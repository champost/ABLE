#!/usr/bin/env python
# Usage 
import sys
from pyfaidx import Fasta

"""
Reads in a fasta file.
"""
def get_raw_fasta(fil):
	return Fasta(fil, sequence_always_upper=True)

"""
Takes a block from a set of individual fasta files given in list l and turns it into a tuple
"""
def block_to_tuple(l, chrom, coord, blocklen):
	a=[]
	for i in range(len(l)):
		a.append(l[i][chrom][coord-blocklen:coord].seq)
	return zip(*a)

"""
Filters out positions containing N (or any other character not in code_set) in any individual in the alignment
"""
def n_filter(t):
	clean_list = []
	code_set = {'A','T','G','C','M','R','W','S','Y','K'}
	for j in t:
		if code_set.issuperset(j):
			 clean_list.append(j)
	return clean_list

"""
Removes monomorphic sites
"""
def remove_mono(t):
	hom_set = {'A','T','G','C'}
	newlist = []
	for j in t:
		if j.count(j[0]) == len(j) and j[0] in hom_set:
			pass
		else:
			newlist.append(j)
	return newlist

"""
Combines the functions above to cut and filter out blocks that contain more than minlen Ns. The output is a dic of block positions and polymorphism info.
"""

def blockcutter_sub(l, chrom, blocklen, minlen, sublen):
	ll={}
	for i in range(blocklen,len(l[0][chrom]),blocklen):
#	for i in range(blocklen,500000,blocklen):
		blockprefltrd = n_filter(block_to_tuple(l,chrom,i,blocklen))
		if len(blockprefltrd) > minlen:
			blockprefltrd_trimmed = blockprefltrd[0:minlen]
			for j in range(sublen,minlen+1,sublen):
				blockprefltrd_trimmed_sub = blockprefltrd_trimmed[j-sublen:j]
				ll[(chrom,i+j-minlen-sublen,i-minlen+j)]=zip(*remove_mono(blockprefltrd_trimmed_sub))
	return ll

"""
Loops above over all chromosomes. The output is still a single dic 
"""
def all_blockcutter_sub(l, chromlist, blocklen, minlen, sub_blocks=1):
	sublen = int(minlen/sub_blocks)
	ll = {}
	for chrom in chromlist:
		ll.update(blockcutter_sub(l, chrom, blocklen, minlen, sublen))		
	return ll

"""
Turns position tuple into a string which can be used as a block label
"""
def blocklabel(l):
	return l[0]+'_'+str(l[1])+'-'+str(l[2])

"""
Filters chromkeys to exclude a list of contigs (e.g. mt or sex chromosomes)
"""
def chrom_filter(chromkeys, flist):
	chromkeysfltrd = []
	for c in chromkeys:
		if c not in flist:
			chromkeysfltrd.append(c)
		else:
			pass
	return chromkeysfltrd

"""
Phases a sequence s containing ambiguity codes for heterozygous sites. The phasing is NOT random but in the order given in phasedict. 
"""
def phase_seqs(s):
	phasedict = {'A':('A','A'),'T':('T','T'),'G':('G','G'),'C':('C','C'),'M':('A','C'),'R':('A','G'),'W':('A','T'),'S':('C','G'),'Y':('C','T'),'K':('G','T')}
	ll = []
	for i in list(s):
		ll.append(phasedict.get(i))
	return zip(*ll)#

def phase_seq_block(s):
	return map(phase_seq,s)
	return zip(*ll)

def main():
	genes = map(get_raw_fasta,(sys.argv)[4:])
	blocklen = int(sys.argv[1])
	# the minimum physical length of each large block 
	minlen = int(sys.argv[2])
	# the # of subblocks:
	sub_blocks = int(sys.argv[3])
	chromkeys = list(genes[0].keys())
	notautosome = ['2a','2b_random','6','MT','X', 'X_random','6_cox_hap1_random']
#	chromkeysfltrd = [chromkeys[8]]
	chromkeysfltrd = chrom_filter(chromkeys, notautosome)
#	print chromkeysfltrd
	outdict = all_blockcutter_sub(genes, chromkeysfltrd, blocklen, minlen, sub_blocks)
	for key, value in outdict.iteritems():
		print 
		print '//'
		print blocklabel(key)
		for i in value:
			#~ for u in phase_seqs(i):
				#~ print ''.join(u)
				#~ print u
			print ''.join(i)
if __name__== "__main__":
	main()
