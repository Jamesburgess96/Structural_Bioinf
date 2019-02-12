#!/usr/bin/python

"""
Structural Bioinformatics assignment 1 - Phi-psi angles

When you finish the assignment, the script should create an
output file containing the amino acid propensities to be buried
in the protein core. Please ONLY modify the code in the two
indicated blocks and do NOT import additional python packages.

To run, make 'Assignment1' your working directory and use:

> python Scripts/readDSSP.py Data/dssp_dirname
"""

# Packages
from sys import argv
import os


def read_AccUnfold():
	""" Read file with maximal possible surface accessibility
	per residue. Returns the values in a dictionary """
	# open file
	try:
		f = open(os.path.join('Data', 'AccUnfold.data'), 'r')
	except:
		print('Error: cannot open AccUnfold.data')
		exit(1)
	
	# create dictionary
	unfolded_acc = {}
	for line in f:
		splitline = line.strip().split()
		# unfolded_acc[amino_acid] = area
		unfolded_acc[splitline[0]] = float(splitline[2])
	f.close()
	print(unfolded_acc)
	return unfolded_acc
	
def read_dir(d, unfolded_acc):
	""" For each DSSP file in the input directory, extract the
	total number of each amino acid type and the number of buried
	amino acids per amino acid type
	"""
	# initialize dictionaries
	# count all amino acids
	#	all_aa_count[aa] = count
	all_aa_count = {}
	# count buried amino acids
	#	buried_aa_count = count
	buried_aa_count = {}
	
	# parse all DSSP files in the directory
	for filename in [fn for fn in os.listdir(d) if fn.endswith('.dssp')]:
		try:
			f = open(os.path.join(d, filename), 'r')
		except:
			print('Error: cannot open', os.path.join(d, filename))
			exit(1)
			
		# start reading at the first line starting with '  # '
		start_reading = False
		for line in f:
			if line.startswith('  # '):
				start_reading = True
			elif start_reading:
				line = line.rstrip()
				
				# amino acid type
				aa_type = line[13]
				
				# skip amino acids marked '!'
				if aa_type == '!':
					continue
					
				if aa_type.islower():
					#print(aa_type)
				### START CODING HERE
				# write conditional statements indicating what
				# needs to happen with 'special' residues. You
				# can use a similar construct as above for 
				# aa_type = '!'
					aa_type = "C"

				if aa_type == "X": #X is any amino acid.
					continue

				if aa_type == "Z":  #Glutamic acid (E) or Glutamine (Q)
					aa_type = "Q"   #CHANGE THIS IN FINAL DRAFT TO POSSIBLY RANDOM E OR Q??
			
				### END CODING HERE
				
				# residue number
				res_num = int(line[5:10])
				# chain ID
				chain = line[11]
				# accessible surface area
				acc = float(line[34:38])
				
				# if the amino acid type is in the dictionaries,
				# add 1 to the total count, else create a new
				# entry in both dictionaries
				if aa_type in all_aa_count:
					all_aa_count[aa_type] += 1
				else:
					all_aa_count[aa_type] = 1
					buried_aa_count[aa_type] = 0
					
				# check if the amino acid is buried
				buried = decide_if_buried(aa_type, acc, unfolded_acc)
				if buried:
					buried_aa_count[aa_type] += 1
				
		# close file
		f.close()
	# return dictionaries
	return all_aa_count, buried_aa_count


def decide_if_buried(aa, acc, unfolded_acc):
	""" Calculate the fraction of buried surface area. If this
	fraction is less than seven percent, the amino acid is 
	considered buried """
	#Goes through each residue in the file, need to determine whether it is buried or not and return true or false for buried.
	buried = False
	### START CODING HERE
	# normal amino acids
	max_area = (unfolded_acc[aa])
	area = (acc/max_area)
	if area < 0.07: #7% is defined as buried.
		buried = True

	### END CODING HERE
	return buried

def print_propensities(all_aa_count, buried_aa_count, outfile):
	""" For each amino acid, calculate the propensity to be
	buried and write it into an output file """
	f = open(outfile, 'w')
	list_aa = sorted(all_aa_count.keys())
	### START CODING HERE
	# you should calculate the propemsity for each amino acid type to be buried
	# you can use all_aa_count and buried_aa_count[aa]  
	# you can use the following loop structure over all amino acids:
	# for aa in list_aa:   
	#     ... all_aa_count[aa] ...
	#     ... buried_aa_count[aa] ...
	# 
	# to print to the output file you can use:
	# print(aa, propensity_buried, file=f) 
	print(buried_aa_count)


	sum_total_residues = 0
	sum_buried_residues = 0

	for aa in list_aa:
		#print(all_aa_count[aa])
		sum_total_residues += all_aa_count[aa] #N(total)
		sum_buried_residues += buried_aa_count[aa] #N(total)s
		propensity_buried = (buried_aa_count[aa]/all_aa_count[aa])/(sum_buried_residues/sum_total_residues)
		print(aa, propensity_buried, file=f)





	#structure type = buried or not.

	### END CODING HERE

def main():
	# check input directory
	d_in = argv[1]
	if not os.path.isdir(d_in):
		print(d_in, 'is not a directory')
		exit(1)
	
	unfolded_acc = read_AccUnfold()
	all_aa_count, buried_aa_count = read_dir(d_in, unfolded_acc)
	print_propensities(all_aa_count, buried_aa_count, os.path.join('Output', 'propensity_buried.txt'))


if __name__ == '__main__':
	main()