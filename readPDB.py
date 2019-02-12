#!/usr/bin/python

"""
Structural Bioinformatics assignment 1 - Phi-psi angles

When you finish the assignment, the script should create an
output file containing the phi and psi angles and secondary
structure assignment of each residue. Please ONLY modify the
code in the three indicated blocks and do NOT import additional
python packages.

To run, make 'Assignment1' your working directory and use:

> python Scripts/readPDB.py Data/pdb_filename
"""

# Packages
from sys import argv
import os
from math import sqrt, atan2, degrees

# Vector functions that we need to calculate the angles
def dot_product(v1, v2):
	""" Calculate the dot product of two vectors """
	return v1[0]*v2[0]+v1[1]*v2[1]+v1[2]*v2[2]

def cross_product(v1, v2):
	""" Calculate the cross product of two vectors """
	i = v1[1]*v2[2] - v1[2]*v2[1]
	j = v1[2]*v2[0] - v1[0]*v2[2]
	k = v1[0]*v2[1] - v1[1]*v2[0]
	return [i,j,k]
	
def magnitude(v):
	""" Calculate the size of a vector """
	return sqrt(v[0]**2 + v[1]**2 + v[2]**2)

# PDB file parser
def readPDB(PDB_file):
	""" Reads a PDB file and stores the atom
	coordinates and amino acid types of the protein """
	# open the file
	try:
		f = open(PDB_file, 'r')
	except:
		print('Error: cannot open', PDB_file)
		
	# dictionaries to store the output
	# pdb atom coordinates:
	#     pdbcoord[chain][residue_number][atom_type] = coordinates
	pdbcoord = {}
	# residue type per chain and residue number (i.e. store the sequence)
	#     pdbseq[chain][resnum] = restype
	pdbseq = {}
	
	# parse each line in the file
	for line in f:
		# remove whitespace at the end of the line
		line = line.strip()
		# only parse the lines containing atom coordinates
		if line[:4] == 'ATOM':
			# ATOM type (e.g. C-alpha)
			atom_type = line[12:16].strip()
			# AMINO ACID type (e.g. alanine)
			aa_type = line[17:20].strip()
			# residue number
			res_num = int(line[22:26])
			# Protein chain
			chain = line[21]
			# coordinates
			xcoord = float(line[30:38])
			ycoord = float(line[38:46])
			zcoord = float(line[46:54])
			
			# if chain does not exists create new entry
			if not chain in pdbcoord:
				pdbcoord[chain] = {}
				pdbseq[chain] = {}
			# if resnum does not exists create new entry
			if not res_num in pdbcoord[chain]:
				pdbcoord[chain][res_num] = {}

			# store coordinates as a vector
			pdbcoord[chain][res_num][atom_type] = [xcoord,ycoord,zcoord]			
			# store sequence
			pdbseq[chain][res_num] = aa_type
			
	# close file
	f.close()
	
	# return dictionaries
	return pdbcoord, pdbseq

### THE FOLLOWING THREE FUNCTIONS ARE THE ONES YOU NEED
### TO EDIT FOR THE ASSIGNMENT. ONLY EDIT THE INDICATED
### BLOCKS
def calculateDihedral(a1, a2, a3, a4):
	""" Calculates the normal vector of the planes
	defined by four atom coordinates """
	### START CODING HERE
	# calculate normal vectors to the planes defined by a1,a2,a3 and a2,a3,a4
	# you may use the functions "cross_product","dot_product" and "magnitude" defined above
	# you can also use the python math functions "atan2" and "degrees"

	#directional vectors.
	b1 = [a2 - a1 for a1, a2 in zip(a1, a2)]
	b2 = [a3 - a2 for a2, a3 in zip(a2, a3)]
	b3 = [a4 - a3 for a3, a4 in zip(a3, a4)]



	#cross b1 x b2 and b2 x b3
	#cross products

	v1 = cross_product(b1, b2)
	v2 = cross_product(b2, b3)


	u1 = []
	u2 = []
	for coordinate in v1:
		u1.append(coordinate/magnitude(v1))
	for coordinate in v2:
		u2.append(coordinate/magnitude(v2))


	cos_angle = (dot_product(v1, v2))/(magnitude(v1)*magnitude(v2))
	print(cos_angle)
	sin_vector = []
	criss_cross = cross_product(v1, v2)
	for coordinate in criss_cross:
		sin_vector.append(coordinate/magnitude(v1)*magnitude(v2))
	sine_angle = dot_product(sin_vector, u1)
	print(sine_angle)

	dihedral = degrees(atan2(sine_angle, cos_angle))



	
	### END CODING HERE
	
	return dihedral

def assign_ss(phi, psi):
	""" Assign a secondary structure type based on the phi
	and psi angles of a residue """
	### START CODING HERE
	
	### END CODING HERE
	return secondary_structure

def print_phi_psi(pdbcoord, pdbseq, outfile):
	""" given the PDB coordinates, calculate the dihedral
	angles of all the residues, assign secondary structure
	types and write them into an output file """
	f = open(outfile, 'w')

	# get the chains from the PDB file
	list_chains = sorted(pdbcoord.keys())


	for chain in list_chains:
		# get the sorted residue numbers from the pdbcoord dictionary
		list_residue_numbers = sorted(pdbcoord[chain].keys())
		for res_num in list_residue_numbers:
			# if certain residues are missing in the PDB file, you will
			# get a KeyError. Make sure your program does not crash, but
			# gives a warning when this happens
			try:
				### START CODING HERE
				#print(pdbcoord[chain][res_num])
				#vectors needed for directional vectors.
				comp_1 = pdbcoord[chain][res_num-1]["C"]
				comp_2 = pdbcoord[chain][res_num]["CA"]
				comp_3 = pdbcoord[chain][res_num]["C"]
				comp_4 = pdbcoord[chain][res_num]["N"]
				comp_5 = pdbcoord[chain][res_num+1]["N"]

				# #phi directional vectors.
				# phi_b1 = [comp_4 - comp_1 for comp_1, comp_4 in zip(comp_1, comp_4)]
				# phi_b2 = [comp_2 - comp_4 for comp_4, comp_2 in zip(comp_4, comp_2)]
				# phi_b3 = [comp_3 - comp_2 for comp_2, comp_3 in zip(comp_2, comp_3)]
				#
				# #psi directional vectors
				# psi_b1 = [comp_2 - comp_4 for comp_4, comp_2 in zip(comp_4, comp_2)]
				# psi_b2 = [comp_3 - comp_2 for comp_2, comp_3 in zip(comp_2, comp_3)]
				# psi_b3 = [comp_5 - comp_3 for comp_3, comp_5 in zip(comp_5, comp_3)]

				phi = calculateDihedral(comp_1, comp_4, comp_2, comp_3)
				psi = calculateDihedral(comp_4, comp_2, comp_3, comp_5)

				print(phi, "phi")
				print(psi, "psi")


				### END CODING HERE

			except KeyError:
				print('WARNING: KeyError:', KeyError, 'in residue', chain, res_num)

			# get amino acid
			aa_type = pdbseq[chain][res_num]
			# write into output file
			#print(chain, res_num, aa_type, phi, psi, ss, file=f)
	f.close()
	print('written:', outfile)

def main():
	# input PDB file
	f_in = argv[1]
	f_out = os.path.join('Output', 'phi_psi.txt')
	
	# read PDB file
	pdbcoord, pdbseq = readPDB(f_in)
	print_phi_psi(pdbcoord, pdbseq, f_out)
	
if __name__ == '__main__':
	main()