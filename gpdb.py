# -*- coding: utf-8 -*-
from prody import *
from Bio import SeqIO
from Bio.Alphabet import IUPAC
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio import AlignIO,SeqIO

oneletter = {
'ASP':'D','GLU':'E','ASN':'N','GLN':'Q',
'ARG':'R','LYS':'K','PRO':'P','GLY':'G',
'CYS':'C','THR':'T','SER':'S','MET':'M',
'TRP':'W','PHE':'F','TYR':'Y','HIS':'H',
'ALA':'A','VAL':'V','LEU':'L','ILE':'I',
'ASX':'B','GLX':'Z','CSO':'C','HIP':'H',
'HSD':'H','HSE':'H','HSP':'H','MSE':'M',
'SEC':'U','SEP':'S','TPO':'T','PTR':'Y',
'XLE':'J','XAA':'X'
}

'''
Non-standard amino acids can be extended in prody with AddNonstdAminoacid()
Must have a CA atom (can set C1 of PYR to CA, but must be done before extending the list)

What is the actual IUPAC one letter code for selenothemionine etc?
This one seems to be complete:
https://raw.github.com/Pymol-Scripts/Pymol-script-repo/master/findseq.py

TODO:
Looks like this will required some changes to prody or subclassing:
- preserve altlocs in output
- preserve anisou in output
- preserve secondary structure in output
'''

def uniq(seq):
    '''
    via Markus Jarderot on Stack Exchange
    Was f7, lifted from http://goo.gl/mQHgi
    '''
    seen = set()
    seen_add = seen.add 
    return [x for x in seq if x not in seen and not seen_add(x)]

def updateAATable(struct,threeletterAbbrv,oneletterAbbrv,CA='CA'):
	'''
	Update table of nonstandard amino acids inherited from prody and oneletterAbbrv dictionary.
	If corresponding residues don't have a CA at set another to this first. Eg set C1 to CA in pyr of adometdc

	'''

	print( "Adding non-standard amino acid %s:%s, %s set to CA"%(threeletterAbbrv,oneletterAbbrv,CA))
	struct.select("resname %s and name %s"%(threeletterAbbrv,CA)).setNames("CA")
	addNonstdAminoacid(threeletterAbbrv)
	oneletter.update({threeletterAbbrv:oneletterAbbrv})
	pass

def seqnum(aln,seqid):
	'''
	Return index of sequence with corresponding sequence id
	'''
	l = [i for i in range(len(aln)) if aln[i].id == seqid]
	return l[0]

def seqbyname(aln,seqid):
	'''
	Return sequence with corresponding sequence id
	'''
	seqRec = [s for s in aln if s.id == seqid][0]
	return seqRec


def renumber_aln(aln,refseq_id,pdbseq_id,first=1):
	'''
	Renumber pdbseq_id in aln according to "resnum" values of refseq_id. 
	Also recorded in "resnum" letter annotation
	
	'''
	pdbseqRec = [s for s in aln if s.id == pdbseq_id][0]
	pdbseqRec.letter_annotations["resnum"]=[None]*len(pdbseqRec)
	refseqRec = [s for s in aln if s.id == refseq_id][0]
	refseqRec.letter_annotations["resnum"]=[None]*len(refseqRec)
	refseqNogap = [i for i in range(len(refseqRec)) if refseqRec[i] != '-']
	resnum = first
	for resindex in refseqNogap:
		refseqRec.letter_annotations["resnum"][resindex]=resnum
		resnum+=1
	
	reslist = [[i,refseqRec.letter_annotations["resnum"][i]] for i in range(len(refseqRec)) if pdbseqRec[i] != '-']
	for [i,r] in reslist:
		pdbseqRec.letter_annotations["resnum"][i]=r

	# print pdbseqRec.letter_annotations["resnum"]

def overlap(aln,seqid1,seqid2):
	'''
	Return list of indices corresponding to aligned residues with assgined residue numbers	
	
	'''
	seq1 = seqbyname(aln,seqid1)
	seq2 = seqbyname(aln,seqid2)

	l = [i for i in range(len(seq1)) if seq1.letter_annotations["resnum"][i] != None and seq2.letter_annotations["resnum"][i] != None]
	return l

def vmdsliceToReslist (vmdslice):
	'''
	Convert resnum selection in vmd to list

	'''
	words = vmdslice.split(" ")
	numlist = []
	last = None
	inrange = False
	for w in words:
		if w != "to" and not inrange:
			last = int(w)
			numlist.append(last)
		elif w != 'to' and inrange:
			for x in range(last,int(w)+1):
				numlist.append(x)
			inrange = False
		elif w == "to":
			inrange = True

	return numlist

def vmdslice(l):
	'''
	Return list of values as compressed slice in vmd selection algebra

	'''
	i = l[0]
	begin = i
	# s = str(i)
	s = ""
	in_interval=0
	for j in l[1:]:
		if j-i == 1 :
			in_interval = 1
			end = j
		elif j -i > 1:			
			if in_interval == 1:
				end = i
				if end - begin >= 2:
					s = s + "%d to %d "%(begin,end)
				else:
					s = s + "%d %d "%(begin,end)		
				begin = j
				in_interval = 0
			else:
				begin = j
				s = s + "%d "%i
		i = j

	if end - begin >= 2:
		s = s + "%d to %d "%(begin,end)
	elif end - begin == 1:
		pass
		s = s + " %d %d"%(begin, end)				
	else:
		s = s + " %d"%(begin)		
	# print l
	return s

def pymol_slice(l):
	'''
	Return list of values as compressed slice in pymol selection algebra

	'''
	i = l[0]
	begin = i
	# s = str(i)
	s = ""
	in_interval=0
	for j in l[1:]:
		if j-i == 1 :
			in_interval = 1
			end = j
		elif j -i > 1:			
			if in_interval == 1:
				end = i
				if end - begin >= 2:
					s = s + "%d-%d,"%(begin,end)
				else:
					s = s + "%d,%d,"%(begin,end)		
				begin = j
				in_interval = 0
			else:
				begin = j
				s = s + "%d,"%i
		i = j

	if end - begin >= 2:
		s = s + "%d-%d,"%(begin,end)
	else:
		s = s + ",%d"%(begin)		
	
	return s

def increment_clash(newresnums,resindices,struct,renumbered_sel="protein"):
	'''
	Fix clashes by incrementing offending residues by the largest previous new residue number.
	This is a bit clumsy, need to still check for clashes after increments and some interleaved 
	residues keep old values due to protein gaps. Also need to iterate over chains.

	First renumber waters/bulk solvent (Cl, Na) from last 

	'''
	# print uniq(struct.getResindices())
	proteinResindices = ' '.join(str(i) for i in resindices)
	nonProteinResnums = uniq(struct.select("not resindex %s"%proteinResindices).getResnums())
	
	intersection = [x for x in newresnums if x in nonProteinResnums]

	if len(intersection) > 0:
		# print intersection
		intersectionids = [uniq(struct.select("resnum %d"%x).getResindices())[0] for x in newresnums if x in nonProteinResnums]
	# 	# nonProteinResindices = struct.select("not protein").getResindices()
		lastProtein = newresnums[-1]
		# print lastProtein
		newclashNums = []
		for i in range(len(intersection)):
		 	atoms = struct.select("resindex %d"%intersectionids[i]).numAtoms()
		 	offsetNum = intersection[i]+lastProtein
		 	newclashNums.append(offsetNum*atoms)
		clashnum_str = ' '.join(str(i) for i in intersectionids)
		# print clashnum_str
		struct.select("resindex %s"%clashnum_str).setResnums(newclashNums)
		# print newclashNums
	# return newclashNums

def sequential_renumber(struct,resindex,resindex_str,selection,number_from):
	newclashNums = []
	last_biggest = number_from
	
	for resn in resindex:
		last_biggest += 1
		natoms = struct.select("resindex %d and (%s)"%(resn,selection)).numAtoms()
		newclashNums += [last_biggest]*natoms
	# print "DEBUG %s"%newclashNums

	struct.select("resindex %s"%resindex_str).setResnums(newclashNums)	
	return last_biggest

def fix_clash(newresnums,resindices,struct,renumbered_selstr="chain A"):
	'''
	This assumes sequential numbering in original pdb. Leave original ligand resnums untouched
	unless they clash, therefore might not be sequential between protein and solvent in final
	output. I tend to rely on these staying the same. Might	do funny things if solvent and ligands
	appear before protein, haven't checked.

	Does not yet deal with multiple chains. Also need to cater for modified amino acids 

	Note: HETATM records that contain standard amino acids not recognised as hetero in selection 
	algebra. Use "hetatm" instead of "hetero"

	'''

	# print struct.select("resindex 810 and ((water or name CL NA) and chain A and not protein)").numAtoms()
	chains_inmodified = [struct.select("resindex %d"%(i)).getChids()[0] for i in resindices]
	# chains_inmodified = [uniq(struct.select("resindex %d"%(i)).getChids())[0] for i in resindices]
	for chain in uniq(chains_inmodified):
		# iterate over chains, clashes between chains are not a problem
		newresnums_ch = [newresnums[i] for i in range(len(newresnums)) if chains_inmodified[i] == chain]
		# print "DEBUG "+str(newresnums_ch)

		solvent_selstr = "(water or name CL NA) and chain %s and not protein"%(chain)
		# not sure what else to treat as bulk solvent, POT?
		hetero_selstr = "(hetatm) and not (%s) and chain %s and not protein"%(solvent_selstr,chain)
		# use hetatm instead of hetero to pick up amino acid ligands
		hetero_intersection = []
		solvent_intersection = []
		last_biggest = newresnums_ch[-1] 

		hetero = struct.select(hetero_selstr)
		if hetero:
			ligands = ', '.join(i for i in uniq(hetero.getResnames()))
			print ("Found the following ligands in chain %s: %s"%(chain,ligands))
			# print uniq(hetero.getChids())
			# print "DEBUG"%uniq(hetero.getResnums())
			hetero_resids = uniq(hetero.getResindices())
			hetero_resids_str = ' '.join(str(i) for i in uniq(hetero.getResindices()))
			hetero_resnums = uniq(hetero.getResnums())
			hetero_intersection = [x for x in newresnums if x in hetero_resnums]

		solvent = struct.select(solvent_selstr)
		if solvent:
			solvent_resids = uniq(solvent.getResindices())
			# print solvent_resids
			solvent_resids_str = ' '.join(str(i) for i in uniq(solvent.getResindices()))
			solvent_resnums = uniq(solvent.getResnums())
			solvent_intersection = [x for x in newresnums if x in solvent_resnums]

		if len (hetero_intersection) > 0:
			print( "WARNING, ligand resnums clash with new resnums")
			print( "Renumbering ligands from %d for %s"%(last_biggest+1,renumbered_selstr))
			# print "DEBUG %s"%hetero_resids
			# print "DEBUG %s"%hetero_resnums
			last_biggest = sequential_renumber(struct,hetero_resids,hetero_resids_str,hetero_selstr,last_biggest)

		if len(solvent_intersection) > 0:
			hetero_newsolvent_intersection = []
			print( "WARNING, solvent resnums clash with new resnums")

			# check for clash between previously ignored ligands and new solvent numbering
			if len(hetero_intersection) > 0 and hetero:
				hetero_newsolvent_intersection = [x for x in hetero_resnums if x in range(last_biggest,last_biggest+len(solvent_resids))]
			if len(hetero_newsolvent_intersection) > 0:
				print ("WARNING, solvent renumbering clashes with ligands")
				print ("Renumbering ligands from %d for %s"%(last_biggest+1,renumbered_selstr))
				last_biggest = sequential_renumber(struct,hetero_resids,hetero_resids_str,hetero_selstr,last_biggest)
			
			print ("Renumbering solvent from %d"%(last_biggest+1))
			sequential_renumber(struct,solvent_resids,solvent_resids_str,solvent_selstr,last_biggest)	

def renumber_struct(struct,seq,selection="chain A"):
	'''
	Renmuber residues in struct according to "resnum" letter annotation of seq. If 
	new numbers conflict with existing resnums of ligands and solvent these must be
	renumbered first. Otherwise, indexing in prody gets confused and only protein 
	resindices can be accessed. Therefore fixing clashes cannot be done after 
	renumbering proteins.

	'''
	# print struct.getTitle()
	currentSelection = struct.select("protein and name CA and (%s) and not hetero"%selection)
	resindices = currentSelection.getResindices()	
	newresnums=[i for i in seq.letter_annotations["resnum"][:] if i != None]
	num_resindices = len(resindices)
	num_newresnums = len(newresnums)

	fix_clash(newresnums,resindices,struct,selection)
	# increment_clash(newresnums,resindices,struct,"all")	

	if num_newresnums == num_resindices and num_resindices > 0:
		resindex_str = ' '.join(str(i) for i in resindices)
		newresnums_long = []
		# [[newresnums[i]]*struct.select("resindex %d"%resindices[i]).numAtoms()  for i in range(len(newresnums))]

		resmatrix = [[newresnums[i],resindices[i]] for i in range(len(newresnums)) ]
		for [newresnum,resindex] in resmatrix:
		# 	# print "%d %d"%(resindex,newresnum)
		# 	struct.select("resindex %d"%resindex).setResnums(newresnum)
		# 	# this is quite slow
			newresnums_long += [newresnum]*struct.select("resindex %d"%resindex).numAtoms() 
		struct.select("resindex %s"%resindex_str).setResnums(newresnums_long) #about 20× faster
	else:
		print( 'ERROR, unequal number of renumbered residues (%d) and \n'\
			'to-be-renumbered residues (%d)\n'%(num_newresnums,num_resindices))
		print( 'Make certain the correct reference sequence was chosen\n'+\
		 		'and that none of the structure sequence is aligned to \n'+\
		 		'gaps. Possible cause: affinity tags etc. Consider using\n'+\
		 		'the "and not sequence XYZ" selector to exclude linkers etc.')

		exit(1)

def write_renumbered(aln,pdbid,pchain):
	nr1 = seqbyname(aln,"%s_%s"%(pdbid,pchain))
	nr1_struct = parsePDB(pdbid,chain=pchain)
	renumber_struct(nr1_struct, "chain %s"%pchain)
	writePDB("%s_%s.renumbered.pdb"%(pdbid,pchain),nr1_struct)


def showoverlap(seqid1,seqid2,aln):
	seq1 = seqbyname(aln, seqid1)
	seq2 = seqbyname(aln, seqid2)
	overlapping = overlap(aln, seqid1, seqid2)
	seq1_overlap = [seq1.letter_annotations["resnum"][i] for i in overlapping] 
	seq2_overlap = [seq2.letter_annotations["resnum"][i] for i in overlapping] 
	seq1_overlap_resname = [seq1.seq[i] for i in overlapping]
	seq2_overlap_resname = [seq2.seq[i] for i in overlapping]

	slice1 = seq1_overlap
	slice1str= ' '.join(str(i) for i in slice1)
	slice2 = seq2_overlap
	slice2str= ' '.join(str(i) for i in slice2)

	print( "%s : %s"%(seqid1,seqid2))

	for i in range(len(seq1_overlap)):
		print ("%s %s : %s %s "%(seq1_overlap[i],seq1_overlap_resname[i],seq2_overlap[i],seq2_overlap_resname[i]))
	


