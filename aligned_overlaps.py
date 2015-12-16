#!env python2
# -*- coding: utf-8 -*-

import argparse
from Bio import AlignIO,SeqIO,ExPASy,SwissProt
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio.Alphabet import IUPAC
from prody.proteins.pdbfile import parsePDB, writePDB
import os
import sys
import re
import tempfile
from gpdb import *
import gpdb

def print_overlaps(pdbid1,pdbid2,slice1,slice2):
	print "%s: %s"%(pdbid1,vmdslice(slice1))
	print "%s: %s"%(pdbid2,vmdslice(slice2))
	pass

def writevmd_script(filename,pdbid1,pdbid2,slice1,slice2):
	vmdscript = open(filename,"w")

	vmdscript.write('''
proc heterofit {{mol1 0} {mol2 1} {sel1 "protein"} {sel2 "protein"} {frame now}} {
# mol2 is stationary
    set list1 [[atomselect $mol1 "$sel1 and name CA"] get resid]
    set list2 [[atomselect $mol2 "$sel2 and name CA"] get resid]

    set fit1 [atomselect $mol1 "$sel1 and name CA"]
    set fit2 [atomselect $mol2 "$sel2 and name CA"]

    set moveby [measure fit $fit1 $fit2]

    [atomselect $mol1 "all" frame $frame] move $moveby
}

''')
	vmdscript.write('proc fit_%s_to_%s {{mol1} {mol2} {sel1 \"protein\"} {sel2 \"protein\"}} {\n'%(pdbid1,pdbid2))
	vmdscript.write("# Fit %s to %s \n"%(pdbid1,pdbid2))
	vmdscript.write("\tputs \"Number of atoms to fit in %s:\"\n"%pdbid1)
	vmdscript.write('\tputs [[atomselect $mol1 \"resid %s and name CA and $sel1\" ] num]\n'%vmdslice(slice1))
	vmdscript.write("\tputs \"Number of atoms to fit in %s:\"\n"%pdbid2)
	vmdscript.write('\tputs [[atomselect $mol2 "resid %s and name CA and $sel2" ] num]\n'%vmdslice(slice2))
	vmdscript.write('\theterofit $mol1 $mol2 "resid %s and $sel1" "resid %s and $sel2"\n'%(vmdslice(slice1),vmdslice(slice2)))
	vmdscript.write("}\n\n")

	# vmdscript.write("proc fit_%s_to_%s {mol1 mol2} {\n"%(pdbid2,pdbid1))
	vmdscript.write('proc fit_%s_to_%s {{mol1} {mol2} {sel1 \"protein\"} {sel2 \"protein\"}} {\n'%(pdbid2,pdbid1))
	vmdscript.write("# Fit %s to %s \n"%(pdbid2,pdbid1))
	vmdscript.write("\tputs \"Number of atoms to fit in %s:\"\n"%pdbid2)
	vmdscript.write("\tputs [[atomselect $mol1 \"resid %s and name CA and $sel1\" ] num]\n"%vmdslice(slice2))
	vmdscript.write("\tputs \"Number of atoms to fit in %s:\"\n"%pdbid1)
	vmdscript.write('\tputs [[atomselect $mol2 "resid %s and name CA and $sel2" ] num]\n'%vmdslice(slice1))
	vmdscript.write('\theterofit $mol1 $mol2 "resid %s and $sel1" "resid %s and $sel2"\n'%(vmdslice(slice2),vmdslice(slice1)))
	vmdscript.write("}\n\n")

	vmdscript.close()
	
	print "VMD script output in %s"%filename	

def write_pymol_script(filename,pdbid1,pdbid2,slice1,slice2):

	pymolscript = open(filename,"w")
	pymolscript.write("import pymol\n")

	pymolscript.write('''
def fit_%s_to_%s(mol1,mol2):
	cmd.do("select tmp1, name CA and (%s) and (i. %s)"%s)
	cmd.do("select tmp2, name CA and (%s) and (i. %s)"%s)
	cmd.do("fit tmp1, tmp2, matchmaker=-1")
	cmd.delete("tmp1")
	cmd.delete("tmp2")

cmd.extend "fit_%s_to_%s",fit_%s_to_%s
'''%(pdbid1,pdbid2,"%s",pymol_slice(slice1),"%mol1","%s",pymol_slice(slice2),"%mol2",\
	pdbid1,pdbid2,pdbid1,pdbid2))

	pymolscript.write('''
def fit_%s_to_%s(mol1,mol2):
	cmd.do("select tmp1, name CA and (%s) and (i. %s)"%s)
	cmd.do("select tmp2, name CA and (%s) and (i. %s)"%s)
	cmd.do("fit tmp1, tmp2, matchmaker=-1")
	cmd.delete("tmp1")
	cmd.delete("tmp2")

cmd.extend" fit_%s_to_%s",fit_%s_to_%s
'''%(pdbid2,pdbid1,"%s",pymol_slice(slice2),"%mol1","%s",pymol_slice(slice1),"%mol2",\
	pdbid2,pdbid1,pdbid2,pdbid1))
	
	
	pymolscript.close()

def overlapping(alnfile,pdbid1,pdbid2,refid1,refid2,writevmd="",writepymol="",filter1=None,filter2=None,first1=1,first2=1):
	tmp=tempfile.gettempdir()

	if os.path.exists(alnfile):
		aln = AlignIO.read(alnfile, "fasta",alphabet=IUPAC.protein)
		# aln = AlignIO.read(alnfile, "clustal",alphabet=IUPAC.protein)
	else:
		print "ERROR, no such alignment: %s"%alnfile
		exit(1)

	aln_ids = [x.id for x in aln]

	for sequence in [pdbid1,pdbid2,refid1,refid2] :
		if not sequence in aln_ids:
			print "No such entry in alignment: %s"%sequence
			return

	if pdbid1 in aln_ids and refid1 in aln_ids and pdbid2 in aln_ids and refid2 in aln_ids:		
		renumber_aln(aln, refid1, pdbid1,first=first1)
		renumber_aln(aln, refid2, pdbid2,first=first2)
		common = overlap(aln, pdbid1, pdbid2)

		pdbSeqRec1 = seqbyname(aln, pdbid1)
		# refSeqRec1 = seqbyname(aln, refid1)
		pdbSeqRec2 = seqbyname(aln, pdbid2)
		# refSeqRec2 = seqbyname(aln, refid2)
		
		resnums1 = pdbSeqRec1.letter_annotations["resnum"] 		
		resnums2 = pdbSeqRec2.letter_annotations["resnum"] 

		#convert filters to vmd-like selection string
		if filter1:
			restrict1 = vmdsliceToReslist(filter1)
		else:	
			restrict1 = resnums1
		if filter2: 
			restrict2 = vmdsliceToReslist(filter2)
		else:
			restrict2 = resnums2

		slice1 = [resnums1[i] for i in common if resnums1[i] in restrict1 and  resnums2[i] in restrict2]
		slice2 = [resnums2[i] for i in common if resnums1[i] in restrict1 and  resnums2[i] in restrict2]

		print "Number of residues to align in %s: %s"%(pdbid1,len(slice1))
		print "Number of residues to align in %s: %s"%(pdbid2,len(slice2))

		if writevmd != "" and writevmd != None:
			writevmd_script(writevmd,pdbid1,pdbid2,slice1,slice2)

		if writepymol != "" and writepymol != None:
			write_pymol_script(writepymol,pdbid1,pdbid2,slice1,slice2)
		
		print_overlaps(pdbid1,pdbid2,slice1,slice2)		
		
	else:
		print "we did not arrive"

def main():
	parser = argparse.ArgumentParser()
 	parser.add_argument("-a","--alignment",type=str,help="Multiple alignment to use for renumbering (fasta format)")	
	parser.add_argument("-p1","--pdbseq1",type=str,help="Sequence id of first sequence renumber.")
	parser.add_argument("-p2","--pdbseq2",type=str,help="Sequence id of second sequence renumber.")
	parser.add_argument("-r1","--refseq1",type=str,help="Reference sequence id for first sequence.")
	parser.add_argument("-r2","--refseq2",type=str,help="Reference sequence id for second sequence.")
	parser.add_argument("-v","--writevmd",type=str,help="Output a vmd .tcl script for fitting")
	parser.add_argument("-py","--writepymol",type=str,help="Output a pymol script for fitting")
	parser.add_argument("-s1","--filter1",type=str,help="VMD Filter for first sequence")
	parser.add_argument("-s2","--filter2",type=str,help="VMD Filter for second sequence")
	parser.add_argument("-f1","--first1",type=int,help="First residue index for refseq1",default=1)
	parser.add_argument("-f2","--first2",type=int,help="First residue index for refseq2",default=1)

	args = parser.parse_args()
	kwargs = {
		'writevmd':args.writevmd,
		'writepymol':args.writepymol,
		"filter1":args.filter1,
		"filter2":args.filter2,
		"first1":args.first1,
		"first2":args.first2
	}

	overlapping(args.alignment, args.pdbseq1, args.pdbseq2, args.refseq1, args.refseq2,**kwargs)

if __name__ == '__main__':
	main()

