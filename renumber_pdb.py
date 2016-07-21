#!/bin/env python3
# -*- coding: utf-8 -*-

from Bio import AlignIO,SeqIO,ExPASy,SwissProt,Seq,SeqRecord
from Bio.Alphabet import IUPAC
from Bio.Emboss.Applications import NeedleCommandline
from prody.proteins import parsePDB, writePDB
import argparse
import os
import sys
import re
import tempfile
from gpdb import *
import gpdb


def renumber_noInputAlign(pdbfile,refseqfile,selection="protein",\
	outfile="renumbered.pdb",newAA=None,first=1):
	'''
	Renumber pdb file (pdbfile) according to reference sequence in refseqfile. 
	Pdb sequence is extracted and aligned with reference sequence using needle 
	from EMBOSS.
	- refseqfile: .fasta file containing the reference sequence by which to 
	renumber
	- selection: atom selection(s) in the the structure file to renumber. 
	Will iterate over comma separated selections to renumber each.
	- pdbfile: original structure file
	- outfile: output structure file
	- newAA: comma separated list of unrepresented amino acids
		XXXYCA: 
		XXX = three letter abbrevation as in pdbfile
		Y = one letter code in the alignment
		CA = atom to use as CA if different from "CA", eg 
		C1 in PVL of 1JEN	

	'''
	# selections = selection.split(",")
	selections = selection
	tmp=tempfile.gettempdir()
	tmp_refseqfile="%s/refseq.fasta"%tmp
	pdbID = re.search("\w+\.\w+", pdbfile).group(0)
	tmp_pdbseqfile="%s/%s.fasta"%(tmp,pdbID)
	tmp_needle="%s/needle.out"%tmp
	if os.path.exists(refseqfile):
		refseqRec = SeqIO.read(refseqfile,"fasta",alphabet=IUPAC.protein )
		refseqRec.id = "refseq"
		SeqIO.write(refseqRec,tmp_refseqfile,"fasta")
	else: 
		print ("ERROR, no such file: %s"%refseqfile)
		exit(1)

	if os.path.exists(pdbfile):
		structure=parsePDB("%s"%pdbfile)
		updateAA(structure,newAA)
	else:
		print ("ERROR, no such file: %s"%pdbfile)
		exit(1)

	modified_selections = []
	for polymer in selections:
		currentSel = structure.select("protein and name CA and %s"%polymer)
		if currentSel:
			pdbseq_str=''.join([oneletter[i] for i in currentSel.getResnames()])
			pdbseqRec=SeqRecord(Seq(pdbseq_str,IUPAC.protein),id=pdbID)
			SeqIO.write(pdbseqRec,tmp_pdbseqfile,"fasta")

			needle_cli = NeedleCommandline(asequence=tmp_pdbseqfile,bsequence=tmp_refseqfile,\
				gapopen=10,gapextend=0.5,outfile=tmp_needle)
			needle_cli()
			aln = AlignIO.read(tmp_needle, "emboss",alphabet=IUPAC.protein )
			# os.remove(tmp_needle)
			# os.remove(tmp_pdbseqfile)		

			gpdb.renumber_aln(aln,"refseq",pdbID,first)
			pdbRenSeq = gpdb.seqbyname(aln, pdbID)
			gpdb.renumber_struct(structure, pdbRenSeq,polymer)
			pdbRenSeq.annotations["resnum"]=str(pdbRenSeq.letter_annotations["resnum"])
			modified_selections.append(polymer)
			# seems to be the only way to store pret residue annotations
			# AlignIO.write(aln,"pdb.outseq","seqxml")		
		else:
			print ('ERROR: Selection \"%s\" has zero CA atoms'%polymer)

	if writePDB(outfile, structure):
		print ("Wrote renumbered %s selections from %s to %s"%\
				(str(modified_selections),pdbfile,outfile))
	os.remove(tmp_refseqfile)

def renumber_InputAlign(alnfile,pdbid,refid,selection="protein"\
	,outfile="renumbered.pdb",pdbfile="",newAA=None,first=1):
	'''
	Renumber input pdb using an exsiting multiple alignment. 
	- alnfile: alignment in .fasta format. Beware of weird 	characters in the 
	sequence ids, eg "|"
	- pdbid: sequence id in the alginment file that corresponds	to the input 
	structure. Must be the same number of residues
	- refid: sequence id corresponding to the reference sequence by which to 
	renumber the pdbid sequence. pdbid musnt' align to any gaps in refid.
	- selection: atom selection(s) in the the structure file to renumber. 
	Will iterate over comma separated selections to renumber each.
	- pdbfile: original structure file
	- outfile: output structure file
	- newAA: comma separated list of unrepresented amino acids
		XXXYCA: 
		XXX = three letter abbrevation as in pdbfile
		Y = one letter code in the alignment
		CA = atom to use as CA if different from "CA", eg 
		C1 in PVL of 1JEN	 

	'''


	selections = selection.split(",")
	tmp=tempfile.gettempdir()

	modified_selections = []

	if os.path.exists(alnfile):
		aln = AlignIO.read(alnfile, "fasta",alphabet=IUPAC.protein)
	else:
		print("ERROR, no such alignment: %s"%alnfile)
		exit(1)

	aln_ids = [x.id for x in aln]
	if pdbid in aln_ids and refid in aln_ids:		
		pdbSeqRec = seqbyname(aln, pdbid)
		if not pdbSeqRec:
			print("ERROR, bad pdbid name")
			exit(1)

		refSeqRec = seqbyname(aln, refid)
		if not refSeqRec:
			print("ERROR, bad refid name")
			exit(1)

		if pdbfile != '':
			if os.path.exists(pdbfile):
				structure = parsePDB(pdbfile)
				updateAA(structure,newAA)
			else:
				print("ERROR, no such pdb file: %s"%pdbfile)
				exit(1)

		renumber_aln(aln, refid, pdbid,first)
		for polymer in selections:
			currentSel = structure.select("not hetero and protein and name CA and %s"%polymer)
			if currentSel:
				renumber_struct(structure, pdbSeqRec, polymer)
				modified_selections.append(polymer)
			else:
				print('ERROR: Selection \"%s\" has zero CA atoms'%polymer)
	else:
		if pdbid not in [x.id for x in aln]:
			print("ERROR, no such sequence to renumber: %s"%pdbid)
		if refid not in [x.id for x in aln]:
			print("ERROR, no such sequence to renumber by: %s"%refid)
		exit(1)

	if writePDB(outfile, structure):
		print("Wrote renumbered %s selections from %s to %s"\
				%(str(modified_selections),pdbfile,outfile))

def updateAA(struct,newAAstr):
	if newAAstr:
		newAAs = newAAstr.split(',')
		for newAA in newAAs:
			# if newAA:
			if len(newAA) > 4:				
				updateAATable(struct,newAA[0:3], newAA[3],CA=newAA[4:])
			else:
				updateAATable(struct,newAA[0:3], newAA[3])		


def main():
	parser = argparse.ArgumentParser()
	parser.add_argument("-s","--structure",type=str,help="PDB file to be renumbered. \
		Defaults to <pdbseq>.pdb when used with \'-a\'")
	parser.add_argument("-a","--alignment",type=str,help="Multiple alignment to use for \
		renumbering (fasta format)")	
	parser.add_argument("-r","--refseq",type=str,help="Reference sequence id (for multiple \
		alginment) or .fasta file of reference sequence by which to renumber. Required.",\
		required=True)
	parser.add_argument("-p","--pdbseq",type=str,help="Sequence id to renumber if using an \
		existing multiple alignment")
	parser.add_argument("-v","--selections",type=str,help="Comma separated list of vmd \
		atom selections in double quotes. Each selection will be renumbered according to the \
		alignment",default="protein",nargs='*')
	parser.add_argument("-o","--outfile",type=str,help="Output .pdb filename. Defaults to -s or -p input +\".renumbered.pdb\"")
	parser.add_argument("-n","--newres",type=str,help="Add a new residue to the table of \
		non-standard amino acids: XXXYZ[Z]. XXX = three-letter abbreviation, Y = one-letter \
		abbreviation, Z[Z] = atom to relabel as CA (if needed)",default=None)
	parser.add_argument("-f","--first",type=int,help="Index of first residue", default=1)

	args = parser.parse_args()
	kwargs = {}

	if args.pdbseq and not args.outfile:
		kwargs['outfile'] = "%s.renumbered.pdb"%args.pdbseq
	elif args.outfile:
		kwargs['outfile'] = args.outfile
	elif not args.outfile:
		# inputStruct = args.structure
		if not args.pdbseq:
			pdbID = re.search("\w+\.\w+", args.structure).group(0)		
			kwargs['outfile'] = re.sub('.pdb','',pdbID)+".renumbered.pdb"
		else:
			kwargs['outfile'] = re.sub('.pdb','',args.pdbseq)+".renumbered.pdb"

	kwargs['newAA'] = args.newres
	kwargs['selection'] = args.selections
	kwargs['first'] = args.first
	if args.alignment:		
		if not args.pdbseq or not args.refseq:
			if not args.pdbseq:
				print("No pdb sequence ID specified")
			if not args.refseq:
				print("No reference sequence ID specified")
			print("ERROR, must specify a reference and structure sequence ID")
			exit(1)
		
		if args.structure:
			kwargs['pdbfile'] = args.structure			
		else:
			kwargs['pdbfile'] = "%s.pdb"%args.pdbseq			
	
		renumber_InputAlign(args.alignment, args.pdbseq, args.refseq,**kwargs)
	elif args.structure and args.refseq:		
		renumber_noInputAlign(args.structure, args.refseq, **kwargs)
	else:
		print("ERROR, must specify a structure or alignment file")
		exit(1)

if __name__ == '__main__':
	main()
