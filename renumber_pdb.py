#!env python
# -*- coding: utf-8 -*-

import argparse
from Bio import AlignIO,SeqIO,ExPASy,SwissProt
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio.Alphabet import IUPAC
from Bio.Emboss.Applications import NeedleCommandline
from prody.proteins.pdbfile import parsePDB, writePDB
import os
import sys
import tempfile
from gpdb import *
import gpdb

def renumber_noInputAlign(pdbfile,refseqfile,selection="protein",outfile="tmp.pdb"):
	'''
	Renumber pdb file (pdbfile) according to reference sequence in refseqfile. 
	Pdb sequence is extracted and aligned with reference sequence using needle 
	 from EMBOSS.

	'''
	selections = selection.split(",")
	tmp=tempfile.gettempdir()
	tmp_refseqfile="%s/refseq.fasta"%tmp
	tmp_pdbseqfile="%s/%s.fasta"%(tmp,pdbfile)
	tmp_needle="%s/needle.out"%tmp
	if os.path.exists(refseqfile):
		refseqRec = SeqIO.parse(refseqfile,"fasta",alphabet=IUPAC.protein ).next()	
		refseqRec.id = "refseq"
		SeqIO.write(refseqRec,tmp_refseqfile,"fasta")
	else: 
		print "ERROR, no such file: %s"%refseqfile
		exit(1)

	if os.path.exists(pdbfile):
		structure=parsePDB("%s"%pdbfile)
	else:
		print "ERROR, no such file: %s"%pdbfile
		exit(1)

	modified_selections = []
	for polymer in selections:
		currentSelection = structure.select("not hetero and protein and name CA and %s"%polymer)
		if currentSelection:
			pdbseq_str=''.join([oneletter[i] for i in currentSelection.getResnames()])
			pdbseqRec=SeqRecord(Seq(pdbseq_str,IUPAC.protein),id=pdbfile)
			SeqIO.write(pdbseqRec,tmp_pdbseqfile,"fasta")

			needle_cli = NeedleCommandline(asequence=tmp_pdbseqfile,bsequence=tmp_refseqfile,gapopen=10,gapextend=0.5,outfile=tmp_needle)
			needle_cli()
			aln = AlignIO.read(tmp_needle, "emboss",alphabet=IUPAC.protein )
			# os.remove(tmp_needle)
			os.remove(tmp_pdbseqfile)		

			gpdb.renumber_aln(aln,"refseq",pdbfile)
			pdbRenSeq = gpdb.seqbyname(aln, pdbfile)		
			gpdb.renumber_struct(structure, pdbRenSeq,polymer)
			# pdbRenSeq = gpdb.seqbyname(aln, pdbfile)
			# seems to be the only way to store per residue annotations
			pdbRenSeq.annotations["resnum"]=str(pdbRenSeq.letter_annotations["resnum"])
			# AlignIO.write(aln,"pdb.outseq","seqxml")		
			modified_selections.append(polymer)
		else:
			print 'ERROR: Selection \"%s\" has zero CA atoms'%polymer

	if writePDB(outfile, structure):
		print "Wrote renumbered %s selections from %s to %s"%(str(modified_selections),pdbfile,outfile)
	os.remove(tmp_refseqfile)

def renumber_InputAlign(alnfile,pdbid,refid,selection="protein",outfile="renumbered.pdb",pdbfile=""):
	selections = selection.split(",")
	tmp=tempfile.gettempdir()

	modified_selections = []

	if os.path.exists(alnfile):
		aln = AlignIO.read(alnfile, "fasta",alphabet=IUPAC.protein)
	else:
		print "ERROR, no such alignment: %s"%alnfile
		exit(1)

	aln_ids = [x.id for x in aln]
	if pdbid in aln_ids and refid in aln_ids:		
		pdbSeqRec = seqbyname(aln, pdbid)
		refSeqRec = seqbyname(aln, refid)

		if pdbfile != '':
			if os.path.exists(pdbfile):
				structure = parsePDB(pdbfile)
			else:
				print "ERROR, no such pdb file: %s"%pdbfile
				exit(1)

		if pdbSeqRec and refSeqRec and structure:
			renumber_aln(aln, refid, pdbid)
			for polymer in selections:
				currentSelection = structure.select("not hetero and protein and name CA and %s"%polymer)
				if currentSelection:
					# renumber_struct(structure, pdbSeqRec)
					renumber_struct(structure, pdbSeqRec, polymer)
					modified_selections.append(polymer)
				else:
					print 'ERROR: Selection \"%s\" has zero CA atoms'%polymer										
	else:
		if pdbid not in [x.id for x in aln]:
			print "ERROR, no such sequence to renumber: %s"%pdbid		
		if refid not in [x.id for x in aln]:
			print "ERROR, no such sequence to renumber by: %s"%refid
		exit(1)

	if writePDB(outfile, structure):
		print "Wrote renumbered %s selections from %s to %s"%(str(modified_selections),pdbfile,outfile)

def main():
	parser = argparse.ArgumentParser()
	parser.add_argument("-s","--structure",type=str,help="PDB file to be renumbered. Defaults to <pdbseq>.pdb when used with \'-a\'")
	parser.add_argument("-a","--alignment",type=str,help="Multiple alignment to use for renumbering (fasta format)")	
	parser.add_argument("-r","--refseq",type=str,help="Reference sequence id (for multiple alginment) or .fasta file of reference sequence"\
		" according to which to renumber. Required.",required=True)
	parser.add_argument("-p","--pdbseq",type=str,help="Sequence id to renumber if using an existing multiple alignment.")
	parser.add_argument("-v","--selections",type=str,help="Comma separated list of vmd atomselections in "\
		"double quotes. Each selection will be renumbered according to the alignment",default="protein")
	parser.add_argument("-o","--outfile",type=str,help="Output .pdb filename",default="renumbered.pdb")

	# combinations: -a -s -p -r / -s -r -v

	args = parser.parse_args()
	kwargs = {}

	# if args.outfile:
	# 		kwargs['outfile'] = args.outfile
	if args.pdbseq and not args.outfile:
		kwargs['outfile'] = "%s.renumbered.pdb"%args.pdbseq
	elif args.outfile:
		kwargs['outfile'] = args.outfile

	if not args.selections:
		args.selections = "protein"

	if args.alignment:		
		if not args.pdbseq or not args.refseq:
			if not args.pdbseq:
				print "No pdb sequence ID specified"
			if not args.refseq:
				print "No reference sequence ID specified"
			print "ERROR, must specify a reference and structure sequence ID"
			exit(1)
		
		if args.structure:
			kwargs['pdbfile'] = args.structure			
		else:
			kwargs['pdbfile'] = "%s.pdb"%args.pdbseq			
		kwargs['selection'] = args.selections
		renumber_InputAlign(args.alignment, args.pdbseq, args.refseq,**kwargs)
	elif args.structure and args.refseq:		
		kwargs['selection'] = args.selections
		renumber_noInputAlign(args.structure, args.refseq, **kwargs)
	else:
		print "ERROR, must specify a structure or alignment file"
		exit(1)
	# doesn't yet cater for reading 
if __name__ == '__main__':
	main()
