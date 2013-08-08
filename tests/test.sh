#!/bin/env zsh
~/bin/align_renumber/renumber_pdb.py -s L812M.B99990003_fit.pdb  -v "chain B","chain D" -r GluN2A_rattusN.fasta -o test.pdb
~/bin/align_renumber/renumber_pdb.py -s test.pdb -v "chain A","chain C" -r GluN1-4A_rattusN.fasta -o L812M.B99990003_fit.renumbered.pdb
diff -q L812M.B99990003_fit.renumbered.pdb outputs/L812M.B99990003_fit.renumbered.pdb
print -l

~/bin/align_renumber/renumber_pdb.py -s 3KG2.pdb -v "chain A","chain B","chain C","chain D" -r AMPA_rat_nosignal.fasta -o 3KG2.renumbered.ns.pdb
~/bin/align_renumber/renumber_pdb.py -s 3KG2.pdb -v "chain A","chain B","chain C","chain D" -r AMPA_rat_nosignal.fasta -o 3KG2.renumbered.ns.pdb
diff -q 3KG2.renumbered.ns.pdb outputs/3KG2.renumbered.ns.pdb
print -l

~/bin/align_renumber/renumber_pdb.py -s 1jen.pdb -v  "chain B","chain A" -r adometdc_hs.fasta -o 1JEN.renumbered.pdb --newres PYROC1
diff -q 1JEN.renumbered.pdb outputs/1JEN.renumbered.pdb
print -l

~/bin/align_renumber/renumber_pdb.py -s 1jen.pdb -v  "chain A B" -r adometdc_hs.fasta -o 1JEN.renumbered2.pdb --newres PYROC1
diff -q 1JEN.renumbered2.pdb outputs/1JEN.renumbered2.pdb
print -l

~/bin/align_renumber/renumber_pdb.py -s 2hfb.pdb -v "chain A","chain B" -r 2hfb.fasta -o 2hfb.renumbered.pdb
diff -q 2hfb.renumbered.pdb outputs/2hfb.renumbered.pdb
print -l
