#!/bin/env zsh
export PYTHONPATH=../:$PYTHONPATH

python align_test.py

diff 1PB7_A.renumbered.pdb outputs/1PB7_A.renumbered.pdb
diff 1PB8_A.renumbered.pdb outputs/1PB8_A.renumbered.pdb
diff 1PB9_A.renumbered.pdb outputs/1PB9_A.renumbered.pdb
diff 1PBQ_A.renumbered.pdb outputs/1PBQ_A.renumbered.pdb
diff 1PBQ_B.renumbered.pdb outputs/1PBQ_B.renumbered.pdb
diff 1Y1M_A.renumbered.pdb outputs/1Y1M_A.renumbered.pdb
diff 1Y1M_B.renumbered.pdb outputs/1Y1M_B.renumbered.pdb
diff 1Y1Z_A.renumbered.pdb outputs/1Y1Z_A.renumbered.pdb
diff 1Y20_A.renumbered.pdb outputs/1Y20_A.renumbered.pdb
diff 2A5S_A.renumbered.pdb outputs/2A5S_A.renumbered.pdb
diff 2A5T_A.renumbered.pdb outputs/2A5T_A.renumbered.pdb
diff 2A5T_B.renumbered.pdb outputs/2A5T_B.renumbered.pdb
diff 3JPW_A.renumbered.pdb outputs/3JPW_A.renumbered.pdb
diff 3JPY_A.renumbered.pdb outputs/3JPY_A.renumbered.pdb
diff 3OEK_A.renumbered.pdb outputs/3OEK_A.renumbered.pdb
diff 3OEL_A.renumbered.pdb outputs/3OEL_A.renumbered.pdb
diff 3OEM_A.renumbered.pdb outputs/3OEM_A.renumbered.pdb
diff 3OEN_A.renumbered.pdb outputs/3OEN_A.renumbered.pdb
diff 3QEL_A.renumbered.pdb outputs/3QEL_A.renumbered.pdb
diff 3QEL_B.renumbered.pdb outputs/3QEL_B.renumbered.pdb
diff 3QEL_C.renumbered.pdb outputs/3QEL_C.renumbered.pdb
diff 3QEL_D.renumbered.pdb outputs/3QEL_D.renumbered.pdb
diff 3QEM_A.renumbered.pdb outputs/3QEM_A.renumbered.pdb
diff 3QEM_B.renumbered.pdb outputs/3QEM_B.renumbered.pdb
diff 3QEM_C.renumbered.pdb outputs/3QEM_C.renumbered.pdb
diff 3QEM_D.renumbered.pdb outputs/3QEM_D.renumbered.pdb

../renumber_pdb.py -s L812M.B99990003_fit.pdb  -v "chain B","chain D" -r GluN2A_rattusN.fasta -o test.pdb
../renumber_pdb.py -s test.pdb -v "chain A","chain C" -r GluN1-4A_rattusN.fasta -o L812M.B99990003_fit.renumbered.pdb
diff -q L812M.B99990003_fit.renumbered.pdb outputs/L812M.B99990003_fit.renumbered.pdb
print -l

../renumber_pdb.py -s 3KG2.pdb -v "chain A","chain B","chain C","chain D" -r AMPA_rat_nosignal.fasta -o 3KG2.renumbered.ns.pdb
../renumber_pdb.py -s 3KG2.pdb -v "chain A","chain B","chain C","chain D" -r AMPA_rat_nosignal.fasta -o 3KG2.renumbered.ns.pdb
diff -q 3KG2.renumbered.ns.pdb outputs/3KG2.renumbered.ns.pdb
print -l

../renumber_pdb.py -s 1jen.pdb -v  "chain B","chain A" -r adometdc_hs.fasta -o 1JEN.renumbered.pdb --newres PYROC1
diff -q 1JEN.renumbered.pdb outputs/1JEN.renumbered.pdb
print -l

../renumber_pdb.py -s 1jen.pdb -v  "chain A B" -r adometdc_hs.fasta -o 1JEN.renumbered2.pdb --newres PYROC1
diff -q 1JEN.renumbered2.pdb outputs/1JEN.renumbered2.pdb
print -l

../renumber_pdb.py -s 2hfb.pdb -v "chain A","chain B" -r 2hfb.fasta -o 2hfb.renumbered.pdb
diff -q 2hfb.renumbered.pdb outputs/2hfb.renumbered.pdb
print -l
