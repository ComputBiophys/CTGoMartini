martinize2 -f 7LJB_clean.pdb -o system.top -x protein_cg.pdb -p backbone -ff martini3001 -cys auto -dssp -go -go-low 0.3 -go-up 1.1 -go-eps 12.0 -name protein > tmp_7LJB
martinize2 -f 4WFE_clean.pdb -o system.top -x protein_cg.pdb -p backbone -ff martini3001 -cys auto -dssp -go -go-low 0.3 -go-up 1.1 -go-eps 12.0 -name protein > tmp_4WFE
