# Undergraduate-graduation-project
Undergraduate graduation project in Soochow University

Python and R are used in data processing.

Description of the py files and R files:
3.py: Python code for data processing.
4.py: Python code for the filtering process of the Reviewed mutation dataset.
5.py: Python code for the filtering process of the Unreviewed mutation dataset.

conservation.py and tongji_seq_len.py: Calculated and analyzed the conservation (using Provean) and residue position (using DSSP and ScanNet) of proteins affected by mutations in the IntAct mutation dataset.
Density distribution and bar graphs were generated.

map.py: Using the Best structures and SIFTS mapping api from PDBe to find protein complex structures for the IntAct mutation dataset. And using PremPS and Mutabind2 to calculate.
map2.py: Using the Best structures and SIFTS mapping api from PDBe to find protein complex structures for the S4191 mutation dataset. And using PremPS and Mutabind2 to calculate.

chayi.R: R Codes for conducting t-tests, Wilcoxon tests, Fisher tests, and chi-square tests .

tongji.py: Python codes to calculate the numbers of 4 mutations and draw .
