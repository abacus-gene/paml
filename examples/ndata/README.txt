README.txt
Ziheng, 17 March 2022


(a) ndata = 3 

(b) ndata = 3 maintree    # this uses the maintree file to generate the subtrees and runs the ML analysis
(b) ndata = 3 maintree 1  # this uses the maintree file to generate the subtrees and runs the ML analysis

(c) ndata = 3 maintree 0  # this generates the genetrees.trees file but does not run ML

This note explains the three options for analyzying multiple datasets.
The multiple alignments are in one sequence file, one after another.
In the example data file 4s-ndata.txt, there are 3 alignments, with 2,
3, and 4 sequences.

In case (a), there will be 3 trees (or tree blocks) for the three
alignments.  The tree file 4s-ndata.trees has the following:

2 1
(S2,S4);
3 1
(S2,S3,S4);
4 1
((S1,S2),S3,S4);

In case (b), the tree file 4s-ndata-maintree.trees includes only one
main tree, which has all species names that occur in all the
alignments.  The "maintree" keyword is used to instruct
codeml/baseml to generate the gene tree for each alignment by pruning
off species that are missing in the alignment.  The program then
conducts the ML analysis.

Case (c) is the same as (b), but the program prints out the subtrees
for the alignments into a file, called genetrees.trees, without doing
the ML iterations.

Options (b) and (c) work with model M0 in codeml and the NSsites
models (you can specify multiple NSsites models in one run).  For
branch and branch-site models, the subtree pruning may cause multiple
branches on the main tree to be merged into one.  If those merged
branches have the same label, the label is retained.  Otherwise the
labels are deleted.  This may cause problems for the analysis as the
resulting subtree may not have any labels.  We have not decided what
to do in this case.
