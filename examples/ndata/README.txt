README.txt
Ziheng, 11 November 2022

This note explains the four options for analyzying multiple datasets implemented in 
baseml and codeml.  The multiple alignments are in one sequence data file, one after 
another.  

[Note for myself: ndata_trees_opt = 0 1 2 3 are for cases a b c d.]

(a) ndata = 3
This is the old option, with one tree block used for all datasets/alignments.
All alignments in the file must have the same set of sequence names.  A common 
use of this option is to simulate replicate datasets (using evolver, say) 
and then analyze them one by one using the same set of trees.  
    cd ndata
	..\..\bin\evolver 5 MCbase.txt
    ..\..\bin\baseml baseml.ctl

The rest of the note is about the new options b, c, d, using the example data files 
The example data file 4s-ndata.txt include 3 alignments, with 2, 3, and 4 sequences, 
respectively.  The tree files are 4s.trees, etc.  

(b) ndata = 3 separate_trees  # each alignment has its own tree block
    baseml baseml-ndata.ctl
    codeml codeml-ndata.ctl
	
(c) ndata = 3 maintree    # this uses the maintree file to generate the subtrees and runs the ML analysis
(c) ndata = 3 maintree 1  # this uses the maintree file to generate the subtrees and runs the ML analysis
    baseml baseml-ndata-maintree-ml.ctl
    codeml codeml-ndata-maintree-ml.ctl

(d) ndata = 3 maintree 0  # this generates the genetrees.trees file but does not run ML
    baseml baseml-ndata-maintree.ctl
    codeml codeml-ndata-maintree.ctl

In case (b), there will be 3 trees (or tree blocks) for the three alignments.  
The tree file 4s-ndata.trees has the following:

2 1
(S2,S4);
3 1
(S2,S3,S4);
4 1
((S1,S2),S3,S4);

In case (c), the tree file 4s-ndata-maintree.trees includes only one
main tree, which has all species names that occur in all the
alignments.  The "maintree" keyword is used to instruct
codeml/baseml to generate the gene tree for each alignment by pruning
off species that are missing in the alignment.  The program then
conducts the ML analysis.

Case (d) is the same as (c), but the program prints out the subtrees
for the alignments into a file, called genetrees.trees, without doing
the ML iterations.

Options (c) and (d) work with model M0 in codeml and the NSsites
models (you can specify multiple NSsites models in one run).  For
branch and branch-site models, the subtree pruning may cause multiple
branches on the main tree to be merged into one.  If those merged
branches have the same label, the label is retained.  Otherwise the
labels are deleted.  This may cause problems for the analysis as the
resulting subtree may not have any labels.  We have not decided what
to do in this case.
