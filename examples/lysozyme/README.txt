(A) 
This folder contains the control file, the sequence data file and
the tree file for demonstrating codon models that use different dN/dS
ratios among lineages (Yang 1998).  The data set is the "small" data
set analyzed in Yang (1998).  The default control file let you
duplicate the results for the small data set in table 1 of Yang
(1998).  Also look at the tree file about specifying the branches of
interest, for which positive selection is tested.

To fix a particular w to 1, arrange the labels so the concerned branch
is the last and then use

       model = 2
   fix_omega = 1
       omega = 1

For example, the tree

 ((1,2) #2, ((3,4) #1, 5), (6,7) );     / * table 1E&J */

fits a model with w0 (background), w1, and w2.  Then the above
specification will force w2 = 1 to be fixed.

Usage:

	codeml lysozymeSmall.ctl

Or you can rename the file lysozyme.ctl as codeml.ctl, and then run 

         codeml


(B) The folder also contains another set of files for the "large" data
set analyzed in Yang (1998).  This data set is also used by Yang and
Nielsen (2002) to test the branch-site models, so I included the files
as well.  The control file lysozymeLarge.ctl specifies model A in that
paper, and should reproduce the results in table 2 of Yang and Nielsen
(2002).

   Model A:  model = 2    NSsites = 2
   Model B:  model = 2    NSsites = 3

Look at the tree file lysozymeLarge.trees for specification of the
"branch of interest" of "foreground" branch.  You can remove the first
line of numbers and the file will be readable from TreeView, which
allows you to show the branch (node) labels as well.

Note that the number of site classes under both models A and B is
fixed, and so the variable ncatG is ignored by codeml.  (To run the
discrete model with only 2 site classes, which is the null model to be
compared with model B, you should specify model = 0, NSsites = 3,
ncatG = 2.)  See the paper for details.  Also please heed the warnings
in the Discussion section of that paper.


References

Yang, Z. 1998. Likelihood ratio tests for detecting positive selection and application to primate lysozyme evolution. Mol. Biol. Evol. 15:568-573.

Yang, Z., and R. Nielsen, 2002 Codon-substitution models for detecting molecular adaptation at individual sites along specific lineages. Mol. Biol. Evol. 19: 908-917.

Ziheng Yang

11 September 2001, last modified on 4 February 2004
