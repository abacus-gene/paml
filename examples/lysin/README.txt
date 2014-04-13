Notes by Ziheng Yang
29 November 2001

   README.txt
   lysin.trees (tree file)
   lysin.nuc  (sequence data file)
   codeml.ctl  (control file)
   lysinResult.txt (results under M0 and M3)
   lysinPosteriorP.txt (posterior probabilities under M3)

This folder contains the control file, the sequence data file and the
tree file for demonstrating codon models that assign different dN/dS
ratios among sites in the sequence (Nielsen & Yang 1998; Yang, Nielsen, 
Goldman & Pedersen 2000).  The included data set is the sperm lysin genes 
from 25 abalone species used in Yang, Swanson & Vacquier (2000).  The default 
control file (with NSsites = 3) lets you duplicate the results in table 1 of 
that paper.  To run the program, try 

	codeml

You can also run several models in one batch by 

	codemlNSsites


The file lysinPosteriorP.txt includes part of the output from the file
rst for model M3 (NSsites=3).  The first 3 columns are the three
probabilities for the three site classes; you can use them to make
figure 1 of Yang, Swanson & Vacquier (2000).  In parentheses are the
most likely class numbers.  The last two columns are the posterior
average w for the site and the probability for the most likely class
(redundant).

Sorry that this example calculation actually takes a long time.  The
M3 model takes 3 hours on my portable (1.2GHz).  I don't have a
smaller data set at hand.

References

Nielsen, R., and Z. Yang. 1998. Likelihood models for detecting positively selected amino acid sites and applications to the HIV-1 envelope gene. Genetics 148:929-936.
Yang, Z., R. Nielsen, N. Goldman and A.-M. K. Pedersen. 2000. Codon-substitution models for heterogeneous selection pressure at amino acid sites. Genetics 155:431-449.
Yang, Z., W. J. Swanson and V. D. Vacquier. 2000. Maximum likelihood analysis of molecular adaptation in abalone sperm lysin reveals variable selective pressures among lineages and sites. Mol. Biol. Evol. 17:1446-1455.
Yang, Z., and J. P. Bielawski. 2000. Statistical methods for detecting molecular adaptation. Trends Ecol. Evol. 15:496-503.
