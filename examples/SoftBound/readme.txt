MCMCTREE
Ziheng Yang

4 October 2005

This program now implements the method of Yang and Rannala (2005) for
calculating species divergence times on a given rooted tree,
accommodating gene sequences from multiple loci and accommodating
multiple fossil calibration points.  The current version enforces the
molecular clock, although work is under way to relax the clock
assumption.

The control file is mcmctree.ctl.  An example is included below, which
you can use to duplicate the results in table 3 of Yang and Rannala
(2005).


          seed = -1234567
       seqfile = mtCDNApri123.txt
      treefile = mtCDNApri.trees
       outfile = out

         ndata = 3
       usedata = 1    * 0: no data; 1:seq like; 2:normal approximation
         clock = 1    * 1: global clock; 2: independent rates; 3: correlated rates

         model = 0    * 0:JC69, 1:K80, 2:F81, 3:F84, 4:HKY85
         alpha = 0    * alpha for gamma rates at sites
         ncatG = 5    * No. categories in discrete gamma

     cleandata = 0    * remove sites with ambiguity data (1:yes, 0:no)?
 BlengthMethod = 0    * 0: arithmetic; 1: geometric; 2: Brownian

       BDparas = 2 2 .1   * birth, death, sampling
   kappa_gamma = 6 2      * gamma prior for kappa
   alpha_gamma = 1 1      * gamma prior for alpha

   rgene_gamma = 2 2     * gamma prior for rate for genes
  sigma2_gamma = 1 1    * gamma prior for sigma^2     (for clock=2 or 3)

      finetune = .2  0.5  0.5  0.1 .9  * times, rates, mixing, paras, RateParas

         print = 1
        burnin = 10000
      sampfreq = 10
       nsample = 10000


Here is a brief description of the variables.  I don't have time to
explain things about the MCMC algorithm.  You might read the document
for MCMCcoal, as some of the notes and warnings are applicable.

The tree and fossil calibration information are specified in the tree
file.

If seed is given a negative number, a random seed is used, determined
from the current time, which means that different runs using the same
control file will start from different places.  If you use a positive
number, that number will be used as the seed, so that running the same
program multiple times using the same control file will produce the
same results.

ndata is the number of loci (or site partitions) in a combined
analysis.  The program allows some species to be missing at some loci.
The example included partitions the three codon positions into 3
different genes.  In the combined analysis of multiple gene loci, the
same substitution model is used, but different parameters are assigned
and estimated for each partition.  You can look at the included
example sequence data file for the format.

Use clock = 1, which means a global clock.  

The 5 finetune parameters are like step lengths in the proposals.
When you run the program, it prints out four or five percentages on
the screen.  These are the acceptance rate for the proposals that
change the species divergence times, for changing the rates, for a
mixing step, for changing subsitution parameters (for models more
complex than JC69 with substitution parameters such as the
transition/transversion rate ratio kappa or gamma shape parameter
alpha), and for changing parameters in the rate-evolution model (only
used if clock = 2 or 3; do not use).  You should change the finetune
variables so that the corresponding acceptance proportions lie between
10% and 90%, or preferably between 15% to 80%.  If the acceptance rate
is too small, make the finetune steplength smaller.

      finetune = .2  0.5  0.5  0.1 .9  * times, rates, mixing, paras, RateParas

The variables rgene_gamma specifies the gamma prior for the rate for
each gene.  The gamma distribution is specified by two parameters a
and b, so that the distribution has mean a/b.  Ignore sigma2_gamma.


   rgene_gamma = 2 2     * gamma prior for rate for genes
  sigma2_gamma = 1 1    * gamma prior for sigma^2     (for clock=2 or 3)


Reference

Yang, Z., and B. Rannala. 2005. Bayesian estimation of species
divergence times under a molecular clock using multiple fossil
calibrations with soft bounds. Mol. Biol. Evol.: 23: 212-226.
