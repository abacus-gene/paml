
    seqfile = mtprim9.nuc  * sequence data file name
    outfile = mcmctree.out * main result file name
   treefile = trees.9s     * tree file for initial rooted tree topology
     LHfile = LHs  * LH file name, read (MCMC=0) or rewritten (MCMC=1)
       MCMC = 1     * 0: read LHs from LHfile, 1: use MCMC to generate LHs
       beta = 0.15  * prob{labeled history change}, used only if MCMC=1
       seed = 1234567 * random number seed

     delta0 = .5     * small number for MCMC, used if MCMC=1 only
     delta1 = .5     * smaller number for comparing candidate LHs

      model = 3      * 0:JC69, 1:K80, 2:F81, 3:F84, 4:HKY85
      kappa = 1.63 2.127  * kappa in K80, F84, or HKY85
      alpha = 0      * alpha for gamma rates at sites
      ncatG = 4      * No. categories in discrete gamma

   hierarch = 0      * 1:hierarchical; 0:empirical Bayes analysis
      birth = 6.7    * birth rate
      death = 2.5    * death rate
     sample = .06   * sampling fraction
     mutate = .24 * mutation rate (# of mutations from root to present)
