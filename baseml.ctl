
      seqfile = brown.nuc * sequence data file name
      outfile = mlb       * main result file
     treefile = trees.5s  * tree structure file name

        noisy = 3    * 0,1,2,3: how much rubbish on the screen
      verbose = 1    * 1: detailed output, 0: concise output
      runmode = 0    * 0: user tree;  1: semi-automatic;  2: automatic

        model = 0    * 0:JC69, 1:K80, 2:F81, 3:F84, 4:HKY85, 5:TN93, 6:REV
        Mgene = 0    * 0:rates, 1:separate; 2:diff pi, 3:diff kapa, 4:all diff
 		 
  given_kappa = 0
        kappa = 6    * initial or given kappa
 		 
  given_alpha = 0
        alpha = 0.2  * initial or given alpha, 0:infinity (constant rate)
       Malpha = 0    * 1: different alpha's for genes, 0: one alpha
        ncatG = 8    * # of categories in the dG, AdG, or nparK models of rates

    given_rho = 1
          rho = 0.   * initial or given rho,   0:no correlation
        nparK = 0    * rate-class models. 1:rK, 2:rK&fK, 3:rK&MK(1/K), 4:rK&MK 

        clock = 0    * 0: no clock, unrooted tree, 1: clock, rooted tree
        nhomo = 0    * 0 & 1: homogeneous, 2: N1, 3: N2
        getSE = 0    * 0: don't want them, 1: want S.E.s of estimates
 RateAncestor = 0    * (1/0): rates (alpha>0) or ancestral states (alpha=0)
