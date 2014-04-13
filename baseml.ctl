      seqfile = brown.nuc * sequence data file name
      outfile = mlb       * main result file
     treefile = brown.trees  * tree structure filename

        noisy = 3   * 0,1,2,3: how much rubbish on the screen
      verbose = 0   * 1: detailed output, 0: concise output
      runmode = 0   * 0: user tree;  1: semi-automatic;  2: automatic
                    * 3: StepwiseAddition; (4,5):PerturbationNNI 

        model = 4   * 0:JC69, 1:K80, 2:F81, 3:F84, 4:HKY85, 5:TN93, 6:REV
        Mgene = 0   * 0:rates, 1:separate; 2:diff pi, 3:diff kapa, 4:all diff

    fix_kappa = 0   * 0: estimate kappa; 1: fix kappa at value below
        kappa = 5  * initial or fixed kappa

    fix_alpha = 0   * 0: estimate alpha; 1: fix alpha at value below
        alpha = 0.3  * initial or fixed alpha, 0:infinity (constant rate)
       Malpha = 0   * 1: different alpha's for genes, 0: one alpha
        ncatG = 8   * # of categories in the dG, AdG, or nparK models of rates

      fix_rho = 1   * 0: estimate rho; 1: fix rho at value below 
          rho = 0.  * initial or fixed rho,   0:no correlation
        nparK = 0   * rate-class models. 1:rK, 2:rK&fK, 3:rK&MK(1/K), 4:rK&MK 

        clock = 0   * 0: no clock, unrooted tree, 1: clock, rooted tree
        nhomo = 0   * 0 & 1: homogeneous, 2: kappa for branches, 3: N1, 4: N2
        getSE = 0   * 0: don't want them, 1: want S.E.s of estimates
 RateAncestor = 1   * (0,1,2): rates (alpha>0) or ancestral states

*   Small_Diff = 4e-7
*    cleandata = 0  * remove sites with ambiguity data (1:yes, 0:no)?
*        ndata = 1

*        icode = 0  * (with RateAncestor=1. try "GC" in data,model=4,Mgene=4)
