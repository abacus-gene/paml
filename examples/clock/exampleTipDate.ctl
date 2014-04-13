      seqfile = exampleTipDate.phy
      outfile = mlb
     treefile = exampleTipDate.trees

        noisy = 3   * 0,1,2,3: how much rubbish on the screen
      verbose = 0   * 1: detailed output, 0: concise output

        model = 4  * 0:JC69, 1:K80, 2:F81, 3:F84, 4:HKY85, 5:TN93, 6:REV
    fix_kappa = 0  * 0: estimate kappa; 1: fix kappa at value below
        kappa = 2.4   * initial or fixed kappa

    fix_alpha = 1  * 0: estimate alpha; 1: fix alpha at value below
        alpha = 0  * initial or fixed alpha, 0:infinity (constant rate)
        ncatG = 8  * # of categories in the dG, AdG, or nparK models of rates

        clock = 3   * 0:no clock, 1:global clock; 2:local clock; 3:TipDate
        getSE = 1   * 0: don't want them, 1: want S.E.s of estimates
 RateAncestor = 0   * (0,1,2): rates (alpha>0) or ancestral states

   Small_Diff = .1e-6
