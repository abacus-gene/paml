
       seqfile = HIV2ge.txt
      treefile = HIV2ge.tre

      outfile = mlb       * main result file
        noisy = 3   * 0,1,2,3: how much rubbish on the screen
      verbose = 1   * 1: detailed output, 0: concise output
      runmode = 0   * 0:user tree;  1:semi-automatic;  2:automatic
                    * 3:StepwiseAddition; (4,5):PerturbationNNI 

        model = 4   * 0:JC69, 1:K80, 2:F81, 3:F84, 4:HKY85, 5:T92, 6:TN93, 7:REV
                    * 8:UNREST, 9:REVu; 10:UNRESTu

*       model = 10  [0]  /* JC69 */
*       model = 10  [1 (TC CT AG GA)]  /* K80 */
*       model = 10  [11 (TA) (TG) (CT) (CA) (CG) (AT) (AC) (AG) (GT) (GC) (GA) ]  /* UNREST */
*       model = 10  [5 (AC CA) (AG GA) (AT TA) (CG GC) (CT TC)]  /* SYM */
*       model = 10  [5 (CT GA) (TC AG) (AT TA) (AC TG) (GT CA) ]
*       model = 9   [2 (TA TG CA CG) (AG)]    /* TN93 */

        Mgene = 0   * 0:rates, 1:separate; 2:diff pi, 3:diff kappa, 4:all diff

      TipDate = 1 100
        clock = 1   * 0:no clock, 1:clock; 2:local clock; 3:CombinedAnalysis
    fix_kappa = 0   * 0: estimate kappa; 1: fix kappa at value below
        kappa = 2   * initial or fixed kappas
 
    fix_alpha = 0   * 0: estimate alpha; 1: fix alpha at value below
        alpha = 0.5   * initial or fixed alpha, 0:infinity (constant rate)
        ncatG = 5   * # of categories in the dG, AdG, or nparK models of rates
      fix_rho = 1   * 0: estimate rho; 1: fix rho at value below
          rho = 0   * initial or fixed rho, 0:no correlation
       Malpha = 0   * 1: different alpha's for genes, 0: one alpha
        nparK = 0   * rate-class models. 1:rK, 2:rK&fK, 3:rK&MK(1/K), 4:rK&MK 

        getSE = 1   * 0: don't want SEs of estimates, 1: want SEs
 RateAncestor = 0   * (0,1,2): rates (alpha>0) or ancestral states
       method = 0 * Optimization method 0: simultaneous; 1: one branch a time

   Small_Diff = 0.5e-6
    cleandata = 0  * remove sites with ambiguity data (1:yes, 0:no)?
  fix_blength = 0 * 0: ignore, -1: random, 1: initial, 2: fixed
