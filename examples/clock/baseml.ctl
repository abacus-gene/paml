      seqfile = mtYoder3.nuc * sequence data file name
     treefile = DatingNUC.trees  * tree structure filename

      outfile = mlb       * main result file
        noisy = 9  * 0,1,2,3: how much rubbish on the screen
      verbose = 1  * 1: detailed output, 0: concise output
      runmode = 0  * 0: user tree;  1: semi-automatic;  2: automatic
                   * 3: StepwiseAddition; (4,5):PerturbationNNI 

        model = 4  * 0:C69, 1:K80, 2:F81, 3:F84, 4:HKY85, 5:TN93, 6:REV

    fix_kappa = 0  * 0: estimate kappa; 1: fix kappa at value below
        kappa = 5  * initial or fixed kappa

    fix_alpha = 1  * 0: estimate alpha; 1: fix alpha at value below
        alpha = 0. * initial or fixed alpha, 0:infinity (constant rate)
        ncatG = 8  * # of categories in the dG, AdG, or nparK models of rates

        clock = 2  * 0:no clock, 1:clock; 2:local clock; 3:TipDate
        getSE = 0  * 0: don't want them, 1: want S.E.s of estimates
 RateAncestor = 0  * (0,1,2): rates (alpha>0) or ancestral states

   Small_Diff = 4e-7

* This is for table 3 "one rate" in Yoder and Yang (2000 MBE 17: 1081-1090).
* Estimation of dates from third codon positions
* (1) At the prompt "node date", input 37 25  and then  -1 -1.
* (2) to run the global clock, chose clock = 1.
