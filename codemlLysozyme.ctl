      seqfile = lysozymeSmall.nuc   * sequence data file name
      outfile = mlc              * main result file name
     treefile = trees.lysozyme      * tree structure file name

        noisy = 9   * 0,1,2,3,9: how much rubbish on the screen
      verbose = 1   * 1: detailed output, 0: concise output
      runmode = 0   * 0: user tree;  1: semi-automatic;  2: automatic
                    * 3: StepwiseAddition; (4,5):PerturbationNNI 

      seqtype = 1   * 1:codons; 2:AAs; 3:codons-->AAs
    CodonFreq = 2   * 0:1/61 each, 1:F1X4, 2:F3X4, 3:codon table
        model = 2
                    * models for codons:
                        * 0:one, 1:b, 2:2 or more dN/dS ratios for branches

      NSsites = 0   * dN/dS among sites. 0:no variation, 1:neutral, 2:positive
        icode = 0   * 0:standard genetic code; 1:mammalian mt; 2-10:see below
        Mgene = 0   * 0:rates, 1:separate; 2:pi, 3:kappa, 4:all

    fix_kappa = 0   * 1: kappa fixed, 0: kappa to be estimated
        kappa = 2   * initial or fixed kappa
    fix_omega = 1   * 1: omega or omega_1 fixed, 0: estimate 
        omega = 1   * initial or fixed omega, for codons or codon-transltd AAs

    fix_alpha = 1   * 0: estimate gamma shape parameter; 1: fix it at alpha
        alpha = .0  * initial or fixed alpha, 0:infinity (constant rate)
       Malpha = 0   * different alphas for genes
        ncatG = 8   * # of categories in the dG or AdG models of rates

        clock = 0   * 0: no clock, unrooted tree, 1: clock, rooted tree
        getSE = 1   * 0: don't want them, 1: want S.E.s of estimates
 RateAncestor = 1   * (1/0): rates (alpha>0) or ancestral states (alpha=0)



* Specifications for duplicating table 1 of Yang (1998 MBE 15:568-573) 
*
*   table 1A: model = 0; fix_omega = 0;
*   table 1B: model = 2; fix_omega = 0;
*       When running the program, input the branch marks as:  
*        0 0 0 0 1 0 0 0 0 0 0
*   table 1C: model = 2; fix_omega = 0;
*       When running the program, input the branch marks as:  
*        1 0 0 0 0 0 0 0 0 0 0
*   table 1D: model = 2; fix_omega = 0;
*       When running the program, input the branch marks as:  
*        1 0 0 0 1 0 0 0 0 0 0
*   table 1E: model = 2; fix_omega = 0;
*       When running the program, input the branch marks as:  
*        1 0 0 0 2 0 0 0 0 0 0
*   table 1I: model = 2, fix_omega = 1; omega = 1;
*       When running the program, input the branch marks as:  
*        1 0 0 0 2 0 0 0 0 0 0
*   table 1J: model = 2; fix_omega = 1; omega = 1;
*       When running the program, input the branch marks as:  
*        2 0 0 0 1 0 0 0 0 0 0
*   Fig 2 (free ratio model, slow): model = 1; 
*
* When you choose model = 2, the program output the branches in the tree as 
* 8..9 9..1 9..2 8..10 10..11 11..3 11..4 10..5 8..12 12..6 12..7
* You need to draw the tree and identify the nodes and branches.  In this 
* case branch h is 8..9 (the first branch) and branch c is 10..11 (the fourth
* branch).  You can prepare the marks in a file and copy and paste them 
* when you run the program.  When fix_omega = 1, the last omega (identified 
* by the greatest mark) is fixed at the value given for omega.
*