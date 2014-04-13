      seqfile = mtYoder.aa * sequence data file name
      outfile = mlc       * main result file
     treefile = DatingAA.trees  * tree structure filename

        noisy = 9   * 0,1,2,3,9: how much rubbish on the screen
      verbose = 1   * 1: detailed output, 0: concise output

      seqtype = 2   * 1:codons; 2:AAs; 3:codons-->AAs
   aaRatefile = ../../mtREV24.dat
        model = 2
                    * models for AAs or codon-translated AAs:
                        * 0:poisson, 1:proportional, 2:Empirical, 3:Empirical+F
                        * 6:FromCodon, 7:AAClasses, 8:REVaa_0, 9:REVaa(nr=189)

    fix_alpha = 0   * 0: estimate gamma shape parameter; 1: fix it at alpha
        alpha = 0.2  * initial or fixed alpha, 0:infinity (constant rate)
        ncatG = 8   * # of categories in the dG or AdG models of rates

        clock = 2   * 0: no clock, unrooted tree, 1: clock, rooted tree
        getSE = 0   * 0: don't want them, 1: want S.E.s of estimates
 RateAncestor = 0   * (0,1,2): rates (alpha>0) or ancestral states (1 or 2)

   Small_Diff = 1e-6

* This is for table 2 Estimation of dates from amino acid sequences
* in Yoder and Yang (2000 MBE 17: 1081-1090).
* (1) At the prompt "node date", input 37 25  and then  -1 -1.
* (2) to run the global clock, chose clock = 1.
