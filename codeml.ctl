
      seqfile = stewart.aa    * sequence data file name
      outfile = mlc          * main result file name
     treefile = stewart.trees  * tree structure file name

        noisy = 3   * 0,1,2,3,9: how much rubbish on the screen
      verbose = 0   * 1: detailed output, 0: concise output
      runmode = 0   * 0: user tree;  1: semi-automatic;  2: automatic
                    * 3: StepwiseAddition; (4,5):PerturbationNNI; -2: pairwise

      seqtype = 2   * 1:codons; 2:AAs; 3:codons-->AAs
    CodonFreq = 2   * 0:1/61 each, 1:F1X4, 2:F3X4, 3:codon table
       aaDist = 0   * 0:equal, +:geometric; -:linear, {1-5:G1974,Miyata,c,p,v}
   aaRatefile = jones.dat * only used for aa seqs with model=empirical(_F)
                           * dayhoff.dat, jones.dat, mtmam.dat, or your own
        model = 2
                    * models for codons:
                        * 0:one, 1:b, 2:2 or more dN/dS ratios for branches
                    * models for AAs or codon-translated AAs:
                        * 0:poisson, 1:proportional, 2:Empirical, 3:Empirical+F
                        * 6:FromCodon, 7:AAClasses, 8:REVaa_0, 9:REVaa(nr=189)

      NSsites = 0   * 0:one w;1:neutral;2:selection; 3:discrete;4:freqs;
                    * 5:gamma;6:2gamma;7:beta;8:beta&w;9:beta&gamma;
                    * 10:beta&gamma+1; 11:beta&normal>1; 12:0&2normal>1;
                    * 13:3normal>0
        icode = 0  * 0:universal code; 1:mammalian mt; 2-10:see below
        Mgene = 0   * 0:rates, 1:separate; 2:pi, 3:kappa, 4:all

    fix_kappa = 0   * 1: kappa fixed, 0: kappa to be estimated
        kappa = 2.3   * initial or fixed kappa
    fix_omega = 0   * 1: omega or omega_1 fixed, 0: estimate 
        omega = 1  * initial or fixed omega, for codons or codon-based AAs

    fix_alpha = 1   * 0: estimate gamma shape parameter; 1: fix it at alpha
        alpha = 0.  * initial or fixed alpha, 0:infinity (constant rate)
       Malpha = 0   * different alphas for genes
        ncatG = 8   * # of categories in the dG or AdG models of rates

      fix_rho = 1   * 0: estimate rho (correlation parameter); 1: fix it at rho
          rho = 0.  * initial or fixed rho,   0:no correlation

        clock = 0   * 0: no clock, unrooted tree, 1: clock, rooted tree
        getSE = 0   * 0: don't want them, 1: want S.E.s of estimates
 RateAncestor = 1   * (0,1,2): rates (alpha>0) or ancestral states (1 or 2)

* Genetic codes: 0:universal, 1:mammalian mt., 2:yeast mt., 3:mold mt.,
* 4: invertebrate mt., 5: ciliate nuclear, 6: echinoderm mt., 
* 7: euplotid mt., 8: alternative yeast nu. 9: ascidian mt., 
* 10: blepharisma nu.
* These codes correspond to transl_table 1 to 11 of GENEBANK.

*   Small_Diff = 8e-7
*        ndata = 5
*    cleandata = 1  * remove sites with ambiguity data (1:yes, 0:no)?