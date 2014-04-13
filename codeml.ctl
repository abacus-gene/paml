
      seqfile = stewart.aa  * sequence data file name
      outfile = mlc          * main result file name
     treefile = trees.6s    * tree structure file name

        noisy = 3   * 0,1,2,3,9: how much rubbish on the screen
      verbose = 0   * 1: detailed output, 0: concise output
      runmode = 0   * 0: user tree;  1: semi-automatic;  2: automatic
                    * 3: StepwiseAddition; (4,5):PerturbationNNI 

      seqtype = 2   * 1:codons; 2:AAs; 3:codons-->AAs
    CodonFreq = 2   * 0:1/61 each, 1:F1X4, 2:F3X4, 3:codon table
       aaDist = 0   * 0:equal aa distance, 1: Grantham1974 (for codon models)
   aaRatefile = jones.dat  * only used for aa seqs and when model=empirical(_F)
                           * dayhoff.dat, jones.dat, mtmam.dat, or your own
        model = 2
                    * models for codons:
                        *  0:one N/S rate, 1: b ratios, 2: 2 or more ratios
                    * models for AAs or codon-translated AAs:
                        *  0:poisson, 1:equal_input, 2:Empirical, 3:Empirical+F
                        *  6:FromCodon, 7:AAClasses, 8:REVaa_0, 9:REVaa(nr=189)

      NSsites = 0   * dN/dS among sites. 0:no variation, 1:neutral, 2:positive
        icode = 0   * 0:standard genetic code; 1:mammalian mt; 2-10:below
        Mgene = 0   * 0:rates, 1:separate; 2:pi, 3:kappa, 4:all

    fix_kappa = 0   * 1: kappa fixed, 0: kappa to be estimated
        kappa = 4.6 * initial or fixed kappa
    fix_omega = 0   * 1: omega or omega_1 fixed, 0: estimate 
        omega = 1   * initial or fixed omega, for codons or codon-transltd AAs
		 
    fix_alpha = 1   * 0: estimate gamma shape parameter; 1: fix it at alpha
        alpha = 0.  * initial or fixed alpha, 0:infinity (constant rate)
       Malpha = 0   * different alphas for genes
        ncatG = 8   * # of categories in the dG or AdG models of rates

      fix_rho = 1   * 0: estimate rho (correlation parameter); 1: fix it at rho
          rho = 0.  * initial or fixed rho,   0:no correlation

        clock = 0   * 0: no clock, unrooted tree, 1: clock, rooted tree
        getSE = 0   * 0: don't want them, 1: want S.E.s of estimates
 RateAncestor = 1   * (1/0): rates (alpha>0) or ancestral states (alpha=0)

* Genetic codes: 0:universal, 1: mammalian mt., 2: yeast mt., 3: mold mt.,
* 4: invertebrate mt., 5: ciliate nuclear, 6: echinoderm mt., 
* 7: euplotid mt., 8: alternative yeast nu. 9: ascidian mt., 
* 10: blepharisma nu.
* These codes correspond to transl_table 1 to 11 of GENEBANK.

