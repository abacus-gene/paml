
      seqfile = stewart.aa   * sequence data file name
      outfile = mlc          * main result file name

        noisy = 3   * 0,1,2,3,9: how much rubbish on the screen
      runmode = 0   * 0: user tree;  1: semi-automatic;  2: automatic

      seqtype = 2   * 1:codons; 2:AAs; 3:codons-->AAs
        model = 3
                    * models for AAs or codon-translated AAs:
                    *  0:poisson, 1:equal_input, 2:Dayhoff(KYM), 3:Jones(KYM)
                    *  4:Dayhoff_Pi, 5:Jones_Pi, 6:FromCodon, 7:NGrantham
                    *  8:REVaa_0, 9:REVaa(nr=189)
                    * models for codons:
                    *  0:1/61 each, 1:1X4, 2:3X4, 3:codon table, 7:NGrantham

        icode = 0   * 0:standard genetic code; 1:mammalian mt; 2-10:below
        Mgene = 0   * 0:rates, 1:separate; 2:pi, 3:kapa, 4:all

   given_kapa = 0   * 1: kapa given, 2: omega given, 3: both given, 0: estimate
         kapa = 2.  * initial or given kappa
        omega = 1.4  * initial or given omega, for codons or codon-transltd AAs
		 
   given_alfa = 1
         alfa = 0.  * initial or given alpha, 0:infinity (constant rate)
        Malfa = 0   * different alphas for genes
        ncatG = 8   * # of categories in the dG or AdG models of rates

    given_rho = 1 
          rho = 0.  * initial or given rho,   0:no correlation

        clock = 0   * 0: no clock, unrooted tree, 1: clock, rooted tree
        getSE = 0   * 0: don't want them, 1: want S.E.s of estimates
 RateAncestor = 1    * (1/0): rates (alpha>0) or ancestral states (alpha=0)


* Genetic codes: 0:standard, 1: mammalian mt., 2: yeast mt., 3: mold mt.,
* 4: invertebrate mt., 5: ciliate nuclear, 6: echinoderm mt., 
* 7: euplotid mt., 8: alternative yeast nu. 9: ascidian mt., 
* 10: blepharisma nu.   
* These codes correspond to transl_table 1 to 11 of GENEBANK.
