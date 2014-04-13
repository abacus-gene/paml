          seed = -1357

        seqfile = HIV2ge.txt
       treefile = HIV2ge.tre
        outfile = out

       seqtype = 0  * 0: nucleotides; 1:codons; 2:AAs
       usedata = 0 in.BV.HKYG5  * 0: no data; 1:seq like; 2:use in.BV; 3: out.BV

         ndata = 1   * 3 for mt2G3P and 6 for 4G6P
         clock = 1    * 1: global clock; 2: independent rates; 3: correlated rates
       TipDate = 1 100  * TipDate (1) & time unit

       RootAge = B(0.5, 2.0, 0.01, 0.01)  * used if no fossil for root
*       RootAge = B(2.9, 3.1, 0.01, 0.01)  * used if no fossil for root
*       RootAge = B(0.2, 1.5, 0.01, 0.01)  * used if no fossil for root
*       RootAge = B(0.8, 2.5, 0.01, 0.01)  * used if no fossil for root

         model = 4    * 0:JC69, 1:K80, 2:F81, 3:F84, 4:HKY85
    aaRatefile = wag.dat
         alpha = 0.5    * alpha for gamma rates at sites
         ncatG = 5    * No. categories in discrete gamma

     cleandata = 0    * remove sites with ambiguity data (1:yes, 0:no)?

       BDparas = 2 1 0 1.8  * lambda, mu, rho, psi for birth-death-sampling model
   kappa_gamma = 6 2      * gamma prior for kappa
   alpha_gamma = 1 1      * gamma prior for alpha

   rgene_gamma = 2 10 2   * gammaDir prior for rate for genes
  sigma2_gamma = 2 2 1   * gammaDir prior for sigma^2     (for clock=2 or 3)
      finetune = 1: 0.5 0.2 .1 .05 .05 .05 * auto (0 or 1): times, musigma2, rates, mixing, paras, FossilErr

         print = 1
        burnin = 20000
      sampfreq = 2
       nsample = 100000
