          seed = 1357
      seqfile = H1.txt
*      seqfile = ../H1.CP123.txt
     treefile = H1.tre

       outfile = out
       seqtype = 0  * 0: nucleotides; 1:codons; 2:AAs
         noisy = 3

       usedata = 2 in.BV.HKYG5  * 0: no data; 1:seq like; 2:use in.BV; 3: out.BV

         ndata = 1   * 
         clock = 1   * 1: global clock; 2: independent rates; 3: correlated rates
       TipDate = 1 100  * TipDate (1) & time unit

*   fossilerror = 0    * no fossil errors
*        RootAge = B(1.1, 2.1, .025, .025)  * used if no fossil for root
        RootAge = B(1, 5, .001, .001)  * used if no fossil for root
*        RootAge = B(4, 6, .001, .001)  * used if no fossil for root

         model = 4    * 0:JC69, 1:K80, 2:F81, 3:F84, 4:HKY85
    aaRatefile = wag.dat
         alpha = 0.5    * alpha for gamma rates at sites
         ncatG = 5    * No. categories in discrete gamma

     cleandata = 0    * remove sites with ambiguity data (1:yes, 0:no)?

*       BDparas = 10 5 .001  * lambda, mu, rho, psi for birth-death-sampling model
*       BDparas = 2   1 0 1.8  * lambda, mu, rho, psi for birth-death-sampling model
*       BDparas = 4   2 0 3.6  * lambda, mu, rho, psi for birth-death-sampling model
*       BDparas = 1 0.5 0 0.9  * lambda, mu, rho, psi for birth-death-sampling model
*       BDparas = 20 10 0 18  * lambda, mu, rho, psi for birth-death-sampling model
*       BDparas = 0.2 0.1 0 0.18  * lambda, mu, rho, psi for birth-death-sampling model
       BDparas = 2   1 0 1.8  * lambda, mu, rho, psi for birth-death-sampling model

   kappa_gamma = 2 1   * gamma prior for kappa
   alpha_gamma = 2 4   * gamma prior for alpha

   rgene_gamma = 2 10   * gamma prior for rate for genes
  sigma2_gamma = 1 20    * gamma prior for sigma^2  (for clock=2)

      finetune = 1: 0.2 0.04 .05 .05 .05 .05  * auto (0 or 1) : times, rates, mixing, paras, RateParas, FossilErr

         print = 1
        burnin = 10000
      sampfreq = 5
       nsample = 20000
