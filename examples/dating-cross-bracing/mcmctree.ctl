          seed = -1
       seqfile = mc.txt
      treefile = 6S.trees
      mcmcfile = mcmc.txt
       outfile = out.txt

         ndata = 1
       seqtype = 0 * 0 : nucleotides; 1:codons; 2:AAs
       usedata = 1 * 0: no data; 1:seq Like; 2:normal approximation; 3:out.BV(in.BV)
         clock = 2 * 1: Global clock; 2: independent rates; 3: correlated rates
*      RootAge = 'B(3.20,4.51,0.025,0.01)'  * safe constraint on root age, used if no fossil for root.

         model = 0 * 0 : JC69, 1: K80, 2: F81, 3: F84, 4: HKY85
         alpha = 0    * alpha for gamma rates at sites
         ncatG = 5    * No. categories In discrete gamma
*   duplication = 1
     cleandata = 0    * remove sites With ambiguity data (1:yes, 0:no)?

       BDparas = 1 1 0.1   * birth, death, sampling
   kappa_gamma = 6 2       * gamma prior for kappa
   alpha_gamma = 1 1       * gamma prior for alpha

   rgene_gamma = 2 20 1    * gammaDir prior for rate for genes
  sigma2_gamma = 1 10 1    * gammaDir prior for sigma^2     (for clock=2 or 3)

      finetune = 1: .1 .1 .1 .1 .1 .1 * auto (0 or 1): times, musigma2, rates, mixing, paras, FossilErr

         print = 1   * 0: no mcmc sample; 1: everything except branch rates 2: everything
        burnin = 20000
      sampfreq = 100
       nsample = 2000
