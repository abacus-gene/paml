          seed = -1
       seqfile = combined.phy
      treefile = species.trees
      mcmcfile = mcmc.txt
       outfile = out.txt

         ndata = 1
       seqtype = 2
       usedata = 2 in.BV 1
         clock = 32 * 3 rate_martingale_AR, 31 lograte_martingale_AR, 32 full GBM
*      RootAge = 'B(3.20,4.51,0.025,0.01)'  * safe constraint on root age, used if no fossil for root.

         model = 0 * 0 : JC69, 1: K80, 2: F81, 3: F84, 4: HKY85
         alpha = 0.5
         ncatG = 4
   duplication = 0
     cleandata = 0    * remove sites With ambiguity data (1:yes, 0:no)?

       BDparas = 0.4 0.2 0.1 M
   kappa_gamma = 6 2       * gamma prior for kappa
   alpha_gamma = 1 1       * gamma prior for alpha

   rgene_gamma = 1 50 1
  sigma2_gamma = 1 10 1
  drift_lnorm  = 0 0.1 1

      finetune = 1: .1 .1 .1 .1 .1 .1 * auto (0 or 1): times, musigma2, rates, mixing, paras, FossilErr

         print = 1
        burnin = 1000
      sampfreq = 50
       nsample = 1000

