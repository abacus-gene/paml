          seed = -1234567

       seqfile = ../papers/DateSoftBound/mtCDNApri123.txt
      treefile = ../papers/DateSoftBound/mtCDNApri.trees

       outfile = out

         ndata = 3
       usedata = 1    * 0: no data; 1:seq like; 2:normal approximation
         clock = 1    * 1: global clock

         model = 0    * 0:JC69, 1:K80, 2:F81, 3:F84, 4:HKY85
         alpha = 0    * alpha for gamma rates at sites
         ncatG = 5    * No. categories in discrete gamma

     cleandata = 0    * remove sites with ambiguity data (1:yes, 0:no)?
 BlengthMethod = 0    * 0: arithmetic; 1: geometric; 2: Brownian

       BDparas = 2 2 .1   * birth, death, sampling
   kappa_gamma = 6 2      * gamma prior for kappa
   alpha_gamma = 1 1      * gamma prior for alpha

   rgene_gamma = 2 2     * gamma prior for rate for genes
  sigma2_gamma = 1 1    * gamma prior for sigma^2     (for clock=2 or 3)

      finetune = .02  0.5  0.5  0.1 .9  * times, rates, mixing, paras, RateParas

         print = 1
        burnin = 10000
      sampfreq = 10
       nsample = 100000
