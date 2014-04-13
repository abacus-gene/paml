README.txt
18 September 2011, Ziheng Yang

This example data file is for duplicating the analysis of Stadler and
Yang (2012).  The method dates divergence times when sequences have
sample dates, as in the case of viral sequences.  The data set was
aligned and analyzed by Lemey et al. (2003).


(A) Data format 

Have a look at the sequence alignment file HIV2ge.txt.  Note that the
end of the sequence name has the sample date.  For example, P03h1995
was isolated in 1995.


(B) ML analysis (Stadler and Yang 2012 table #)

To run the ML analysis, cd to the current folder and then run baseml.
This uses the control file baseml.ctl.

On Windows:

   cd paml4.6
   cd examples\TipDate
   ..\..\bin\baseml

On Unix/Linux/Mac OSX:

   cd paml4.6
   cd examples/TipDate
   ../../bin/baseml

The TipDate model is specified by the following line in baseml.ctl,
where the second number is the time unit (100 years in this case):

       TipDate = 1 100  * TipDate (1) & time unit


(C) Bayesian analysis (Stadler and Yang 2012 table #)

The control file is mcmctree.ctl is for running mcmctree.  Right now
the seqfile and treefile are specified as follows.  These are set up
so that you can run two copies of mcmctree in two different folders
(r1 and r2 inside the folder) at the same time.

        seqfile = ../HIV2ge.txt
       treefile = ../HIV2ge.tre

Start two command terminals, and cd to r1 and r2 respectively.  Then in each window type

   ..\..\..\bin\mcmctree ..\mcmctree.ctl (on Windows)

   ..\..\..\bin\mcmctree ..\mcmctree.ctl (on Unix/Linux/MacOSX)

You confirm that the two runs produce very similar results.  If not,
you need increase nsample or burnin etc.

Have a look at the control file mcmctree.ctl, the paml/mcmctree
documentation in the doc/ folder.

The following line in the control file specifies the TipDate model.
The second number is the time unit.

       TipDate = 1 100  * TipDate (1) & time unit

When you run the default analysis, you will see the following printout on the monitor.

TipDate model
Date range: (1995.00, 1982.00) => (0, 0.13). TimeUnit = 100.00.

The program scans the sequence names, and take the last field in the
sequence name as the sampling date.  The most recent sample date will
be time zero, and the other times will be rescaled using the time
unit.  For example the sequences in the example dataset are from 1982
to 1995, with 1995 becoming time 0, and 1982 becoming 0.13 as one time
unit is specified to be 100 years.

You should specify a prior on the age of the root of the tree.

       RootAge = B(0.5, 2.0, 0.01, 0.02)

The above means that the age of root is between 0.5 and 2.0 time units
(that is, between 1945 and 1795), but those bounds are soft, so that
the minimum- and maximum-age bounds are violated with probabilities 1%
and 2%, respectively.

Other things that are important include the prior on the substitution rate 

   rgene_gamma = 2 10   * gamma prior G(alpha, beta) for rate for genes

Note that the gamma distribution G(alpha, beta) has shape parameter
alpha and scale parameter beta.  You specify alpha depending on how
much confidence you have (with alpha = 1 or 2 to be diffuse priors and
alpha = 5 or 10 to be informative priors, say), and then specify the
scale parameter by having the mean alpha/beta in the right range for
the data.  Here alpha = 2 means a fairly diffuse prior, while
alpha/beta = 0.5 means 0.5 changes per site per time unit, or 0.005
changes per site per year (since the time unit is 100 years).

The clock model may also be important.

      clock = 1    * 1: global clock; 2: independent rates; 3: correlated rates

References

|P, Pybus OG, Wang B, Saksena NK, Salemi M, Vandamme AM. 2003. Tracing the origin and history of the HIV-2 epidemic. Proc Natl Acad Sci USA 100:6588-6592.

Stadler, T. and Z. Yang.  2012  Dating phylogenies with sequentially sampled tips.  Syst Biol, submitted


//end of file
