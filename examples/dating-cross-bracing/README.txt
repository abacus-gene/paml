Notes by Sandra Álvarez-Carretero
Last modified: 20 January 2026

Variable "duplication" is a switch ("0" for no and "1" for yes)
to turn on the option for duplicated genes/proteins.
There are two types of constraints that can be enabled, equality
and inequality constraints, which will be triggered based on the
notation used in the tree file to identify specific nodes.

For instance, when using paralogous genes, there may be some nodes
in the tree that correspond to the same speciation event, which we
may want to fix to the same age. We can identify such nodes in the
tree file by using node labels such as #1, #2, etc.

The following Newick trees show how this notation would be
incorporated to specify this type of node age constraints in
MCMCtree:


(((A1, A2) [#1 B{0.2, 0.4}], A3) #2, ((B1, B2) #1, B3) #2 ) 'B(0.9,1.1)';
(((A1, A2) [#1 B{0.2, 0.4}], A3) #2, ((B1, B2) [#1], B3) [#2] ) 'B(0.9,1.1)';
(((A1, A2) [#1 B{0.2, 0.4}], A3) #2, ((B1, B2) [#1 B{0.2, 0.4}], B3) #2 ) 'B(0.9,1.1)';

(((A1, A2) [#1 >0.2 <0.4], A3) #2, ((B1, B2) #1, B3) #2 ) >0.9<1.1;
(((A1, A2) [#1 >0.2 <0.4], A3) #2, ((B1, B2) [#1], B3) [#2] ) >0.9<1.1;
(((A1, A2) [#1 >0.2 <0.4], A3) #2, ((B1, B2) [#1 >0.2 <0.4], B3) #2 ) >0.9<1.1;


All the trees above are equivalent ways of specifying equality
constraints (also known as cross-bracing) in MCMCtree. Among
the nodes on the tree with the same label, one is chosen as
the "driver" node while the others are labelled as "mirror"
nodes. If calibration information is provided on one of the
shared nodes, the same information will therefore apply to all
shared nodes. If calibration information is provided on multiple
shared nodes, that information must be the same.
The time prior (or the prior on all node ages on the tree) is
constructed by using a density at the root of the tree
(specified by the user via the input tree file), while the
ages of all non-calibration nodes are given by the uniform
density. This time prior is similar to that used by
Thorne et al. (1998). The parameters in the birth-death process
with species sampling (i.e., lambda, mu, rho; which are specified
using variable "BDparas") are ignored.

==================================================================
EXAMPLE 
==================================================================

You can run an example using the trees in "6s-bracing.trees":

```
(((A1,A2) #2, A3) [#1], ((B1,B2) #2, B3) [#1 B{0.5,0.7}]) >0.9<1.1;
```

> NOTE: we pay attention to the first tree, which is the one that enables
cross-bracing. The second tree shows another example in which the
node age constraints are considered independent (not cross-braced). This
tree is used in another example when duplication is not enabled (i.e.,
"duplication = 0" in "mcmctree.ctl", see input tree file "6s.trees").

There are two different cross-bracing events specified in the 
Newick tree in "6s-bracing.trees":

* Cross-bracing 1: the MRCA of B1 and B2 and the MRCA of A1 and A2.
They are not only cross-braced, but also a node age constraint is
being enforced: a soft-bound calibration with 0.5 as the minimum
age and 0.7 as the maximum age. The same posterior time densities 
will be inferred for both node ages.

* Cross-bracing 2: the MRCA of A1 and A2 is cross-braced with the 
MRCA of B1 and B2. There are no additional node age constraints. 
The same posterior time densities will be inferred for both node ages.

You can run MCMCtree as follows within this directory:

```
mcmctree mcmctree-bracing.ctl
```

The screen output will look more or less like the following, which
may change depending on the PAML version. The output below has
been generated with MCMCtree in PAML v4.10.10:

```
MCMCTREE in paml version 4.10.10, 20 Jan 2026

Reading options from mcmctree-bracing.ctl..

Reading main tree.
(((A1, A2), A3), ((B1, B2), B3));
   0   1   2   1   2
node  7 is mirrored by node 9.
node  8 is mirrored by node 10.
   0   1   2   0   0
   0  -1  -1   7   8

duplication dating:  2 node ages are mirrored

Reading sequence data..  1 loci


*** Locus 1 ***
ns = 6          ls = 1000
Reading sequences, sequential format..
Reading seq # 6: B3
Sequences read..
Counting site patterns..  0:00
Compressing,    592 patterns at   1000 /   1000 sites (100.0%),  0:00
Collecting fpatt[] & pose[],    592 patterns at   1000 /   1000 sites (100.0%),  0:00
Counting frequencies..
592 patterns, clean


Fossil calibration information used.
Node   7:   B (  0.9000,  1.1000,  0.0250,  0.0250 )
Node   8:   B (  0.5000,  0.7000,  0.0250,  0.0250 )

189440 bytes for conP
88 bytes for rates.

20000 burnin, sampled every 100, 2000 samples
Approximating posterior
(Settings: cleandata=1  print=1  saveconP=1)

getting initial values to start MCMC.

************
Species tree
ns = 6  nnode = 11
 father   node  name              time     sons          fossil
      9      1  A1                0.00000
      9      2  A2                0.00000
      8      3  A3                0.00000
     11      4  B1                0.00000
     11      5  B2                0.00000
     10      6  B3                0.00000
      0      7                    1.04006  ( 8 10)  B ( 0.9000, 1.1000, 0.0250, 0.0250 )
      7      8                    0.39348  ( 9  3)  B ( 0.5000, 0.7000, 0.0250, 0.0250 )
      8      9                    0.31393  ( 1  2)
      7     10                    0.39348  (11  6)
     10     11                    0.31393  ( 4  5)

(((1, 2), 3), ((4, 5), 6));
(((A1, A2), A3), ((B1, B2), B3));
(((A1: 0.313929, A2: 0.313929): 0.079549, A3: 0.393477): 0.646583, ((B1: 0.313929, B2: 0.313929): 0.079549, B3: 0.393477): 0.646583);

priors:
        mu_i ~ gammaDir(2.000, 20.000, 1.000), SD(mu_i) =  0.07071, corr(mu_i,mu_j) = -0.50000
        sigma2 ~ gammaDir(1.0000, 10.0000, 1.0000)

Initial parameters (np = 7):
  1.040060  0.393477  0.313929  0.393477  0.313929  0.145800  0.044323

lnL0 = -10313.16

Starting MCMC (np = 7) . . .
paras: 5 times, 1 mu, 1 sigma2 (& rates, kappa, alpha)
 -8% 0.20 0.39 0.91 0.00 0.00  1.037 0.596 0.215 0.596 0.215 0.499 0.457 -7433.91
(nsteps = 18)
Current Pjump:   0.20225  0.39275  0.91175  0.00000  0.00000  0.98225  0.88100  0.77725  0.92800  0.61775  0.84050  0.67600  0.69225  0.98400  0.98700  0.83050  0.96575  0.33250
Current steps:   0.09943  0.07014  0.00787  0.02631  0.05150  0.00529  0.07640  0.06771  0.01514  0.09704  0.04518  0.07917  0.06183  0.01001  0.00589  0.05570  0.01011  0.09027
New     steps:   0.06417  0.09763  0.11065  0.00026  0.00052  0.37204  0.79276  0.36419  0.26162  0.27811  0.34647  0.27847  0.23117  0.78185  0.56616  0.40084  0.36833  0.10198

 -6% 0.40 0.18 0.26 0.00 0.00  1.018 0.662 0.382 0.662 0.382 0.707 0.035 -7408.59
(nsteps = 18)
Current Pjump:   0.39725  0.18050  0.25650  0.00000  0.00000  0.06200  0.46300  0.13350  0.20150  0.15600  0.14900  0.18575  0.21050  0.07300  0.11700  0.23200  0.22875  0.25825
Current steps:   0.06417  0.09763  0.11065  0.00026  0.00052  0.37204  0.79276  0.36419  0.26162  0.27811  0.34647  0.27847  0.23117  0.78185  0.56616  0.40084  0.36833  0.10198
New     steps:   0.09067  0.05583  0.09256  0.00000  0.00001  0.07134  1.38478  0.15212  0.16817  0.13649  0.16212  0.16415  0.15573  0.17673  0.20654  0.30009  0.27154  0.08596

 -5% 0.26 0.39 0.34 0.00 0.00  1.016 0.659 0.381 0.659 0.381 0.709 0.038 -7408.65  0:05
 -4% 0.26 0.38 0.34 0.00 0.00  1.014 0.661 0.382 0.661 0.382 0.712 0.033 -7408.67
(nsteps = 18)
Current Pjump:   0.26150  0.37850  0.33575  0.00000  0.00000  0.49200  0.23075  0.42475  0.34650  0.41325  0.40225  0.35525  0.36125  0.51450  0.41225  0.30850  0.32125  0.32175
Current steps:   0.09067  0.05583  0.09256  0.00000  0.00001  0.07134  1.38478  0.15212  0.16817  0.13649  0.16212  0.16415  0.15573  0.17673  0.20654  0.30009  0.27154  0.08596
New     steps:   0.07751  0.07409  0.10580  0.00000  0.00000  0.13653  1.03062  0.23517  0.19977  0.20328  0.23289  0.20109  0.19481  0.36302  0.30660  0.31007  0.29434  0.09335

 -2% 0.33 0.26 0.28 0.00 0.00  1.016 0.663 0.383 0.663 0.383 0.711 0.032 -7408.72
(nsteps = 18)
Current Pjump:   0.32650  0.25925  0.27750  0.00000  0.00000  0.23950  0.36175  0.26050  0.27425  0.25700  0.25400  0.26975  0.26275  0.24300  0.26950  0.29000  0.29275  0.27825
Current steps:   0.07751  0.07409  0.10580  0.00000  0.00000  0.13653  1.03062  0.23517  0.19977  0.20328  0.23289  0.20109  0.19481  0.36302  0.30660  0.31007  0.29434  0.09335
New     steps:   0.08566  0.06272  0.09672  0.00000  0.00000  0.10585  1.29146  0.20016  0.18019  0.17042  0.19270  0.17801  0.16741  0.28598  0.27113  0.29812  0.28610  0.08559

  0% 0.30 0.32 0.31 0.00 0.00  1.012 0.665 0.382 0.665 0.382 0.716 0.032 -7408.67  0:12

(nsteps = 18)
Current Pjump:   0.29625  0.31650  0.31425  0.00000  0.00000  0.33375  0.26725  0.30350  0.31775  0.32275  0.32600  0.31525  0.34375  0.33225  0.31550  0.29625  0.28475  0.30100
Current steps:   0.08566  0.06272  0.09672  0.00000  0.00000  0.10585  1.29146  0.20016  0.18019  0.17042  0.19270  0.17801  0.16741  0.28598  0.27113  0.29812  0.28610  0.08559
New     steps:   0.08442  0.06680  0.10213  0.00000  0.00000  0.12012  1.13128  0.20289  0.19279  0.18576  0.21258  0.18869  0.19694  0.32277  0.28765  0.29379  0.26936  0.08593

  5% 0.30 0.31 0.28 0.00 0.00  1.013 0.662 0.380 0.662 0.380 0.711 0.034 -7408.67  0:20
 10% 0.30 0.31 0.29 0.00 0.00  1.015 0.663 0.381 0.663 0.381 0.710 0.035 -7408.68  0:28
 15% 0.30 0.30 0.29 0.00 0.00  1.015 0.663 0.381 0.663 0.381 0.710 0.034 -7408.70  0:34
 20% 0.30 0.31 0.29 0.00 0.00  1.015 0.663 0.381 0.663 0.381 0.711 0.034 -7408.69  0:40
 25% 0.30 0.31 0.28 0.00 0.00  1.015 0.663 0.382 0.663 0.382 0.710 0.035 -7408.67  0:47
 30% 0.30 0.31 0.28 0.00 0.00  1.015 0.663 0.382 0.663 0.382 0.710 0.035 -7408.66  0:53
 35% 0.30 0.31 0.28 0.00 0.00  1.015 0.663 0.381 0.663 0.381 0.710 0.034 -7408.69  1:04
 40% 0.30 0.31 0.29 0.00 0.00  1.015 0.663 0.381 0.663 0.381 0.710 0.034 -7408.69  1:14
 45% 0.30 0.31 0.29 0.00 0.00  1.015 0.663 0.381 0.663 0.381 0.711 0.034 -7408.69  1:20
 50% 0.30 0.31 0.28 0.00 0.00  1.015 0.663 0.382 0.663 0.382 0.711 0.034 -7408.69  1:27
 55% 0.30 0.31 0.28 0.00 0.00  1.015 0.663 0.382 0.663 0.382 0.711 0.033 -7408.69  1:33
 60% 0.30 0.31 0.29 0.00 0.00  1.015 0.663 0.382 0.663 0.382 0.711 0.033 -7408.70  1:39
 65% 0.30 0.31 0.29 0.00 0.00  1.015 0.663 0.381 0.663 0.381 0.711 0.033 -7408.69  1:45
 70% 0.30 0.31 0.29 0.00 0.00  1.015 0.663 0.382 0.663 0.382 0.711 0.033 -7408.69  1:51
 75% 0.30 0.31 0.29 0.00 0.00  1.015 0.663 0.382 0.663 0.382 0.711 0.033 -7408.68  2:00
 80% 0.30 0.31 0.29 0.00 0.00  1.015 0.663 0.382 0.663 0.382 0.711 0.033 -7408.69  2:07
 85% 0.30 0.31 0.29 0.00 0.00  1.015 0.663 0.382 0.663 0.382 0.711 0.033 -7408.69  2:13
 90% 0.30 0.31 0.29 0.00 0.00  1.015 0.663 0.382 0.663 0.382 0.711 0.033 -7408.69  2:21
 95% 0.30 0.31 0.29 0.00 0.00  1.015 0.663 0.381 0.663 0.381 0.711 0.033 -7408.69  2:30
100% 0.30 0.31 0.29 0.00 0.00  1.015 0.663 0.382 0.663 0.382 0.711 0.033 -7408.69  2:41

Time used:  2:41
Summarizing MCMC samples . ..

Data file has a header line.
2001 records, 9 variables
Collecting mean, median, min, max, percentiles, etc.
                             9/     9 done   2:42


Posterior means (95% Equal-tail CI) (95% HPD CI) HPD-CI-width

t_n7           1.0155 ( 0.9106,  1.0992) ( 0.9178,  1.1038)  0.1860
t_n8           0.6626 ( 0.5839,  0.7074) ( 0.5984,  0.7133)  0.1149
t_n9           0.3820 ( 0.3138,  0.4578) ( 0.3120,  0.4555)  0.1435
t_n10          0.6626 ( 0.5839,  0.7074) ( 0.5984,  0.7133)  0.1149
t_n11          0.3820 ( 0.3138,  0.4578) ( 0.3120,  0.4555)  0.1435
mu             0.7111 ( 0.5910,  0.8201) ( 0.5910,  0.8201)  0.2290
sigma2         0.0327 ( 0.0011,  0.1247) ( 0.0000,  0.0977)  0.0976
lnL        -7408.7696 (-7412.8990, -7405.8580) (-7412.3570, -7405.6170)  6.7400

mean    1.0155  0.6626  0.3820  0.6626  0.3820  0.7111  0.0327  -7408.7696
Eff     1.1389  1.5496  0.9180  1.5496  0.9180  0.9534  0.8147  1.0434
time prior: Birth-Death-Sampling
rate prior: Log-Normal

Time used:  2:42
```

You can see the following section in the screen output...


```
Reading main tree.
(((A1, A2), A3), ((B1, B2), B3));
   0   1   2   1   2
node  7 is mirrored by node 9.
node  8 is mirrored by node 10.
   0   1   2   0   0
   0  -1  -1   7   8

duplication dating:  2 node ages are mirrored
```

... which shows which nodes are being cross-braced and 
confirms that variable "duplication" has been enabled.

In addition, you can see that the inferred mean posterior
times and CIs for cross-braced nodes are the same: nodes 8
and 10 were cross-braced and so were 9 and 11:

```
Posterior means (95% Equal-tail CI) (95% HPD CI) HPD-CI-width

t_n7           1.0155 ( 0.9106,  1.0992) ( 0.9178,  1.1038)  0.1860
t_n8           0.6626 ( 0.5839,  0.7074) ( 0.5984,  0.7133)  0.1149
t_n9           0.3820 ( 0.3138,  0.4578) ( 0.3120,  0.4555)  0.1435
t_n10          0.6626 ( 0.5839,  0.7074) ( 0.5984,  0.7133)  0.1149
t_n11          0.3820 ( 0.3138,  0.4578) ( 0.3120,  0.4555)  0.1435
mu             0.7111 ( 0.5910,  0.8201) ( 0.5910,  0.8201)  0.2290
sigma2         0.0327 ( 0.0011,  0.1247) ( 0.0000,  0.0977)  0.0976
lnL        -7408.7696 (-7412.8990, -7405.8580) (-7412.3570, -7405.6170)  6.7400
```

You can check that the posterior time densities are also the same
in the "mcmc.txt" file that shall be generated.

You can compare the inferred divergence times when enabling cross-bracing
to those estimated when nodes are calibrated independently. You can run
MCMCtree using the control file "mcmctree.ctl", which uses the input tree
file "6s.trees":

```
mcmctree mcmctree.ctl
```

You will see that no screen output refers to cross-braced nodes and
the mean posterior time densities for nodes 8 and 10 differ as the 
node age constraints had not been cross-braced.
