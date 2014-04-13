Codon Simulation Programs

README.txt 

July 2004


(A) evolver: Versions for simulating codon sequences

This folder contains a few programs (windows executables or the .exe
files) for simulating codon sequences under codon models.  If you use
unix/linux/osx, you delete the .exe files, and compile the programs
yourself, as follows.  Note that evolver.c and associated programs are
in the folder src/.

evolverNSbranches:
     cc -O2 -o evolverNSbranches -DCodonNSbranches evolver.c tools.c -lm

evolverNSsites:
     cc -O2 -o evolverNSsites -DCodonNSsites evolver.c tools.c -lm

evolverNSbranches takes input from MCcodonNSbranches.dat and simulate
codon sequences with different omega ratios (dN/dS) for different
branches on the tree, that is, under the model of Yang (1998).  Check
that file for details, a copy of which is included in this folder.
Note that you specify the branch lengths in the tree using the symbol
":", as in other programs, and you specify the omega ratio for the
branch using the symbol "#", a notation not used by other programs.

evolverNSsites takes input from MCcodonNSsites.dat and simulates codon
sequences with a few classes of sites with different omega ratios
(dN/dS).  This is the discrete (M3) model of Yang et al. (2000).  You
can of course use this format to simulate other models as well,
because they are its special cases.  For example, model M2a
(PositiveSelecion) uses three site classes with omega < 1, = 1, and >
1, so you can use a omega = 1 in the data file.  Check the file for
details, a copy of which is included in this folder.

The branch length under a codon model is defined as the expected
number of nucleotide substitutions per codon.


(B) PositiveSites: Program for summarizing simulation results
concerning inference of sites under positive selection under site
models.

When you use evolverNSsites (see above) to simulate data sets, the
program will generate a file called siterates, listing the sites that
have omega > 1 in the simulation and are thus truly under positive
selection.  When the data sets are analyzed using codeml, codeml will
print results in a file called mlc, which may include a list of sites
inferred to be under positive selection with calculated empirical
Bayes posterior probabilities.  PositiveSites compares those siterates
and mlc to calculate a few measures of performance of the codeml
analysis.  This is the kind of analyses used in the simulation studies
of Anisimova et al. (2002) and Wong et al. (2004).

Currently codeml uses two procedures to calculate the posterior
probabilities: the Naive Empirical Bayes (NEB: Nielsen & Yang 1998; Yang et
al. 2000) and Bayes empirical Bayes (BEB: Yang, Wong, Nielsen
submitted).  NEB uses maximum likelihood estimates of parameters
without accounting for their errors, while BEB uses a prior to correct
for uncertainties in the parameter estimates.  codeml prints out the
NEB results, followed by the BEB results in both mlc and rst, while
PositiveSites processes mlc only.  The included windows executables
PositiveSitesNEB and PositiveSitesBEB are for processing the NEB and
BEB results, respectively.  To compile for unix, use the following
command to compile the included source file.

   cc -DNEB -o PositiveSitesNEB -O2 PositiveSites.c -lm
   cc -DBEB -o PositiveSitesBEB -O2 PositiveSites.c -lm

The program relies on a unique string in the codeml output mlc to
identify results for each simulated replicate.  NEB relies on the
string "Naive Empirical".  BEB relies on the string "Bayes Empirical".
The source file PositiveSites.c is rather short so you can check the
source file.

PositiveSitesNEB and PositiveSitesBEB work with codeml in paml 3.14.

You run the programs in the following order.  The first two commands
obviously simulate data sets and analyze them.  The third line
summarizes the results, and you have to tell PositiveSites the number
of codons in the sequence (100) and the number of simulated replicates
on the command line.

	evolverNSsites
	codeml
	PositiveSitesNEB 100 500


Now enter the measurements.  PositiveSites calculates the accuracy,
power, and false positive rate of codeml inference of positive
selection sites.  The measures are defined as follows (Anisimova et
al. 2002; Wong et al. 2004).

                             codeml inference
                              +        -         Total
        evolver    +          N++      N+-       N+.
                   -          N-+      N--       N-.
                   Total      N.+      N.-       N

   Accuracy      = N++/N.+
   Power         = N++/N+.
   FalsePositive = N-+/N-.

The program collects N++ (NmatchB & NmatchC), N+. (NEvolver), and N.+
(NCodemlB & NcodemlC), and then calculates the three measures as
above.  Note that codeml inference depends on a cutoff P, hence the B
(for binned) and C (for cumulative) difference.  All proportions are
calculated as the ratio of averages, taking the ratio after counting
sites over replicate data sets.  Output is on the screen.  There are
more descriptions of those measures in Yang, Wong & Nielsen
(submitted).


(C) Testing the NEB algorithm.  You can confirm the calculation of
codeml and also the simulation program by fixing parameters at their
true values (values used in the simulation).  For example, you can use
a small tree, simulated 100 codons in sequence and 500 datasets, using
fixed parameters for the w distribution (say under model 2A).  They
analyze the data using codeml but with all parameters including branch
lengths, kappa, and parameters in the w distribution at their true
values.  You can do this by using in.codeml.  See the paml
Documentation about how to set this up.  Then when codeml says a site
is under positive selection with probability 0.9, that site is under
positive selection with probability 0.9.  The following is a sample
output from PositiveSitesNEB in one such simulation.  The second line
means that 1645 sites were listed by codeml as being under positive
selection with posterior probability in the bin 0.55-0.60.  Among
them, 55.1% of them are truly under positive selection (on the evolver
list).  This proportion is in the column "AccuracyBin", accuracy
binned.  If both evolver and codeml are correct, we should expect this
proportion to be between 0.55 and 0.60, and in this case the match
seems to be good enough.  There are 65422 sites with posterior
probability exceeding 0.55 according to codeml, and 93.2% of them are
truly under positive selection.  This proportion is in the column
"AccuracyCum", accuracy accumulated.  There is no theory to predict
what this proportion is except that it should exceed 0.55.  The last
column "Power" is the power of the codeml analysis.  At the 55%
cutoff, 88.2% of the 69098 sites truly under positive selection were
picked up by codeml.  Another measure, used by Wong et al. (2004), is
the false positive rate.  This is defined as the number of sites
falsely identified to be under posiitve selection by codeml among all
sites not under positive selection.

 P              AccuracyBin      Pcut     AccuracyCum     Power

 0.50 - 0.55:   0.519 ( 1863)    >0.50:   0.920 (67285)   0.896 (69098)
 0.55 - 0.60:   0.551 ( 1645)    >0.55:   0.932 (65422)   0.882 (69098)
 0.60 - 0.65:   0.611 ( 1742)    >0.60:   0.941 (63777)   0.869 (69098)
 0.65 - 0.70:   0.660 ( 1757)    >0.65:   0.951 (62035)   0.854 (69098)
 0.70 - 0.75:   0.713 ( 1813)    >0.70:   0.959 (60278)   0.837 (69098)
 0.75 - 0.80:   0.775 ( 1995)    >0.75:   0.967 (58465)   0.818 (69098)
 0.80 - 0.85:   0.821 ( 2466)    >0.80:   0.974 (56470)   0.796 (69098)
 0.85 - 0.90:   0.876 ( 3057)    >0.85:   0.981 (54004)   0.766 (69098)
 0.90 - 0.95:   0.925 ( 4785)    >0.90:   0.987 (50947)   0.728 (69098)
 0.95 - 0.99:   0.975 (10202)    >0.95:   0.993 (46162)   0.664 (69098)
 0.99 - 1.00:   0.998 (35960)    >0.99:   0.998 (35960)   0.520 (69098)

A test like this only confirms that the programs are likely to be
correctly coded.  It does not tell much about the performance of the
codeml analysis of real data sets, as in real data analysis, the model
used by codeml might be wrong and also the true values of parameters
are unknown.


A similar test can be conducted to confirm the BEB calculation.  In
this case the data should be simulated by drawing parameter values
(such as the proportions and omega ratios) from the prior.  We (Yang,
Wong, Nielsen submitted) inlcuded an example of this test, which is
also used to illustrate the measures of performance.  There is
probably no need for you to simulate under the prior, but if you want
a that simulation program, you can contact me.


References

Anisimova, M., J. P. Bielawski, and Z. Yang. 2002. Accuracy and power
of Bayes prediction of amino acid sites under positive
selection. Molecular Biology and Evolution 19:950-958.

Nielsen, R., and Z. Yang. 1998. Likelihood models for detecting
positively selected amino acid sites and applications to the HIV-1
envelope gene. Genetics 148:929-936.

Wong, W. S. W., Z. Yang, N. Goldman, and R. Nielsen. 2004. Accuracy
and power of statistical methods for detecting adaptive evolution in
protein coding sequences and for identifying positively selected
sites. Genetics:in press.

Yang, Z. 1998. Likelihood ratio tests for detecting positive selection
and application to primate lysozyme evolution. Molecular Biology and
Evolution 15:568-573.

Yang, Z., R. Nielsen, N. Goldman, and
A.-M. K. Pedersen. 2000. Codon-substitution models for heterogeneous
selection pressure at amino acid sites. Genetics 155:431-449.
