MCMCTREE
Ziheng Yang

4 October 2005, last modified October 2013


(A) You can use the files in this folder to run mcmctree to duplicate the
results of Yang and Rannala (2006: table 3) and Rannala and Yang
(2007, table table 2).  

    ..\..\bin\mcmctree

Make your window wider (100 columns) before running the program.  See
the section on mcmctree in the paml manual for a description of the
program.  The default control file with usedata = 1 means the use of
the exact likelihood calculation.


(B) To use the approximate likelihood calculation, do the following.

(b1) Copy baseml.exe from the paml/bin folder to the current folder.

(b2) Edit mcmctree.ctl to use usedata = 3.  Run mcmctree.  

     ..\..\bin\mcmctree

     This generates a file named out.BV.  Rename it in.BV

(b3) Edit mcmctree.ctl to use usedata = 2.  Run mcmctree again.  

     ..\..\bin\mcmctree

The exact and approximate likelihood calculations generate very
similar results for this dataset.


(C)  To run the InfiniteSites program, type one of the following

     ..\..\bin\InfiniteSites

     ..\..\bin\InfiniteSites mcmctree.Infinitesites.ctl


InfiniteSites reads the same control file as mcmctree.  In addition it
reads the file FixedDsClock23.txt if clock = 2 or 3.  It calculates
the posterior according to the theory in the 2006 and 2007 papers,
assuming that the branch lengths are fixed.  See also dos Reis & Yang
(2013).


References

dos Reis M, Yang Z (2013) The unbearable uncertainty of Bayesian divergence 
time estimation. J. Syst. Evol. 51:30-43.

Rannala, B., and Z. Yang. 2007. Inferring speciation times under an
episodic molecular clock. Syst. Biol. 56: 453-466.

Yang, Z., and B. Rannala. 2006. Bayesian estimation of species
divergence times under a molecular clock using multiple fossil
calibrations with soft bounds. Mol. Biol. Evol.: 23: 212-226.
