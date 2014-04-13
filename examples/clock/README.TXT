readme.txt
November 2000

This directory contains example files for global and local molecular
clock analysis.  You can print out the control (.ctl) files, have a
look at the tree files (.trees) and result files (.rst).  Note that
the data and result files might be very large, so do not print them
without checking first.

The following files are for duplicating the results in table 3 "one
rate" in Yoder and Yang (2000 MBE 17: 1081-1090): Estimation of dates
from third codon positions.  You run the program by the command
baseml.  When the iteration finishes, the program will ask you to
input "reference node & date", input 37 25 and then -1 -1.

To run the global clock, chose clock = 1.

 baseml.ctl
 mtYoder3.nuc
 DatingNUC.trees
 mtYoder3.rst


The following files are for duplicating the amino acid based analysis
in the Yoder and Yang paper.

 codeml.ctl
 mtYoder.aa
 DatingAA.trees
 mtYoderAA.rst

The following set of files are for the TipDate analysis.  Look at the result file exampleTipdate.rst for more information.

 exampleTipDate.ctl
 exampleTipDate.phy
 exampleTipDate.trees
 exampleTipDate.rst

