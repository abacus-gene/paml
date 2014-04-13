Notes by Ziheng, 12 August 2003

This is about codon models for McDonald & Kreithman-style test using
the ape mt data used in Hasegawa et al. (1998).

mtCDNAape.nuc
mtCDNAape.trees

The following are model specifications and log likelihood values, as
well as MLEs.  "W" and "B" refer to within- and between-species
branches.  "R" and "C" refer to radical and conserved amino acid
changes.


np=11 model=0 aaDist=0: lnL=-20486.034301  k=20.74839  w=0.04414

np=12 model=2 aaDist=0: lnL=-20444.099676  k=21.59077  wW=0.28638  wB=0.03693

np=12 model=0 aaDist=7: lnL=-20482.229434  k=20.52018  wR=0.02745  wC=0.04658

np=14 model=2 aaDist=7: lnL=-20440.382774  k=21.36004  wWR=0.15012  wWC=0.30470  wBR=0.02380  wBC=0.03885



Reference

Hasegawa, M., Y. Cao and Z. Yang. 1998. Preponderance of slightly
deleterious polymorphism in mitochondrial DNA: replacement/synonymous
rate ratio is much higher within species than between
species. Molecular Biology and Evolution 15:1499-1505.

Yang, Z., R. Nielsen and M. Hasegawa. 1998. Models of amino acid
substitution and applications to mitochondrial protein
evolution. Molecular Biology and Evolution 15:1600-1611.
