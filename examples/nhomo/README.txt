README.txt
Ziheng, 29 January 2022

This dataset is used to illustrate the nonhomogeneous models implemented 
in baseml (nhomo = 4 and 5).  Matsumoto and Akashi (2018) simulated 10 
replicate datasets under different parameter settings to emulate base 
composition evolution at synonymous sites in protein-coding genes in six 
species from the Drosophila melanogaster subgroup.  There are two 
sequences for each species (the so-called collapsed sequences), with 12 
sequences in total.  The dataset used here is replicate 1 of the first 
scenario of stationary evolution, and we used the dataset for MCU0.95, 
with high GC.

Two model settings were used in the analyses: 
(X) model=7, fix_kappa=0, nhomo=4 and 
(Y) model=7, fix_kappa=2, nhomo=5, 

Setting X (nhomo=4) assumes that every branch in the tree has its own vector 
of base-composition parameters.  
In setting Y (nhomo=5), sequences 1 and 2 are assumed to have the same 
frequency parameters, as are sequences 3 & 4, 5 & 6, 7 & 8, 9 & 10, and 11
& 12.  The model is specified by labelling branches on the tree.  Please see
the tree file tree-nhomo5.txt.

For more details, search for 'nhomo' in pamlDOC.pdf.

Thanks to Tomo for providing the example data file and notes.

Note (for TM):
dataset_1: 16_collapse_MCU=0.5_stationary_10samples (10 reps)
dataset_2: 16_collapse_MCU=0.7_stationary_10samples (10 reps)
dataset_3: 16_collapse_MCU=0.95_stationary_10samples (10 reps)  [replicate 1 is used here.]


References

Matsumoto T, Akashi H, Yang Z. 2015. Evaluation of ancestral sequence
reconstruction methods to infer nonstationary patterns of nucleotide
substitution. Genetics 200:873-890.

Matsumoto T, Akashi H. 2018. Distinguishing among evolutionary forces
acting on genome-wide base composition: computer simulation analysis
of approximate methods for inferring site frequency spectra of derived
mutations. G3 (Bethesda) 8:1755-1769.
