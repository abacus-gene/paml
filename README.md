# Phylogenetic Analysis by Maximum Likelihood

**PAML** (for **Phylogenetic Analysis by Maximum Likelihood**) is a package of programs for phylogenetic analyses of DNA or protein sequences using maximum likelihood. It is maintained by [**Ziheng Yang**](http://abacus.gene.ucl.ac.uk/) and distributed under the GNU GPL v3.

Before running any of the programs included in the `PAML` software, please go through the different sections of this Wiki, where a detailed explanation of how to install and run the programs is given. One of the main problems users face is related to **data formatting**. In that way, please make sure that you go through the **section [`Data formatting` in the PAML Wiki](https://github.com/abacus-gene/paml/wiki/Data-formatting)** before running any `PAML` program and format your input data files and control files accordingly.

> [!IMPORTANT]
>
> * Problems with input data, control files, error/warning messages output by the `PAML` programs (which inform users about formatting issues in their input files or wrong settings in their control files), and general questions should be posted in the **[`PAML` discussion group](https://groups.google.com/g/pamlsoftware?pli=1)**. **Before posting a message, please [use the search tool in the `PAML` discussion group](https://groups.google.com/g/pamlsoftware) to check whether your question/s have already been asked by other `PAML` users. You should also refer to [the FAQs document](https://github.com/abacus-gene/paml/blob/master/doc/pamlFAQs.pdf)**. If you still cannot find an answer to your question/s, please post them in the [**`PAML` discussion group**](https://groups.google.com/g/pamlsoftware).
> * Please, **do not** paste the screen output and the error/warning that you get without any reference when you are requesting help. Instead, make sure that you attach the input data and the control files you are using, explain how you ran `PAML`, which `PAML` version you are using, etc. In that way, your issue will be **much easier and faster to troubleshoot**. You can find more information about **how to report your warning/error message [on this website](https://uxwritinghub.com/error-message-examples/)**.
> * **Please, only raise an issue on this GitHub repository when you experience technical problems such as compiling issues, programs aborting or not running at all, etc.**

## Installation

ANSI C source codes and executable files are distributed for UNIX, Linux, and Mac OS X (see the [latest stable release available on this GitHub repository](https://github.com/abacus-gene/paml/releases)). `PAML` is not good for tree making, although it may be used to estimate parameters and test hypotheses to study the evolutionary process once you have reconstructed trees using other programs such as [`RAxML-NG`](https://github.com/amkozlov/raxml-ng), [`IQ-TREE`](http://www.iqtree.org/), [`PAUP*`](https://paup.phylosolutions.com/), [`PHYLIP`](https://evolution.genetics.washington.edu/phylip/doc/main.html), [`PhyML`](http://www.atgc-montpellier.fr/phyml/), etc.

To download and install `PAML` software, please follow the next links in the wiki:

* [Downloading and installing `PAML`](https://github.com/abacus-gene/paml/wiki/Installation)
  * [`PAML` for Linux](https://github.com/abacus-gene/paml/wiki/Installation#paml-for-linux)
    * [How to export paths in Linux](https://github.com/abacus-gene/paml/wiki/Installation#exporting-paths-linux)
  * [`PAML` for Mac OS X](https://github.com/abacus-gene/paml/wiki/Installation#paml-for-mac-os-x)
    * [How to export paths in Mac OS X](https://github.com/abacus-gene/paml/wiki/Installation#exporting-paths-mac-os-x)
  * [`PAML` for Windows 10 or later](https://github.com/abacus-gene/paml/wiki/Installation#paml-for-windows-10-or-later)
    * [How to export paths using the Windows Subsystem for Linux](https://github.com/abacus-gene/paml/wiki/Installation#exporting-paths-wsl)
  * [Old compiled versions of `PAML`](https://github.com/abacus-gene/paml/wiki/Installation#old-compiled-versions-of-paml)
    * [`PAML` for Mac OS X](https://github.com/abacus-gene/paml/wiki/Installation#old-paml-versions-for-mac-os-x)
    * [`PAML` for Windows 9x/NT/2000/XP/Vista 7](https://github.com/abacus-gene/paml/wiki/Installation#old-paml-versions-for-windows-9xnt2000xpvista-7)
  * [Graphical user interface: PAML-X](https://github.com/abacus-gene/paml/wiki/Installation#graphical-user-interface-paml-x)
 
> If you are a macOS user and you are experiencing issues such as "trace trap" or "segmentation fault" when running `PAML` programs, please download the `dev` branch and follow the installation instructions given in the section [`PAML` for Mac OS X in the PAML Wiki](https://github.com/abacus-gene/paml/wiki/Installation#paml-for-mac-os-x). 

## Documentation

The [PAML Wiki](https://github.com/abacus-gene/paml/wiki/) is still under construction :wrench: -- we are slowly migrating the [PAML documentation in PDF format](https://github.com/abacus-gene/paml/blob/master/doc/pamlDOC.pdf) to a more interactive and engaging tutorial!

In the meantime, you can access the following sections, which are also detailed in PAML Wiki home section:

* [Overview](https://github.com/abacus-gene/paml/wiki/Overview)
  * [`PAML` documentation](https://github.com/abacus-gene/paml/wiki/Overview#paml-documentation)
  * [What `PAML` programs can do](https://github.com/abacus-gene/paml/wiki/Overview#what-paml-programs-can-do)
  * [What `PAML` programs cannot do](https://github.com/abacus-gene/paml/wiki/Overview#what-paml-programs-cannot-do)
* [Data formatting](https://github.com/abacus-gene/paml/wiki/Data-formatting)
  * [Sequence data file format](https://github.com/abacus-gene/paml/wiki/Data-formatting#sequence-data-file-format)
    * [Sequential and interleaved formats](https://github.com/abacus-gene/paml/wiki/Data-formatting#sequential-and-interleaved-formats)
      * [Species/sequence names](https://github.com/abacus-gene/paml/wiki/Data-formatting#speciessequence-names)
      * [Option `G`](https://github.com/abacus-gene/paml/wiki/Data-formatting#option-g)
      * [Option `G` for codon sequences (`codeml` with `seqtype = 1`)](https://github.com/abacus-gene/paml/wiki/Data-formatting#option-g-for-codon-sequences-codeml-with-seqtype--1)
    * [Site pattern counts](https://github.com/abacus-gene/paml/wiki/Data-formatting#site-pattern-counts)
  * [Tree file format and representations of tree topology](https://github.com/abacus-gene/paml/wiki/Data-formatting#tree-file-format-and-representations-of-tree-topology)
    * [Parenthesis notation](https://github.com/abacus-gene/paml/wiki/Data-formatting#parenthesis-notation)
    * [Tree files produced by `PAUP` and `MacClade`](https://github.com/abacus-gene/paml/wiki/Data-formatting#tree-files-produced-by-paup-and-macclade)
    * [Branch or node labels](https://github.com/abacus-gene/paml/wiki/Data-formatting#branch-or-node-labels)
    * [Divergence date symbol `@`](https://github.com/abacus-gene/paml/wiki/Data-formatting#divergence-date-symbol-)
    * [Branch representation of tree topology](https://github.com/abacus-gene/paml/wiki/Data-formatting#branch-representation-of-tree-topology)
* [Substitution models](https://github.com/abacus-gene/paml/wiki/Substitution-models)
  * [Nucleotide substitution models](https://github.com/abacus-gene/paml/wiki/Substitution-models#nucleotide-substitution-models)
  * [Codon substitution models](https://github.com/abacus-gene/paml/wiki/Substitution-models#codon-substitution-models)
    * [Branch models](https://github.com/abacus-gene/paml/wiki/Substitution-models#branch-models)
    * [Site models](https://github.com/abacus-gene/paml/wiki/Substitution-models#site-models)
      * [Suzuki and Gojobori’s (1999) method for detecting sites under positive selection](https://github.com/abacus-gene/paml/wiki/Substitution-models#suzuki-and-gojoboris-1999-method-for-detecting-sites-under-positive-selection)
    * [Branch-site models](https://github.com/abacus-gene/paml/wiki/Substitution-models#branch-site-models)
    * [Clade models](https://github.com/abacus-gene/paml/wiki/Substitution-models#clade-models)
    * [Mutation-selection models](https://github.com/abacus-gene/paml/wiki/Substitution-models#mutation-selection-model)
  * [Amino acid substitution models](https://github.com/abacus-gene/paml/wiki/Substitution-models#amino-acid-substitution-models)
* **`PAML` programs**
  * [BASEML](https://github.com/abacus-gene/paml/wiki/BASEML)
  * [BASEMLG](https://github.com/abacus-gene/paml/wiki/BASEMLG)
  * [CODEML](https://github.com/abacus-gene/paml/wiki/CODEML)
  * [Evolver]([Evolver](https://github.com/abacus-gene/paml/wiki/Evolver))
  * [Infinitesites](https://github.com/abacus-gene/paml/wiki/Infinitesites)
  * [MCMCtree]([MCMCtree](https://github.com/abacus-gene/paml/wiki/MCMCtree))
  * [yn00](https://github.com/abacus-gene/paml/wiki/yn00)

## Citing `PAML`

**If you use `PAML`**, please cite the following:

* [Yang, Z (1997). PAML: a program package for phylogenetic analysis by maximum likelihood. *Comput. Appl. Biosci.* 13, 555-556](http://abacus.gene.ucl.ac.uk/ziheng/pdf/1997YangCABIOSv13p555.pdf).
* [Yang, Z (2007). PAML 4: Phylogenetic Analysis by Maximum Likelihood. *Mol. Biol. Evol.* 24, 1586-1591](https://academic.oup.com/mbe/article-pdf/24/8/1586/3853532/msm088.pdf).

### In addition...

**If you use the `PAML` program `MCMCtree`**, please cite the following papers if you have used/run...

* ... the **approximate likelihood calculation** to speed up analyses with phylogenomic datasets (calculating branch lengths, Hessian, and gradient):
  * [dos Reis M and Yang Z (2011). Approximate likelihood calculation for Bayesian estimation of divergence times. *Mol. Biol. Evol.* 28, 2161-2172](http://abacus.gene.ucl.ac.uk/ziheng/pdf/2011dosReisYangMBEv28p2161.pdf).
* ... **Bayesian model selection** analyses:
  * [dos Reis M, et al. (2018). Using phylogenomic data to explore the effects of relaxed clocks and calibration strategies on divergence time estimation: Primates as a test case. *Syst. Biol.* 67, 594–615](http://abacus.gene.ucl.ac.uk/ziheng/pdf/2018dosReis.Primates.pdf). Also, please cite the [`mcmc3r` R package](https://github.com/dosreislab/mcmc3r) if you use it ([see tutorial via this link](https://dosreislab.github.io/2017/10/24/marginal-likelihood-mcmc3r.html)).
* ... the models for **continuous morphological characters** implemented in `MCMCtree`:
  * [Álvarez-Carretero S, et al. (2019). Bayesian estimation of species divergence times using correlated quantitative characters. *Syst. Biol.* 68, 967–986](http://abacus.gene.ucl.ac.uk/ziheng/pdf/2019Alvarez-CarreteroSB.pdf). Also, please cite the [`mcmc3r` R package](https://github.com/dosreislab/mcmc3r) if you use it ([see tutorial via this link](https://github.com/dosreislab/mcmc3r/blob/master/vignettes/Reproduce_Carnivora_analysis.Rmd)).
* ... the protocol **Bayesian Molecular Clock Dating Using Genome-Scale Datasets**:
  * [dos Reis M and Yang, Z (2019). Bayesian Molecular Clock Dating Using Genome-Scale Datasets. *In*: Anisimova, M. (eds) Evolutionary Genomics. Methods in Molecular Biology, vol 1910. Humana, New York, NY](https://link.springer.com/protocol/10.1007/978-1-4939-9074-0_10). You can also [access this chapter and the code used throughout the protocol on the `divtime` GitHub repository maintained by Mario dos Reis](https://github.com/mariodosreis/divtime).
* ... the protocol **Environmental Microbial Evolution: Methods and Protocols**:
  * [dos Reis M (2022) *In*: Haiwei Luo (ed.) Environmental Microbial Evolution: Methods and Protocols. Methods in Molecular Biology, vol 2569. Humana, New York, NY](https://link.springer.com/protocol/10.1007/978-1-0716-2691-7_1). You can [access the code used throughout the protocol on the `microdiv` GitHub repository maintained by Mario dos Reis](https://github.com/dosreislab/microdiv).
* ... the **Bayesian sequential-subtree (BSS) approach** and/or the scripts to **fit skew-*t* distributions to fossil calibrations**:
  * [Álvarez-Carretero S, et al. (2022) A species-level timeline of mammal evolution integrating phylogenomic data. *Nature* 602, 263–267](https://rdcu.be/cDHW7). You can find [the tutorial to reproduce these analyses on the `mammals_dating` GitHub repository maintained by Sandra Álvarez-Carretero](https://github.com/sabifo4/mammals_dating).

**If you use the `PAML` program `CODEML`**, please cite the following papers if you have used/run...

* ... the protocol **A Beginners guide to estimating the non-synonymous to synonymous rate ratio of all protein-coding genes in a genome**:
  * [Jeffares DC, Tomiczek B, Sojo V, dos Reis M (2015). A Beginners Guide to Estimating the Non-synonymous to Synonymous Rate Ratio of all Protein-Coding Genes in a Genome. *In*: Peacock, C. (eds) Parasite Genomics Protocols. Methods in Molecular Biology, vol 1201. Humana Press, New York, NY.](https://link.springer.com/protocol/10.1007/978-1-4939-1438-8_4).
* ... the protocol **Beginner's guide on the use of PAML to detect positive selection** and/or the [corresponding GitHub tutorial on the `positive-selection` repository](https://github.com/abacus-gene/paml-tutorial/tree/main/positive-selection):
  * [Álvarez-Carretero S, Kapli P, Yang Z (2023). Beginner's guide on the use of PAML to detect positive selection, Mol Biol Evol, 40(4):msad041.](https://doi.org/10.1093/molbev/msad041). Remember to read the [supplementary material where we discuss (i) analyses and checks you should carry out before running tests of positive selection with CODEML, (ii) gene tree VS species tree, and (iii) the usage of rooted and unrooted trees](http://abacus.gene.ucl.ac.uk/ziheng/pdf/2023Alvarez-Carretero-codeml-SI.pdf).

## Additional information

Changes and bug fixes until v4.10.1 were documented in the file [doc/pamlHistory.txt](https://github.com/abacus-gene/paml/blob/master/doc/pamlHistory.txt). Changes in later versions have been documented for each release in [the `releases` section of this GitHub repository](https://github.com/abacus-gene/paml/releases).
