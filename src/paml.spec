Summary: Phylogenetic analyses of DNA or protein sequences using maximum likelihood
Name: paml
Version: 3.13
Release: 2
License: free for academic use only
Group: Applications/Scientific
URL: http://abacus.gene.ucl.ac.uk/software/paml.html
Source: ftp://abacus.gene.ucl.ac.uk/pub/paml/%{name}%{version}.tar.Z
BuildRoot: %{_tmppath}/%{name}-%{version}-%{release}-root


%description
PAML (for Phylogentic Analysis by Maximum Likelihood) contains a few
programs for model fitting and phylogenetic tree reconstruction using
nucleotide or amino-acid sequence data.  PAML is distributed free of
charge for academic use only.

#*** WARNING ***
#When using codeml or evolver, they may complain "XXXX.dat not found".
#The missing file should be found in %{_docdir}/%{name}-%{version}.
#Copy it to the directory where your input data exists.  Or make a sym-link.
#See the bottom of p.5 of %{_docdir}/%{name}-%{version}/doc/pamlDOC.pdf
#for details.


%prep
[ "$RPM_BUILD_ROOT" != "/" ] && rm -rf $RPM_BUILD_ROOT
%setup -q -n %{name}%{version}


%build
export CFLAGS=$RPM_OPT_FLAGS
cd src; make -e


%install
EXECUTABLES='baseml basemlg chi2 codeml codemlsites evolver mcmctree pamp yn00'
install --mode=755 -d $RPM_BUILD_ROOT/%{_prefix}/bin
cd src; install --mode=755 $EXECUTABLES  $RPM_BUILD_ROOT/%{_prefix}/bin


%clean
[ "$RPM_BUILD_ROOT" != "/" ] && rm -rf $RPM_BUILD_ROOT


%files
%defattr(-,root,root)
%doc *.trees GeneticCode.txt *.dat README.txt *.nuc *.ctl doc examples paupblock paupend paupstart *.aa

%{_prefix}/bin/*


%changelog
* Tue May 28 2002 Hunter Matthews <thm@duke.edu>
- Spec file tweaks, updates.

* Tue Apr  9 2002 Naoki Takebayashi <ntakebay@bio.indiana.edu> [3.12-1]
- updated to 3.12
- slight change in %file

* Tue Mar 5 2002 Naoki Takebayashi <ntakebay@bio.indiana.edu> [3.1-1]
- updated to 3.1

* Sat May 5 2001 Naoki Takebayashi <ntakebay@bio.indiana.edu> [3.0d-1]
- first release of RPM package
