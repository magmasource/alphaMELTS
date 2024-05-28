alphaMELTS
==========

alphaMELTS started as a text menu-driven interface to subroutine versions of the MELTS ([Ghiorso & Sack, 1995](http://dx.doi.org/10.1007/s004100050036)), pMELTS ([Ghiorso et al., 2002](http://dx.doi.org/10.1029/2001GC000217)), and pHMELTS ([Asimow et al., 2004](http://dx.doi.org/10.1029/2003GC000568)) models of thermodynamic equilibrium in silicate systems. Formerly known as 'Adiabat_1ph', and described in a software brief in G-cubed ([Smith & Asimow, 2005](http://dx.doi.org/10.1029/2004GC000816)), it evolved to 'alphaMELTS' via Adiabat_1ph &rarr; A1ph &rarr; Alph &rarr; alpha &rarr; alphaMELTS.

alphaMELTS 2 was a complete rewrite incorporating the latest source code for Rhyolite-MELTS 1.0.2, 1.1.0 and 1.2.0 ([Gualda et al. 2012](https://doi.org/10.1093/petrology/egr080); [Ghiorso & Gualda, 2015](https://doi.org/10.1007/s00410-015-1141-8)), and the pMELTS model forked from the ENKI-portal on GitLab (see [xMELTS](https://gitlab.com/ENKI-portal/xMELTS)). The program is split into libraries of alphaMELTS functions and a text-based front-end to access them. The external interface (i.e. the menu options and files) of the alphamelts executable resembles the older version, but is more flexible and easier to use, with tab completion and command history / logging. alphaMELTS can be automated and called from the command line, or from scripts in MATLAB, Python or R.

Version 2.3.0 brings together standalone, MATLAB and Python versions that were previously posted and versioned in separate repositories. Source code for the underlying libalphaMELTS C library is available for browsing and an old Makefile is provided for reference. Workflows to build and package alphaMELTS interfaces on multiple platforms will be added soon, and comprehensive documentation and a software brief are upcoming.

How does alphaMELTS compare to other versions of MELTS?
=======================================================

THe alphaMELTS suite of software packages is intended to bridge the middle ground in terms of ease of use and power. All versions have a small footprint and are relatively easy to install locally and natively on all platforms. It is simple to automate repetative tasks compared to graphical user interfaces (including the original (Rhyolite-)MELTS GUI and easyMelts) and does not require an internet connection (unlike MELTS for Excel). It has a gentler learning curve than ENKI's Thermoengine package for those new to Python, but not the same level of power and scalability that Thermoengine can offer.

alphaMELTS forms the base for MELTS calculations in [PetThermoTools](https://github.com/gleesonm1/PetThermoTools) - an open-source Python3 tool for performing phase equilibria calculations using the MELTS family or [Holland et al. (2018)](https://doi.org/10.1093/petrology/egy048) thermodynamic models. PetThermoTools is suitable all users whether they are new to Python or advanced users and can easy integrate with other open-source Python packages, such as [Thermobar](https://github.com/PennyWieser/Thermobar) and [PySulfSat](https://github.com/PennyWieser/PySulfSat).

Will I get the same results with alphaMELTS as other MELTS software?
====================================================================

You can see the forked xMELTS code and early (alpha)MELTS for MATLAB development on GitHub. Separate Git repositories for alphaMELTS 2 and alphaMELTS for MATLAB/Python are archived on the MAGMA Source GitList server.
* https://github.com/magmasource/MAGMA
* https://magmasource.caltech.edu/gitlist/alphaMELTS2.git/
* https://magmasource.caltech.edu/gitlist/MELTS_Matlab.git/

The most noticeable difference between alphaMELTS and other MELTS software is that the feldspar phase has been separated into plagiocalse and alkali-feldspar to ease bookkeeping. All models use the pMELTS names for solid phases (i.e. no spaces in names). Some other changes include:

    Tridymite and cristobalite were tuned to match experiments on natural compositions.
    Rutile fixed (error in Tfus inherited from Samsonov), and adjusted for pMELTS.
    Corundum adjusted for pMELTS, taking account of pMELTS adjustments to alumina liquid.

More details will be listed in the GitHub Wiki and Changelog.

Installation
============

alphaMELTS has been repackaged so that users need only download the files needed for the interface (standalone app, MATLAB, or Python) and platform that they intend to use (Windows 10/11 PC, macOS Intel or Apple Silicon, or Linux PC). For now, the installation process for each interface is inherited from the separate packages / repositories and is documented in the GitHub Wiki.
