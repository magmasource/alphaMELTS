# alphaMELTS@CIT

<img src="https://raw.githubusercontent.com/magmasource/alphaMELTS/main/docs/logo_final%20ellipse.png" alt="New alphaMELTS Logo" width="40%" align="right">

alphaMELTS started as a text menu-driven interface to subroutine versions of the MELTS ([Ghiorso & Sack, 1995](http://dx.doi.org/10.1007/s004100050036)), pMELTS ([Ghiorso et al., 2002](http://dx.doi.org/10.1029/2001GC000217)), and pHMELTS ([Asimow et al., 2004](http://dx.doi.org/10.1029/2003GC000568)) models of thermodynamic equilibrium in silicate systems. Formerly known as 'Adiabat_1ph', and described in a software brief in G-cubed ([Smith & Asimow, 2005](http://dx.doi.org/10.1029/2004GC000816)), it evolved to 'alphaMELTS' via Adiabat_1ph &rarr; A1ph &rarr; Alph &rarr; alpha &rarr; alphaMELTS.

alphaMELTS 2 was a complete rewrite incorporating the latest source code for Rhyolite-MELTS 1.0.2, 1.1.0 and 1.2.0 ([Gualda et al. 2012](https://doi.org/10.1093/petrology/egr080); [Ghiorso & Gualda, 2015](https://doi.org/10.1007/s00410-015-1141-8)) and the pMELTS model, as forked from the ENKI-portal on GitLab (see [xMELTS](https://gitlab.com/ENKI-portal/xMELTS)). The program is split into libraries of alphaMELTS functions and a text-based front-end to access them. The external interface (i.e. the menu options and files) of the alphamelts executable resembles the older version, but is more flexible and easier to use, with tab completion and command history / logging. alphaMELTS can be automated and called from the command line, or from scripts in MATLAB, Python or R.

Version 2.3.0 brings together standalone, MATLAB and Python versions that were previously posted and versioned in separate repositories. Source code for the underlying libalphaMELTS C library is available for browsing and an old Makefile is provided for reference. Workflows to build and package alphaMELTS interfaces on multiple platforms will be added soon, and comprehensive documentation and a software brief are upcoming.

## How does alphaMELTS compare to other versions of MELTS?

The alphaMELTS suite of software packages is intended to bridge the middle ground in terms of ease of use and power. All versions have a small footprint and are relatively easy to install locally and natively on all platforms. It is simple to automate repetative tasks compared to graphical user interfaces (including the original (Rhyolite-)MELTS GUI and easyMelts) and does not require an internet connection (unlike MELTS for Excel). It has a gentler learning curve than ENKI's [Thermoengine](https://gitlab.com/ENKI-portal/ThermoEngine) package for those new to Python, but not the same level of power and scalability that Thermoengine can offer.

alphaMELTS forms the base for MELTS calculations in [PetThermoTools](https://github.com/gleesonm1/PetThermoTools), an open-source Python3 tool for performing phase equilibria calculations using the MELTS family or [Holland et al. (2018)](https://doi.org/10.1093/petrology/egy048) thermodynamic models. PetThermoTools is suitable for all users whether they are new to scripting or advanced users of Python and can easy integrate with other open-source Python packages, such as [Thermobar](https://github.com/PennyWieser/Thermobar) and [PySulfSat](https://github.com/PennyWieser/PySulfSat).

## Will I get the same results with alphaMELTS as other MELTS software?

You can see the forked xMELTS code and early (alpha)MELTS for MATLAB development on GitHub. Separate Git repositories for alphaMELTS 2 and alphaMELTS for MATLAB/Python are archived on the MAGMA Source GitList server.
* https://github.com/magmasource/MAGMA
* https://magmasource.caltech.edu/gitlist/alphaMELTS2.git/
* https://magmasource.caltech.edu/gitlist/MELTS_Matlab.git/

The most noticeable difference between alphaMELTS and other MELTS software is that the feldspar phase has been separated into plagiocalse and alkali-feldspar to ease bookkeeping. All models use the pMELTS names for solid phases (i.e. no spaces in names). Some other changes include:

    Tridymite and cristobalite were tuned to match experiments on natural compositions.
    Rutile fixed (error in Tfus inherited from Samsonov), and adjusted for pMELTS.
    Corundum adjusted for pMELTS, taking account of pMELTS adjustments to alumina liquid.

All of these phases tended not to be stabilized when they should have been. These properties of these phases are now closer to correct, but they may need to be suppressed if recreating past simulations or performing calculations outside the narrow composition range for which there are constraints.

As of alphaMELTS 2.3.1 (released August 3, 2024) the fluid phase is always called 'fluid', regardless of which Rhyolite-MELTS or pMELTS model is being used. This differs from easyMelts and other MELTS software where the fluid phase is 'water' for models with pure-H<sub>2</sub>O fluid (Rhyolite-MELTS v1.0.2, pMELTS) and 'fluid' for the models that can include CO<sub>2</sub>.

Details of other changes will be listed in the GitHub Wiki and Changelog.

## Installation

alphaMELTS has been repackaged so that users need only download the files needed for the interface (standalone app, MATLAB, or Python) and the platform (Windows 10/11, macOS, Linux) and processor they intend to use (e.g. x86_64 for PC or Intel Mac; arm64 / aarch64 for Apple Silicon Macs). For now, the installation process for each interface is inherited from the separate packages / repositories and is documented in the [GitHub Wiki](https://github.com/magmasource/alphaMELTS/wiki). Click the drop down in the sidebar to the right to navigate to the correct set of instructions for installation and testing.

## Notifications and Support

You can get email or GitHub app notifications when a new version of alphaMELTS is released by selecting 'Watch' above, then choosing 'Custom' and ticking the 'Releases' checkbox. See [Choosing your notification settings](https://docs.github.com/en/account-and-profile/managing-subscriptions-and-notifications-on-github/setting-up-notifications/configuring-notifications#choosing-your-notification-settings) for more details on configuring the style of notification you receive.

Limited documentation is available in the 'docs' folder, and the 'tutorial' example is provided for each interface in the 'examples' folder. You can also check the [MELTS Software Users Forum](https://magmasource.caltech.edu/forum/) both in the dedicated alphaMELTS 2 / alphaMELTS for MATLAB/Python area and elsewhere, as the answers to many queries apply to more than version of MELTS. You are welcome to email the alphaMELTS developers with queries about calculations you are trying do - please find our contact details on GitHub or see the old [alphaMELTS 1.9 Support Page](https://magmasource.caltech.edu/alphamelts/support.php).

If you have found a bug or want to make a feature request or suggestion, please check the [GitHub issues](https://github.com/magmasource/alphaMELTS/issues) page and open a new issue if needed.
