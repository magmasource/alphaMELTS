alphaMELTS app version of tutorial: 'Quick start using MELTS'
========================================================================

This README will show you how to work through one page of the MELTS tutorial, using the new version of the standalone alphaMELTS program:

[Testing MELTS installation using the sample input file - Try this first](https://melts.ofm-research.org/Manual/UnixManHtml/examples-test.html)

The page is part of the 'Getting Started: Testing and Examples' section of the [MELTS manual](https://melts.ofm-research.org/manual.html), written by Mark Ghiorso for the original graphical-user-interface (GUI) version. The idea of this topic is to show the similarities and differences between his GUI and our text-menu versions of MELTS and how to reproduce rhyolite-MELTS GUI MELTS calculations in alphaMELTS version 2 onwards. It should be useful for newcomers to MELTS / alphaMELTS users, or for alphaMELTS users who want to use the rhyolite-MELTS and mixed H~2~O-CO~2~ fluid models, and for those who are familiar with the the rhyolite-MELTS GUI version. This was posted for the release version of alphaMELTS 2.0 but should apply to alphaMELTS 2.3.0.

This example uses an input file that has been modified from the [Morb.melts](https://melts.ofm-research.org/Melts-ftp/Morb.melts) file used in the original tutorial. The modifications simply bring it up to date with the latest GUI versions of MELTS and do not affect software results. Early versions of the MELTS software used a 'Mode: Fractionation' line in the input file to select fractional crystallization, whereas in later versions this is replaced by 'Mode: Fractionate Solids'. Also some phases were suppressed in the original tutorial because the solid solution models were still being developed / tested. None of the phases are stable for the run conditions of the tutorial example, so suppressing these phases has no effect. However, the 'Suppress: ' lines can cause read errors as the names of some phases have changed slightly since (see [here](https://magmasource.caltech.edu/forum/index.php/topic,85.0.html)).

A version of the modified Morb.melts is included here. Note that this file is identical to the one linked above, except that it has Windows file endings, but the [morb.melts](https://melts.ofm-research.org/Macosx/morb.melts) version linked on the [Mac OS X download site](https://melts.ofm-research.org/macosx.html) differs slightly in Cr~2~O~3~ content and in the *P-T* conditions.

Assuming alphaMELTS is correctly installed, open Terminal (on Mac), cmd.exe (on native Windows), the WSL or similar. Navigate to the tutorial folder by typing 'cd MELTS/tutorial'. alphaMELTS is started by typing 'run-alphamelts.pl' (note the dash, '-'; 'run_alphamelts.pl' with an underscore, '_',  is for alphaMELTS 1.X). The script will ask whether you would like to add some command line switches. Unlike alphaMELTS 1.X, the *settings_file* is optional in alphaMELTS 2+. For now we won't use one, so just press return for the default settings.

```
Mac:tutorial psmith$ run-alphamelts.pl
No command line switches specified! Please enter switches now (or press return for defaults).
```

The alphamelts program will detect that there are no settings and print the same set of questions as the GUI for choosing the thermodynamic model. There is more explanation of the various MELTS versions in the [MELTS decision tree](https://melts.ofm-research.org/MELTS-decision-tree.html). However, rhyolite-MELTS 1.0.2 should give identical results to the original MELTS model for the bulk composition used here. So, choose 'n' and press enter.

```
Checking for updates ... 2.02 ... 2.02.01
ALPHAMELTS_VERSION 2.0201


ALPHAMELTS_CALC_MODE not set!

---> Default calculation mode is rhyolite-MELTS (v. 1.0.2).  Change this? (y or n): [b]n[/b]
---> Calculation mode is rhyolite-MELTS (public release v 1.0.2).

*** alphaMELTS 2.2.1 (Aug 22 2023 14:02:19) -- rhyolite-MELTS P-T path w/ or w/o liquid ***

... etc ...

alphaMELTS main menu (enter 'x' when done):
 1. Read MELTS file to set composition and system properties
 2. Twiddle starting or continuation parameters
 3. Single (batch) calculation
 4. Execute follow path, mineral isograd or melt contour)
 5. Set fO2 buffer
 6. Set H2O (ppm) or aH2O
 7. Impose isenthalpic, isentropic or isochoric conditions
 8. Adjust solid phase settings (including fractionation)
 9. Adjust liquid or fluid settings (including melt extraction)
10. Turn phase diagram mode on / off (includes wet liquidus)
11. Turn geothermal, grid or PT-path file mode on / off
12. Source mixer / Remixer (includes AFC and flux melting)
13. Update system properties, or I/O settings
14. Write out one or more MELTS files
15. Write thermodynamic output for all phases
16. Calculate integrated melt and output file(s)
17. Fit parental melt composition with amoeba
-1. Turn off menu display for options 1-17
 X. QUIT
Your choice: 
```

Several menu items are missing even in the release version of alphaMELTS 2+. Some of these will be available soon; others will take a little longer. Choose option 1 to read in 'Morb.melts'. Note that you only need type the first letter(s) of Morb.melts, and then press the 'tab' button to complete the name:

```
Your choice: 1
MELTS filename: Morb.melts
Input file read. Waiting for command or user input.
```

In alphaMELTS 1.X, some lines of Morb.melts were ignored so only those lines that were used were echoed back. alphaMELTS 2+ uses all lines in Morb.melts, so it does not print them unless there is a problem (or you use the '-d' command switch with run-alphamelts.pl). Although the front-end program looks very different to the one in the [original tutorial](https://melts.ofm-research.org/Manual/UnixManHtml/examples-test.html), the settings for the underlying MELTS algorithm are now exactly the same as described in the 'Four things have been accomplished:' list there. You can verify that the initial pressure and temperature are set correctly with option 2:

```
Your choice: 2
Parameters to twiddle (enter 'x' to return):
 1. Initial Temperature (1200.00 C)
 2. Final Temperature (1000.00 C)
 3. Temperature Increment (-3.00 C)
 4. Initial Pressure (500.00 bars)
 5. Final Pressure (500.00 bars)
 6. Pressure Increment (-0.00 bars) 
Choose: x
```

The phase diagram mode that was a robust way to find the liquidus in alphaMELTS 1.X is not yet hooked up to the menu in this version. However, simple 'Find Liquidus' options have been built into the next step of the calculation. So go straight to option 3, the 'Single (batch) calculation' and choose the option to use the liquidus temperature:

```
Your choice: 3
Starting guess to use (enter 'x' to return):
0. Norm-calculated subsolidus assemblage
1. Superliquidus starting guess at current temperature (1200.00 C).
2. Superliquidus starting guess at liquidus temperature.
3. Superliquidus starting guess at wet liquidus temperature.
Choose: 2
<> Found the liquidus at T = 1220.31 (C).
...Checking saturation state of potential solids.
...Projecting equality constraints.
...Minimizing the thermodynamic potential.
...Checking saturation state of potential solids.
<> Stable liquid solid assemblage achieved.
Initial alphaMELTS calculation at: P 500.000000 (bars), T 1220.312500 (C)
liquid:    SiO2 TiO2 Al2O3 Fe2O3 Cr2O3  FeO  MnO  MgO  NiO  CoO  CaO Na2O  K2O P2O5  H2O 
100.350 g 48.51 1.01 17.58  0.89  0.03 7.56 0.00 9.07 0.00 0.00 12.41 2.64 0.03 0.08 0.20 
Activity of H2O = 0.0179039  Melt fraction = 1
```

Note, this uses the 'Find Liquidus' algorithm from the GUI. It is only really reliable if used before any other phases are added to the assemblage, which is why is has been incorporated into the initial calculation step. This 'Find Liquidus' routine may give incorrect results if used later in the run, and is also unsuitable for water saturated systems. For water saturated systems use the (slower) option 3, but note that this should also only be used at the start of the run.

Now choose option 4 to continue the calculation from the current conditions, incrementing the temperature by 3 &deg;C until the temperature hits the Final Temperature (1000 &deg;C). At each stage crystallized solids will be removed.

Two questions replace the functionality of the ALPHAMELTS_SAVE_ALL and ALPHAMELTS_SKIP_FAILURE environment variables in alphaMELTS 1.9. As this is the first time calling menu option 4, the answer to the first question does not matter. Choose '0' for the second question to make the run as similar to the original GUI run as possible:

```
Your choice: 4
Record output / fractionations from scratch (1) or append to previous (0)? 0
Maximum number of steps to take, or 0 for 'unlimited' (i.e. limited by other criteria);
enter negative number to limit failed steps instead (e.g. -4 to skip failure 4 times): 0
...Checking saturation state of potential solids.
...Projecting equality constraints.
...Minimizing the thermodynamic potential.
...Checking saturation state of potential solids.
<> Stable liquid solid assemblage achieved.
Initial alphaMELTS calculation at: P 500.000000 (bars), T 1220.312500 (C)
liquid:    SiO2 TiO2 Al2O3 Fe2O3 Cr2O3  FeO  MnO  MgO  NiO  CoO  CaO Na2O  K2O P2O5  H2O 
100.350 g 48.51 1.01 17.58  0.89  0.03 7.56 0.00 9.07 0.00 0.00 12.41 2.64 0.03 0.08 0.20 
Activity of H2O = N/A  Melt fraction = 1

... etc ...

...Checking saturation state of potential solids.
...Adding the solid phase apatite to the assemblage.
...Projecting equality constraints.
...Minimizing the thermodynamic potential.
...Checking saturation state of potential solids.
<> Stable liquid solid assemblage achieved.
alphaMELTS at: P 500.000000 (bars), T 1001.312500 (C)
liquid:    SiO2 TiO2 Al2O3 Fe2O3 Cr2O3  FeO  MnO  MgO  NiO  CoO  CaO Na2O  K2O P2O5  H2O 
12.5185 g 58.75 0.53  8.91  1.50  0.00 14.11 0.00 0.72 0.00 0.00 5.08 7.98 0.20 0.62 1.60 
Activity of H2O = N/A  Melt fraction = 0.981517
olivine: 0.059366 g, composition (Ca0.01Mg0.25Fe''0.75Mn0.00Co0.00Ni0.00)2SiO4
clinopyroxene: 0.025007 g, composition cpx Na0.03Ca0.81Fe''0.55Mg0.50Fe'''0.05Ti0.02Al0.16Si1.89O6
feldspar: 0.127634 g, composition K0.00Na0.59Ca0.41Al1.41Si2.59O8
spinel: 0.018452 g, composition Fe''1.57Mg0.10Fe'''0.55Al0.10Cr0.00Ti0.68O4
apatite: 0 g, composition Ca5(PO4)3OH
Minimum Temperature reached
Successful return from alphamelts...

alphaMELTS main menu (enter 'x' when done):
```

Once it has finished executing, choose option 0 to quit (as opposed to shutting the terminal window). You should find 7 new text files have been generated, which are described in the alphaMELTS 1.X documentation. (Note that in alphaMELTS 2 the 'Phase_mass_tbl.txt' and 'Phase_vol_table.txt' files are generated after the run by the run-alphamelts.pl script.)

alphaMELTS 2 has another optional line for the *melts_file* that can take one of four values:

* 'Output: txt' - alphaMELTS-style output (see the alphaMELTS documentation file)
* 'Output: tbl' - GUI-style output, including melts.out and .tbl files that can be used with [Excel-based tools](https://magmasource.caltech.edu/forum/index.php/topic,828.0.html) for tabulation and plotting (though these tools are old and not guaranteed to work in newer versions of Excel)
* 'Output: both' - alphaMELTS and GUI output
* 'Output: none' - no text output

By default the output is alphaMELTS-style only. As of November 18th, 2020, multiple instances of phases are listed as 'phase1', 'phase2' etc. for consistency with alphaMELTS for MATLAB/Python. If you want the old alphaMELTS format (i.e. 'phase_0', 'phase_1') then use 'Output: both_0', or 'Output: txt_0' instead.

The output will differ from the original slightly from the as [GUI Liquid table](https://melts.ofm-research.org/Melts-ftp/melts-liquid.tbl-MORB) gives thermodynamic quantities, such as entropy, for the liquid whereas the main output table in alphaMELTS, 'System_main_tbl.txt', gives system values instead (including a small amount of solids prior to fractionation). The most striking difference is that alphaMELTS 2+ separates feldspars into plagioclase and alkali-feldspar phases, where the GUI (and other versions, including easyMelts) have a single felspar solution phase.

For GUI-style output, there will also be slight discrepancies between the newly generated melts.out file and the [original output file](https://melts.ofm-research.org/Melts-ftp/melts.out-MORB) due to developments since the original was made. The most obvious are the addition of Mn, Co and Ni to the olivine solid solution model, though none of these components are in the Morb.melts bulk composition, and splitting of the pyroxene phase into separate clinopyroxene and orthopyroxene phases (see [here](https://magmasource.caltech.edu/forum/index.php/topic,85.0.html)). Either way, the phase compositions are essentially the same though and differences in the phase masses are a few hundredths of a gram, at most.
