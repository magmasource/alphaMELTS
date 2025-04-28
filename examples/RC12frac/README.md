Exercise: isobaric fractional crystallization of a back-arc basalt

This calculation can be done in easyMELTS or in alphaMELTS. These are the alphaMELTS instructions.

For this calculation, which involves H2O-CO2 fluid in a mafic system, we will use rhyoliteMELTS 1.2.0.

The input file is KM0417_RC12_frac.melts
	This is a moderately primitive back-arc basin basalt glass, discussed in Bezos et al. 2009 (10.1029/2008JB005924)
	We'll do a conventional fractional crystallization run, down temperature at constant pressure in 1 °C steps
	The input file bravely claims it will go all the way to 600 °C, but it will likely exit on failure before this
	The option "mode: fractionate solids" makes this a fractional rather than batch crystallization run
	The initial temperature is already at the liquidus, so you don't have to search for it.

Procedure:

1. run-alphamelts.pl -p output_directory to start alphamelts2
 
2. Say "y n y" to select rhyoliteMELTS 1.2.0 calibration

3. Option 1 to read file, use the filename KM0317_RC12_frac.melts

4. Execute the path (menu option 4, suboption 1 for superliquidus initial guess at given T,
		suboption 1 to start over with saving output, suboption 0 for unlimited iterations)

5. Scroll back up and look at the screen output that was generated to quickly answer some basic questions:
* What was the liquidus phase? 
* What was the composition of this mineral at the liquidus?
* Does the model predict sensible zoning of this phase as crystallization continues?
* What was the second phase to join the crystallization assemblage?
* What happens from 1182 to 1177 °C? Does this seem right to you?
* Look at the spinel compositions; how would you name the various spinels and describe their changes?
* Look at the liquid compositions; what would you call the last liquid that we see? Are all the oxides in the liquid changing in a way that comports with your intuition?

6. X to exit. Your output is automatically saved and moved to the output directory you declared.
	
7. Things to plot and examine:
* This is a good case to use Melts_Excel combine_tbl to plot up your output.
* Plot (or examine) the melt fraction vs. temperature and describe the shape of the curve.
* Plot olivine forsterite content against plagioclase anorthite content; this can be quite diagnostic of different suites of igneous rocks, especially cumulates.
* Try looking at MgO variation diagrams; do you understand why all the oxides are doing what they do?
* Plot FeO* vs. MgO in the melt; identify the changes in slope and what crystallizing assemblages they relate to.

8. Variations worth trying:
* The "thing that happens at 1182 to 1177 °C" is undesirable. Is is not the correct behavior. We will do a workaround for this in the next exercise.
* Perhaps MELTS could get past the point where this run stopped, if we let it. You could run again and give a negative number to the "how many iterations" question in menu option 4 to tell it to skip over the first (or first few, or all) failed calculations.

9. Related calculations:
* Try changing the fO2 buffer or turning it off; how is the result different?
* Try changing the initial H2O content; how is the result different?
* Try changing the pressure; how is the result different?
