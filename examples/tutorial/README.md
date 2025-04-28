alphaMELTS Examples
===================

This directory contains examples of the 'tutorial' 'Quick start using MELTS' for various interfaces. It is intended for comparison with other versions of the tutorial that are hosted on the MELTS Software Users forum:

* [easyMelts version of tutorial: 'Quick start using MELTS'](https://magmasource.caltech.edu/forum/index.php/topic,992.0.html)
* [other pinned versions of tutorial: 'Quick start using MELTS'](https://magmasource.caltech.edu/forum/index.php/board,12.0.html)

The MATLAB and Python versions contain comments in the files, whereas the standalone app has a README (taken from the [forum example](https://magmasource.caltech.edu/forum/index.php/topic,975.0.html)). The R version is a proof of concept. The plan is to provide alphaMELTS for Python wrappers for R and Julia that handle the reticulate and PyCall libraries respecively.

The [Generic Mapping Tools](https://www.generic-mapping-tools.org/) (GMT) program can be used to automate plotting on the command line and works well with standard alphaMELTS output. It also has Python and Julia wrappers (see [PyGMT](https://www.pygmt.org/latest/) and [GMT.jl](https://github.com/GenericMappingTools/GMT.jl)), as well as an older MATLAB interface ([GMT/MEX](https://github.com/GenericMappingTools/gmtmex)). The GMT_Plot.pdf figure should be reproducable with the standalone alphaMELTS app, the MATLAB version and Python / Jupyter notebook versions and provides a quick way to verify visually that the example has run correctly (with the caveat that the MATLAB version does not rotate / crop the figure correctly and has only been tested on Windows).

More examples can be found on the MAGMA Source GitList site: https://magmasource.caltech.edu/gitlist/Workshops.git/

These will eventually be copied to a dedicated site, to be hosted at alphamelts.caltech.edu.

