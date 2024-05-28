#!/usr/bin/env python
# -*- coding: utf-8 -*-
# Filename: tasplot.py
 
############################################################################
#
#  COPYRIGHT:  (C) 2015 John A Stevenson / @volcan01010
#                       Joaquin Cortés
#  WEBSITE: http://all-geo.org/volcan01010
#
#  This program is free software; you can redistribute it and/or modify
#  it under the terms of the GNU General Public License as published by
#  the Free Software Foundation; either version 3 of the License, or
#  (at your option) any later version.
#
#  This program is distributed in the hope that it will be useful,
#  but WITHOUT ANY WARRANTY; without even the implied warranty of
#  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#  GNU General Public License for more details.
#
#  http://www.gnu.org/licenses/gpl-3.0.html
#
#############################################################################/

__all__ = ['add_LeMaitre_fields']

class MissingModuleException(Exception):
    pass

 
# Plot LeMaitre lines
def add_LeMaitre_fields(plot_axes, fontsize=8, color=(0.6, 0.6, 0.6)):
    """Add fields for geochemical classifications from LeMaitre et al (2002)
    to pre-existing axes.  If necessary, the axes object can be retrieved via
    plt.gca() command. e.g.
    
    ax1 = plt.gca()
    add_LeMaitre_fields(ax1)
    ax1.plot(silica, total_alkalis, 'o')
    
    Fontsize and color options can be used to change from the defaults.
    
    It may be necessary to follow the command with plt.draw() to update
    the plot.
    
    Le Maitre RW (2002) Igneous rocks : IUGS classification and glossary of
        terms : recommendations of the International Union of Geological 
        Sciences Subcommission on the Systematics of igneous rocks, 2nd ed. 
        Cambridge University Press, Cambridge
"""

    # Check matplotlib is imported
    import sys
    if 'matplotlib.pyplot' not in sys.modules:
        raise MissingModuleException("""Matplotlib not imported.
        Matplotlib is installed as part of many scientific packages and is
        required to create plots.""")
    
    # Check that plot_axis can plot
    if 'plot' not in dir(plot_axes):
        raise TypeError('plot_axes is not a matplotlib axes instance.')
    
    # Prepare the field information
    from collections import namedtuple
    FieldLine = namedtuple('FieldLine', 'x1 y1 x2 y2')
    lines = (FieldLine(x1=41, y1=0, x2=41, y2=7),
             FieldLine(x1=41, y1=7, x2=52.5, y2=14),
             FieldLine(x1=45, y1=0, x2=45, y2=5),
             FieldLine(x1=41, y1=3, x2=45, y2=3),
             FieldLine(x1=45, y1=5, x2=61, y2=13.5),
             FieldLine(x1=45, y1=5, x2=52, y2=5),
             FieldLine(x1=52, y1=5, x2=69, y2=8),
             FieldLine(x1=49.4, y1=7.3, x2=52, y2=5),
             FieldLine(x1=52, y1=5, x2=52, y2=0),
             FieldLine(x1=48.4, y1=11.5, x2=53, y2=9.3),
             FieldLine(x1=53, y1=9.3, x2=57, y2=5.9),
             FieldLine(x1=57, y1=5.9, x2=57, y2=0),
             FieldLine(x1=52.5, y1=14, x2=57.6, y2=11.7),
             FieldLine(x1=57.6, y1=11.7, x2=63, y2=7),
             FieldLine(x1=63, y1=7, x2=63, y2=0),
             FieldLine(x1=69, y1=12, x2=69, y2=8),
             FieldLine(x1=45, y1=9.4, x2=49.4, y2=7.3),
             FieldLine(x1=69, y1=8, x2=77, y2=0))

    FieldName = namedtuple('FieldName', 'name x y rotation')
    names = (FieldName('Picro\nbasalt', 43, 2, 0),
             FieldName('Basalt', 48.5, 2, 0),
             FieldName('Basaltic\nandesite', 54.5, 2, 0),
             FieldName('Andesite', 60, 2, 0),
             FieldName('Dacite', 68.5, 2, 0),
             FieldName('Rhyolite', 76, 9, 0),
             FieldName('Trachyte\n(Q < 20%)\n\nTrachydacite\n(Q > 20%)',
                       64.5, 11.5, 0),
             FieldName('Basaltic\ntrachyandesite', 53, 8, -20),
             FieldName('Trachy-\nbasalt', 49, 6.2, 0),
             FieldName('Trachyandesite', 57.2, 9, 0),
             FieldName('Phonotephrite', 49, 9.6, 0),
             FieldName('Tephriphonolite', 53.0, 11.8, 0),
             FieldName('Phonolite', 57.5, 13.5, 0),
             FieldName('Tephrite\n(Ol < 10%)', 45, 8, 0),
             FieldName('Foidite', 44, 11.5, 0),
             FieldName('Basanite\n(Ol > 10%)', 43.5, 6.5, 0))

    # Plot the lines and fields
    for line in lines:
        plot_axes.plot([line.x1, line.x2], [line.y1, line.y2],
                       '-', color=color, zorder=0)
    for name in names:
        plot_axes.text(name.x, name.y, name.name, color=color, size=fontsize,
                 horizontalalignment='center', verticalalignment='top',
                 rotation=name.rotation, zorder=0)