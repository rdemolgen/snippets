# PM1_plots

created by Stuart Cannon s.cannon@exeter.ac.uk

This code is designed to aggregate data to contribute to the interpretation of novel variants and deduce whether said variant is in a pathogenic variation hotspot, with an absence of benign variation.
It was written as part of a trainee project and although it produces a useful plot, there is much scope for development in terms of code structure,
the data included as well as visualisation.

It will plot and svg using gnuplot with the amino acid residue along the x-axis and varying y axis scales. Data sources include:
* Harvard Exac API 
* HGMD
* UniProt 

If consurf data are present within the 'consurf_grades' directory they will also be plotted. 
note they should be saved as is output by consurf ("consurf.grades") in a directory  named as the user input gene name.
e.g. if you are trying to plot ABCC8, call the directory containing the ABCC8 consurf data "ABCC8".

You can generate consurf conservation scores here:

http://consurf.tau.ac.il/2016/

Prerequisites:
* Developed with python 3.6
* Gnuplot version 5
* Mechanicalsoup
* bs4 (Beautifulsoup 4)
* numpy
* pandas
* collections

## usage 

    python PM1_plotter.py [gene] [amino acid residue number]

## example usage

    python PM1_plotter.py ABCC8 123

