# snippets
## Small snippets of useful code

# exac_coverage_20x.py

The following method returns the number of individuals covered at 20x for the genomic coordinated submitted by the user
from the Exac API:

http://exac.hms.harvard.edu/

This is best utilised when a variant is not reported in Exac.

prerequisites: 
* python 3
* sys, requests

## usage:

    python exac_coverage_20x.py [Chromosome] [genomic coordinate]

example usage:

    python exac_coverage_20x.py 1 123456

To elaborate, the absence of a variant in Exac implies it is not identified in any of the cohort i.e. 60706 individuals.
However, not all regions are covered equally and therefore reporting based on individuals covered at a certain depth enables
a more informative figure.

As a caveat to this method, the allele number reported in the Exac browser includes those from individuals with a depth of coverage (DP) >= 10 
and a genotype quality >= 20. Because there is not a genotype quality score available through this method (as it is queries the depth metrics only)
this method will provide a greater allele number (number of individuals x 2) than the Exac browser for any variant in the Exac database.

# PM1_plots

This code is designed to aggregate data to contribute to the interpretation of novel variants and deduce whether said variant is in a pathogenic variation hotspot, with an absence of benign variation.

It will plot and svg using gnuplot with the amino acid residue along the x-axis. Data sources include:
* Harvard Exac API 
* HGMD
* UniProt 

If consurf data are present they will also be plotted. 

Prerequisites:
* Developed with python 3.6
* Gnuplot version 5
* Mechanicalsoup
* bs4 (beautiful soup 4)
* numpy
* pandas
* collections

## usage 

    python PM1_plotter.py [gene] [amino acid residue number]

## example usage

    python PM1_plotter.py ABCC8 123

