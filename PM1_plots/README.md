# PM1_plots

### Update 26/01/18 ###
 - Plots of ExAC data will now include all allele frequencies seperated by zygosity (heterozygous, homozygous, hemizygous).
 - When entering a password for HGMD log-in, this will not appear on the screen.
 - The amino acid residue number will appear in the X-axis of the plot.

## Introduction ##
This code is designed to aggregate data to contribute to the interpretation of novel variants and deduce whether said variant is in a pathogenic variation hotspot, with an absence of benign variation.
It was written as part of a trainee project and although it produces a useful plot, there is much scope for development in terms of code structure,
the data included as well as visualisation.

It will plot and svg using gnuplot with the amino acid residue along the x-axis and varying y axis scales. Data sources include:
* Harvard Exac API 
* HGMD
* UniProt 
* Consurf

If Consurf data are present within the 'consurf_grades' directory they will also be plotted. 
N.B. Consurf data should be saved as output by Consurf ("consurf.grades") in a directory  named as the user input gene name.
e.g. if you are trying to plot ABCC8, call the directory containing the ABCC8 consurf data "ABCC8".

You can generate Consurf conservation scores here:

http://consurf.tau.ac.il/2016/

## Installing ##
For CentOs users, there is a script (install.bash) that you may run to install the needed dependencies basing in the Conda package, 
dependency and enviromnent manager (https://conda.io). Alternately, as a standalone installer, you should use the clone_and_install.bash script,
which will clone the full snippets project, change to this directory and install the required dependencies.

## Prerequisites: ##
* Developed with python 3.6
* Gnuplot version 5
* Mechanicalsoup
* bs4 (Beautifulsoup 4)
* numpy
* pandas
* collections
# Chrome driver

## Usage ##

    python PM1_plotter.py [gene] [amino acid residue number]

## Example usage ##

    python PM1_plotter.py ABCC8 123

## Caveats ##

Please be aware that there are some caveats to bear in mind when using PM1 plots regarding transcripts.

The UniProt section of plot uses the canonical transcript as identified by UniProt, given the following criteria defined by Uniprot:
 - It is the most prevalent.
 - It is the most similar to orthologous sequences found in other species.
 - By virtue of its length or amino acid composition, it allows the clearest description of domains, isoforms, polymorphisms, post-translational modifications, etc.
 - In the absence of any information, we choose the longest sequence.

The HGMD plot uses the transcript identified by a RefSeq number on their gene page.

The ExAC plot uses the canonical Ensembl transcript, which Ensembl set according to the following hierarchy: 
  - (1) Longest CCDS translation with no stop codons.
  - (2) If no (1), choose the longest Ensembl/Havana merged translation with no stop codons.
  - (3) If no (2), choose the longest translation with no stop codons.
  - (4) If no translation, choose the longest non-protein-coding transcript. 
  
So it is possible that the transcripts used to create the plots may be different. For example, when  CASR is run, it retrieves Ensembl transcript ENST00000498619, which is used by ExAC and is 1088 amino acids long. This transcript maps to NM_001178065 rather than NM_000388, which is the transcript that HGMD uses. The Uniprot transcript is 1078 amino acids long (the same length as NM_000388).

The transcripts used for each section of the graph will be printed to screen.
**The graph plotted will be based on the length of the canonical transcript as defined by Uniprot.**

Updated and maintained by Verity Fryer (verity.fryer@nhs.net)
Created by Stuart Cannon (s.cannon@exeter.ac.uk)
