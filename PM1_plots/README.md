# PM1_plots

Created by Stuart Cannon s.cannon@exeter.ac.uk

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

## Prerequisites:
* Developed with python 3.6
* Gnuplot version 5
* Mechanicalsoup
* bs4 (Beautifulsoup 4)
* numpy
* pandas
* collections

## Usage 

    python PM1_plotter.py [gene] [amino acid residue number]

## Example usage

    python PM1_plotter.py ABCC8 123

## Caveats


Please be aware that there are some caveats to bear in mind when using PM1 plots regarding transcripts.

The UniProt section of plot uses the canonical transcript as identified by UniProt, given the following criteria defined by Uniprot:
 - It is the most prevalent.
  - It is the most similar to orthologous sequences found in other species.By virtue of its length or amino acid composition, it allows the clearest description of domains, isoforms, polymorphisms, post-translational modifications, etc.
 - In the absence of any information, we choose the longest sequence.

The HGMD plot uses the transcript identified by a RefSeq number on their gene page.

The ExAC plot uses the canonical Ensembl transcript, which Ensembl set according to the following hierarchy: 
  (1) Longest CCDS translation with no stop codons.
  (2) If no (1), choose the longest Ensembl/Havana merged translation with no stop codons.
  (3) If no (2), choose the longest translation with no stop codons.
  (4) If no translation, choose the longest non-protein-coding transcript. 
  
So it is possible that the transcripts used to create the plots may be different. For example, when  CASR is run, it
retrieves Ensembl transcript ENST00000498619, which is used by ExAC and is 1088 amino acids long. This transcript maps to NM_001178065 rather than NM_000388, which is the transcript that HGMD uses. The Uniprot transcript is 1078 amino acids long (the
same length as NM_000388).

**The graph plotted will be based on the length of the canonical transcript as defined by Uniprot.**
