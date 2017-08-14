# snippets

created by Stuart Cannon s.cannon@exeter.ac.uk

## This repo has been created to contain scripts which have been found useful for a clinical diagnostic setting.
Each directory contains the code required for a specific purpose; the details of which are included in separated READMEs.

To utilise these tools you can pull the entire repo, meet the dependencies and then use methods as needed. 
alternatively, if you have a git version > 1.7.0 you can use the sparse checkout feature:

    mkdir <dir_name>
    cd <dir_name>
    git init
    git remote add -f origin <url>

Now you need to configure which directories to pull from the remote

    git config core.sparseCheckout true
    echo "target/dir_or_file/" >> .git/info/sparse-checkout
    echo "other/target/dir_or_file/" >> .git/info/sparse-checkout

And finally, pull the target to your machine

    git pull origin master

# Exac_coverage

The following method returns the number of individuals covered at 20x for the genomic coordinated submitted by the user
from the Exac REST API:

http://exac.hms.harvard.edu/

This is best utilised when a variant is not reported in Exac and you are interested in how many individuals it was not reported in.

# PM1_plots

This code is designed to aggregate data to provide a visual aid, and therefore piece of evidence, to support or refute the presence of a mutational hotspot.
This is in line with the ACGS variant interpretation guidlines description of point PM1. 

Note - It was written as part of a trainee project and although it produces a useful plot, there is much scope for development in terms of code structure,
the data included as well as visualisation.

