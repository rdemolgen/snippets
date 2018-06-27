#!/bon/bash
#Script to install the required dependencies for the python scripts in this folders
#If you need some programs running in Python 3.6 and others in other Python version (like Python 2.7), you may 
#take adevantages of the conda enviromments. Conda environments are sets of Python executables and modules compatible to them.
#If you need to create an environment specific for Python 3.6, you may do it with the following command:
#conda create --name python36 python=3.6
#Then you, should activate the enviromnent with the command
#source activate python36
#each subsequent time you need to use this modules, you only need to activate the module
#When you finnish, you may deaactivate this enviromnent (to use the base environment) with the command:
#source deactivate
#ot you may activate other environment (which will replace the current one)
#Pysam for the exac browser(it install samtools, htslib and the like)
conda install -c bioconda gnuplot numpy pandas pysam
#MehanicalSoup 5 install as dependency beautifulsoup4-4.6.0
conda install -c conda-forge mechanicalsoup=5
conda install pymongo flask
pip install flask-runner
pip install flask-errormail
