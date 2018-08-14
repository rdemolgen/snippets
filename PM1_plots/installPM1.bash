#!/bin/bash
#You may install this software by downloading the compressed file via Github webpage or by its URL via wget <URL>, or using git with the command 
git clone https://github.com/Lucioric2000/snippets
cd snippets/PM1_plots

#For the first installs in the Centos server you should execute:
echo -e "[google-chrome]\\nname=google-chrome\\nbaseurl=http://dl.google.com/linux/chrome/rpm/stable/x86_64\\nenabled=1\\ngpgcheck=1\\ngpgkey=https://dl-ssl.google.com/linux/linux_signing_key.pub"|sudo tee /etc/yum.repos.d/google-chrome.repo
sudo yum -y install git nano wget bzip2 gcc libX11 libX11-devel xclock xorg-x11-drivers xorg-x11-docs xorg-x11-xinit unzip google-chrome-stable

#Installation of conda and conda packages
function conda_install(){
	#Install the Miniconda Python pachages manager
	#As is is a complex procedure, it is packed in a function. If you need to re-run this script without reinstalling Anaconda,
	#comment out the conda_install() function call several lines below.
	conda_home=/opt/conda
	python_version=3
	echo "Next, the Miniconda package will be downloaded and installed"
	wget https://repo.continuum.io/miniconda/Miniconda${python_version}-latest-Linux-x86_64.sh
	chmod +x Miniconda${python_version}-latest-Linux-x86_64.sh
	echo "You should install Miniconda to the default path there appears"
	sudo sh Miniconda${python_version}-latest-Linux-x86_64.sh -p ${conda_home}
	rm Miniconda${python_version}-latest-Linux-x86_64.sh
	#Make the updated shell path available in this session:
	export PATH="$PATH:${conda_home}"
    #Output the contents of ~/.bashrc plus the content enclosed in qoutes (which is in a string representation that handles newline characters) to the file ~/.bashrc.new.qiagen
    echo -e "\n#Shell environment for qiagen\nexport LD_LIBRARY_PATH=\$LD_LIBRARY_PATH:${srv_qiagen}/bin/ssw/src/"|cat ~/.bashrc ->~/.bashrc.new.qiagen
    #Move the file ~/.bashrc.new.qiagen to ~/.bashrc (overwriting the existent ~/.bashrc withouk asking for confirmation)
    mv -f ~/.bashrc.new.qiagen ~/.bashrc
    #Make the updated shell path available in this session:
    source ~/.bashrc
}
condabin=`which conda`
if [ -z $condabin ]
then
	conda_home=/opt/conda
	conda_install
else 
    conda_home=${condabin%/bin/conda}
    echo "Conda installation found at $conda_home. Script will use tht installation."
fi
source activate base
udo ${conda_home}/bin/conda install -c bioconda numpy pandas pysam
sudo ${conda_home}/bin/conda install -c conda-forge mechanicalsoup selenium
sudo ${conda_home}/bin/conda install pymongo flask cython lxml
sudo ${conda_home}/bin/pip install flask-runner flask-errormail
#Gnuplot installation
wget https://cytranet.dl.sourceforge.net/project/gnuplot/gnuplot/5.2.4/gnuplot-5.2.4.tar.gz
tar -xvzf gnuplot-5.2.4.tar.gz
cd gnuplot-5.2.4
./configure
make
make check
sudo make install
cd ..
sudo rm -rf gnuplot-5.2.4*

#Chrome driver (For Selenium)
wget https://chromedriver.storage.googleapis.com/2.40/chromedriver_linux64.zip
unzip chromedriver_linux64.zip
sudo mv chromedriver /usr/local/bin

#echo "Now you are in the folder `pwd`"

echo "You may execute the PM1_plotter script using, for example, the code 'python PM1_plotter.py ABCC8 123' (for the ABCC1 gene) on the subdirectory snippets/PM1_plotter."
