#FROM ubuntu:18.04
#FROM ubuntu:14.04
FROM centos:7
#LABEL Remarks=”This is a dockerfile for snippets and PM1_plotter" 
MAINTAINER Lucio Montero <Lucioric2000@hotmail.com>
USER root
#--centos settings
RUN echo -e "[google-chrome]\\nname=google-chrome\\nbaseurl=http://dl.google.com/linux/chrome/rpm/stable/x86_64\\nenabled=1\\ngpgcheck=1\\ngpgkey=https://dl-ssl.google.com/linux/linux_signing_key.pub"|tee /etc/yum.repos.d/google-chrome.repo
RUN yum -y update && yum -y install git nano wget bzip2 gcc libX11 libX11-devel xclock xorg-x11-drivers xorg-x11-docs xorg-x11-xinit zip unzip google-chrome-stable
#--Ubuntu settings
#RUN  apt-get -m update && apt-get -y install git nano wget bzip2 gcc x11-apps unzip zip unzip gnupg
#RUN echo "deb [arch=amd64] http://dl.google.com/linux/chrome/deb/ stable main"|tee  /etc/apt/sources.list.d/google-chrome.list
#RUN wget https://dl.google.com/linux/linux_signing_key.pub
#RUN apt-key add linux_signing_key.pub
#RUN  apt-get -y install google-chrome-stable
#FROM dockerfile/chrome
RUN git clone https://github.com/Lucioric2000/snippets
FROM continuumio/miniconda3:latest
RUN conda install -y -c bioconda numpy pandas pysam
RUN conda install -y -c conda-forge mechanicalsoup selenium
RUN conda install -y pymongo flask cython lxml
RUN pip install --upgrade pip
RUN pip install flask-runner flask-errormail
#RUN wget https://chromedriver.storage.googleapis.com/2.40/chromedriver_linux64.zip -P /usr/bin 
#RUN unzip /usr/bin/chromedriver_linux64.zip
COPY chromedriver /usr/bin

# by default /bin/bash is executed
#CMD ["/bin/bash"]
