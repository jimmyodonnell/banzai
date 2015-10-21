#!/bin/bash
# Idempotent shell script to install system level prerequisites for 
# the Banzai NGS pipeline.  Designed to be executed from Vagrantfile.
# Usage: provision.sh centos vagrant (default)
if [ "$EUID" -ne 0 ]
then echo "Please run as root"
    exit 1
fi

if [ $1 ]
then
    OS=$1
else
    OS='centos7'
fi      

if [ $2 ] 
then
    USER=$2
else
    USER='vagrant'
fi
if id -u "$USER" >/dev/null 2>&1; 
then
    echo "user $USER exists"
else
    echo "user $USER does not exist"
    exit 1
fi

cd /home/$USER
mkdir Downloads && cd Downloads

# OS specific provisioning
# TODO: Add stanza for other OSes, e.g. 'ubuntu'
if [ $OS = 'centos7' ]
then
    echo Disable SELinux
    sed -i 's/SELINUX=enforcing/SELINUX=disabled/' /etc/selinux/config
    mkdir /selinux
    echo 0 > /selinux/enforce

    echo Add epel, remi, and required system libraries and packages
    yum makecache fast
    yum -y install wget git
    wget -q -N http://dl.fedoraproject.org/pub/epel/7/x86_64/e/epel-release-7-5.noarch.rpm
    wget -q -N http://rpms.famillecollet.com/enterprise/remi-release-7.rpm
    rpm -Uvh remi-release-7*.rpm epel-release-7*.rpm

    yum install -y xz-libs gcc perl-Archive-Any perl-Digest-MD5 perl-List-MoreUtils R
    yum install -y autoconf automake python-devel python-pip zlib-devel perl-Archive-Tar

    yum -y groups install "GNOME Desktop"
fi

# Commands that work on any *nix

echo Install PEAR
wget -q -N http://sco.h-its.org/exelixis/web/software/pear/files/pear-0.9.6-bin-64.tar.gz
tar xzf pear-0.9.6-bin-64.tar.gz
cp pear-0.9.6-bin-64/pear-0.9.6-bin-64 /usr/local/bin

echo Build and install SWARM
git clone https://github.com/torognes/swarm.git
cd swarm/src/
make
cd ../../
cp swarm/bin/swarm /usr/local/bin

echo Build and install VSEARCH
git clone https://github.com/torognes/vsearch.git
cd vsearch
./autogen.sh
./configure
make
make install
cd ..
make

echo Install CUTADAPT
pip install cutadapt

echo Build and install SEQTK
git clone https://github.com/lh3/seqtk.git
cd seqtk
make
cp seqtk /usr/local/bin
cp trimadap /usr/local/bin
cd ..

echo Install BLAST+
wget -q -N ftp://ftp.ncbi.nlm.nih.gov/blast/executables/LATEST/ncbi-blast-2.2.31+-2.x86_64.rpm
rpm -ivh ncbi-blast-2.2.31+-2.x86_64.rpm

echo Install MEGAN
wget -q -N http://ab.inf.uni-tuebingen.de/data/software/megan5/download/MEGAN_unix_5_10_6.sh
unset DISPLAY
bash MEGAN_unix_5_10_6.sh -q

echo Install VEGAN and GTOOLS
wget -q -N https://cran.rstudio.com/src/contrib/gtools_3.5.0.tar.gz
cat <<EOR > install-packages.r
install.packages("vegan", repos="http://r-forge.r-project.org/")
install.packages("gtools_3.5.0.tar.gz", repos=NULL, type="source")
EOR
Rscript install-packages.r

echo Install FASTQC
pushd /opt
wget -q -N http://www.bioinformatics.babraham.ac.uk/projects/fastqc/fastqc_v0.11.3.zip
unzip fastqc_v0.11.3.zip
chmod +x /opt/FastQC/fastqc
ln -s /opt/FastQC/fastqc /usr/local/bin/fastqc
popd

# Appears to be installed by R
##echo Install RSTUDIO SERVER
##wget -q -N http://download2.rstudio.org/rstudio-server-rhel-0.99.486-x86_64.rpm
##rpm -ivh rstudio-server-rhel-0.99.486-x86_64.rpm

echo Build database for locate command
updatedb

##echo Configure and start services
##/usr/bin/systemctl enable httpd.service
##/usr/bin/systemctl start httpd.service

##echo Modifying local firewall to allow incoming connections on port 80
##firewall-cmd --zone=public --add-port=80/tcp --permanent
##firewall-cmd --reload

echo Configuring vim edit environment
cd /home/$USER
cat <<EOT > .vimrc
:set tabstop=4
:set expandtab
:set shiftwidth=4
EOT

echo Cloning Banzai repo from https://github.com/MBARIMike/banzai...
mkdir dev && cd dev
git clone https://github.com/MBARIMike/banzai banzai

echo Giving user $USER ownership of everything in /home/$USER
chown -R $USER /home/$USER

echo -------------------------------------
echo Banzai pipeline provisioning finished
echo -------------------------------------

