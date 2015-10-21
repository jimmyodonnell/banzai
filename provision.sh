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

    echo Add epel, remi, and postgres repositories
    yum makecache fast
    yum -y install wget git
    wget -q -N http://dl.fedoraproject.org/pub/epel/7/x86_64/e/epel-release-7-5.noarch.rpm
    wget -q -N http://rpms.famillecollet.com/enterprise/remi-release-7.rpm
    rpm -Uvh remi-release-7*.rpm epel-release-7*.rpm


    yum -y groups install "GNOME Desktop"
fi

# Commands that work on any *nix


echo Build and install gdal
wget -q -N http://download.osgeo.org/gdal/2.0.0/gdal-2.0.0.tar.gz        
tar xzf gdal-2.0.0.tar.gz
cd gdal-2.0.0
export PATH=$(pwd):$PATH
./configure --with-python
gmake && gmake install
cd ..


# Required to install the netCDF4 python module
echo "Need to sudo to install hdf5 packages..."
sudo yum -y install hdf5 hdf5-devel
if [ $? -ne 0 ] ; then
    echo "Exiting $0"
    exit 1
fi

# Required to install the netCDF4 python module
wget ftp://ftp.unidata.ucar.edu/pub/netcdf/netcdf-4.3.3.tar.gz
tar -xzf netcdf-4.3.3.tar.gz
cd netcdf-4.3.3
./configure --enable-hl --enable-shared
make; sudo make install
cd ..


echo Build database for locate command
updatedb

echo Configure and start services
/usr/bin/systemctl enable httpd.service
/usr/bin/systemctl start httpd.service

echo Modifying local firewall to allow incoming connections on port 80
firewall-cmd --zone=public --add-port=80/tcp --permanent
firewall-cmd --reload

echo Configuring vim edit environment
cd /home/$USER
cat <<EOT > .vimrc
:set tabstop=4
:set expandtab
:set shiftwidth=4
EOT

echo Cloning Banzai repo from https://github.com/MBARIMike/banzai...
mkdir dev && cd dev
git clone https://github.com/MBARIMike/banzai banzaigit

echo Giving user $USER ownership of everything in /home/$USER
chown -R $USER /home/$USER

echo Banzai provisioning finished
echo ----------------------------

