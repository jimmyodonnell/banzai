Building a Banzai development system 
====================================

First, install [Vagrant](https://www.vagrantup.com/) and and [VirtualBox](VirtualBox.md)
&mdash; there are standard installers for Mac, Windows, and Linux. Then create an empty folder 
off your home directory such as `Vagrants/banzaivm`, open a command prompt window, cd to that 
folder, and enter these commands:

    curl "https://raw.githubusercontent.com/jimmyodonnell/banzai/master/Vagrantfile.sh" -o Vagrantfile
    curl "https://raw.githubusercontent.com/jimmyodonnell/banzai/master/provision.sh" -o provision.sh
    vagrant up --provider virtualbox

The `vagrant up` command takes about an hour to build and install all the required software
needed to execute the Banzai pipeline.

All connections to this virtual machine are performed from the the directory you installed 
it in; you must cd to it (e.g. `cd ~/Vagrants/banzaivm`) before logging in with the 
`vagrant ssh -- -X` command.  After installation finishes log into your new virtual machine 
and test it:

    vagrant ssh -- -X
