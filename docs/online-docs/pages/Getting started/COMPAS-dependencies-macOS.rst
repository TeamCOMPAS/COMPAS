macOS installation
==================

It is strongly recommended to update to the latest version of macOS through the App Store. You can find what macOS version you are 
using by clicking on the Apple symbol on the top left of your screen and clicking "About This Mac".

The next step is to install or update Xcode. You can find it directly in the App Store or at `Xcode <https://developer.apple.com/xcode/>`__\ . 
Note: Xcode installation requires around 20 GB of disk space. If you are low on disk space, you may consider installing a ``C++`` 
compiler directly.
 
Once Xcode is installed, open a Terminal, and execute the following command to install the required command line developer tools::
 
    xcode-select --install

Next, you need to install several extra libraries and python modules. Popular ways of installing them are via package managers MacPorts and Homebrew. 
We give instructions for installing ``gsl``, ``boost``, and ``hdf5`` with Homebrew. To install Homebrew, run::

    /usr/bin/ruby -e "$(curl -fsSL https://raw.githubusercontent.com/Homebrew/install/master/install)"

If the installation was successful, the following should run without error::

    brew --version

Now install ``gsl``, ``boost``, and ``hdf5`` using Homebrew by running::

    brew install gsl
    brew install boost
    brew install hdf5

