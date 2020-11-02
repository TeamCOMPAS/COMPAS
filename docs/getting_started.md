[//]: ## (grip -b getting_started.md)

# Getting Started

  * [1. Installing COMPAS and Dependencies](#1-installing-compas-and-dependencies)
    + [1.1 Instructions for Linux](#11-instructions-for-linux)
    + [1.2 Instructions for macOS](#12-instructions-for-macos)
    + [1.3 Setting up the Makefile and Compiling](#13-setting-up-the-makefile-and-compiling)
    + [1.4 Installing Python](#14-installing-python)
  * [2. Evolving your first binary](#2-evolving-your-first-binary)
    + [2.1 Running COMPAS from a grid file](#21-running-compas-from-a-grid-file)
    + [2.2 Examining detailed output](#22-examining-detailed-output)
  * [3. Further queries](#3-further-queries)

## 1. Installing COMPAS and Dependencies
First change to the directory within which you wish to store your copy of COMPAS. In the rest of this document, we use as an example `~/codes`:

    cd ~/codes

Use `git clone` to download the COMPAS repository. If you do not have git installed, you may follow the instructions on https://www.atlassian.com/git/tutorials/install-git.

If you have not yet configured GitHub with SSH, you can clone over HTTPS:

    git clone https://github.com/TeamCOMPAS/COMPAS.git

With SSH configured, you can clone with

    git clone git@github.com:TeamCOMPAS/COMPAS.git

COMPAS requires a C++ compiler, and the libraries gsl and boost. We include installation instructions for Ubuntu/Linux OS and macOS.

### 1.1 Instructions for Linux
You will need to install the following packages (and their prerequisites) using your package manager:

| Package | Ubuntu (apt)     | Fedora (dnf)<br />CentOS (yum)<br />RHEL (yum) |
|---------|------------------|---------------|
| g++     | g++              | gcc           |
| Boost   | libboost-all-dev | boost-devel   |
| GSL     | libgsl-dev       | gsl gsl-devel |

So, in Ubuntu, type

    sudo apt-get install g++ libboost-all-dev libgsl-dev

In Fedora,

    sudo dnf install gcc boost-devel gsl gsl-devel

In RHEL or CentOS,

    sudo yum install gcc boost-devel gsl gsl-devel

### 1.2 Instructions for macOS
We suggest you first update to the latest version of macOS through the App Store. You can find what macOS version you are using by clicking on the Apple symbol on the top left of your screen and clicking "About This Mac".

The next step is to install or update Xcode. You can find it directly in the App Store or at https://developer.apple.com/xcode/. Note: Xcode installation requires around 20 GB of disk space. If you are low on disk space, you may consider installing a C++ compiler directly.
 
Once you have Xcode installed, open a Terminal, and run the following to install the required command line developer tools:
 
    xcode-select --install

Next, you need to install several extra libraries and python modules. Popular ways of installing them are via package managers MacPorts and Homebrew. We give instructions for installing boost and gsl with Homebrew. To install Homebrew, run

    /usr/bin/ruby -e "$(curl -fsSL https://raw.githubusercontent.com/Homebrew/install/master/install)"

 If the installation was successful, the following should run without error:

    brew --version

Now install gsl and boost using Homebrew by running

    brew install gsl
    brew install boost

### 1.3 Setting up the Makefile and Compiling
Time to actually install COMPAS. We first need to define an environment variable for the root directory of COMPAS in your shell start-up file for COMPAS to run properly. For example, if you use bash as your shell, open `~/.bashrc` with a text editor and put in the following:

    export COMPAS_ROOT_DIR=~/codes/COMPAS

where `~/codes` should be replaced with the path to the directory where you cloned the COMPAS repository. For this to take effect, either restart your bash session or run

    source ~/.bashrc

If your shell is zsh (which is the default of macOS 10.15), set the environment variable as above in `~/.zshrc` instead of `~/.bashrc`. If your shell is csh, set the environment variable in `~/.cshrc` using

    setenv COMPAS_ROOT_DIR ~/codes/COMPAS
    
Now go to the COMPAS source code directory:

    cd $COMPAS_ROOT_DIR/src

In this directory you will find the file `Makefile`, which you need to edit to point to your boost library. If you installed boost using Homebrew, find the path to the boost libraries using

    $ brew info boost
    boost: stable 1.72.0 (bottled), HEAD
    Collection of portable C++ source libraries
    https://www.boost.org/
    /usr/local/Cellar/boost/1.72.0 (14,466 files, 648.5MB) *
    ...

Copy the path, which in this case is `/usr/local/Cellar/boost/1.72.0`, and add it to the appropriate line of the Makefile:

    BOOST = /usr/local/Cellar/boost/1.72.0
 
Then, compile again by running `make -f Makefile`.

### 1.4 Installing Python
Python and some selected libraries are required for interfacing with the code, and also for post-processing. We recommend using python3. The matplotlib and numpy libraries should also be installed. The libraries scipy, astropy, and pandas are also used in some other scripts.

First check if you have python3 installed. If you do, the following should give you the version number:

    python3 --version

If you do not have python3 installed, install it by following the instructions below for your OS:

* For macOS: We recommend installing python and its libraries using MacPorts. You can follow the instructions on https://astrofrog.github.io/macports-python/
* For Linux, install `python3` using your package manager (e.g. in Ubuntu, run `sudo apt-get install python3`). We recommend installing the required python libraries using the package installer pip. E.g. To install numpy, run `pip install numpy`.


## 2. Evolving your first binary
The python script pythonSubmit.py specifies all the program options (physics assumptions, output types) and runs COMPAS in the terminal. Although the primary functionality of COMPAS is to evolve a whole population of binary stars rapidly, for now, let's focus on evolving a single stellar system and examining the detailed output.

To start, clone the default file, pythonSubmit.py in the `defaults` directory, into the `demo` directory. Then change to the `demo` directory:

    cp $COMPAS_ROOT_DIR/defaults/pythonSubmit.py $COMPAS_ROOT_DIR/docs/demo/
    cd $COMPAS_ROOT_DIR/docs/demo

### 2.1 Running COMPAS from a grid file
In population synthesis, the initial stellar population is usually generated by drawing the primary mass, secondary mass, semi-major axis, and eccentricity from their respective distributions specified in the program options. However, we illustrate COMPAS's ability to specify a grid of initial values for single and binary star evolution using COMPAS's grid functionality.

An example grid file, `Grid_demo.txt`, has been included in the `demo` directory. Open it with a text editor to view it:

    # Demo BSE Grid file
  
    Mass_1, Mass_2, Metallicity_1, Metallicity_2, Eccentricity, Separation

    63.6,   27.8,   0.001,         0.001,         0.0,          3.40


It should be clear that this grid file specifies a binary of zero-age main sequence stars with primary mass 63.6 Msol, secondary mass 27.8 Msol, primary metallicity 0.001, secondary metallicity 0.001, zero eccentricity, and semi-major axis of 3.40 AU. For more detailed documentation of COMPAS's grid functionality for both single and binary stars, please see [Specifications](./COMPAS_Doc.pdf).

We now have to tell our python submit to take its input from this grid file. To do this, open `pythonSubmit.py` with a text editor, and specify the grid filename:

    grid_filename = 'Grid_demo.txt'

To print the time evolution of binary properties, we need to turn on detailed output in the python submit:

    detailed_output = True                         # WARNING: this creates a data heavy file

Now let's run COMPAS!

    $ python3 pythonSubmit.py
    ...
    COMPAS v02.08.01

    Compact Object Mergers: Population Astrophysics and Statistics 
    by Team COMPAS (http://compas.science/index.html)
    A binary star simulator

    Start generating binaries at Mon Sep 14 15:09:14 2020

    Evolution of current binary stopped: Double compact object
    0: Evolution stopped: Black_Hole + Black_Hole

    Simulation completed

    End generating binaries at Mon Sep 14 15:09:14 2020

    Clock time = 0.168678 CPU seconds
    Wall time  = 0:0:0 (hh:mm:ss)

Congratulations! You just made a binary black hole. And it didn't even take a second.


### 2.2 Examining detailed output
The COMPAS run just now produces a new directory `COMPAS_Output`, inside which you will find the following files/directories:

* `Detailed_Output`: This directory contains `BSE_Detailed_Output_0.csv`, which records the detailed time evolution of binary. This is only produced if `detailed_output = True` in the python submit. 
* `BSE_Double_Compact_Objects.csv`: Record of binaries that formed double compact objects (in our example, we evolved one binary that became a binary black hole, which was recorded in this file.)
* `BSE_Supernovae.csv`: Record of supernovae.
* `BSE_System_Parameters.csv`: Record of all evolved binaries.
* `BSE_Common_Envelopes.csv`: Record of common envelope events.
* `Run_Details`: Record of run details.

We examine `BSE_Detailed_Output_0.csv` to look at the evolution of the two stars. A default python plotting script has been included to visualise the data. Let's run the script:

    python3 $COMPAS_ROOT_DIR/docs/demo/demo_singlesys_plotter.py

This should produce the following plot:  
![demo_plot](demo/demo_plot.png)


COMPAS provides many tools for analysing and post-processing the data. Please view the post-processing documentation in `COMPAS/postProcessing`. 


## 3. Further queries
If you have any queries unanswered by this document, your best bet is to consult our more detailed [Specifications Document](./COMPAS_Doc.pdf), which is included in your installation at `$COMPAS_ROOT_DIR/docs/COMPAS_Doc.pdf`.

If this still doesn't answer your question, you can join the [COMPAS User Google Group](https://groups.google.com/forum/#!members/compas-user) to engage in discussions with COMPAS users and developers, or email your queries to compas-user@googlegroups.com.
