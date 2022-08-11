Linux installation
==================

You will need to install the following packages (and their prerequisites) using your package manager:


                .. list-table::
                   :widths: 25 40 35
                   :header-rows: 1
                   :class: bordered
                
                   * - Package
                     - Ubuntu (apt)
                     - RHEL (yum)
                   * - g++
                     - g++
                     - gcc
                   * - GSL
                     - libgsl-dev
                     - gsl gsl-devel
                   * - Boost
                     - libboost-all-dev
                     - boost-devel
                   * - hdf5
                     - libhdf5-serial-dev
                     - hdf5-devel

For Ubuntu, type::

    sudo apt-get install g++ libboost-all-dev libgsl-dev libhdf5-serial-dev

For Fedora::

    sudo dnf install gcc boost-devel gsl gsl-devel hdf5-devel

For RHEL or CentOS::

    sudo yum install gcc boost-devel gsl gsl-devel hdf5-devel

