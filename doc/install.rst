.. index:: installation

.. highlight:: bash

.. _installation:

Installation 
============

TRIQS and its applications are provided *a la carte*:
after you have installed the TRIQS library, you will be able to easily install
various TRIQS-based applications: impurity solvers, tools for DFT+DMFT
calculations, etc.

This page describes the installation of the TRIQS library itself. The
installation procedure of the applications is described on their respective
sites, under 'Install'.

Philosophy
----------

The TRIQS project is in perpetual dynamic evolution such that we can get our
developments directly to users straight out of the oven.

However, we also understand that some users may not wish to constantly update
their codes, and are happy to use an older but perhaps more stable version. 

To this end, we propose two options to the user:

#. You follow the master branch in the git repository of TRIQS and all applications.
   This will guarantee that you are using the latest stable release including essential
   bug-fixes. Note that new releases might occasionally include changes of the API, which
   we summarize in a set of release notes.
   We use continuous integration (`travis <https://travis-ci.org/TRIQS/triqs/branches>`_ and `jenkins <https://triqs.ipht.cnrs.fr/jenkins/view/TRIQS>`_) 
   to ensure that the master branch always compiles and passes all tests. This is checked
   for both the TRIQS library and several public (and private) applications.

#. You use a version tag, e.g. version 1.4, for TRIQS and all applications.
   This guarantees complete reproducibility, while you might be missing out on the latest
   features and bug-fixes.

Prerequisites
-------------

The TRIQS library relies on a certain number of standard libraries and tools
described in the :ref:`list of requirements <requirements>`. Please pay
particular attention to the :ref:`C++ compilers<require_cxx_compilers>` and to
:ref:`Python virtual environments<python_virtualenv>`.  Here are instructions to
install these necessary libraries on two standard systems:

.. toctree::
   :maxdepth: 1

   installation/ubuntu
   installation/osx_install

.. note:: Please ensure that you have all the dependencies installed before proceeding further!

Installation steps
------------------

You need to install first Cpp2Py and then TRIQS.

We provide hereafter the build instructions in the form of a documented bash script.
You can adapt INSTALL_PREFIX, NCORES for your local settings.
Note that, contrary to previous versions of TRIQS,  
the installation directory CMAKE_INSTALL_PREFIX is now mandatory in the cmake configuration.


.. literalinclude:: install.sh
   :language: bash


Environment setup
^^^^^^^^^^^^^^^^^^^

Cpp2Py and TRIQS both provide a small script (`cpp2pyvars.sh` and `triqsvars.sh`)
to load their respective installation into your :ref:`environment variables <environment_vars>`.
Please source them with the proper replacement of INSTALL_PREFIX::

        source INSTALL_PREFIX/share/cpp2pyvars.sh
        source INSTALL_PREFIX/share/triqsvars.sh

To automate this process, please add these two lines to your `~/.bash_profile <https://en.wikipedia.org/wiki/Bash_(Unix_shell)#Startup_scripts>`_
(or `~/.zprofile <http://zsh.sourceforge.net/FAQ/zshfaq03.html#l19>`_) 

Further reading
------------------
.. toctree::
   :maxdepth: 1

   installation/install_options
   installation/environment_vars
   installation/python
   installation/clang
