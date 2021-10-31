COMPAS code repository
======================

The public COMPAS code resides in a ``GitHub``\ [#f1]_ repository.  You will need access to ``GitHub``, and the ``git`` version 
control system installed.

If you do not have ``git`` installed, visit `Install Git <https://www.atlassian.com/git/tutorials/install-git>`__ and follow the instructions there.


Getting your first copy of the COMPAS source code
-------------------------------------------------

While you don't need a ``GitHub`` account to read the COMPAS ``GitHub`` repository, you will need an account to push your code 
changes to the repository. If you plan to contribute to the COMPAS code base, you will need to create a ``GitHub`` account.


To Fork or Not to Fork?
-----------------------

You can clone the COMPAS repository directly, or first create a fork of the repository and clone the fork. 

Forking a ``GitHub`` repository creates a copy of the repository on ``GitHub``, in your account.  The original repository, and the fork 
created, are linked on ``GitHub``.

Cloning a repository creates a copy of the repository on your local system.

See `The difference between forking and cloning a repository <https://github.community/t/the-difference-between-forking-and-cloning-a-repository/10189>`__ 
to learn more about the pros and cons of forking a ``GitHub`` repository.


If you choose to fork the COMPAS repository, use the “fork” button at the top right of the ``GitHub`` repository page.


Creating a clone of the GitHub repository
-----------------------------------------

Whether you forked the COMPAS repository or chose to work directly with the repository, you will need to clone (copy) the repository to 
your local system.

First, change to the directory on your local system within which you wish to store your local copy of the COMPAS source code
(throughout this documentation we use the example directory `~/codes`):

::

    cd ~/codes


If you have configured your ``GitHub`` account for ``SSH``, you can clone with:

::

    git clone git@github.com:USERNAME/REPONAME.git


If you have not yet configured your ``GitHub`` account with ``SSH``, you can clone over ``HTTPS``:

::

    git clone https://github.com/USERNAME/REPONAME.git


(Subsititute your GitHub username for "USERNAME", and the GitHub repository name for "REPONAME"
("COMPAS" if you did not create a fork))

At the completion of the command you will have a copy (clone) of the COMPAS ``GitHub`` repository in the `~/codes` directory on your 
local system.

You can also use the green "Code" dropdown on the ``GitHub`` repository page to clone the repository.



Setting your username and email address
---------------------------------------

Before you can push changes, you must ensure that your name and email address are set:

::

   cd ~/codes
   git config --global user.name "Fred Nurk"
   git config --global user.email "fred.nurk@aplace.adomain"


You should  now be ready to start using COMPAS!


.. rubric:: Footnotes

.. [#f1] https://github.com/
