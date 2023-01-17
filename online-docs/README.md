COMPAS online documentation
===========================

The COMPAS online documentation is available at this URL:

    https://compas.science/docs

which redirects to the underlying readthedocs page at:

    https://compas.readthedocs.io/en/latest/index.html


The source files for the online documentation are in the 'online-docs' directory at:

    https://github.com/TeamCOMPAS/COMPAS


The online documentation is provided via readthedocs.  To learn more about readthedocs, visit:

    https://readthedocs.org/


readthedocs uses the Sphinx document generator.  To learn more about Sphinx, visit:

https://www.sphinx-doc.org/en/master/index.html


The file 'requirements.txt', in the 'docs' directory, lists the software/python modules required
to build the documentation.  If you plan to build the documentation locally (either so you can 
test your changes before you push them to the repository, or just to view the documentation locally
rather than online), you will need to install these requirements.

The requirements.txt file informs readthedocs what dependencies need to be installed in order for
readthedocs to build the online documentation. Some modules are commented ('#') in requirements.txt
because readthedocs either doesn't require them, or installs them by default - you will need to install
the commented modules to build the documentation locally.

Python is required, so install that if it is not already installed.  Then install 'sphinx', and the
modules listed in requirements.txt (including the ones commented).


Updating the documentation
--------------------------

The documentation source file are ReST (Restructured text) files - similar to markdown.  With readthedocs,
eash .rst file is compiled into a .html file, so viewable with a web browser. The documentation source
files (the .rst files) are in 'docs/online-docs' (only the 'requirements.txt' file is required by
readthedocs to be in the 'docs' directory - if that file is not in the 'docs' directory, or top-level
directory of the repo, readthedocs will not build the online documentation).

The 'online-docs' directory contains files and sub-directories that are structured to match the structured
of the online documentation. Find the .rst file you want to modify, and make the changes in your favourite
text editor. If you need to add new .rst files, follow the existing structure, and make sure you link the
files into an index (toctree) somewhere.


Building the documentation locally
----------------------------------

Once you have Sphinx and the dependencies (from 'requirements.txt') installed, navigate to the 'docs/online-docs'
directory in you local COMPAS repo, and type:

    make clean
    make html

If everything has been installed correctly this will first remove the existing .html files for the documentation
(make clean), and then recreate them (make html).  During the build process you may see some warnings that some
documents (.rst files) are not included in any toctree - that's ok, not all our pages are accessible via a table
of contents (toctree).

If the build completes successfully, there will be .html files in the 'docs/online-docs/\_build/html' directory.
To view the newly build documentation, open 'docs/online-docs/\_build/html/index.html' with your web browser
(e.g. type 'file://path-to-compas/docs/online-docs/\_build/html/index.html' into your web browser address bar, 
where 'path-to-compas' is the path to your local COMPAS repo).

Aside: the '\_build' directory is not required on the remote repo (and will only bloat the repo), so you should
add it to your .gitignore.


Pushing the changes online
--------------------------

Once you are satisfied with your changes, push the updated source files to the COMPAS repo as you would any source
changes.  If things work properly, that's all you need to do: readthedocs has webhooks that will notice the change
and automatically rebuild the online documentation.  The process (noticing the change, rebuilding the documentation,
then deploying the updated web pages) could take up to 15 minutes or so. If something goes wrong and the changes are
not noticed by readthedocs, or you just don't want to wait 15 minutes for your changes to appear on the web, you can
log into the readthedocs project and initiate a rebuild manually.

The COMPAS readthedocs project name is 'compas', and the project page is:

    https://readthedocs.org/projects/compas/

You need to login to the readthedocs project to do anything other than look at the dashboard. Once you are logged in 
you can rebuild the project. To manually start the rebuild, make sure you are on the 'Overview' page (select the
'Overview' button if you are not) and select the 'Build' button. The build may sometimes fail with an "environment"
error - the solution is to wait a few minutes and try again (readthedocs has some concurrency and timing limits).
Logon details can be found on the COMPAS slack workspace (devel\_compas\_documentation channel).

