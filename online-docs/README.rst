COMPAS online documentation
===========================

- COMPAS homepage: https://compas.science/docs
- Code documentation: https://compas.readthedocs.io/
- COMPAS Github: https://github.com/TeamCOMPAS/COMPAS


Tools
-----

- ReadTheDocs: https://readthedocs.org/
- Sphinx: https://www.sphinx-doc.org/en/master/index.html

Updating the documentation
--------------------------

1. Edit .rst files in 'docs/online-docs'
2. Follow existing structure for new files
3. Link new files to an index (toctree)


Building the documentation locally
----------------------------------

In the repository root directory:

.. code-block:: bash

    pip install -e '.[dev]'
    cd online-docs
    make clean
    make html




View results in 'docs/online-docs/_build/html/index.html'

Make sure there arnt any broken links! See the build logs:

```
(pages/User guide/docker: line   22) ok        https://stackoverflow.com/questions/23735149/what-is-the-difference-between-a-docker-image-and-a-container
(pages/Developer guide/Developer build/docker-developer: line   23) ok        https://www.atlassian.com/continuous-delivery/principles/continuous-integration-vs-delivery-vs-deployment
(pages/User guide/Running COMPAS/running-via-docker: line   41) broken    https://stackoverflow.com/questions/23735149/what-is-the-difference-between-a-docker-image-and-a-container#:~:text=An%20instance%20of%20an%20image,of%20layers%20as%20you%20describe.&text=You%20can%20see%20all%20your,an%20image%20is%20a%20container - Anchor '%3A~%3Atext%3DAn%20instance%20of%20an%20image%2Cof%20layers%20as%20you%20describe.%26text%3DYou%20can%20see%20all%20your%2Can%20image%20is%20a%20container' not found
(pages/User guide/docker: line   49) ok        https://www.docker.com/

```


Pushing the changes online
--------------------------

1. Push updated (.rst) source files to COMPAS repo (not the _build dir)
2. ReadTheDocs automatically rebuilds (takes up to 15 minutes)
3. For manual rebuild:
   - Log in to https://readthedocs.org/projects/compas/
   - Go to 'Overview' page
   - Click 'Build' button

Note: If build fails with "environment" error, wait and retry.

Logon details can be found on the COMPAS slack workspace (devel\_compas\_documentation channel).
