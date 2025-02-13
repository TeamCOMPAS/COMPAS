Git Workflow
============

If you are not familiar with ``git``, please refer to a tutorial such as `this one <https://github.com/git-guides>`_.

In ``COMPAS`` we work in the ``dev`` branch, and use the ``master`` branch for the major version changes. The main workflow is as follows:

1. Create a git-issue describing the tasks you want to work on (this is optional, but recommended).
2. Create a new branch from ``dev`` with a descriptive name (e.g. ``feature/issue-123``), and make a draft pull request to ``dev``. This will allow you to work on the code collaboratively and get feedback from the team.
3. Make your changes, commit them, and push them to the remote repository. Every commit will trigger the continuous integration (CI) tests.
4. Once you are happy with your changes, check that CI tests are passing, and switch the PR to ``ready for review``. This will make it clear to the team that your changes are ready for review.

The team will review your changes, and may ask you to make some modifications. You can make these changes in the same branch, and push them to the remote repository. Once the changes are approved, the PR will be merged into ``dev``.

CI Tests
--------

There are a few tests that are run automatically when you push your changes to the remote repository. These tests are:

1. `spell-checking <https://github.com/TeamCOMPAS/COMPAS/blob/dev/.github/workflows/precommit-checks.yml>`_
   This ensures that docstrings and comments are correctly spelled.
2. `COMPAS compile test <https://github.com/TeamCOMPAS/COMPAS/blob/dev/.github/workflows/compas-compile-ci.yml>`_
   This ensures that COMPAS C++ and python utilities can be correctly compiled (and COMPAS can run on a fiducial binary system).
3. `COMPAS py-utils unit tests <https://github.com/TeamCOMPAS/COMPAS/blob/dev/.github/workflows/compas-compile-ci.yml>`_
   This ensures that some of the python utilities are working as expected.


The tests will fail if the fiducial binary system does not lead to a binary black hole merger.

.. literalinclude:: ../../../py_tests/test_data/run.sh
   :linenos:
   :language: bash