Python Developer Guide
=====================

Install the python post processing code with the test suite::

    pip install -e .[test]


Run the test suite::

    pytest py_tests

When adding new features, please add tests to the test suite.


Python styling
--------------

We use the `black <https://github.com/psf/black>`_ code formatter.
Please run `black .` before committing python files.



