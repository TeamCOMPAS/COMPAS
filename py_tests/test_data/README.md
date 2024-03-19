# Test Data

This folder contains COMPAS configs used for testing. 
COMPAS is run from within the pytest environment and the output is cached for the remainder of the tests.

The pytest runs the following command:
```bash
compas_run_submit fiducial_bbh_configs.yaml
```

The pytest cache is _not_ cleared locally (but will be cleared on the CI server).
Hence, this may need to be rerun locally if changes are made to COMPAS.