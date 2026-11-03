Focused COMPAS merger regressions

Run the tests with the normal compiler and GSL/Boost/HDF5 settings used for COMPAS:

```sh
make -C src test-main-sequence-mergers ODIR=/tmp/compas-merger-test-obj BDIR=/tmp/compas-merger-test-bin
```

Choose a fresh object directory when headers have changed. The test target does not update the `src/COMPAS` compatibility link. Supply the usual `GSLINCDIR`, `GSLLIBDIR`, `BOOSTINCDIR`, `BOOSTLIBDIR`, `HDF5INCDIR`, and `HDF5LIBDIR` overrides if needed on your machine.

The test executable links the production stellar classes and checks fully mixed remnant composition, hydrogen mass, age bounds, zero-time recalculation, and subsequent burning. It covers ZAMS masses 0.7–60 solar masses, including the MIST radiative/convective interpolation, and fresh/hydrogen-free limits. It also checks the 20-solar-mass review example and the separate CHE age mapping. Six runs cover LINEAR and MIST, MANDEL/HURLEY/BRCEK cores, and main-sequence magnetic-flux decay.

These tests establish consistency of the implemented prescription. They do not calibrate the MIST extrapolation at high masses, validate full merger mixing against detailed stellar models, or reproduce the magnetar population fit.

The [recorded validation results](validation/compas_merger_validation_2026-11-03.json) and [provenance notes](validation/README.md) are included in this repository for review with release `03.55.02`.
