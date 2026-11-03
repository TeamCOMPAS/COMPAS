Merger validation evidence

[compas_merger_validation_2026-11-03.json](compas_merger_validation_2026-11-03.json) contains the 3 November baseline comparison, its run settings, six binary configurations, and the focused C++ regression counts. The implementation was subsequently committed as `24b5161283543d1065ca18c7b1c20bec481c3b31`; its earlier patch digest is retained as historical provenance.

The `release_verification` section records fresh checks of `03.55.02`. The tag identifies the release commit; use that tag and commit when citing a rerun. These small checks validate implementation consistency and control tracks, not the magnetar population fit. Existing paper results remain results of their originally cited builds until the populations are rerun.

The [log definitions](merger_logdefs.txt) add the hydrogen and magnetic-field columns used in the comparisons. Runs use HDF5 detailed output with `--logfile-detailed-output-record-types -1`; masses, periods, seed and other options are recorded in the JSON. The six cases compare MIST and LINEAR mergers, an explicit BRCEK core model, and a wide non-merging system.

The scientific build tag is `03.55.02` without a leading `v`: this repository automatically publishes `v*` tags to PyPI. Fetch it with `git fetch origin tag 03.55.02` and check out the tag for a reproducible rerun.
