# Research version with circumbinary-disk (CBD) evolution

This fork contains the CBD implementation used in Grichener et al. (2026).
The CBD implementation is a research extension of [COMPAS](https://github.com/TeamCOMPAS/COMPAS)
and is not currently part of the official COMPAS release, and is in processes of being integrated into COMPAS. 

## Summary
This PR adds an optional post-common-envelope circumbinary disk (CBD) framework to COMPAS. The new framework checks whether a CBD can form in post-common-envelope systems containing a neutron star or black hole companion, based on angular-momentum considerations. If a CBD forms, the framework evolves the post-CEE binary under disk-driven accretion and angular-momentum exchange based on an external library called [spindler-c](https://github.com/ruggero-valli/spindler-c). When the compact object accretion rate exceeds the Eddington limit, the framework applies an Eddington-limited accretion correction, reducing the actually accreted mass and ejecting the difference as winds, which carry away angular momentum. 

The framework is disabled by default, so existing COMPAS behavior is preserved unless the new CBD options are explicitly enabled by the user supplied boolean option --circumbinary-disk.

I used ChatGPT as a coding assistant while translating my research prototype into main-branch-worthy (hopefully) C++ code. I've continuously reviewed the implementation and tested it against my prototype. 

## Motivation
Post-common-envelope CBDs affect the orbital separations, eccentricities, and masses of binary systems, and therefore their survival probabilities, merger times, and observable properties. This PR provides a first implementation of such a framework in COMPAS.

## Main changes
- Adds a new 'CircumbinaryDisk' module.
- Adds user-facing options controlling CBD formation and evolution.
- Adds CBD-related logging quantities, outputs to BSE_Circumbinary_Disks.
- Adds documentation of the physical assumptions, units, and limitations.

## Dependencies
Using the CBD module requires installing [spindler-c](https://github.com/ruggero-valli/spindler-c), which itself requires [libinterpolate](https://gitlab.com/rob.izzard/librinterpolate). While working on the CBD module I encountered a bug in spindler-c. I've opened a PR to fix it, but for now using [my version](https://github.com/AldanaGrichener/spindler-c) would be better.
By default, COMPAS will compile without spindler-c, making the CBD module unavailable. COMPAS users who want to use the CBD module will have to supply paths to spindler-c at compile time. 

## Modified files
- Adds: CircumbinaryDiskFormation.cpp, CircumbinaryDiskEvolutionSpindler.cpp, CircumbinaryDiskEvolutionUnavailable.cpp, and CircumbinaryDisk.h. 
- Modifies: BaseBinaryStar.cpp, BaseBinaryStar.h, Log.cpp, Log.h, LogTypedefs.h, Makefile, Options.cpp, Options.h, typedefs.h, yaml.h. 

## CBD user-supplied options and their defaults 
I've added quite a few user-supplied options as part of the CBD framework. Some of these should perhaps be hardcoded, while others are legacy options which ultimately shouldn't be part of the main branch. Here's a table summarizing the added options:
 
| Option | Default | Allowed values / range | Symbol | Purpose |
|---|---:|---|---|---|
| `--circumbinary-disk` | `false` | Boolean | --- | Enables post-CEE circumbinary disk formation and evolution. |
| `--circumbinary-disk-envelope-mass-fraction-supplied` | `0.05` | `[0, 0.2]`; also requires $\beta f_M <= 1$ | $f_M$ | Fraction of the common envelope mass available for CBD formation. |
| `--circumbinary-disk-beta` | `3.0` | `>= 0`; also requires $\beta f_M <= 1$ | $\beta$ | Ratio between fraction of binary's AM and fraction of envelope's mass available for CBD formation. |
| `--circumbinary-disk-lifetime` | `1.0e5 yr` | `> 0 yr` | $\tau_{CBD}$ | Duration over which the supplied CBD mass is delivered. |
| `--circumbinary-disk-structure-factor` | `4.0` | `> 1` | $\zeta$ | Disk structure factor multiplying the CBD angular-momentum threshold. |
| `--circumbinary-disk-inner-radius-over-separation` | `2.5` | `>= 1` | $R_{CBD,in}$ | Inner radius of the CBD in units of the post-CE binary separation. |
| `--circumbinary-disk-evolution-mode` | `INSTANTANEOUS` | `INSTANTANEOUS`, `TIMESTEPPED` | --- | Selects whether the CBD-binary interaction is evolved immediately after CE or over later COMPAS timesteps. |
| `--circumbinary-disk-eddington-mode` | `GE23` | `CAP`, `GE23`, `TM23` | --- | Selects the Eddington-limited CBD accretion prescription. |
| `--circumbinary-disk-eddington-cap-factor` | `1.0` | `> 0` | $c_{Edd}$ | Multiplier applied to the Eddington accretion rate in `CAP` mode. |
| `--circumbinary-disk-eddington-ge23-factor` | `2.0 / 3.0` | `(0, 1]` | $f_{GE23}$ | Retained fraction of the requested super-Eddington accretion rate in `GE23` mode. |
| `--post-common-envelope-eccentricity-prescription` | `CIRCULAR` | `CIRCULAR`, `PRE_RLOF_CAP`, `PRE_RLOF_FRACTION` | --- | Controls the eccentricity assigned immediately after CE. |
| `--post-common-envelope-eccentricity-cap` | `0.1` | `[0, 1)` | $e_{CE,max}$ | Used by `PRE_RLOF_CAP`: $e_{postCE} = min(e_{preRLOF}, e_{CE,max})$ |
| `--post-common-envelope-eccentricity-fraction` | `0.1` | `[0, 1]` | $f_e$ | Used by `PRE_RLOF_FRACTION`: $e_{postCE} = f_e e_{preRLOF}$. |
| `--circumbinary-disk-use-synchronized-envelope-angular-momentum` | `false` | Boolean | --- | Legacy option for compatibility with [Unger+2026](https://ui.adsabs.harvard.edu/abs/2026ApJ...998...79U/abstract). Uses the synchronized-envelope angular-momentum estimate in the CBD formation criterion instead of COMPAS total angular momentum. |

## CBD logic, options, and assumptions

CBD formation: we assume a CBD of mass $f_M M_{env}$ is formed whenever the system could retain enough angular momentum in the CE process to form one. The AM available for CBD formation is parametrized as $J_{avail} = \beta f_M J_{tot}$ where $J_{tot}$ is the binary's pre-RLOF total AM. The AM required to support a CBD with an inner radius $R_{CBD,in}$ is given by the Keplerian AM times the structure factor $\zeta$.

Evolution of the orbit and mass of the binary given a CBD: if a CBD is formed, the binary's total mass, mass ratio, semi-major axis, and eccentricity all evolve according to derivatives provided by spindler-c. The combination of total mass and mass ratio dictates the mass supplied to each binary component. A non-compact star accretes all the supplied mass, while the accretion towards a neutron star or a black hole might be suppressed in case of super Eddington accretion, which isn't accounted for in spindler-c. We thus add to spindler-c's derivatives terms which account for the Eddington limit. We allow three different suppression modes in case of super-Eddington accretion:


1. GE23 – $`\dot{M}_{actual} = f_{GE23}\dot{M}_{supplied}`$
Allows for strongly super-Eddington accretion. Based on [Ghodla & Eldridge 2023](https://ui.adsabs.harvard.edu/abs/2023MNRAS.523.1711G/abstract). 
2. TM23 –  $`\dot{M}_{actual} \propto \dot{M}_{supplied} ^ {\frac{2}{5}}`$
Allows for mildly super-Eddington accretion. Based on [Tuna & Metzger 2023.](https://ui.adsabs.harvard.edu/abs/2023ApJ...955..125T/abstract) 
3. CAP – $`\dot{M}_{actual} = min(\dot{M}_{supplied}, c_{Edd} \dot{M}_{Edd})`$
Ad hoc option.
 
In all cases, the difference between the supplied and retained masses is assumed to be launched in winds from the accretor, carrying away its specific orbital angular momentum.

CBD lifetime: the supplied CBD mass is delivered over a user-specified lifetime at a constant rate. By default, we assume that the lifetime is shorter than any other post-CEE timescale, and therefore the entire CBD-binary interaction can be carried out "instantaneously" in the aftermath of the CEE. I've also included a "timestepped" mode, where the CBD-binary interaction is concurrent with whatever else is going on in the system (winds, tides, etc.), but as it requires more significant changes in the code and its benefit for realistic CBD lifetimes is small, perhaps it'd be better to cut it. 

## Additional post-CEE eccentricity functionality
COMPAS currently sets the post-CEE eccentricity to zero, which makes sense since it is a big unknown, which usually doesn't impact the results much. However, since the CBD-driven orbital evolution, and particularly eccentricity pumping, is sensitive to the eccentricity at disk formation, I added a user-supplied option for '--post-common-envelope-eccentricity-prescription'. The options are

1. 'CIRCULAR' (default) – the regular COMPAS behavior.
2. 'PRE_RLOF_CAP' -  the user provides a maximal post-CEE eccentricity through '--post-common-envelope-eccentricity', the pre-RLOF eccentricity is retained if it was below the cap.
3. 'PRE_RLOF_FRACTION'  - the post-CEE eccentricity is set to a user-provided fraction of the pre-RLOF eccentricity through '--post-common-envelope-eccentricity-fraction'.



