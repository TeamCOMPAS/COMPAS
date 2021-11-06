---
title: 'COMPAS: A rapid binary population synthesis suite'
tags:
  - Python
  - C++
  - astronomy
  - gravitational waves
  - binary evolution
authors:
  - name: Team COMPAS
    affiliation: "1" # (Multiple affiliations must be quoted)
  - name: Jeff Riley 
    affiliation: "2,3"
  - name: Poojan Agrawal
    affiliation: "4,3"
  - name: Jim W. Barrett
    affiliation: 5
  - name: Kristian N. K. Boyett
    affiliation: 6
  - name: Floor S. Broekgaarden
    affiliation: 7
  - name: Debatri Chattopadhyay
    affiliation: "8,4,3"
  - name: Sebastian M. Gaebel
    affiliation: 9
  - name: Fabian Gittins
    affiliation: 10
  - name: Ryosuke Hirai
    affiliation: "2,3"
  - name: George Howitt
    affiliation: 11
  - name: Stephen Justham
    affiliation: "12,13,14"
  - name: Lokesh Khandelwal
    affiliation: 12
  - name: Floris Kummer
    affiliation: 12
  - name: Mike Y. M. Lau
    affiliation: "2,3"
  - name: Ilya Mandel
    affiliation: "2,3,5"
  - name: Selma E. de Mink
    affiliation: "14,12,7"
  - name: Coenraad Neijssel
    affiliation: "5,3"
  - name: Tim Riley
    affiliation: "2,3"
  - name: Lieke van Son
    affiliation: "7,12,14"
  - name: Simon Stevenson
    affiliation: "4,3"
  - name: Alejandro Vigna-Gómez
    affiliation: "15,16"
  - name: Serena Vinciguerra
    affiliation: 12
  - name: Tom Wagg
    affiliation: "7,14"
  - name: Reinhold Willcox
    affiliation: "2,3"
affiliations:
 - name: The public COMPAS code is a product of work by the entire COMPAS collaboration over many years; we therefore kindly request that, in recognition of this team effort, the paper is cited as Team COMPAS - J. Riley et al.
   index: 1
 - name: School of Physics and Astronomy, Monash University, Clayton, Victoria 3800, Australia
   index: 2
 - name: OzGrav, Australian Research Council Centre of Excellence for Gravitational Wave Discovery, Australia
   index: 3
 - name: Centre for Astrophysics and Supercomputing, Swinburne University of Technology, Hawthorn, VIC 3122, Australia
   index: 4
 - name: Institute of Gravitational Wave Astronomy and School of Physics and Astronomy, University of Birmingham, Birmingham, B15 2TT
   index: 5
 - name: Department of Physics, University of Oxford, Denys Wilkinson Building, Keble Road, Oxford OX1 3RH, UK
   index: 6
 - name: Center for Astrophysics |Harvard & Smithsonian, 60 Garden St., Cambridge, MA 02138, USA
   index: 7
 - name: School of Physics and Astronomy, Cardiff University, Cardiff, CF24 3AA, United Kingdom
   index: 8
 - name: Max Planck Institute for Gravitational Physics (Albert Einstein Institute), Callinstrasse 38, D-30167 Hannover, Germany
   index: 9
 - name: Mathematical Sciences and STAG Research Centre, University of Southampton, Southampton SO17 1BJ, UK
   index: 10
 - name: School of Physics, University of Melbourne, Parkville, Victoria, 3010, Australia
   index: 11
 - name: Anton Pannekoek Institute of Astronomy and GRAPPA, Science Park 904, University of Amsterdam, 1098XH Amsterdam, The Netherlands
   index: 12
 - name: School of Astronomy & Space Science, University of the Chinese Academy of Sciences, Beijing 100012, China
   index: 13
 - name: Max Planck Institute for Astrophysics, Karl-Schwarzschild-Str. 1, 85748 Garching, Germany
   index: 14
 - name: DARK, Niels Bohr Institute, University of Copenhagen, Jagtvej 128, 2200, Copenhagen, Denmark
   index: 15
 - name: Niels Bohr International Academy, The Niels Bohr Institute, Blegdamsvej 17, 2100 Copenhagen, Denmark
   index: 16
date: xx Month 2021
bibliography: paper.bib

# Optional fields if submitting to a AAS journal too, see this blog post:
# https://blog.joss.theoj.org/2018/12/a-new-collaboration-with-aas-publishing
aas-doi: 10.3847/xxxxx <- update this with the DOI from AAS once you know it.
aas-journal: Astrophysical Journal Supplements <- The name of the AAS journal.
---

# Summary

Most massive stars---those with initial masses greater than 8 $M_\odot$---are born with another massive star as a companion [@Sana:2012Sci;@Moe:2017ApJS]. Massive binary stars are responsible for producing many exotic astrophysical phenomena, such as the observed diversity of supernovae, binary pulsars, X-ray binaries and merging compact objects. The latter are now regularly observed by the ground-based gravitational wave observatories Advanced LIGO and Virgo [@abbott2016observation;@GWTC2].  Population models of massive binary evolution make it possible to interpret existing observations and to make predictions for future observing campaigns.  

# Statement of need

Binary population synthesis generates population models of isolated stellar binaries under a set of parametrized assumptions.  These models permit comparisons against observational data sets, such as X-ray binaries of gravitational-wave mergers.   

In particular, rapid binary population synthesis is needed in order to efficiently explore a broad parameter space of uncertain assumptions about the physics of stellar and binary evolution, including supernova remnant masses and natal kicks, mass transfer efficiency and stability, and the outcome of common-envelope events.  

A range of binary population synthesis codes have been developed over the last three decades.  These include the Scenario Machine [@Scenario], IBiS [@IBiS], SeBa [@SeBa],  BSE [@Hurley:2002rf], StarTrack [@Belczynski:2008],  binary$\_$c [@BinaryC],  MOBSE [@2018MNRAS.474.2959G]  and COSMIC [@2019arXiv191100903B].  These codes range from private to semi-public to fully public, and differ in the range of available tools, computational complexity, and speed of execution.

[COMPAS](https://compas.science) is a rapid binary population synthesis suite. It parametrizes complex astrophysical processes with prescriptions calibrated to detailed models.  COMPAS is designed to allow for flexible modifications as evolutionary models improve.  All code is fully public and, including pre-processing and post-processing tools.  COMPAS is computationally efficient, with a focus on the statistical analysis of large populations, particularly but not exclusively in the context of gravitational-wave astronomy.  


# Details

The core engine of COMPAS---responsible for calculating the evolution of single [@Hurley:2000pk] and binary [@Hurley:2002rf] stars---is written in object oriented C++ for speed and flexibility. COMPAS is able to simulate the evolution of a typical binary over 10 Gyr in approximately 10 milliseconds.

A detailed description of the implementation of the COMPAS suite can be found in @COMPAS:2021methodsPaper.

In addition to the core stellar and binary evolution engine, we provide Python scripts for both pre- and post-processing COMPAS outputs. Post-processing can account for integrating populations formed throughout cosmic history [@2019MNRAS.490.3740N] and methods to account for gravitational-wave selection effects [@Barrett:2017fcw]. A set of examples is also provided.

COMPAS is *embarrassingly* parallel and can be trivially run on high performance computers and distributed on cloud computing.

COMPAS was initially designed to focus on studies of merging binaries containing neutron stars and black holes that are being observed through gravitational waves [@Stevenson2017FormationEvolution;@2018MNRAS.481.4009V]. 
In recent years, the scope of systems investigated with COMPAS has expanded to incorporate, e.g., Be X-ray binaries [@Vinciguerra:2020] and luminous red novae [@Howitt:2020] (see @COMPAS:2021methodsPaper or [the COMPAS collaboration website](https://compas.science) for a summary of COMPAS publications to date.)

COMPAS development happens on [Github](https://github.com/TeamCOMPAS/COMPAS). We maintain a [Zenodo community](https://zenodo.org/communities/compas/) where data from many publications using COMPAS is publicly available. 


# Acknowledgements

Multiple authors are supported by the Australian Research Council Centre of Excellence for Gravitational Wave Discovery (OzGrav), through project number CE170100004. Multiple authors were funded in part by the National Science Foundation under Grant No. (NSF grant number 2009131), the Netherlands Organization for Scientific Research (NWO) as part of the Vidi research program BinWaves with project number 639.042.728 and by the European Union’s Horizon 2020 research and innovation program from the European Research Council (ERC, Grant agreement No. 715063).  FSB is supported in part by the Prins Bernard Cultuurfonds studiebeurs. IM is a recipient of an Australian Research Council Future Fellowship (FT190100574).  AVG acknowledges funding support by the Danish National Research Foundation (DNRF132)


# References
