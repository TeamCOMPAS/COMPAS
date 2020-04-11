[//]: ## (grip -b README.md)

![COMPASlogo](docs/COMPASlogo.png)

# Compact Object Mergers: Population Astrophysics & Statistics

[//]: ## (Outline features)
COMPAS is a publicly available rapid binary population synthesis code (https://compas.science/) that is developed for easily adjustable evolution prescriptions and model parameters. COMPAS draws binary properties from a set of initial distributions, and evolves a binary star from zero-age main sequence to the end of its life as two compact remnants. It has been used for inference from observations of gravitational wave mergers, Galactic double neutron stars, and radio pulsars.

## Documentation
[Getting started](./docs/getting_started.md)  
[Specifications](./docs/COMPAS_Doc.pdf)  
[Post-processing](./postProcessing/README.txt)  
[git workflow for developers](./docs/dev-git-workflow.md)  

## Contact
Please email your queries to compas-user@googlegroups.com. You are also welcome to join the [COMPAS User Google Group](https://groups.google.com/forum/#!members/compas-user) to engage in discussions with COMPAS users and developers.

## Acknowledgements
If you use this code or parts of this code for results presented in a scientific publication, we would greatly appreciate if you send us your paper reference and make your input settings and output data publicly available by uploading it to the [COMPAS Zenodo community](https://zenodo.org/communities/compas/). Please also kindly include citations to the two following papers:

1. Stevenson S., Vigna-Gómez A., Mandel I., Barrett J. W., Neijssel C. J., Perkins D., de Mink S. E., 2017, [Nature Communications, 8, 14906](https://ui.adsabs.harvard.edu/abs/2017NatCo...814906S/abstract)
2. Vigna-Gómez A., Neijssel C. J., Stevenson S., Barrett J. W., Belczynski K., Justham S., de Mink S., M&uuml;ller B., Podsiadlowski Ph., Renzo M., Szécsi D., Mandel I., 2018, [MNRAS, 481, 4009](https://ui.adsabs.harvard.edu/abs/2018MNRAS.481.4009V/abstract)

We anticipate releasing a more detailed and comprehensive methods paper in the future. We also greatly appreciate an acknowledgement of the form: 

>_Simulations in this paper made use of the COMPAS rapid binary population synthesis code which is freely available at http://github.com/TeamCOMPAS/COMPAS_.

Furthermore,
  * If you use COMPAS's importance sampling algorithm STROOPWAFEL, please cite 

     Broekgaarden F. S., Justham S., de Mink S. E., Gair J., Mandel I., Stevenson S., Barrett J. W., Vigna-Gómez A., Neijssel C. J., 2019, [MNRAS, 490, 5228](https://ui.adsabs.harvard.edu/abs/2019MNRAS.490.5228B/abstract)

  * If using the COMPAS model of gravitational wave selection effects, please cite

     Barrett J. W., Gaebel S. M., Neijssel C. J., Vigna-Gómez A., Stevenson S., Berry C. P. L., Farr W. M., Mandel I., 2018, [MNRAS, 477, 4685](https://ui.adsabs.harvard.edu/abs/2018MNRAS.477.4685B/abstract)

  * If using COMPAS's integration over cosmic star formation history, please cite 

     Neijssel C. J., Vigna-Gómez A., Stevenson S., Barrett J. W., Gaebel S. M., Broekgaarden F. S., de Mink S. E., Szécsi D., Vinciguerra S., Mandel I., 2019, [MNRAS, 490, 3740](https://ui.adsabs.harvard.edu/abs/2019MNRAS.490.3740N/abstract)

  * If using the COMPAS model of (pulsational) pair instability supernova, please cite 

     Stevenson S., Sampson M., Powell J., Vigna-Gómez A., Neijssel C. J., Szécsi D., Mandel I., 2019, [ApJ, 882, 121](https://ui.adsabs.harvard.edu/abs/2019ApJ...882..121S/abstract)
     
  * If evolving pulsar spins and magnetic fields with COMPAS, please cite
  
     Chattopadhyay D., Stevenson S., Hurley J. R., Rossi L. J., Flynn C., 2020,  [MNRAS](https://ui.adsabs.harvard.edu/abs/2020MNRAS.tmp..697C/abstract)

## License
[MIT](https://choosealicense.com/licenses/mit/)
