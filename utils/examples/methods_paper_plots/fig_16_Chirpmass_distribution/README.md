# How to recreate/edit the BBH chirp mass distribution from the COMPAS methods paper

This folder contains everything you need to reproduce the BBH chirpmass distribution from the COMPAS methods paper (_arXiv_: https://arxiv.org/abs/2109.10352;  Team COMPAS et al. 2021). 


The data used for this figure is publicly available at: https://zenodo.org/record/5655483

This data set contains the output of 10,000,000 binaries evolved using COMPAS 02.21.00, using adaptive importance sampling (STROOPWAFEL, Broekgaarden et al. 2019), sampling from a metallicity uniform in $\log(Z) \in [10^{-4},0.03]$. More details can be found in `Run_Details.txt`.

### Data reporduction
The data can be reproduced by running version `02.21.00` of COMPAS, 

1. Run `stroopwafel_interface.py`, that reads in the `Fig16_pythonSubmit.py.py` (both contained in this folder).
2. Calculate the rates by running ```FastCosmicIntegration.py```  from COMPAS's post-processing tools, with the following flags altered from their default values:


```:::bash
    python FastCosmicIntegration.py  --mu0 0.035 --muz -0.23 --sigma0 0.39 --sigmaz 0.0 --alpha 0.0 --weight mixture_weight --zstep 0.01 --sens O3 --m1min 10. --aSF 0.01 --bSF 2.77 --cSF 2.9 --dSF 4.7 
```

