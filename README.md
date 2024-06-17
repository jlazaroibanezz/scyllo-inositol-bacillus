# scyllo-inositol-bacillus

This repository contains the code necessary to run the simulations performed in the conference paper "Flexible Nets to Improve GEM Cell Factories by
Combining Kinetic and Proteomics Data" submitted to the 22nd conference on Computational Methods in Systems Biology (CMSB 2024).

The generate_sc_iYO844_def.py file contains the code necessary to implement the reactions necessary for the production of the metabolite 
of interest, the scyllo-inositol. The file needs the previous installation of the COBRA package and needs to be in the same directory as
the constraint-based model that shall be modified, in this case, the iYO844 model from Bacillus subtilis.

## COBRApy installation 

```
pip install cobra
```
## Generation of scyllo-inositol producing model

Run in a terminal the file generate_sc_iYO844_def.py

```
python3 generate_sc_iYO844_def.py
```

The result is a modified constraint-based model named sc_iYO844 with the capacity to synthesize scyllo-inositol in the working directory.


