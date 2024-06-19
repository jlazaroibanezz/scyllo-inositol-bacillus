# scyllo-inositol-bacillus

This repository contains the code necessary to run the simulations performed in the conference paper "Flexible Nets to Improve GEM Cell Factories by
Combining Kinetic and Proteomics Data" submitted to the 22nd conference on Computational Methods in Systems Biology (CMSB 2024).

The generate_sc_iYO844_def.py file contains the code required to implement the reactions necessary for the production of the metabolite 
of interest, the scyllo-inositol. The file needs the previous installation of the COBRA package and needs to be in the same directory as
the constraint-based model that shall be modified, in this case, the iYO844 model from Bacillus subtilis.

## COBRApy installation 

Open a terminal and run the following command:

```
$ pip install cobra
```
## Generation of scyllo-inositol producing model

Run in a terminal the file generate_sc_iYO844_def.py using Python:

```
$ python3 generate_sc_iYO844_def.py
```

The result is a modified constraint-based model named sc_iYO844 with the capacity to synthesize scyllo-inositol in the working directory.

## Simulations using Flexible Nets

Install the fnyzer package: 

```
$ pip install fnyzer
```
To run the simulations with fnyzer it is necessary to download a solver that will solve the linear programming problem. We suggest GLPK, Gurobi or
CPLEX. In our case, we used [CPLEX](https://www.ibm.com/es-es/products/ilog-cplex-optimization-studio).

Once installed the solver and the fnyzer Python package, run in a terminal:

```
$ python3 gecko_smom_FNs_def.py
```
The parameter growth (0.3 h-1 by default) can be adjusted to obtain the different scyllo-inositol production values generated in the simulations 
that constitute the graphs in Figures 7 and 8 in the paper "Flexible Nets to Improve GEM Cell Factories by Combining Kinetic and Proteomics Data".

For instance:

```
$ python3 gecko_smom_FNs_def.py -g 0
```
The expected output would be:

```
GECKO+sMOMENT constraints
Growth rate 0 h-1
Scyllo-inositol production:  5.303801878612387 mmol gDW-1 h-1
```



