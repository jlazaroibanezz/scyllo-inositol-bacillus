# scyllo-inositol-bacillus

This repository contains the code necessary to run the simulations performed in the conference paper "Flexible Nets to Improve GEM Cell Factories by
Combining Kinetic and Proteomics Data" submitted to the 22nd conference on Computational Methods in Systems Biology (CMSB 2024).

The generate_sc_iYO844_def.py file contains the code required to implement the reactions necessary for the production of the metabolite 
of interest, the scyllo-inositol. The file needs the previous installation of the COBRA package and needs to be in the same directory as
the constraint-based model that shall be modified, in this case, the [iYO844.xml](http://bigg.ucsd.edu/models/iYO844) model from Bacillus subtilis.

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
CPLEX. In our case, we used [CPLEX](https://www.ibm.com/es-es/products/ilog-cplex-optimization-studio). However, GLPK is the most accessible one.

Once installed the solver and the fnyzer Python package, run in a terminal:

```
$ python3 gecko_smom_FNs_def.py -h 
```
This command shows an overview of the tunable parameters:

```
usage: gecko_smom_FNs_def_jj.py [-h] [-g GROWTH] [-u GLC_UPTAKE_BOUND] [-t CONS_TYPE]

options:
  -h, --help            show this help message and exit
  -g GROWTH, --growth GROWTH
                        growth rate (h-1)
  -u GLC_UPTAKE_BOUND, --glc_uptake_bound GLC_UPTAKE_BOUND
                        upper bound for the glucose uptake rate (mmol gDW-1 h-1)
  -t CONS_TYPE, --cons_type CONS_TYPE
                        Type of constraints added: g+s --> GECKO+sMOMENT constraints; g
                        --> GECKO constraints; n --> No enzymatic constraints; s -->
                        sMOMENT constraints
```

Three parameters can be modified in order to obtain the scyllo-inositol production values generated in the simulations that constitute the graphs in Figures 7 and 8 
in the paper "Flexible Nets to Improve GEM Cell Factories by Combining Kinetic and Proteomics Data":

1. The parameter growth (0.3 h-1 by default): Sets the fixed growth rate to perform the simulations while the scyllo-inositol production is maximized.
2. The parameter glc_uptake (100 mmol gDW-1 h-1 by default): The glucose uptake rate is unconstraint by default, but it can be adjusted to the experimental value 7.71  mmol gDW-1 h-1.
3. The parameter cons_type: Sets the type of constraints that are integrated in the sc_iYO844 model.

For instance, consider that we want to compute the scyllo-inositol production and the glucose uptake rate when the growth rate is 0.032 and applying the only the GECKO constraints: 

```
$ python3 gecko_smom_FNs_def.py -g 0.032 -t g
```
The expected output would be:

```
GECKO constraints
Growth rate 0.032 h-1
Glucose uptake rate 24.79710095198935 mmol gDW-1 h-1
Scyllo-inositol production:  19.74652927198935 mmol gDW-1 h-1
```


