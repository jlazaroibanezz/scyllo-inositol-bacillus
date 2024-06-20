### This script transforms the constraint based model sc_iYO844.xml to a
### Flexible Net, adds enzymatic constraints as described in the paper
### "Flexible Nets to Improve GEM Cell Factories by Combining Kinetic
###   and Proteomics Data" J. Lazaro et al., and optimizes the production
### of scyllo-inositol

import argparse
from fnyzer import FNFactory
from fnyzer import cobra2fn
import cobra

# Constants
def_growth_rt = 0.3    # Default growth rate (h-1)
def_glc_up_bound = 100 # Default upper bound for the glucose uptake rate (mmol gDW-1 h-1)
def_cons = 'g+s'       # Default model constraints: GECKO+sMOMENT constraints
P = 0.01067            # Total enzyme abundance (g gDW-1)

# Parse parameters
parser = argparse.ArgumentParser()
parser.add_argument("-g", "--growth", default=def_growth_rt, help="growth rate (h-1)")
parser.add_argument("-u", "--glc_uptake_bound", default=def_glc_up_bound, help="upper bound for the glucose uptake rate (mmol gDW-1 h-1)")
parser.add_argument("-t", "--cons_type", default=def_cons, help="Type of constraints added:\n g+s --> GECKO+sMOMENT constraints; g --> GECKO constraints; n --> No enzymatic constraints; s --> sMOMENT constraints")
args = parser.parse_args()
growth_rt = args.growth
glc_uptake_bound = args.glc_uptake_bound
cons_type = args.cons_type

# Read model and transform to Flexible Net
model = cobra.io.read_sbml_model('sc_iYO844.xml')
fnet = cobra2fn(model) 
fnet['solver'] = 'cplex_direct' # 'glpk' # cplex_direct, gurobi

# Reactions with enzymatic constraints
#    typ: 'fonly':  Irreversible and forward only reaction
#         'fcat':   Reversible and forward-catalyzed reaction 
#         'bcat':   Reversible and backward-catalyzed reaction 
#         'scillo': Scyllo-inositol synthesis reactions
#   kcat: kcat (h-1)
#     ea: Enzyme abundance (mmol gDW-1)
#     mw: Molecular weight (kDa)
rxns = {'G6PDH2r': {'typ': 'fonly',  'kcat': 626400  , 'ea': 8.05e-06, 'mw': 58   },
        'CS':      {'typ': 'fonly',  'kcat': 176400  , 'ea': 2.51e-05, 'mw': 42   },
        'GAPD':    {'typ': 'fonly',  'kcat': 252000  , 'ea': 5.77e-05, 'mw': 35   },
        'OXADC':   {'typ': 'fonly',  'kcat': 212400  , 'ea': 6.21e-07, 'mw': 43   },
        'OXGDC':   {'typ': 'fonly',  'kcat': 720     , 'ea': 6.80e-08, 'mw': 63.92},
        'PGI':     {'typ': 'fcat',   'kcat': 453600  , 'ea': 1.55e-05, 'mw': 61.53},
        'TPI':     {'typ': 'fcat',   'kcat': 540000  , 'ea': 1.28e-05, 'mw': 26.97},
        'PGK':     {'typ': 'bcat',   'kcat': 1184400 , 'ea': 3.60e-05, 'mw': 41.12},
        'PGM':     {'typ': 'bcat',   'kcat': 2757240 , 'ea': 8.85e-06, 'mw': 74   },
        'ENO':     {'typ': 'fcat',   'kcat': 4.6943e5, 'ea': 3.17e-05, 'mw': 46.58},
        'ICDHyr':  {'typ': 'fcat',   'kcat': 295200  , 'ea': 1.10e-04, 'mw': 46   },
        'FUM':     {'typ': 'fcat',   'kcat': 614988  , 'ea': 7.29e-06, 'mw': 50   },
        'LDH_L':   {'typ': 'fcat',   'kcat': 23099976, 'ea': 3.60e-06, 'mw': 35   },
        'PGCD':    {'typ': 'fcat',   'kcat': 52416   , 'ea': 1.90e-05, 'mw': 56   },
        'MCITL2':  {'typ': 'fcat',   'kcat': 68400   , 'ea': 6.80e-08, 'mw': 32   },
        'PTAr':    {'typ': 'fcat',   'kcat': 2345760 , 'ea': 8.49e-06, 'mw': 34   },
        'MDH':     {'typ': 'fcat',   'kcat': 637560  , 'ea': 1.06e-04, 'mw': 33   },
        'MI1PP':   {'typ': 'scillo', 'kcat': 23400   , 'ea': 1.06e-06, 'mw': 29.61},
        'INS2D':   {'typ': 'scillo', 'kcat': 78624   , 'ea': 3.18e-06, 'mw': 38.2 },
 'scino_synthesis':{'typ': 'scillo', 'kcat': 304668  , 'ea': 2.08e-06, 'mw': 39.95}
}

### Add enzymatic constraints to each reaction
fnet['shandlers'] = {}
fnet['extracons'] = []

fnet['places']['E_total'] = P  # Total enzyme abundance

if cons_type == 'g+s': # GECKO+sMOMENT constraints
    for r in rxns:
        typ, kcat, ea, mw = rxns[r]['typ'], rxns[r]['kcat'], rxns[r]['ea'], rxns[r]['mw']
        if typ in ['fonly', 'fcat']:
            fnet['places']['E_'+r+'_f'] = ea
            fnet['trans']['t_'+r+'_f']['l0'] = 0
            fnet['shandlers']['s_'+r+'_f'] = [{'e':('E_'+r+'_f','s_'+r+'_f'),
                                               'r':('s_'+r+'_f','t_'+r+'_f'),
                                               'p':('E_total','s_'+r+'_f')},
                                              'r<='+str(kcat)+'*'+str(ea),
                                              'r<='+str(kcat)+'*p'+'/'+str(mw)]
            if typ == 'fcat':
                fnet['extracons'].append("l0['t_"+str(r)+"_b'] == 0")
			
        elif typ == 'bcat':
            fnet['places']['E_'+r+'_b'] = ea
            fnet['trans']['t_'+r+'_b']['l0'] = 0
            fnet['extracons'].append("l0['t_"+r+"_f'] == 0")
            fnet['shandlers']['s_'+r+'_b'] = [{'e':('E_'+r+'_b','s_'+r+'_b'),
                                               'r':('s_'+r+'_b','t_'+r+'_b'),
                                               'p':('E_total','s_'+r+'_b')},
                                              'r<='+str(kcat)+'*'+str(ea),
                                              'r<='+str(kcat)+'*p'+'/' + str(mw)]
	
        elif typ == 'scillo':
            fnet['places']['E_'+r+'_f'] = ea
            fnet['trans']['t_'+r+'_f']['l0'] = 0
            fnet['shandlers']['s_'+r+'_f'] = [{'e':('E_'+r+'_f','s_'+r+'_f'),
                                               'r':('s_'+r+'_f','t_'+r+'_f'),
                                               'p':('E_total','s_'+r+'_f')},
                                              'r>='+str(kcat)+'*'+str(ea),
                                              'r<='+str(kcat)+'*p'+'/'+str(mw)]
                           
    if glc_uptake_bound != def_glc_up_bound:
        print("GECKO+sMOMENT+glucose uptake rate constraints")
    else:
        print("GECKO+sMOMENT constraints")
                           
elif cons_type == 'g': # GECKO constraints
    for r in rxns:
        typ, kcat, ea = rxns[r]['typ'], rxns[r]['kcat'], rxns[r]['ea']
        if typ in ['fonly', 'fcat']:
            fnet['places']['E_'+r+'_f'] = ea
            fnet['trans']['t_'+r+'_f']['l0'] = 0
            fnet['shandlers']['s_'+r+'_f'] = [{'e':('E_'+r+'_f','s_'+r+'_f'),
                                               'r':('s_'+r+'_f','t_'+r+'_f')},
                                              'r<='+str(kcat)+'*'+str(ea)]
            if typ == 'fcat':
                fnet['extracons'].append("l0['t_"+str(r)+"_b'] == 0")
			
        elif typ == 'bcat':
            fnet['places']['E_'+r+'_b'] = ea
            fnet['trans']['t_'+r+'_b']['l0'] = 0
            fnet['extracons'].append("l0['t_"+r+"_f'] == 0")
            fnet['shandlers']['s_'+r+'_b'] = [{'e':('E_'+r+'_b','s_'+r+'_b'),
                                               'r':('s_'+r+'_b','t_'+r+'_b')},
                                              'r<='+str(kcat)+'*'+str(ea)]
	
        elif typ == 'scillo':
            fnet['places']['E_'+r+'_f'] = ea
            fnet['trans']['t_'+r+'_f']['l0'] = 0
            fnet['shandlers']['s_'+r+'_f'] = [{'e':('E_'+r+'_f','s_'+r+'_f'),
                                               'r':('s_'+r+'_f','t_'+r+'_f')},
                                              'r>='+str(kcat)+'*'+str(ea)]
                           
    if glc_uptake_bound != def_glc_up_bound:
        print("GECKO+glucose uptake rate constraints")
    else:
        print("GECKO constraints")

elif cons_type == 's': # sMOMENT constraints
    for r in rxns:
        typ, kcat, mw = rxns[r]['typ'], rxns[r]['kcat'], rxns[r]['mw']
        if typ in ['fonly', 'fcat']:
            fnet['trans']['t_'+r+'_f']['l0'] = 0
            fnet['shandlers']['s_'+r+'_f'] = [{'r':('s_'+r+'_f','t_'+r+'_f'),
                                               'p':('E_total','s_'+r+'_f')},
                                              'r<='+str(kcat)+'*p/'+str(mw)]
            if typ == 'fcat':
                fnet['extracons'].append("l0['t_"+str(r)+"_b'] == 0")
			
        elif typ == 'bcat':
            fnet['trans']['t_'+r+'_b']['l0'] = 0
            fnet['extracons'].append("l0['t_"+r+"_f'] == 0")
            fnet['shandlers']['s_'+r+'_b'] = [{'r':('s_'+r+'_b','t_'+r+'_b'),
                                               'p':('E_total','s_'+r+'_b')},
                                              'r<='+str(kcat)+'*p/' + str(mw)]
	
        elif typ == 'scillo':
            fnet['trans']['t_'+r+'_f']['l0'] = 0
            fnet['shandlers']['s_'+r+'_f'] = [{'r':('s_'+r+'_f','t_'+r+'_f'),
                                               'p':('E_total','s_'+r+'_f')},
                                              'r<='+str(kcat)+'*p/'+str(mw)]
                           
    if glc_uptake_bound != def_glc_up_bound:
        print("sMOMENT+glucose uptake rate constraints")
    else:
        print("sMOMENT constraints")
        
elif cons_type == 'n': # No constraints
    print("No enzymatic constraints")

else:
    print("Incorrect constraints type. Use -h for help.")
    exit()
	
### Set parameters and optimize
fnet['obj'] = {'f': "l0['t_EX_scino_f']", 'sense': 'max'}  # Define objective function
fnet['extracons'].append("l0['t_BIOMASS_BS_10_f']== "+str(growth_rt)) # Fix growth rate
fnet['extracons'].append("l0['t_EX_glc__D_e_b']<= "+str(glc_uptake_bound))  # Upper bound glucose uptake rate 
fnet['options'] = {'antype': 'cst'}  # Constant steady-state
netobj = FNFactory(fnet) 
netobj.optimize()
solution = netobj.objval 

print("Growth rate:", growth_rt, "h-1")
print("Glucose uptake rate:", netobj.trans['t_EX_glc__D_e_b'].avl, 'mmol gDW-1 h-1')
print('Scyllo-inositol production:', solution, 'mmol gDW-1 h-1')
