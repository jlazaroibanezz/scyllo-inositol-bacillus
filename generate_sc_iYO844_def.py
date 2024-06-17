### This script generates the contraint based model sc_iYO844.xml
### which is the result of adding the reactions described in
### "Flexible Nets to Improve GEM Cell Factories by Combining Kinetic
###   and Proteomics Data" J. Lazaro et al. for the production
### of scyllo-inositol to the constraint based iYO844.xml

import cobra
from cobra import Model, Reaction, Metabolite

model = cobra.io.read_sbml_model("iYO844.xml")  # Load model

# Transform glucose-6-phosphate into myo-inositol-1-phosphate
reaction = Reaction('MI1P_synthesis')
reaction.name = 'Synthesis of myo-inositol-1-phosphate'
reaction.lower_bound = 0  
reaction.upper_bound = 1000

g6p_c = model.metabolites.get_by_id("g6p_c") 
mi1p__D_c = model.metabolites.get_by_id("mi1p__D_c")

reaction.add_metabolites({
    g6p_c: -1.0,
    mi1p__D_c: 1.0
})

# Add scyllo-inositol and synthesis of it by using 2-inosose
reaction2 = Reaction('scino_synthesis')
reaction2.name = 'Synthesis of scyllo-inositol'
reaction2.lower_bound = 0  
reaction2.upper_bound = 1000  

ins_c = model.metabolites.get_by_id("2ins_c") 
nadph_c = model.metabolites.get_by_id("nadph_c") 
nadp_c = model.metabolites.get_by_id("nadp_c") 
h_c = model.metabolites.get_by_id("h_c") 
scino_c = Metabolite(
    'scino_c',
    name='scyllo-inositol',
    compartment='c')

reaction2.add_metabolites({
    ins_c: -1.0,
    scino_c: 1.0,
    nadph_c: -1.0,
    nadp_c : 1.0,
    h_c : -1.0
})

# Add exchange reaction for the scyllo-inositol
reaccEX = Reaction('EX_scino')
reaccEX.name = 'Exchange reaction to allow scyllo-inositol to leave the system'
reaccEX.lower_bound = 0.0
reaccEX.upper_bound = 1000.0
reaccEX.add_metabolites({scino_c: -1.0})
model.add_reactions([reaction, reaction2, reaccEX])

# Knock out the reaction that is knocked out in the original paper
model.reactions.get_by_id('INSCR').lower_bound = 0
model.reactions.get_by_id('INSCR').upper_bound = 0
model.reactions.get_by_id('INSCR')
# Set the glucose maximum uptake rate to the one in the enzimatically constraint model to compare results
model.reactions.get_by_id('EX_glc__D_e').lower_bound = -100

# Set the scyllo-inositol synthesis as the objective and solve
model.objective = "scino_synthesis"
solution = model.optimize()
print('The scyllo-inositol production is: ', solution.objective_value, 'mmol gDW-1 h-1')

cobra.io.write_sbml_model(model, 'sc_iYO844.xml')

