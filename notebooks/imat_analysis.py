#!/usr/bin/env python
# coding: utf-8

# # Kinase IMAT Modeling

# This notebook explores the impact that Serine-Threonine protein kinases have on *Mycobacterium tuberculosis* metabolism by integration of gene expression data with Genome Scale Metabolic Models (GSMM). This is done through use of the IMAT algorithm, described by [Shlomi et al., 2008](https://www.nature.com/articles/nbt.1487). This algorithm attempts to find an optimal trade off between having high flux through high expression reactions, and low flux through low expression reactions in order to determine a flux state which matches provided gene expression data. This flux state can then be used to generate a condition specific model used for futher analysis and modeling. 

# ## IMAT Algorithm

# The IMAT algorithm:
# 1. Categorizes gene expression into low expression, high expression, and inbetween (normally expressed as -1, 1, and 0 respectively)
# 2. Adds binary indicator variables ($y_i^+$ and $y_i^-$) to the metabolic model to represent whether a reaction with a high expression value, has a flux above a cutoff ($\epsilon$) in the forward or reverse direction (+ for forward, - for reverse)
# 3. Adds binary indicator variables ($y_i^+$) to the metabolic model to represent whether a reaction with a low expression value, has a flux below a cutoff ($\epsilon$)
# 4. Maximizes the sum of high expression reactions active in either direction, plus the sum of low expression reactions that are inactive
# 
# This can be expressed as:  
# 
# $max_{v,y^+,y^-}(\sum_{i \in R_H} (y_i^+ + y_i^-) + \sum_{i \in R_L}y_i^+) $  
#   
# Such That  
# 
#   
# $S \cdot v =0 $    
# <br>
# $v_{min}\leq v \leq v_{max}$    
# <br>
# $v_i + y_i^+(v_{min,i}-\epsilon) \geq v_{min,i}; i\in R_H$  
# <br>
# $v_i + y_i^-(v_{max,i}+\epsilon) \leq v_{max,i}; i \in R_H$  
# <br>
# $v_{min,i}(1-y_i^+) \leq v_i \leq v_{max,i}(1-y_i^+);i \in R_L$  
# <br>
# $v \in R^m$    
# <br>
# $y_i^+, y_i^- \in [0,1]$

# ## Imports

# In[1]:


# Standard Library imports
import pathlib

# External Imports
import cobra
import metworkpy
import numpy as np
import pandas as pd

# Local Imports


# ## Setup
# Set the default solver to Gurobi to significantly speed up the analysis

# In[2]:


cobra.Configuration.solver = "gurobi"


# ## Data Preprocessing

# The first step is to read in the expression data, and convert it into trinarized reaction information.

# ### Read in Expression Data

# In[3]:


gene_expression = pd.read_csv(pathlib.Path('..') / 'data' / 'gene-expression' / 'kinase_RNA_seq.csv', index_col=0).drop( "GrowthPhase", axis=1)
# Convert the RPKM measurements into TPM
gene_expression = metworkpy.utils.rpkm_to_tpm(gene_expression)


# In[4]:


# Separate the expression data into biological replicates
condition_dict = {
    "WT":None,
    "PknB_KD":None,
    "PknB_IND":None,
    "PknD_KO":None,
    "PknD_IND":None,
    "PknE_KO":None,
    "PknE_IND":None,
    "PknF_KO":None,
    "PknF_IND":None,
    "PknG_KO":None,
    "PknG_IND":None,
    "PknH_KO":None,
    "PknH_IND":None,
    "PknJ_KO":None,
    "PknJ_IND":None,
    "PknK_KO":None,
    "PknK_IND":None,
    "PknL_KO":None,
    "PknL_IND":None,
}
for condition in condition_dict:
    condition_dict[condition] = gene_expression.loc[gene_expression.index.str.startswith(condition)]


# ### Convert Gene Expression into Trinarized Reactions
# Use the Gene-Protein-Reaction rules from the genome scale metabolic model to convert the gene expression data into trinarized reaction data.

# In[5]:


# Convert gene expression into qualitative weights
for condition in condition_dict:
    condition_dict[condition] = metworkpy.utils.expr_to_weights(condition_dict[condition], quantile = 0.15, 
                                                                aggregator=np.median, sample_axis=0)


# In[6]:


# Convert qualitative weights into trinarized reactions

# Load the iEK1011 model
iek1011 = metworkpy.read_model(pathlib.Path('..') / 'data' / 'Models' / 'iEK1011_m7H10_media.json')

# Convert the gene expression data into trinarized reaction data
# WARNING: This step may take a while, it is not implemented very efficiently currently
for condition in condition_dict:
    p = pathlib.Path('..') / 'cache' / f'{condition}_rxn_weights.json'
    if p.exists():
        condition_dict[condition] = pd.read_json(p, typ = "series")
    else:
        condition_dict[condition] = metworkpy.parse.gene_to_rxn_weights(iek1011, condition_dict[condition])
        condition_dict[condition].to_json(p)


# ## IMAT
# Now that the data has been processed, IMAT can be used to generate condition specific models in several ways. 

# First, we can create a dataframe to hold the biomass growth for all methods across all conditions

# In[7]:


biomass_growth_df = pd.DataFrame(0., index = condition_dict.keys(), columns=["simple_bounds", "subset_model", "fva_model", "milp_model"])


# **WARNING**: This section takes a *very* long time to run! 

# ### Simple Bounds

# In[8]:


METHOD = "simple_bounds"
model_dict = dict()

print(METHOD)
for condition in condition_dict:
    print(condition)
    p = pathlib.Path('..') / 'results' /  METHOD / f'iek1011_{condition}_model.json'
    if p.exists():
        model_dict[condition] = metworkpy.read_model(str(p))
    else:
        model_dict[condition] = metworkpy.imat.generate_model(model=iek1011.copy(), rxn_weights=condition_dict[condition], method=METHOD)
        metworkpy.write_model(model_dict[condition], model_path=p)   
    biomass_growth_df.loc[condition, METHOD] = model_dict[condition].slim_optimize()


# ### Subset Model

# In[9]:


METHOD = "subset_model"
model_dict = dict()
print(METHOD)
for condition in condition_dict:
    print(condition)
    p = pathlib.Path('..') / 'results' /  METHOD / f'iek1011_{condition}_model.json'
    if p.exists():
        model_dict[condition] = metworkpy.read_model(str(p))
    else:
        model_dict[condition] = metworkpy.imat.generate_model(model=iek1011.copy(), rxn_weights=condition_dict[condition], method=METHOD)
        metworkpy.write_model(model_dict[condition], model_path=p)
    biomass_growth_df.loc[condition, METHOD] = model_dict[condition].slim_optimize()


# ### FVA Model

# In[10]:


METHOD = "fva_model"
model_dict = dict()
print(METHOD)
for condition in condition_dict:
    print(condition)
    p = pathlib.Path('..') / 'results' /  METHOD / f'iek1011_{condition}_model.json'
    if p.exists():
        model_dict[condition] = metworkpy.read_model(str(p))
    else:
        model_dict[condition] = metworkpy.imat.generate_model(model=iek1011.copy(), rxn_weights=condition_dict[condition], method=METHOD, loopless=False)
        metworkpy.write_model(model_dict[condition], model_path=p)
    biomass_growth_df.loc[condition, METHOD] = model_dict[condition].slim_optimize()


# ### MILP Model

# In[ ]:


METHOD = "milp_model"
model_dict = dict()
print(METHOD)
for condition in condition_dict:
    print(condition)
    p = pathlib.Path('..') / 'results' /  METHOD / f'iek1011_{condition}_model.json'
    if p.exists():
        model_dict[condition] = metworkpy.read_model(str(p))
    else:
        model_dict[condition] = metworkpy.imat.generate_model(model=iek1011.copy(), rxn_weights=condition_dict[condition], method=METHOD)
        metworkpy.write_model(model_dict[condition], model_path=p)
    biomass_growth_df.loc[condition, METHOD] = model_dict[condition].slim_optimize()


# ### Biomass Results

# In[ ]:


base_path = pathlib.Path("..") / "results"
biomass_growth_df.sort_index().to_csv(base_path / "biomass_raw.csv")

biomass_prop = biomass_growth_df.drop("WT", axis=0) / biomass_growth_df.loc["WT"]
biomass_prop.sort_index().to_csv(base_path / "biomass_proportion.csv")

