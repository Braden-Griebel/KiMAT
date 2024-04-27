#!/usr/bin/env python
# coding: utf-8

# ## Imports

# In[5]:


# Standard Library imports
import os
import pathlib

# External Imports
import cobra
from cobra.flux_analysis import single_gene_deletion
import metworkpy
import numpy as np
import pandas as pd

# Local Imports


# ## Setup

# In[6]:


cobra.Configuration.solver = "cplex"
## Local Developement
# BASE_PATH = pathlib.Path("..")
## HPC
HOME = os.getenv("HOME")
BASE_PATH = pathlib.Path(HOME) / "projects" / "KiMAT"


# ## Data Preparation

# In[9]:


# Read in gene expression
gene_expression = pd.read_csv(BASE_PATH / 'data' / 'gene-expression' / 'kinase_RNA_seq.csv', index_col=0).drop( "GrowthPhase", axis=1)
# Convert the RPKM measurements into TPM
gene_expression = metworkpy.utils.rpkm_to_tpm(gene_expression)
# Read in quantile normalized gene expression
normalized_expression = pd.read_csv(BASE_PATH / "data" / "quant_normed" / "xprs_norm.tsv", sep="\t", index_col=0).transpose()


# In[30]:


# Separate the expression data into biological replicates
expr_dict = {
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
norm_expr_dict = {
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


for condition in expr_dict:
    expr_dict[condition] = gene_expression.loc[gene_expression.index.str.startswith(condition)]
    norm_expr_dict[condition] = normalized_expression.loc[normalized_expression.index.str.startswith(condition)]


# ## Read in Base Model

# In[65]:


iek1011 = metworkpy.read_model(BASE_PATH / "data" / "Models" / "iEK1011_m7H10_media.json")
for method in ["simple_bounds", "subset_model", "fva_model", "milp_model"]:
    metworkpy.write_model(iek1011, 
                          BASE_PATH / "results" / "diff_expr_models" / method / "iek1011_WT_model.json")


# ## Convert Gene Expression into Trinarized Reactions

# In[33]:


rxn_weight_dict = {
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
norm_rxn_weight_dict = {
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

diff_rxn_weight_dict = {
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

wt_expression = gene_expression.loc[gene_expression.index.str.startswith("WT")].median(axis=0)+1e-10

for condition in rxn_weight_dict:
    p = BASE_PATH / 'cache' / 'rxn_weights' / f'{condition}_rxn_weights.json'
    norm_p = BASE_PATH / 'cache' /'norm_rxn_weights' / f'{condition}_rxn_weights.json'
    if p.exists():
        rxn_weight_dict[condition] = pd.read_json(p, typ = "series")
    else: 
        gene_weights = metworkpy.utils.expr_to_weights(expression = expr_dict[condition],
                                                      quantile = 0.15, 
                                                      subset = iek1011.genes.list_attr("id"),
                                                      sample_axis=0)
        rxn_weight_dict[condition] = metworkpy.parse.gene_to_rxn_weights(iek1011, gene_weights)
        rxn_weight_dict[condition].to_json(p)
    if norm_p.exists():
        norm_rxn_weight_dict[condition] = pd.read_json(norm_p, typ = "series")
    else: 
        norm_gene_weights = metworkpy.utils.expr_to_weights(expression = norm_expr_dict[condition],
                                                      quantile = 0.15, 
                                                      subset = iek1011.genes.list_attr("id"),
                                                      sample_axis=0)
        norm_rxn_weight_dict[condition] = metworkpy.parse.gene_to_rxn_weights(iek1011, norm_gene_weights)
        norm_rxn_weight_dict[condition].to_json(norm_p)
    if condition != "WT":
        diff_p = BASE_PATH / "cache" / "diff_rxn_weights" / f"{condition}_rxn_weights.json"
        if diff_p.exists():
            diff_rxn_weight_dict[condition] = pd.read_json(diff_p, typ="series")
        else:
            log2fc = np.log2(gene_expression.loc[gene_expression.index.str.startswith(condition)].median(axis=0).div(wt_expression)+1e-10)
            diff_gene_weights = log2fc.apply(lambda x: -1 if x<=-1 else (1 if x>=1 else 0))
            diff_rxn_weight_dict[condition] = metworkpy.parse.gene_to_rxn_weights(iek1011, diff_gene_weights)
            diff_rxn_weight_dict[condition].to_json(diff_p)


# ## IMAT

# In[37]:


biomass_growth_df = pd.DataFrame(np.NaN, index = rxn_weight_dict.keys(), 
                                 columns=["simple_bounds", "subset_model", "fva_model", "milp_model"],
                                dtype="float")
norm_biomass_growth_df = pd.DataFrame(np.NaN, index = norm_rxn_weight_dict.keys(), 
                                 columns=["simple_bounds", "subset_model", "fva_model", "milp_model"],
                                dtype="float")
diff_biomass_growth_df = pd.DataFrame(np.NaN, index = norm_rxn_weight_dict.keys(), 
                                     columns = ["simple_bounds", "subset_model", "fva_model", "milp_model"],
                                     dtype="float")
wt_objective = iek1011.slim_optimize()
diff_biomass_growth_df.loc['WT'] = wt_objective


# In[29]:


def imat_per_method(method, weights_dict, results_path, biomass_growth_df, skip_list=None):
    print(f"METHOD: {method}")
    print("**********************")
    for condition in weights_dict:
        print(f"\tCondition: {condition}")
        print("\t--------------")
        if condition in skip_list:
            print("SKIPPED")
            continue
        p = results_path / method / f"iek1011_{condition}_model.json"
        if p.exists():
            model = metworkpy.read_model(str(p))
        else:
            try:
                model = metworkpy.imat.generate_model(model = iek1011.copy(), 
                                                      rxn_weights=weights_dict[condition],
                                                     method = method)
                metworkpy.write_model(model, model_path=p)
            except Exception as exp:
                continue
        try:
            biomass_growth_df.loc[condition, method] = model.slim_optimize()
        except Exception as exp:
            pass


# ### Simple Bounds

# In[ ]:


imat_per_method(method="simple_bounds", weights_dict=rxn_weight_dict, 
               results_path=BASE_PATH / "results" / "expr_models", 
                biomass_growth_df=biomass_growth_df)
imat_per_method(method="simple_bounds", weights_dict=norm_rxn_weight_dict, 
               results_path=BASE_PATH / "results" / "norm_expr_models", 
                biomass_growth_df=norm_biomass_growth_df)
imat_per_method(method="simple_bounds", weights_dict = diff_rxn_weight_dict, 
               results_path = BASE_PATH / "results" / "diff_expr_models",
               biomass_growth_df = diff_biomass_growth_df)


# ### Subset Model

# In[ ]:


imat_per_method(method="subset_model", weights_dict=rxn_weight_dict, 
               results_path=BASE_PATH / "results" / "expr_models", 
                biomass_growth_df=biomass_growth_df)
imat_per_method(method="subset_model", weights_dict=norm_rxn_weight_dict, 
               results_path=BASE_PATH / "results" / "norm_expr_models", 
                biomass_growth_df=norm_biomass_growth_df)
imat_per_method(method="subset_model", weights_dict = diff_rxn_weight_dict, 
               results_path = BASE_PATH / "results" / "diff_expr_models",
               biomass_growth_df = diff_biomass_growth_df)


# ### FVA Model

# In[ ]:


imat_per_method(method="fva_model", weights_dict=rxn_weight_dict, 
               results_path=BASE_PATH / "results" / "expr_models", 
                biomass_growth_df=biomass_growth_df)
imat_per_method(method="fva_model", weights_dict=norm_rxn_weight_dict, 
               results_path=BASE_PATH / "results" / "norm_expr_models", 
                biomass_growth_df=norm_biomass_growth_df)
imat_per_method(method="fva_model", weights_dict = diff_rxn_weight_dict, 
               results_path = BASE_PATH / "results" / "diff_expr_models",
               biomass_growth_df = diff_biomass_growth_df)


# ### MILP Model

# In[ ]:


imat_per_method(method="milp_model", weights_dict=rxn_weight_dict, 
               results_path=BASE_PATH / "results" / "expr_models", 
                biomass_growth_df=biomass_growth_df)
imat_per_method(method="milp_model", weights_dict=norm_rxn_weight_dict, 
               results_path=BASE_PATH / "results" / "norm_expr_models", 
                biomass_growth_df=norm_biomass_growth_df)
imat_per_method(method="milp_model", weights_dict = diff_rxn_weight_dict, 
               results_path = BASE_PATH / "results" / "diff_expr_models",
               biomass_growth_df = diff_biomass_growth_df)


# ### Save Biomass DataFrames

# In[ ]:


# Basic results
biomass_growth_df.to_csv(BASE_PATH / "results" / "expr_models" / "biomass_growth.csv")
biomass_growth_df[~(biomass_growth_df.index=="WT")].div(biomass_growth_df.loc["WT"], axis=1).to_csv(BASE_PATH / "results" / "expr_models" / "biomass_growth_prop.csv")
# Qnorm results
norm_biomass_growth_df.to_csv(BASE_PATH / "results" / "norm_expr_models" / "norm_biomass_growth.csv")
norm_biomass_growth_df[~(norm_biomass_growth_df.index=="WT")].div(norm_biomass_growth_df.loc["WT"], axis=1).to_csv(
    BASE_PATH / "results" / "norm_expr_models" / "norm_biomass_growth_prop.csv")
# Differential Expression Results
diff_biomass_growth_df.to_csv(BASE_PATH / "results" / "diff_expr_models" / "diff_biomass_growth.csv")
diff_biomass_growth_df[~(diff_biomass_growth_df.index=="WT")].div(diff_biomass_growth_df.loc["WT"], axis=1).to_csv(
    BASE_PATH / "results" / "diff_expr_models" / "diff_biomass_growth_prop.csv")


# ## Essentiality Analysis

# In[3]:


def essentiality_analysis(model_dir: pathlib.Path, 
                          cond_list:list[str], 
                          wt:str,
                          out_dir: pathlib.Path)->None: 
    wt_model = metworkpy.read_model(str(model_dir / f"iek1011_{wt}_model.json"))
    flux_res_df = pd.DataFrame(np.NaN, 
                             dtype="float", 
                             index = wt_model.genes.list_attr("id"),
                             columns = [wt]+cond_list
                            )
    prop_res_df = pd.DataFrame(np.NaN, 
                             dtype="float", 
                             index = wt_model.genes.list_attr("id"),
                             columns = [wt]+cond_list
                            )
    ess_res_df = pd.DataFrame(np.NaN, 
                             dtype="float", 
                             index = wt_model.genes.list_attr("id"),
                             columns = [wt]+cond_list
                            )
    for cond in [wt]+cond_list:
        print("\t-----------------")
        print(f"\tCondition: {cond}")
        print("\t-----------------")
        cond_model = metworkpy.read_model(str(model_dir / f"iek1011_{cond}_model.json"))
        try:
            max_obj = cond_model.slim_optimize()
            ko_res = single_gene_deletion(cond_model)
            ko_res["ids"] = ko_res["ids"].apply(lambda x: next(iter(x)))
            ko_res = ko_res.set_index("ids")
            flux_res_df.loc[ko_res.index, cond] = ko_res["growth"]
            prop_res_df.loc[flux_res_df.index, cond] = flux_res_df.loc[:, cond] / max_obj
            ess_res_df.loc[(prop_res_df.loc[:,cond]<=5e-2).index, cond] = 1
            ess_res_df.loc[(~(prop_res_df.loc[:,cond]>5e-2)).index, cond] = 0
        except Exception as exp:
            continue
    flux_res_df.to_csv(out_dir / "ko_flux.csv")
    prop_res_df.to_csv(out_dir / "ko_prop.csv")
    ess_res_df.to_csv(out_dir / "ko_essentiality.csv")        


# In[4]:


# Normal Expression Essentiality Analysis
EXPRESSION_METHOD = "expr"
print(f"Starting: {EXPRESSION_METHOD}")
for method in ["simple_bounds", "subset_model", "fva_model", "milp_model"]:
    essentiality_analysis(
        model_dir = BASE_PATH / "results" / f"{EXPRESSION_METHOD}_models" / method,
        cond_list = ["PknB_KD",
                     "PknB_IND",
                     "PknD_KO",
                     "PknD_IND",
                     "PknE_KO",
                     "PknE_IND",
                     "PknF_KO",
                     "PknF_IND",
                     "PknG_KO",
                     "PknG_IND",
                     "PknH_KO",
                     "PknH_IND",
                     "PknJ_KO",
                     "PknJ_IND",
                     "PknK_KO",
                     "PknK_IND",
                     "PknL_KO",
                     "PknL_IND",],
        wt = "WT",
        out_dir = BASE_PATH / "results" / "essentiality" /  f"{EXPRESSION_METHOD}_essentiality" / method
       )
# Quantile Normalized Expression Essentiality Analysis
EXPRESSION_METHOD = "norm_expr"
print(f"Starting: {EXPRESSION_METHOD}")
for method in ["simple_bounds", "subset_model", "fva_model", "milp_model"]:
    essentiality_analysis(
        model_dir = BASE_PATH / "results" / f"{EXPRESSION_METHOD}_models" / method,
        cond_list = ["PknB_KD",
                     "PknB_IND",
                     "PknD_KO",
                     "PknD_IND",
                     "PknE_KO",
                     "PknE_IND",
                     "PknF_KO",
                     "PknF_IND",
                     "PknG_KO",
                     "PknG_IND",
                     "PknH_KO",
                     "PknH_IND",
                     "PknJ_KO",
                     "PknJ_IND",
                     "PknK_KO",
                     "PknK_IND",
                     "PknL_KO",
                     "PknL_IND",],
        wt = "WT",
        out_dir = BASE_PATH / "results" / "essentiality" /  f"{EXPRESSION_METHOD}_essentiality" / method
       )
# Differential Expression Essentiality Analysis
EXPRESSION_METHOD = "diff_expr"
print(f"Starting: {EXPRESSION_METHOD}")
for method in ["simple_bounds", "subset_model", "fva_model", "milp_model"]:
    essentiality_analysis(
        model_dir = BASE_PATH / "results" / f"{EXPRESSION_METHOD}_models" / method,
        cond_list = ["PknB_KD",
                     "PknB_IND",
                     "PknD_KO",
                     "PknD_IND",
                     "PknE_KO",
                     "PknE_IND",
                     "PknF_KO",
                     "PknF_IND",
                     "PknG_KO",
                     "PknG_IND",
                     "PknH_KO",
                     "PknH_IND",
                     "PknJ_KO",
                     "PknJ_IND",
                     "PknK_KO",
                     "PknK_IND",
                     "PknL_KO",
                     "PknL_IND",],
        wt = "WT",
        out_dir = BASE_PATH / "results" / "essentiality" /  f"{EXPRESSION_METHOD}_essentiality" / method
       )


# In[ ]:




