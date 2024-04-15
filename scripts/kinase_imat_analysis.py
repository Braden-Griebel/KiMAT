# Standard Library imports
import os
import pathlib

# External Imports
import cobra
import metworkpy
import numpy as np
import pandas as pd

# Local Imports

# SETUP
cobra.Configuration.solver = "cplex"
HOME = os.getenv("HOME")
BASE_PATH = pathlib.Path(HOME) / "projects" / "KiMAT" 

gene_expression = pd.read_csv(BASE_PATH / 'data' / 'gene-expression' / 'kinase_RNA_seq.csv', index_col=0).drop( "GrowthPhase", axis=1)
# Convert the RPKM measurements into TPM
gene_expression = metworkpy.utils.rpkm_to_tpm(gene_expression)

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
    
# Convert gene expression into qualitative weights
for condition in condition_dict:
    condition_dict[condition] = metworkpy.utils.expr_to_weights(condition_dict[condition], quantile = 0.15, 
                                                                aggregator=np.median, sample_axis=0)
                                                                
# Convert qualitative weights into trinarized reactions

# Load the iEK1011 model
iek1011 = metworkpy.read_model(BASE_PATH / 'data' / 'Models' / 'iEK1011_m7H10_media.json')

# Convert the gene expression data into trinarized reaction data
# WARNING: This step may take a while, it is not implemented very efficiently currently
for condition in condition_dict:
    p = BASE_PATH / 'cache' / f'{condition}_rxn_weights.json'
    if p.exists():
        condition_dict[condition] = pd.read_json(p, typ = "series")
    else:
        condition_dict[condition] = metworkpy.parse.gene_to_rxn_weights(iek1011, condition_dict[condition])
        condition_dict[condition].to_json(p)
        
biomass_growth_df = pd.DataFrame(0., index = condition_dict.keys(), columns=["simple_bounds", "subset_model", "fva_model", "milp_model"])

                                                              
                                                                                                                            
# Function for performing the IMAT analysis in each condition

def imat_per_method(method, condition_dict, base_path, biomass_growth_df):
    print(f"Method: {method}")
    print("***************************")
    for condition in condition_dict:
        print(condition)
        print("---------------")
        p = base_path / 'results' / method / f"iek1011_{condition}_model.json"
        if p.exists():
            model = metworkpy.read_model(str(p))
        else:
            try:
                model = metworkpy.imat.generate_model(model = iek1011.copy(), rxn_weights=condition_dict[condition], method=method)
                metworkpy.write_model(model, model_path=p)
            except:
                continue
        biomass_growth_df.loc[condition, method] = model.slim_optimize()
    
        
imat_per_method(method = "simple_bounds", condition_dict=condition_dict, base_path=BASE_PATH, biomass_growth_df=biomass_growth_df)
imat_per_method(method = "subset_model", condition_dict=condition_dict, base_path=BASE_PATH, biomass_growth_df=biomass_growth_df)
imat_per_method(method = "fva_model", condition_dict=condition_dict, base_path=BASE_PATH, biomass_growth_df=biomass_growth_df)
imat_per_method(method = "milp_model", condition_dict=condition_dict, base_path=BASE_PATH, biomass_growth_df=biomass_growth_df)        

biomass_growth_df.sort_index().to_csv(BASE_PATH / "results" / "biomass_raw.csv")

biomass_prop = biomass_growth_df.drop("WT", axis=0) / biomass_growth_df.loc["WT"]
biomass_prop.sort_index().to_csv(BASE_PATH/ "results" / "biomass_proportion.csv")                                                     