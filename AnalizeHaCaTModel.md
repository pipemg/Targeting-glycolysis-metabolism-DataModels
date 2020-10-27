

```python
import cobra
import math
#!wget http://www.ebi.ac.uk/biomodels-main/download?mid=MODEL1603150001 -O recon2.2.xml
Recon2 = cobra.io.read_sbml_model("Models/recon2.2.xml")
```



```python
hacat_model=cobra.io.read_sbml_model("Models/hgu133APlus2_hacat_model_0418_DEMEM6429_n5.sbml")

```


```python
#Biomass bounds
for rxex in hacat_model.exchanges:
    rxex.bounds=(-0.01,1)
    
### DMEM 6429 medium
#### AminoAcids

hacat_model.reactions.EX_ala_L_LPAREN_e_RPAREN_.lower_bound=-1
hacat_model.reactions.EX_arg_L_LPAREN_e_RPAREN_.lower_bound=-1
hacat_model.reactions.EX_asn_L_LPAREN_e_RPAREN_.lower_bound=-1
hacat_model.reactions.EX_asp_L_LPAREN_e_RPAREN_.lower_bound=-1
hacat_model.reactions.EX_cys_L_LPAREN_e_RPAREN_.lower_bound=-1
hacat_model.reactions.EX_gln_L_LPAREN_e_RPAREN_.lower_bound=-1
hacat_model.reactions.EX_glu_L_LPAREN_e_RPAREN_.lower_bound=-1
hacat_model.reactions.EX_gly_LPAREN_e_RPAREN_.lower_bound=-1
hacat_model.reactions.EX_his_L_LPAREN_e_RPAREN_.lower_bound=-1
hacat_model.reactions.EX_ile_L_LPAREN_e_RPAREN_.lower_bound=-1
hacat_model.reactions.EX_leu_L_LPAREN_e_RPAREN_.lower_bound=-1
hacat_model.reactions.EX_lys_L_LPAREN_e_RPAREN_.lower_bound=-1 
hacat_model.reactions.EX_met_L_LPAREN_e_RPAREN_.lower_bound=-1
hacat_model.reactions.EX_phe_L_LPAREN_e_RPAREN_.lower_bound=-1
hacat_model.reactions.EX_pro_L_LPAREN_e_RPAREN_.lower_bound=-1
hacat_model.reactions.EX_ser_L_LPAREN_e_RPAREN_.lower_bound=-1
hacat_model.reactions.EX_thr_L_LPAREN_e_RPAREN_.lower_bound=-1
hacat_model.reactions.EX_trp_L_LPAREN_e_RPAREN_.lower_bound=-1
hacat_model.reactions.EX_tyr_L_LPAREN_e_RPAREN_.lower_bound=-1
hacat_model.reactions.EX_val_L_LPAREN_e_RPAREN_.lower_bound=-1

### DMEM 6429 medium
#Carbon Sources
hacat_model.reactions.EX_glc_LPAREN_e_RPAREN_.lower_bound=-4.5
hacat_model.reactions.EX_pyr_LPAREN_e_RPAREN_.lower_bound= -1

### DMEM 6429 medium
#Minerals, vitamins, other
#hela_model.reactions.EX_ca2_LPAREN_e_RPAREN_.lower_bound=-1
hacat_model.reactions.EX_chol_LPAREN_e_RPAREN_.lower_bound=-1
#hela_model.reactions.EX_cl_LPAREN_e_RPAREN_.lower_bound=-1
#cobalt2
#cu2
hacat_model.reactions.EX_fe2_LPAREN_e_RPAREN_.lower_bound=-1
#hela_model.reactions.EX_fe3_LPAREN_e_RPAREN_.lower_bound=-1
hacat_model.reactions.EX_fol_LPAREN_e_RPAREN_.lower_single_gene_deletionbound=-1
#hela_model.reactions.EX_k_LPAREN_e_RPAREN_.lower_bound=-1
hacat_model.reactions.EX_h2o_LPAREN_e_RPAREN_.lower_bound=-10
#h2s
hacat_model.reactions.EX_inost_LPAREN_e_RPAREN_.lower_bound=-1
#magnesium
#manganase
#hela_model.reactions.EX_ncam_LPAREN_e_RPAREN_.lower_bound=-1
hacat_model.reactions.EX_pi_LPAREN_e_RPAREN_.lower_bound=-10
#hela_model.reactions.EX_pnto_R_LPAREN_e_RPAREN_.lower_bound=-1
#hela_model.reactions.EX_pydxn_LPAREN_e_RPAREN_.lower_bound=-1
hacat_model.reactions.EX_ribflv_LPAREN_e_RPAREN_.lower_bound=-1
hacat_model.reactions.EX_so4_LPAREN_e_RPAREN_.lower_bound=-1
#hela_model.reactions.EX_thm_LPAREN_e_RPAREN_.lower_bound=-1
#zn2

#EXTRA
hacat_model.reactions.EX_o2_LPAREN_e_RPAREN_.bounds=(-1000,0)
#hacat_model.reactions.EX_h2o2_LPAREN_e_RPAREN_.bounds=(0,100)
#hacat_model.reactions.EX_o2s_LPAREN_e_RPAREN_.bounds=(0,100)



hacat_model.reactions.biomass_reaction.lower_bound=0.008
hacat_model.reactions.biomass_reaction.upper_bound=0.037


fba_solution = hacat_model.optimize()
hacat_model.optimize()
```


    ---------------------------------------------------------------------------

    KeyError                                  Traceback (most recent call last)

    /opt/conda/lib/python3.5/site-packages/cobra/core/dictlist.py in __getattr__(self, attr)
        447         try:
    --> 448             return DictList.get_by_id(self, attr)
        449         except KeyError:


    /opt/conda/lib/python3.5/site-packages/cobra/core/dictlist.py in get_by_id(self, id)
         53         """return the element with a matching id"""
    ---> 54         return list.__getitem__(self, self._dict[id])
         55 


    KeyError: 'EX_glc_LPAREN_e_RPAREN_'

    
    During handling of the above exception, another exception occurred:


    AttributeError                            Traceback (most recent call last)

    <ipython-input-9-0f1e34e7386a> in <module>()
         29 ### DMEM 6429 medium
         30 #Carbon Sources
    ---> 31 hacat_model.reactions.EX_glc_LPAREN_e_RPAREN_.lower_bound=-4.5
         32 hacat_model.reactions.EX_pyr_LPAREN_e_RPAREN_.lower_bound= -1
         33 


    /opt/conda/lib/python3.5/site-packages/cobra/core/dictlist.py in __getattr__(self, attr)
        449         except KeyError:
        450             raise AttributeError("DictList has no attribute or entry %s" %
    --> 451                                  attr)
        452 
        453     def __dir__(self):


    AttributeError: DictList has no attribute or entry EX_glc_LPAREN_e_RPAREN_



```python
fba_solution.fluxes.to_csv('hacat_fba_fluxes_221117_DEMEM6429_n5.csv')
```


```python
pfba_solution = cobra.flux_analysis.pfba(hacat_model)
```


```python
pfba_solution.fluxes.to_csv('hacat_pfba_fluxes_221117_DEMEM6429_n5.csv')
```


```python
hacat_model.summary(fva=1.00)
```


```python
hacat_model.reactions.HMGCOASim.reaction
```


```python
hacat_model.metabolites.hmgcoa_m.summary(fva=1.0)
```


```python
import pandas
from time import time

import cobra.test

from cobra.flux_analysis import (
    single_gene_deletion, single_reaction_deletion, double_gene_deletion,
    double_reaction_deletion)

deletion_results=single_reaction_deletion(hacat_model, hacat_model.reactions)
```


```python
deletion_results[deletion_results.status != 'optimal']
```

Doubling time of human cancerous cell lines
16.2 Hrs 
http://bionumbers.hms.harvard.edu/bionumber.aspx?&id=100685&ver=21


```python
closedModel=hacat_model.copy()
closedModel.add_reaction(Recon2.reactions.DM_atp_c_)
closedModel.objective="DM_atp_c_"

for rx in closedModel.exchanges:
    rx.upper_bound= 100
    rx.lower_bound= -0.00000001
    

#Biomass
closedModel.reactions.biomass_other.lower_bound=0
closedModel.reactions.biomass_carbohydrate.lower_bound=0
closedModel.reactions.biomass_DNA.lower_bound=0
closedModel.reactions.biomass_lipid.lower_bound=0
closedModel.reactions.biomass_protein.lower_bound=0
closedModel.reactions.biomass_RNA.lower_bound=0
closedModel.reactions.biomass_reaction.lower_bound=0
closedModel.reactions.biomass_reaction.upper_bound=0


#closedModel.reactions.r1147.bounds=(0,0)
#closedModel.reactions.r2375.bounds=(0,0)
# test for max ATP hydrolysis flux from only o2 and the defined carbon
# source
# glucose aerobic

closedModel.reactions.EX_o2_LPAREN_e_RPAREN_.bounds=(-1000,0.0)
closedModel.reactions.EX_h2o_LPAREN_e_RPAREN_.bounds=(-1000,1000)
closedModel.reactions.EX_co2_LPAREN_e_RPAREN_.bounds=(-0.0,1000)
closedModel.reactions.EX_glc_LPAREN_e_RPAREN_.bounds=(-1,-1)
closedModel.optimize()
```


    ---------------------------------------------------------------------------

    KeyError                                  Traceback (most recent call last)

    /opt/conda/lib/python3.5/site-packages/cobra/core/dictlist.py in __getattr__(self, attr)
        447         try:
    --> 448             return DictList.get_by_id(self, attr)
        449         except KeyError:


    /opt/conda/lib/python3.5/site-packages/cobra/core/dictlist.py in get_by_id(self, id)
         53         """return the element with a matching id"""
    ---> 54         return list.__getitem__(self, self._dict[id])
         55 


    KeyError: 'EX_glc_LPAREN_e_RPAREN_'

    
    During handling of the above exception, another exception occurred:


    AttributeError                            Traceback (most recent call last)

    <ipython-input-8-b3b984bf7ac6> in <module>()
         28 closedModel.reactions.EX_h2o_LPAREN_e_RPAREN_.bounds=(-1000,1000)
         29 closedModel.reactions.EX_co2_LPAREN_e_RPAREN_.bounds=(-0.0,1000)
    ---> 30 closedModel.reactions.EX_glc_LPAREN_e_RPAREN_.bounds=(-1,-1)
         31 closedModel.optimize()


    /opt/conda/lib/python3.5/site-packages/cobra/core/dictlist.py in __getattr__(self, attr)
        449         except KeyError:
        450             raise AttributeError("DictList has no attribute or entry %s" %
    --> 451                                  attr)
        452 
        453     def __dir__(self):


    AttributeError: DictList has no attribute or entry EX_glc_LPAREN_e_RPAREN_



```python
cobra.flux_analysis.pfba(closedModel).fluxes.to_csv('hacat_pfba_fluxes_221117_DEMEM6429_n5_glucose_aerobic.csv')
```


```python
closedModel=hacat_model.copy()
closedModel.add_reaction(Recon2.reactions.DM_atp_c_)
closedModel.objective="DM_atp_c_"

for rx in closedModel.exchanges:
    rx.upper_bound= 100
    rx.lower_bound= -0.00000001
    

#Biomass
closedModel.reactions.biomass_other.lower_bound=0
closedModel.reactions.biomass_carbohydrate.lower_bound=0
closedModel.reactions.biomass_DNA.lower_bound=0
closedModel.reactions.biomass_lipid.lower_bound=0
closedModel.reactions.biomass_protein.lower_bound=0
closedModel.reactions.biomass_RNA.lower_bound=0
closedModel.reactions.biomass_reaction.lower_bound=0
closedModel.reactions.biomass_reaction.upper_bound=0

# test for max ATP hydrolysis flux from only o2 and the defined carbon
# source
# glucose anaerobic

closedModel.reactions.EX_o2_LPAREN_e_RPAREN_.bounds=(0,0)
closedModel.reactions.EX_h2o_LPAREN_e_RPAREN_.bounds=(-1000,1000)
closedModel.reactions.EX_co2_LPAREN_e_RPAREN_.bounds=(0.0,1000)
closedModel.reactions.EX_glc_LPAREN_e_RPAREN_.bounds=(-1,-1)
closedModel.optimize()
```


```python
cobra.flux_analysis.pfba(closedModel).fluxes.to_csv('hacat_pfba_fluxes_221117_DEMEM6429_n5_glucose_anaerobic.csv')
```


```python
closedModel=hacat_model.copy()
closedModel.add_reaction(Recon2.reactions.DM_atp_c_)
closedModel.objective="DM_atp_c_"

for rx in closedModel.exchanges:
    rx.upper_bound= 100
    rx.lower_bound= -0.00000001
    

#Biomass
closedModel.reactions.biomass_other.lower_bound=0
closedModel.reactions.biomass_carbohydrate.lower_bound=0
closedModel.reactions.biomass_DNA.lower_bound=0
closedModel.reactions.biomass_lipid.lower_bound=0
closedModel.reactions.biomass_protein.lower_bound=0
closedModel.reactions.biomass_RNA.lower_bound=0
closedModel.reactions.biomass_reaction.lower_bound=0
closedModel.reactions.biomass_reaction.upper_bound=0

# glutamine aerobic
closedModel.reactions.EX_o2_LPAREN_e_RPAREN_.bounds=(-1000,0)
closedModel.reactions.EX_h2o_LPAREN_e_RPAREN_.bounds=(-1000,1000)
closedModel.reactions.EX_co2_LPAREN_e_RPAREN_.bounds=(0,1000)
closedModel.reactions.EX_gln_L_LPAREN_e_RPAREN_.bounds=(-1,-1)

closedModel.optimize()
```


```python
cobra.flux_analysis.pfba(closedModel).fluxes.to_csv('hacat_pfba_fluxes_221117_DEMEM6429_n5_glutamine_aerobic.csv')
```


```python
closedModel=hacat_model.copy()
closedModel.add_reaction(Recon2.reactions.DM_atp_c_)
closedModel.objective="DM_atp_c_"

for rx in closedModel.exchanges:
    rx.upper_bound= 100
    rx.lower_bound= -0.00000001 #-0.00474 #this is the minimum value
    

#Biomass
closedModel.reactions.biomass_other.lower_bound=0
closedModel.reactions.biomass_carbohydrate.lower_bound=0
closedModel.reactions.biomass_DNA.lower_bound=0
closedModel.reactions.biomass_lipid.lower_bound=0
closedModel.reactions.biomass_protein.lower_bound=0
closedModel.reactions.biomass_RNA.lower_bound=0
closedModel.reactions.biomass_reaction.lower_bound=0
closedModel.reactions.biomass_reaction.upper_bound=0

# glutamine anaerobic
closedModel.reactions.EX_o2_LPAREN_e_RPAREN_.bounds=(0,0)
closedModel.reactions.EX_h2o_LPAREN_e_RPAREN_.bounds=(-1000,1000)
closedModel.reactions.EX_co2_LPAREN_e_RPAREN_.bounds=(0,1000)
closedModel.reactions.EX_gln_L_LPAREN_e_RPAREN_.bounds=(-1,-1)

closedModel.optimize()
```


```python
cobra.flux_analysis.pfba(closedModel).fluxes.to_csv('hacat_pfba_fluxes_221117_DEMEM6429_n5_glutamine_anaerobic.csv')
```


```python
closedModel=hacat_model.copy()
closedModel.add_reaction(Recon2.reactions.DM_atp_c_)
closedModel.objective="DM_atp_c_"

for rx in closedModel.exchanges:
    rx.upper_bound= 100
    rx.lower_bound= -0.00000001 #-0.00474
    
    
#Biomass
closedModel.reactions.biomass_other.lower_bound=0
closedModel.reactions.biomass_carbohydrate.lower_bound=0
closedModel.reactions.biomass_DNA.lower_bound=0
closedModel.reactions.biomass_lipid.lower_bound=0
closedModel.reactions.biomass_protein.lower_bound=0
closedModel.reactions.biomass_RNA.lower_bound=0
closedModel.reactions.biomass_reaction.lower_bound=0
closedModel.reactions.biomass_reaction.upper_bound=0


## fru aerobic
closedModel.reactions.EX_o2_LPAREN_e_RPAREN_.bounds=(-1000,0)
closedModel.reactions.EX_h2o_LPAREN_e_RPAREN_.bounds=(-1000,1000)
closedModel.reactions.EX_co2_LPAREN_e_RPAREN_.bounds=(0,1000)
closedModel.reactions.EX_fru_LPAREN_e_RPAREN_.bounds=(-1,-1)

closedModel.optimize()
```


```python
cobra.flux_analysis.pfba(closedModel).fluxes.to_csv('hacat_pfba_fluxes_221117_DEMEM6429_n5_fructose_aerobic.csv')
```


```python
closedModel=hacat_model.copy()
closedModel.add_reaction(Recon2.reactions.DM_atp_c_)
closedModel.objective="DM_atp_c_"

for rx in closedModel.exchanges:
    rx.upper_bound= 100
    rx.lower_bound= -0.00000001
    
    
#Biomass
closedModel.reactions.biomass_other.lower_bound=0
closedModel.reactions.biomass_carbohydrate.lower_bound=0
closedModel.reactions.biomass_DNA.lower_bound=0
closedModel.reactions.biomass_lipid.lower_bound=0
closedModel.reactions.biomass_protein.lower_bound=0
closedModel.reactions.biomass_RNA.lower_bound=0
closedModel.reactions.biomass_reaction.lower_bound=0
closedModel.reactions.biomass_reaction.upper_bound=0


## fru aerobic
closedModel.reactions.EX_o2_LPAREN_e_RPAREN_.bounds=(0,0)
closedModel.reactions.EX_h2o_LPAREN_e_RPAREN_.bounds=(-1000,1000)
closedModel.reactions.EX_co2_LPAREN_e_RPAREN_.bounds=(0,1000)
closedModel.reactions.EX_fru_LPAREN_e_RPAREN_.bounds=(-1,-1)

closedModel.optimize()
```


```python
cobra.flux_analysis.pfba(closedModel).fluxes.to_csv('hacat_pfba_fluxes_221117_DEMEM6429_n5_fructose_anaerobic.csv')
```


```python
import pandas
from time import time

import cobra.test

from cobra.flux_analysis import (
    single_gene_deletion, single_reaction_deletion, double_gene_deletion,
    double_reaction_deletion)

deletion_results=single_reaction_deletion(hacat_model, hacat_model.reactions)
```


```python
deletion_results[deletion_results["status"]!="optimal"]
```


```python
closedModel=hacat_model.copy()
for rx in closedModel.reactions:
    rx.bounds=(-0.001,0.001)
    fba=closedModel.optimize()
    fba2=hacat_model.optimize()
    if(fba.f != fba2.f/2 ):
        print(rx.id)
```


```python
#Glycolysis
parsFBA=
```


```python

```


```python

```
