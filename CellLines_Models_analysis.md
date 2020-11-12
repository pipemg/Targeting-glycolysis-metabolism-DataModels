
### Reactions, Genes and confidences

In order to start a reconstruction CORDA requires you to assign a confidence score to each reaction in your base model. This can be done by a variety of methods, and even by hand, but the most common way is to assign confidence based on proteome or gene expression data.
CORDA manages a total of 5 confidence levels:

+ -1 for reactions that are not present and should not be included in the model
+ 0 for reactions with unknown confidence which may be included in the model if necessary
+ 1 for low confidence reactions that should be included if necessary
+ 2 for medium confidence reactions that should be included if necessary
+ 3 for high confidence reactions that must be included if possible in any way

The most tedious step here is usaully mapping the confidence for genes or proteins to the distinct reactions. Many of the larger models come with gene-reaction rules in the form
gene1 and gene2 or (gene3 and gene4)
and the individual confidence values have to be mapped from the gene confidence levels. Here "and" is evaluated by the minimum confidence and "or" by the maximum confidence. The Python package includes a handy function to do this for you automatically in a safe manner. For that you will require the gene-reaction rule (Recon 1 and 2 include them in their model for instance) and a dictionary mapping genes/proteins to their confidence values. 


### Load the ubiquitin scores of the microarrays

#### Definitions
   If we have a matrix M where rows are unique genes and columns are samples (cancer samples for example) we define g[i,j] as the gene i in the sample j where g[i,j] in {0,1}  where 0 means a gene that is "Off" in the sample and 1 represents "On".
   Whith this we calculate the percentil (0 to 1) of each row (genes) and this will be the ubiquitin score for gene i in the model, this score represents the number of times a gene is considered on or off in that specific data, and we do this for each model we whant to do. 


```python
import cobra
#!wget http://www.ebi.ac.uk/biomodels-main/download?mid=MODEL1603150001 -O recon2.2.xml
Recon2 = cobra.io.read_sbml_model("Models/recon2.2.xml")
#Recon204 = cobra.io.load_matlab_model("Models/Recon2.v04.mat")
Recon2.reactions.OIVD3m.gene_reaction_rule='HGNC:2698 and HGNC:987 and HGNC:986 and HGNC:2898 or HGNC:2698 and HGNC:2898 and HGNC:987 and HGNC:986'
Recon2.reactions.OIVD1m.gene_reaction_rule='HGNC:2698 and HGNC:987 and HGNC:986 and HGNC:2898 or HGNC:2698 and HGNC:2898 and HGNC:987 and HGNC:986'
Recon2.reactions.OIVD2m.gene_reaction_rule='HGNC:2698 and HGNC:987 and HGNC:986 and HGNC:2898 or HGNC:2698 and HGNC:2898 and HGNC:987 and HGNC:986'
```

    cobra/io/sbml.py:235 [1;31mUserWarning[0m: M_h_x appears as a reactant and product FAOXC220200x



```python
for rxex in Recon2.exchanges:
    rxex.bounds=(0.001,10)

### DMEM 6429

Recon2.reactions.EX_ala_L_LPAREN_e_RPAREN_.lower_bound=-0.01
Recon2.reactions.EX_arg_L_LPAREN_e_RPAREN_.lower_bound=-0.084
Recon2.reactions.EX_asn_L_LPAREN_e_RPAREN_.lower_bound=-0.01
Recon2.reactions.EX_asp_L_LPAREN_e_RPAREN_.lower_bound=-0.01
Recon2.reactions.EX_ca2_LPAREN_e_RPAREN_.lower_bound=-1
Recon2.reactions.EX_chol_LPAREN_e_RPAREN_.lower_bound=-0.05
Recon2.reactions.EX_cl_LPAREN_e_RPAREN_.lower_bound=-1
Recon2.reactions.EX_cys_L_LPAREN_e_RPAREN_.lower_bound=-0.0626
Recon2.reactions.EX_fe2_LPAREN_e_RPAREN_.lower_bound=-1
Recon2.reactions.EX_fe3_LPAREN_e_RPAREN_.lower_bound=-1
Recon2.reactions.EX_fol_LPAREN_e_RPAREN_.lower_bound=-0.05
Recon2.reactions.EX_glc_LPAREN_e_RPAREN_.bounds=(-4.5,0)
Recon2.reactions.EX_gln_L_LPAREN_e_RPAREN_.lower_bound=-0.584
Recon2.reactions.EX_glu_L_LPAREN_e_RPAREN_.lower_bound=-0.01
Recon2.reactions.EX_gly_LPAREN_e_RPAREN_.lower_bound=-0.03
Recon2.reactions.EX_h_LPAREN_e_RPAREN_.lower_bound=-10
Recon2.reactions.EX_h2o_LPAREN_e_RPAREN_.lower_bound=-10
Recon2.reactions.EX_his_L_LPAREN_e_RPAREN_.lower_bound=-0.042
Recon2.reactions.EX_ile_L_LPAREN_e_RPAREN_.lower_bound=-0.105
Recon2.reactions.EX_inost_LPAREN_e_RPAREN_.lower_bound=-0.05
Recon2.reactions.EX_k_LPAREN_e_RPAREN_.lower_bound=-1
Recon2.reactions.EX_leu_L_LPAREN_e_RPAREN_.lower_bound=-0.105
Recon2.reactions.EX_lys_L_LPAREN_e_RPAREN_.lower_bound=-0.146
Recon2.reactions.EX_met_L_LPAREN_e_RPAREN_.lower_bound=-0.03
Recon2.reactions.EX_na1_LPAREN_e_RPAREN_.lower_bound=-1
Recon2.reactions.EX_ncam_LPAREN_e_RPAREN_.lower_bound=-0.05
Recon2.reactions.EX_o2_LPAREN_e_RPAREN_.lower_bound=-10
Recon2.reactions.EX_phe_L_LPAREN_e_RPAREN_.lower_bound=-0.066
Recon2.reactions.EX_pi_LPAREN_e_RPAREN_.lower_bound=-1
Recon2.reactions.EX_pnto_R_LPAREN_e_RPAREN_.lower_bound=-0.05
Recon2.reactions.EX_pro_L_LPAREN_e_RPAREN_.lower_bound=-0.01
Recon2.reactions.EX_pydxn_LPAREN_e_RPAREN_.lower_bound=-0.05
Recon2.reactions.EX_pyr_LPAREN_e_RPAREN_.lower_bound=-0.11
Recon2.reactions.EX_ribflv_LPAREN_e_RPAREN_.lower_bound=-0.05
Recon2.reactions.EX_ser_L_LPAREN_e_RPAREN_.lower_bound=-0.042
Recon2.reactions.EX_so4_LPAREN_e_RPAREN_.lower_bound=-1
Recon2.reactions.EX_thm_LPAREN_e_RPAREN_.lower_bound=-0.05
Recon2.reactions.EX_thr_L_LPAREN_e_RPAREN_.lower_bound=-0.095
Recon2.reactions.EX_trp_L_LPAREN_e_RPAREN_.lower_bound=-0.016
Recon2.reactions.EX_tyr_L_LPAREN_e_RPAREN_.lower_bound=-0.10379
Recon2.reactions.EX_val_L_LPAREN_e_RPAREN_.lower_bound=-0.094
Recon2.reactions.EX_lac_L_LPAREN_e_RPAREN_.bounds=(-1,10)
Recon2.reactions.EX_lac_D_LPAREN_e_RPAREN_.bounds=(-0.001,0.001)

Recon2.reactions.biomass_reaction.bounds=(0.008, 0.05)

```


```python
Recon2.reactions.get_by_id("biomass_DNA").reaction
```




    '0.941642857142857 datp_n + 0.674428571428572 dctp_n + 0.707 dgtp_n + 0.935071428571429 dttp_n --> biomass_DNA_c'




```python
import pandas as pd
ub_scores_matrix =pd.read_csv("Binary_CellLines_RECON2_2.tsv",sep=",",index_col=0,header=0)
confidence_scores_matrix=ub_scores_matrix.copy()
```


```python
for model in ub_scores_matrix.columns:
    #High Confidence rate 
       
    HC=ub_scores_matrix[model][ub_scores_matrix[model] >= 0.75].index.tolist()
    confidence_scores_matrix[model][HC]=3
        
    #Medium Confidence Rate 
    MC=ub_scores_matrix[model][(ub_scores_matrix[model] < 0.75) & (ub_scores_matrix[model] >=0.5) ].index.tolist()
    confidence_scores_matrix[model][MC]=2
        
    #Low
    LC=ub_scores_matrix[model][(ub_scores_matrix[model] < 0.5) & (ub_scores_matrix[model] >=0.25) ].index.tolist()
    confidence_scores_matrix[model][LC]=1
    
    #unknown
    ZC=ub_scores_matrix[model][(ub_scores_matrix[model] <.25) & (ub_scores_matrix[model] > 0) ].index.tolist()
    confidence_scores_matrix[model][ZC]=0
    
    #Negative
    NC=ub_scores_matrix[model][(ub_scores_matrix[model] == 0)].index.tolist()
    confidence_scores_matrix[model][NC]=-1
 
```


```python
from corda import reaction_confidence


conf_HeLa = {}
conf_keratinocytes= {}

for r in Recon2.reactions:
    if(r.gene_reaction_rule!=''):
        conf_HeLa[r.id]=reaction_confidence(r.gene_reaction_rule,confidence_scores_matrix["MaxHeLa"])
        conf_keratinocytes[r.id]=reaction_confidence(r.gene_reaction_rule,confidence_scores_matrix["Maxkeratinocytes"])
    else:
        conf_HeLa[r.id]=1
        conf_keratinocytes[r.id]=1
```


```python
%matplotlib inline
import pandas as pd

df=pd.DataFrame({'Cancer': conf_CancerBiopsy})
df.Cancer.value_counts()

```


    ---------------------------------------------------------------------------

    NameError                                 Traceback (most recent call last)

    <ipython-input-7-ddcefe5cc071> in <module>()
          2 import pandas as pd
          3 
    ----> 4 df=pd.DataFrame({'Cancer': conf_CancerBiopsy})
          5 df.Cancer.value_counts()


    NameError: name 'conf_CancerBiopsy' is not defined



```python
conf_HeLa["EX_arg_L_LPAREN_e_RPAREN_"]=3
conf_HeLa["EX_cys_L_LPAREN_e_RPAREN_"]=3
conf_HeLa["EX_gln_L_LPAREN_e_RPAREN_"]=3
conf_HeLa["EX_gly_LPAREN_e_RPAREN_"]=3
conf_HeLa["EX_his_L_LPAREN_e_RPAREN_"]=3
conf_HeLa["EX_ile_L_LPAREN_e_RPAREN_"]=3
conf_HeLa["EX_leu_L_LPAREN_e_RPAREN_"]=3
conf_HeLa["EX_lys_L_LPAREN_e_RPAREN_"]=3
conf_HeLa["EX_met_L_LPAREN_e_RPAREN_"]=3
conf_HeLa["EX_phe_L_LPAREN_e_RPAREN_"]=3
conf_HeLa["EX_ser_L_LPAREN_e_RPAREN_"]=3
conf_HeLa["EX_thr_L_LPAREN_e_RPAREN_"]=3
conf_HeLa["EX_trp_L_LPAREN_e_RPAREN_"]=3
conf_HeLa["EX_tyr_L_LPAREN_e_RPAREN_"]=3
conf_HeLa["EX_thm_LPAREN_e_RPAREN_"]=3
conf_HeLa["EX_val_L_LPAREN_e_RPAREN_"]=3
conf_HeLa["EX_glc_LPAREN_e_RPAREN_"]=3
conf_HeLa["EX_pyr_LPAREN_e_RPAREN_"]=3
conf_HeLa["EX_ca2_LPAREN_e_RPAREN_"]=3
conf_HeLa["EX_cl_LPAREN_e_RPAREN_"]=3
conf_HeLa["EX_k_LPAREN_e_RPAREN_"]=3
conf_HeLa["EX_na1_LPAREN_e_RPAREN_"]=3
conf_HeLa["EX_fe2_LPAREN_e_RPAREN_"]=3
conf_HeLa["EX_fe3_LPAREN_e_RPAREN_"]=3
conf_HeLa["EX_so4_LPAREN_e_RPAREN_"]=3
conf_HeLa["EX_pi_LPAREN_e_RPAREN_"]=3
conf_HeLa["EX_h_LPAREN_e_RPAREN_"]=3
conf_HeLa["EX_h2o_LPAREN_e_RPAREN_"]=3
conf_HeLa["EX_o2_LPAREN_e_RPAREN_"]=3
conf_HeLa["EX_co2_LPAREN_e_RPAREN_"]=3
conf_HeLa["EX_lac_L_LPAREN_e_RPAREN_"]=3
conf_HeLa["EX_lac_D_LPAREN_e_RPAREN_"]=-1
conf_HeLa["EX_pro_L_LPAREN_e_RPAREN_"]=3
conf_HeLa["EX_ala_L_LPAREN_e_RPAREN_"]=3
conf_HeLa["EX_asn_L_LPAREN_e_RPAREN_"]=3
conf_HeLa["EX_asp_L_LPAREN_e_RPAREN_"]=3
conf_HeLa["EX_glu_L_LPAREN_e_RPAREN_"]=3
conf_HeLa["EX_chol_LPAREN_e_RPAREN_"]=3
conf_HeLa["EX_fol_LPAREN_e_RPAREN_"]=3
conf_HeLa["EX_inost_LPAREN_e_RPAREN_"]=3
conf_HeLa["EX_ncam_LPAREN_e_RPAREN_"]=3
conf_HeLa["EX_pnto_R_LPAREN_e_RPAREN_"]=3
conf_HeLa["EX_pydxn_LPAREN_e_RPAREN_"]=3
conf_HeLa["EX_ribflv_LPAREN_e_RPAREN_"]=3



```


```python
conf_keratinocytes["EX_arg_L_LPAREN_e_RPAREN_"]=3
conf_keratinocytes["EX_cys_L_LPAREN_e_RPAREN_"]=3
conf_keratinocytes["EX_gln_L_LPAREN_e_RPAREN_"]=3
conf_keratinocytes["EX_gly_LPAREN_e_RPAREN_"]=3
conf_keratinocytes["EX_his_L_LPAREN_e_RPAREN_"]=3
conf_keratinocytes["EX_ile_L_LPAREN_e_RPAREN_"]=3
conf_keratinocytes["EX_leu_L_LPAREN_e_RPAREN_"]=3
conf_keratinocytes["EX_lys_L_LPAREN_e_RPAREN_"]=3
conf_keratinocytes["EX_met_L_LPAREN_e_RPAREN_"]=3
conf_keratinocytes["EX_phe_L_LPAREN_e_RPAREN_"]=3
conf_keratinocytes["EX_ser_L_LPAREN_e_RPAREN_"]=3
conf_keratinocytes["EX_thr_L_LPAREN_e_RPAREN_"]=3
conf_keratinocytes["EX_trp_L_LPAREN_e_RPAREN_"]=3
conf_keratinocytes["EX_tyr_L_LPAREN_e_RPAREN_"]=3
conf_keratinocytes["EX_thm_LPAREN_e_RPAREN_"]=3
conf_keratinocytes["EX_val_L_LPAREN_e_RPAREN_"]=3
conf_keratinocytes["EX_glc_LPAREN_e_RPAREN_"]=3
conf_keratinocytes["EX_pyr_LPAREN_e_RPAREN_"]=3
conf_keratinocytes["EX_ca2_LPAREN_e_RPAREN_"]=3
conf_keratinocytes["EX_cl_LPAREN_e_RPAREN_"]=3
conf_keratinocytes["EX_k_LPAREN_e_RPAREN_"]=3
conf_keratinocytes["EX_na1_LPAREN_e_RPAREN_"]=3
conf_keratinocytes["EX_fe2_LPAREN_e_RPAREN_"]=3
conf_keratinocytes["EX_fe3_LPAREN_e_RPAREN_"]=3
conf_keratinocytes["EX_so4_LPAREN_e_RPAREN_"]=3
conf_keratinocytes["EX_pi_LPAREN_e_RPAREN_"]=3
conf_keratinocytes["EX_h_LPAREN_e_RPAREN_"]=3
conf_keratinocytes["EX_h2o_LPAREN_e_RPAREN_"]=3
conf_keratinocytes["EX_o2_LPAREN_e_RPAREN_"]=3
conf_keratinocytes["EX_co2_LPAREN_e_RPAREN_"]=3
conf_keratinocytes["EX_lac_L_LPAREN_e_RPAREN_"]=3
conf_keratinocytes["EX_lac_D_LPAREN_e_RPAREN_"]=-1
conf_keratinocytes["EX_pro_L_LPAREN_e_RPAREN_"]=3
conf_keratinocytes["EX_ala_L_LPAREN_e_RPAREN_"]=3
conf_keratinocytes["EX_asn_L_LPAREN_e_RPAREN_"]=3
conf_keratinocytes["EX_asp_L_LPAREN_e_RPAREN_"]=3
conf_keratinocytes["EX_glu_L_LPAREN_e_RPAREN_"]=3
conf_keratinocytes["EX_chol_LPAREN_e_RPAREN_"]=3
conf_keratinocytes["EX_fol_LPAREN_e_RPAREN_"]=3
conf_keratinocytes["EX_inost_LPAREN_e_RPAREN_"]=3
conf_keratinocytes["EX_ncam_LPAREN_e_RPAREN_"]=3
conf_keratinocytes["EX_pnto_R_LPAREN_e_RPAREN_"]=3
conf_keratinocytes["EX_pydxn_LPAREN_e_RPAREN_"]=3
conf_keratinocytes["EX_ribflv_LPAREN_e_RPAREN_"]=3
```


```python
conf_HeLa["biomass_reaction"]=3
conf_HeLa["biomass_DNA"]=3
conf_HeLa["biomass_RNA"]=3
conf_HeLa["biomass_carbohydrate"]=3
conf_HeLa["biomass_lipid"]=3
conf_HeLa["biomass_other"]=3
conf_HeLa["biomass_protein"]=3
conf_HeLa["DM_atp_c_"]=3


conf_keratinocytes["biomass_reaction"]=3
conf_keratinocytes["biomass_DNA"]=3
conf_keratinocytes["biomass_RNA"]=3
conf_keratinocytes["biomass_carbohydrate"]=3
conf_keratinocytes["biomass_lipid"]=3
conf_keratinocytes["biomass_other"]=3
conf_keratinocytes["biomass_protein"]=3
conf_keratinocytes["DM_atp_c_"]=3
```


```python
metas = ['0.014 biomass_DNA_c + 0.058 biomass_RNA_c + 0.071 biomass_carbohydrate_c + 0.097 biomass_lipid_c + 0.054 biomass_other_c + 0.706 biomass_protein_c --> ', '3pg_c', '4abut_c', '4hpro_LT_c', 'accoa_m', 'accoa_m --> coa_m', 'ade_c', 'adn_c', 'adp_c', 'akg_m', 'ala_B_c', 'ala_L_c', 'amet_c', 'amp_c', 'arg_L_c', 'asn_L_c', 'asp_D_c', 'asp_L_c', 'atp_c', 'bhb_c', 'cdp_c', 'CE1936_c', 'chol_c', 'chsterol_c', 'cit_c', 'citr_L_c', 'cmp_c', 'creat_c', 'crm_hs_c', 'crtn_c', 'ctp_c', 'cys_L_c', 'dag_hs_c', 'dhap_c', 'e4p_c', 'f6p_c', 'fdp_c', 'fru_c', 'fum_c', 'g1p_c', 'g3p_c', 'g6p_c', 'gdp_c', 'glc_D_c', 'glc_D_e', 'glc_D_e --> glc_D_c', 'gln_L_m', 'gln_L_c', 'gln_L_e --> gln_L_c', 'glu_L_c', 'glyb_c', 'gly_c', 'gmp_c', 'gthox_c', 'gthrd_c', 'gua_c', 'HC00342_c', 'his_L_c', 'hxan_c', 'icit_c', 'ile_L_c', 'lac_L_c', 'leu_L_c', 'leu_L_c',  'lys_L_c', 'mag_hs_c', 'mal_L_c', 'met_L_c', 'nad_c', 'nadh_c', 'nadh_m', 'nad_m', 'nadp_c', 'nadph_c', 'nadph_m', 'nadph_m', 'nadp_m', 'oaa_m', 'orn_c', 'pa_hs_c', 'pe_hs_c', 'pep_c', 'phe_L_c', 'pmtcoa_c --> coa_c', 'pro_D_c', 'pro_L_c', 'ps_hs_c', 'ptrc_c', 'pyr_c', 'pyr_m', 'r5p_c', 'ru5p_D_c', 's7p_c', 'ser_L_c', 'spmd_c', 'succ_c', 'succoa_m --> coa_m', 'tag_hs_c', 'thr_L_c', 'trp_L_c', 'tym_c', 'tyr_L_c', 'udp_c', 'ump_c', 'utp_c', 'val_L_c']
```


```python
%%time


from corda import CORDA

opt_HeLa = CORDA(model=Recon2, confidence=conf_HeLa, n=5,met_prod=metas,  penalty_factor=1000) 
opt_HeLa.build()
print(opt_HeLa)
model_HeLa=opt_HeLa.cobra_model(name="HeLaModel")
print(model_HeLa.optimize())

```


```python
%%time

from corda import CORDA

opt_Keratinocytes = CORDA(model=Recon2, confidence=conf_keratinocytes, n=5, met_prod=metas,  penalty_factor=1000) 
opt_Keratinocytes.build()
print(opt_Keratinocytes)

model_Keratinocytes=opt_Keratinocytes.cobra_model(name="Keratinocytes")
print(model_Keratinocytes.optimize())
```


```python
cp=model_CancerBiopsy.copy()
cp.optimize()
cp.summary(fva=True)
```


```python
cp=model_Keratinocytes.copy()
cp.optimize()
cp.summary(fva=True)
```


```python
import pandas as pd

def get_drugable_targets(normal_Model, disease_Model, model_name,  eps=0.01):
    Nids = [r.id for r in normal_Model.reactions]
    Dids = [r.id for r in disease_Model.reactions]
        
    nmodel=normal_Model.copy()
    dmodel=disease_Model.copy()

    common_rxs = list(set(Nids) & set(Dids))
    print("Common reactions size",len(common_rxs))
    unique_Nrx = list(set(Nids) - set(Dids))
    print("Normal unique reactions size",len(unique_Nrx))
    unique_Drx = list(set(Dids) - set(Nids))
    print("Disease unique reactions size",len(unique_Drx))

    nflx0=normal_Model.optimize().f
    dflx0=disease_Model.optimize().f   

        
    results={}

    for rx in common_rxs:
        #print(rx)

        nbounds=nmodel.reactions.get_by_id(rx).bounds
        dbounds=dmodel.reactions.get_by_id(rx).bounds
        
        nmodel.reactions.get_by_id(rx).bounds=(-eps,eps)
        dmodel.reactions.get_by_id(rx).bounds=(-eps,eps)
                
        nfba=nmodel.optimize()    
        dfba=dmodel.optimize()
        
        nflx1=nfba.f
        dflx1=dfba.f
        
        results[rx]={}
        
        results[rx]["model"]=model_name
        results[rx]["gene_rule"]=nmodel.reactions.get_by_id(rx).gene_reaction_rule
         
        results[rx]["norm_flux"]=nflx0
        results[rx]["dise_flux"]=dflx0 
        
        results[rx]["del_norm_flux"]=nflx1
        results[rx]["del_dise_flux"]=dflx1 
            
       # results[rx]["norm_prolif_ratio"]=nflx1/nflx0
        #results[rx]["dise_prolif_ratio"]=dflx1/dflx0
        
        #results[rx]["norm_dise_ratio"]=(nflx1/nflx0)/(dflx1/dflx0)
        
        nmodel.reactions.get_by_id(rx).bounds=nbounds
        dmodel.reactions.get_by_id(rx).bounds=dbounds

        
    for rx in unique_Nrx:
        #print(rx)
        
        nbounds=nmodel.reactions.get_by_id(rx).bounds
        
        nmodel.reactions.get_by_id(rx).bounds=(-eps,eps)      
        
        nfba=nmodel.optimize()    
        
        nflx1=nfba.f
        
        results[rx]={}
        
        results[rx]["model"]=model_name
        results[rx]["gene_rule"]=nmodel.reactions.get_by_id(rx).gene_reaction_rule

        results[rx]["norm_flux"]=nflx0
        results[rx]["dise_flux"]=dflx0 
        
        results[rx]["del_norm_flux"]=nflx1
        results[rx]["del_dise_flux"]=dflx0

            
       # results[rx]["norm_prolif_ratio"]=nflx1/nflx0
       # results[rx]["dise_prolif_ratio"]=dflx0
        
       # results[rx]["norm_dise_ratio"]=nflx1/nflx0
            
        nmodel.reactions.get_by_id(rx).bounds=nbounds
        
    for rx in unique_Drx:
        #print(rx)
        
        dbounds=dmodel.reactions.get_by_id(rx).bounds
    
        dmodel.reactions.get_by_id(rx).bounds=(-eps,eps)
        dfba=dmodel.optimize()
        dflx1=dfba.f
        
        results[rx]={}

        results[rx]["model"]=model_name
        results[rx]["gene_rule"]=dmodel.reactions.get_by_id(rx).gene_reaction_rule

        results[rx]["norm_flux"]=nflx0
        results[rx]["dise_flux"]=dflx0 
        
        results[rx]["del_norm_flux"]=nflx0
        results[rx]["del_dise_flux"]=dflx1 

        #results[rx]["norm_prolif_ratio"]=nflx0
       # results[rx]["dise_prolif_ratio"]=dflx1/dflx0

        
       # results[rx]["norm_dise_ratio"]=1/(dflx1/dflx0)
        
        dmodel.reactions.get_by_id(rx).bounds=dbounds


    return(pd.DataFrame(results).transpose())

```


```python
R_CLines=get_drugable_targets(
    normal_Model=model_KeratinocytesMean, 
    disease_Model=model_HeLaMean, 
    model_name="CellLines")

R_CLines["cancer_rate"]=R_CLines["del_dise_flux"]/R_CLines["dise_flux"]
R_CLines["normal_rate"]=R_CLines["del_norm_flux"]/R_CLines["norm_flux"]
R_CLines["Total_rate"]=R_CLines["cancer_rate"]/R_CLines["normal_rate"]
R_CLines[(R_CLines["cancer_rate"]<0.999) & (R_CLines["normal_rate"]>.999) ].sort_values(by="Total_rate", ascending=True)


```


```python
%%time

from multiprocessing import Pool
from corda import CORDA

pool = Pool()

res1=pool.apply_async(model_create, [Recon2, conf_HeLaMean, 5, metas ,100] )
res2=pool.apply_async(model_create, [Recon2, conf_HeLaMean, 5, metas ,500] )
res3=pool.apply_async(model_create, [Recon2, conf_HeLaMean, 5, metas ,1000] )
res4=pool.apply_async(model_create, [Recon2, conf_HeLaMean, 5, metas ,2500] )
res5=pool.apply_async(model_create, [Recon2, conf_HeLaMean, 5, metas ,5000] )

res6=pool.apply_async(model_create, [Recon2, conf_keratinocytesMean, 5, metas ,100] )
res7=pool.apply_async(model_create, [Recon2, conf_keratinocytesMean, 5, metas ,500] )
res8=pool.apply_async(model_create, [Recon2, conf_keratinocytesMean, 5, metas ,1000] )
res9=pool.apply_async(model_create, [Recon2, conf_keratinocytesMean, 5, metas ,2500] )
res0=pool.apply_async(model_create, [Recon2, conf_keratinocytesMean, 5, metas ,5000] )

opt_HeLa1=res1.get() 
opt_HeLa2=res2.get() 
opt_HeLa3=res3.get() 
opt_HeLa4=res4.get() 
opt_HeLa5=res5.get() 

opt_Keratinocytes1=res6.get() 
opt_Keratinocytes2=res7.get() 
opt_Keratinocytes3=res8.get() 
opt_Keratinocytes4=res9.get() 
opt_Keratinocytes5=res0.get() 
```


```python
FBA=cobra.flux_analysis.pfba(model_HeLa)
FBA.fluxes.to_csv("HeLa_flux.csv")

FBA=cobra.flux_analysis.pfba(model_Keratinocytes)
FBA.fluxes.to_csv("Keratinocytes_flux.csv")
```


```python
cobra.io.write_sbml_model(model_HeLa, "Full_HeLa_model_0519_DEMEM6429_n5_FINAL.sbml")
cobra.io.write_sbml_model(model_Keratinocytes, "Full_Kerat_model_0519_DEMEM6429_n5_FINAL.sbml")
```


```python
closedModel=model_HeLaMean.copy()
closedModel.objective="DM_atp_c_"

for rx in closedModel.exchanges:
    rx.upper_bound= 100
    rx.lower_bound= 0
    

#Biomass
closedModel.reactions.biomass_other.lower_bound=0
closedModel.reactions.biomass_carbohydrate.lower_bound=0
closedModel.reactions.biomass_DNA.lower_bound=0
closedModel.reactions.biomass_lipid.lower_bound=0
closedModel.reactions.biomass_protein.lower_bound=0
closedModel.reactions.biomass_RNA.lower_bound=0
closedModel.reactions.biomass_reaction.lower_bound=0
closedModel.reactions.biomass_reaction.upper_bound=0



######################################################################
## Glucose aerobic
######################################################################

closedModel.reactions.EX_o2_LPAREN_e_RPAREN_.bounds=(-1000,-1)
closedModel.reactions.EX_h2o_LPAREN_e_RPAREN_.bounds=(-1000,1000)
closedModel.reactions.EX_co2_LPAREN_e_RPAREN_.bounds=(0,1000)
closedModel.reactions.EX_glc_LPAREN_e_RPAREN_.bounds=(-1,-1)
closedModel.reactions.EX_gln_L_LPAREN_e_RPAREN_.bounds=(0,0)

FBA = cobra.flux_analysis.pfba(closedModel)
print("===========================")
print("Glucose aerobic")
print("===========================")
print("Oxigen use",closedModel.reactions.EX_o2_LPAREN_e_RPAREN_.reaction,FBA["EX_o2_LPAREN_e_RPAREN_"])
print("H2O",closedModel.reactions.EX_h2o_LPAREN_e_RPAREN_.reaction,FBA["EX_h2o_LPAREN_e_RPAREN_"])
print("CO2 production",closedModel.reactions.EX_co2_LPAREN_e_RPAREN_.reaction,FBA["EX_co2_LPAREN_e_RPAREN_"])
print("Glucose consumption",closedModel.reactions.EX_glc_LPAREN_e_RPAREN_.reaction,FBA["EX_glc_LPAREN_e_RPAREN_"])
print("Glutamine consumption",closedModel.reactions.EX_gln_L_LPAREN_e_RPAREN_.reaction,FBA["EX_gln_L_LPAREN_e_RPAREN_"])
print("ATP production",FBA["DM_atp_c_"])

closedModel.summary(fva=True)
```


```python


######################################################################
## Glucose anaerobic
######################################################################

closedModel.reactions.EX_o2_LPAREN_e_RPAREN_.bounds=(0,0)
closedModel.reactions.EX_h2o_LPAREN_e_RPAREN_.bounds=(-1000,1000)
closedModel.reactions.EX_co2_LPAREN_e_RPAREN_.bounds=(0,1000)
closedModel.reactions.EX_glc_LPAREN_e_RPAREN_.bounds=(-1,-1)
closedModel.reactions.EX_gln_L_LPAREN_e_RPAREN_.bounds=(0,0)

FBA = cobra.flux_analysis.pfba(closedModel)
print("===========================")
print("Glucose anaerobic")
print("===========================")
print("Oxigen use",closedModel.reactions.EX_o2_LPAREN_e_RPAREN_.reaction,FBA["EX_o2_LPAREN_e_RPAREN_"])
print("H2O",closedModel.reactions.EX_h2o_LPAREN_e_RPAREN_.reaction,FBA["EX_h2o_LPAREN_e_RPAREN_"])
print("CO2 production",closedModel.reactions.EX_co2_LPAREN_e_RPAREN_.reaction,FBA["EX_co2_LPAREN_e_RPAREN_"])
print("Glucose consumption",closedModel.reactions.EX_glc_LPAREN_e_RPAREN_.reaction,FBA["EX_glc_LPAREN_e_RPAREN_"])
print("Glutamine consumption",closedModel.reactions.EX_gln_L_LPAREN_e_RPAREN_.reaction,FBA["EX_gln_L_LPAREN_e_RPAREN_"])
print("ATP production",FBA["DM_atp_c_"])

closedModel.summary(fva=True)
```


```python

```


```python
model_HeLaMean1=opt_HeLa1.cobra_model(name="HeLaMean2")
print(model_HeLaMean1.optimize())
model_HeLaMean2=opt_HeLa2.cobra_model(name="HeLaMean2")
print(model_HeLaMean2.optimize())
model_HeLaMean3=opt_HeLa3.cobra_model(name="HeLaMean3")
print(model_HeLaMean3.optimize())
model_HeLaMean4=opt_HeLa4.cobra_model(name="HeLaMean4")
print(model_HeLaMean4.optimize())
model_HeLaMean5=opt_HeLa5.cobra_model(name="HeLaMean5")
print(model_HeLaMean5.optimize())
```


```python
model_KeratinocytesMean1=opt_Keratinocytes1.cobra_model(name="Keratinocytes1")
print(model_KeratinocytesMean1.optimize())
model_KeratinocytesMean2=opt_Keratinocytes2.cobra_model(name="Keratinocytes2")
print(model_KeratinocytesMean2.optimize())
model_KeratinocytesMean3=opt_Keratinocytes3.cobra_model(name="Keratinocytes3")
print(model_KeratinocytesMean3.optimize())
model_KeratinocytesMean4=opt_Keratinocytes4.cobra_model(name="Keratinocytes4")
print(model_KeratinocytesMean4.optimize())
model_KeratinocytesMean5=opt_Keratinocytes5.cobra_model(name="Keratinocytes5")
print(model_KeratinocytesMean5.optimize())
```


```python
R_CLines=get_drugable_targets(
    normal_Model=model_KeratinocytesMean1, 
    disease_Model=model_HeLaMean1, 
    model_name="CellLines")
R_CLines["Dise_ratio"]=R_CLines["del_dise_flux"]/R_CLines["dise_flux"]
R_CLines["Norm_ratio"]=R_CLines["del_norm_flux"]/R_CLines["norm_flux"]
R_CLines[R_CLines["Dise_ratio"]< 0.999].sort_values(by="Dise_ratio", ascending=True)
```


```python
R_CLines=get_drugable_targets(
    normal_Model=model_KeratinocytesMean2, 
    disease_Model=model_HeLaMean2, 
    model_name="CellLines")
R_CLines["Dise_ratio"]=R_CLines["del_dise_flux"]/R_CLines["dise_flux"]
R_CLines["Norm_ratio"]=R_CLines["del_norm_flux"]/R_CLines["norm_flux"]
R_CLines[R_CLines["Dise_ratio"]< 0.999].sort_values(by="Dise_ratio", ascending=True)
```


```python
R_CLines=get_drugable_targets(
    normal_Model=model_KeratinocytesMean3, 
    disease_Model=model_HeLaMean3, 
    model_name="CellLines")
R_CLines["Dise_ratio"]=R_CLines["del_dise_flux"]/R_CLines["dise_flux"]
R_CLines["Norm_ratio"]=R_CLines["del_norm_flux"]/R_CLines["norm_flux"]
R_CLines[R_CLines["Dise_ratio"]< 0.999].sort_values(by="Dise_ratio", ascending=True)
```


```python
R_CLines=get_drugable_targets(
    normal_Model=model_KeratinocytesMean4, 
    disease_Model=model_HeLaMean4, 
    model_name="CellLines")
R_CLines["Dise_ratio"]=R_CLines["del_dise_flux"]/R_CLines["dise_flux"]
R_CLines["Norm_ratio"]=R_CLines["del_norm_flux"]/R_CLines["norm_flux"]
R_CLines[R_CLines["Dise_ratio"]< 0.999].sort_values(by="Dise_ratio", ascending=True)
```


```python
R_CLines=get_drugable_targets(
    normal_Model=model_KeratinocytesMean5, 
    disease_Model=model_HeLaMean5, 
    model_name="CellLines")
R_CLines["Dise_ratio"]=R_CLines["del_dise_flux"]/R_CLines["dise_flux"]
R_CLines["Norm_ratio"]=R_CLines["del_norm_flux"]/R_CLines["norm_flux"]
R_CLines[R_CLines["Dise_ratio"]< 0.999].sort_values(by="Dise_ratio", ascending=True)
```


```python
model_KeratinocytesMean1.summary(fva=True)
```


```python

```
