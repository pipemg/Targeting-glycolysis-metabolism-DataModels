

```python
import cobra
#!wget http://www.ebi.ac.uk/biomodels-main/download?mid=MODEL1603150001 -O recon2.2.xml
Recon2 = cobra.io.read_sbml_model("Models/recon2.2.xml")
#Recon2 = cobra.io.load_matlab_model("Models/Recon2.v04.mat")
```

    cobra/io/sbml.py:235 [1;31mUserWarning[0m: M_h_x appears as a reactant and product FAOXC220200x



```python
print("Reactions:", len(Recon2.reactions))
print("Metabolites:", len(Recon2.metabolites))
print(Recon2.optimize())
```

    Reactions: 7785
    Metabolites: 5324
    <Solution 555.786 at 0x7f8d8c337f60>



```python
for rxex in Recon2.exchanges:
    rxex.bounds=(-0.01,10)

### Western Diet
#### Sugars
Recon2.reactions.get_by_id("EX_arab_L_LPAREN_e_RPAREN_").lower_bound=-0.17878295 
#Cellobiose not in Humana
#Recon2.reactions.get_by_id("EX_drib_LPAREN_e_RPAREN_").lower_bound=-0.17878295 #????
Recon2.reactions.get_by_id("EX_fru_LPAREN_e_RPAREN_").lower_bound=-0.14898579
Recon2.reactions.get_by_id("EX_fuc_L_LPAREN_e_RPAREN_").lower_bound=-0.14898579 
Recon2.reactions.get_by_id("EX_gal_LPAREN_e_RPAREN_").lower_bound=-0.14898579 
Recon2.reactions.get_by_id("EX_glc_LPAREN_e_RPAREN_").lower_bound=-0.14898579
#Falta Ex_glcn en Recon2.2
Recon2.reactions.get_by_id("EX_lcts_LPAREN_e_RPAREN_").lower_bound=-0.07449289 
Recon2.reactions.get_by_id("EX_malt_LPAREN_e_RPAREN_").lower_bound=-0.07449289 
Recon2.reactions.get_by_id("EX_man_LPAREN_e_RPAREN_").lower_bound=-0.14898579 
#Melibiose  not in Human
#D-Mannitol not in Human
Recon2.reactions.get_by_id("EX_oxa_LPAREN_e_RPAREN_").lower_bound=-0.44695737 
Recon2.reactions.get_by_id("EX_rib_D_LPAREN_e_RPAREN_").lower_bound=-0.17878295 
#L-Rhamnose not in Human
Recon2.reactions.get_by_id("EX_sucr_LPAREN_e_RPAREN_").lower_bound=-0.07449289 

#Estos no estan en sangre
#Recon2.reactions.get_by_id("EX_tre_LPAREN_e_RPAREN_").lower_bound=-0.07449289 
#Recon2.reactions.get_by_id("EX_xyl_D_LPAREN_e_RPAREN_").lower_bound=-0.17878295
#Recon2.reactions.get_by_id("EX_strch1_LPAREN_e_RPAREN_").lower_bound=-0.25733909 

#fiber not in human

### Western Diet
#### Fat
Recon2.reactions.get_by_id("EX_arachd_LPAREN_e_RPAREN_").lower_bound=-0.00332813 
Recon2.reactions.get_by_id("EX_chsterol_LPAREN_e_RPAREN_").lower_bound=-0.00495795
Recon2.reactions.get_by_id("EX_glyc_LPAREN_e_RPAREN_").lower_bound=-1.79965486 
Recon2.reactions.get_by_id("EX_hdca_LPAREN_e_RPAREN_").lower_bound=-0.39637090
Recon2.reactions.get_by_id("EX_hdcea_LPAREN_e_RPAREN_").lower_bound=-0.03651697
Recon2.reactions.get_by_id("EX_lnlc_LPAREN_e_RPAREN_").lower_bound=-0.35910921
Recon2.reactions.get_by_id("EX_lnlnca_LPAREN_e_RPAREN_").lower_bound=-0.01756512
Recon2.reactions.get_by_id("EX_lnlncg_LPAREN_e_RPAREN_").lower_bound=-0.01756512
Recon2.reactions.get_by_id("EX_ocdca_LPAREN_e_RPAREN_").lower_bound=-0.16928260
Recon2.reactions.get_by_id("EX_ocdcea_LPAREN_e_RPAREN_").lower_bound=-0.68144465
#No en sangre
#Recon2.reactions.get_by_id("EX_octa_LPAREN_e_RPAREN_").lower_bound=-0.01294272
#Recon2.reactions.get_by_id("EX_ttdca_LPAREN_e_RPAREN_").lower_bound=-0.06867567

### Western Diet
#### Protein
Recon2.reactions.get_by_id("EX_ala_L_LPAREN_e_RPAREN_").lower_bound=-1
Recon2.reactions.get_by_id("EX_arg_L_LPAREN_e_RPAREN_").lower_bound=-0.15
Recon2.reactions.get_by_id("EX_asn_L_LPAREN_e_RPAREN_").lower_bound=-0.225
Recon2.reactions.get_by_id("EX_asp_L_LPAREN_e_RPAREN_").lower_bound=-0.225
Recon2.reactions.get_by_id("EX_cys_L_LPAREN_e_RPAREN_").lower_bound=-1
Recon2.reactions.get_by_id("EX_gln_L_LPAREN_e_RPAREN_").lower_bound=-0.18
Recon2.reactions.get_by_id("EX_glu_L_LPAREN_e_RPAREN_").lower_bound=-0.18
Recon2.reactions.get_by_id("EX_gly_LPAREN_e_RPAREN_").lower_bound=-0.45
Recon2.reactions.get_by_id("EX_his_L_LPAREN_e_RPAREN_").lower_bound=-0.15
Recon2.reactions.get_by_id("EX_ile_L_LPAREN_e_RPAREN_").lower_bound=-0.15
Recon2.reactions.get_by_id("EX_leu_L_LPAREN_e_RPAREN_").lower_bound=-0.15
Recon2.reactions.get_by_id("EX_lys_L_LPAREN_e_RPAREN_").lower_bound=-0.15 
Recon2.reactions.get_by_id("EX_met_L_LPAREN_e_RPAREN_").lower_bound=-0.18
Recon2.reactions.get_by_id("EX_phe_L_LPAREN_e_RPAREN_").lower_bound=-1
Recon2.reactions.get_by_id("EX_pro_L_LPAREN_e_RPAREN_").lower_bound=-0.18
Recon2.reactions.get_by_id("EX_ser_L_LPAREN_e_RPAREN_").lower_bound=-1
Recon2.reactions.get_by_id("EX_thr_L_LPAREN_e_RPAREN_").lower_bound=-0.225
Recon2.reactions.get_by_id("EX_trp_L_LPAREN_e_RPAREN_").lower_bound=-0.08181818
Recon2.reactions.get_by_id("EX_tyr_L_LPAREN_e_RPAREN_").lower_bound=-1
Recon2.reactions.get_by_id("EX_val_L_LPAREN_e_RPAREN_").lower_bound=-0.18

### Western Diet
#### Minterals, Vitamins, others

#12dgr180
#26dap_M
#2dmmq8
#2obut
Recon2.reactions.get_by_id("EX_3mop_LPAREN_e_RPAREN_").lower_bound=-1
#4abz
#4hbz
Recon2.reactions.get_by_id("EX_ac_LPAREN_e_RPAREN_").lower_bound=-1
Recon2.reactions.get_by_id("EX_acgam_LPAREN_e_RPAREN_").lower_bound=-1
#No en sangre acmana
#Recon2.reactions.get_by_id("EX_acmana_LPAREN_e_RPAREN_").lower_bound=-1
#acnam
Recon2.reactions.get_by_id("EX_ade_LPAREN_e_RPAREN_").lower_bound=-1
Recon2.reactions.get_by_id("EX_adn_LPAREN_e_RPAREN_").lower_bound=-1
#adocbl
Recon2.reactions.get_by_id("EX_ala_D_LPAREN_e_RPAREN_").lower_bound=-1
Recon2.reactions.get_by_id("EX_amp_LPAREN_e_RPAREN_").lower_bound=-1
#arab_D
#No en sangre btn
#Recon2.reactions.get_by_id("EX_btn_LPAREN_e_RPAREN_").lower_bound=-1
Recon2.reactions.get_by_id("EX_ca2_LPAREN_e_RPAREN_").lower_bound=-1
#cbl1
Recon2.reactions.get_by_id("EX_cgly_LPAREN_e_RPAREN_").lower_bound=-1
Recon2.reactions.get_by_id("EX_chol_LPAREN_e_RPAREN_").lower_bound=-1
Recon2.reactions.get_by_id("EX_cit_LPAREN_e_RPAREN_").lower_bound=-1
Recon2.reactions.get_by_id("EX_cl_LPAREN_e_RPAREN_").lower_bound=-1
#cobalt2
Recon2.reactions.get_by_id("EX_csn_LPAREN_e_RPAREN_").lower_bound=-1
#cu2
Recon2.reactions.get_by_id("EX_dad_2_LPAREN_e_RPAREN_").lower_bound=-1
Recon2.reactions.get_by_id("EX_dcyt_LPAREN_e_RPAREN_").lower_bound=-1
Recon2.reactions.get_by_id("EX_ddca_LPAREN_e_RPAREN_").lower_bound=-1
Recon2.reactions.get_by_id("EX_dgsn_LPAREN_e_RPAREN_").lower_bound=-1
Recon2.reactions.get_by_id("EX_fe2_LPAREN_e_RPAREN_").lower_bound=-1
# Methemoglobin normalmente no en sangre
#Recon2.reactions.get_by_id("EX_fe3_LPAREN_e_RPAREN_").lower_bound=-1
#fe3dcit
Recon2.reactions.get_by_id("EX_fald_LPAREN_e_RPAREN_").lower_bound=-1
Recon2.reactions.get_by_id("EX_fol_LPAREN_e_RPAREN_").lower_bound=-1
Recon2.reactions.get_by_id("EX_for_LPAREN_e_RPAREN_").lower_bound=-1
Recon2.reactions.get_by_id("EX_fum_LPAREN_e_RPAREN_").lower_bound=-1
Recon2.reactions.get_by_id("EX_gam_LPAREN_e_RPAREN_").lower_bound=-1
Recon2.reactions.get_by_id("EX_glu_L_LPAREN_e_RPAREN_").lower_bound=-1
Recon2.reactions.get_by_id("EX_glyc3p_LPAREN_e_RPAREN_").lower_bound=-1
Recon2.reactions.get_by_id("EX_gthox_LPAREN_e_RPAREN_").lower_bound=-1
Recon2.reactions.get_by_id("EX_gthrd_LPAREN_e_RPAREN_").lower_bound=-1
Recon2.reactions.get_by_id("EX_gua_LPAREN_e_RPAREN_").lower_bound=-1
Recon2.reactions.get_by_id("EX_h_LPAREN_e_RPAREN_").lower_bound=-1
Recon2.reactions.get_by_id("EX_h2o_LPAREN_e_RPAREN_").lower_bound=-10
Recon2.reactions.get_by_id("EX_h2o2_LPAREN_e_RPAREN_").lower_bound=-0.01

#h2s
Recon2.reactions.get_by_id("EX_hxan_LPAREN_e_RPAREN_").lower_bound=-1
#indole
Recon2.reactions.get_by_id("EX_k_LPAREN_e_RPAREN_").lower_bound=-1
#lanost
Recon2.reactions.get_by_id("EX_meoh_LPAREN_e_RPAREN_").lower_bound=-10
#metsox_S_L
#mg2
#mn2
#mobd
#mqn8
Recon2.reactions.get_by_id("EX_na1_LPAREN_e_RPAREN_").lower_bound=-1
Recon2.reactions.get_by_id("EX_nac_LPAREN_e_RPAREN_").lower_bound=-1
Recon2.reactions.get_by_id("EX_ncam_LPAREN_e_RPAREN_").lower_bound=-1
#nmn nmn_e doesn't exists in Recon 2.2 
Recon2.reactions.get_by_id("EX_no2_LPAREN_e_RPAREN_").lower_bound=-1
Recon2.reactions.get_by_id("EX_orn_LPAREN_e_RPAREN_").lower_bound=-1
# No en sangre
#Recon2.reactions.get_by_id("EX_pheme_LPAREN_e_RPAREN_").lower_bound=-1

Recon2.reactions.get_by_id("EX_pi_LPAREN_e_RPAREN_").lower_bound=-1
#Pimelate  not in Human
Recon2.reactions.get_by_id("EX_pnto_R_LPAREN_e_RPAREN_").lower_bound=-1
Recon2.reactions.get_by_id("EX_ptrc_LPAREN_e_RPAREN_").lower_bound=-1
# No en sangre
#Recon2.reactions.get_by_id("EX_pydam_LPAREN_e_RPAREN_").lower_bound=-1
Recon2.reactions.get_by_id("EX_pydx_LPAREN_e_RPAREN_").lower_bound=-1
Recon2.reactions.get_by_id("EX_pydx5p_LPAREN_e_RPAREN_").lower_bound=-1
Recon2.reactions.get_by_id("EX_pydxn_LPAREN_e_RPAREN_").lower_bound=-1
#ubiquinone-8 not in Human
Recon2.reactions.get_by_id("EX_ribflv_LPAREN_e_RPAREN_").lower_bound=-1
#No en sangre
Recon2.reactions.get_by_id("EX_sel_LPAREN_e_RPAREN_").lower_bound=-1
#Siroheme not in Human
Recon2.reactions.get_by_id("EX_so4_LPAREN_e_RPAREN_").lower_bound=-1
Recon2.reactions.get_by_id("EX_spmd_LPAREN_e_RPAREN_").lower_bound=-1
Recon2.reactions.get_by_id("EX_thm_LPAREN_e_RPAREN_").lower_bound=-1
Recon2.reactions.get_by_id("EX_thymd_LPAREN_e_RPAREN_").lower_bound=-1
Recon2.reactions.get_by_id("EX_ura_LPAREN_e_RPAREN_").lower_bound=-1
Recon2.reactions.get_by_id("EX_uri_LPAREN_e_RPAREN_").lower_bound=-1
#Xanthine not in Human
#zn2 not in Human

Recon2.reactions.biomass_reaction.lower_bound=0.008
Recon2.reactions.biomass_reaction.upper_bound=0.1


```


```python
Recon2.optimize()
```




    <Solution 0.100 at 0x7f8d8c1a1a20>




```python
import pandas as pd
ub_scores_matrix =pd.read_csv("Binary_Byposys_RECON2_2.tsv",sep=",",index_col=0,header=0)
confidence_scores_matrix=ub_scores_matrix.copy()
confidence_scores_matrix.head()
```




<div>
<style>
    .dataframe thead tr:only-child th {
        text-align: right;
    }

    .dataframe thead th {
        text-align: left;
    }

    .dataframe tbody tr th {
        vertical-align: top;
    }
</style>
<table border="1" class="dataframe">
  <thead>
    <tr style="text-align: right;">
      <th></th>
      <th>CancerBiopsyHGU133A</th>
      <th>NormalBiopsyHGU133A</th>
      <th>CancerBiopsyHGU133Plus2</th>
      <th>NormalBiopsyHGU133APlus2</th>
      <th>MaxCancerBiopsy</th>
      <th>MaxNormalBiopsy</th>
      <th>MeanCancerBiopsy</th>
      <th>MeanNormalBiopsy</th>
    </tr>
  </thead>
  <tbody>
    <tr>
      <th>HGNC:10006</th>
      <td>0.000000</td>
      <td>0.0000</td>
      <td>0.000000</td>
      <td>0.0</td>
      <td>0.000000</td>
      <td>0.0000</td>
      <td>0.000000</td>
      <td>0.00000</td>
    </tr>
    <tr>
      <th>HGNC:1027</th>
      <td>0.015625</td>
      <td>0.0625</td>
      <td>0.116667</td>
      <td>0.0</td>
      <td>0.116667</td>
      <td>0.0625</td>
      <td>0.066146</td>
      <td>0.03125</td>
    </tr>
    <tr>
      <th>HGNC:10293</th>
      <td>0.000000</td>
      <td>0.0000</td>
      <td>0.483333</td>
      <td>0.0</td>
      <td>0.483333</td>
      <td>0.0000</td>
      <td>0.241667</td>
      <td>0.00000</td>
    </tr>
    <tr>
      <th>HGNC:10297</th>
      <td>0.000000</td>
      <td>0.0000</td>
      <td>0.541667</td>
      <td>0.0</td>
      <td>0.541667</td>
      <td>0.0000</td>
      <td>0.270833</td>
      <td>0.00000</td>
    </tr>
    <tr>
      <th>HGNC:10451</th>
      <td>0.878472</td>
      <td>0.6000</td>
      <td>1.000000</td>
      <td>0.5</td>
      <td>1.000000</td>
      <td>0.6000</td>
      <td>0.939236</td>
      <td>0.55000</td>
    </tr>
  </tbody>
</table>
</div>




```python
for model in ub_scores_matrix.columns:
    #High Confidence rate 
       
    HC=ub_scores_matrix[model][ub_scores_matrix[model] >= .70].index.tolist()
    confidence_scores_matrix[model][HC]=3
        
    #Medium Confidence Rate 
    MC=ub_scores_matrix[model][(ub_scores_matrix[model] <  .70) & (ub_scores_matrix[model] >= 0.50) ].index.tolist()
    confidence_scores_matrix[model][MC]=2
        
    #Low
    LC=ub_scores_matrix[model][(ub_scores_matrix[model] < 0.40) & (ub_scores_matrix[model] >=0.3) ].index.tolist()
    confidence_scores_matrix[model][LC]=1
    
    #unknown
    ZC=ub_scores_matrix[model][(ub_scores_matrix[model] <.3) & (ub_scores_matrix[model] >= 0.1) ].index.tolist()
    confidence_scores_matrix[model][ZC]=0
    
    #Negative
    NC=ub_scores_matrix[model][(ub_scores_matrix[model] < 0.1)].index.tolist()
    confidence_scores_matrix[model][NC]=-1
 
```


```python
from corda import reaction_confidence


conf_MeanCancerBiopsy = {}
conf_MeanNormalBiopsy = {}

for r in Recon2.reactions:
    if(r.gene_reaction_rule!=''):
        conf_MeanCancerBiopsy[r.id]=reaction_confidence(r.gene_reaction_rule,confidence_scores_matrix["MeanCancerBiopsy"])
        conf_MeanNormalBiopsy[r.id]=reaction_confidence(r.gene_reaction_rule,confidence_scores_matrix["MeanNormalBiopsy"])
    else:
        conf_MeanCancerBiopsy[r.id]=0
        conf_MeanNormalBiopsy[r.id]=0
```


```python
conf_MeanCancerBiopsy["biomass_reaction"]=3
conf_MeanCancerBiopsy["biomass_DNA"]=3
conf_MeanCancerBiopsy["biomass_RNA"]=3
conf_MeanCancerBiopsy["biomass_carbohydrate"]=3
conf_MeanCancerBiopsy["biomass_lipid"]=3
conf_MeanCancerBiopsy["biomass_other"]=3
conf_MeanCancerBiopsy["biomass_protein"]=3
conf_MeanCancerBiopsy["DM_atp_c_"]=3

conf_MeanNormalBiopsy["biomass_reaction"]=3
conf_MeanNormalBiopsy["biomass_DNA"]=3
conf_MeanNormalBiopsy["biomass_RNA"]=3
conf_MeanNormalBiopsy["biomass_carbohydrate"]=3
conf_MeanNormalBiopsy["biomass_lipid"]=3
conf_MeanNormalBiopsy["biomass_other"]=3
conf_MeanNormalBiopsy["biomass_protein"]=3
conf_MeanNormalBiopsy["DM_atp_c_"]=3

Recon2.objective="biomass_reaction"

```


```python
metas = ['adp_c', 'atp_c','glc_D_c',  'fru_c', 'nad_c', 'nadh_c','nad_m', 'nadh_m','nadph_c', 'nadph_m', 'nadp_c', 'nadp_m', 'cmp_c', 'HC00342_c', 'glcn_c', 'citr_L_c', 'glyb_c', 'icit_c', '3pg_c', 'accoa_m ->coa_m', 'akg_m', 'e4p_c', 'f6p_c', 'g3p_c', 'g6p_c', 'oaa_m', 'pep_c', 'pyr_c', 'r5p_c', 'succoa_m ->coa_m', 'ala_L_c', 'arg_L_c', 'asn_L_c', 'asp_L_c', 'val_L_c', 'adp_c', 'thr_L_c', 'leu_L_c', 'gln_L_c', 'glu_L_c', 'gly_c', 'pro_D_c', 'pro_L_c', 'ser_L_c', 'ctp_c', 'fdp_c','utp_c', 'pmtcoa_c -> coa_c', 'chsterol_c', 'dag_hs_c', 'tag_hs_c', 'mag_hs_c', 'gthox_c','gthrd_c','ru5p_D_c','crm_hs_c', 'pa_hs_c', 'pe_hs_c', 'ps_hs_c', 'hxan_c', 'creat_c', 'crtn_c', 'chol_c', 'orn_c', '4hpro_LT_c', 's7p_c', 'amp_c', 'udp_c', 'ade_c', 'asp_D_c', 'adn_c', 'ump_c', 'his_L_c', 'tyr_L_c', 'phe_L_c', 'gdp_c', 'ctp_c', 'cys_L_c', 'amet_c', 'cit_c', 'tym_c', 'succ_c', 'CE1936_c', 'gua_c', 'cdp_c', 'g1p_c', 'lys_L_c', 'bhb_c', 'ile_L_c', 'gmp_c', 'dhap_c', 'fum_c', 'mal_L_c', 'spmd_c', 'ala_B_c', 'trp_L_c', 'lac_L_c', 'met_L_c', 'ptrc_c', '4abut_c','0.014 biomass_DNA_c + 0.058 biomass_RNA_c + 0.071 biomass_carbohydrate_c + 0.097 biomass_lipid_c + 0.054 biomass_other_c + 0.706 biomass_protein_c --> ']
```


```python
%%time

from corda import CORDA

opt_MeanCancerBiopsy = CORDA(model=Recon2, confidence=conf_MeanCancerBiopsy, n=10, met_prod=metas,  penalty_factor=100 ) 
opt_MeanCancerBiopsy.build()
print(opt_MeanCancerBiopsy)
model_MeanCancerBiopsy=opt_MeanCancerBiopsy.cobra_model(name="MeanCancerBiopsy")
print(model_MeanCancerBiopsy.optimize())
```


    ---------------------------------------------------------------------------

    ValueError                                Traceback (most recent call last)

    <timed exec> in <module>()


    /opt/conda/lib/python3.5/site-packages/corda/corda.py in __init__(self, model, confidence, met_prod, n, penalty_factor, support)
        114             if r.id in confidence:
        115                 if confidence[r.id] not in [-1, 0, 1, 2, 3]:
    --> 116                     raise ValueError("Not a valid confidence value!")
        117                 else:
        118                     self.conf[r.id] = confidence[r.id]


    ValueError: Not a valid confidence value!



```python
model_MeanCancerBiopsy.summary(fva=True)
```


    ---------------------------------------------------------------------------

    NameError                                 Traceback (most recent call last)

    <ipython-input-15-df9b175e60c0> in <module>()
    ----> 1 model_MeanCancerBiopsy.summary(fva=True)
    

    NameError: name 'model_MeanCancerBiopsy' is not defined



```python
%%time

from corda import CORDA


opt_MeanNormalBiopsy = CORDA(model=Recon2, confidence=conf_MeanNormalBiopsy, n=10, met_prod=metas,  penalty_factor=100) 
opt_MeanNormalBiopsy.build()
print(opt_MeanNormalBiopsy)
model_MeanNormalBiopsy=opt_MeanNormalBiopsy.cobra_model(name="MeanNormalBiopsy")
print(model_MeanNormalBiopsy.optimize())
```


```python
model_MeanNormalBiopsy.summary(fva=True)
```


```python
cobra.io.write_sbml_model(model_MeanCancerBiopsy, "Full_SCC_model_1218_western_diet_n10_pf100.sbml")
cobra.io.write_sbml_model(model_MeanNormalBiopsy, "Full_Normal_model_1218_western_diet_n10_pf100.sbml")
```


```python
import pandas as pd

def get_drugable_targets(normal_Model, disease_Model, model_name,  eps=0.1):
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
        
        results[rx]["del_norm_flux"]=nflx1
        results[rx]["del_dise_flux"]=dflx1 

        #results[rx]["norm_prolif_ratio"]=nflx0
       # results[rx]["dise_prolif_ratio"]=dflx1/dflx0

        
       # results[rx]["norm_dise_ratio"]=1/(dflx1/dflx0)
        
        dmodel.reactions.get_by_id(rx).bounds=dbounds


    return(pd.DataFrame(results).transpose())

```


```python
R_Biopsies=get_drugable_targets(
    normal_Model=model_MeanNormalBiopsy, 
    disease_Model=model_MeanCancerBiopsy, model_name="Biopsys", eps=0.1)
```


    ---------------------------------------------------------------------------

    NameError                                 Traceback (most recent call last)

    <ipython-input-3-d8a72369d17a> in <module>()
    ----> 1 R_Biopsies=get_drugable_targets(
          2     normal_Model=model_MeanNormalBiopsy,
          3     disease_Model=model_MeanCancerBiopsy, model_name="Biopsys", eps=0.1)


    NameError: name 'get_drugable_targets' is not defined



```python
R_Biopsies["dise_prolif_ratio"]=R_Biopsies["del_dise_flux"]/R_Biopsies["dise_flux"]
R_Biopsies["norm_prolif_ratio"]=R_Biopsies["del_norm_flux"]/R_Biopsies["norm_flux"]
R_Biopsies["total_ratio"]=R_Biopsies["dise_prolif_ratio"]/R_Biopsies["norm_prolif_ratio"]

```


    ---------------------------------------------------------------------------

    NameError                                 Traceback (most recent call last)

    <ipython-input-4-c48cb8ff9e37> in <module>()
    ----> 1 R_Biopsies["dise_prolif_ratio"]=R_Biopsies["del_dise_flux"]/R_Biopsies["dise_flux"]
          2 R_Biopsies["norm_prolif_ratio"]=R_Biopsies["del_norm_flux"]/R_Biopsies["norm_flux"]
          3 R_Biopsies["total_ratio"]=R_Biopsies["dise_prolif_ratio"]/R_Biopsies["norm_prolif_ratio"]


    NameError: name 'R_Biopsies' is not defined



```python
R_Biopsies[R_Biopsies["dise_prolif_ratio"]<.999].sort_values(by="total_ratio", ascending=True)

```




<div>
<style>
    .dataframe thead tr:only-child th {
        text-align: right;
    }

    .dataframe thead th {
        text-align: left;
    }

    .dataframe tbody tr th {
        vertical-align: top;
    }
</style>
<table border="1" class="dataframe">
  <thead>
    <tr style="text-align: right;">
      <th></th>
      <th>del_dise_flux</th>
      <th>del_norm_flux</th>
      <th>dise_flux</th>
      <th>gene_rule</th>
      <th>model</th>
      <th>norm_flux</th>
      <th>dise_prolif_ratio</th>
      <th>norm_prolif_ratio</th>
      <th>total_ratio</th>
    </tr>
  </thead>
  <tbody>
    <tr>
      <th>GAPD</th>
      <td>0.00914891</td>
      <td>0.0223156</td>
      <td>0.0562509</td>
      <td>HGNC:24864 or HGNC:4141</td>
      <td>Biopsys</td>
      <td>0.0223156</td>
      <td>0.162645</td>
      <td>1</td>
      <td>0.162645</td>
    </tr>
    <tr>
      <th>ENO</th>
      <td>0.00914891</td>
      <td>0.0223156</td>
      <td>0.0562509</td>
      <td>HGNC:3350 or HGNC:3354 or HGNC:3354 or HGNC:3353</td>
      <td>Biopsys</td>
      <td>0.0223156</td>
      <td>0.162645</td>
      <td>1</td>
      <td>0.162645</td>
    </tr>
    <tr>
      <th>PGM</th>
      <td>0.00914891</td>
      <td>0.0223156</td>
      <td>0.0562509</td>
      <td>HGNC:1093 or HGNC:8888 or HGNC:8889</td>
      <td>Biopsys</td>
      <td>0.0223156</td>
      <td>0.162645</td>
      <td>1</td>
      <td>0.162645</td>
    </tr>
    <tr>
      <th>PGK</th>
      <td>0.00914891</td>
      <td>0.0223156</td>
      <td>0.0562509</td>
      <td>HGNC:8896 or HGNC:8898</td>
      <td>Biopsys</td>
      <td>0.0223156</td>
      <td>0.162645</td>
      <td>1</td>
      <td>0.162645</td>
    </tr>
  </tbody>
</table>
</div>




```python
cp_model_MeanCancerBiopsy=model_MeanCancerBiopsy.copy()
cp_model_MeanCancerBiopsy.reactions.CMPSAS.bounds=(-0.1,0.1)
print(cp_model_MeanCancerBiopsy.optimize())
print(model_MeanCancerBiopsy.optimize())
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


    KeyError: 'CMPSAS'

    
    During handling of the above exception, another exception occurred:


    AttributeError                            Traceback (most recent call last)

    <ipython-input-100-8c2bb1473fe9> in <module>()
          1 cp_model_MeanCancerBiopsy=model_MeanCancerBiopsy.copy()
    ----> 2 cp_model_MeanCancerBiopsy.reactions.CMPSAS.bounds=(-0.1,0.1)
          3 print(cp_model_MeanCancerBiopsy.optimize())
          4 print(model_MeanCancerBiopsy.optimize())


    /opt/conda/lib/python3.5/site-packages/cobra/core/dictlist.py in __getattr__(self, attr)
        449         except KeyError:
        450             raise AttributeError("DictList has no attribute or entry %s" %
    --> 451                                  attr)
        452 
        453     def __dir__(self):


    AttributeError: DictList has no attribute or entry CMPSAS



```python
cp_model_MeanNormalBiopsy=model_MeanNormalBiopsy.copy()
cp_model_MeanNormalBiopsy.reactions.CMPSAS.bounds=(-0.1,0.1)
print(cp_model_MeanNormalBiopsy.optimize())
print(model_MeanNormalBiopsy.optimize())
```

    <Solution 0.100 at 0x7f3c046d66d8>
    <Solution 0.071 at 0x7f3c046d6668>



```python

```


```python
import random
import numpy


def shuffle(conf):
    lst = [x for x in range(0,len(conf.items()))]

    random.shuffle(lst)
    return lst

def randomize_model(conf_scores):
    lst=shuffle(conf_scores)
    a=numpy.array(lst)
    (rxindex,rxvalues) = zip(*conf_scores.items())
    rxvalues=numpy.array(rxvalues)

    nrxval=rxvalues[a]

    dictionary = dict(zip(list(rxindex), list(nrxval)))
    return(dictionary)



```


```python
for model in ub_scores_matrix.columns:
    #High Confidence rate 
       
    HC=ub_scores_matrix[model][ub_scores_matrix[model] >= 0.9].index.tolist()
    confidence_scores_matrix[model][HC]=3
        
    #Medium Confidence Rate 
    MC=ub_scores_matrix[model][(ub_scores_matrix[model] <  0.9) & (ub_scores_matrix[model] >= 0.7) ].index.tolist()
    confidence_scores_matrix[model][MC]=2
        
    #Low
    LC=ub_scores_matrix[model][(ub_scores_matrix[model] < 0.7) & (ub_scores_matrix[model] >=0.5) ].index.tolist()
    confidence_scores_matrix[model][LC]=1
    
    #unknown
    ZC=ub_scores_matrix[model][(ub_scores_matrix[model] <.5) & (ub_scores_matrix[model] >= 0.3) ].index.tolist()
    confidence_scores_matrix[model][ZC]=0
    
    #Negative
    NC=ub_scores_matrix[model][(ub_scores_matrix[model] < 0.3)].index.tolist()
    confidence_scores_matrix[model][NC]=-1
 
```


```python
import multiprocessing as mp
from corda import CORDA
import string

output = mp.Queue(maxsize=24)

# define a example function
def multiple_randomized_analysis(modelC, modelN, Recon, pos):
    cConf_scores=randomize_model(modelC)
    nConf_scores=randomize_model(modelN)


    cConf_scores["biomass_reaction"]=3
    cConf_scores["biomass_DNA"]=3
    cConf_scores["biomass_RNA"]=3
    cConf_scores["biomass_carbohydrate"]=3
    cConf_scores["biomass_lipid"]=3
    cConf_scores["biomass_other"]=3
    cConf_scores["biomass_protein"]=3
    cConf_scores["DM_atp_c_"]=3

    nConf_scores["biomass_reaction"]=3
    nConf_scores["biomass_DNA"]=3
    nConf_scores["biomass_RNA"]=3
    nConf_scores["biomass_carbohydrate"]=3
    nConf_scores["biomass_lipid"]=3
    nConf_scores["biomass_other"]=3
    nConf_scores["biomass_protein"]=3
    nConf_scores["DM_atp_c_"]=3
    
    Recon.objective="biomass_reaction"
    
    metas = ['adp_c', 'atp_c','glc_D_c',  'fru_c', 'nad_c', 'nadh_c','nad_m', 'nadh_m','nadph_c', 'nadph_m', 'nadp_c', 'nadp_m', 'cmp_c', 'HC00342_c', 'glcn_c', 'citr_L_c', 'glyb_c', 'icit_c', '3pg_c', 'accoa_m ->coa_m', 'akg_m', 'e4p_c', 'f6p_c', 'g3p_c', 'g6p_c', 'oaa_m', 'pep_c', 'pyr_c', 'r5p_c', 'succoa_m ->coa_m', 'ala_L_c', 'arg_L_c', 'asn_L_c', 'asp_L_c', 'val_L_c', 'adp_c', 'thr_L_c', 'leu_L_c', 'gln_L_c', 'glu_L_c', 'gly_c', 'pro_D_c', 'pro_L_c', 'ser_L_c', 'ctp_c', 'fdp_c','utp_c', 'pmtcoa_c -> coa_c', 'chsterol_c', 'dag_hs_c', 'tag_hs_c', 'mag_hs_c', 'gthox_c','gthrd_c','ru5p_D_c','crm_hs_c', 'pa_hs_c', 'pe_hs_c', 'ps_hs_c', 'hxan_c', 'creat_c', 'crtn_c', 'chol_c', 'orn_c', '4hpro_LT_c', 's7p_c', 'amp_c', 'udp_c', 'ade_c', 'asp_D_c', 'adn_c', 'ump_c', 'his_L_c', 'tyr_L_c', 'phe_L_c', 'gdp_c', 'ctp_c', 'cys_L_c', 'amet_c', 'cit_c', 'tym_c', 'succ_c', 'CE1936_c', 'gua_c', 'cdp_c', 'g1p_c', 'lys_L_c', 'bhb_c', 'ile_L_c', 'gmp_c', 'dhap_c', 'fum_c', 'mal_L_c', 'spmd_c', 'ala_B_c', 'trp_L_c', 'lac_L_c', 'met_L_c', 'ptrc_c', '4abut_c','0.014 biomass_DNA_c + 0.058 biomass_RNA_c + 0.071 biomass_carbohydrate_c + 0.097 biomass_lipid_c + 0.054 biomass_other_c + 0.706 biomass_protein_c --> ']    

    opt_c= CORDA(model=Recon, confidence=cConf_scores, n=5, met_prod=metas,  penalty_factor=1000 ) 
    opt_n= CORDA(model=Recon, confidence=nConf_scores, n=5, met_prod=metas,  penalty_factor=1000 ) 
    opt_c.build()
    opt_n.build()
    model_c=opt_c.cobra_model(name="C")
    model_n=opt_n.cobra_model(name="N")
     
    targets=get_drugable_targets(model_n, model_c, "biopsys" )
    f="Targets_pf_1000_Biopsys_" + str(pos) + ".tsv"                                                                  
    targets.to_csv(f)

# Setup a list of processes that we want to run
pool = mp.Pool(processes=24)

for i in range(1,101):
    pool.apply_async(multiple_randomized_analysis, args=(conf_MeanCancerBiopsy,conf_MeanNormalBiopsy,Recon2,i))

pool.close()
pool.join()

```

from corda import CORDA

def randomized_analysis(modelC, modelN, Recon):
    cConf_scores=randomize_model(modelC)
    nConf_scores=randomize_model(modelN)
    
    metas = ['adp_c', 'atp_c','glc_D_c',  'fru_c', 'nad_c', 'nadh_c','nad_m', 'nadh_m','nadph_c', 'nadph_m', 'nadp_c', 'nadp_m', 'cmp_c', 'HC00342_c', 'glcn_c', 'citr_L_c', 'glyb_c', 'icit_c', '3pg_c', 'accoa_m ->coa_m', 'akg_m', 'e4p_c', 'f6p_c', 'g3p_c', 'g6p_c', 'oaa_m', 'pep_c', 'pyr_c', 'r5p_c', 'succoa_m ->coa_m', 'ala_L_c', 'arg_L_c', 'asn_L_c', 'asp_L_c', 'val_L_c', 'adp_c', 'thr_L_c', 'leu_L_c', 'gln_L_c', 'glu_L_c', 'gly_c', 'pro_D_c', 'pro_L_c', 'ser_L_c', 'ctp_c', 'fdp_c','utp_c', 'pmtcoa_c -> coa_c', 'chsterol_c', 'dag_hs_c', 'tag_hs_c', 'mag_hs_c', 'gthox_c','gthrd_c','ru5p_D_c','crm_hs_c', 'pa_hs_c', 'pe_hs_c', 'ps_hs_c', 'hxan_c', 'creat_c', 'crtn_c', 'chol_c', 'orn_c', '4hpro_LT_c', 's7p_c', 'amp_c', 'udp_c', 'ade_c', 'asp_D_c', 'adn_c', 'ump_c', 'his_L_c', 'tyr_L_c', 'phe_L_c', 'gdp_c', 'ctp_c', 'cys_L_c', 'amet_c', 'cit_c', 'tym_c', 'succ_c', 'CE1936_c', 'gua_c', 'cdp_c', 'g1p_c', 'lys_L_c', 'bhb_c', 'ile_L_c', 'gmp_c', 'dhap_c', 'fum_c', 'mal_L_c', 'spmd_c', 'ala_B_c', 'trp_L_c', 'lac_L_c', 'met_L_c', 'ptrc_c', '4abut_c','0.014 biomass_DNA_c + 0.058 biomass_RNA_c + 0.071 biomass_carbohydrate_c + 0.097 biomass_lipid_c + 0.054 biomass_other_c + 0.706 biomass_protein_c --> ']    
 
    opt_c= CORDA(model=Recon, confidence=cConf_scores, n=5, met_prod=metas,  penalty_factor=2500 ) 
    opt_n= CORDA(model=Recon, confidence=nConf_scores, n=5, met_prod=metas,  penalty_factor=2500 ) 
    opt_c.build()
    opt_n.build()
    model_c=opt_c.cobra_model(name="C")
    model_n=opt_n.cobra_model(name="N")

    targets=get_drugable_targets(model_n, model_c, "biopsys" )
    targets.to_csv("Targets.tsv",sep="\t")
    targets=targets[ (targets["norm_dise_ratio"] > 1.000001) ].sort_values(by="norm_dise_ratio", ascending=False)
    targets.to_csv("Targets2_pf_2500.tsv",sep="\t")

    
randomized_analysis(modelC=conf_CancerBiopsy,  modelN=conf_NormalBiopsy, Recon=Recon2)



```python
%matplotlib inline

from mpl_toolkits import mplot3d
import numpy as np
import matplotlib.pyplot as plt

fig = plt.figure()
ax = plt.axes(projection='3d')

ax = plt.axes(projection='3d')

# Data for a three-dimensional line
zline = np.linspace(0, 15, 1000)
xline = np.sin(zline)
yline = np.cos(zline)
ax.plot3D(xline, yline, zline, 'gray')

# Data for three-dimensional scattered points
zdata = 15 * np.random.random(100)
xdata = np.sin(zdata) + 0.1 * np.random.randn(100)
ydata = np.cos(zdata) + 0.1 * np.random.randn(100)
ax.scatter3D(xdata, ydata, zdata, c=zdata, cmap='Greens');
```


```python
def f(x, y):
    return np.sin(np.sqrt(x ** 2 + y ** 2))

x = np.linspace(-6, 6, 30)
y = np.linspace(-6, 6, 30)

X, Y = np.meshgrid(x, y)
Z = f(X, Y)
fig = plt.figure()
ax = plt.axes(projection='3d')
ax.contour3D(X, Y, Z, 50, cmap='binary')
ax.set_xlabel('x')
ax.set_ylabel('y')
ax.set_zlabel('z');

```


```python
import pandas as pd
 
```


```python

```
