

```python
import cobra
Recon2 = cobra.io.read_sbml_model("Models/recon2.2.xml")

```

    cobra/io/sbml.py:235 [1;31mUserWarning[0m: M_h_x appears as a reactant and product FAOXC220200x



```python
Recon2.reactions.OIVD3m.gene_reaction_rule='HGNC:2698 and HGNC:987 and HGNC:986 and HGNC:2898 or HGNC:2698 and HGNC:2898 and HGNC:987 and HGNC:986'
Recon2.reactions.OIVD1m.gene_reaction_rule='HGNC:2698 and HGNC:987 and HGNC:986 and HGNC:2898 or HGNC:2698 and HGNC:2898 and HGNC:987 and HGNC:986'
Recon2.reactions.OIVD2m.gene_reaction_rule='HGNC:2698 and HGNC:987 and HGNC:986 and HGNC:2898 or HGNC:2698 and HGNC:2898 and HGNC:987 and HGNC:986'
```


```python
print("Reactions:", len(Recon2.reactions))
print("Metabolites:", len(Recon2.metabolites))
print(Recon2.optimize())
```

    Reactions: 7785
    Metabolites: 5324
    <Solution 555.786 at 0x7fbadfdbe780>



```python
bm = Recon2.reactions.get_by_id("biomass_reaction")
print(bm.build_reaction_string())
print(bm.notes)
```

    0.014 biomass_DNA_c + 0.058 biomass_RNA_c + 0.071 biomass_carbohydrate_c + 0.097 biomass_lipid_c + 0.054 biomass_other_c + 0.706 biomass_protein_c --> 
    {'GENE ASSOCIATION': ['']}



```python
Recon2.optimize()
```




    <Solution 555.786 at 0x7fbadfb810b8>




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




    <Solution 0.100 at 0x7fbadfdbe048>




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
       
    HC=ub_scores_matrix[model][ub_scores_matrix[model] >= .7].index.tolist()
    confidence_scores_matrix[model][HC]=3
        
    #Medium Confidence Rate 
    MC=ub_scores_matrix[model][(ub_scores_matrix[model] <  .7) & (ub_scores_matrix[model] >= 0.5) ].index.tolist()
    confidence_scores_matrix[model][MC]=2
        
    #Low
    LC=ub_scores_matrix[model][(ub_scores_matrix[model] < 0.5) & (ub_scores_matrix[model] >=0.4) ].index.tolist()
    confidence_scores_matrix[model][LC]=1
    
    #unknown
    ZC=ub_scores_matrix[model][(ub_scores_matrix[model] <.4) & (ub_scores_matrix[model] >= 0.2) ].index.tolist()
    confidence_scores_matrix[model][ZC]=0
    
    #Negative
    NC=ub_scores_matrix[model][(ub_scores_matrix[model] < 0.2 )].index.tolist()
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
#n=5  penalty_factor=900  <Solution 0.056 at 0x7fabebb0f278>


from corda import CORDA

opt_MeanCancerBiopsy = CORDA(model=Recon2, confidence=conf_MeanCancerBiopsy, n=5, met_prod=metas,  penalty_factor=1000 ) 
opt_MeanCancerBiopsy.build()
print(opt_MeanCancerBiopsy)
model_MeanCancerBiopsy=opt_MeanCancerBiopsy.cobra_model(name="MeanCancerBiopsy")
print(model_MeanCancerBiopsy.optimize())
```

    build status: reconstruction complete
    Inc. reactions: 2041/7886
     - unclear: 912/3621
     - exclude: 218/2637
     - low and medium: 226/845
     - high: 685/783
    
    <Solution 0.056 at 0x7fbade0569e8>
    CPU times: user 10min 23s, sys: 128 ms, total: 10min 23s
    Wall time: 10min 23s



```python
model_MeanCancerBiopsy.summary(fva=True)
```

    IN FLUXES                                      OUT FLUXES                                      OBJECTIVES
    ---------------------------------------------  ----------------------------------------------  -----------------------
    id                   Flux  Range               id                   Flux  Range                biomass_reac...  0.0563
    ---------------  --------  ------------------  ---------------  --------  -------------------
    h2o_e            0.499     [-3.46, 7.85]       h_e              1.17      [-1, 6.3]
    glyc_e           0.32      [0, 1.8]            lac_L_e          0.777     [-0.01, 8.1]
    cit_e            0.31      [-0.845, 1]         ac_e             0.353     [0, 4.72]
    glc_D_e          0.149     [0.149, -0.899]     mal_L_e          0.32      [-0.01, 3.39]
    ser_L_e          0.0495    [0.0124, 1]         lac_D_e          0.297     [0, 8.11]
    gam_e            0.0384    [0, 1]              pi_e             0.127     [-0.106, 0.614]
    lys_L_e          0.0333    [0.0333, 0.0333]    nh4_e            0.0736    [-0.01, 4.42]
    glu_L_e          0.0314    [0, 1]              urate_e          0.0462    [0, 0.0506]
    leu_L_e          0.0307    [0.0207, 0.15]      glyc_S_e         0.037     [0, 0.988]
    gly_e            0.0303    [-0.256, 0.45]      co2_e            0.027     [-0.01, 4.8]
    pro_L_e          0.0232    [0.0232, 0.0232]    abt_e            0.02      [0, 0.929]
    arg_L_e          0.0202    [0.0102, 0.0202]    ribflv_e         0.01      [0, 0.02]
    asp_L_e          0.0198    [0, 0.225]          glygn4_e         0.01      [0, 0.01]
    val_L_e          0.0198    [0.00983, 0.18]     strch2_e         0.01      [0, 0.01]
    gln_L_e          0.0183    [0.18, -1.93]       ascb_L_e         0.01      [-0.01, 0.01]
    thr_L_e          0.0176    [0.0176, 0.0176]    adn_e            0.01      [0.01, -1]
    3mop_e           0.0161    [-0.134, 1]         ala_L_e          0.00861   [-0.3, 1.06]
    asn_L_e          0.0157    [0.0157, 0.0157]    dcmp_e           0.00727   [0.00727, -0.01]
    phe_L_e          0.0146    [0.0146, 0.0146]    Rtotal_e         0.00672   [-0.01, 0.019]
    glyc3p_e         0.0108    [-0.01, 0.434]      c51crn_e         0.00548   [0, 0.01]
    ps_hs_e          0.01      [0.01, 0.01]        am4n9cs_e        0.00475   [0, 0.01]
    fad_e            0.01      [0, 0.01]           ala_B_e          0.00452   [-0.01, 0.982]
    glygn2_e         0.01      [0, 0.01]           q10h2_e          0.00452   [0.01, -0.01]
    o2s_e            0.01      [0, 0.01]           dttp_m           0.00315   [0, 0.0163]
    strch1_e         0.01      [0, 0.01]           acetone_e        0.00312   [0, 0.0607]
    dctp_n           0.01      [-0.00727, 0.01]    vacc_e           0.00135   [-0.01, 0.124]
    dhdascb_e        0.01      [0.01, -0.01]       dag_hs_e         0.000984  [-0.000984, 0.0159]
    datp_n           0.01      [0.01, -0.0109]     fum_e            0         [-1, 2.4]
    dgtp_n           0.01      [0.01, -0.0109]     5oxpro_e         0         [0, 2.38]
    atp_e            0.01      [0.01, -0.0232]     akg_e            0         [-0.01, 2.25]
    xmp_e            0.01      [0.01, -0.0332]     fald_e           0         [-1, 1.96]
    acac_e           0.01      [0.01, -0.0507]     HC00342_e        0         [-0.01, 1.83]
    o2_e             0.01      [0.01, -0.114]      slfcys_e         0         [-0.01, 1.01]
    h2o2_e           0.01      [0.01, -0.125]      ile_L_e          0         [-0.15, 0.984]
    ha_e             0.01      [0.01, -0.199]      xylt_e           0         [-0.01, 0.919]
    ppa_e            0.01      [0.01, -0.361]      gthrd_e          0         [-0.281, 0.617]
    ins_e            0.01      [0.01, -1]          Lcystin_c        0         [-0.01, 0.504]
    pyr_e            0.01      [0.01, -2.95]       fe3_e            0         [-0.01, 0.496]
    tyr_L_e          0.00898   [0.00898, 0.00898]  ha_pre1_e        0         [-0.01, 0.395]
    met_L_e          0.00861   [0.00861, 0.0418]   gal_e            0         [0, 0.315]
    lpchol_hs_e      0.00771   [0, 0.01]           oxa_e            0         [0, 0.265]
    his_L_e          0.00711   [0, 0.00711]        3mob_e           0         [-0.01, 0.16]
    utp_e            0.00616   [0.00301, 0.01]     HC00250_e        0         [0, 0.143]
    crn_e            0.00548   [0, 0.01]           gthox_e          0         [0, 0.135]
    pe_hs_e          0.00508   [-0.00688, 0.01]    taur_c           0         [-0.01, 0.134]
    am9csa_e         0.00475   [0, 0.01]           hdcea_e          0         [0, 0.12]
    q10_e            0.00452   [-0.01, 0.01]       4mop_e           0         [-0.01, 0.119]
    biomass_other_c  0.00304   [0.00304, 0.00304]  ocdcea_e         0         [0, 0.108]
    cys_L_e          0.00262   [-0.0274, 1]        fuc_L_e          0         [0, 0.0949]
    datp_m           0.00256   [0.00699, -0.0209]  nrvnc_e          0         [0, 0.0948]
    for_e            0.002     [0.00585, -1.48]    anth_c           0         [0, 0.0811]
    strdnc_e         0.00135   [0, 0.01]           pheme_e          0         [0, 0.0567]
    inost_e          0.00131   [0.00131, 0.00131]  Ser_Thr_l        0         [-0.01, 0.05]
    sph1p_e          0.000984  [0.000984, 0.01]    ahcys_e          0         [0, 0.0332]
    pglyc_hs_e       0.00082   [0.00082, 0.00082]  acald_e          0         [-0.01, 0.0309]
    trp_L_e          0.000748  [0.000748, 0.0818]  adrn_e           0         [0, 0.0301]
    dttp_n           0.000736  [-0.00626, 0.01]    clpnd_e          0         [0, 0.03]
    meoh_e           0         [-1, 1.96]          cholate_e        0         [0, 0.0256]
    so4_e            0         [-0.13, 1]          bilglcur_e       0         [0, 0.024]
    fe2_e            0         [-0.01, 0.499]      co_e             0         [0, 0.024]
    arab_L_e         0         [0, 0.179]          prostgh2_e       0         [0, 0.024]
    hdca_e           0         [0, 0.131]          dgtp_m           0         [0, 0.0209]
    ocdca_e          0         [0, 0.125]          1mncam_e         0         [0, 0.02]
    lnlc_e           0         [0, 0.0538]         nac_e            0         [0, 0.02]
    3hpvs_e          0         [0, 0.01]           tmndnc_e         0         [-0.01, 0.02]
    4hpro_LT_m       0         [0, 0.01]           dctp_m           0         [0, 0.0173]
    C02470_e         0         [0, 0.01]           glyb_e           0         [0, 0.0172]
    HC00822_e        0         [0, 0.01]           gchola_e         0         [-0.01, 0.0156]
    T_antigen_g      0         [0, 0.01]           tchola_e         0         [-0.01, 0.0154]
    citr_L_c         0         [0, 0.01]           core5_g          0         [0, 0.012]
    crvnc_e          0         [0, 0.01]           core7_g          0         [0, 0.012]
    csa_e            0         [0, 0.01]           3ddcrn_e         0         [0, 0.01]
    dcsptn1_e        0         [0, 0.01]           3deccrn_e        0         [0, 0.01]
    ddca_e           0         [0, 0.01]           3hpvstet_e       0         [0, 0.01]
    etoh_e           0         [0, 0.01]           3ivcrn_e         0         [0, 0.01]
    fmn_e            0         [0, 0.01]           3octdeccrn_e     0         [0, 0.01]
    gluala_e         0         [0, 0.01]           3tdcrn_e         0         [0, 0.01]
    leuktrA4_e       0         [0, 0.01]           3ttetddcoacrn_e  0         [0, 0.01]
    octa_e           0         [0, 0.01]           am4ncs_e         0         [0, 0.01]
    retinol_e        0         [0, 0.01]           c6crn_e          0         [0, 0.01]
    sbt_D_e          0         [0, 0.01]           c8crn_e          0         [0, 0.01]
    tacr_e           0         [0, 0.01]           ddeccrn_e        0         [0, 0.01]
    whhdca_e         0         [0, 0.01]           epoxtac_e        0         [0, 0.01]
    pe_hs_r          0         [-0.00688, 0.01]    ivcrn_e          0         [0, 0.01]
    Asn_X_Ser_Thr_l  0         [-0.01, 0.01]       leuktrB4_e       0         [0, 0.01]
    eicostet_e       0         [-0.01, 0.01]       leuktrE4_e       0         [0, 0.01]
    carn_e           0         [0, 0.00711]        retfa_e          0         [0, 0.01]
                                                   retn_e           0         [0, 0.01]
                                                   HC01609_e        0         [-0.01, 0.01]
                                                   HC01610_e        0         [-0.01, 0.01]
                                                   cmp_e            0         [-0.01, 0.01]
                                                   ctp_e            0         [-0.01, 0.01]
                                                   dtmp_e           0         [-0.01, 0.01]
                                                   dttp_e           0         [-0.01, 0.01]
                                                   n2m2nmasn_e      0         [-0.01, 0.01]
                                                   ncam_c           0         [-0.01, 0.01]
                                                   fuc14galacgl...  0         [0, 0.00902]
                                                   sphs1p_e         0         [0, 0.00902]
                                                   man_e            0         [0.03, -0.114]



```python
%%time
#n=5  penalty_factor=900  <Solution 0.051 at 0x7fabebf29160>
#n=5  penalty_factor=500  <Solution 0.050 at 0x7fabebf29160>

from corda import CORDA


opt_MeanNormalBiopsy = CORDA(model=Recon2, confidence=conf_MeanNormalBiopsy, n=5, met_prod=metas,  penalty_factor=1000) 
opt_MeanNormalBiopsy.build()
print(opt_MeanNormalBiopsy)
model_MeanNormalBiopsy=opt_MeanNormalBiopsy.cobra_model(name="MeanNormalBiopsy")
print(model_MeanNormalBiopsy.optimize())
```


```python
model_MeanNormalBiopsy.summary(fva=True)
```


```python
cobra.io.write_sbml_model(model_MeanCancerBiopsy, "Full_SCC_model_1218_western_diet_n5.sbml")
cobra.io.write_sbml_model(model_MeanNormalBiopsy, "Full_Normal_model_1218_western_diet_n5.sbml")
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
R_Biopsies["dise_prolif_ratio"]=R_Biopsies["del_dise_flux"]/R_Biopsies["dise_flux"]
R_Biopsies["norm_prolif_ratio"]=R_Biopsies["del_norm_flux"]/R_Biopsies["norm_flux"]
R_Biopsies["total_ratio"]=R_Biopsies["dise_prolif_ratio"]/R_Biopsies["norm_prolif_ratio"]

```


```python
R_Biopsies[R_Biopsies["dise_prolif_ratio"]<.999].sort_values(by="total_ratio", ascending=True)

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

for i in range(1,49):
    pool.apply_async(multiple_randomized_analysis, args=(conf_MeanCancerBiopsy,conf_MeanNormalBiopsy,Recon2,i))

pool.close()
pool.join()

```

    Common reactions size 1120
    Normal unique reactions size 1059
    Disease unique reactions size 1100


    cobra/util/solver.py:404 [1;31mUserWarning[0m: solver status is 'infeasible'


    Common reactions size 1169
    Normal unique reactions size 960
    Disease unique reactions size 1116


    cobra/util/solver.py:404 [1;31mUserWarning[0m: solver status is 'infeasible'


    Common reactions size 1102
    Normal unique reactions size 1091
    Disease unique reactions size 1226


    cobra/util/solver.py:404 [1;31mUserWarning[0m: solver status is 'infeasible'


    Common reactions size 1131
    Normal unique reactions size 1085
    Disease unique reactions size 1162


    cobra/util/solver.py:404 [1;31mUserWarning[0m: solver status is 'infeasible'


    Common reactions size 1019
    Normal unique reactions size 987
    Disease unique reactions size 1279


    cobra/util/solver.py:404 [1;31mUserWarning[0m: solver status is 'infeasible'


    Common reactions size 1121
    Normal unique reactions size 1024
    Disease unique reactions size 1206


    cobra/util/solver.py:404 [1;31mUserWarning[0m: solver status is 'infeasible'


    Common reactions size 1154
    Normal unique reactions size 1044
    Disease unique reactions size 1248


    cobra/util/solver.py:404 [1;31mUserWarning[0m: solver status is 'infeasible'


    Common reactions size 1121
    Normal unique reactions size 1013
    Disease unique reactions size 1323


    cobra/util/solver.py:404 [1;31mUserWarning[0m: solver status is 'infeasible'


    Common reactions size 1123
    Normal unique reactions size 1089
    Disease unique reactions size 1278


    cobra/util/solver.py:404 [1;31mUserWarning[0m: solver status is 'infeasible'


    Common reactions size 1104
    Normal unique reactions size 1095
    Disease unique reactions size 1279
    Common reactions size 1092
    Normal unique reactions size 1077
    Disease unique reactions size 1283


    cobra/util/solver.py:404 [1;31mUserWarning[0m: solver status is 'infeasible'


    Common reactions size 1150
    Normal unique reactions size 1120
    Disease unique reactions size 1185


    cobra/util/solver.py:404 [1;31mUserWarning[0m: solver status is 'infeasible'


    Common reactions size 1080
    Normal unique reactions size 967
    Disease unique reactions size 1237
    Common reactions size 1101
    Normal unique reactions size 1090
    Disease unique reactions size 1258
    Common reactions size 1266
    Normal unique reactions size 1116
    Disease unique reactions size 1139
    Common reactions size 1092
    Normal unique reactions size 1032
    Disease unique reactions size 1298


    cobra/util/solver.py:404 [1;31mUserWarning[0m: solver status is 'infeasible'
    cobra/util/solver.py:404 [1;31mUserWarning[0m: solver status is 'infeasible'


    Common reactions size 1149
    Normal unique reactions size 1148
    Disease unique reactions size 1165


    cobra/util/solver.py:404 [1;31mUserWarning[0m: solver status is 'infeasible'


    Common reactions size 1136
    Normal unique reactions size 1119
    Disease unique reactions size 1323


    cobra/util/solver.py:404 [1;31mUserWarning[0m: solver status is 'infeasible'


    Common reactions size 1114
    Normal unique reactions size 939
    Disease unique reactions size 1315


    cobra/util/solver.py:404 [1;31mUserWarning[0m: solver status is 'infeasible'


    Common reactions size 1186
    Normal unique reactions size 1075
    Disease unique reactions size 1189


    cobra/util/solver.py:404 [1;31mUserWarning[0m: solver status is 'infeasible'


    Common reactions size 1014
    Normal unique reactions size 1130
    Disease unique reactions size 1294


    cobra/util/solver.py:404 [1;31mUserWarning[0m: solver status is 'infeasible'


    Common reactions size 1151
    Normal unique reactions size 1065
    Disease unique reactions size 1299


    cobra/util/solver.py:404 [1;31mUserWarning[0m: solver status is 'infeasible'
    cobra/util/solver.py:404 [1;31mUserWarning[0m: solver status is 'infeasible'


    Common reactions size 1053
    Normal unique reactions size 1113
    Disease unique reactions size 1126


    cobra/util/solver.py:404 [1;31mUserWarning[0m: solver status is 'infeasible'


    Common reactions size 993
    Normal unique reactions size 1139
    Disease unique reactions size 1254


    cobra/util/solver.py:404 [1;31mUserWarning[0m: solver status is 'infeasible'


    Common reactions size 1041
    Normal unique reactions size 1112
    Disease unique reactions size 1212


    cobra/util/solver.py:404 [1;31mUserWarning[0m: solver status is 'infeasible'


    Common reactions size 1054
    Normal unique reactions size 1137
    Disease unique reactions size 1266


    Process ForkPoolWorker-13:
    Process ForkPoolWorker-8:
    Process ForkPoolWorker-22:
    Process ForkPoolWorker-3:
    Process ForkPoolWorker-20:
    Process ForkPoolWorker-21:
    Process ForkPoolWorker-19:
    Process ForkPoolWorker-10:
    Process ForkPoolWorker-18:
    Process ForkPoolWorker-2:
    Process ForkPoolWorker-16:
    Process ForkPoolWorker-7:
    Process ForkPoolWorker-6:
    Process ForkPoolWorker-1:
    Process ForkPoolWorker-4:
    Process ForkPoolWorker-14:
    Process ForkPoolWorker-24:
    Process ForkPoolWorker-12:
    Process ForkPoolWorker-11:
    Traceback (most recent call last):
      File "/opt/conda/lib/python3.5/multiprocessing/process.py", line 249, in _bootstrap
        self.run()
      File "/opt/conda/lib/python3.5/multiprocessing/process.py", line 93, in run
        self._target(*self._args, **self._kwargs)
      File "/opt/conda/lib/python3.5/multiprocessing/pool.py", line 119, in worker
        result = (True, func(*args, **kwds))
      File "<ipython-input-23-9847a66c3d0b>", line 38, in multiple_randomized_analysis
        opt_n.build()
      File "/opt/conda/lib/python3.5/site-packages/corda/corda.py", line 236, in build
        need = self.associated(include)
      File "/opt/conda/lib/python3.5/site-packages/corda/corda.py", line 182, in associated
        self.__zero_objective()
      File "/opt/conda/lib/python3.5/site-packages/corda/corda.py", line 138, in __zero_objective
        {v: 0 for v in self.model.variables})
    Traceback (most recent call last):
    Traceback (most recent call last):
    Traceback (most recent call last):
    Traceback (most recent call last):
    Traceback (most recent call last):
      File "/opt/conda/lib/python3.5/site-packages/corda/corda.py", line 138, in <dictcomp>
        {v: 0 for v in self.model.variables})
      File "/opt/conda/lib/python3.5/multiprocessing/process.py", line 249, in _bootstrap
        self.run()
    Traceback (most recent call last):
    Traceback (most recent call last):
    Traceback (most recent call last):
    Traceback (most recent call last):
      File "/opt/conda/lib/python3.5/multiprocessing/process.py", line 249, in _bootstrap
        self.run()
      File "/opt/conda/lib/python3.5/multiprocessing/process.py", line 249, in _bootstrap
        self.run()
      File "/opt/conda/lib/python3.5/multiprocessing/process.py", line 249, in _bootstrap
        self.run()
      File "/opt/conda/lib/python3.5/multiprocessing/process.py", line 249, in _bootstrap
        self.run()
    KeyboardInterrupt
      File "/opt/conda/lib/python3.5/multiprocessing/process.py", line 249, in _bootstrap
        self.run()
      File "/opt/conda/lib/python3.5/multiprocessing/process.py", line 249, in _bootstrap
        self.run()
      File "/opt/conda/lib/python3.5/multiprocessing/process.py", line 93, in run
        self._target(*self._args, **self._kwargs)
      File "/opt/conda/lib/python3.5/multiprocessing/process.py", line 249, in _bootstrap
        self.run()
      File "/opt/conda/lib/python3.5/multiprocessing/process.py", line 249, in _bootstrap
        self.run()
      File "/opt/conda/lib/python3.5/multiprocessing/process.py", line 93, in run
        self._target(*self._args, **self._kwargs)
      File "/opt/conda/lib/python3.5/multiprocessing/process.py", line 93, in run
        self._target(*self._args, **self._kwargs)
      File "/opt/conda/lib/python3.5/multiprocessing/process.py", line 93, in run
        self._target(*self._args, **self._kwargs)
      File "/opt/conda/lib/python3.5/multiprocessing/process.py", line 93, in run
        self._target(*self._args, **self._kwargs)
      File "/opt/conda/lib/python3.5/multiprocessing/pool.py", line 119, in worker
        result = (True, func(*args, **kwds))
      File "/opt/conda/lib/python3.5/multiprocessing/process.py", line 93, in run
        self._target(*self._args, **self._kwargs)
      File "/opt/conda/lib/python3.5/multiprocessing/process.py", line 93, in run
        self._target(*self._args, **self._kwargs)
      File "/opt/conda/lib/python3.5/multiprocessing/process.py", line 93, in run
        self._target(*self._args, **self._kwargs)
      File "/opt/conda/lib/python3.5/multiprocessing/process.py", line 93, in run
        self._target(*self._args, **self._kwargs)
      File "/opt/conda/lib/python3.5/multiprocessing/pool.py", line 119, in worker
        result = (True, func(*args, **kwds))
      File "/opt/conda/lib/python3.5/multiprocessing/pool.py", line 119, in worker
        result = (True, func(*args, **kwds))
      File "/opt/conda/lib/python3.5/multiprocessing/pool.py", line 119, in worker
        result = (True, func(*args, **kwds))
      File "/opt/conda/lib/python3.5/multiprocessing/pool.py", line 119, in worker
        result = (True, func(*args, **kwds))
      File "/opt/conda/lib/python3.5/multiprocessing/pool.py", line 119, in worker
        result = (True, func(*args, **kwds))
      File "/opt/conda/lib/python3.5/multiprocessing/pool.py", line 119, in worker
        result = (True, func(*args, **kwds))
      File "/opt/conda/lib/python3.5/multiprocessing/pool.py", line 119, in worker
        result = (True, func(*args, **kwds))
      File "<ipython-input-23-9847a66c3d0b>", line 38, in multiple_randomized_analysis
        opt_n.build()
      File "<ipython-input-23-9847a66c3d0b>", line 38, in multiple_randomized_analysis
        opt_n.build()
      File "/opt/conda/lib/python3.5/multiprocessing/pool.py", line 119, in worker
        result = (True, func(*args, **kwds))
      File "<ipython-input-23-9847a66c3d0b>", line 38, in multiple_randomized_analysis
        opt_n.build()
      File "<ipython-input-23-9847a66c3d0b>", line 38, in multiple_randomized_analysis
        opt_n.build()
      File "<ipython-input-23-9847a66c3d0b>", line 38, in multiple_randomized_analysis
        opt_n.build()
      File "<ipython-input-23-9847a66c3d0b>", line 38, in multiple_randomized_analysis
        opt_n.build()
      File "<ipython-input-23-9847a66c3d0b>", line 38, in multiple_randomized_analysis
        opt_n.build()
      File "<ipython-input-23-9847a66c3d0b>", line 38, in multiple_randomized_analysis
        opt_n.build()
      File "/opt/conda/lib/python3.5/site-packages/corda/corda.py", line 272, in build
        penalize_medium=False)
      File "/opt/conda/lib/python3.5/site-packages/corda/corda.py", line 272, in build
        penalize_medium=False)
      File "<ipython-input-23-9847a66c3d0b>", line 38, in multiple_randomized_analysis
        opt_n.build()
      File "/opt/conda/lib/python3.5/site-packages/corda/corda.py", line 272, in build
        penalize_medium=False)
      File "/opt/conda/lib/python3.5/site-packages/corda/corda.py", line 272, in build
        penalize_medium=False)
      File "/opt/conda/lib/python3.5/site-packages/corda/corda.py", line 272, in build
        penalize_medium=False)
      File "/opt/conda/lib/python3.5/site-packages/corda/corda.py", line 272, in build
        penalize_medium=False)
      File "/opt/conda/lib/python3.5/site-packages/corda/corda.py", line 272, in build
        penalize_medium=False)
      File "/opt/conda/lib/python3.5/site-packages/corda/corda.py", line 182, in associated
        self.__zero_objective()
      File "/opt/conda/lib/python3.5/site-packages/corda/corda.py", line 272, in build
        penalize_medium=False)
      File "/opt/conda/lib/python3.5/site-packages/corda/corda.py", line 244, in build
        need = self.associated(include, penalize_medium=False)
      File "/opt/conda/lib/python3.5/site-packages/corda/corda.py", line 182, in associated
        self.__zero_objective()
      File "/opt/conda/lib/python3.5/site-packages/corda/corda.py", line 182, in associated
        self.__zero_objective()
      File "/opt/conda/lib/python3.5/site-packages/corda/corda.py", line 182, in associated
        self.__zero_objective()
      File "/opt/conda/lib/python3.5/site-packages/corda/corda.py", line 182, in associated
        self.__zero_objective()
      File "/opt/conda/lib/python3.5/site-packages/corda/corda.py", line 182, in associated
        self.__zero_objective()
      File "/opt/conda/lib/python3.5/site-packages/corda/corda.py", line 195, in associated
        self.__corda_objective(pen)
      File "/opt/conda/lib/python3.5/site-packages/corda/corda.py", line 138, in __zero_objective
        {v: 0 for v in self.model.variables})
      File "/opt/conda/lib/python3.5/site-packages/corda/corda.py", line 195, in associated
        self.__corda_objective(pen)
      File "/opt/conda/lib/python3.5/site-packages/corda/corda.py", line 138, in __zero_objective
        {v: 0 for v in self.model.variables})
      File "/opt/conda/lib/python3.5/site-packages/corda/corda.py", line 195, in associated
        self.__corda_objective(pen)
      File "/opt/conda/lib/python3.5/site-packages/corda/corda.py", line 138, in __zero_objective
        {v: 0 for v in self.model.variables})
      File "/opt/conda/lib/python3.5/site-packages/corda/corda.py", line 138, in __zero_objective
        {v: 0 for v in self.model.variables})
      File "/opt/conda/lib/python3.5/site-packages/corda/corda.py", line 138, in __zero_objective
        {v: 0 for v in self.model.variables})
      File "/opt/conda/lib/python3.5/site-packages/corda/corda.py", line 138, in __zero_objective
        {v: 0 for v in self.model.variables})
      File "/opt/conda/lib/python3.5/site-packages/optlang/glpk_interface.py", line 352, in set_linear_coefficients
        glp_set_obj_coef(self.problem.problem, variable._index, float(coefficient))
      File "/opt/conda/lib/python3.5/site-packages/corda/corda.py", line 134, in __corda_objective
        self.model.objective.set_linear_coefficients(pen)
      File "/opt/conda/lib/python3.5/site-packages/optlang/glpk_interface.py", line 352, in set_linear_coefficients
        glp_set_obj_coef(self.problem.problem, variable._index, float(coefficient))
      File "/opt/conda/lib/python3.5/site-packages/corda/corda.py", line 134, in __corda_objective
        self.model.objective.set_linear_coefficients(pen)
      File "/opt/conda/lib/python3.5/site-packages/corda/corda.py", line 134, in __corda_objective
        self.model.objective.set_linear_coefficients(pen)
      File "/opt/conda/lib/python3.5/site-packages/optlang/glpk_interface.py", line 352, in set_linear_coefficients
        glp_set_obj_coef(self.problem.problem, variable._index, float(coefficient))
      File "/opt/conda/lib/python3.5/site-packages/optlang/glpk_interface.py", line 87, in _index
        i = glp_find_col(self.problem.problem, str(self.name))
    KeyboardInterrupt
      File "/opt/conda/lib/python3.5/site-packages/optlang/glpk_interface.py", line 352, in set_linear_coefficients
        glp_set_obj_coef(self.problem.problem, variable._index, float(coefficient))
      File "/opt/conda/lib/python3.5/site-packages/optlang/glpk_interface.py", line 352, in set_linear_coefficients
        glp_set_obj_coef(self.problem.problem, variable._index, float(coefficient))
      File "/opt/conda/lib/python3.5/site-packages/optlang/glpk_interface.py", line 352, in set_linear_coefficients
        glp_set_obj_coef(self.problem.problem, variable._index, float(coefficient))
    KeyboardInterrupt
    Traceback (most recent call last):
    Traceback (most recent call last):
      File "/opt/conda/lib/python3.5/site-packages/optlang/glpk_interface.py", line 352, in set_linear_coefficients
        glp_set_obj_coef(self.problem.problem, variable._index, float(coefficient))
      File "/opt/conda/lib/python3.5/site-packages/optlang/glpk_interface.py", line 352, in set_linear_coefficients
        glp_set_obj_coef(self.problem.problem, variable._index, float(coefficient))
    KeyboardInterrupt
      File "/opt/conda/lib/python3.5/site-packages/optlang/glpk_interface.py", line 86, in _index
        if self.problem is not None:
      File "/opt/conda/lib/python3.5/site-packages/optlang/glpk_interface.py", line 87, in _index
        i = glp_find_col(self.problem.problem, str(self.name))
      File "/opt/conda/lib/python3.5/site-packages/optlang/glpk_interface.py", line 87, in _index
        i = glp_find_col(self.problem.problem, str(self.name))
    Traceback (most recent call last):
    Traceback (most recent call last):
    Traceback (most recent call last):
    Traceback (most recent call last):
    Traceback (most recent call last):
      File "/opt/conda/lib/python3.5/multiprocessing/process.py", line 249, in _bootstrap
        self.run()
      File "/opt/conda/lib/python3.5/multiprocessing/process.py", line 249, in _bootstrap
        self.run()
      File "/opt/conda/lib/python3.5/site-packages/optlang/glpk_interface.py", line 87, in _index
        i = glp_find_col(self.problem.problem, str(self.name))
      File "/opt/conda/lib/python3.5/site-packages/optlang/glpk_interface.py", line 86, in _index
        if self.problem is not None:
    KeyboardInterrupt
    KeyboardInterrupt
      File "/opt/conda/lib/python3.5/site-packages/optlang/glpk_interface.py", line 352, in set_linear_coefficients
        glp_set_obj_coef(self.problem.problem, variable._index, float(coefficient))
    KeyboardInterrupt
      File "/opt/conda/lib/python3.5/multiprocessing/process.py", line 249, in _bootstrap
        self.run()
      File "/opt/conda/lib/python3.5/multiprocessing/process.py", line 249, in _bootstrap
        self.run()
      File "/opt/conda/lib/python3.5/multiprocessing/process.py", line 249, in _bootstrap
        self.run()
      File "/opt/conda/lib/python3.5/multiprocessing/process.py", line 93, in run
        self._target(*self._args, **self._kwargs)
      File "/opt/conda/lib/python3.5/multiprocessing/process.py", line 249, in _bootstrap
        self.run()
      File "/opt/conda/lib/python3.5/multiprocessing/process.py", line 249, in _bootstrap
        self.run()
    Traceback (most recent call last):
      File "/opt/conda/lib/python3.5/multiprocessing/process.py", line 93, in run
        self._target(*self._args, **self._kwargs)
    KeyboardInterrupt
      File "/opt/conda/lib/python3.5/site-packages/optlang/glpk_interface.py", line 87, in _index
        i = glp_find_col(self.problem.problem, str(self.name))
    KeyboardInterrupt
      File "/opt/conda/lib/python3.5/multiprocessing/process.py", line 93, in run
        self._target(*self._args, **self._kwargs)
      File "/opt/conda/lib/python3.5/multiprocessing/process.py", line 93, in run
        self._target(*self._args, **self._kwargs)
      File "/opt/conda/lib/python3.5/multiprocessing/pool.py", line 119, in worker
        result = (True, func(*args, **kwds))
      File "/opt/conda/lib/python3.5/multiprocessing/process.py", line 93, in run
        self._target(*self._args, **self._kwargs)
      File "/opt/conda/lib/python3.5/multiprocessing/process.py", line 93, in run
        self._target(*self._args, **self._kwargs)
      File "/opt/conda/lib/python3.5/multiprocessing/pool.py", line 119, in worker
        result = (True, func(*args, **kwds))
      File "/opt/conda/lib/python3.5/multiprocessing/process.py", line 249, in _bootstrap
        self.run()
    KeyboardInterrupt
      File "/opt/conda/lib/python3.5/multiprocessing/pool.py", line 119, in worker
        result = (True, func(*args, **kwds))
      File "/opt/conda/lib/python3.5/multiprocessing/pool.py", line 119, in worker
        result = (True, func(*args, **kwds))
      File "<ipython-input-23-9847a66c3d0b>", line 38, in multiple_randomized_analysis
        opt_n.build()
      File "/opt/conda/lib/python3.5/multiprocessing/pool.py", line 119, in worker
        result = (True, func(*args, **kwds))
      File "<ipython-input-23-9847a66c3d0b>", line 42, in multiple_randomized_analysis
        targets=get_drugable_targets(model_n, model_c, "biopsys" )
      File "<ipython-input-23-9847a66c3d0b>", line 38, in multiple_randomized_analysis
        opt_n.build()
      File "/opt/conda/lib/python3.5/multiprocessing/process.py", line 93, in run
        self._target(*self._args, **self._kwargs)
      File "/opt/conda/lib/python3.5/site-packages/corda/corda.py", line 272, in build
        penalize_medium=False)
      File "<ipython-input-18-b299aa21be4a>", line 32, in get_drugable_targets
        nfba=nmodel.optimize()
      File "<ipython-input-23-9847a66c3d0b>", line 42, in multiple_randomized_analysis
        targets=get_drugable_targets(model_n, model_c, "biopsys" )
      File "/opt/conda/lib/python3.5/site-packages/corda/corda.py", line 236, in build
        need = self.associated(include)
      File "/opt/conda/lib/python3.5/multiprocessing/pool.py", line 119, in worker
        result = (True, func(*args, **kwds))
      File "/opt/conda/lib/python3.5/site-packages/cobra/core/model.py", line 840, in optimize
        self.solver.optimize()
      File "/opt/conda/lib/python3.5/site-packages/corda/corda.py", line 196, in associated
        sol = self.model.solver.optimize()
      File "/opt/conda/lib/python3.5/site-packages/corda/corda.py", line 196, in associated
        sol = self.model.solver.optimize()
    Process ForkPoolWorker-9:
      File "<ipython-input-18-b299aa21be4a>", line 32, in get_drugable_targets
        nfba=nmodel.optimize()
      File "/opt/conda/lib/python3.5/site-packages/optlang/interface.py", line 1471, in optimize
        status = self._optimize()
      File "<ipython-input-23-9847a66c3d0b>", line 38, in multiple_randomized_analysis
        opt_n.build()
      File "/opt/conda/lib/python3.5/site-packages/optlang/interface.py", line 1468, in optimize
        status = self._optimize()
      File "/opt/conda/lib/python3.5/multiprocessing/process.py", line 93, in run
        self._target(*self._args, **self._kwargs)
      File "/opt/conda/lib/python3.5/multiprocessing/pool.py", line 119, in worker
        result = (True, func(*args, **kwds))
      File "/opt/conda/lib/python3.5/site-packages/cobra/core/model.py", line 840, in optimize
        self.solver.optimize()
      File "/opt/conda/lib/python3.5/site-packages/optlang/interface.py", line 1468, in optimize
        status = self._optimize()
      File "/opt/conda/lib/python3.5/site-packages/optlang/glpk_interface.py", line 669, in _optimize
        status = self._run_glp_simplex()
      File "/opt/conda/lib/python3.5/site-packages/optlang/glpk_interface.py", line 669, in _optimize
        status = self._run_glp_simplex()
      File "/opt/conda/lib/python3.5/site-packages/optlang/glpk_interface.py", line 643, in _run_glp_simplex
        return_value = glp_simplex(self.problem, self.configuration._smcp)
      File "/opt/conda/lib/python3.5/multiprocessing/pool.py", line 119, in worker
        result = (True, func(*args, **kwds))
      File "<ipython-input-23-9847a66c3d0b>", line 38, in multiple_randomized_analysis
        opt_n.build()
      File "/opt/conda/lib/python3.5/site-packages/optlang/glpk_interface.py", line 669, in _optimize
        status = self._run_glp_simplex()
      File "/opt/conda/lib/python3.5/site-packages/optlang/interface.py", line 1471, in optimize
        status = self._optimize()
    Traceback (most recent call last):
      File "/opt/conda/lib/python3.5/site-packages/optlang/glpk_interface.py", line 643, in _run_glp_simplex
        return_value = glp_simplex(self.problem, self.configuration._smcp)
    KeyboardInterrupt
      File "/opt/conda/lib/python3.5/site-packages/optlang/glpk_interface.py", line 643, in _run_glp_simplex
        return_value = glp_simplex(self.problem, self.configuration._smcp)
      File "<ipython-input-23-9847a66c3d0b>", line 38, in multiple_randomized_analysis
        opt_n.build()
      File "/opt/conda/lib/python3.5/site-packages/corda/corda.py", line 272, in build
        penalize_medium=False)
      File "/opt/conda/lib/python3.5/site-packages/corda/corda.py", line 236, in build
        need = self.associated(include)
      File "<ipython-input-23-9847a66c3d0b>", line 38, in multiple_randomized_analysis
        opt_n.build()
    KeyboardInterrupt
      File "/opt/conda/lib/python3.5/site-packages/optlang/glpk_interface.py", line 669, in _optimize
        status = self._run_glp_simplex()
      File "/opt/conda/lib/python3.5/multiprocessing/process.py", line 249, in _bootstrap
        self.run()
      File "/opt/conda/lib/python3.5/site-packages/corda/corda.py", line 272, in build
        penalize_medium=False)
      File "/opt/conda/lib/python3.5/site-packages/optlang/glpk_interface.py", line 643, in _run_glp_simplex
        return_value = glp_simplex(self.problem, self.configuration._smcp)
      File "/opt/conda/lib/python3.5/site-packages/corda/corda.py", line 202, in associated
        sol = self.model.solver.primal_values
      File "/opt/conda/lib/python3.5/site-packages/corda/corda.py", line 182, in associated
        self.__zero_objective()
      File "/opt/conda/lib/python3.5/site-packages/corda/corda.py", line 272, in build
        penalize_medium=False)
    KeyboardInterrupt
      File "/opt/conda/lib/python3.5/site-packages/corda/corda.py", line 196, in associated
        sol = self.model.solver.optimize()
      File "/opt/conda/lib/python3.5/multiprocessing/process.py", line 93, in run
        self._target(*self._args, **self._kwargs)
    KeyboardInterrupt
      File "/opt/conda/lib/python3.5/site-packages/corda/corda.py", line 138, in __zero_objective
        {v: 0 for v in self.model.variables})
      File "/opt/conda/lib/python3.5/multiprocessing/pool.py", line 119, in worker
        result = (True, func(*args, **kwds))
      File "/opt/conda/lib/python3.5/site-packages/optlang/interface.py", line 1234, in primal_values
        zip(self._get_variables_names(), self._get_primal_values())
      File "/opt/conda/lib/python3.5/site-packages/optlang/glpk_interface.py", line 352, in set_linear_coefficients
        glp_set_obj_coef(self.problem.problem, variable._index, float(coefficient))
      File "<ipython-input-23-9847a66c3d0b>", line 38, in multiple_randomized_analysis
        opt_n.build()
      File "/opt/conda/lib/python3.5/site-packages/optlang/interface.py", line 1222, in _get_variables_names
        return [variable.name for variable in self.variables]
      File "/opt/conda/lib/python3.5/site-packages/corda/corda.py", line 272, in build
        penalize_medium=False)
      File "/opt/conda/lib/python3.5/site-packages/optlang/interface.py", line 1222, in <listcomp>
        return [variable.name for variable in self.variables]
      File "/opt/conda/lib/python3.5/site-packages/optlang/container.py", line 76, in __iter__
        if original_length != len(self._object_list):
    Process ForkPoolWorker-15:
    KeyboardInterrupt
      File "/opt/conda/lib/python3.5/site-packages/corda/corda.py", line 202, in associated
        sol = self.model.solver.primal_values
      File "/opt/conda/lib/python3.5/site-packages/corda/corda.py", line 196, in associated
        sol = self.model.solver.optimize()
      File "/opt/conda/lib/python3.5/site-packages/optlang/interface.py", line 1468, in optimize
        status = self._optimize()
      File "/opt/conda/lib/python3.5/site-packages/optlang/interface.py", line 443, in problem
        return getattr(self, '_problem', None)
      File "/opt/conda/lib/python3.5/site-packages/optlang/interface.py", line 1234, in primal_values
        zip(self._get_variables_names(), self._get_primal_values())
    Traceback (most recent call last):
      File "/opt/conda/lib/python3.5/site-packages/optlang/interface.py", line 1468, in optimize
        status = self._optimize()
    KeyboardInterrupt
      File "/opt/conda/lib/python3.5/site-packages/optlang/glpk_interface.py", line 669, in _optimize
        status = self._run_glp_simplex()
      File "/opt/conda/lib/python3.5/multiprocessing/process.py", line 249, in _bootstrap
        self.run()
      File "/opt/conda/lib/python3.5/site-packages/optlang/glpk_interface.py", line 669, in _optimize
        status = self._run_glp_simplex()
      File "/opt/conda/lib/python3.5/multiprocessing/process.py", line 93, in run
        self._target(*self._args, **self._kwargs)
      File "/opt/conda/lib/python3.5/site-packages/optlang/glpk_interface.py", line 643, in _run_glp_simplex
        return_value = glp_simplex(self.problem, self.configuration._smcp)
      File "/opt/conda/lib/python3.5/site-packages/optlang/glpk_interface.py", line 643, in _run_glp_simplex
        return_value = glp_simplex(self.problem, self.configuration._smcp)
    KeyboardInterrupt
    KeyboardInterrupt
    Process ForkPoolWorker-17:
      File "/opt/conda/lib/python3.5/multiprocessing/pool.py", line 119, in worker
        result = (True, func(*args, **kwds))
      File "<ipython-input-23-9847a66c3d0b>", line 38, in multiple_randomized_analysis
        opt_n.build()
      File "/opt/conda/lib/python3.5/site-packages/corda/corda.py", line 272, in build
        penalize_medium=False)
      File "/opt/conda/lib/python3.5/site-packages/corda/corda.py", line 196, in associated
        sol = self.model.solver.optimize()
    Traceback (most recent call last):
      File "/opt/conda/lib/python3.5/site-packages/optlang/interface.py", line 1468, in optimize
        status = self._optimize()
      File "/opt/conda/lib/python3.5/multiprocessing/process.py", line 249, in _bootstrap
        self.run()
      File "/opt/conda/lib/python3.5/site-packages/optlang/glpk_interface.py", line 669, in _optimize
        status = self._run_glp_simplex()
    Traceback (most recent call last):
      File "/opt/conda/lib/python3.5/multiprocessing/process.py", line 93, in run
        self._target(*self._args, **self._kwargs)
      File "/opt/conda/lib/python3.5/site-packages/optlang/glpk_interface.py", line 643, in _run_glp_simplex
        return_value = glp_simplex(self.problem, self.configuration._smcp)
      File "/opt/conda/lib/python3.5/multiprocessing/process.py", line 249, in _bootstrap
        self.run()
      File "/opt/conda/lib/python3.5/multiprocessing/pool.py", line 119, in worker
        result = (True, func(*args, **kwds))
    KeyboardInterrupt
      File "/opt/conda/lib/python3.5/site-packages/optlang/interface.py", line 1222, in _get_variables_names
        return [variable.name for variable in self.variables]
      File "<ipython-input-23-9847a66c3d0b>", line 42, in multiple_randomized_analysis
        targets=get_drugable_targets(model_n, model_c, "biopsys" )
      File "/opt/conda/lib/python3.5/site-packages/optlang/interface.py", line 1222, in <listcomp>
        return [variable.name for variable in self.variables]
      File "<ipython-input-18-b299aa21be4a>", line 29, in get_drugable_targets
        nmodel.reactions.get_by_id(rx).bounds=(-eps,eps)
      File "/opt/conda/lib/python3.5/site-packages/optlang/interface.py", line 206, in name
        return self._name
      File "/opt/conda/lib/python3.5/site-packages/cobra/util/context.py", line 69, in wrapper
        f(self, new_value)
    KeyboardInterrupt
      File "/opt/conda/lib/python3.5/site-packages/cobra/core/reaction.py", line 265, in bounds
        update_forward_and_reverse_bounds(self)
      File "/opt/conda/lib/python3.5/site-packages/cobra/core/reaction.py", line 1083, in update_forward_and_reverse_bounds
        reaction.reverse_variable._lb = None
      File "/opt/conda/lib/python3.5/site-packages/cobra/core/reaction.py", line 162, in reverse_variable
        self.reverse_id]
      File "/opt/conda/lib/python3.5/site-packages/cobra/core/reaction.py", line 111, in reverse_id
        self.id.encode('utf-8')).hexdigest()[0:5]))
    KeyboardInterrupt
      File "/opt/conda/lib/python3.5/multiprocessing/process.py", line 93, in run
        self._target(*self._args, **self._kwargs)
    Process ForkPoolWorker-23:
      File "/opt/conda/lib/python3.5/multiprocessing/pool.py", line 119, in worker
        result = (True, func(*args, **kwds))
    Traceback (most recent call last):
      File "<ipython-input-23-9847a66c3d0b>", line 38, in multiple_randomized_analysis
        opt_n.build()
      File "/opt/conda/lib/python3.5/multiprocessing/process.py", line 249, in _bootstrap
        self.run()
      File "/opt/conda/lib/python3.5/site-packages/corda/corda.py", line 272, in build
        penalize_medium=False)
      File "/opt/conda/lib/python3.5/multiprocessing/process.py", line 93, in run
        self._target(*self._args, **self._kwargs)
      File "/opt/conda/lib/python3.5/site-packages/corda/corda.py", line 196, in associated
        sol = self.model.solver.optimize()
      File "/opt/conda/lib/python3.5/multiprocessing/pool.py", line 119, in worker
        result = (True, func(*args, **kwds))
      File "/opt/conda/lib/python3.5/site-packages/optlang/interface.py", line 1468, in optimize
        status = self._optimize()
      File "/opt/conda/lib/python3.5/site-packages/optlang/glpk_interface.py", line 669, in _optimize
        status = self._run_glp_simplex()
      File "<ipython-input-23-9847a66c3d0b>", line 38, in multiple_randomized_analysis
        opt_n.build()
      File "/opt/conda/lib/python3.5/site-packages/optlang/glpk_interface.py", line 643, in _run_glp_simplex
        return_value = glp_simplex(self.problem, self.configuration._smcp)
      File "/opt/conda/lib/python3.5/site-packages/corda/corda.py", line 272, in build
        penalize_medium=False)
    KeyboardInterrupt
      File "/opt/conda/lib/python3.5/site-packages/corda/corda.py", line 196, in associated
        sol = self.model.solver.optimize()
      File "/opt/conda/lib/python3.5/site-packages/optlang/interface.py", line 1468, in optimize
        status = self._optimize()
      File "/opt/conda/lib/python3.5/site-packages/optlang/glpk_interface.py", line 669, in _optimize
        status = self._run_glp_simplex()
      File "/opt/conda/lib/python3.5/site-packages/optlang/glpk_interface.py", line 643, in _run_glp_simplex
        return_value = glp_simplex(self.problem, self.configuration._smcp)
    KeyboardInterrupt



    ---------------------------------------------------------------------------

    KeyboardInterrupt                         Traceback (most recent call last)

    <ipython-input-23-9847a66c3d0b> in <module>()
         51 
         52 pool.close()
    ---> 53 pool.join()
    

    /opt/conda/lib/python3.5/multiprocessing/pool.py in join(self)
        508         util.debug('joining pool')
        509         assert self._state in (CLOSE, TERMINATE)
    --> 510         self._worker_handler.join()
        511         self._task_handler.join()
        512         self._result_handler.join()


    /opt/conda/lib/python3.5/threading.py in join(self, timeout)
       1052 
       1053         if timeout is None:
    -> 1054             self._wait_for_tstate_lock()
       1055         else:
       1056             # the behavior of a negative timeout isn't documented, but


    /opt/conda/lib/python3.5/threading.py in _wait_for_tstate_lock(self, block, timeout)
       1068         if lock is None:  # already determined that the C code is done
       1069             assert self._is_stopped
    -> 1070         elif lock.acquire(block, timeout):
       1071             lock.release()
       1072             self._stop()


    KeyboardInterrupt: 


    Common reactions size 1106
    Normal unique reactions size 1210
    Disease unique reactions size 1165


    cobra/util/solver.py:404 [1;31mUserWarning[0m: solver status is 'infeasible'


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
