
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
and the individual confidence values have to be mapped from the gene confidence levels. Here "and" is evaluated by the minimum confidence and "or" by the maximum confidence. The Python package includes a handy function to do this for you automatically in a safe manner. For that you will require the gene-reaction rule (Recon 1 and 2 include them in their model for instance) and a dictionary mapping genes/proteins to their confidence values. For examples:
In [1]:


### Load the ubiquitin scores of the microarrays

#### Definitions
   If we have a matrix M where rows are unique genes and columns are samples (cancer samples for example) we define g[i,j] as the gene i in the sample j where g[i,j] in {0,1}  where 0 means a gene that is "Off" in the sample and 1 represents "On".
   Whith this we calculate the percentil (0 to 1) of each row (genes) and this will be the ubiquitin score for gene i in the model, this score represents the number of times a gene is considered on or off in that specific data, and we do this for each model we whant to do. 


```python
import cobra
#!wget http://www.ebi.ac.uk/biomodels-main/download?mid=MODEL1603150001 -O recon2.2.xml
Recon2 = cobra.io.read_sbml_model("Models/recon2.2.xml")
#Recon204 = cobra.io.load_matlab_model("Models/Recon2.v04.mat")
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
```

    Reactions: 7785
    Metabolites: 5324



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




    <Solution 555.786 at 0x7f3de2fd3630>




```python
for rxex in Recon2.exchanges:
    rxex.bounds=(-0.01,10)

### Western Diet
#### Sugars
Recon2.reactions.EX_arab_L_LPAREN_e_RPAREN_.lower_bound=-0.17878295 
#Cellobiose not in Human
Recon2.reactions.EX_drib_LPAREN_e_RPAREN_.lower_bound=-0.17878295 
Recon2.reactions.EX_fru_LPAREN_e_RPAREN_.lower_bound=-0.14898579
Recon2.reactions.EX_fuc_L_LPAREN_e_RPAREN_.lower_bound=-0.14898579 
Recon2.reactions.EX_gal_LPAREN_e_RPAREN_.lower_bound=-0.14898579 
Recon2.reactions.EX_glc_LPAREN_e_RPAREN_.lower_bound=-0.14898579
#Falta Ex_glcn en Recon2.2
Recon2.reactions.EX_lcts_LPAREN_e_RPAREN_.lower_bound=-0.07449289 
Recon2.reactions.EX_malt_LPAREN_e_RPAREN_.lower_bound=-0.07449289 
Recon2.reactions.EX_man_LPAREN_e_RPAREN_.lower_bound=-0.14898579 
#Melibiose  not in Human
#D-Mannitol not in Human
Recon2.reactions.EX_oxa_LPAREN_e_RPAREN_.lower_bound=-0.44695737 
Recon2.reactions.EX_rib_D_LPAREN_e_RPAREN_.lower_bound=-0.17878295 
#L-Rhamnose not in Human
Recon2.reactions.EX_sucr_LPAREN_e_RPAREN_.lower_bound=-0.07449289 
Recon2.reactions.EX_tre_LPAREN_e_RPAREN_.lower_bound=-0.07449289 
Recon2.reactions.EX_xyl_D_LPAREN_e_RPAREN_.lower_bound=-0.17878295
Recon2.reactions.EX_strch1_LPAREN_e_RPAREN_.lower_bound=-0.25733909 

#fiber not in human

### Western Diet
#### Fat
Recon2.reactions.EX_arachd_LPAREN_e_RPAREN_.lower_bound=-0.00332813 
Recon2.reactions.EX_chsterol_LPAREN_e_RPAREN_.lower_bound=-0.00495795
Recon2.reactions.EX_glyc_LPAREN_e_RPAREN_.lower_bound=-1.79965486 
Recon2.reactions.EX_hdca_LPAREN_e_RPAREN_.lower_bound=-0.39637090
Recon2.reactions.EX_hdcea_LPAREN_e_RPAREN_.lower_bound=-0.03651697
Recon2.reactions.EX_lnlc_LPAREN_e_RPAREN_.lower_bound=-0.35910921
Recon2.reactions.EX_lnlnca_LPAREN_e_RPAREN_.lower_bound=-0.01756512
Recon2.reactions.EX_lnlncg_LPAREN_e_RPAREN_.lower_bound=-0.01756512
Recon2.reactions.EX_ocdca_LPAREN_e_RPAREN_.lower_bound=-0.16928260
Recon2.reactions.EX_ocdcea_LPAREN_e_RPAREN_.lower_bound=-0.68144465
Recon2.reactions.EX_octa_LPAREN_e_RPAREN_.lower_bound=-0.01294272
Recon2.reactions.EX_ttdca_LPAREN_e_RPAREN_.lower_bound=-0.06867567

### Western Diet
#### Protein
Recon2.reactions.EX_ala_L_LPAREN_e_RPAREN_.lower_bound=-1
Recon2.reactions.EX_arg_L_LPAREN_e_RPAREN_.lower_bound=-0.15
Recon2.reactions.EX_asn_L_LPAREN_e_RPAREN_.lower_bound=-0.225
Recon2.reactions.EX_asp_L_LPAREN_e_RPAREN_.lower_bound=-0.225
Recon2.reactions.EX_cys_L_LPAREN_e_RPAREN_.lower_bound=-1
Recon2.reactions.EX_gln_L_LPAREN_e_RPAREN_.lower_bound=-0.18
Recon2.reactions.EX_glu_L_LPAREN_e_RPAREN_.lower_bound=-0.18
Recon2.reactions.EX_gly_LPAREN_e_RPAREN_.lower_bound=-0.45
Recon2.reactions.EX_his_L_LPAREN_e_RPAREN_.lower_bound=-0.15
Recon2.reactions.EX_ile_L_LPAREN_e_RPAREN_.lower_bound=-0.15
Recon2.reactions.EX_leu_L_LPAREN_e_RPAREN_.lower_bound=-0.15
Recon2.reactions.EX_lys_L_LPAREN_e_RPAREN_.lower_bound=-0.15 
Recon2.reactions.EX_met_L_LPAREN_e_RPAREN_.lower_bound=-0.18
Recon2.reactions.EX_phe_L_LPAREN_e_RPAREN_.lower_bound=-1
Recon2.reactions.EX_pro_L_LPAREN_e_RPAREN_.lower_bound=-0.18
Recon2.reactions.EX_ser_L_LPAREN_e_RPAREN_.lower_bound=-1
Recon2.reactions.EX_thr_L_LPAREN_e_RPAREN_.lower_bound=-0.225
Recon2.reactions.EX_trp_L_LPAREN_e_RPAREN_.lower_bound=-0.08181818
Recon2.reactions.EX_tyr_L_LPAREN_e_RPAREN_.lower_bound=-1
Recon2.reactions.EX_val_L_LPAREN_e_RPAREN_.lower_bound=-0.18

### Western Diet
#### Minterals, Vitamins, others

#12dgr180
#26dap_M
#2dmmq8
#2obut
Recon2.reactions.EX_3mop_LPAREN_e_RPAREN_.lower_bound=-1
#4abz
#4hbz
Recon2.reactions.EX_ac_LPAREN_e_RPAREN_.lower_bound=-1
Recon2.reactions.EX_acgam_LPAREN_e_RPAREN_.lower_bound=-1
Recon2.reactions.EX_acmana_LPAREN_e_RPAREN_.lower_bound=-1
#acnam
Recon2.reactions.EX_ade_LPAREN_e_RPAREN_.lower_bound=-1
Recon2.reactions.EX_adn_LPAREN_e_RPAREN_.lower_bound=-1
#adocbl
Recon2.reactions.EX_ala_D_LPAREN_e_RPAREN_.lower_bound=-1
Recon2.reactions.EX_amp_LPAREN_e_RPAREN_.lower_bound=-1
#arab_D
Recon2.reactions.EX_btn_LPAREN_e_RPAREN_.lower_bound=-1
Recon2.reactions.EX_ca2_LPAREN_e_RPAREN_.lower_bound=-1
#cbl1
Recon2.reactions.EX_cgly_LPAREN_e_RPAREN_.lower_bound=-1
Recon2.reactions.EX_chol_LPAREN_e_RPAREN_.lower_bound=-1
#chor
Recon2.reactions.EX_cit_LPAREN_e_RPAREN_.lower_bound=-1
Recon2.reactions.EX_cl_LPAREN_e_RPAREN_.lower_bound=-1
#cobalt2
Recon2.reactions.EX_csn_LPAREN_e_RPAREN_.lower_bound=-1
#cu2
Recon2.reactions.EX_dad_2_LPAREN_e_RPAREN_.lower_bound=-1
Recon2.reactions.EX_dcyt_LPAREN_e_RPAREN_.lower_bound=-1
Recon2.reactions.EX_ddca_LPAREN_e_RPAREN_.lower_bound=-1
Recon2.reactions.EX_dgsn_LPAREN_e_RPAREN_.lower_bound=-1
Recon2.reactions.EX_fe2_LPAREN_e_RPAREN_.lower_bound=-1
Recon2.reactions.EX_fe3_LPAREN_e_RPAREN_.lower_bound=-1
#fe3dcit
Recon2.reactions.EX_fald_LPAREN_e_RPAREN_.lower_bound=-1
Recon2.reactions.EX_fol_LPAREN_e_RPAREN_.lower_bound=-1
Recon2.reactions.EX_for_LPAREN_e_RPAREN_.lower_bound=-1
Recon2.reactions.EX_fum_LPAREN_e_RPAREN_.lower_bound=-1
Recon2.reactions.EX_gam_LPAREN_e_RPAREN_.lower_bound=-1
Recon2.reactions.EX_glu_L_LPAREN_e_RPAREN_=-1
Recon2.reactions.EX_glyc3p_LPAREN_e_RPAREN_.lower_bound=-1
Recon2.reactions.EX_gthox_LPAREN_e_RPAREN_.lower_bound=-1
Recon2.reactions.EX_gthrd_LPAREN_e_RPAREN_.lower_bound=-1
Recon2.reactions.EX_gua_LPAREN_e_RPAREN_.lower_bound=-1
Recon2.reactions.EX_h_LPAREN_e_RPAREN_.lower_bound=-1
Recon2.reactions.EX_h2o2_LPAREN_e_RPAREN_.lower_bound=-10
#h2
#h2s
Recon2.reactions.EX_hxan_LPAREN_e_RPAREN_.lower_bound=-1
#indole
Recon2.reactions.EX_k_LPAREN_e_RPAREN_.lower_bound=-1
#lanost
Recon2.reactions.EX_meoh_LPAREN_e_RPAREN_.lower_bound=-10
#metsox_S_L
#mg2
#mn2
#mobd
#mqn8
Recon2.reactions.EX_na1_LPAREN_e_RPAREN_.lower_bound=-1
Recon2.reactions.EX_nac_LPAREN_e_RPAREN_.lower_bound=-1
Recon2.reactions.EX_ncam_LPAREN_e_RPAREN_.lower_bound=-1
#nmn nmn_e doesn't exists in Recon 2.2 
Recon2.reactions.EX_no2_LPAREN_e_RPAREN_.lower_bound=-1
Recon2.reactions.EX_orn_LPAREN_e_RPAREN_.lower_bound=-1
Recon2.reactions.EX_pheme_LPAREN_e_RPAREN_.lower_bound=-1

Recon2.reactions.EX_pi_LPAREN_e_RPAREN_.lower_bound=-1
#Pimelate  not in Human
Recon2.reactions.EX_pnto_R_LPAREN_e_RPAREN_.lower_bound=-1
Recon2.reactions.EX_ptrc_LPAREN_e_RPAREN_.lower_bound=-1
Recon2.reactions.EX_pydam_LPAREN_e_RPAREN_.lower_bound=-1
Recon2.reactions.EX_pydx_LPAREN_e_RPAREN_.lower_bound=-1
Recon2.reactions.EX_pydx5p_LPAREN_e_RPAREN_.lower_bound=-1
Recon2.reactions.EX_pydxn_LPAREN_e_RPAREN_.lower_bound=-1
#ubiquinone-8 not in Human
Recon2.reactions.EX_ribflv_LPAREN_e_RPAREN_.lower_bound=-1
Recon2.reactions.EX_sel_LPAREN_e_RPAREN_.lower_bound=-1
#Siroheme not in Human
Recon2.reactions.EX_so4_LPAREN_e_RPAREN_.lower_bound=-1
Recon2.reactions.EX_spmd_LPAREN_e_RPAREN_.lower_bound=-1
Recon2.reactions.EX_thm_LPAREN_e_RPAREN_.lower_bound=-1
Recon2.reactions.EX_thymd_LPAREN_e_RPAREN_.lower_bound=-1
Recon2.reactions.EX_ura_LPAREN_e_RPAREN_.lower_bound=-1
Recon2.reactions.EX_uri_LPAREN_e_RPAREN_.lower_bound=-1
#Xanthine not in Human
#zn2 not in Human

#Recon2.reactions.biomass_reaction.lower_bound=0.008
#Recon2.reactions.biomass_reaction.upper_bound=0.1
```

#extra #HMDB
Recon2.reactions.EX_xylt_LPAREN_e_RPAREN_.lower_bound=-0.1 #0.7
Recon2.reactions.r0408.upper_bound=0
Recon2.reactions.r0408.lower_bound=0
Recon2.reactions.EX_HC00250_LPAREN_e_RPAREN_.lower_bound=0
Recon2.reactions.EX_HC00250_LPAREN_e_RPAREN_.upper_bound=0



```python
Recon2.optimize()
```




    <Solution 0.270 at 0x7f3e151f6ac8>




```python
import pandas as pd
ub_scores_matrix =pd.read_csv("Binary_CellLines_RECON2_2.tsv",sep=",",index_col=0,header=0)
confidence_scores_matrix=ub_scores_matrix
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
      <th>HeLaHGU133A</th>
      <th>KeratinocytesHGU133A</th>
      <th>HeLaHGU133Plus2</th>
      <th>KeratinocytesHGU133Plus2</th>
      <th>MaxHeLa</th>
      <th>Maxkeratinocytes</th>
      <th>MeanHeLa</th>
      <th>Meankeratinocytes</th>
    </tr>
  </thead>
  <tbody>
    <tr>
      <th>HGNC:10006</th>
      <td>0.0</td>
      <td>0.0</td>
      <td>0.00</td>
      <td>0.0000</td>
      <td>0.00</td>
      <td>0.0000</td>
      <td>0.00</td>
      <td>0.00000</td>
    </tr>
    <tr>
      <th>HGNC:1027</th>
      <td>0.0</td>
      <td>0.0</td>
      <td>0.18</td>
      <td>0.0000</td>
      <td>0.18</td>
      <td>0.0000</td>
      <td>0.09</td>
      <td>0.00000</td>
    </tr>
    <tr>
      <th>HGNC:10293</th>
      <td>0.0</td>
      <td>0.0</td>
      <td>1.00</td>
      <td>0.8750</td>
      <td>1.00</td>
      <td>0.8750</td>
      <td>0.50</td>
      <td>0.43750</td>
    </tr>
    <tr>
      <th>HGNC:10297</th>
      <td>0.0</td>
      <td>0.0</td>
      <td>0.54</td>
      <td>0.5625</td>
      <td>0.54</td>
      <td>0.5625</td>
      <td>0.27</td>
      <td>0.28125</td>
    </tr>
    <tr>
      <th>HGNC:10451</th>
      <td>1.0</td>
      <td>1.0</td>
      <td>1.00</td>
      <td>1.0000</td>
      <td>1.00</td>
      <td>1.0000</td>
      <td>1.00</td>
      <td>1.00000</td>
    </tr>
  </tbody>
</table>
</div>




```python
for model in ub_scores_matrix.columns:
    #High Confidence rate 
       
    HC=ub_scores_matrix[model][ub_scores_matrix[model] >= 0.75].index.tolist()
    confidence_scores_matrix[model][HC]=3
        
    #Medium Confidence Rate 
    MC=ub_scores_matrix[model][(ub_scores_matrix[model] <  0.75) & (ub_scores_matrix[model] >= 0.5) ].index.tolist()
    confidence_scores_matrix[model][MC]=2
        
    #Low
    LC=ub_scores_matrix[model][(ub_scores_matrix[model] < 0.5) & (ub_scores_matrix[model] >=0.25) ].index.tolist()
    confidence_scores_matrix[model][LC]=1
    
    #unknown
    ZC=ub_scores_matrix[model][(ub_scores_matrix[model] <.25) & (ub_scores_matrix[model] >= 0.1) ].index.tolist()
    confidence_scores_matrix[model][ZC]=0
    
    #Negative
    NC=ub_scores_matrix[model][(ub_scores_matrix[model] < 0.1)].index.tolist()
    confidence_scores_matrix[model][NC]=-1
 
```


```python
list(confidence_scores_matrix)
```




    ['HeLaHGU133A',
     'KeratinocytesHGU133A',
     'HeLaHGU133Plus2',
     'KeratinocytesHGU133Plus2',
     'MaxHeLa',
     'Maxkeratinocytes',
     'MeanHeLa',
     'Meankeratinocytes']




```python
from corda import reaction_confidence



conf_HeLaHGU133A = {}
conf_HeLaPlus2 = {}
conf_keratinocytesHGU133A = {}
conf_keratinocytesPlus2 = {}

conf_HeLaMax = {}
conf_HeLaMean = {}
conf_keratinocytesMean= {}
conf_keratinocytesMax= {}


for r in Recon2.reactions:
        conf_HeLaHGU133A[r.id]=reaction_confidence(r.gene_reaction_rule,confidence_scores_matrix["HeLaHGU133A"])
        conf_HeLaPlus2[r.id]=reaction_confidence(r.gene_reaction_rule,confidence_scores_matrix["HeLaHGU133Plus2"])
        conf_HeLaMax[r.id]=reaction_confidence(r.gene_reaction_rule,confidence_scores_matrix["MaxHeLa"])
        conf_HeLaMean[r.id]=reaction_confidence(r.gene_reaction_rule,confidence_scores_matrix["MeanHeLa"])

        conf_keratinocytesPlus2[r.id]=reaction_confidence(r.gene_reaction_rule,confidence_scores_matrix["KeratinocytesHGU133Plus2"])
        conf_keratinocytesHGU133A[r.id]=reaction_confidence(r.gene_reaction_rule,confidence_scores_matrix["KeratinocytesHGU133A"])
        conf_keratinocytesMax[r.id]=reaction_confidence(r.gene_reaction_rule,confidence_scores_matrix["Maxkeratinocytes"])
        conf_keratinocytesMean[r.id]=reaction_confidence(r.gene_reaction_rule,confidence_scores_matrix["Meankeratinocytes"])

```


```python
conf_HeLaHGU133A["biomass_reaction"]=3
conf_HeLaHGU133A["biomass_DNA"]=3
conf_HeLaHGU133A["biomass_RNA"]=3
conf_HeLaHGU133A["biomass_carbohydrate"]=3
conf_HeLaHGU133A["biomass_lipid"]=3
conf_HeLaHGU133A["biomass_other"]=3
conf_HeLaHGU133A["biomass_protein"]=3

conf_HeLaPlus2["biomass_reaction"]=3
conf_HeLaPlus2["biomass_DNA"]=3
conf_HeLaPlus2["biomass_RNA"]=3
conf_HeLaPlus2["biomass_carbohydrate"]=3
conf_HeLaPlus2["biomass_lipid"]=3
conf_HeLaPlus2["biomass_other"]=3
conf_HeLaPlus2["biomass_protein"]=3


conf_keratinocytesHGU133A["biomass_reaction"]=3
conf_keratinocytesHGU133A["biomass_DNA"]=3
conf_keratinocytesHGU133A["biomass_RNA"]=3
conf_keratinocytesHGU133A["biomass_carbohydrate"]=3
conf_keratinocytesHGU133A["biomass_lipid"]=3
conf_keratinocytesHGU133A["biomass_other"]=3
conf_keratinocytesHGU133A["biomass_protein"]=3

conf_keratinocytesPlus2["biomass_reaction"]=3
conf_keratinocytesPlus2["biomass_DNA"]=3
conf_keratinocytesPlus2["biomass_RNA"]=3
conf_keratinocytesPlus2["biomass_carbohydrate"]=3
conf_keratinocytesPlus2["biomass_lipid"]=3
conf_keratinocytesPlus2["biomass_other"]=3
conf_keratinocytesPlus2["biomass_protein"]=3



Recon2.objective="biomass_reaction"
```


```python
metas = ['adp_c','atp_c', 'glc_D_c',  'fru_c','nad_c', 'nadh_c','nad_m', 'nadh_m','nadph_c', 'nadph_m', 'nadp_c', 'nadp_m', 'cmp_c', 'HC00342_c', 'glcn_c', 'citr_L_c', 'glyb_c', 'icit_c', '3pg_c', 'accoa_m ->coa_m', 'akg_m', 'e4p_c', 'f6p_c', 'g3p_c', 'g6p_c', 'oaa_m', 'pep_c', 'pyr_c', 'r5p_c', 'succoa_m ->coa_m', 'ala_L_c', 'arg_L_c', 'asn_L_c', 'asp_L_c', 'val_L_c', 'adp_c', 'thr_L_c', 'leu_L_c', 'gln_L_c', 'glu_L_c', 'gly_c', 'pro_D_c', 'pro_L_c', 'ser_L_c', 'ctp_c', 'fdp_c','utp_c', 'pmtcoa_c -> coa_c', 'chsterol_c', 'dag_hs_c', 'tag_hs_c', 'mag_hs_c', 'gthox_c','gthrd_c','ru5p_D_c','crm_hs_c', 'pa_hs_c', 'pe_hs_c', 'ps_hs_c', 'hxan_c', 'creat_c', 'crtn_c', 'chol_c', 'orn_c', '4hpro_LT_c', 's7p_c', 'amp_c', 'udp_c', 'ade_c', 'asp_D_c', 'adn_c', 'ump_c', 'his_L_c', 'tyr_L_c', 'phe_L_c', 'gdp_c', 'ctp_c', 'cys_L_c', 'amet_c', 'cit_c', 'tym_c', 'succ_c', 'CE1936_c', 'gua_c', 'cdp_c', 'g1p_c', 'lys_L_c', 'bhb_c', 'ile_L_c', 'gmp_c', 'dhap_c', 'fum_c', 'mal_L_c', 'spmd_c', 'ala_B_c', 'trp_L_c', 'lac_L_c', 'met_L_c', 'ptrc_c', '4abut_c']
```


```python
%%time

from corda import CORDA



opt_HeLaHGU133A = CORDA(model=Recon2, confidence=conf_HeLaHGU133A, n=5, met_prod=metas,  penalty_factor=1000) 
opt_HeLaHGU133A.build()
print(opt_HeLaHGU133A)
model_HeLaHGU133A=opt_HeLaHGU133A.cobra_model(name="HeLaHGU133A")
print(model_HeLaHGU133A.optimize())


```

    build status: reconstruction complete
    Inc. reactions: 1963/7885
     - unclear: 782/3337
     - exclude: 383/3311
     - low and medium: 55/423
     - high: 743/814
    
    <Solution 0.108 at 0x7f3de0ad0518>
    CPU times: user 10min 18s, sys: 444 ms, total: 10min 18s
    Wall time: 10min 18s



```python
%%time

from corda import CORDA


opt_HeLaPlus2 = CORDA(model=Recon2, confidence=conf_HeLaPlus2, n=5, met_prod=metas,  penalty_factor=1000) 
opt_HeLaPlus2.build()
print(opt_HeLaPlus2)
model_HeLaPlus2=opt_HeLaPlus2.cobra_model(name="HeLaPlus2")
print(model_HeLaPlus2.optimize())

```

    build status: reconstruction complete
    Inc. reactions: 3019/7885
     - unclear: 945/3258
     - exclude: 349/2010
     - low and medium: 113/818
     - high: 1612/1799
    
    <Solution 0.216 at 0x7f3de07585c0>
    CPU times: user 19min 51s, sys: 884 ms, total: 19min 52s
    Wall time: 19min 52s



```python
%%time

from corda import CORDA

opt_keratinocytesHGU133A = CORDA(model=Recon2, confidence=conf_keratinocytesHGU133A, n=5, met_prod=metas,  penalty_factor=1000) 
opt_keratinocytesHGU133A.build()
print(opt_keratinocytesHGU133A)
model_keratinocytesHGU133A=opt_keratinocytesHGU133A.cobra_model(name="keratinocytesHGU133A")

print(model_keratinocytesHGU133A.optimize())
```

    build status: reconstruction complete
    Inc. reactions: 1517/7885
     - unclear: 695/3317
     - exclude: 231/3619
     - low and medium: 90/391
     - high: 501/558
    
    <Solution 0.100 at 0x7f3d97699eb8>
    CPU times: user 7min 53s, sys: 432 ms, total: 7min 53s
    Wall time: 7min 53s



```python
%%time

from corda import CORDA

opt_keratinocytesPlus2 = CORDA(model=Recon2, confidence=conf_keratinocytesPlus2, n=5, met_prod=metas,  penalty_factor=1000) 
opt_keratinocytesPlus2.build()
print(opt_keratinocytesPlus2)
model_keratinocytesPlus2=opt_keratinocytesPlus2.cobra_model(name="keratinocytesPlus2")

print(model_keratinocytesPlus2.optimize())

```

    build status: reconstruction complete
    Inc. reactions: 3363/7885
     - unclear: 1014/3247
     - exclude: 307/1960
     - low and medium: 145/600
     - high: 1897/2078
    
    <Solution 0.253 at 0x7f3d97beddd8>
    CPU times: user 18min 28s, sys: 772 ms, total: 18min 29s
    Wall time: 18min 29s



```python

```

    <Solution 0.058 at 0x7f3de12f2e80>
    <Solution 0.176 at 0x7f3de12f2da0>
    <Solution 0.061 at 0x7f3de12f2c18>
    <Solution 0.253 at 0x7f3de12f2b00>



```python

```


```python
model_KeratinocytesMean.reactions.biomass_reaction.boundaryd=0.008
model_HeLaMean.reactions.biomass_reaction.upper_bound=0.043

print(model_HeLaMean.optimize())


model_KeratinocytesMean.reactions.biomass_reaction.lower_bound=0.008
model_KeratinocytesMean.reactions.biomass_reaction.upper_bound=0.037

print(model_KeratinocytesMean.optimize())

```

    <Solution 0.043 at 0x7f1738583390>
    <Solution 0.037 at 0x7f1738583160>



```python
cobra.io.write_sbml_model(model_HeLaMean, "Full_HeLa_model_1118_DEMEM6429_n3.sbml")
cobra.io.write_sbml_model(model_KeratinocytesMean, "Full_Kerat_model_1118_DEMEM6429_n3.sbml")
```


```python
%matplotlib inline  

from matplotlib_venn import venn2, venn3
 
# Make the diagram
venn2([set([ k.id for k in model_HeLaMean.reactions ]),set([ k.id for k in model_KeratinocytesMean.reactions ])],("HeLa","Keratinocytes"))

```




    <matplotlib_venn._common.VennDiagram at 0x7f1739b76630>




![png](output_24_1.png)



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
        
        results[rx]["norm_flux"]=nflx1
        results[rx]["dise_flux"]=dflx1
 
        results[rx]["norm_prolif_ratio"]=nflx1/nflx0
        results[rx]["dise_prolif_ratio"]=dflx1/dflx0
        
        results[rx]["norm_dise_ratio"]=(nflx1/nflx0)/(dflx1/dflx0)
        
        results[rx]["genes"]=normal_Model.reactions.get_by_id(rx).gene_name_reaction_rule

        
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
        
        results[rx]["norm_flux"]=nflx1
        results[rx]["dise_flux"]=dflx0
 
        results[rx]["norm_prolif_ratio"]=nflx1/nflx0
        results[rx]["dise_prolif_ratio"]=(nflx1/nflx0)
        
        results[rx]["norm_dise_ratio"]=nflx1/nflx0
        
        results[rx]["genes"]=normal_Model.reactions.get_by_id(rx).gene_name_reaction_rule

            
        nmodel.reactions.get_by_id(rx).bounds=nbounds
        
    for rx in unique_Drx:
        #print(rx)
        
        dbounds=dmodel.reactions.get_by_id(rx).bounds
    
        dmodel.reactions.get_by_id(rx).bounds=(-eps,eps)
        dfba=dmodel.optimize()
        dflx1=dfba.f
        
        results[rx]={}

        results[rx]["model"]=model_name

        results[rx]["norm_flux"]=nflx0
        results[rx]["dise_flux"]=dflx1
 
        results[rx]["norm_prolif_ratio"]=nflx0
        results[rx]["dise_prolif_ratio"]=dflx1/dflx0
        
        results[rx]["norm_dise_ratio"]=1/(dflx1/dflx0)
        
        results[rx]["genes"]=dmodel.reactions.get_by_id(rx).gene_name_reaction_rule

        dmodel.reactions.get_by_id(rx).bounds=dbounds


    return(pd.DataFrame(results).transpose())

```


```python
R_CLines=get_drugable_targets(
    normal_Model=model_KeratinocytesMean, 
    disease_Model=model_HeLaMean, 
    model_name="CellLines")

#Glyceraldehyde-3-phosphate dehydrogenase GAPD
R_CLines[ (R_CLines["norm_dise_ratio"] > 1.0000001) ].sort_values(by="norm_dise_ratio", ascending=False)
#R_Biopsy.sort_values(by="norm_dise_ratio", ascending=False)
```

    Common reactions size 1203
    Normal unique reactions size 212
    Disease unique reactions size 777





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
      <th>dise_flux</th>
      <th>dise_prolif_ratio</th>
      <th>genes</th>
      <th>model</th>
      <th>norm_dise_ratio</th>
      <th>norm_flux</th>
      <th>norm_prolif_ratio</th>
    </tr>
  </thead>
  <tbody>
    <tr>
      <th>biomass_reaction</th>
      <td>0.0581578</td>
      <td>1.35251</td>
      <td></td>
      <td>CellLines</td>
      <td>1.21002</td>
      <td>0.0605526</td>
      <td>1.63656</td>
    </tr>
  </tbody>
</table>
</div>




```python
%%time

from corda import CORDA


opt_MeanKeratinocytes = CORDA(model=Recon2, confidence=conf_Meankeratinocytes, n=2, met_prod=metas) 
opt_MeanKeratinocytes.build()
print(opt_MeanKeratinocytes)
model_MeanKeratinocytes=opt_MeanKeratinocytes.cobra_model(name="MeanKeratinocytes")
print(model_MeanKeratinocytes.optimize())

```

    build status: reconstruction complete
    Inc. reactions: 3155/7885
     - unclear: 1007/3256
     - exclude: 403/2659
     - low and medium: 0/0
     - high: 1745/1970
    
    <Solution 0.253 at 0x7fb139a68b38>
    CPU times: user 18min 33s, sys: 1.01 s, total: 18min 34s
    Wall time: 18min 34s



```python
print(model_MeanHeLa.reactions.count)

print(model_MeanHeLa.optimize(),"3629")
print(model_MeanKeratinocytes.optimize(), "3196")
```

    <built-in method count of DictList object at 0x7fb13f699548>
    <Solution 0.270 at 0x7fb14115b320> 3629
    <Solution 0.253 at 0x7fb14115be10> 3196



```python
R_Biopsy=get_drugable_targets(
    normal_Model=model_MeanKeratinocytes, 
    disease_Model=model_MeanHeLa, 
    model_name="CellLines")

#Glyceraldehyde-3-phosphate dehydrogenase GAPD
#R_Biopsy[ (R_Biopsy["norm_dise_ratio"] > 1.0000001) ].sort_values(by="norm_dise_ratio", ascending=False)
R_Biopsy.sort_values(by="norm_dise_ratio", ascending=False)
```

    Common reactions size 2737
    Normal unique reactions size 318
    Disease unique reactions size 783





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
      <th>dise_flux</th>
      <th>dise_prolif_ratio</th>
      <th>genes</th>
      <th>model</th>
      <th>norm_dise_ratio</th>
      <th>norm_flux</th>
      <th>norm_prolif_ratio</th>
    </tr>
  </thead>
  <tbody>
    <tr>
      <th>EX_leu_L_LPAREN_e_RPAREN_</th>
      <td>0.183305</td>
      <td>0.678353</td>
      <td></td>
      <td>CellLines</td>
      <td>1.17333</td>
      <td>0.201635</td>
      <td>0.795934</td>
    </tr>
    <tr>
      <th>EX_pro_L_LPAREN_e_RPAREN_</th>
      <td>0.242436</td>
      <td>0.89718</td>
      <td></td>
      <td>CellLines</td>
      <td>1.06667</td>
      <td>0.242436</td>
      <td>0.956992</td>
    </tr>
    <tr>
      <th>biomass_reaction</th>
      <td>0.1</td>
      <td>0.370069</td>
      <td></td>
      <td>CellLines</td>
      <td>1.06667</td>
      <td>0.1</td>
      <td>0.39474</td>
    </tr>
    <tr>
      <th>biomass_protein</th>
      <td>0.141643</td>
      <td>0.524177</td>
      <td></td>
      <td>CellLines</td>
      <td>1.06667</td>
      <td>0.141643</td>
      <td>0.559122</td>
    </tr>
    <tr>
      <th>EX_gam_LPAREN_e_RPAREN_</th>
      <td>0.27022</td>
      <td>1</td>
      <td></td>
      <td>CellLines</td>
      <td>1</td>
      <td>0.253331</td>
      <td>1</td>
    </tr>
    <tr>
      <th>EX_sbt_DASH_d_LPAREN_e_RPAREN_</th>
      <td>0.27022</td>
      <td>1</td>
      <td></td>
      <td>CellLines</td>
      <td>1</td>
      <td>0.253331</td>
      <td>1</td>
    </tr>
    <tr>
      <th>r0497</th>
      <td>0.27022</td>
      <td>1</td>
      <td>or  or</td>
      <td>CellLines</td>
      <td>1</td>
      <td>0.253331</td>
      <td>1</td>
    </tr>
    <tr>
      <th>r0782</th>
      <td>0.27022</td>
      <td>1</td>
      <td></td>
      <td>CellLines</td>
      <td>1</td>
      <td>0.253331</td>
      <td>1</td>
    </tr>
    <tr>
      <th>PROSTGE1t3</th>
      <td>0.27022</td>
      <td>1</td>
      <td></td>
      <td>CellLines</td>
      <td>1</td>
      <td>0.253331</td>
      <td>1</td>
    </tr>
    <tr>
      <th>SCP21cx</th>
      <td>0.27022</td>
      <td>1</td>
      <td></td>
      <td>CellLines</td>
      <td>1</td>
      <td>0.253331</td>
      <td>1</td>
    </tr>
    <tr>
      <th>r0653</th>
      <td>0.27022</td>
      <td>1</td>
      <td>( and ) or ( and )</td>
      <td>CellLines</td>
      <td>1</td>
      <td>0.253331</td>
      <td>1</td>
    </tr>
    <tr>
      <th>B3GALTg</th>
      <td>0.27022</td>
      <td>1</td>
      <td></td>
      <td>CellLines</td>
      <td>1</td>
      <td>0.253331</td>
      <td>1</td>
    </tr>
    <tr>
      <th>EX_malt_LPAREN_e_RPAREN_</th>
      <td>0.27022</td>
      <td>1</td>
      <td></td>
      <td>CellLines</td>
      <td>1</td>
      <td>0.253331</td>
      <td>1</td>
    </tr>
    <tr>
      <th>RE3250X</th>
      <td>0.27022</td>
      <td>1</td>
      <td>or</td>
      <td>CellLines</td>
      <td>1</td>
      <td>0.253331</td>
      <td>1</td>
    </tr>
    <tr>
      <th>S6TASE12ly</th>
      <td>0.27022</td>
      <td>1</td>
      <td></td>
      <td>CellLines</td>
      <td>1</td>
      <td>0.253331</td>
      <td>1</td>
    </tr>
    <tr>
      <th>CLPNDt</th>
      <td>0.27022</td>
      <td>1</td>
      <td></td>
      <td>CellLines</td>
      <td>1</td>
      <td>0.253331</td>
      <td>1</td>
    </tr>
    <tr>
      <th>SMVGLUChep</th>
      <td>0.27022</td>
      <td>1</td>
      <td>or</td>
      <td>CellLines</td>
      <td>1</td>
      <td>0.253331</td>
      <td>1</td>
    </tr>
    <tr>
      <th>SIAASE4ly</th>
      <td>0.27022</td>
      <td>1</td>
      <td>and  and  and</td>
      <td>CellLines</td>
      <td>1</td>
      <td>0.253331</td>
      <td>1</td>
    </tr>
    <tr>
      <th>G14T5g</th>
      <td>0.27022</td>
      <td>1</td>
      <td>or  or  or</td>
      <td>CellLines</td>
      <td>1</td>
      <td>0.253331</td>
      <td>1</td>
    </tr>
    <tr>
      <th>EX_pglyc_hs_LPAREN_e_RPAREN_</th>
      <td>0.27022</td>
      <td>1</td>
      <td></td>
      <td>CellLines</td>
      <td>1</td>
      <td>0.253331</td>
      <td>1</td>
    </tr>
    <tr>
      <th>DSAT</th>
      <td>0.27022</td>
      <td>1</td>
      <td></td>
      <td>CellLines</td>
      <td>1</td>
      <td>0.253331</td>
      <td>1</td>
    </tr>
    <tr>
      <th>GUAPRT</th>
      <td>0.27022</td>
      <td>1</td>
      <td></td>
      <td>CellLines</td>
      <td>1</td>
      <td>0.253331</td>
      <td>1</td>
    </tr>
    <tr>
      <th>SCP21x</th>
      <td>0.27022</td>
      <td>1</td>
      <td></td>
      <td>CellLines</td>
      <td>1</td>
      <td>0.253331</td>
      <td>1</td>
    </tr>
    <tr>
      <th>GULNter</th>
      <td>0.27022</td>
      <td>1</td>
      <td></td>
      <td>CellLines</td>
      <td>1</td>
      <td>0.253331</td>
      <td>1</td>
    </tr>
    <tr>
      <th>r1817</th>
      <td>0.27022</td>
      <td>1</td>
      <td></td>
      <td>CellLines</td>
      <td>1</td>
      <td>0.253331</td>
      <td>0.253331</td>
    </tr>
    <tr>
      <th>TDCHOLAte</th>
      <td>0.27022</td>
      <td>1</td>
      <td></td>
      <td>CellLines</td>
      <td>1</td>
      <td>0.253331</td>
      <td>0.253331</td>
    </tr>
    <tr>
      <th>r0907</th>
      <td>0.27022</td>
      <td>1</td>
      <td></td>
      <td>CellLines</td>
      <td>1</td>
      <td>0.253331</td>
      <td>0.253331</td>
    </tr>
    <tr>
      <th>r1873</th>
      <td>0.27022</td>
      <td>1</td>
      <td></td>
      <td>CellLines</td>
      <td>1</td>
      <td>0.253331</td>
      <td>0.253331</td>
    </tr>
    <tr>
      <th>r1954</th>
      <td>0.27022</td>
      <td>1</td>
      <td></td>
      <td>CellLines</td>
      <td>1</td>
      <td>0.253331</td>
      <td>0.253331</td>
    </tr>
    <tr>
      <th>UREAt</th>
      <td>0.27022</td>
      <td>1</td>
      <td>or  or  or</td>
      <td>CellLines</td>
      <td>1</td>
      <td>0.253331</td>
      <td>0.253331</td>
    </tr>
    <tr>
      <th>...</th>
      <td>...</td>
      <td>...</td>
      <td>...</td>
      <td>...</td>
      <td>...</td>
      <td>...</td>
      <td>...</td>
    </tr>
    <tr>
      <th>r0783</th>
      <td>0.27022</td>
      <td>1</td>
      <td></td>
      <td>CellLines</td>
      <td>1</td>
      <td>0.253331</td>
      <td>1</td>
    </tr>
    <tr>
      <th>HDCEAtr</th>
      <td>0.27022</td>
      <td>1</td>
      <td></td>
      <td>CellLines</td>
      <td>1</td>
      <td>0.253331</td>
      <td>1</td>
    </tr>
    <tr>
      <th>r0779</th>
      <td>0.27022</td>
      <td>1</td>
      <td></td>
      <td>CellLines</td>
      <td>1</td>
      <td>0.253331</td>
      <td>1</td>
    </tr>
    <tr>
      <th>RE3383M</th>
      <td>0.27022</td>
      <td>1</td>
      <td>or  or</td>
      <td>CellLines</td>
      <td>1</td>
      <td>0.253331</td>
      <td>1</td>
    </tr>
    <tr>
      <th>GALASE16ly</th>
      <td>0.27022</td>
      <td>1</td>
      <td>( and  and  and ) or ( and )</td>
      <td>CellLines</td>
      <td>1</td>
      <td>0.253331</td>
      <td>1</td>
    </tr>
    <tr>
      <th>r0580</th>
      <td>0.27022</td>
      <td>1</td>
      <td></td>
      <td>CellLines</td>
      <td>1</td>
      <td>0.253331</td>
      <td>1</td>
    </tr>
    <tr>
      <th>r1547</th>
      <td>0.27022</td>
      <td>1</td>
      <td></td>
      <td>CellLines</td>
      <td>1</td>
      <td>0.253331</td>
      <td>1</td>
    </tr>
    <tr>
      <th>EX_h_LPAREN_e_RPAREN_</th>
      <td>0.27022</td>
      <td>1</td>
      <td></td>
      <td>CellLines</td>
      <td>1</td>
      <td>0.253331</td>
      <td>1</td>
    </tr>
    <tr>
      <th>RE3040R</th>
      <td>0.27022</td>
      <td>1</td>
      <td></td>
      <td>CellLines</td>
      <td>1</td>
      <td>0.253331</td>
      <td>1</td>
    </tr>
    <tr>
      <th>NACHEXA11ly</th>
      <td>0.27022</td>
      <td>1</td>
      <td>( and ) or</td>
      <td>CellLines</td>
      <td>1</td>
      <td>0.253331</td>
      <td>1</td>
    </tr>
    <tr>
      <th>HSD3B12r</th>
      <td>0.27022</td>
      <td>1</td>
      <td></td>
      <td>CellLines</td>
      <td>1</td>
      <td>0.253331</td>
      <td>1</td>
    </tr>
    <tr>
      <th>RE3564M</th>
      <td>0.27022</td>
      <td>1</td>
      <td>or  or  or</td>
      <td>CellLines</td>
      <td>1</td>
      <td>0.253331</td>
      <td>1</td>
    </tr>
    <tr>
      <th>r2538</th>
      <td>0.27022</td>
      <td>1</td>
      <td></td>
      <td>CellLines</td>
      <td>1</td>
      <td>0.253331</td>
      <td>1</td>
    </tr>
    <tr>
      <th>EX_q10h2_LPAREN_e_RPAREN_</th>
      <td>0.27022</td>
      <td>1</td>
      <td></td>
      <td>CellLines</td>
      <td>1</td>
      <td>0.253331</td>
      <td>1</td>
    </tr>
    <tr>
      <th>HMBS</th>
      <td>0.27022</td>
      <td>1</td>
      <td></td>
      <td>CellLines</td>
      <td>1</td>
      <td>0.253331</td>
      <td>1</td>
    </tr>
    <tr>
      <th>SACCD4m</th>
      <td>0.27022</td>
      <td>1</td>
      <td></td>
      <td>CellLines</td>
      <td>1</td>
      <td>0.253331</td>
      <td>0.253331</td>
    </tr>
    <tr>
      <th>UNK3</th>
      <td>0.27022</td>
      <td>1</td>
      <td></td>
      <td>CellLines</td>
      <td>1</td>
      <td>0.253331</td>
      <td>0.253331</td>
    </tr>
    <tr>
      <th>r1788</th>
      <td>0.27022</td>
      <td>1</td>
      <td></td>
      <td>CellLines</td>
      <td>1</td>
      <td>0.253331</td>
      <td>0.253331</td>
    </tr>
    <tr>
      <th>r1679</th>
      <td>0.27022</td>
      <td>1</td>
      <td></td>
      <td>CellLines</td>
      <td>1</td>
      <td>0.253331</td>
      <td>0.253331</td>
    </tr>
    <tr>
      <th>FAOXCPRIST3x</th>
      <td>0.27022</td>
      <td>1</td>
      <td>( and  and  and ( or )) or ( and  and ( or ) a...</td>
      <td>CellLines</td>
      <td>1</td>
      <td>0.253331</td>
      <td>0.253331</td>
    </tr>
    <tr>
      <th>EX_glyc3p_LPAREN_e_RPAREN_</th>
      <td>0.27022</td>
      <td>1</td>
      <td></td>
      <td>CellLines</td>
      <td>1</td>
      <td>0.253331</td>
      <td>1</td>
    </tr>
    <tr>
      <th>EX_biocyt_LPAREN_e_RPAREN_</th>
      <td>0.274957</td>
      <td>1.01753</td>
      <td></td>
      <td>CellLines</td>
      <td>0.982772</td>
      <td>0.253331</td>
      <td>0.253331</td>
    </tr>
    <tr>
      <th>EX_lys_L_LPAREN_e_RPAREN_</th>
      <td>0.185776</td>
      <td>0.6875</td>
      <td></td>
      <td>CellLines</td>
      <td>0.969697</td>
      <td>0.168888</td>
      <td>0.666667</td>
    </tr>
    <tr>
      <th>EX_h2o2_LPAREN_e_RPAREN_</th>
      <td>0.27022</td>
      <td>1</td>
      <td></td>
      <td>CellLines</td>
      <td>0.962822</td>
      <td>0.243913</td>
      <td>0.962822</td>
    </tr>
    <tr>
      <th>H2O2t</th>
      <td>0.27022</td>
      <td>1</td>
      <td></td>
      <td>CellLines</td>
      <td>0.962822</td>
      <td>0.243913</td>
      <td>0.962822</td>
    </tr>
    <tr>
      <th>ENO</th>
      <td>0.248565</td>
      <td>0.919862</td>
      <td>or  or  or</td>
      <td>CellLines</td>
      <td>0.906837</td>
      <td>0.21132</td>
      <td>0.834165</td>
    </tr>
    <tr>
      <th>PGM</th>
      <td>0.248565</td>
      <td>0.919862</td>
      <td>or  or</td>
      <td>CellLines</td>
      <td>0.906837</td>
      <td>0.21132</td>
      <td>0.834165</td>
    </tr>
    <tr>
      <th>PGK</th>
      <td>0.247178</td>
      <td>0.914729</td>
      <td>or</td>
      <td>CellLines</td>
      <td>0.892236</td>
      <td>0.206757</td>
      <td>0.816154</td>
    </tr>
    <tr>
      <th>GAPD</th>
      <td>0.247178</td>
      <td>0.914729</td>
      <td>or</td>
      <td>CellLines</td>
      <td>0.812597</td>
      <td>0.188303</td>
      <td>0.743306</td>
    </tr>
    <tr>
      <th>LYStiDF</th>
      <td>0.27022</td>
      <td>1</td>
      <td>or  or  or</td>
      <td>CellLines</td>
      <td>0.666667</td>
      <td>0.168888</td>
      <td>0.666667</td>
    </tr>
  </tbody>
</table>
<p>3838 rows × 7 columns</p>
</div>




```python
R_Biopsy[ (R_Biopsy["norm_dise_ratio"] > 1.0000001) ].sort_values(by="norm_dise_ratio", ascending=False)

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
      <th>dise_flux</th>
      <th>dise_prolif_ratio</th>
      <th>genes</th>
      <th>model</th>
      <th>norm_dise_ratio</th>
      <th>norm_flux</th>
      <th>norm_prolif_ratio</th>
    </tr>
  </thead>
  <tbody>
    <tr>
      <th>EX_leu_L_LPAREN_e_RPAREN_</th>
      <td>0.183305</td>
      <td>0.678353</td>
      <td></td>
      <td>CellLines</td>
      <td>1.17333</td>
      <td>0.201635</td>
      <td>0.795934</td>
    </tr>
    <tr>
      <th>EX_pro_L_LPAREN_e_RPAREN_</th>
      <td>0.242436</td>
      <td>0.89718</td>
      <td></td>
      <td>CellLines</td>
      <td>1.06667</td>
      <td>0.242436</td>
      <td>0.956992</td>
    </tr>
    <tr>
      <th>biomass_protein</th>
      <td>0.141643</td>
      <td>0.524177</td>
      <td></td>
      <td>CellLines</td>
      <td>1.06667</td>
      <td>0.141643</td>
      <td>0.559122</td>
    </tr>
    <tr>
      <th>biomass_reaction</th>
      <td>0.1</td>
      <td>0.370069</td>
      <td></td>
      <td>CellLines</td>
      <td>1.06667</td>
      <td>0.1</td>
      <td>0.39474</td>
    </tr>
  </tbody>
</table>
</div>




```python
import pandas
from time import time

import cobra.test

from cobra.flux_analysis import (
    single_gene_deletion, single_reaction_deletion, double_gene_deletion,
    double_reaction_deletion)

deletion_results=single_reaction_deletion(model_MeanHeLa, model_MeanHeLa.reactions)

deletion_results
```


```python
import pandas
from time import time

import cobra.test

from cobra.flux_analysis import (
    single_gene_deletion, single_reaction_deletion, double_gene_deletion,
    double_reaction_deletion)

deletion_results=single_reaction_deletion(model_MeanKeratinocytes, model_MeanKeratinocytes.reactions)

deletion_results[deletion_results.status != 'optimal']
```


```python
closedModel=model_MaxCancerBiopsy.copy()
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
# glucose aerobic

closedModel.reactions.EX_o2_LPAREN_e_RPAREN_.bounds=(-1000,0.0)
closedModel.reactions.EX_h2o_LPAREN_e_RPAREN_.bounds=(-1000,1000)
closedModel.reactions.EX_co2_LPAREN_e_RPAREN_.bounds=(-0.0,1000)
closedModel.reactions.EX_glc_LPAREN_e_RPAREN_.bounds=(-1,-1)
closedModel.optimize()
```


```python
closedModel=model_MaxNormalBiopsy.copy()
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
# glucose aerobic

closedModel.reactions.EX_o2_LPAREN_e_RPAREN_.bounds=(-1000,0.0)
closedModel.reactions.EX_h2o_LPAREN_e_RPAREN_.bounds=(-1000,1000)
closedModel.reactions.EX_co2_LPAREN_e_RPAREN_.bounds=(-0.0,1000)
closedModel.reactions.EX_glc_LPAREN_e_RPAREN_.bounds=(-1,-1)
closedModel.optimize()
```


```python
closedModel=model_MaxCancerBiopsy.copy()
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
closedModel=model_MaxNormalBiopsy.copy()
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
closedModel=model_MaxCancerBiopsy.copy()
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
closedModel=model_MaxNormalBiopsy.copy()
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
closedModel=model_MaxCancerBiopsy.copy()
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

# glutamine anaerobic
closedModel.reactions.EX_o2_LPAREN_e_RPAREN_.bounds=(0,0)
closedModel.reactions.EX_h2o_LPAREN_e_RPAREN_.bounds=(-1000,1000)
closedModel.reactions.EX_co2_LPAREN_e_RPAREN_.bounds=(0,1000)
closedModel.reactions.EX_gln_L_LPAREN_e_RPAREN_.bounds=(-1,-1)

closedModel.optimize()
```


```python
closedModel=model_MaxNormalBiopsy.copy()
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

# glutamine anaerobic
closedModel.reactions.EX_o2_LPAREN_e_RPAREN_.bounds=(0,0)
closedModel.reactions.EX_h2o_LPAREN_e_RPAREN_.bounds=(-1000,1000)
closedModel.reactions.EX_co2_LPAREN_e_RPAREN_.bounds=(0,1000)
closedModel.reactions.EX_gln_L_LPAREN_e_RPAREN_.bounds=(-1,-1)

closedModel.optimize()
```


```python

```


```python

```

    Common reactions size 2762
    Normal unique reactions size 334
    Disease unique reactions size 767





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
      <th>dise_flux</th>
      <th>dise_prolif_ratio</th>
      <th>genes</th>
      <th>model</th>
      <th>norm_dise_ratio</th>
      <th>norm_flux</th>
      <th>norm_prolif_ratio</th>
    </tr>
  </thead>
  <tbody>
    <tr>
      <th>EX_leu_L_LPAREN_e_RPAREN_</th>
      <td>0.183305</td>
      <td>0.678353</td>
      <td></td>
      <td>CellLines</td>
      <td>1.17333</td>
      <td>0.201635</td>
      <td>0.795934</td>
    </tr>
    <tr>
      <th>EX_pro_L_LPAREN_e_RPAREN_</th>
      <td>0.242436</td>
      <td>0.89718</td>
      <td></td>
      <td>CellLines</td>
      <td>1.1146</td>
      <td>0.253331</td>
      <td>1</td>
    </tr>
    <tr>
      <th>biomass_protein</th>
      <td>0.141643</td>
      <td>0.524177</td>
      <td></td>
      <td>CellLines</td>
      <td>1.06667</td>
      <td>0.141643</td>
      <td>0.559122</td>
    </tr>
    <tr>
      <th>biomass_reaction</th>
      <td>0.1</td>
      <td>0.370069</td>
      <td></td>
      <td>CellLines</td>
      <td>1.06667</td>
      <td>0.1</td>
      <td>0.39474</td>
    </tr>
  </tbody>
</table>
</div>




```python

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
      <th>dise_flux</th>
      <th>dise_prolif_ratio</th>
      <th>genes</th>
      <th>model</th>
      <th>norm_dise_ratio</th>
      <th>norm_flux</th>
      <th>norm_prolif_ratio</th>
    </tr>
  </thead>
  <tbody>
    <tr>
      <th>EX_leu_L_LPAREN_e_RPAREN_</th>
      <td>0.183305</td>
      <td>0.678353</td>
      <td></td>
      <td>CellLines</td>
      <td>1.17333</td>
      <td>0.201635</td>
      <td>0.795934</td>
    </tr>
    <tr>
      <th>EX_pro_L_LPAREN_e_RPAREN_</th>
      <td>0.242436</td>
      <td>0.89718</td>
      <td></td>
      <td>CellLines</td>
      <td>1.1146</td>
      <td>0.253331</td>
      <td>1</td>
    </tr>
    <tr>
      <th>biomass_protein</th>
      <td>0.141643</td>
      <td>0.524177</td>
      <td></td>
      <td>CellLines</td>
      <td>1.06667</td>
      <td>0.141643</td>
      <td>0.559122</td>
    </tr>
    <tr>
      <th>biomass_reaction</th>
      <td>0.1</td>
      <td>0.370069</td>
      <td></td>
      <td>CellLines</td>
      <td>1.06667</td>
      <td>0.1</td>
      <td>0.39474</td>
    </tr>
    <tr>
      <th>GLUPROASCT1</th>
      <td>0.27022</td>
      <td>1</td>
      <td></td>
      <td>CellLines</td>
      <td>1</td>
      <td>0.253331</td>
      <td>0.253331</td>
    </tr>
    <tr>
      <th>PUNP3</th>
      <td>0.27022</td>
      <td>1</td>
      <td></td>
      <td>CellLines</td>
      <td>1</td>
      <td>0.253331</td>
      <td>1</td>
    </tr>
    <tr>
      <th>r1400</th>
      <td>0.27022</td>
      <td>1</td>
      <td></td>
      <td>CellLines</td>
      <td>1</td>
      <td>0.253331</td>
      <td>1</td>
    </tr>
    <tr>
      <th>RDH3a</th>
      <td>0.27022</td>
      <td>1</td>
      <td>or  or  or  or</td>
      <td>CellLines</td>
      <td>1</td>
      <td>0.253331</td>
      <td>1</td>
    </tr>
    <tr>
      <th>CO2tp</th>
      <td>0.27022</td>
      <td>1</td>
      <td></td>
      <td>CellLines</td>
      <td>1</td>
      <td>0.253331</td>
      <td>1</td>
    </tr>
    <tr>
      <th>NACHEXA12ly</th>
      <td>0.27022</td>
      <td>1</td>
      <td>( and ) or</td>
      <td>CellLines</td>
      <td>1</td>
      <td>0.253331</td>
      <td>1</td>
    </tr>
    <tr>
      <th>GTHOm</th>
      <td>0.27022</td>
      <td>1</td>
      <td></td>
      <td>CellLines</td>
      <td>1</td>
      <td>0.253331</td>
      <td>1</td>
    </tr>
    <tr>
      <th>DDCRNe</th>
      <td>0.27022</td>
      <td>1</td>
      <td></td>
      <td>CellLines</td>
      <td>1</td>
      <td>0.253331</td>
      <td>1</td>
    </tr>
    <tr>
      <th>CMPACNAtg</th>
      <td>0.27022</td>
      <td>1</td>
      <td></td>
      <td>CellLines</td>
      <td>1</td>
      <td>0.253331</td>
      <td>1</td>
    </tr>
    <tr>
      <th>OMPDC</th>
      <td>0.27022</td>
      <td>1</td>
      <td></td>
      <td>CellLines</td>
      <td>1</td>
      <td>0.253331</td>
      <td>1</td>
    </tr>
    <tr>
      <th>NDPK6</th>
      <td>0.27022</td>
      <td>1</td>
      <td>or  or ( and ) or</td>
      <td>CellLines</td>
      <td>1</td>
      <td>0.253331</td>
      <td>1</td>
    </tr>
    <tr>
      <th>r1135</th>
      <td>0.27022</td>
      <td>1</td>
      <td></td>
      <td>CellLines</td>
      <td>1</td>
      <td>0.253331</td>
      <td>1</td>
    </tr>
    <tr>
      <th>RE1521X</th>
      <td>0.27022</td>
      <td>1</td>
      <td></td>
      <td>CellLines</td>
      <td>1</td>
      <td>0.253331</td>
      <td>1</td>
    </tr>
    <tr>
      <th>r0682</th>
      <td>0.27022</td>
      <td>1</td>
      <td></td>
      <td>CellLines</td>
      <td>1</td>
      <td>0.253331</td>
      <td>1</td>
    </tr>
    <tr>
      <th>MM7Cbg</th>
      <td>0.27022</td>
      <td>1</td>
      <td>or  or</td>
      <td>CellLines</td>
      <td>1</td>
      <td>0.253331</td>
      <td>1</td>
    </tr>
    <tr>
      <th>EBP2r</th>
      <td>0.27022</td>
      <td>1</td>
      <td></td>
      <td>CellLines</td>
      <td>1</td>
      <td>0.253331</td>
      <td>1</td>
    </tr>
    <tr>
      <th>HCO3_NAt</th>
      <td>0.27022</td>
      <td>1</td>
      <td></td>
      <td>CellLines</td>
      <td>1</td>
      <td>0.253331</td>
      <td>0.253331</td>
    </tr>
    <tr>
      <th>NTD7e</th>
      <td>0.27022</td>
      <td>0.27022</td>
      <td>or</td>
      <td>CellLines</td>
      <td>1</td>
      <td>0.253331</td>
      <td>1</td>
    </tr>
    <tr>
      <th>FADH2ETC</th>
      <td>0.27022</td>
      <td>1</td>
      <td>( and  and ) or ( and )</td>
      <td>CellLines</td>
      <td>1</td>
      <td>0.253331</td>
      <td>1</td>
    </tr>
    <tr>
      <th>r1584</th>
      <td>0.27022</td>
      <td>1</td>
      <td></td>
      <td>CellLines</td>
      <td>1</td>
      <td>0.253331</td>
      <td>1</td>
    </tr>
    <tr>
      <th>IVCOAACBP</th>
      <td>0.27022</td>
      <td>1</td>
      <td></td>
      <td>CellLines</td>
      <td>1</td>
      <td>0.253331</td>
      <td>1</td>
    </tr>
    <tr>
      <th>r2511</th>
      <td>0.27022</td>
      <td>1</td>
      <td></td>
      <td>CellLines</td>
      <td>1</td>
      <td>0.253331</td>
      <td>1</td>
    </tr>
    <tr>
      <th>RE3231C</th>
      <td>0.27022</td>
      <td>1</td>
      <td></td>
      <td>CellLines</td>
      <td>1</td>
      <td>0.253331</td>
      <td>1</td>
    </tr>
    <tr>
      <th>MMSAD3m</th>
      <td>0.27022</td>
      <td>1</td>
      <td></td>
      <td>CellLines</td>
      <td>1</td>
      <td>0.253331</td>
      <td>1</td>
    </tr>
    <tr>
      <th>PCHOLP_hs</th>
      <td>0.27022</td>
      <td>1</td>
      <td>or</td>
      <td>CellLines</td>
      <td>1</td>
      <td>0.253331</td>
      <td>1</td>
    </tr>
    <tr>
      <th>EX_chol_LPAREN_e_RPAREN_</th>
      <td>0.27022</td>
      <td>1</td>
      <td></td>
      <td>CellLines</td>
      <td>1</td>
      <td>0.253331</td>
      <td>1</td>
    </tr>
    <tr>
      <th>...</th>
      <td>...</td>
      <td>...</td>
      <td>...</td>
      <td>...</td>
      <td>...</td>
      <td>...</td>
      <td>...</td>
    </tr>
    <tr>
      <th>H6ET3er</th>
      <td>0.27022</td>
      <td>1</td>
      <td>and</td>
      <td>CellLines</td>
      <td>1</td>
      <td>0.253331</td>
      <td>1</td>
    </tr>
    <tr>
      <th>r2108</th>
      <td>0.27022</td>
      <td>0.27022</td>
      <td></td>
      <td>CellLines</td>
      <td>1</td>
      <td>0.253331</td>
      <td>1</td>
    </tr>
    <tr>
      <th>EX_nrpphr_LPAREN_e_RPAREN_</th>
      <td>0.27022</td>
      <td>1</td>
      <td></td>
      <td>CellLines</td>
      <td>1</td>
      <td>0.253331</td>
      <td>1</td>
    </tr>
    <tr>
      <th>RE1100G</th>
      <td>0.27022</td>
      <td>1</td>
      <td></td>
      <td>CellLines</td>
      <td>1</td>
      <td>0.253331</td>
      <td>0.253331</td>
    </tr>
    <tr>
      <th>NACHORCTL3le</th>
      <td>0.27022</td>
      <td>1</td>
      <td></td>
      <td>CellLines</td>
      <td>1</td>
      <td>0.253331</td>
      <td>0.253331</td>
    </tr>
    <tr>
      <th>B3GALT3g</th>
      <td>0.27022</td>
      <td>1</td>
      <td></td>
      <td>CellLines</td>
      <td>1</td>
      <td>0.253331</td>
      <td>0.253331</td>
    </tr>
    <tr>
      <th>r1692</th>
      <td>0.27022</td>
      <td>1</td>
      <td></td>
      <td>CellLines</td>
      <td>1</td>
      <td>0.253331</td>
      <td>0.253331</td>
    </tr>
    <tr>
      <th>TRPO2</th>
      <td>0.27022</td>
      <td>1</td>
      <td>or  or</td>
      <td>CellLines</td>
      <td>1</td>
      <td>0.253331</td>
      <td>0.253331</td>
    </tr>
    <tr>
      <th>r1823</th>
      <td>0.27022</td>
      <td>1</td>
      <td></td>
      <td>CellLines</td>
      <td>1</td>
      <td>0.253331</td>
      <td>0.253331</td>
    </tr>
    <tr>
      <th>r1898</th>
      <td>0.27022</td>
      <td>1</td>
      <td></td>
      <td>CellLines</td>
      <td>1</td>
      <td>0.253331</td>
      <td>0.253331</td>
    </tr>
    <tr>
      <th>r0193</th>
      <td>0.27022</td>
      <td>1</td>
      <td></td>
      <td>CellLines</td>
      <td>1</td>
      <td>0.253331</td>
      <td>0.253331</td>
    </tr>
    <tr>
      <th>r1682</th>
      <td>0.27022</td>
      <td>1</td>
      <td></td>
      <td>CellLines</td>
      <td>1</td>
      <td>0.253331</td>
      <td>0.253331</td>
    </tr>
    <tr>
      <th>HSD3B2r</th>
      <td>0.27022</td>
      <td>1</td>
      <td></td>
      <td>CellLines</td>
      <td>1</td>
      <td>0.253331</td>
      <td>0.253331</td>
    </tr>
    <tr>
      <th>r1949</th>
      <td>0.27022</td>
      <td>1</td>
      <td></td>
      <td>CellLines</td>
      <td>1</td>
      <td>0.253331</td>
      <td>0.253331</td>
    </tr>
    <tr>
      <th>r1952</th>
      <td>0.27022</td>
      <td>1</td>
      <td></td>
      <td>CellLines</td>
      <td>1</td>
      <td>0.253331</td>
      <td>0.253331</td>
    </tr>
    <tr>
      <th>FAOXC15NADPx</th>
      <td>0.27022</td>
      <td>1</td>
      <td>and</td>
      <td>CellLines</td>
      <td>1</td>
      <td>0.253331</td>
      <td>0.253331</td>
    </tr>
    <tr>
      <th>r1668</th>
      <td>0.27022</td>
      <td>1</td>
      <td></td>
      <td>CellLines</td>
      <td>1</td>
      <td>0.253331</td>
      <td>0.253331</td>
    </tr>
    <tr>
      <th>r1990</th>
      <td>0.27022</td>
      <td>1</td>
      <td></td>
      <td>CellLines</td>
      <td>1</td>
      <td>0.253331</td>
      <td>0.253331</td>
    </tr>
    <tr>
      <th>FPGS8</th>
      <td>0.27022</td>
      <td>1</td>
      <td></td>
      <td>CellLines</td>
      <td>1</td>
      <td>0.253331</td>
      <td>1</td>
    </tr>
    <tr>
      <th>CYTK5n</th>
      <td>0.27022</td>
      <td>1</td>
      <td></td>
      <td>CellLines</td>
      <td>1</td>
      <td>0.253331</td>
      <td>1</td>
    </tr>
    <tr>
      <th>B3GNT12g</th>
      <td>0.27022</td>
      <td>1</td>
      <td></td>
      <td>CellLines</td>
      <td>1</td>
      <td>0.253331</td>
      <td>1</td>
    </tr>
    <tr>
      <th>EX_biocyt_LPAREN_e_RPAREN_</th>
      <td>0.274957</td>
      <td>1.01753</td>
      <td></td>
      <td>CellLines</td>
      <td>0.982772</td>
      <td>0.253331</td>
      <td>0.253331</td>
    </tr>
    <tr>
      <th>EX_lys_L_LPAREN_e_RPAREN_</th>
      <td>0.185776</td>
      <td>0.6875</td>
      <td></td>
      <td>CellLines</td>
      <td>0.969697</td>
      <td>0.168888</td>
      <td>0.666667</td>
    </tr>
    <tr>
      <th>EX_h2o2_LPAREN_e_RPAREN_</th>
      <td>0.27022</td>
      <td>1</td>
      <td></td>
      <td>CellLines</td>
      <td>0.961206</td>
      <td>0.243504</td>
      <td>0.961206</td>
    </tr>
    <tr>
      <th>H2O2t</th>
      <td>0.27022</td>
      <td>1</td>
      <td></td>
      <td>CellLines</td>
      <td>0.961206</td>
      <td>0.243504</td>
      <td>0.961206</td>
    </tr>
    <tr>
      <th>PGM</th>
      <td>0.248517</td>
      <td>0.919685</td>
      <td>or  or</td>
      <td>CellLines</td>
      <td>0.907012</td>
      <td>0.21132</td>
      <td>0.834165</td>
    </tr>
    <tr>
      <th>ENO</th>
      <td>0.248517</td>
      <td>0.919685</td>
      <td>or  or  or</td>
      <td>CellLines</td>
      <td>0.907012</td>
      <td>0.21132</td>
      <td>0.834165</td>
    </tr>
    <tr>
      <th>PGK</th>
      <td>0.24713</td>
      <td>0.91455</td>
      <td>or</td>
      <td>CellLines</td>
      <td>0.89241</td>
      <td>0.206757</td>
      <td>0.816154</td>
    </tr>
    <tr>
      <th>GAPD</th>
      <td>0.24713</td>
      <td>0.91455</td>
      <td>or</td>
      <td>CellLines</td>
      <td>0.812756</td>
      <td>0.188303</td>
      <td>0.743306</td>
    </tr>
    <tr>
      <th>LYStiDF</th>
      <td>0.27022</td>
      <td>1</td>
      <td>or  or  or</td>
      <td>CellLines</td>
      <td>0.666667</td>
      <td>0.168888</td>
      <td>0.666667</td>
    </tr>
  </tbody>
</table>
<p>3863 rows × 7 columns</p>
</div>




```python
cbmodel=model_MeanNormalBiopsy.copy()
pbounds=cbmodel.reactions.SPHMDAc.bounds
cbmodel.reactions.PCHOLHSTDe.bounds=[-0.1,0.1]
print(model_MeanNormalBiopsy.optimize())
print(cbmodel.optimize())

```


```python
crx=model_MaxNormalBiopsy.reactions
df=pandas.DataFrame(columns=['Name','ID','Reaction','Genes','RxRule'])
for idex in R_Biopsy.index.tolist():
    rx=crx.get_by_id(idex)
    df2=pandas.DataFrame(
        {
         'Name': [rx.name],
         'ID': [idex],
         'Reaction': [rx.reaction],
         'Genes' : [rx.genes],
         'RxRule' : [rx.gene_reaction_rule],
        }, columns=['Name','ID','Reaction','Genes','RxRule'])
    df=df.append(df2)
    
df
```


```python
import cobra.test
from cobra.flux_analysis import production_envelope

prod_env = production_envelope(model_MaxNormalBiopsy, ["EX_glc_LPAREN_e_RPAREN_", "EX_o2_LPAREN_e_RPAREN_"])
```


```python
prod_env
```


```python
model_MaxNormalBiopsy.reactions.EX_ac_LPAREN_e_RPAREN_.reaction
```


```python
model_MaxNormalBiopsy.objective.name
```


```python
%matplotlib inline

prod_env = production_envelope(
    model_MaxNormalBiopsy, ["EX_co2_LPAREN_e_RPAREN_","EX_glc_LPAREN_e_RPAREN_","EX_gly_LPAREN_e_RPAREN_",])

```


```python
prod_env
```


```python
prod_env.plot(
  kind='line', x='EX_glc_LPAREN_e_RPAREN_', y='EX_glc_LPAREN_e_RPAREN_')
```
