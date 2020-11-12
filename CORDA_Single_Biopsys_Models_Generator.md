
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
Recon1 = cobra.io.read_sbml_model("Models/RECON1.xml")
```




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

    <Solution 555.786 at 0x7f9132f6bf98>

```python
for rxex in Recon2.exchanges:
    rxex.bounds=(-0.01,10)

### Western Diet
#### Sugars
Recon2.reactions.get_by_id("EX_arab_L_LPAREN_e_RPAREN_").lower_bound=-0.17878295 
Recon2.reactions.get_by_id("EX_drib_LPAREN_e_RPAREN_").lower_bound=-0.17878295 
Recon2.reactions.get_by_id("EX_fru_LPAREN_e_RPAREN_").lower_bound=-0.14898579
Recon2.reactions.get_by_id("EX_fuc_L_LPAREN_e_RPAREN_").lower_bound=-0.14898579 
Recon2.reactions.get_by_id("EX_gal_LPAREN_e_RPAREN_").lower_bound=-0.14898579 
Recon2.reactions.get_by_id("EX_glc_LPAREN_e_RPAREN_").lower_bound=-0.14898579
Recon2.reactions.get_by_id("EX_lcts_LPAREN_e_RPAREN_").lower_bound=-0.07449289 
Recon2.reactions.get_by_id("EX_malt_LPAREN_e_RPAREN_").lower_bound=-0.07449289 
Recon2.reactions.get_by_id("EX_man_LPAREN_e_RPAREN_").lower_bound=-0.14898579 
Recon2.reactions.get_by_id("EX_oxa_LPAREN_e_RPAREN_").lower_bound=-0.44695737 
Recon2.reactions.get_by_id("EX_rib_D_LPAREN_e_RPAREN_").lower_bound=-0.17878295 
Recon2.reactions.get_by_id("EX_sucr_LPAREN_e_RPAREN_").lower_bound=-0.07449289 
Recon2.reactions.get_by_id("EX_tre_LPAREN_e_RPAREN_").lower_bound=-0.07449289 
Recon2.reactions.get_by_id("EX_xyl_D_LPAREN_e_RPAREN_").lower_bound=-0.17878295
Recon2.reactions.get_by_id("EX_strch1_LPAREN_e_RPAREN_").lower_bound=-0.25733909 


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
Recon2.reactions.get_by_id("EX_octa_LPAREN_e_RPAREN_").lower_bound=-0.01294272
Recon2.reactions.get_by_id("EX_ttdca_LPAREN_e_RPAREN_").lower_bound=-0.06867567

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

Recon2.reactions.get_by_id("EX_3mop_LPAREN_e_RPAREN_").lower_bound=-1
Recon2.reactions.get_by_id("EX_ac_LPAREN_e_RPAREN_").lower_bound=-1
Recon2.reactions.get_by_id("EX_acgam_LPAREN_e_RPAREN_").lower_bound=-1
Recon2.reactions.get_by_id("EX_acmana_LPAREN_e_RPAREN_").lower_bound=-1
Recon2.reactions.get_by_id("EX_ade_LPAREN_e_RPAREN_").lower_bound=-1
Recon2.reactions.get_by_id("EX_adn_LPAREN_e_RPAREN_").lower_bound=-1
Recon2.reactions.get_by_id("EX_ala_D_LPAREN_e_RPAREN_").lower_bound=-1
Recon2.reactions.get_by_id("EX_amp_LPAREN_e_RPAREN_").lower_bound=-1
Recon2.reactions.get_by_id("EX_btn_LPAREN_e_RPAREN_").lower_bound=-1
Recon2.reactions.get_by_id("EX_ca2_LPAREN_e_RPAREN_").lower_bound=-1
Recon2.reactions.get_by_id("EX_cgly_LPAREN_e_RPAREN_").lower_bound=-1
Recon2.reactions.get_by_id("EX_chol_LPAREN_e_RPAREN_").lower_bound=-1
Recon2.reactions.get_by_id("EX_cit_LPAREN_e_RPAREN_").lower_bound=-1
Recon2.reactions.get_by_id("EX_cl_LPAREN_e_RPAREN_").lower_bound=-1
Recon2.reactions.get_by_id("EX_csn_LPAREN_e_RPAREN_").lower_bound=-1
Recon2.reactions.get_by_id("EX_dad_2_LPAREN_e_RPAREN_").lower_bound=-1
Recon2.reactions.get_by_id("EX_dcyt_LPAREN_e_RPAREN_").lower_bound=-1
Recon2.reactions.get_by_id("EX_ddca_LPAREN_e_RPAREN_").lower_bound=-1
Recon2.reactions.get_by_id("EX_dgsn_LPAREN_e_RPAREN_").lower_bound=-1
Recon2.reactions.get_by_id("EX_fe2_LPAREN_e_RPAREN_").lower_bound=-1
Recon2.reactions.get_by_id("EX_fe3_LPAREN_e_RPAREN_").lower_bound=-1
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
Recon2.reactions.get_by_id("EX_h2o2_LPAREN_e_RPAREN_").lower_bound=-10
Recon2.reactions.get_by_id("EX_hxan_LPAREN_e_RPAREN_").lower_bound=-1
Recon2.reactions.get_by_id("EX_k_LPAREN_e_RPAREN_").lower_bound=-1
Recon2.reactions.get_by_id("EX_meoh_LPAREN_e_RPAREN_").lower_bound=-10
Recon2.reactions.get_by_id("EX_na1_LPAREN_e_RPAREN_").lower_bound=-1
Recon2.reactions.get_by_id("EX_nac_LPAREN_e_RPAREN_").lower_bound=-1
Recon2.reactions.get_by_id("EX_ncam_LPAREN_e_RPAREN_").lower_bound=-1
Recon2.reactions.get_by_id("EX_no2_LPAREN_e_RPAREN_").lower_bound=-1
Recon2.reactions.get_by_id("EX_orn_LPAREN_e_RPAREN_").lower_bound=-1
Recon2.reactions.get_by_id("EX_pheme_LPAREN_e_RPAREN_").lower_bound=-1

Recon2.reactions.get_by_id("EX_pi_LPAREN_e_RPAREN_").lower_bound=-1
Recon2.reactions.get_by_id("EX_pnto_R_LPAREN_e_RPAREN_").lower_bound=-1
Recon2.reactions.get_by_id("EX_ptrc_LPAREN_e_RPAREN_").lower_bound=-1
Recon2.reactions.get_by_id("EX_pydam_LPAREN_e_RPAREN_").lower_bound=-1
Recon2.reactions.get_by_id("EX_pydx_LPAREN_e_RPAREN_").lower_bound=-1
Recon2.reactions.get_by_id("EX_pydx5p_LPAREN_e_RPAREN_").lower_bound=-1
Recon2.reactions.get_by_id("EX_pydxn_LPAREN_e_RPAREN_").lower_bound=-1
Recon2.reactions.get_by_id("EX_ribflv_LPAREN_e_RPAREN_").lower_bound=-1
Recon2.reactions.get_by_id("EX_sel_LPAREN_e_RPAREN_").lower_bound=-1
Recon2.reactions.get_by_id("EX_so4_LPAREN_e_RPAREN_").lower_bound=-1
Recon2.reactions.get_by_id("EX_spmd_LPAREN_e_RPAREN_").lower_bound=-1
Recon2.reactions.get_by_id("EX_thm_LPAREN_e_RPAREN_").lower_bound=-1
Recon2.reactions.get_by_id("EX_thymd_LPAREN_e_RPAREN_").lower_bound=-1
Recon2.reactions.get_by_id("EX_ura_LPAREN_e_RPAREN_").lower_bound=-1
Recon2.reactions.get_by_id("EX_uri_LPAREN_e_RPAREN_").lower_bound=-1

#Recon2.reactions.biomass_reaction.lower_bound=0.008
#Recon2.reactions.biomass_reaction.upper_bound=0.1


```

New optimal value with the medium bounds.

```python
Recon2.optimize()
```

 <Solution 0.270 at 0x7f91331aae10>


```python
import pandas as pd
ub_scores_matrix =pd.read_csv("Binary_Byposys_RECON2_2.tsv",sep=",",index_col=0,header=0)
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
from corda import reaction_confidence



conf_CancerBiopsyHGU133A = {}
conf_NormalBiopsyHGU133A = {}
conf_CancerBiopsyHGU133Plus2 = {}
conf_NormalBiopsyHGU133APlus2 = {}


for r in Recon2.reactions:
        conf_CancerBiopsyHGU133A[r.id]=reaction_confidence(r.gene_reaction_rule,confidence_scores_matrix["CancerBiopsyHGU133A"])
        conf_NormalBiopsyHGU133A[r.id]=reaction_confidence(r.gene_reaction_rule,confidence_scores_matrix["NormalBiopsyHGU133A"])
        conf_CancerBiopsyHGU133Plus2[r.id]=reaction_confidence(r.gene_reaction_rule,confidence_scores_matrix["CancerBiopsyHGU133Plus2"])
        conf_NormalBiopsyHGU133APlus2[r.id]=reaction_confidence(r.gene_reaction_rule,confidence_scores_matrix["NormalBiopsyHGU133APlus2"])

```


```python
conf_CancerBiopsyHGU133A["biomass_reaction"]=3
conf_CancerBiopsyHGU133A["biomass_DNA"]=3
conf_CancerBiopsyHGU133A["biomass_RNA"]=3
conf_CancerBiopsyHGU133A["biomass_carbohydrate"]=3
conf_CancerBiopsyHGU133A["biomass_lipid"]=3
conf_CancerBiopsyHGU133A["biomass_other"]=3
conf_CancerBiopsyHGU133A["biomass_protein"]=3

conf_NormalBiopsyHGU133A["biomass_reaction"]=3
conf_NormalBiopsyHGU133A["biomass_DNA"]=3
conf_NormalBiopsyHGU133A["biomass_RNA"]=3
conf_NormalBiopsyHGU133A["biomass_carbohydrate"]=3
conf_NormalBiopsyHGU133A["biomass_lipid"]=3
conf_NormalBiopsyHGU133A["biomass_other"]=3
conf_NormalBiopsyHGU133A["biomass_protein"]=3


conf_CancerBiopsyHGU133Plus2["biomass_reaction"]=3
conf_CancerBiopsyHGU133Plus2["biomass_DNA"]=3
conf_CancerBiopsyHGU133Plus2["biomass_RNA"]=3
conf_CancerBiopsyHGU133Plus2["biomass_carbohydrate"]=3
conf_CancerBiopsyHGU133Plus2["biomass_lipid"]=3
conf_CancerBiopsyHGU133Plus2["biomass_other"]=3
conf_CancerBiopsyHGU133Plus2["biomass_protein"]=3

conf_NormalBiopsyHGU133APlus2["biomass_reaction"]=3
conf_NormalBiopsyHGU133APlus2["biomass_DNA"]=3
conf_NormalBiopsyHGU133APlus2["biomass_RNA"]=3
conf_NormalBiopsyHGU133APlus2["biomass_carbohydrate"]=3
conf_NormalBiopsyHGU133APlus2["biomass_lipid"]=3
conf_NormalBiopsyHGU133APlus2["biomass_other"]=3
conf_NormalBiopsyHGU133APlus2["biomass_protein"]=3


Recon2.objective="biomass_reaction"
```


```python
metas = ['adp_c', 'atp_c','glc_D_c',  'fru_c', 'nad_c', 'nadh_c','nad_m', 'nadh_m','nadph_c', 'nadph_m', 'nadp_c', 'nadp_m', 'cmp_c', 'HC00342_c', 'glcn_c', 'citr_L_c', 'glyb_c', 'icit_c', '3pg_c', 'accoa_m ->coa_m', 'akg_m', 'e4p_c', 'f6p_c', 'g3p_c', 'g6p_c', 'oaa_m', 'pep_c', 'pyr_c', 'r5p_c', 'succoa_m ->coa_m', 'ala_L_c', 'arg_L_c', 'asn_L_c', 'asp_L_c', 'val_L_c', 'adp_c', 'thr_L_c', 'leu_L_c', 'gln_L_c', 'glu_L_c', 'gly_c', 'pro_D_c', 'pro_L_c', 'ser_L_c', 'ctp_c', 'fdp_c','utp_c', 'pmtcoa_c -> coa_c', 'chsterol_c', 'dag_hs_c', 'tag_hs_c', 'mag_hs_c', 'gthox_c','gthrd_c','ru5p_D_c','crm_hs_c', 'pa_hs_c', 'pe_hs_c', 'ps_hs_c', 'hxan_c', 'creat_c', 'crtn_c', 'chol_c', 'orn_c', '4hpro_LT_c', 's7p_c', 'amp_c', 'udp_c', 'ade_c', 'asp_D_c', 'adn_c', 'ump_c', 'his_L_c', 'tyr_L_c', 'phe_L_c', 'gdp_c', 'ctp_c', 'cys_L_c', 'amet_c', 'cit_c', 'tym_c', 'succ_c', 'CE1936_c', 'gua_c', 'cdp_c', 'g1p_c', 'lys_L_c', 'bhb_c', 'ile_L_c', 'gmp_c', 'dhap_c', 'fum_c', 'mal_L_c', 'spmd_c', 'ala_B_c', 'trp_L_c', 'lac_L_c', 'met_L_c', 'ptrc_c', '4abut_c','0.014 biomass_DNA_c + 0.058 biomass_RNA_c + 0.071 biomass_carbohydrate_c + 0.097 biomass_lipid_c + 0.054 biomass_other_c + 0.706 biomass_protein_c --> ']
```


```python
%%time

from corda import CORDA

opt_NormalBiopsyHGU133APlus2 = CORDA(model=Recon2, confidence=conf_NormalBiopsyHGU133APlus2, n=5, met_prod=metas, penalty_factor=500 ) 
opt_NormalBiopsyHGU133APlus2.build()
print(opt_NormalBiopsyHGU133APlus2)
model_NormalBiopsyHGU133APlus2=opt_NormalBiopsyHGU133APlus2.cobra_model(name="NormalBiopsyHGU133APlus2")
print(model_NormalBiopsyHGU133APlus2.optimize())
```

    build status: reconstruction complete
    Inc. reactions: 2508/7886
     - unclear: 993/3312
     - exclude: 310/2794
     - low and medium: 189/638
     - high: 1016/1142
    
    <Solution 0.129 at 0x7f911de87390>
    CPU times: user 13min 37s, sys: 1.25 s, total: 13min 39s
    Wall time: 13min 39s



```python
%%time

from corda import CORDA


opt_CancerBiopsyHGU133Plus2 = CORDA(model=Recon2, confidence=conf_CancerBiopsyHGU133Plus2, n=5, met_prod=metas, penalty_factor=500) 
opt_CancerBiopsyHGU133Plus2.build()
print(opt_CancerBiopsyHGU133Plus2)
model_CancerBiopsyHGU133Plus2=opt_CancerBiopsyHGU133Plus2.cobra_model(name="CancerBiopsyHGU133Plus2")
print(model_CancerBiopsyHGU133Plus2.optimize())
```

    build status: reconstruction complete
    Inc. reactions: 2742/7886
     - unclear: 945/3247
     - exclude: 279/2447
     - low and medium: 142/645
     - high: 1376/1547
    
    <Solution 0.132 at 0x7f911d5ecbe0>
    CPU times: user 15min 52s, sys: 1.7 s, total: 15min 54s
    Wall time: 15min 54s



```python
%%time

from corda import CORDA

opt_CancerBiopsyHGU133A = CORDA(model=Recon2, confidence=conf_CancerBiopsyHGU133A, n=5, met_prod=metas,  penalty_factor=500) 
opt_CancerBiopsyHGU133A.build()
print(opt_CancerBiopsyHGU133A)
model_CancerBiopsyHGU133A=opt_CancerBiopsyHGU133A.cobra_model(name="CancerBiopsyHGU133A")
print(model_CancerBiopsyHGU133A.optimize())
```

    build status: reconstruction complete
    Inc. reactions: 1338/7886
     - unclear: 662/3316
     - exclude: 194/3454
     - low and medium: 88/662
     - high: 394/454
    
    <Solution 0.154 at 0x7f9130dca5f8>
    CPU times: user 6min 55s, sys: 552 ms, total: 6min 56s
    Wall time: 6min 56s



```python
%%time

from corda import CORDA

opt_NormalBiopsyHGU133A = CORDA(model=Recon2, confidence=conf_NormalBiopsyHGU133A, n=5, met_prod=metas,  penalty_factor=500 ) 
opt_NormalBiopsyHGU133A.build()
print(opt_NormalBiopsyHGU133APlus2)
model_NormalBiopsyHGU133A=opt_NormalBiopsyHGU133A.cobra_model(name="NormalBiopsyHGU133A")
print(model_NormalBiopsyHGU133A.optimize())
```

    build status: reconstruction complete
    Inc. reactions: 2407/2407
     - unclear: 993/993
     - exclude: 310/310
     - low and medium: 189/189
     - high: 915/915
    
    <Solution 0.088 at 0x7f911d19fb70>
    CPU times: user 7min 37s, sys: 928 ms, total: 7min 38s
    Wall time: 7min 38s



```python
%matplotlib inline  

from matplotlib_venn import venn2, venn3
 
# Make the diagram
venn2([set([ k.id for k in model_CancerBiopsyHGU133APlus2.reactions ]),set([ k.id for k in model_CancerBiopsyHGU133A.reactions ])],("Plus2","133A"))




```


    ---------------------------------------------------------------------------

    NameError                                 Traceback (most recent call last)

    <ipython-input-72-77d4292e842b> in <module>()
          4 
          5 # Make the diagram
    ----> 6 venn2([set([ k.id for k in model_CancerBiopsyHGU133APlus2.reactions ]),set([ k.id for k in model_CancerBiopsyHGU133A.reactions ])],("Plus2","133A"))
          7 
          8 


    NameError: name 'model_CancerBiopsyHGU133APlus2' is not defined



```python
%matplotlib inline  

from matplotlib_venn import venn2, venn3
 
# Make the diagram
venn2([set([ k.id for k in model_NormalBiopsyHGU133APlus2.reactions ]),set([ k.id for k in model_NormalBiopsyHGU133A.reactions ])],("Plus2","133A"))




```




    <matplotlib_venn._common.VennDiagram at 0x7f911ca65828>




![png](output_19_1.png)



```python
cobra.io.write_sbml_model(model_MeanCancerBiopsy, "Full_SCC_model_1118_western_diet_n6.sbml")
cobra.io.write_sbml_model(model_MeanNormalBiopsy, "Full_Normal_model_1118_western_diet_n6.sbml")
```


    ---------------------------------------------------------------------------

    NameError                                 Traceback (most recent call last)

    <ipython-input-31-19ae9cc94a82> in <module>()
    ----> 1 cobra.io.write_sbml_model(model_MeanCancerBiopsy, "Full_SCC_model_1118_western_diet_n6.sbml")
          2 cobra.io.write_sbml_model(model_MeanNormalBiopsy, "Full_Normal_model_1118_western_diet_n6.sbml")


    NameError: name 'model_MeanCancerBiopsy' is not defined



```python
model_MeanCancerBiopsy.reactions.biomass_reaction.lower_bound=0.008
model_MeanNormalBiopsy.reactions.biomass_reaction.lower_bound=0.008
model_MeanCancerBiopsy.reactions.biomass_reaction.upper_bound=1
model_MeanNormalBiopsy.reactions.biomass_reaction.upper_bound=1


print(model_MeanNormalBiopsy.optimize(), 1705)
print(model_MeanCancerBiopsy.optimize(), 1813)

```

    <Solution 0.050 at 0x7f5c6df69c88> 1705
    <Solution 0.056 at 0x7f5c6df69240> 1813



```python
%matplotlib inline  

from matplotlib_venn import venn2, venn3
 
# Make the diagram
venn2([set([ k.id for k in model_MeanCancerBiopsy.reactions ]),set([ k.id for k in model_MeanNormalBiopsy.reactions ])],("Cancer","Normal"))

```




    <matplotlib_venn._common.VennDiagram at 0x7f5c6decf7f0>




![png](output_22_1.png)



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
        results[rx]["dise_prolif_ratio"]=dflx0
        
        results[rx]["norm_dise_ratio"]=nflx1/nflx0
            
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
        
        dmodel.reactions.get_by_id(rx).bounds=dbounds


    return(pd.DataFrame(results).transpose())

```


```python
Targets_biopsys=get_drugable_targets(model_MeanNormalBiopsy, model_MeanCancerBiopsy, "biopsys" )
#Targets_biopsys.to_csv("Drugable_targets_FULL_biopsys.tsv",sep="\t")

```


    ---------------------------------------------------------------------------

    NameError                                 Traceback (most recent call last)

    <ipython-input-76-b6fe64716c84> in <module>()
    ----> 1 Targets_biopsys=get_drugable_targets(model_MeanNormalBiopsy, model_MeanCancerBiopsy, "biopsys" )
          2 #Targets_biopsys.to_csv("Drugable_targets_FULL_biopsys.tsv",sep="\t")


    NameError: name 'model_MeanNormalBiopsy' is not defined



```python
Targets_biopsys[ (Targets_biopsys["norm_dise_ratio"] > 1.0000001) ].sort_values(by="norm_dise_ratio", ascending=False)

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
      <th>model</th>
      <th>norm_dise_ratio</th>
      <th>norm_flux</th>
      <th>norm_prolif_ratio</th>
    </tr>
  </thead>
  <tbody>
    <tr>
      <th>AKGDm</th>
      <td>0.008</td>
      <td>0.14222</td>
      <td>biopsys</td>
      <td>7.03136</td>
      <td>0.0504142</td>
      <td>1</td>
    </tr>
    <tr>
      <th>SUCD1m</th>
      <td>0.008</td>
      <td>0.14222</td>
      <td>biopsys</td>
      <td>7.03136</td>
      <td>0.0504142</td>
      <td>1</td>
    </tr>
    <tr>
      <th>CO2t</th>
      <td>0.0322384</td>
      <td>0.573119</td>
      <td>biopsys</td>
      <td>1.74484</td>
      <td>0.0504142</td>
      <td>1</td>
    </tr>
    <tr>
      <th>r0354</th>
      <td>0.0562509</td>
      <td>1</td>
      <td>biopsys</td>
      <td>1.18692</td>
      <td>0.0598375</td>
      <td>1.18692</td>
    </tr>
    <tr>
      <th>r0357</th>
      <td>0.0562509</td>
      <td>1</td>
      <td>biopsys</td>
      <td>1.18692</td>
      <td>0.0598375</td>
      <td>1.18692</td>
    </tr>
    <tr>
      <th>HEX7</th>
      <td>0.0562509</td>
      <td>1</td>
      <td>biopsys</td>
      <td>1.18692</td>
      <td>0.0598375</td>
      <td>1.18692</td>
    </tr>
    <tr>
      <th>DUTPDPm</th>
      <td>0.0562509</td>
      <td>1</td>
      <td>biopsys</td>
      <td>1.18692</td>
      <td>0.0598375</td>
      <td>1.18692</td>
    </tr>
    <tr>
      <th>HEX4</th>
      <td>0.0562509</td>
      <td>1</td>
      <td>biopsys</td>
      <td>1.18692</td>
      <td>0.0598375</td>
      <td>1.18692</td>
    </tr>
    <tr>
      <th>r0360</th>
      <td>0.0562509</td>
      <td>1</td>
      <td>biopsys</td>
      <td>1.18692</td>
      <td>0.0598375</td>
      <td>1.18692</td>
    </tr>
    <tr>
      <th>r0361</th>
      <td>0.0562509</td>
      <td>1</td>
      <td>biopsys</td>
      <td>1.18692</td>
      <td>0.0598375</td>
      <td>1.18692</td>
    </tr>
    <tr>
      <th>r0358</th>
      <td>0.0562509</td>
      <td>1</td>
      <td>biopsys</td>
      <td>1.18692</td>
      <td>0.0598375</td>
      <td>1.18692</td>
    </tr>
    <tr>
      <th>HEX1</th>
      <td>0.0562509</td>
      <td>1</td>
      <td>biopsys</td>
      <td>1.18692</td>
      <td>0.0598375</td>
      <td>1.18692</td>
    </tr>
    <tr>
      <th>r0355</th>
      <td>0.0562509</td>
      <td>1</td>
      <td>biopsys</td>
      <td>1.18692</td>
      <td>0.0598375</td>
      <td>1.18692</td>
    </tr>
    <tr>
      <th>DM_dctp_n_</th>
      <td>0.0562509</td>
      <td>1</td>
      <td>biopsys</td>
      <td>1.16822</td>
      <td>0.0588951</td>
      <td>1.16822</td>
    </tr>
    <tr>
      <th>EX_utp_LPAREN_e_RPAREN_</th>
      <td>0.0562509</td>
      <td>1</td>
      <td>biopsys</td>
      <td>1.16822</td>
      <td>0.0588951</td>
      <td>1.16822</td>
    </tr>
    <tr>
      <th>TRIOK</th>
      <td>0.0562509</td>
      <td>1</td>
      <td>biopsys</td>
      <td>1.09346</td>
      <td>0.0551258</td>
      <td>1.09346</td>
    </tr>
    <tr>
      <th>r1527</th>
      <td>0.0562509</td>
      <td>0.0562509</td>
      <td>biopsys</td>
      <td>1.09346</td>
      <td>0.0551258</td>
      <td>1.09346</td>
    </tr>
    <tr>
      <th>PPAP</th>
      <td>0.0562509</td>
      <td>0.0562509</td>
      <td>biopsys</td>
      <td>1.09346</td>
      <td>0.0551258</td>
      <td>1.09346</td>
    </tr>
    <tr>
      <th>KHK2</th>
      <td>0.0562509</td>
      <td>1</td>
      <td>biopsys</td>
      <td>1.09346</td>
      <td>0.0551258</td>
      <td>1.09346</td>
    </tr>
    <tr>
      <th>NTD1</th>
      <td>0.0562509</td>
      <td>1</td>
      <td>biopsys</td>
      <td>1.09346</td>
      <td>0.0551258</td>
      <td>1.09346</td>
    </tr>
    <tr>
      <th>PIK5n</th>
      <td>0.0562509</td>
      <td>1</td>
      <td>biopsys</td>
      <td>1.09346</td>
      <td>0.0551258</td>
      <td>1.09346</td>
    </tr>
    <tr>
      <th>GLNS</th>
      <td>0.0562509</td>
      <td>1</td>
      <td>biopsys</td>
      <td>1.09346</td>
      <td>0.0551258</td>
      <td>1.09346</td>
    </tr>
    <tr>
      <th>NTD4</th>
      <td>0.0562509</td>
      <td>1</td>
      <td>biopsys</td>
      <td>1.09346</td>
      <td>0.0551258</td>
      <td>1.09346</td>
    </tr>
    <tr>
      <th>CHOLATEt3</th>
      <td>0.0562509</td>
      <td>0.0562509</td>
      <td>biopsys</td>
      <td>1.09346</td>
      <td>0.0551258</td>
      <td>1.09346</td>
    </tr>
    <tr>
      <th>PIK4n</th>
      <td>0.0562509</td>
      <td>1</td>
      <td>biopsys</td>
      <td>1.09346</td>
      <td>0.0551258</td>
      <td>1.09346</td>
    </tr>
    <tr>
      <th>r1162</th>
      <td>0.0562509</td>
      <td>1</td>
      <td>biopsys</td>
      <td>1.09346</td>
      <td>0.0551258</td>
      <td>1.09346</td>
    </tr>
    <tr>
      <th>EX_q10h2_LPAREN_e_RPAREN_</th>
      <td>0.0562509</td>
      <td>0.0562509</td>
      <td>biopsys</td>
      <td>1.08411</td>
      <td>0.0546547</td>
      <td>1.08411</td>
    </tr>
    <tr>
      <th>EX_ppa_LPAREN_e_RPAREN_</th>
      <td>0.0562509</td>
      <td>1</td>
      <td>biopsys</td>
      <td>1.08411</td>
      <td>0.0546547</td>
      <td>1.08411</td>
    </tr>
    <tr>
      <th>EX_atp_LPAREN_e_RPAREN_</th>
      <td>0.0562509</td>
      <td>1</td>
      <td>biopsys</td>
      <td>1.08411</td>
      <td>0.0546547</td>
      <td>1.08411</td>
    </tr>
    <tr>
      <th>RE0912C</th>
      <td>0.0562509</td>
      <td>1</td>
      <td>biopsys</td>
      <td>1.02869</td>
      <td>0.0518607</td>
      <td>1.02869</td>
    </tr>
    <tr>
      <th>TRDR</th>
      <td>0.0562509</td>
      <td>1</td>
      <td>biopsys</td>
      <td>1.02869</td>
      <td>0.0518607</td>
      <td>1.02869</td>
    </tr>
    <tr>
      <th>r0377</th>
      <td>0.0562509</td>
      <td>1</td>
      <td>biopsys</td>
      <td>1.02869</td>
      <td>0.0518607</td>
      <td>1.02869</td>
    </tr>
    <tr>
      <th>VITD3t</th>
      <td>0.0562509</td>
      <td>1</td>
      <td>biopsys</td>
      <td>1.0206</td>
      <td>0.0514528</td>
      <td>1.0206</td>
    </tr>
    <tr>
      <th>DNDPt55m</th>
      <td>0.0562509</td>
      <td>1</td>
      <td>biopsys</td>
      <td>1.01869</td>
      <td>0.0513566</td>
      <td>1.01869</td>
    </tr>
    <tr>
      <th>DNDPt56m</th>
      <td>0.0562509</td>
      <td>1</td>
      <td>biopsys</td>
      <td>1.01869</td>
      <td>0.0513566</td>
      <td>1.01869</td>
    </tr>
    <tr>
      <th>GMPS2</th>
      <td>0.0562509</td>
      <td>1</td>
      <td>biopsys</td>
      <td>1.01429</td>
      <td>0.0511348</td>
      <td>1.01429</td>
    </tr>
    <tr>
      <th>NTD7</th>
      <td>0.0562509</td>
      <td>1</td>
      <td>biopsys</td>
      <td>1.01175</td>
      <td>0.0510067</td>
      <td>1.01175</td>
    </tr>
    <tr>
      <th>NT5C</th>
      <td>0.0562509</td>
      <td>1</td>
      <td>biopsys</td>
      <td>1.01175</td>
      <td>0.0510067</td>
      <td>1.01175</td>
    </tr>
    <tr>
      <th>ADNK1m</th>
      <td>0.0562509</td>
      <td>1</td>
      <td>biopsys</td>
      <td>1.01175</td>
      <td>0.0510067</td>
      <td>1.01175</td>
    </tr>
    <tr>
      <th>ACITL</th>
      <td>0.0562509</td>
      <td>1</td>
      <td>biopsys</td>
      <td>1.00935</td>
      <td>0.0508854</td>
      <td>1.00935</td>
    </tr>
    <tr>
      <th>NTD11</th>
      <td>0.0562509</td>
      <td>1</td>
      <td>biopsys</td>
      <td>1.00935</td>
      <td>0.0508854</td>
      <td>1.00935</td>
    </tr>
    <tr>
      <th>TIGCRNe</th>
      <td>0.0562509</td>
      <td>1</td>
      <td>biopsys</td>
      <td>1.00935</td>
      <td>0.0508854</td>
      <td>1.00935</td>
    </tr>
    <tr>
      <th>r1516</th>
      <td>0.0562509</td>
      <td>1</td>
      <td>biopsys</td>
      <td>1.00935</td>
      <td>0.0508854</td>
      <td>1.00935</td>
    </tr>
    <tr>
      <th>NTD10</th>
      <td>0.0562509</td>
      <td>1</td>
      <td>biopsys</td>
      <td>1.00716</td>
      <td>0.0507753</td>
      <td>1.00716</td>
    </tr>
    <tr>
      <th>r0283</th>
      <td>0.0562509</td>
      <td>1</td>
      <td>biopsys</td>
      <td>1.00599</td>
      <td>0.0507163</td>
      <td>1.00599</td>
    </tr>
    <tr>
      <th>CDS</th>
      <td>0.0562509</td>
      <td>1</td>
      <td>biopsys</td>
      <td>1.00458</td>
      <td>0.0506449</td>
      <td>1.00458</td>
    </tr>
    <tr>
      <th>HMGLm</th>
      <td>0.0562509</td>
      <td>1</td>
      <td>biopsys</td>
      <td>1.00289</td>
      <td>0.05056</td>
      <td>1.00289</td>
    </tr>
    <tr>
      <th>r1515</th>
      <td>0.0562509</td>
      <td>1</td>
      <td>biopsys</td>
      <td>1.00082</td>
      <td>0.0504558</td>
      <td>1.00082</td>
    </tr>
  </tbody>
</table>
</div>



Aqui se analizan los modelos


```python
import pandas
from time import time

import cobra.test

from cobra.flux_analysis import (
    single_gene_deletion, single_reaction_deletion, double_gene_deletion,
    double_reaction_deletion)

deletion_results=single_reaction_deletion(model_MeanCancerBiopsy, model_MeanCancerBiopsy.reactions)

cancer_del = deletion_results[deletion_results.status != 'optimal']
deletion_results[deletion_results.status != 'optimal']

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
      <th>flux</th>
      <th>status</th>
    </tr>
  </thead>
  <tbody>
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

deletion_results=single_reaction_deletion(model_MeanNormalBiopsy, model_MeanNormalBiopsy.reactions)

normal_del = deletion_results[deletion_results.status != 'optimal']
deletion_results[deletion_results.status != 'optimal']

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
      <th>flux</th>
      <th>status</th>
    </tr>
  </thead>
  <tbody>
  </tbody>
</table>
</div>




```python
Targets_biopsys=get_drugable_targets(model_MeanNormalBiopsy, model_MeanCancerBiopsy, "biopsys" )
Targets_biopsys.to_csv("Drugable_targets_FULL_biopsys.tsv",sep="\t")
Targets_biopsys[ (Targets_biopsys["norm_dise_ratio"] > 1.0000001) ].sort_values(by="norm_dise_ratio", ascending=False)

```

    Common reactions size 1352
    Normal unique reactions size 252
    Disease unique reactions size 360





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
      <th>SUCD1m</th>
      <td>0.00779375</td>
      <td>0.138553</td>
      <td>and  and  and</td>
      <td>biopsys</td>
      <td>7.21744</td>
      <td>0.0504142</td>
      <td>1</td>
    </tr>
    <tr>
      <th>AKGDm</th>
      <td>0.00781568</td>
      <td>0.138943</td>
      <td>and  and  and</td>
      <td>biopsys</td>
      <td>7.19719</td>
      <td>0.0504142</td>
      <td>1</td>
    </tr>
    <tr>
      <th>CO2t</th>
      <td>0.0322384</td>
      <td>0.573119</td>
      <td></td>
      <td>biopsys</td>
      <td>1.74484</td>
      <td>0.0504142</td>
      <td>1</td>
    </tr>
    <tr>
      <th>HEX1</th>
      <td>0.0562509</td>
      <td>1</td>
      <td>or  or  or  or</td>
      <td>biopsys</td>
      <td>1.18692</td>
      <td>0.0598375</td>
      <td>1.18692</td>
    </tr>
    <tr>
      <th>r0354</th>
      <td>0.0562509</td>
      <td>1</td>
      <td>or  or  or</td>
      <td>biopsys</td>
      <td>1.18692</td>
      <td>0.0598375</td>
      <td>1.18692</td>
    </tr>
    <tr>
      <th>r0357</th>
      <td>0.0562509</td>
      <td>1</td>
      <td>or  or</td>
      <td>biopsys</td>
      <td>1.18692</td>
      <td>0.0598375</td>
      <td>1.18692</td>
    </tr>
    <tr>
      <th>DUTPDPm</th>
      <td>0.0562509</td>
      <td>1</td>
      <td>or</td>
      <td>biopsys</td>
      <td>1.18692</td>
      <td>0.0598375</td>
      <td>1.18692</td>
    </tr>
    <tr>
      <th>HEX7</th>
      <td>0.0562509</td>
      <td>1</td>
      <td>or  or  or  or</td>
      <td>biopsys</td>
      <td>1.18692</td>
      <td>0.0598375</td>
      <td>1.18692</td>
    </tr>
    <tr>
      <th>r0355</th>
      <td>0.0562509</td>
      <td>1</td>
      <td>or  or  or</td>
      <td>biopsys</td>
      <td>1.18692</td>
      <td>0.0598375</td>
      <td>1.18692</td>
    </tr>
    <tr>
      <th>HEX4</th>
      <td>0.0562509</td>
      <td>1</td>
      <td>or  or  or  or</td>
      <td>biopsys</td>
      <td>1.18692</td>
      <td>0.0598375</td>
      <td>1.18692</td>
    </tr>
    <tr>
      <th>r0360</th>
      <td>0.0562509</td>
      <td>1</td>
      <td>or  or</td>
      <td>biopsys</td>
      <td>1.18692</td>
      <td>0.0598375</td>
      <td>1.18692</td>
    </tr>
    <tr>
      <th>r0358</th>
      <td>0.0562509</td>
      <td>1</td>
      <td>or  or</td>
      <td>biopsys</td>
      <td>1.18692</td>
      <td>0.0598375</td>
      <td>1.18692</td>
    </tr>
    <tr>
      <th>r0361</th>
      <td>0.0562509</td>
      <td>1</td>
      <td>or  or</td>
      <td>biopsys</td>
      <td>1.18692</td>
      <td>0.0598375</td>
      <td>1.18692</td>
    </tr>
    <tr>
      <th>DM_dctp_n_</th>
      <td>0.0562509</td>
      <td>1</td>
      <td></td>
      <td>biopsys</td>
      <td>1.16822</td>
      <td>0.0588951</td>
      <td>1.16822</td>
    </tr>
    <tr>
      <th>EX_utp_LPAREN_e_RPAREN_</th>
      <td>0.0562509</td>
      <td>1</td>
      <td></td>
      <td>biopsys</td>
      <td>1.16822</td>
      <td>0.0588951</td>
      <td>1.16822</td>
    </tr>
    <tr>
      <th>KHK2</th>
      <td>0.0562509</td>
      <td>1</td>
      <td></td>
      <td>biopsys</td>
      <td>1.09346</td>
      <td>0.0551258</td>
      <td>1.09346</td>
    </tr>
    <tr>
      <th>PIK4n</th>
      <td>0.0562509</td>
      <td>1</td>
      <td></td>
      <td>biopsys</td>
      <td>1.09346</td>
      <td>0.0551258</td>
      <td>1.09346</td>
    </tr>
    <tr>
      <th>r1162</th>
      <td>0.0562509</td>
      <td>1</td>
      <td></td>
      <td>biopsys</td>
      <td>1.09346</td>
      <td>0.0551258</td>
      <td>1.09346</td>
    </tr>
    <tr>
      <th>PPAP</th>
      <td>0.0562509</td>
      <td>0.0562509</td>
      <td>or  or</td>
      <td>biopsys</td>
      <td>1.09346</td>
      <td>0.0551258</td>
      <td>1.09346</td>
    </tr>
    <tr>
      <th>TRIOK</th>
      <td>0.0562509</td>
      <td>1</td>
      <td></td>
      <td>biopsys</td>
      <td>1.09346</td>
      <td>0.0551258</td>
      <td>1.09346</td>
    </tr>
    <tr>
      <th>CHOLATEt3</th>
      <td>0.0562509</td>
      <td>0.0562509</td>
      <td>or</td>
      <td>biopsys</td>
      <td>1.09346</td>
      <td>0.0551258</td>
      <td>1.09346</td>
    </tr>
    <tr>
      <th>r1527</th>
      <td>0.0562509</td>
      <td>0.0562509</td>
      <td></td>
      <td>biopsys</td>
      <td>1.09346</td>
      <td>0.0551258</td>
      <td>1.09346</td>
    </tr>
    <tr>
      <th>NTD4</th>
      <td>0.0562509</td>
      <td>1</td>
      <td>or  or  or  or</td>
      <td>biopsys</td>
      <td>1.09346</td>
      <td>0.0551258</td>
      <td>1.09346</td>
    </tr>
    <tr>
      <th>PIK5n</th>
      <td>0.0562509</td>
      <td>1</td>
      <td></td>
      <td>biopsys</td>
      <td>1.09346</td>
      <td>0.0551258</td>
      <td>1.09346</td>
    </tr>
    <tr>
      <th>GLNS</th>
      <td>0.0562509</td>
      <td>1</td>
      <td>or</td>
      <td>biopsys</td>
      <td>1.09346</td>
      <td>0.0551258</td>
      <td>1.09346</td>
    </tr>
    <tr>
      <th>NTD1</th>
      <td>0.0562509</td>
      <td>1</td>
      <td>or</td>
      <td>biopsys</td>
      <td>1.09346</td>
      <td>0.0551258</td>
      <td>1.09346</td>
    </tr>
    <tr>
      <th>EX_atp_LPAREN_e_RPAREN_</th>
      <td>0.0562509</td>
      <td>1</td>
      <td></td>
      <td>biopsys</td>
      <td>1.08411</td>
      <td>0.0546547</td>
      <td>1.08411</td>
    </tr>
    <tr>
      <th>EX_q10h2_LPAREN_e_RPAREN_</th>
      <td>0.0562509</td>
      <td>0.0562509</td>
      <td></td>
      <td>biopsys</td>
      <td>1.08411</td>
      <td>0.0546547</td>
      <td>1.08411</td>
    </tr>
    <tr>
      <th>EX_ppa_LPAREN_e_RPAREN_</th>
      <td>0.0562509</td>
      <td>1</td>
      <td></td>
      <td>biopsys</td>
      <td>1.08411</td>
      <td>0.0546547</td>
      <td>1.08411</td>
    </tr>
    <tr>
      <th>TRDR</th>
      <td>0.0562509</td>
      <td>1</td>
      <td></td>
      <td>biopsys</td>
      <td>1.02869</td>
      <td>0.0518607</td>
      <td>1.02869</td>
    </tr>
    <tr>
      <th>r0377</th>
      <td>0.0562509</td>
      <td>1</td>
      <td></td>
      <td>biopsys</td>
      <td>1.02869</td>
      <td>0.0518607</td>
      <td>1.02869</td>
    </tr>
    <tr>
      <th>RE0912C</th>
      <td>0.0562509</td>
      <td>1</td>
      <td>or</td>
      <td>biopsys</td>
      <td>1.02869</td>
      <td>0.0518607</td>
      <td>1.02869</td>
    </tr>
    <tr>
      <th>VITD3t</th>
      <td>0.0562509</td>
      <td>1</td>
      <td></td>
      <td>biopsys</td>
      <td>1.0206</td>
      <td>0.0514528</td>
      <td>1.0206</td>
    </tr>
    <tr>
      <th>DNDPt55m</th>
      <td>0.0562509</td>
      <td>1</td>
      <td></td>
      <td>biopsys</td>
      <td>1.01869</td>
      <td>0.0513566</td>
      <td>1.01869</td>
    </tr>
    <tr>
      <th>DNDPt56m</th>
      <td>0.0562509</td>
      <td>1</td>
      <td></td>
      <td>biopsys</td>
      <td>1.01869</td>
      <td>0.0513566</td>
      <td>1.01869</td>
    </tr>
    <tr>
      <th>GMPS2</th>
      <td>0.0562509</td>
      <td>1</td>
      <td></td>
      <td>biopsys</td>
      <td>1.01429</td>
      <td>0.0511348</td>
      <td>1.01429</td>
    </tr>
    <tr>
      <th>ADNK1m</th>
      <td>0.0562509</td>
      <td>1</td>
      <td></td>
      <td>biopsys</td>
      <td>1.01175</td>
      <td>0.0510067</td>
      <td>1.01175</td>
    </tr>
    <tr>
      <th>NT5C</th>
      <td>0.0562509</td>
      <td>1</td>
      <td>or  or  or  or</td>
      <td>biopsys</td>
      <td>1.01175</td>
      <td>0.0510067</td>
      <td>1.01175</td>
    </tr>
    <tr>
      <th>NTD7</th>
      <td>0.0562509</td>
      <td>1</td>
      <td>or  or  or</td>
      <td>biopsys</td>
      <td>1.01175</td>
      <td>0.0510067</td>
      <td>1.01175</td>
    </tr>
    <tr>
      <th>NTD11</th>
      <td>0.0562509</td>
      <td>1</td>
      <td>or  or  or</td>
      <td>biopsys</td>
      <td>1.00935</td>
      <td>0.0508854</td>
      <td>1.00935</td>
    </tr>
    <tr>
      <th>TIGCRNe</th>
      <td>0.0562509</td>
      <td>1</td>
      <td></td>
      <td>biopsys</td>
      <td>1.00935</td>
      <td>0.0508854</td>
      <td>1.00935</td>
    </tr>
    <tr>
      <th>ACITL</th>
      <td>0.0562509</td>
      <td>1</td>
      <td></td>
      <td>biopsys</td>
      <td>1.00935</td>
      <td>0.0508854</td>
      <td>1.00935</td>
    </tr>
    <tr>
      <th>r1516</th>
      <td>0.0562509</td>
      <td>1</td>
      <td></td>
      <td>biopsys</td>
      <td>1.00935</td>
      <td>0.0508854</td>
      <td>1.00935</td>
    </tr>
    <tr>
      <th>NTD10</th>
      <td>0.0562509</td>
      <td>1</td>
      <td>or  or  or  or</td>
      <td>biopsys</td>
      <td>1.00716</td>
      <td>0.0507753</td>
      <td>1.00716</td>
    </tr>
    <tr>
      <th>r0283</th>
      <td>0.0562509</td>
      <td>1</td>
      <td></td>
      <td>biopsys</td>
      <td>1.00599</td>
      <td>0.0507163</td>
      <td>1.00599</td>
    </tr>
    <tr>
      <th>CDS</th>
      <td>0.0562509</td>
      <td>1</td>
      <td>or</td>
      <td>biopsys</td>
      <td>1.00458</td>
      <td>0.0506449</td>
      <td>1.00458</td>
    </tr>
    <tr>
      <th>HMGLm</th>
      <td>0.0562509</td>
      <td>1</td>
      <td>or</td>
      <td>biopsys</td>
      <td>1.00289</td>
      <td>0.05056</td>
      <td>1.00289</td>
    </tr>
    <tr>
      <th>r1515</th>
      <td>0.0562509</td>
      <td>1</td>
      <td></td>
      <td>biopsys</td>
      <td>1.00082</td>
      <td>0.0504558</td>
      <td>1.00082</td>
    </tr>
  </tbody>
</table>
</div>




```python
cpNorm=model_MeanNormalBiopsy
cpCanc=model_MeanCancerBiopsy

print(cpNorm.optimize())
print(cpCanc.optimize())

cpNorm.reactions.AKGDm.bounds=(-0.0,0.0)
cpNorm.reactions.AKGDm.bounds=(-0.0,0.0)

print(cpNorm.optimize())
print(cpCanc.optimize())

```

    <Solution 0.050 at 0x7f7ccded3048>
    <Solution 0.056 at 0x7f7ccded30b8>
    <Solution 0.050 at 0x7f7ccded39b0>
    <Solution 0.056 at 0x7f7ccded3898>



```python

```


```python
cbmodel=model_MeanNormalBiopsy.copy()
pbounds=cbmodel.reactions.SPHMDAc.bounds
cbmodel.reactions.PCHOLHSTDe.bounds=[-0.1,0.1]
print(model_MeanNormalBiopsy.optimize())
print(cbmodel.optimize())

```

Con la finalidad de probar la capacidad fisiologica del modelo se verificaron las posibles fuentes de obtención de energia. 
Incluir más contexto fisiologico

Incluir mayor contexto de normoxia e hipoxia


```python
# glucose aerobic CancerBiopsy


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

closedModel.reactions.EX_o2_LPAREN_e_RPAREN_.bounds=(-1000,0.0)
closedModel.reactions.EX_h2o_LPAREN_e_RPAREN_.bounds=(-1000,1000)
closedModel.reactions.EX_co2_LPAREN_e_RPAREN_.bounds=(-0.0,1000)
closedModel.reactions.EX_glc_LPAREN_e_RPAREN_.bounds=(-1,-1)
closedModel.optimize()
```


```python
# glucose anaerobic CancerBiopsy

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
    model_MaxNormalBiopsy, ["EX_co2_LPAREN_e_RPAREN_","EX_glc_LPAREN_e_RPAREN_","EX_glc_LPAREN_e_RPAREN_",])
```

Glucose - Lactate (exp)
Glucose - CO2
glutamine - Lactate
Glutamine - CO2
Glutamate trasferase - ROS
phospatidil Choline-**



```python
prod_env.head()
```


```python
%matplotlib inline

from mpl_toolkits.mplot3d import axes3d
import matplotlib.pyplot as plt

fig = plt.figure()
ax = fig.add_subplot(111, projection='3d')

X=prod_env["EX_co2_LPAREN_e_RPAREN_"]
Y=prod_env["EX_glc_LPAREN_e_RPAREN_"]
Z=prod_env["EX_gly_LPAREN_e_RPAREN_"]

surf=ax.plot_wireframe(X, Y, Z, rstride=10, cstride=10)

plt.show()
```


```python

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
        results[rx]["dise_prolif_ratio"]=dflx0
        
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
R_Biopsy=get_drugable_targets(
    normal_Model=model_MaxNormalBiopsy, 
    disease_Model=model_MaxCancerBiopsy, 
    model_name="CancerBiopsy")

#Glyceraldehyde-3-phosphate dehydrogenase GAPD
R_Biopsy = R_Biopsy[ (R_Biopsy["dise_prolif_ratio"] < 1.0) ].sort_values(by="norm_dise_ratio", ascending=False)
```


```python
R_Biopsy
```


```python
model_MaxCancerBiopsy.optimize()
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
model_MaxNormalBiopsy.reactions.PCHOLPr_hs.flux
```


```python
model_MaxNormalBiopsy.reactions.LYStiDF.name
```


```python
def metabolite_summary_dataframes(model=None, met="", solution=None, threshold=0.01, fva=False,
                       names=False, floatfmt='.3g'):
    """
    Get a summary of the production and consumption fluxes.

    This method requires the model for which this metabolite is a part
    to be solved.

    Parameters
    ----------
    solution : cobra.Solution, optional
        A previously solved model solution to use for generating the
        summary. If none provided (default), the summary method will
        resolve the model. Note that the solution object must match the
        model, i.e., changes to the model such as changed bounds,
        added or removed reactions are not taken into account by this
        method.
    threshold : float, optional
        Threshold below which fluxes are not reported.
    fva : pandas.DataFrame, float or None, optional
        Whether or not to include flux variability analysis in the output.
        If given, fva should either be a previous FVA solution matching
        the model or a float between 0 and 1 representing the
        fraction of the optimum objective to be searched.
    names : bool, optional
        Emit reaction and metabolite names rather than identifiers (default
        False).
    floatfmt : string, optional
        Format string for floats (default '.3g').

    """
    if names:
        emit = attrgetter('name')
    else:
        emit = attrgetter('id')
        
    if solution is None:
        model.slim_optimize(error_value=None)
        solution = get_solution(model, reactions=model.metabolites.get_by_id(met).reactions)

    rxns = sorted(met.reactions, key=attrgetter("id"))
    rxn_id = list()
    rxn_name = list()
    flux = list()
    reaction = list()
    for rxn in rxns:
        rxn_id.append(rxn.id)
        rxn_name.append(rxn.name)
        flux.append(solution[rxn.id] * rxn.metabolites[met])
      #  txt = rxn.build_reaction_string(use_metabolite_names=names)
       # reaction.append((txt, 40 if fva is not None else 50))
        reaction.append(rxn.reaction)
        
    flux_summary = pd.DataFrame({
        "id": rxn_name,
        "flux": flux,
        "reaction": reaction
    }, index=rxn_id)
    

    if fva is not None:
        if hasattr(fva, 'columns'):
            fva_results = fva
        else:
            fva_results = flux_variability_analysis(
                met.model, list(met.reactions), fraction_of_optimum=fva)
  
        for rxn in rxns:
            fmax = rxn.metabolites[met] * fva_results.at[rxn.id, "maximum"]
            fmin = rxn.metabolites[met] * fva_results.at[rxn.id, "minimum"]
            if abs(fmin) <= abs(fmax):
                flux_summary.at[rxn.id, "fmax"] = fmax
                flux_summary.at[rxn.id, "fmin"] = fmin
            else:
                # Reverse fluxes.
                flux_summary.at[rxn.id, "fmax"] = fmin
                flux_summary.at[rxn.id, "fmin"] = fmax

    assert flux_summary["flux"].sum() < 1E-6, "Error in flux balance"

    flux_summary = cobra.flux_analysis.summary._process_flux_dataframe(flux_summary, fva, threshold,
                                           floatfmt)

    flux_summary['percent'] = 0
    total_flux = flux_summary.loc[flux_summary.is_input, "flux"].sum()

    flux_summary.loc[flux_summary.is_input, 'percent'] = \
        flux_summary.loc[flux_summary.is_input, 'flux'] / total_flux
    flux_summary.loc[~flux_summary.is_input, 'percent'] = \
        flux_summary.loc[~flux_summary.is_input, 'flux'] / total_flux

    flux_summary['percent'] = flux_summary.percent.apply(
        lambda x: '{:.0%}'.format(x))

    if fva is not None:
        flux_table = tabulate(
            flux_summary.loc[:, ['percent', 'flux', 'fva_fmt', 'id',
                                 'reaction']].values, floatfmt=floatfmt,
            headers=['%', 'FLUX', 'RANGE', 'RXN ID', 'REACTION']).split('\n')
    else:
        flux_table = tabulate(
            flux_summary.loc[:, ['percent', 'flux', 'id', 'reaction']].values,
            floatfmt=floatfmt, headers=['%', 'FLUX', 'RXN ID', 'REACTION']
        ).split('\n')

    flux_table_head = flux_table[:2]


    return(flux_summary)
          # [flux_summary.is_input.values],pd.np.array(flux_table[2:])[~flux_summary.is_input.values]
```


```python
fba_solution = cobra.flux_analysis.pfba(model_MaxCancerBiopsy)
fba_solution.metabolites.pchol_hs_c.summary()
```


```python
rxs=model_MaxCancerBiopsy.metabolites.get_by_id("pchol_hs_c").reactions

for rx in rxs:
    metbols=model_MaxCancerBiopsy.reactions.get_by_id(rx.id).metabolites
    for met in metbols:
        if(met.formula_weight>100):
            print(met.id, met.formula_weight)
```


```python
mat=cobra.util.array.create_stoichiometric_matrix(model_MaxCancerBiopsy, array_type='DataFrame', dtype=None)

mat.values.T*mat.values
```


```python
start=metabolite_summary_dataframes(model,model_MaxCancerBiopsy.metabolites.get_by_id("pchol_hs_c"), solution=fba_solution, threshold=0.01, fva=True, names=True)
start=start[start["flux"]!=0]
up_rx=start.index[start["is_input"]==True]
down_rx=start.index[start["is_input"]==False]

met_matrix=[[0]]
names=["pchol_hs_c"]

for rx in up_rx:
    metbols = model_MaxCancerBiopsy.reactions.get_by_id(rx).metabolites
    for met in metbols:
        next_met=[]
        if(met.formula_weight>100):
            names.append(met.id)
            
    df = pd.DataFrame([[1, 2], [3, 4]], columns=list('AB'))
    
    
    for met in metbols:
        met_matrix.ad
    
    
for rx in down_rx:
    metbols = model_MaxCancerBiopsy.reactions.get_by_id(rx).metabolites
    

```


```python
def get_pathway(model, met)
    start=metabolite_summary_dataframes(model,model_MaxCancerBiopsy.metabolites.get_by_id(met), solution=fba_solution, threshold=0.01, fva=True, names=True)
    up_rx=start.index["is_input"==True]
    down_rx=start.index["is_input"==True]
```


```python

```


```python
print(model_MaxNormalBiopsy.reactions.PCHOLPr_hs.reaction, model_MaxNormalBiopsy.reactions.PCHOLPr_hs.flux )
print(model_MaxNormalBiopsy.reactions.PCHOLP_hs.reaction)
print(model_MaxNormalBiopsy.reactions.PCHOLPg_hs.reaction)
print(model_MaxNormalBiopsy.metabolites.chol_r.reactions)
print(model_MaxNormalBiopsy.metabolites.pa_hs_r.reactions)


```


```python
#Análisis of choline phosphatase
Cpfba_solution = cobra.flux_analysis.pfba(model_MaxCancerBiopsy)
Npfba_solution = cobra.flux_analysis.pfba(model_MaxNormalBiopsy)


print("\nCancerModel\n",)
Cpfba_solution.metabolites.pchol_hs_r.summary()
print("\nNormalModel\n",)
Npfba_solution.metabolites.pchol_hs_r.summary()

```


```python
for i in model_MaxNormalBiopsy.exchanges:
    print(i)
```


```python
def get_pathway(model, metabolite,heavy=True, parsimonius=True):
    if(parsimonius):
        fba_solution = cobra.flux_analysis.pfba(model)
    else:
        fba_solution = model.optimize()

    #return(model.metabolites.get_by_id(metabolite).summary())

    print(fba_solution.metabolites.get_by_id(metabolite).summary(fva=True, names=True))
    
    #sum=cobra.flux_analysis.summary.metabolite_summary(met=metabolite, solution=fba_solution, threshold=0.01, fva=True, names=True))

get_pathway(model=model_MaxNormalBiopsy,metabolite="pchol_hs_r")
```


```python
#Análisis of choline phosphatase
print("\nCancerModel\n",)
model_MaxCancerBiopsy.metabolites.pchol_hs_r.summary()
print("\nNormalModel\n",)
model_MaxNormalBiopsy.metabolites.pchol_hs_r.summary()
```


```python
rx=model_MaxNormalBiopsy.metabolites.lac_D_e.reactions
for r in rx:
    print(r.id,r.name,r.reaction)
```


```python
#Comparación de modelos
model_MaxCancerBiopsy.reactions.get_by_id("biomass_reaction").reaction
```
