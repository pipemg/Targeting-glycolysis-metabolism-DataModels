
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


In order to create the models of cervix biopsis samples we use
the MODEL1603150001 of biomodels, recon2.2. with a rule correction

```python
import cobra


#!wget http://www.ebi.ac.uk/biomodels-main/download?mid=MODEL1603150001 -O recon2.2.xml
Recon2 = cobra.io.read_sbml_model("Models/recon2.2.xml")
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
for rxex in Recon2.exchanges:
    rxex.bounds=(-0.001,10)
    
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
Recon2.reactions.EX_lac_D_LPAREN_e_RPAREN_.bounds=(0,0.001)


Recon2.reactions.biomass_reaction.bounds=(0.008, 0.05)
```


```python
Recon2.optimize()
```

    <Solution 0.050 at 0x7f18e3147278>



```python
Recon2.summary()
```

    IN FLUXES                  OUT FLUXES                 OBJECTIVES
    -------------------------  -------------------------  ---------------------
    o2_e             0.301     hco3_e           0.0857    biomass_reac...  0.05
    ile_L_e          0.0869    ala_B_e          0.0736
    so4_e            0.0673    slfcys_e         0.0673
    cys_L_e          0.0626    acac_e           0.0536
    arg_L_e          0.0466    nh4_e            0.0461
    lys_L_e          0.0296    ac_e             0.0306
    gln_L_e          0.0236    urea_e           0.0286
    leu_L_e          0.0213    co2_e            0.0235
    val_L_e          0.0176    succ_e           0.0176
    his_L_e          0.0174    hxan_e           0.0138
    thr_L_e          0.0156    akg_e            0.0109
    phe_L_e          0.013     drib_e           0.00647
    tyr_L_e          0.0103    urate_e          0.0057
    ala_L_e          0.01      ura_e            0.00403
    asn_L_e          0.01      acetone_e        0.00332
    asp_L_e          0.01      pan4p_e          0.003
    pro_L_e          0.01      pnto_R_e         0.003
    chol_e           0.00706   taur_c           0.003
    met_L_e          0.00665   HC02191_e        0.002
    glc_D_e          0.0056    dttp_e           0.002
    biomass_other_c  0.0027    man_e            0.00113
    inost_e          0.00117   13_cis_retng...  0.001
    3mop_e           0.001     23cump_e         0.001
    3ump_e           0.001     5oxpro_e         0.001
    4hpro_LT_e       0.001     CE1950_e         0.001
    4hpro_LT_m       0.001     CE2250_e         0.001
    4mop_e           0.001     CE5798_e         0.001
    CE0074_e         0.001     HC01610_e        0.001
    CE5560_e         0.001     Tyr_ggn_e        0.001
    CE5797_e         0.001     ach_e            0.001
    HC00250_e        0.001     c81crn_e         0.001
    HC00342_e        0.001     ctp_e            0.001
    HC01609_e        0.001     dgmp_e           0.001
    HC02192_e        0.001     epoxtac_e        0.001
    Lcystin_c        0.001     itp_e            0.001
    Lcystin_e        0.001     leuktrE4_e       0.001
    Rtotal_e         0.001     meoh_e           0.001
    adn_e            0.001     ribflv_e         0.001
    adprib_e         0.001     tchola_e         0.001
    adrn_e           0.001     ethamp_r         0.000889
    ahcys_e          0.001     Asn_X_Ser_Thr_l  0.000491
    aicar_e          0.001     1glyc_hs_e       0.000271
    alaala_e         0.001     dttp_m           0.000247
    amp_e            0.001
    arachcoa_e       0.001
    arachd_e         0.001
    asp_D_e          0.001
    atp_e            0.001
    c81coa_c         0.001
    cdp_e            0.001
    cgly_e           0.001
    cit_e            0.001
    cmp_e            0.001
    coa_e            0.001
    creat_e          0.001
    crn_e            0.001
    cynt_e           0.001
    cytd_e           0.001
    dad_2_e          0.001
    dag_hs_e         0.001
    datp_m           0.001
    datp_n           0.001
    dca_e            0.001
    dcmp_e           0.001
    dcsptn1_e        0.001
    dctp_n           0.001
    ddca_e           0.001
    dgsn_e           0.001
    dgtp_e           0.001
    dgtp_m           0.001
    dgtp_n           0.001
    dlnlcg_e         0.001
    dpcoa_e          0.001
    dtdp_e           0.001
    dtmp_e           0.001
    dttp_n           0.001
    fad_e            0.001
    fum_e            0.001
    gal_e            0.001
    gchola_e         0.001
    gdp_e            0.001
    gluala_e         0.001
    glyb_e           0.001
    glyc3p_e         0.001
    glygly_e         0.001
    glygn2_e         0.001
    glyleu_e         0.001
    glypro_e         0.001
    glysar_e         0.001
    gmp_e            0.001
    gsn_e            0.001
    gthrd_e          0.001
    hdca_e           0.001
    hdcea_e          0.001
    idp_e            0.001
    ins_e            0.001
    leugly_e         0.001
    leuktrC4_e       0.001
    leuleu_e         0.001
    lnlncg_e         0.001
    lpchol_hs_e      0.001
    mag_hs_e         0.001
    mal_L_e          0.001
    malcoa_e         0.001
    o2s_e            0.001
    ocdca_e          0.001
    ocdcea_e         0.001
    octdececoa_c     0.001
    orn_e            0.001
    orot_e           0.001
    pchol_hs_e       0.001
    pe_hs_e          0.001
    pe_hs_r          0.001
    pglyc_hs_e       0.001
    pmtcoa_r         0.001
    ppi_e            0.001
    progly_e         0.001
    prpp_e           0.001
    ps_hs_e          0.001
    retfa_e          0.001
    sarcs_e          0.001
    so3_e            0.001
    tacr_e           0.001
    tag_hs_e         0.001
    thym_e           0.001
    ttdca_e          0.001
    udp_e            0.001
    ump_e            0.001
    uri_e            0.001
    utp_e            0.001
    xmp_e            0.001
    sphs1p_e         0.000889
    sph1p_e          0.000874
    strch1_e         0.000812
    docosac_e        0.000704
    trp_L_e          0.000665
    n2m2nmasn_e      0.000491
    octa_e           0.000222
    cbasp_e          0.000129
    adp_e            0.00011



```python
import pandas as pd
ub_scores_matrix =pd.read_csv("Binary_CellLines_RECON2_2.tsv",sep=",",index_col=0,header=0)
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
%matplotlib inline
ub_scores_matrix["MeanHeLa"].hist(bins=20)

```




    <matplotlib.axes._subplots.AxesSubplot at 0x7f18e1e62e80>




![png](output_8_1.png)



```python
%matplotlib inline

ub_scores_matrix["Meankeratinocytes"].hist(bins=10)
```

![png](output_9_1.png)



```python
for model in ub_scores_matrix.columns:
    #High Confidence rate 
       
    HC=ub_scores_matrix[model][ub_scores_matrix[model] >= .75].index.tolist()
    confidence_scores_matrix[model][HC]=3
        
    ##Medium Confidence Rate 
    MC=ub_scores_matrix[model][(ub_scores_matrix[model] < 0.75) & (ub_scores_matrix[model] >= .5) ].index.tolist()
    confidence_scores_matrix[model][MC]=2
        
    #Low
    LC=ub_scores_matrix[model][(ub_scores_matrix[model] < 0.5) & (ub_scores_matrix[model] >= 0.25) ].index.tolist()
    confidence_scores_matrix[model][LC]=1
    
    #unknown
    ZC=ub_scores_matrix[model][(ub_scores_matrix[model] <.25) & (ub_scores_matrix[model] > 0) ].index.tolist()
    confidence_scores_matrix[model][ZC]=0
    
    #Negative
    NC=ub_scores_matrix[model][(ub_scores_matrix[model] ==0)].index.tolist()
    confidence_scores_matrix[model][NC]=-1
 
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
      <td>-1.0</td>
      <td>-1.0</td>
      <td>-1.0</td>
      <td>-1.0</td>
      <td>-1.0</td>
      <td>-1.0</td>
      <td>-1.0</td>
      <td>-1.0</td>
    </tr>
    <tr>
      <th>HGNC:1027</th>
      <td>-1.0</td>
      <td>-1.0</td>
      <td>0.0</td>
      <td>-1.0</td>
      <td>0.0</td>
      <td>-1.0</td>
      <td>0.0</td>
      <td>-1.0</td>
    </tr>
    <tr>
      <th>HGNC:10293</th>
      <td>-1.0</td>
      <td>-1.0</td>
      <td>3.0</td>
      <td>3.0</td>
      <td>3.0</td>
      <td>3.0</td>
      <td>2.0</td>
      <td>1.0</td>
    </tr>
    <tr>
      <th>HGNC:10297</th>
      <td>-1.0</td>
      <td>-1.0</td>
      <td>2.0</td>
      <td>2.0</td>
      <td>2.0</td>
      <td>2.0</td>
      <td>1.0</td>
      <td>1.0</td>
    </tr>
    <tr>
      <th>HGNC:10451</th>
      <td>3.0</td>
      <td>3.0</td>
      <td>3.0</td>
      <td>3.0</td>
      <td>3.0</td>
      <td>3.0</td>
      <td>3.0</td>
      <td>3.0</td>
    </tr>
    <tr>
      <th>HGNC:10452</th>
      <td>3.0</td>
      <td>3.0</td>
      <td>3.0</td>
      <td>3.0</td>
      <td>3.0</td>
      <td>3.0</td>
      <td>3.0</td>
      <td>3.0</td>
    </tr>
    <tr>
      <th>HGNC:1047</th>
      <td>-1.0</td>
      <td>-1.0</td>
      <td>-1.0</td>
      <td>-1.0</td>
      <td>-1.0</td>
      <td>-1.0</td>
      <td>-1.0</td>
      <td>-1.0</td>
    </tr>
    <tr>
      <th>HGNC:1048</th>
      <td>-1.0</td>
      <td>-1.0</td>
      <td>-1.0</td>
      <td>-1.0</td>
      <td>-1.0</td>
      <td>-1.0</td>
      <td>-1.0</td>
      <td>-1.0</td>
    </tr>
    <tr>
      <th>HGNC:10536</th>
      <td>-1.0</td>
      <td>-1.0</td>
      <td>-1.0</td>
      <td>-1.0</td>
      <td>-1.0</td>
      <td>-1.0</td>
      <td>-1.0</td>
      <td>-1.0</td>
    </tr>
    <tr>
      <th>HGNC:10540</th>
      <td>-1.0</td>
      <td>3.0</td>
      <td>3.0</td>
      <td>3.0</td>
      <td>3.0</td>
      <td>3.0</td>
      <td>1.0</td>
      <td>3.0</td>
    </tr>
    <tr>
      <th>HGNC:10545</th>
      <td>2.0</td>
      <td>3.0</td>
      <td>3.0</td>
      <td>3.0</td>
      <td>3.0</td>
      <td>3.0</td>
      <td>3.0</td>
      <td>3.0</td>
    </tr>
    <tr>
      <th>HGNC:10547</th>
      <td>-1.0</td>
      <td>1.0</td>
      <td>3.0</td>
      <td>3.0</td>
      <td>3.0</td>
      <td>3.0</td>
      <td>1.0</td>
      <td>2.0</td>
    </tr>
    <tr>
      <th>HGNC:10571</th>
      <td>3.0</td>
      <td>3.0</td>
      <td>3.0</td>
      <td>3.0</td>
      <td>3.0</td>
      <td>3.0</td>
      <td>3.0</td>
      <td>3.0</td>
    </tr>
    <tr>
      <th>HGNC:10606</th>
      <td>1.0</td>
      <td>-1.0</td>
      <td>3.0</td>
      <td>3.0</td>
      <td>3.0</td>
      <td>3.0</td>
      <td>2.0</td>
      <td>2.0</td>
    </tr>
    <tr>
      <th>HGNC:1062</th>
      <td>-1.0</td>
      <td>-1.0</td>
      <td>3.0</td>
      <td>3.0</td>
      <td>3.0</td>
      <td>3.0</td>
      <td>1.0</td>
      <td>2.0</td>
    </tr>
    <tr>
      <th>HGNC:1063</th>
      <td>-1.0</td>
      <td>-1.0</td>
      <td>3.0</td>
      <td>3.0</td>
      <td>3.0</td>
      <td>3.0</td>
      <td>2.0</td>
      <td>2.0</td>
    </tr>
    <tr>
      <th>HGNC:10680</th>
      <td>3.0</td>
      <td>3.0</td>
      <td>3.0</td>
      <td>3.0</td>
      <td>3.0</td>
      <td>3.0</td>
      <td>3.0</td>
      <td>3.0</td>
    </tr>
    <tr>
      <th>HGNC:10681</th>
      <td>0.0</td>
      <td>-1.0</td>
      <td>3.0</td>
      <td>3.0</td>
      <td>3.0</td>
      <td>3.0</td>
      <td>2.0</td>
      <td>2.0</td>
    </tr>
    <tr>
      <th>HGNC:10682</th>
      <td>3.0</td>
      <td>-1.0</td>
      <td>3.0</td>
      <td>3.0</td>
      <td>3.0</td>
      <td>3.0</td>
      <td>3.0</td>
      <td>2.0</td>
    </tr>
    <tr>
      <th>HGNC:10683</th>
      <td>-1.0</td>
      <td>-1.0</td>
      <td>3.0</td>
      <td>3.0</td>
      <td>3.0</td>
      <td>3.0</td>
      <td>2.0</td>
      <td>2.0</td>
    </tr>
    <tr>
      <th>HGNC:10691</th>
      <td>-1.0</td>
      <td>-1.0</td>
      <td>-1.0</td>
      <td>-1.0</td>
      <td>-1.0</td>
      <td>-1.0</td>
      <td>-1.0</td>
      <td>-1.0</td>
    </tr>
    <tr>
      <th>HGNC:10761</th>
      <td>-1.0</td>
      <td>-1.0</td>
      <td>-1.0</td>
      <td>-1.0</td>
      <td>-1.0</td>
      <td>-1.0</td>
      <td>-1.0</td>
      <td>-1.0</td>
    </tr>
    <tr>
      <th>HGNC:108</th>
      <td>-1.0</td>
      <td>-1.0</td>
      <td>-1.0</td>
      <td>-1.0</td>
      <td>-1.0</td>
      <td>-1.0</td>
      <td>-1.0</td>
      <td>-1.0</td>
    </tr>
    <tr>
      <th>HGNC:10817</th>
      <td>-1.0</td>
      <td>-1.0</td>
      <td>1.0</td>
      <td>3.0</td>
      <td>1.0</td>
      <td>3.0</td>
      <td>0.0</td>
      <td>2.0</td>
    </tr>
    <tr>
      <th>HGNC:10818</th>
      <td>-1.0</td>
      <td>-1.0</td>
      <td>2.0</td>
      <td>0.0</td>
      <td>2.0</td>
      <td>0.0</td>
      <td>1.0</td>
      <td>0.0</td>
    </tr>
    <tr>
      <th>HGNC:10850</th>
      <td>-1.0</td>
      <td>-1.0</td>
      <td>1.0</td>
      <td>1.0</td>
      <td>1.0</td>
      <td>1.0</td>
      <td>0.0</td>
      <td>0.0</td>
    </tr>
    <tr>
      <th>HGNC:10852</th>
      <td>3.0</td>
      <td>-1.0</td>
      <td>3.0</td>
      <td>3.0</td>
      <td>3.0</td>
      <td>3.0</td>
      <td>3.0</td>
      <td>2.0</td>
    </tr>
    <tr>
      <th>HGNC:10856</th>
      <td>-1.0</td>
      <td>-1.0</td>
      <td>-1.0</td>
      <td>-1.0</td>
      <td>-1.0</td>
      <td>-1.0</td>
      <td>-1.0</td>
      <td>-1.0</td>
    </tr>
    <tr>
      <th>HGNC:10860</th>
      <td>-1.0</td>
      <td>-1.0</td>
      <td>0.0</td>
      <td>2.0</td>
      <td>0.0</td>
      <td>2.0</td>
      <td>0.0</td>
      <td>1.0</td>
    </tr>
    <tr>
      <th>HGNC:10862</th>
      <td>-1.0</td>
      <td>-1.0</td>
      <td>1.0</td>
      <td>-1.0</td>
      <td>1.0</td>
      <td>-1.0</td>
      <td>0.0</td>
      <td>-1.0</td>
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
      <td>...</td>
    </tr>
    <tr>
      <th>HGNC:9462</th>
      <td>0.0</td>
      <td>-1.0</td>
      <td>3.0</td>
      <td>2.0</td>
      <td>3.0</td>
      <td>2.0</td>
      <td>2.0</td>
      <td>1.0</td>
    </tr>
    <tr>
      <th>HGNC:9463</th>
      <td>-1.0</td>
      <td>-1.0</td>
      <td>-1.0</td>
      <td>-1.0</td>
      <td>-1.0</td>
      <td>-1.0</td>
      <td>-1.0</td>
      <td>-1.0</td>
    </tr>
    <tr>
      <th>HGNC:9465</th>
      <td>2.0</td>
      <td>1.0</td>
      <td>3.0</td>
      <td>3.0</td>
      <td>3.0</td>
      <td>3.0</td>
      <td>3.0</td>
      <td>2.0</td>
    </tr>
    <tr>
      <th>HGNC:95</th>
      <td>-1.0</td>
      <td>-1.0</td>
      <td>3.0</td>
      <td>3.0</td>
      <td>3.0</td>
      <td>3.0</td>
      <td>1.0</td>
      <td>2.0</td>
    </tr>
    <tr>
      <th>HGNC:9577</th>
      <td>-1.0</td>
      <td>-1.0</td>
      <td>3.0</td>
      <td>3.0</td>
      <td>3.0</td>
      <td>3.0</td>
      <td>1.0</td>
      <td>1.0</td>
    </tr>
    <tr>
      <th>HGNC:9587</th>
      <td>3.0</td>
      <td>1.0</td>
      <td>3.0</td>
      <td>3.0</td>
      <td>3.0</td>
      <td>3.0</td>
      <td>3.0</td>
      <td>2.0</td>
    </tr>
    <tr>
      <th>HGNC:9588</th>
      <td>-1.0</td>
      <td>-1.0</td>
      <td>3.0</td>
      <td>3.0</td>
      <td>3.0</td>
      <td>3.0</td>
      <td>2.0</td>
      <td>2.0</td>
    </tr>
    <tr>
      <th>HGNC:9592</th>
      <td>-1.0</td>
      <td>-1.0</td>
      <td>-1.0</td>
      <td>-1.0</td>
      <td>-1.0</td>
      <td>-1.0</td>
      <td>-1.0</td>
      <td>-1.0</td>
    </tr>
    <tr>
      <th>HGNC:9599</th>
      <td>2.0</td>
      <td>-1.0</td>
      <td>2.0</td>
      <td>3.0</td>
      <td>2.0</td>
      <td>3.0</td>
      <td>2.0</td>
      <td>2.0</td>
    </tr>
    <tr>
      <th>HGNC:9603</th>
      <td>-1.0</td>
      <td>-1.0</td>
      <td>-1.0</td>
      <td>0.0</td>
      <td>-1.0</td>
      <td>0.0</td>
      <td>-1.0</td>
      <td>0.0</td>
    </tr>
    <tr>
      <th>HGNC:9604</th>
      <td>-1.0</td>
      <td>-1.0</td>
      <td>1.0</td>
      <td>1.0</td>
      <td>1.0</td>
      <td>1.0</td>
      <td>0.0</td>
      <td>0.0</td>
    </tr>
    <tr>
      <th>HGNC:9605</th>
      <td>0.0</td>
      <td>3.0</td>
      <td>2.0</td>
      <td>3.0</td>
      <td>2.0</td>
      <td>3.0</td>
      <td>1.0</td>
      <td>3.0</td>
    </tr>
    <tr>
      <th>HGNC:964</th>
      <td>-1.0</td>
      <td>-1.0</td>
      <td>-1.0</td>
      <td>-1.0</td>
      <td>-1.0</td>
      <td>-1.0</td>
      <td>-1.0</td>
      <td>-1.0</td>
    </tr>
    <tr>
      <th>HGNC:9689</th>
      <td>2.0</td>
      <td>-1.0</td>
      <td>3.0</td>
      <td>3.0</td>
      <td>3.0</td>
      <td>3.0</td>
      <td>3.0</td>
      <td>1.0</td>
    </tr>
    <tr>
      <th>HGNC:9721</th>
      <td>-1.0</td>
      <td>-1.0</td>
      <td>0.0</td>
      <td>-1.0</td>
      <td>0.0</td>
      <td>-1.0</td>
      <td>0.0</td>
      <td>-1.0</td>
    </tr>
    <tr>
      <th>HGNC:9722</th>
      <td>-1.0</td>
      <td>-1.0</td>
      <td>3.0</td>
      <td>3.0</td>
      <td>3.0</td>
      <td>3.0</td>
      <td>2.0</td>
      <td>2.0</td>
    </tr>
    <tr>
      <th>HGNC:9723</th>
      <td>-1.0</td>
      <td>-1.0</td>
      <td>0.0</td>
      <td>2.0</td>
      <td>0.0</td>
      <td>2.0</td>
      <td>0.0</td>
      <td>1.0</td>
    </tr>
    <tr>
      <th>HGNC:9725</th>
      <td>3.0</td>
      <td>2.0</td>
      <td>3.0</td>
      <td>3.0</td>
      <td>3.0</td>
      <td>3.0</td>
      <td>3.0</td>
      <td>3.0</td>
    </tr>
    <tr>
      <th>HGNC:9726</th>
      <td>-1.0</td>
      <td>-1.0</td>
      <td>-1.0</td>
      <td>-1.0</td>
      <td>-1.0</td>
      <td>-1.0</td>
      <td>-1.0</td>
      <td>-1.0</td>
    </tr>
    <tr>
      <th>HGNC:9752</th>
      <td>-1.0</td>
      <td>-1.0</td>
      <td>3.0</td>
      <td>3.0</td>
      <td>3.0</td>
      <td>3.0</td>
      <td>2.0</td>
      <td>2.0</td>
    </tr>
    <tr>
      <th>HGNC:9755</th>
      <td>-1.0</td>
      <td>-1.0</td>
      <td>0.0</td>
      <td>-1.0</td>
      <td>0.0</td>
      <td>-1.0</td>
      <td>0.0</td>
      <td>-1.0</td>
    </tr>
    <tr>
      <th>HGNC:976</th>
      <td>-1.0</td>
      <td>-1.0</td>
      <td>2.0</td>
      <td>3.0</td>
      <td>2.0</td>
      <td>3.0</td>
      <td>1.0</td>
      <td>2.0</td>
    </tr>
    <tr>
      <th>HGNC:977</th>
      <td>-1.0</td>
      <td>-1.0</td>
      <td>0.0</td>
      <td>-1.0</td>
      <td>0.0</td>
      <td>-1.0</td>
      <td>0.0</td>
      <td>-1.0</td>
    </tr>
    <tr>
      <th>HGNC:986</th>
      <td>-1.0</td>
      <td>-1.0</td>
      <td>-1.0</td>
      <td>-1.0</td>
      <td>-1.0</td>
      <td>-1.0</td>
      <td>-1.0</td>
      <td>-1.0</td>
    </tr>
    <tr>
      <th>HGNC:987</th>
      <td>-1.0</td>
      <td>-1.0</td>
      <td>2.0</td>
      <td>3.0</td>
      <td>2.0</td>
      <td>3.0</td>
      <td>1.0</td>
      <td>1.0</td>
    </tr>
    <tr>
      <th>HGNC:9919</th>
      <td>-1.0</td>
      <td>-1.0</td>
      <td>-1.0</td>
      <td>1.0</td>
      <td>-1.0</td>
      <td>1.0</td>
      <td>-1.0</td>
      <td>0.0</td>
    </tr>
    <tr>
      <th>HGNC:9920</th>
      <td>-1.0</td>
      <td>-1.0</td>
      <td>-1.0</td>
      <td>-1.0</td>
      <td>-1.0</td>
      <td>-1.0</td>
      <td>-1.0</td>
      <td>-1.0</td>
    </tr>
    <tr>
      <th>HGNC:9922</th>
      <td>-1.0</td>
      <td>-1.0</td>
      <td>-1.0</td>
      <td>-1.0</td>
      <td>-1.0</td>
      <td>-1.0</td>
      <td>-1.0</td>
      <td>-1.0</td>
    </tr>
    <tr>
      <th>HGNC:9940</th>
      <td>-1.0</td>
      <td>-1.0</td>
      <td>-1.0</td>
      <td>-1.0</td>
      <td>-1.0</td>
      <td>-1.0</td>
      <td>-1.0</td>
      <td>-1.0</td>
    </tr>
    <tr>
      <th>HGNC:9959</th>
      <td>-1.0</td>
      <td>-1.0</td>
      <td>-1.0</td>
      <td>-1.0</td>
      <td>-1.0</td>
      <td>-1.0</td>
      <td>-1.0</td>
      <td>-1.0</td>
    </tr>
  </tbody>
</table>
<p>1595 rows × 8 columns</p>
</div>




```python
%matplotlib inline

confidence_scores_matrix["MaxHeLa"].hist(bins=10)
```

![png](output_12_1.png)



```python
from corda import reaction_confidence


conf_HeLaMean = {}
conf_keratinocytesMean= {}

for r in Recon2.reactions:
    if(r.gene_reaction_rule!=''):
        conf_HeLaMean[r.id]=reaction_confidence(r.gene_reaction_rule,confidence_scores_matrix["MaxHeLa"])
        conf_keratinocytesMean[r.id]=reaction_confidence(r.gene_reaction_rule,confidence_scores_matrix["Maxkeratinocytes"])
    else:
        conf_HeLaMean[r.id]=1
        conf_keratinocytesMean[r.id]=1
```


```python

 
conf_HeLaMean["EX_arg_L_LPAREN_e_RPAREN_"]=3
conf_HeLaMean["EX_cys_L_LPAREN_e_RPAREN_"]=3
conf_HeLaMean["EX_gln_L_LPAREN_e_RPAREN_"]=3
conf_HeLaMean["EX_gly_LPAREN_e_RPAREN_"]=3
conf_HeLaMean["EX_his_L_LPAREN_e_RPAREN_"]=3
conf_HeLaMean["EX_ile_L_LPAREN_e_RPAREN_"]=3
conf_HeLaMean["EX_leu_L_LPAREN_e_RPAREN_"]=3
conf_HeLaMean["EX_lys_L_LPAREN_e_RPAREN_"]=3
conf_HeLaMean["EX_met_L_LPAREN_e_RPAREN_"]=3
conf_HeLaMean["EX_phe_L_LPAREN_e_RPAREN_"]=3
conf_HeLaMean["EX_ser_L_LPAREN_e_RPAREN_"]=3
conf_HeLaMean["EX_thr_L_LPAREN_e_RPAREN_"]=3
conf_HeLaMean["EX_trp_L_LPAREN_e_RPAREN_"]=3
conf_HeLaMean["EX_tyr_L_LPAREN_e_RPAREN_"]=3
conf_HeLaMean["EX_thm_LPAREN_e_RPAREN_"]=3
conf_HeLaMean["EX_val_L_LPAREN_e_RPAREN_"]=3
conf_HeLaMean["EX_glc_LPAREN_e_RPAREN_"]=3
conf_HeLaMean["EX_pyr_LPAREN_e_RPAREN_"]=3
conf_HeLaMean["EX_ca2_LPAREN_e_RPAREN_"]=3
conf_HeLaMean["EX_cl_LPAREN_e_RPAREN_"]=3
conf_HeLaMean["EX_k_LPAREN_e_RPAREN_"]=3
conf_HeLaMean["EX_na1_LPAREN_e_RPAREN_"]=3
conf_HeLaMean["EX_fe2_LPAREN_e_RPAREN_"]=3
conf_HeLaMean["EX_fe3_LPAREN_e_RPAREN_"]=3
conf_HeLaMean["EX_so4_LPAREN_e_RPAREN_"]=3
conf_HeLaMean["EX_pi_LPAREN_e_RPAREN_"]=3
conf_HeLaMean["EX_h_LPAREN_e_RPAREN_"]=3
conf_HeLaMean["EX_h2o_LPAREN_e_RPAREN_"]=3
conf_HeLaMean["EX_o2_LPAREN_e_RPAREN_"]=3
conf_HeLaMean["EX_co2_LPAREN_e_RPAREN_"]=3
conf_HeLaMean["EX_lac_L_LPAREN_e_RPAREN_"]=3
conf_HeLaMean["EX_lac_D_LPAREN_e_RPAREN_"]=-1
conf_HeLaMean["EX_pro_L_LPAREN_e_RPAREN_"]=3
conf_HeLaMean["EX_ala_L_LPAREN_e_RPAREN_"]=3
conf_HeLaMean["EX_asn_L_LPAREN_e_RPAREN_"]=3
conf_HeLaMean["EX_asp_L_LPAREN_e_RPAREN_"]=3
conf_HeLaMean["EX_glu_L_LPAREN_e_RPAREN_"]=3
conf_HeLaMean["EX_chol_LPAREN_e_RPAREN_"]=3
conf_HeLaMean["EX_fol_LPAREN_e_RPAREN_"]=3
conf_HeLaMean["EX_inost_LPAREN_e_RPAREN_"]=3
conf_HeLaMean["EX_ncam_LPAREN_e_RPAREN_"]=3
conf_HeLaMean["EX_pnto_R_LPAREN_e_RPAREN_"]=3
conf_HeLaMean["EX_pydxn_LPAREN_e_RPAREN_"]=3
conf_HeLaMean["EX_ribflv_LPAREN_e_RPAREN_"]=3


conf_keratinocytesMean["EX_arg_L_LPAREN_e_RPAREN_"]=3
conf_keratinocytesMean["EX_cys_L_LPAREN_e_RPAREN_"]=3
conf_keratinocytesMean["EX_gln_L_LPAREN_e_RPAREN_"]=3
conf_keratinocytesMean["EX_gly_LPAREN_e_RPAREN_"]=3
conf_keratinocytesMean["EX_his_L_LPAREN_e_RPAREN_"]=3
conf_keratinocytesMean["EX_ile_L_LPAREN_e_RPAREN_"]=3
conf_keratinocytesMean["EX_leu_L_LPAREN_e_RPAREN_"]=3
conf_keratinocytesMean["EX_lys_L_LPAREN_e_RPAREN_"]=3
conf_keratinocytesMean["EX_met_L_LPAREN_e_RPAREN_"]=3
conf_keratinocytesMean["EX_phe_L_LPAREN_e_RPAREN_"]=3
conf_keratinocytesMean["EX_ser_L_LPAREN_e_RPAREN_"]=3
conf_keratinocytesMean["EX_thr_L_LPAREN_e_RPAREN_"]=3
conf_keratinocytesMean["EX_trp_L_LPAREN_e_RPAREN_"]=3
conf_keratinocytesMean["EX_tyr_L_LPAREN_e_RPAREN_"]=3
conf_keratinocytesMean["EX_thm_LPAREN_e_RPAREN_"]=3
conf_keratinocytesMean["EX_val_L_LPAREN_e_RPAREN_"]=3
conf_keratinocytesMean["EX_glc_LPAREN_e_RPAREN_"]=3
conf_keratinocytesMean["EX_pyr_LPAREN_e_RPAREN_"]=3
conf_keratinocytesMean["EX_ca2_LPAREN_e_RPAREN_"]=3
conf_keratinocytesMean["EX_cl_LPAREN_e_RPAREN_"]=3
conf_keratinocytesMean["EX_k_LPAREN_e_RPAREN_"]=3
conf_keratinocytesMean["EX_na1_LPAREN_e_RPAREN_"]=3
conf_keratinocytesMean["EX_fe2_LPAREN_e_RPAREN_"]=3
conf_keratinocytesMean["EX_fe3_LPAREN_e_RPAREN_"]=3
conf_keratinocytesMean["EX_so4_LPAREN_e_RPAREN_"]=3
conf_keratinocytesMean["EX_pi_LPAREN_e_RPAREN_"]=3
conf_keratinocytesMean["EX_h_LPAREN_e_RPAREN_"]=3
conf_keratinocytesMean["EX_h2o_LPAREN_e_RPAREN_"]=3
conf_keratinocytesMean["EX_o2_LPAREN_e_RPAREN_"]=3
conf_keratinocytesMean["EX_co2_LPAREN_e_RPAREN_"]=3
conf_keratinocytesMean["EX_lac_L_LPAREN_e_RPAREN_"]=3
conf_keratinocytesMean["EX_lac_D_LPAREN_e_RPAREN_"]=-1
conf_keratinocytesMean["EX_pro_L_LPAREN_e_RPAREN_"]=3
conf_keratinocytesMean["EX_ala_L_LPAREN_e_RPAREN_"]=3
conf_keratinocytesMean["EX_asn_L_LPAREN_e_RPAREN_"]=3
conf_keratinocytesMean["EX_asp_L_LPAREN_e_RPAREN_"]=3
conf_keratinocytesMean["EX_glu_L_LPAREN_e_RPAREN_"]=3
conf_keratinocytesMean["EX_chol_LPAREN_e_RPAREN_"]=3
conf_keratinocytesMean["EX_fol_LPAREN_e_RPAREN_"]=3
conf_keratinocytesMean["EX_inost_LPAREN_e_RPAREN_"]=3
conf_keratinocytesMean["EX_ncam_LPAREN_e_RPAREN_"]=3
conf_keratinocytesMean["EX_pnto_R_LPAREN_e_RPAREN_"]=3
conf_keratinocytesMean["EX_pydxn_LPAREN_e_RPAREN_"]=3
conf_keratinocytesMean["EX_ribflv_LPAREN_e_RPAREN_"]=3


```


```python
%matplotlib inline
import pandas as pd

df=pd.DataFrame({'HeLa': conf_HeLaMean})
df.HeLa.value_counts()

```
     1.0    3257
     3.0    1744
    -1.0    1547
     0.0     677
     2.0     560
    Name: HeLa, dtype: int64


```python
%matplotlib inline
import pandas as pd

df=pd.DataFrame({'Keratinocytes': conf_keratinocytesMean})
df.Keratinocytes.value_counts()

```

     1.0    3373
     3.0    2021
    -1.0    1724
     0.0     440
     2.0     227
    Name: Keratinocytes, dtype: int64

```python
conf_HeLaMean["biomass_reaction"]=3
conf_HeLaMean["biomass_DNA"]=3
conf_HeLaMean["biomass_RNA"]=3
conf_HeLaMean["biomass_carbohydrate"]=3
conf_HeLaMean["biomass_lipid"]=3
conf_HeLaMean["biomass_other"]=3
conf_HeLaMean["biomass_protein"]=3
conf_HeLaMean["DM_atp_c_"]=3

conf_keratinocytesMean["biomass_reaction"]=3
conf_keratinocytesMean["biomass_DNA"]=3
conf_keratinocytesMean["biomass_RNA"]=3
conf_keratinocytesMean["biomass_carbohydrate"]=3
conf_keratinocytesMean["biomass_lipid"]=3
conf_keratinocytesMean["biomass_other"]=3
conf_keratinocytesMean["biomass_protein"]=3
conf_keratinocytesMean["DM_atp_c_"]=3

Recon2.objective="biomass_reaction"
```


```python
metas = ['0.014 biomass_DNA_c + 0.058 biomass_RNA_c + 0.071 biomass_carbohydrate_c + 0.097 biomass_lipid_c + 0.054 biomass_other_c + 0.706 biomass_protein_c --> ', '3pg_c', '4abut_c', '4hpro_LT_c', 'accoa_m', 'accoa_m -> coa_m', 'ade_c', 'adn_c', 'adp_c', 'akg_m', 'ala_B_c', 'ala_L_c', 'amet_c', 'amp_c', 'arg_L_c', 'asn_L_c', 'asp_D_c', 'asp_L_c', 'atp_c', 'bhb_c', 'cdp_c', 'CE1936_c', 'chol_c', 'chsterol_c', 'cit_c', 'citr_L_c', 'cmp_c', 'creat_c', 'crm_hs_c', 'crtn_c', 'ctp_c', 'cys_L_c', 'dag_hs_c', 'dhap_c', 'e4p_c', 'f6p_c', 'fdp_c', 'fru_c', 'fum_c', 'g1p_c', 'g3p_c', 'g6p_c', 'gdp_c', 'glc_D_c', 'glc_D_e', 'g6p_c',  'f6p_c','glc_D_e -> glc_D_c', 'gln_L_c','gln_L_e', 'gln_L_e -> gln_L_c', 'glu_L_c', 'glu_L_e', 'glyb_c', 'gly_c', 'gmp_c', 'gthox_c', 'gthrd_c', 'gua_c', 'HC00342_c', 'his_L_c', 'hxan_c', 'icit_c', 'ile_L_c', 'lac_L_c', 'leu_L_c', 'leu_L_c',  'lys_L_c', 'mag_hs_c', 'mal_L_c', 'met_L_c', 'nad_c', 'nadh_c', 'nadh_m', 'nad_m', 'nadp_c', 'nadph_c', 'nadph_m', 'nadph_m', 'nadp_m', 'oaa_m', 'orn_c', 'pa_hs_c', 'pe_hs_c', 'pep_c', 'phe_L_c', 'pmtcoa_c -> coa_c', 'pro_D_c', 'pro_L_c', 'ps_hs_c', 'ptrc_c', 'pyr_c', 'pyr_m', 'r5p_c', 'ru5p_D_c', 's7p_c', 'ser_L_c', 'spmd_c', 'succ_c', 'succoa_m -> coa_m', 'tag_hs_c', 'thr_L_c', 'trp_L_c', 'tym_c', 'tyr_L_c', 'udp_c', 'ump_c', 'utp_c', 'val_L_c']
```

```python
%%time


from corda import CORDA

opt_HeLaMean = CORDA(model=Recon2, confidence=conf_HeLaMean, n=5,met_prod=metas,  penalty_factor=1000) 
opt_HeLaMean.build()
print(opt_HeLaMean)
model_HeLaMean=opt_HeLaMean.cobra_model(name="HeLaMean")
print(model_HeLaMean.optimize())

```

    build status: reconstruction complete
    Inc. reactions: 3006/7894
     - unclear: 204/677
     - exclude: 232/1547
     - low and medium: 896/3809
     - high: 1674/1861
    
    <Solution 0.027 at 0x7f188e592d30>
    CPU times: user 22min 55s, sys: 884 ms, total: 22min 56s
    Wall time: 22min 56s



```python
cp = model_HeLaMean.copy()
cp.optimize()
cp.summary(fva=True)
FBA=cp.optimize()
FBA.fluxes.to_csv("HeLaPrueba.tsv")
```

    IN FLUXES                                        OUT FLUXES                                       OBJECTIVES
    -----------------------------------------------  -----------------------------------------------  -----------------------
    id                   Flux  Range                 id                   Flux  Range                 biomass_reac...  0.0266
    ---------------  --------  --------------------  ---------------  --------  --------------------
    glc_D_e          0.306     [0.26, 4.5]           lac_L_e          0.601     [-1, 9.62]
    ser_L_e          0.023     [0.042, -1.22]        h_e              0.597     [-0.715, 10]
    o2_e             0.0172    [0.00272, 10]         co2_e            0.0254    [-0.001, 6.6]
    lys_L_e          0.0158    [0.0158, 0.0168]      glyc_e           0.0107    [-0.001, 1.26]
    leu_L_e          0.0145    [0.0145, 0.0155]      urate_e          0.00375   [0, 0.411]
    arg_L_e          0.0112    [-0.0263, 0.084]      ura_e            0.00375   [-0.001, 0.00629]
    glu_L_e          0.01      [-0.001, 0.01]        ade_e            0.00357   [-0.001, 0.32]
    pro_L_e          0.01      [0.01, -0.0945]       ac_e             0.00308   [-0.001, 0.0635]
    val_L_e          0.00938   [0.00938, 0.094]      vacc_e           0.002     [-0.001, 0.505]
    asp_L_e          0.00938   [0, 0.01]             ptrc_e           0.00168   [-0.001, 0.181]
    gly_e            0.00878   [0.03, -1.69]         urea_e           0.00168   [0, 0.11]
    gln_L_e          0.00867   [-0.277, 0.584]       for_e            0.00154   [-0.001, 0.252]
    thr_L_e          0.00832   [0.00832, 0.00832]    hdcea_e          0.001     [-0.001, 0.569]
    ile_L_e          0.00761   [0.00761, 0.105]      cholate_e        0.001     [0, 0.219]
    asn_L_e          0.00743   [0.01, -0.851]        fol_c            0.001     [-0.001, 0.052]
    phe_L_e          0.0069    [0.0069, 0.066]       thym_e           0.001     [-0.001, 0.00861]
    pi_e             0.00562   [-0.000892, 1]        1mncam_e         0.001     [0, 0.002]
    tyr_L_e          0.00425   [-0.0549, 0.104]      ribflv_e         0.001     [0, 0.002]
    his_L_e          0.00336   [0.00236, 0.00336]    leuktrD4_e       0.001     [-0.001, 0.002]
    chol_e           0.00305   [0.00305, 0.05]       h2o_e            0.000747  [-5.26, 10]
    met_L_e          0.00207   [0.00207, 0.03]       dttp_m           0.000652  [0, 0.00994]
    pyr_e            0.00155   [0.11, -10]           so4_e            0.000618  [0.0933, -0.164]
    biomass_other_c  0.00144   [0.00144, 0.00144]    3mlda_e          0.000419  [0, 0.001]
    dag_hs_e         0.001     [0.001, 0.001]        dgpi_prot_hs_r   0.000333  [0.000333, 0.000333]
    lpchol_hs_e      0.001     [0.001, 0.001]        gpi_sig_r        0.000333  [0.000333, 0.000333]
    pchol_hs_e       0.001     [0.001, 0.001]        lnlc_e           0.00025   [-0.001, 0.508]
    pe_hs_e          0.001     [0.001, 0.001]        cys_L_e          0.000144  [0.0287, -0.0626]
    pe_hs_r          0.001     [0.001, 0.001]        h2o2_e           0         [-0.001, 10]
    ps_hs_e          0.001     [0.001, 0.001]        gal_e            0         [0, 2.12]
    tag_hs_e         0.001     [0.001, 0.001]        abt_e            0         [0, 1.7]
    4hpro_LT_m       0.001     [0, 0.001]            oxa_e            0         [0, 1.58]
    5mthf_e          0.001     [0, 0.001]            ala_L_e          0         [-0.01, 1.35]
    Lcystin_e        0.001     [0, 0.001]            glyc_S_e         0         [0, 1.2]
    arab_L_e         0.001     [0, 0.001]            fe2_e            0         [-1, 1]
    dcsptn1_e        0.001     [0, 0.001]            fe3_e            0         [-1, 1]
    fad_e            0.001     [0, 0.001]            5oxpro_e         0         [0, 0.936]
    gtp_e            0.001     [0, 0.001]            ppi_e            0         [0, 0.5]
    o2s_e            0.001     [0, 0.001]            nrvnc_e          0         [0, 0.367]
    strdnc_e         0.001     [0, 0.001]            gua_e            0         [-0.001, 0.354]
    thymd_e          0.001     [0, 0.001]            chsterol_e       0         [0, 0.225]
    gdp_e            0.001     [0.001, -0.001]       C02528_e         0         [0, 0.217]
    ncam_c           0.001     [0.001, -0.001]       fuc_L_e          0         [0, 0.188]
    gmp_e            0.001     [0.001, -0.002]       xolest2_hs_e     0         [-0.001, 0.169]
    leuktrC4_e       0.001     [0.001, -0.002]       4hphac_e         0         [0, 0.159]
    cytd_e           0.001     [0.001, -0.00629]     tymsf_e          0         [0, 0.159]
    dctp_n           0.001     [0.001, -0.00629]     fucfucgalacg...  0         [0, 0.152]
    udp_e            0.001     [0.001, -0.00629]     orn_e            0         [-0.001, 0.105]
    ump_e            0.001     [0.001, -0.00629]     galgalfucfuc...  0         [0, 0.096]
    uri_e            0.001     [0.001, -0.00629]     glyb_e           0         [0, 0.047]
    utp_e            0.001     [0.001, -0.00629]     gthox_e          0         [0, 0.0456]
    ahcys_e          0.001     [0.001, -0.00815]     citr_L_c         0         [-0.001, 0.0348]
    atp_e            0.001     [0.001, -0.00815]     2hb_e            0         [0, 0.0279]
    gsn_e            0.001     [0.001, -0.00815]     q10h2_e          0         [0, 0.0261]
    dttp_n           0.001     [0.001, -0.00894]     anth_c           0         [0, 0.0156]
    5mta_e           0.001     [0.001, -0.00915]     bildglcur_e      0         [0, 0.0151]
    dad_2_e          0.001     [0.001, -0.0142]      din_e            0         [-0.001, 0.0142]
    Rtotal_e         0.001     [0.001, -0.172]       bilglcur_e       0         [-0.001, 0.0141]
    gchola_e         0.001     [0.001, -0.216]       co_e             0         [0, 0.0131]
    ocdcea_e         0.001     [0.001, -0.488]       pheme_e          0         [0, 0.0131]
    man_e            0.001     [0.001, -1.21]        octdececoa_c     0         [-0.000333, 0.0102]
    xylt_e           0.001     [0.001, -1.7]         35cgmp_e         0         [0, 0.00915]
    hco3_e           0.001     [0.001, -3.96]        camp_e           0         [0, 0.00915]
    inost_e          0.000953  [0.000953, -2.83]     spmd_e           0         [-0.001, 0.00915]
    ins_e            0.000712  [0.001, -0.00815]     ala_B_e          0         [0, 0.00829]
    whhdca_e         0.000611  [0, 0.001]            adn_e            0         [-0.001, 0.00815]
    sph1p_e          0.000465  [0.001, -0.498]       dcmp_e           0         [-0.001, 0.00629]
    orot_e           0.00046   [0, 0.001]            dcyt_e           0         [-0.001, 0.00629]
    hista_e          0.000419  [0, 0.001]            sprm_c           0         [0, 0.00508]
    pglyc_hs_e       0.000388  [0.000388, 0.000388]  Ser_Thr_l        0         [-0.001, 0.005]
    trp_L_e          0.000354  [0.000354, 0.016]     adrn_e           0         [-0.001, 0.003]
    datp_n           0.000351  [0.000351, 0.000351]  5adtststeron...  0         [0, 0.002]
    pre_prot_r       0.000333  [0.000333, 0.000333]  andrstrnglc_e    0         [0, 0.002]
    dgtp_n           0.000263  [0.000263, 0.000263]  estroneglc_e     0         [0, 0.002]
    ha_e             0.000167  [0.001, -0.417]       imp_e            0         [0, 0.002]
    mag_hs_e         0.000134  [0, 0.001]            n5m2masn_g       0         [0, 0.002]
    acac_e           0.000117  [0, 0.001]            nac_e            0         [0, 0.002]
    fol_e            0         [-0.003, 0.05]        arachd_e         0         [-0.001, 0.002]
    pnto_R_e         0         [0, 0.00915]          core8_g          0         [0, 0.0012]
    3hpvs_e          0         [0, 0.001]            3bcrn_e          0         [0, 0.001]
    4nph_e           0         [0, 0.001]            3hdececrn_e      0         [0, 0.001]
    5fthf_e          0         [0, 0.001]            3hpvstet_e       0         [0, 0.001]
    HC00822_e        0         [0, 0.001]            3ivcrn_e         0         [0, 0.001]
    acmp_e           0         [0, 0.001]            3octdec2crn_e    0         [0, 0.001]
    bilirub_e        0         [0, 0.001]            3octdeccrn_e     0         [0, 0.001]
    carn_e           0         [0, 0.001]            3octdece1crn_e   0         [0, 0.001]
    crn_e            0         [0, 0.001]            3tetd7ecoacrn_e  0         [0, 0.001]
    crvnc_e          0         [0, 0.001]            3thexddcoacrn_e  0         [0, 0.001]
    datp_m           0         [0, 0.001]            3ttetddcoacrn_e  0         [0, 0.001]
    dca_e            0         [0, 0.001]            4nphsf_e         0         [0, 0.001]
    dheas_e          0         [0, 0.001]            aprgstrn_e       0         [0, 0.001]
    doco13ac_e       0         [0, 0.001]            c4crn_e          0         [0, 0.001]
    dopa_e           0         [0, 0.001]            c5dc_e           0         [0, 0.001]
    etoh_e           0         [0, 0.001]            c6crn_e          0         [0, 0.001]
    fmn_e            0         [0, 0.001]            clpnd_e          0         [0, 0.001]
    gam_e            0         [0, 0.001]            dmhptcrn_e       0         [0, 0.001]
    gluala_e         0         [0, 0.001]            dopasf_e         0         [0, 0.001]
    hdca_e           0         [0, 0.001]            epoxtac_e        0         [0, 0.001]
    hpdca_e          0         [0, 0.001]            ivcrn_e          0         [0, 0.001]
    leuktrA4_e       0         [0, 0.001]            lac_D_e          0         [0, 0.001]
    lnlnca_e         0         [0, 0.001]            leuktrB4_e       0         [0, 0.001]
    lnlncg_e         0         [0, 0.001]            prostgf2_e       0         [0, 0.001]
    n2m2nmasn_e      0         [0, 0.001]            retn_e           0         [0, 0.001]
    ocdca_e          0         [0, 0.001]            sulpacmp_e       0         [0, 0.001]
    octa_e           0         [0, 0.001]            triodthysuf_e    0         [0, 0.001]
    phyt_e           0         [0, 0.001]            tststeroneglc_e  0         [0, 0.001]
    prgstrn_e        0         [0, 0.001]            5adtststerone_e  0         [-0.001, 0.001]
    prostge2_e       0         [0, 0.001]            Asn_X_Ser_Thr_l  0         [-0.001, 0.001]
    ptdca_e          0         [0, 0.001]            andrstrn_e       0         [-0.001, 0.001]
    retinol_e        0         [0, 0.001]            dlnlcg_e         0         [-0.001, 0.001]
    tacr_e           0         [0, 0.001]            estradiol_e      0         [-0.001, 0.001]
    tmndnc_e         0         [0, 0.001]            idp_e            0         [-0.001, 0.001]
    triodthy_e       0         [0, 0.001]            itp_e            0         [-0.001, 0.001]
    tststerone_e     0         [0, 0.001]            meoh_e           0         [-0.001, 0.001]
    ttdca_e          0         [0, 0.001]            taur_c           0         [-0.001, 0.001]
    T_antigen_g      0         [-0.0002, 0.001]      tchola_e         0         [-0.001, 0.001]
    pmtcoa_r         0         [-0.000333, 0.001]    C02470_e         0         [0, -0.001]
    eicostet_e       0         [-0.001, 0.001]       arach_e          0         [0, -0.001]
    estrones_e       0         [-0.001, 0.001]
    fald_e           0         [-0.001, 0.001]



```python
cp = model_HeLaMean.copy()
cp.reactions.EX_ser_L_LPAREN_e_RPAREN_.bounds=(0,0)
cp.reactions.EX_gly_LPAREN_e_RPAREN_.bounds=(0,0)
cp.optimize()
cp.summary(fva=True)
FBA=cp.optimize()
FBA.fluxes.to_csv("HeLaPrueba2.tsv")
```

    IN FLUXES                                        OUT FLUXES                                       OBJECTIVES
    -----------------------------------------------  -----------------------------------------------  -----------------------
    id                   Flux  Range                 id                   Flux  Range                 biomass_reac...  0.0266
    ---------------  --------  --------------------  ---------------  --------  --------------------
    glc_D_e          0.357     [0.273, 4.5]          h_e              0.702     [-0.657, 10]
    o2_e             0.0337    [0.00272, 10]         lac_L_e          0.699     [-1, 9.56]
    pi_e             0.0168    [-0.000892, 1]        co2_e            0.0271    [-0.001, 6.56]
    lys_L_e          0.0158    [0.0158, 0.0168]      fol_c            0.0111    [-0.001, 0.052]
    arg_L_e          0.0152    [-0.0263, 0.084]      urate_e          0.00761   [0, 0.399]
    gln_L_e          0.0148    [-0.241, 0.584]       so4_e            0.00659   [0.0933, -0.164]
    leu_L_e          0.0145    [0.0145, 0.0155]      ac_e             0.00657   [-0.001, 0.0635]
    fol_e            0.0111    [-0.003, 0.05]        ptrc_e           0.00568   [-0.001, 0.181]
    asp_L_e          0.01      [0, 0.01]             urea_e           0.00568   [0, 0.11]
    glu_L_e          0.01      [-0.001, 0.01]        thym_e           0.00423   [-0.001, 0.00861]
    pro_L_e          0.01      [0.01, -0.0945]       dttp_m           0.00338   [0, 0.00994]
    val_L_e          0.00938   [0.00938, 0.094]      2hb_e            0.00235   [0, 0.0279]
    thr_L_e          0.00832   [0.00832, 0.00832]    vacc_e           0.00233   [-0.001, 0.504]
    ala_L_e          0.00781   [0.01, -1.31]         h2o_e            0.00229   [-5.14, 10]
    ile_L_e          0.00761   [0.00761, 0.105]      hdcea_e          0.00133   [-0.001, 0.568]
    asn_L_e          0.00743   [0.01, -0.815]        cholate_e        0.001     [0, 0.219]
    phe_L_e          0.0069    [0.0069, 0.066]       1mncam_e         0.001     [0, 0.002]
    met_L_e          0.00542   [0.00207, 0.03]       ribflv_e         0.001     [0, 0.002]
    tyr_L_e          0.0045    [-0.0549, 0.104]      leuktrD4_e       0.001     [-0.001, 0.002]
    cys_L_e          0.00348   [-0.0287, 0.0626]     3mlda_e          0.001     [0, 0.001]
    his_L_e          0.00336   [0.00236, 0.00336]    for_e            0.000543  [-0.001, 0.252]
    chol_e           0.00305   [0.00305, 0.05]       dgpi_prot_hs_r   0.000333  [0.000333, 0.000333]
    biomass_other_c  0.00144   [0.00144, 0.00144]    gpi_sig_r        0.000333  [0.000333, 0.000333]
    dag_hs_e         0.001     [0.001, 0.001]        pmtcoa_r         0.000329  [0.000333, -0.001]
    lpchol_hs_e      0.001     [0.001, 0.001]        q10h2_e          0.000256  [0, 0.0261]
    pchol_hs_e       0.001     [0.001, 0.001]        lnlc_e           0.00025   [-0.001, 0.507]
    pe_hs_e          0.001     [0.001, 0.001]        pyr_e            0         [-0.11, 10]
    pe_hs_r          0.001     [0.001, 0.001]        h2o2_e           0         [-0.001, 10]
    ps_hs_e          0.001     [0.001, 0.001]        gal_e            0         [0, 2.11]
    tag_hs_e         0.001     [0.001, 0.001]        abt_e            0         [0, 1.69]
    4hpro_LT_m       0.001     [0, 0.001]            oxa_e            0         [0, 1.54]
    Lcystin_e        0.001     [0, 0.001]            glyc_S_e         0         [0, 1.1]
    acac_e           0.001     [0, 0.001]            fe2_e            0         [-1, 1]
    arab_L_e         0.001     [0, 0.001]            fe3_e            0         [-1, 1]
    dcsptn1_e        0.001     [0, 0.001]            5oxpro_e         0         [0, 0.932]
    fad_e            0.001     [0, 0.001]            ppi_e            0         [0, 0.5]
    hista_e          0.001     [0, 0.001]            nrvnc_e          0         [0, 0.366]
    o2s_e            0.001     [0, 0.001]            gua_e            0         [-0.001, 0.339]
    orot_e           0.001     [0, 0.001]            ade_e            0         [-0.001, 0.31]
    strdnc_e         0.001     [0, 0.001]            chsterol_e       0         [0, 0.224]
    thymd_e          0.001     [0, 0.001]            C02528_e         0         [0, 0.217]
    gdp_e            0.001     [0.001, -0.001]       fuc_L_e          0         [0, 0.188]
    ncam_c           0.001     [0.001, -0.001]       xolest2_hs_e     0         [-0.001, 0.169]
    gmp_e            0.001     [0.001, -0.002]       4hphac_e         0         [0, 0.159]
    leuktrC4_e       0.001     [0.001, -0.002]       tymsf_e          0         [0, 0.159]
    cytd_e           0.001     [0.001, -0.00629]     fucfucgalacg...  0         [0, 0.151]
    dctp_n           0.001     [0.001, -0.00629]     orn_e            0         [-0.001, 0.105]
    udp_e            0.001     [0.001, -0.00629]     galgalfucfuc...  0         [0, 0.0957]
    ump_e            0.001     [0.001, -0.00629]     glyb_e           0         [0, 0.047]
    ura_e            0.001     [0.001, -0.00629]     gthox_e          0         [0, 0.0456]
    uri_e            0.001     [0.001, -0.00629]     citr_L_c         0         [-0.001, 0.0348]
    utp_e            0.001     [0.001, -0.00629]     anth_c           0         [0, 0.0156]
    adn_e            0.001     [0.001, -0.00815]     bildglcur_e      0         [0, 0.0151]
    ahcys_e          0.001     [0.001, -0.00815]     bilglcur_e       0         [-0.001, 0.0141]
    atp_e            0.001     [0.001, -0.00815]     co_e             0         [0, 0.0131]
    gsn_e            0.001     [0.001, -0.00815]     pheme_e          0         [0, 0.0131]
    ins_e            0.001     [0.001, -0.00815]     35cgmp_e         0         [0, 0.00915]
    dttp_n           0.001     [0.001, -0.00894]     camp_e           0         [0, 0.00915]
    dad_2_e          0.001     [0.001, -0.0142]      5mta_e           0         [-0.001, 0.00915]
    din_e            0.001     [0.001, -0.0142]      spmd_e           0         [-0.001, 0.00915]
    Rtotal_e         0.001     [0.001, -0.171]       ala_B_e          0         [0, 0.00829]
    gchola_e         0.001     [0.001, -0.214]       dcmp_e           0         [-0.001, 0.00629]
    ocdcea_e         0.001     [0.001, -0.487]       sprm_c           0         [0, 0.00508]
    man_e            0.001     [0.001, -1.21]        Ser_Thr_l        0         [-0.001, 0.005]
    xylt_e           0.001     [0.001, -1.69]        adrn_e           0         [-0.001, 0.003]
    hco3_e           0.001     [0.001, -3.93]        5adtststeron...  0         [0, 0.002]
    inost_e          0.000953  [0.000953, -2.82]     andrstrnglc_e    0         [0, 0.002]
    crvnc_e          0.000863  [0, 0.001]            estroneglc_e     0         [0, 0.002]
    glyc_e           0.000768  [0.001, -1.2]         imp_e            0         [0, 0.002]
    dcyt_e           0.00067   [0.001, -0.00629]     n5m2masn_g       0         [0, 0.002]
    sph1p_e          0.000465  [0.001, -0.488]       nac_e            0         [0, 0.002]
    pglyc_hs_e       0.000388  [0.000388, 0.000388]  arachd_e         0         [-0.001, 0.002]
    trp_L_e          0.000354  [0.000354, 0.016]     core8_g          0         [0, 0.0012]
    datp_n           0.000351  [0.000351, 0.000351]  3bcrn_e          0         [0, 0.001]
    pre_prot_r       0.000333  [0.000333, 0.000333]  3hdececrn_e      0         [0, 0.001]
    octdececoa_c     0.000329  [0.000333, -0.0102]   3hpvstet_e       0         [0, 0.001]
    dgtp_n           0.000263  [0.000263, 0.000263]  3ivcrn_e         0         [0, 0.001]
    ha_e             0.000167  [0.001, -0.416]       3octdec2crn_e    0         [0, 0.001]
    pnto_R_e         0         [0, 0.00915]          3octdeccrn_e     0         [0, 0.001]
    3hpvs_e          0         [0, 0.001]            3octdece1crn_e   0         [0, 0.001]
    4nph_e           0         [0, 0.001]            3tetd7ecoacrn_e  0         [0, 0.001]
    5fthf_e          0         [0, 0.001]            3thexddcoacrn_e  0         [0, 0.001]
    5mthf_e          0         [0, 0.001]            3ttetddcoacrn_e  0         [0, 0.001]
    C02470_e         0         [0, 0.001]            4nphsf_e         0         [0, 0.001]
    HC00822_e        0         [0, 0.001]            aprgstrn_e       0         [0, 0.001]
    acmp_e           0         [0, 0.001]            c4crn_e          0         [0, 0.001]
    arach_e          0         [0, 0.001]            c5dc_e           0         [0, 0.001]
    bilirub_e        0         [0, 0.001]            c6crn_e          0         [0, 0.001]
    carn_e           0         [0, 0.001]            clpnd_e          0         [0, 0.001]
    crn_e            0         [0, 0.001]            dmhptcrn_e       0         [0, 0.001]
    datp_m           0         [0, 0.001]            dopasf_e         0         [0, 0.001]
    dca_e            0         [0, 0.001]            epoxtac_e        0         [0, 0.001]
    dheas_e          0         [0, 0.001]            ivcrn_e          0         [0, 0.001]
    doco13ac_e       0         [0, 0.001]            lac_D_e          0         [0, 0.001]
    dopa_e           0         [0, 0.001]            leuktrB4_e       0         [0, 0.001]
    etoh_e           0         [0, 0.001]            prostgf2_e       0         [0, 0.001]
    fmn_e            0         [0, 0.001]            retn_e           0         [0, 0.001]
    gam_e            0         [0, 0.001]            sulpacmp_e       0         [0, 0.001]
    gluala_e         0         [0, 0.001]            triodthysuf_e    0         [0, 0.001]
    gtp_e            0         [0, 0.001]            tststeroneglc_e  0         [0, 0.001]
    hdca_e           0         [0, 0.001]            5adtststerone_e  0         [-0.001, 0.001]
    hpdca_e          0         [0, 0.001]            Asn_X_Ser_Thr_l  0         [-0.001, 0.001]
    leuktrA4_e       0         [0, 0.001]            andrstrn_e       0         [-0.001, 0.001]
    lnlnca_e         0         [0, 0.001]            dlnlcg_e         0         [-0.001, 0.001]
    lnlncg_e         0         [0, 0.001]            eicostet_e       0         [-0.001, 0.001]
    mag_hs_e         0         [0, 0.001]            estradiol_e      0         [-0.001, 0.001]
    n2m2nmasn_e      0         [0, 0.001]            idp_e            0         [-0.001, 0.001]
    ocdca_e          0         [0, 0.001]            itp_e            0         [-0.001, 0.001]
    octa_e           0         [0, 0.001]            meoh_e           0         [-0.001, 0.001]
    phyt_e           0         [0, 0.001]            taur_c           0         [-0.001, 0.001]
    prgstrn_e        0         [0, 0.001]            tchola_e         0         [-0.001, 0.001]
    prostge2_e       0         [0, 0.001]
    ptdca_e          0         [0, 0.001]
    retinol_e        0         [0, 0.001]
    tacr_e           0         [0, 0.001]
    tmndnc_e         0         [0, 0.001]
    triodthy_e       0         [0, 0.001]
    tststerone_e     0         [0, 0.001]
    ttdca_e          0         [0, 0.001]
    whhdca_e         0         [0, 0.001]
    T_antigen_g      0         [-0.0002, 0.001]
    estrones_e       0         [-0.001, 0.001]
    fald_e           0         [-0.001, 0.001]



```python
FBA.metabolites.hpyr_c.summary()
```

    PRODUCING REACTIONS -- 3-hydroxypyruvate (hpyr_c)
    -------------------------------------------------
    %       FLUX  RXN ID    REACTION
    ----  ------  --------  ------------------------------------
    100%   0.333  r0160     pyr_c + ser_L_c --> ala_L_c + hpyr_c
    
    CONSUMING REACTIONS -- 3-hydroxypyruvate (hpyr_c)
    -------------------------------------------------
    %       FLUX  RXN ID    REACTION
    ----  ------  --------  ------------------------------------
    100%   0.333  HPYRDC    h_c + hpyr_c --> co2_c + gcald_c



```python
FBA.metabolites.gly_c.summary()
```

    PRODUCING REACTIONS -- glycine (gly_c)
    --------------------------------------
    %      FLUX  RXN ID      REACTION
    ---  ------  ----------  --------------------------------------------------
    84%  0.352   GLYtm       gly_c <=> gly_m
    16%  0.0663  GLYt4       gly_e + na1_e --> gly_c + na1_c
    
    CONSUMING REACTIONS -- glycine (gly_c)
    --------------------------------------
    %      FLUX  RXN ID      REACTION
    ---  ------  ----------  --------------------------------------------------
    42%  0.177   GHMT2r      ser_L_c + thf_c <=> gly_c + h2o_c + mlthf_c
    37%  0.154   r0629       cholcoa_c + gly_c <=> coa_c + gchola_c
    17%  0.0726  r1552       gly_e + pro_L_c <=> gly_c + pro_L_e
    3%   0.0143  biomass...  0.716189801699717 ala_L_c + 0.508866855524079 a...



```python
%%time

from corda import CORDA

opt_KeratinocytesMean = CORDA(model=Recon2, confidence=conf_keratinocytesMean, n=5, met_prod=metas,  penalty_factor=1000) 
opt_KeratinocytesMean.build()
print(opt_KeratinocytesMean)

model_KeratinocytesMean=opt_KeratinocytesMean.cobra_model(name="KeratinocytesMean")
print(model_KeratinocytesMean.optimize())
```

    build status: reconstruction complete
    Inc. reactions: 3340/7894
     - unclear: 109/440
     - exclude: 259/1724
     - low and medium: 1015/3592
     - high: 1957/2138
    
    <Solution 0.024 at 0x7f188ce6d668>
    CPU times: user 28min 1s, sys: 1.51 s, total: 28min 3s
    Wall time: 28min 3s



```python
cp = model_KeratinocytesMean.copy()

cp.optimize()
cp.summary(fva=True)

```

    IN FLUXES                                        OUT FLUXES                                     OBJECTIVES
    -----------------------------------------------  ---------------------------------------------  -----------------------
    id                   Flux  Range                 id                   Flux  Range               biomass_reac...  0.0242
    ---------------  --------  --------------------  ---------------  --------  ------------------
    glc_D_e          0.272     [0.203, 4.5]          h_e              0.54      [-0.619, 10]
    o2_e             0.0545    [0.00246, 8.68]       lac_L_e          0.506     [-1, 9.73]
    phe_L_e          0.0261    [0.00629, 0.066]      co2_e            0.0369    [-0.001, 10]
    ser_L_e          0.0159    [0.042, -1.12]        pyr_e            0.0283    [-0.11, 9.86]
    lys_L_e          0.0144    [0.0144, 0.0154]      h2o_e            0.0249    [-1.54, 10]
    leu_L_e          0.0132    [0.0132, 0.0142]      tyr_L_e          0.0159    [0.0558, -0.104]
    pro_L_e          0.01      [0.01, 0.01]          ala_B_e          0.00853   [-0.001, 1.31]
    glu_L_e          0.01      [0.01, -1.08]         pi_e             0.00363   [0.0162, -0.896]
    val_L_e          0.00896   [0.00855, 0.094]      dgsn_e           0.00321   [-0.001, 0.0188]
    arg_L_e          0.00871   [0.00671, 0.00871]    Rtotal_e         0.00192   [-0.001, 0.247]
    asp_L_e          0.00814   [0, 0.01]             gua_e            0.00167   [-0.001, 0.342]
    gln_L_e          0.0079    [-0.251, 0.584]       cit_e            0.00164   [-0.001, 1.09]
    thr_L_e          0.00758   [0.00758, 0.00758]    vacc_e           0.001     [-0.001, 0.461]
    ile_L_e          0.00694   [0.00694, 0.105]      cgly_e           0.001     [-0.001, 0.0898]
    asn_L_e          0.00677   [0.00677, -0.828]     3ivcrn_e         0.001     [0, 0.001]
    ala_L_e          0.00539   [0.01, -1.16]         epoxtac_e        0.001     [0, 0.001]
    met_L_e          0.00464   [0.00371, 0.03]       ndersv_e         0.001     [0, 0.001]
    gly_e            0.00419   [0.03, -1.14]         2hb_e            0.000927  [0, 0.0273]
    his_L_e          0.00307   [0.00206, 0.00307]    ethamp_r         0.000576  [0, 0.852]
    biomass_other_c  0.00131   [0.00131, 0.00131]    ade_e            0.000379  [-0.001, 0.347]
    4hpro_LT_m       0.001     [0, 0.001]            xolest2_hs_e     8.4e-05   [-0.001, 0.245]
    crn_e            0.001     [0, 0.001]            oxa_e            0         [0, 2.16]
    orot_e           0.001     [0, 0.001]            gal_e            0         [0, 2.15]
    pchol_hs_e       0.001     [0, 0.001]            glyc_S_e         0         [0, 1.6]
    rsv_e            0.001     [0, 0.001]            glyc_e           0         [-0.001, 1.6]
    tacr_e           0.001     [0, 0.001]            acmana_e         0         [0, 1.32]
    tag_hs_e         0.001     [0, 0.001]            taur_c           0         [-0.001, 1.09]
    dag_hs_e         0.001     [-0.000924, 0.001]    taur_e           0         [-0.001, 1.09]
    glyc3p_e         0.001     [-0.000924, 0.001]    5oxpro_e         0         [0, 1.09]
    lpchol_hs_e      0.001     [-0.000924, 0.001]    HC00342_e        0         [-0.001, 1.09]
    pe_hs_e          0.001     [-0.000924, 0.001]    fe2_e            0         [-1, 1]
    pe_hs_r          0.001     [-0.000924, 0.001]    ha_pre1_e        0         [-0.001, 0.848]
    ps_hs_e          0.001     [-0.000924, 0.001]    hdcea_e          0         [-0.001, 0.522]
    dcsptn1_e        0.001     [0.001, -0.001]       hdca_e           0         [-0.001, 0.511]
    cdp_e            0.001     [0.001, -0.002]       elaid_e          0         [-0.001, 0.461]
    cmp_e            0.001     [0.001, -0.002]       urate_e          0         [0, 0.429]
    ctp_e            0.001     [0.001, -0.002]       ha_e             0         [-0.001, 0.424]
    cytd_e           0.001     [0.001, -0.00753]     chsterol_e       0         [0, 0.391]
    dctp_n           0.001     [0.001, -0.00753]     cholate_e        0         [0, 0.378]
    udp_e            0.001     [0.001, -0.00753]     sph1p_e          0         [-0.001, 0.367]
    ump_e            0.001     [0.001, -0.00753]     tdchola_e        0         [0, 0.354]
    ura_e            0.001     [0.001, -0.00753]     eicostet_e       0         [-0.001, 0.352]
    uri_e            0.001     [0.001, -0.00753]     nrvnc_e          0         [0, 0.329]
    utp_e            0.001     [0.001, -0.00753]     4hphac_e         0         [0, 0.16]
    aicar_e          0.001     [0.001, -0.0136]      acac_e           0         [-0.001, 0.16]
    atp_e            0.001     [0.001, -0.0136]      glyb_e           0         [0, 0.128]
    ins_e            0.001     [0.001, -0.0136]      spc_hs_e         0         [0, 0.102]
    xmp_e            0.001     [0.001, -0.0136]      q10h2_e          0         [0, 0.0904]
    datp_n           0.001     [0.001, -0.0143]      q10_e            0         [-0.001, 0.089]
    dgtp_n           0.001     [0.001, -0.0188]      orn_e            0         [0, 0.0865]
    din_e            0.001     [0.001, -0.0188]      pheacgln_e       0         [0, 0.0597]
    gthrd_e          0.001     [0.001, -0.0898]      dad_2_e          0         [-0.001, 0.0188]
    sphs1p_e         0.001     [0.001, -0.367]       35cgmp_e         0         [0, 0.0146]
    strdnc_e         0.001     [0.001, -0.385]       camp_e           0         [0, 0.0146]
    hxan_e           0.001     [0.001, -0.428]       adn_e            0         [-0.001, 0.0136]
    xylt_e           0.001     [0.001, -1.72]        ahcys_e          0         [-0.001, 0.0136]
    inost_e          0.000565  [0.000565, 0.00156]   gsn_e            0         [-0.001, 0.0136]
    for_e            0.000422  [0.001, -0.395]       bilglcur_e       0         [0, 0.0108]
    pglyc_hs_e       0.000353  [0.000353, 0.000353]  co_e             0         [0, 0.0108]
    trp_L_e          0.000323  [0.000323, 0.000323]  pheme_e          0         [0, 0.0108]
    dttp_n           0.000317  [0.001, -0.00821]     dttp_m           0         [0, 0.00921]
    adrn_e           0.000219  [0, 0.001]            Ser_Thr_l        0         [-0.001, 0.005]
    cys_L_e          0.000202  [-0.0262, 0.0626]     fuc14galacgl...  0         [0, 0.004]
    so4_e            0         [-0.0898, 1]          fuc_L_e          0         [0, 0.004]
    fe3_e            0         [-1, 1]               galfuc12gal1...  0         [0, 0.004]
    3hpvs_e          0         [0, 0.001]            galfucgalacg...  0         [0, 0.004]
    5HPET_c          0         [0, 0.001]            fald_e           0         [-0.001, 0.003]
    5fthf_e          0         [0, 0.001]            man_e            0         [-0.001, 0.003]
    HC00822_e        0         [0, 0.001]            meoh_e           0         [-0.001, 0.003]
    acmp_e           0         [0, 0.001]            34dhoxpeg_e      0         [0, 0.002]
    allop_e          0         [0, 0.001]            5adtststeron...  0         [0, 0.002]
    amp_e            0         [0, 0.001]            andrstrnglc_e    0         [0, 0.002]
    arab_L_e         0         [0, 0.001]            imp_e            0         [0, 0.002]
    arach_e          0         [0, 0.001]            leuktrB4_e       0         [0, 0.002]
    c81coa_c         0         [0, 0.001]            leuktrC4_e       0         [0, 0.002]
    carn_e           0         [0, 0.001]            prostgh2_e       0         [0, 0.002]
    crvnc_e          0         [0, 0.001]            gmp_e            0         [-0.001, 0.002]
    dca_e            0         [0, 0.001]            octdececoa_c     0         [-0.001, 0.002]
    dd2coa_c         0         [0, 0.001]            pmtcoa_r         0         [-0.001, 0.002]
    decdicoa_c       0         [0, 0.001]            fucfucfucgal...  0         [0, 0.00133]
    dheas_e          0         [0, 0.001]            n5m2masn_g       0         [0, 0.00133]
    doco13ac_e       0         [0, 0.001]            core5_g          0         [0, 0.0012]
    dopa_e           0         [0, 0.001]            core7_g          0         [0, 0.0012]
    etoh_e           0         [0, 0.001]            core8_g          0         [0, 0.0012]
    fmn_e            0         [0, 0.001]            3bcrn_e          0         [0, 0.001]
    gam_e            0         [0, 0.001]            3ddcrn_e         0         [0, 0.001]
    gluala_e         0         [0, 0.001]            3deccrn_e        0         [0, 0.001]
    gtp_e            0         [0, 0.001]            3hdececrn_e      0         [0, 0.001]
    h2o2_e           0         [0, 0.001]            3hexdcrn_e       0         [0, 0.001]
    hista_e          0         [0, 0.001]            3hpvstet_e       0         [0, 0.001]
    hpdca_e          0         [0, 0.001]            3mlda_e          0         [0, 0.001]
    lgnc_e           0         [0, 0.001]            3octdec2crn_e    0         [0, 0.001]
    lnlc_e           0         [0, 0.001]            3octdeccrn_e     0         [0, 0.001]
    n2m2nmasn_e      0         [0, 0.001]            3octdece1crn_e   0         [0, 0.001]
    o2s_e            0         [0, 0.001]            3tdcrn_e         0         [0, 0.001]
    ocdca_e          0         [0, 0.001]            3tetd7ecoacrn_e  0         [0, 0.001]
    octa_e           0         [0, 0.001]            3thexddcoacrn_e  0         [0, 0.001]
    ppa_e            0         [0, 0.001]            3ttetddcoacrn_e  0         [0, 0.001]
    pre_prot_r       0         [0, 0.001]            5mta_e           0         [0, 0.001]
    prgstrn_e        0         [0, 0.001]            abt_e            0         [0, 0.001]
    prostge2_e       0         [0, 0.001]            aprgstrn_e       0         [0, 0.001]
    ptdca_e          0         [0, 0.001]            c4crn_e          0         [0, 0.001]
    ptvst_e          0         [0, 0.001]            c51crn_e         0         [0, 0.001]
    retinol_e        0         [0, 0.001]            c5dc_e           0         [0, 0.001]
    spmd_e           0         [0, 0.001]            c6crn_e          0         [0, 0.001]
    tetdece1coa_c    0         [0, 0.001]            c81crn_e         0         [0, 0.001]
    tmndnc_e         0         [0, 0.001]            c8crn_e          0         [0, 0.001]
    tststerone_e     0         [0, 0.001]            clpnd_e          0         [0, 0.001]
    ttdca_e          0         [0, 0.001]            ddece1crn_e      0         [0, 0.001]
    whhdca_e         0         [0, 0.001]            decdicrn_e       0         [0, 0.001]
    T_antigen_g      0         [-0.0002, 0.001]      dem2emgacpai...  0         [0, 0.001]
    andrstrn_e       0         [-0.001, 0.001]       dgpi_prot_hs_r   0         [0, 0.001]
    citr_L_c         0         [-0.001, 0.001]       fol_e            0         [0, 0.001]
    nrpphr_e         0         [-0.001, 0.001]       gpi_sig_r        0         [0, 0.001]
                                                     gthox_e          0         [0, 0.001]
                                                     ivcrn_e          0         [0, 0.001]
                                                     lac_D_e          0         [0, 0.001]
                                                     mem2emgacpai...  0         [0, 0.001]
                                                     oxyp_e           0         [0, 0.001]
                                                     prostgi2_e       0         [0, 0.001]
                                                     ptvstm3_e        0         [0, 0.001]
                                                     retn_e           0         [0, 0.001]
                                                     ribflv_e         0         [0, 0.001]
                                                     rsvlac_e         0         [0, 0.001]
                                                     sprm_c           0         [0, 0.001]
                                                     sulpacmp_e       0         [0, 0.001]
                                                     tetdece1crn_e    0         [0, 0.001]
                                                     tststeroneglc_e  0         [0, 0.001]
                                                     ttdcrn_e         0         [0, 0.001]
                                                     Asn_X_Ser_Thr_l  0         [-0.000333, 0.001]
                                                     5adtststerone_e  0         [-0.001, 0.001]
                                                     citr_L_e         0         [-0.001, 0.001]
                                                     gdp_e            0         [-0.001, 0.001]
                                                     idp_e            0         [-0.001, 0.001]
                                                     itp_e            0         [-0.001, 0.001]
                                                     leuktrA4_e       0         [-0.001, 0.001]



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

    Common reactions size 2511
    Normal unique reactions size 720
    Disease unique reactions size 386



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
      <th>cancer_rate</th>
      <th>normal_rate</th>
      <th>Total_rate</th>
    </tr>
  </thead>
  <tbody>
    <tr>
      <th>PIt8</th>
      <td>0.024931</td>
      <td>0.0242436</td>
      <td>0.026599</td>
      <td>HGNC:10946 or HGNC:10947</td>
      <td>CellLines</td>
      <td>0.0242436</td>
      <td>0.937289</td>
      <td>1</td>
      <td>0.937289</td>
    </tr>
  </tbody>
</table>
</div>





```python
FBA=cobra.flux_analysis.pfba(model_HeLaMean)
FBA.fluxes.to_csv("HeLa_flux.csv")

FBA=cobra.flux_analysis.pfba(model_KeratinocytesMean)
FBA.fluxes.to_csv("Keratinocytes_flux.csv")
```


```python
cobra.io.write_sbml_model(model_HeLaMean, "Full_HeLa_model_0319_DEMEM6429_n5_FINAL.sbml")
cobra.io.write_sbml_model(model_KeratinocytesMean, "Full_Kerat_model_0319_DEMEM6429_n5_FINAL.sbml")
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
%matplotlib inline  

from matplotlib_venn import venn2, venn3
 
# Make the diagram
venn2([set([ k.id for k in model_HeLaMean.reactions ]),set([ k.id for k in model_KeratinocytesMean.reactions ])],("HeLa","Keratinocytes"))

```


```python
import pandas
from time import time

import cobra.test

from cobra.flux_analysis import (
    single_gene_deletion, single_reaction_deletion)

deletion_results=single_reaction_deletion(cobra_model=model_HeLaMean, reaction_list=model_HeLaMean.reactions )

deletion_results[deletion_results.flux<0.01].sort_values(ascending=True, by='flux')


```


```python
import pandas
from time import time

import cobra.test

from cobra.flux_analysis import (
    single_gene_deletion, single_reaction_deletion, double_gene_deletion,
    double_reaction_deletion)

deletion_results=single_reaction_deletion(model_KeratinocytesMean, model_KeratinocytesMean.reactions)

deletion_results[deletion_results.flux<0.01].sort_values(ascending=True, by='flux')
```

## Sanity Checks


```python
%matplotlib inline
from cobra.flux_analysis import production_envelope

prod_env = production_envelope(
    model_HeLaMean, ["EX_co2_LPAREN_e_RPAREN_","EX_glc_LPAREN_e_RPAREN_","EX_gly_LPAREN_e_RPAREN_",])

```



```python
def model_create(rec, conf, nsize, m, pf):
    opt=CORDA(model=rec, confidence=conf, n=nsize, met_prod=m,  penalty_factor=pf ) 
    opt.build()
    return(opt)
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

    CPU times: user 2min 5s, sys: 28.3 s, total: 2min 34s
    Wall time: 45min 29s



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

    <Solution 0.026 at 0x7fd3789b3668>
    <Solution 0.027 at 0x7fd378b33908>
    <Solution 0.027 at 0x7fd3e82a1c88>
    <Solution 0.027 at 0x7fd35f613208>
    <Solution 0.027 at 0x7fd3e9396ba8>



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

    <Solution 0.024 at 0x7fd3479ebe48>
    <Solution 0.024 at 0x7fd370208358>
    <Solution 0.024 at 0x7fd34fa3e320>
    <Solution 0.024 at 0x7fd34f996080>
    <Solution 0.024 at 0x7fd346884d68>



```python
R_CLines=get_drugable_targets(
    normal_Model=model_KeratinocytesMean1, 
    disease_Model=model_HeLaMean1, 
    model_name="CellLines")
R_CLines["Dise_ratio"]=R_CLines["del_dise_flux"]/R_CLines["dise_flux"]
R_CLines["Norm_ratio"]=R_CLines["del_norm_flux"]/R_CLines["norm_flux"]
R_CLines[R_CLines["Dise_ratio"]< 0.999].sort_values(by="Dise_ratio", ascending=True)
```

    Common reactions size 2451
    Normal unique reactions size 687
    Disease unique reactions size 401







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
      <th>Dise_ratio</th>
      <th>Norm_ratio</th>
    </tr>
  </thead>
  <tbody>
    <tr>
      <th>ENO</th>
      <td>0.008</td>
      <td>0.008</td>
      <td>0.025934</td>
      <td>HGNC:3350 or HGNC:3354 or HGNC:3354 or HGNC:3353</td>
      <td>CellLines</td>
      <td>0.0242436</td>
      <td>0.308475</td>
      <td>0.329984</td>
    </tr>
    <tr>
      <th>EX_glc_LPAREN_e_RPAREN_</th>
      <td>0.008</td>
      <td>0.00897042</td>
      <td>0.025934</td>
      <td></td>
      <td>CellLines</td>
      <td>0.0242436</td>
      <td>0.308475</td>
      <td>0.370012</td>
    </tr>
    <tr>
      <th>GAPD</th>
      <td>0.008</td>
      <td>0.008</td>
      <td>0.025934</td>
      <td>HGNC:24864 or HGNC:4141</td>
      <td>CellLines</td>
      <td>0.0242436</td>
      <td>0.308475</td>
      <td>0.329984</td>
    </tr>
    <tr>
      <th>GLCt4</th>
      <td>0.008</td>
      <td>0.00897042</td>
      <td>0.025934</td>
      <td>HGNC:11037 or HGNC:11038 or HGNC:22146 or HGNC...</td>
      <td>CellLines</td>
      <td>0.0242436</td>
      <td>0.308475</td>
      <td>0.370012</td>
    </tr>
    <tr>
      <th>PGK</th>
      <td>0.008</td>
      <td>0.008</td>
      <td>0.025934</td>
      <td>HGNC:8896 or HGNC:8898</td>
      <td>CellLines</td>
      <td>0.0242436</td>
      <td>0.308475</td>
      <td>0.329984</td>
    </tr>
    <tr>
      <th>PGM</th>
      <td>0.008</td>
      <td>0.008</td>
      <td>0.025934</td>
      <td>HGNC:1093 or HGNC:8888 or HGNC:8889</td>
      <td>CellLines</td>
      <td>0.0242436</td>
      <td>0.308475</td>
      <td>0.329984</td>
    </tr>
    <tr>
      <th>TPI</th>
      <td>0.008</td>
      <td>0.00884086</td>
      <td>0.025934</td>
      <td>HGNC:12009</td>
      <td>CellLines</td>
      <td>0.0242436</td>
      <td>0.308475</td>
      <td>0.364668</td>
    </tr>
    <tr>
      <th>biomass_reaction</th>
      <td>0.01</td>
      <td>0.01</td>
      <td>0.025934</td>
      <td></td>
      <td>CellLines</td>
      <td>0.0242436</td>
      <td>0.385594</td>
      <td>0.41248</td>
    </tr>
    <tr>
      <th>biomass_protein</th>
      <td>0.0141643</td>
      <td>0.0141643</td>
      <td>0.025934</td>
      <td></td>
      <td>CellLines</td>
      <td>0.0242436</td>
      <td>0.546167</td>
      <td>0.584249</td>
    </tr>
    <tr>
      <th>LYStiDF</th>
      <td>0.0168888</td>
      <td>0.0242436</td>
      <td>0.025934</td>
      <td>HGNC:11057 or HGNC:11060 or HGNC:11061 or HGNC...</td>
      <td>CellLines</td>
      <td>0.0242436</td>
      <td>0.65122</td>
      <td>1</td>
    </tr>
    <tr>
      <th>EX_lys_L_LPAREN_e_RPAREN_</th>
      <td>0.0168888</td>
      <td>0.0168888</td>
      <td>0.025934</td>
      <td></td>
      <td>CellLines</td>
      <td>0.0242436</td>
      <td>0.65122</td>
      <td>0.696627</td>
    </tr>
    <tr>
      <th>EX_leu_L_LPAREN_e_RPAREN_</th>
      <td>0.0183305</td>
      <td>0.0183305</td>
      <td>0.025934</td>
      <td></td>
      <td>CellLines</td>
      <td>0.0242436</td>
      <td>0.706811</td>
      <td>0.756095</td>
    </tr>
  </tbody>
</table>
</div>




```python
R_CLines=get_drugable_targets(
    normal_Model=model_KeratinocytesMean2, 
    disease_Model=model_HeLaMean2, 
    model_name="CellLines")
R_CLines["Dise_ratio"]=R_CLines["del_dise_flux"]/R_CLines["dise_flux"]
R_CLines["Norm_ratio"]=R_CLines["del_norm_flux"]/R_CLines["norm_flux"]
R_CLines[R_CLines["Dise_ratio"]< 0.999].sort_values(by="Dise_ratio", ascending=True)
```

    Common reactions size 2491
    Normal unique reactions size 713
    Disease unique reactions size 409


    cobra/util/solver.py:404 [1;31mUserWarning[0m: solver status is 'infeasible'





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
      <th>Dise_ratio</th>
      <th>Norm_ratio</th>
    </tr>
  </thead>
  <tbody>
    <tr>
      <th>ENO</th>
      <td>0.008</td>
      <td>0.008</td>
      <td>0.026599</td>
      <td>HGNC:3350 or HGNC:3354 or HGNC:3354 or HGNC:3353</td>
      <td>CellLines</td>
      <td>0.0242436</td>
      <td>0.300763</td>
      <td>0.329984</td>
    </tr>
    <tr>
      <th>EX_glc_LPAREN_e_RPAREN_</th>
      <td>0.008</td>
      <td>0.00901176</td>
      <td>0.026599</td>
      <td></td>
      <td>CellLines</td>
      <td>0.0242436</td>
      <td>0.300763</td>
      <td>0.371717</td>
    </tr>
    <tr>
      <th>GAPD</th>
      <td>0.008</td>
      <td>0.008</td>
      <td>0.026599</td>
      <td>HGNC:24864 or HGNC:4141</td>
      <td>CellLines</td>
      <td>0.0242436</td>
      <td>0.300763</td>
      <td>0.329984</td>
    </tr>
    <tr>
      <th>GLCt4</th>
      <td>0.008</td>
      <td>0.00901176</td>
      <td>0.026599</td>
      <td>HGNC:11037 or HGNC:11038 or HGNC:22146 or HGNC...</td>
      <td>CellLines</td>
      <td>0.0242436</td>
      <td>0.300763</td>
      <td>0.371717</td>
    </tr>
    <tr>
      <th>PGK</th>
      <td>0.008</td>
      <td>0.008</td>
      <td>0.026599</td>
      <td>HGNC:8896 or HGNC:8898</td>
      <td>CellLines</td>
      <td>0.0242436</td>
      <td>0.300763</td>
      <td>0.329984</td>
    </tr>
    <tr>
      <th>PGM</th>
      <td>0.008</td>
      <td>0.008</td>
      <td>0.026599</td>
      <td>HGNC:1093 or HGNC:8888 or HGNC:8889</td>
      <td>CellLines</td>
      <td>0.0242436</td>
      <td>0.300763</td>
      <td>0.329984</td>
    </tr>
    <tr>
      <th>TPI</th>
      <td>0.008</td>
      <td>0.00884086</td>
      <td>0.026599</td>
      <td>HGNC:12009</td>
      <td>CellLines</td>
      <td>0.0242436</td>
      <td>0.300763</td>
      <td>0.364668</td>
    </tr>
    <tr>
      <th>biomass_reaction</th>
      <td>0.01</td>
      <td>0.01</td>
      <td>0.026599</td>
      <td></td>
      <td>CellLines</td>
      <td>0.0242436</td>
      <td>0.375954</td>
      <td>0.41248</td>
    </tr>
    <tr>
      <th>biomass_protein</th>
      <td>0.0141643</td>
      <td>0.0141643</td>
      <td>0.026599</td>
      <td></td>
      <td>CellLines</td>
      <td>0.0242436</td>
      <td>0.532513</td>
      <td>0.584249</td>
    </tr>
    <tr>
      <th>LYStiDF</th>
      <td>0.0168888</td>
      <td>0.0242436</td>
      <td>0.026599</td>
      <td>HGNC:11057 or HGNC:11060 or HGNC:11061 or HGNC...</td>
      <td>CellLines</td>
      <td>0.0242436</td>
      <td>0.634939</td>
      <td>1</td>
    </tr>
    <tr>
      <th>EX_lys_L_LPAREN_e_RPAREN_</th>
      <td>0.0168888</td>
      <td>0.0168888</td>
      <td>0.026599</td>
      <td></td>
      <td>CellLines</td>
      <td>0.0242436</td>
      <td>0.634939</td>
      <td>0.696627</td>
    </tr>
    <tr>
      <th>EX_leu_L_LPAREN_e_RPAREN_</th>
      <td>0.0183305</td>
      <td>0.0183305</td>
      <td>0.026599</td>
      <td></td>
      <td>CellLines</td>
      <td>0.0242436</td>
      <td>0.689141</td>
      <td>0.756095</td>
    </tr>
  </tbody>
</table>
</div>




```python
R_CLines=get_drugable_targets(
    normal_Model=model_KeratinocytesMean3, 
    disease_Model=model_HeLaMean3, 
    model_name="CellLines")
R_CLines["Dise_ratio"]=R_CLines["del_dise_flux"]/R_CLines["dise_flux"]
R_CLines["Norm_ratio"]=R_CLines["del_norm_flux"]/R_CLines["norm_flux"]
R_CLines[R_CLines["Dise_ratio"]< 0.999].sort_values(by="Dise_ratio", ascending=True)
```

    Common reactions size 2509
    Normal unique reactions size 722
    Disease unique reactions size 386





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
      <th>Dise_ratio</th>
      <th>Norm_ratio</th>
    </tr>
  </thead>
  <tbody>
    <tr>
      <th>ENO</th>
      <td>0.008</td>
      <td>0.008</td>
      <td>0.026599</td>
      <td>HGNC:3350 or HGNC:3354 or HGNC:3354 or HGNC:3353</td>
      <td>CellLines</td>
      <td>0.0242436</td>
      <td>0.300763</td>
      <td>0.329984</td>
    </tr>
    <tr>
      <th>EX_glc_LPAREN_e_RPAREN_</th>
      <td>0.008</td>
      <td>0.00908105</td>
      <td>0.026599</td>
      <td></td>
      <td>CellLines</td>
      <td>0.0242436</td>
      <td>0.300763</td>
      <td>0.374575</td>
    </tr>
    <tr>
      <th>GAPD</th>
      <td>0.008</td>
      <td>0.008</td>
      <td>0.026599</td>
      <td>HGNC:24864 or HGNC:4141</td>
      <td>CellLines</td>
      <td>0.0242436</td>
      <td>0.300763</td>
      <td>0.329984</td>
    </tr>
    <tr>
      <th>GLCt4</th>
      <td>0.008</td>
      <td>0.00908105</td>
      <td>0.026599</td>
      <td>HGNC:11037 or HGNC:11038 or HGNC:22146 or HGNC...</td>
      <td>CellLines</td>
      <td>0.0242436</td>
      <td>0.300763</td>
      <td>0.374575</td>
    </tr>
    <tr>
      <th>PGK</th>
      <td>0.008</td>
      <td>0.008</td>
      <td>0.026599</td>
      <td>HGNC:8896 or HGNC:8898</td>
      <td>CellLines</td>
      <td>0.0242436</td>
      <td>0.300763</td>
      <td>0.329984</td>
    </tr>
    <tr>
      <th>PGM</th>
      <td>0.008</td>
      <td>0.008</td>
      <td>0.026599</td>
      <td>HGNC:1093 or HGNC:8888 or HGNC:8889</td>
      <td>CellLines</td>
      <td>0.0242436</td>
      <td>0.300763</td>
      <td>0.329984</td>
    </tr>
    <tr>
      <th>TPI</th>
      <td>0.008</td>
      <td>0.00887755</td>
      <td>0.026599</td>
      <td>HGNC:12009</td>
      <td>CellLines</td>
      <td>0.0242436</td>
      <td>0.300763</td>
      <td>0.366181</td>
    </tr>
    <tr>
      <th>biomass_reaction</th>
      <td>0.01</td>
      <td>0.01</td>
      <td>0.026599</td>
      <td></td>
      <td>CellLines</td>
      <td>0.0242436</td>
      <td>0.375954</td>
      <td>0.41248</td>
    </tr>
    <tr>
      <th>biomass_protein</th>
      <td>0.0141643</td>
      <td>0.0141643</td>
      <td>0.026599</td>
      <td></td>
      <td>CellLines</td>
      <td>0.0242436</td>
      <td>0.532513</td>
      <td>0.584249</td>
    </tr>
    <tr>
      <th>EX_lys_L_LPAREN_e_RPAREN_</th>
      <td>0.0168888</td>
      <td>0.0168888</td>
      <td>0.026599</td>
      <td></td>
      <td>CellLines</td>
      <td>0.0242436</td>
      <td>0.634939</td>
      <td>0.696627</td>
    </tr>
    <tr>
      <th>EX_leu_L_LPAREN_e_RPAREN_</th>
      <td>0.0183305</td>
      <td>0.0183305</td>
      <td>0.026599</td>
      <td></td>
      <td>CellLines</td>
      <td>0.0242436</td>
      <td>0.689141</td>
      <td>0.756095</td>
    </tr>
    <tr>
      <th>PIt8</th>
      <td>0.024931</td>
      <td>0.0242436</td>
      <td>0.026599</td>
      <td>HGNC:10946 or HGNC:10947</td>
      <td>CellLines</td>
      <td>0.0242436</td>
      <td>0.937289</td>
      <td>1</td>
    </tr>
  </tbody>
</table>
</div>




```python
R_CLines=get_drugable_targets(
    normal_Model=model_KeratinocytesMean4, 
    disease_Model=model_HeLaMean4, 
    model_name="CellLines")
R_CLines["Dise_ratio"]=R_CLines["del_dise_flux"]/R_CLines["dise_flux"]
R_CLines["Norm_ratio"]=R_CLines["del_norm_flux"]/R_CLines["norm_flux"]
R_CLines[R_CLines["Dise_ratio"]< 0.999].sort_values(by="Dise_ratio", ascending=True)
```

    Common reactions size 2528
    Normal unique reactions size 713
    Disease unique reactions size 386





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
      <th>Dise_ratio</th>
      <th>Norm_ratio</th>
    </tr>
  </thead>
  <tbody>
    <tr>
      <th>ENO</th>
      <td>0.008</td>
      <td>0.008</td>
      <td>0.026599</td>
      <td>HGNC:3350 or HGNC:3354 or HGNC:3354 or HGNC:3353</td>
      <td>CellLines</td>
      <td>0.0242436</td>
      <td>0.300763</td>
      <td>0.329984</td>
    </tr>
    <tr>
      <th>EX_glc_LPAREN_e_RPAREN_</th>
      <td>0.008</td>
      <td>0.00908105</td>
      <td>0.026599</td>
      <td></td>
      <td>CellLines</td>
      <td>0.0242436</td>
      <td>0.300763</td>
      <td>0.374575</td>
    </tr>
    <tr>
      <th>GAPD</th>
      <td>0.008</td>
      <td>0.008</td>
      <td>0.026599</td>
      <td>HGNC:24864 or HGNC:4141</td>
      <td>CellLines</td>
      <td>0.0242436</td>
      <td>0.300763</td>
      <td>0.329984</td>
    </tr>
    <tr>
      <th>GLCt4</th>
      <td>0.008</td>
      <td>0.00908105</td>
      <td>0.026599</td>
      <td>HGNC:11037 or HGNC:11038 or HGNC:22146 or HGNC...</td>
      <td>CellLines</td>
      <td>0.0242436</td>
      <td>0.300763</td>
      <td>0.374575</td>
    </tr>
    <tr>
      <th>PGK</th>
      <td>0.008</td>
      <td>0.008</td>
      <td>0.026599</td>
      <td>HGNC:8896 or HGNC:8898</td>
      <td>CellLines</td>
      <td>0.0242436</td>
      <td>0.300763</td>
      <td>0.329984</td>
    </tr>
    <tr>
      <th>PGM</th>
      <td>0.008</td>
      <td>0.008</td>
      <td>0.026599</td>
      <td>HGNC:1093 or HGNC:8888 or HGNC:8889</td>
      <td>CellLines</td>
      <td>0.0242436</td>
      <td>0.300763</td>
      <td>0.329984</td>
    </tr>
    <tr>
      <th>TPI</th>
      <td>0.008</td>
      <td>0.00887755</td>
      <td>0.026599</td>
      <td>HGNC:12009</td>
      <td>CellLines</td>
      <td>0.0242436</td>
      <td>0.300763</td>
      <td>0.366181</td>
    </tr>
    <tr>
      <th>biomass_reaction</th>
      <td>0.01</td>
      <td>0.01</td>
      <td>0.026599</td>
      <td></td>
      <td>CellLines</td>
      <td>0.0242436</td>
      <td>0.375954</td>
      <td>0.41248</td>
    </tr>
    <tr>
      <th>biomass_protein</th>
      <td>0.0141643</td>
      <td>0.0141643</td>
      <td>0.026599</td>
      <td></td>
      <td>CellLines</td>
      <td>0.0242436</td>
      <td>0.532513</td>
      <td>0.584249</td>
    </tr>
    <tr>
      <th>EX_lys_L_LPAREN_e_RPAREN_</th>
      <td>0.0168888</td>
      <td>0.0168888</td>
      <td>0.026599</td>
      <td></td>
      <td>CellLines</td>
      <td>0.0242436</td>
      <td>0.634939</td>
      <td>0.696627</td>
    </tr>
    <tr>
      <th>EX_leu_L_LPAREN_e_RPAREN_</th>
      <td>0.0183305</td>
      <td>0.0183305</td>
      <td>0.026599</td>
      <td></td>
      <td>CellLines</td>
      <td>0.0242436</td>
      <td>0.689141</td>
      <td>0.756095</td>
    </tr>
  </tbody>
</table>
</div>




```python
R_CLines=get_drugable_targets(
    normal_Model=model_KeratinocytesMean5, 
    disease_Model=model_HeLaMean5, 
    model_name="CellLines")
R_CLines["Dise_ratio"]=R_CLines["del_dise_flux"]/R_CLines["dise_flux"]
R_CLines["Norm_ratio"]=R_CLines["del_norm_flux"]/R_CLines["norm_flux"]
R_CLines[R_CLines["Dise_ratio"]< 0.999].sort_values(by="Dise_ratio", ascending=True)
```

    Common reactions size 2520
    Normal unique reactions size 720
    Disease unique reactions size 384


    cobra/util/solver.py:404 [1;31mUserWarning[0m: solver status is 'infeasible'





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
      <th>Dise_ratio</th>
      <th>Norm_ratio</th>
    </tr>
  </thead>
  <tbody>
    <tr>
      <th>ENO</th>
      <td>0.008</td>
      <td>0.008</td>
      <td>0.026599</td>
      <td>HGNC:3350 or HGNC:3354 or HGNC:3354 or HGNC:3353</td>
      <td>CellLines</td>
      <td>0.0242436</td>
      <td>0.300763</td>
      <td>0.329984</td>
    </tr>
    <tr>
      <th>EX_glc_LPAREN_e_RPAREN_</th>
      <td>0.008</td>
      <td>0.00908105</td>
      <td>0.026599</td>
      <td></td>
      <td>CellLines</td>
      <td>0.0242436</td>
      <td>0.300763</td>
      <td>0.374575</td>
    </tr>
    <tr>
      <th>GAPD</th>
      <td>0.008</td>
      <td>0.008</td>
      <td>0.026599</td>
      <td>HGNC:24864 or HGNC:4141</td>
      <td>CellLines</td>
      <td>0.0242436</td>
      <td>0.300763</td>
      <td>0.329984</td>
    </tr>
    <tr>
      <th>GLCt4</th>
      <td>0.008</td>
      <td>0.00908105</td>
      <td>0.026599</td>
      <td>HGNC:11037 or HGNC:11038 or HGNC:22146 or HGNC...</td>
      <td>CellLines</td>
      <td>0.0242436</td>
      <td>0.300763</td>
      <td>0.374575</td>
    </tr>
    <tr>
      <th>PGK</th>
      <td>0.008</td>
      <td>0.008</td>
      <td>0.026599</td>
      <td>HGNC:8896 or HGNC:8898</td>
      <td>CellLines</td>
      <td>0.0242436</td>
      <td>0.300763</td>
      <td>0.329984</td>
    </tr>
    <tr>
      <th>PGM</th>
      <td>0.008</td>
      <td>0.008</td>
      <td>0.026599</td>
      <td>HGNC:1093 or HGNC:8888 or HGNC:8889</td>
      <td>CellLines</td>
      <td>0.0242436</td>
      <td>0.300763</td>
      <td>0.329984</td>
    </tr>
    <tr>
      <th>TPI</th>
      <td>0.008</td>
      <td>0.00887755</td>
      <td>0.026599</td>
      <td>HGNC:12009</td>
      <td>CellLines</td>
      <td>0.0242436</td>
      <td>0.300763</td>
      <td>0.366181</td>
    </tr>
    <tr>
      <th>biomass_reaction</th>
      <td>0.01</td>
      <td>0.01</td>
      <td>0.026599</td>
      <td></td>
      <td>CellLines</td>
      <td>0.0242436</td>
      <td>0.375954</td>
      <td>0.41248</td>
    </tr>
    <tr>
      <th>biomass_protein</th>
      <td>0.0141643</td>
      <td>0.0141643</td>
      <td>0.026599</td>
      <td></td>
      <td>CellLines</td>
      <td>0.0242436</td>
      <td>0.532513</td>
      <td>0.584249</td>
    </tr>
    <tr>
      <th>EX_lys_L_LPAREN_e_RPAREN_</th>
      <td>0.0168888</td>
      <td>0.0168888</td>
      <td>0.026599</td>
      <td></td>
      <td>CellLines</td>
      <td>0.0242436</td>
      <td>0.634939</td>
      <td>0.696627</td>
    </tr>
    <tr>
      <th>EX_leu_L_LPAREN_e_RPAREN_</th>
      <td>0.0183305</td>
      <td>0.0183305</td>
      <td>0.026599</td>
      <td></td>
      <td>CellLines</td>
      <td>0.0242436</td>
      <td>0.689141</td>
      <td>0.756095</td>
    </tr>
    <tr>
      <th>PIt8</th>
      <td>0.0257091</td>
      <td>0.0242436</td>
      <td>0.026599</td>
      <td>HGNC:10946 or HGNC:10947</td>
      <td>CellLines</td>
      <td>0.0242436</td>
      <td>0.966545</td>
      <td>1</td>
    </tr>
  </tbody>
</table>
</div>




```python
model_KeratinocytesMean1.summary(fva=True)
```

    IN FLUXES                                        OUT FLUXES                                     OBJECTIVES
    -----------------------------------------------  ---------------------------------------------  -----------------------
    id                   Flux  Range                 id                   Flux  Range               biomass_reac...  0.0242
    ---------------  --------  --------------------  ---------------  --------  ------------------
    glc_D_e          0.262     [0.199, 4.5]          h_e              0.512     [-0.632, 10]
    o2_e             0.041     [0.00246, 8.47]       lac_L_e          0.486     [-1, 9.73]
    ser_L_e          0.0167    [0.042, -1.12]        co2_e            0.0395    [0, 10]
    lys_L_e          0.0144    [0.0144, 0.0154]      pyr_e            0.03      [-0.11, 9.86]
    ile_L_e          0.0141    [0.00694, 0.105]      h2o_e            0.0292    [-1.54, 10]
    leu_L_e          0.0132    [0.0132, 0.0142]      ala_B_e          0.0156    [-0.001, 1.3]
    gly_e            0.0121    [0.03, -1.14]         dgsn_e           0.00353   [-0.001, 0.0198]
    pro_L_e          0.01      [0.01, 0.01]          gua_e            0.00281   [-0.001, 0.343]
    arg_L_e          0.00871   [0.00671, 0.00871]    Rtotal_e         0.002     [-0.001, 0.247]
    val_L_e          0.00865   [0.00855, 0.094]      pi_e             0.00127   [0.0155, -0.895]
    asp_L_e          0.00855   [0, 0.01]             ethamp_r         0.001     [0, 0.846]
    gln_L_e          0.0079    [-0.159, 0.584]       hdca_e           0.001     [-0.001, 0.51]
    thr_L_e          0.00758   [0.00758, 0.00758]    prostgh2_e       0.001     [0, 0.002]
    glu_L_e          0.00736   [0.01, -0.733]        pmtcoa_r         0.001     [-0.001, 0.002]
    asn_L_e          0.00677   [0.00677, -0.736]     lac_D_e          0.001     [0, 0.001]
    phe_L_e          0.00629   [0.00629, 0.066]      5adtststerone_e  0.001     [-0.001, 0.001]
    tyr_L_e          0.00387   [-0.0558, 0.104]      andrstrn_e       0.001     [0.001, -0.001]
    met_L_e          0.00371   [0.00371, 0.03]       so4_e            0.001     [0.0898, -1]
    his_L_e          0.00307   [0.00307, 0.00307]    dad_2_e          0.00068   [-0.001, 0.0198]
    chol_e           0.00217   [0.05, -0.128]        for_e            0.000495  [-0.001, 0.394]
    biomass_other_c  0.00131   [0.00131, 0.00131]    5mta_e           0.000495  [0, 0.001]
    cys_L_e          0.00113   [-0.0262, 0.0626]     sprm_c           0.000495  [0, 0.001]
    pe_hs_e          0.001     [0.000342, 0.001]     elaid_e          0.000108  [-0.001, 0.46]
    dheas_e          0.001     [0, 0.001]            orn_e            0.000105  [0, 0.0865]
    o2s_e            0.001     [0, 0.001]            glyc_S_e         0         [0, 2.73]
    orot_e           0.001     [0, 0.001]            oxa_e            0         [0, 2.17]
    pchol_hs_e       0.001     [0, 0.001]            gal_e            0         [0, 2.15]
    prostge2_e       0.001     [0, 0.001]            xylt_e           0         [0, 1.72]
    tag_hs_e         0.001     [0, 0.001]            acmana_e         0         [0, 1.32]
    tststerone_e     0.001     [0, 0.001]            ala_L_e          0         [-0.01, 1.16]
    dag_hs_e         0.001     [-0.000595, 0.001]    taur_c           0         [-0.001, 1.09]
    glyc3p_e         0.001     [-0.000595, 0.001]    taur_e           0         [-0.001, 1.09]
    lpchol_hs_e      0.001     [-0.000595, 0.001]    fe2_e            0         [-1, 1]
    ps_hs_e          0.001     [-0.000595, 0.001]    fe3_e            0         [-1, 1]
    cdp_e            0.001     [0.001, -0.002]       ha_pre1_e        0         [-0.001, 0.849]
    cmp_e            0.001     [0.001, -0.002]       5oxpro_e         0         [0, 0.744]
    gmp_e            0.001     [0.001, -0.002]       vacc_e           0         [-0.001, 0.46]
    octdececoa_c     0.001     [0.001, -0.002]       urate_e          0         [0, 0.429]
    cytd_e           0.001     [0.001, -0.00853]     ha_e             0         [-0.001, 0.425]
    dctp_n           0.001     [0.001, -0.00853]     chsterol_e       0         [0, 0.391]
    dcyt_e           0.001     [0.001, -0.00853]     cholate_e        0         [0, 0.378]
    udp_e            0.001     [0.001, -0.00853]     tdchola_e        0         [0, 0.354]
    ump_e            0.001     [0.001, -0.00853]     eicostet_e       0         [-0.001, 0.352]
    uri_e            0.001     [0.001, -0.00853]     ade_e            0         [-0.001, 0.348]
    utp_e            0.001     [0.001, -0.00853]     nrvnc_e          0         [0, 0.328]
    aicar_e          0.001     [0.001, -0.0136]      xolest2_hs_e     0         [-0.001, 0.245]
    atp_e            0.001     [0.001, -0.0136]      4hphac_e         0         [0, 0.16]
    ins_e            0.001     [0.001, -0.0136]      acac_e           0         [-0.001, 0.16]
    xmp_e            0.001     [0.001, -0.0136]      spc_hs_e         0         [0, 0.14]
    datp_n           0.001     [0.001, -0.0143]      q10h2_e          0         [0, 0.0917]
    dgtp_n           0.001     [0.001, -0.0198]      q10_e            0         [-0.001, 0.0903]
    din_e            0.001     [0.001, -0.0198]      cgly_e           0         [-0.001, 0.0898]
    sphs1p_e         0.001     [0.001, -0.366]       gthrd_e          0         [-0.001, 0.0898]
    strdnc_e         0.001     [0.001, -0.385]       pheacgln_e       0         [0, 0.0597]
    hxan_e           0.001     [0.001, -0.428]       2hb_e            0         [0, 0.0273]
    hdcea_e          0.001     [0.001, -0.522]       35cgmp_e         0         [0, 0.0146]
    HC00342_e        0.001     [0.001, -0.742]       camp_e           0         [0, 0.0146]
    cit_e            0.001     [0.001, -0.742]       adn_e            0         [-0.001, 0.0136]
    glyc_e           0.001     [0.001, -2.73]        bilglcur_e       0         [0, 0.0108]
    ura_e            0.000971  [0.001, -0.00853]     co_e             0         [0, 0.0108]
    inost_e          0.000565  [0.000565, 0.000894]  pheme_e          0         [0, 0.0108]
    spmd_e           0.000495  [0, 0.001]            dttp_m           0         [0, 0.0102]
    ahcys_e          0.000495  [0.001, -0.0136]      Ser_Thr_l        0         [-0.001, 0.005]
    gsn_e            0.00046   [0.001, -0.0136]      fuc14galacgl...  0         [0, 0.004]
    sph1p_e          0.000424  [0.001, -0.367]       fuc_L_e          0         [0, 0.004]
    pglyc_hs_e       0.000353  [0.000353, 0.000353]  galfucgalacg...  0         [0, 0.004]
    pe_hs_r          0.000342  [0.000342, 0.001]     man_e            0         [-0.001, 0.003]
    trp_L_e          0.000323  [0.000323, 0.000323]  meoh_e           0         [-0.001, 0.003]
    dttp_n           0.000317  [0.001, -0.00921]     34dhoxpeg_e      0         [0, 0.002]
    amp_e            0.000302  [0, 0.001]            5adtststeron...  0         [0, 0.002]
    3hpvs_e          0         [0, 0.001]            andrstrnglc_e    0         [0, 0.002]
    4hpro_LT_m       0         [0, 0.001]            imp_e            0         [0, 0.002]
    5HPET_c          0         [0, 0.001]            leuktrB4_e       0         [0, 0.002]
    5fthf_e          0         [0, 0.001]            leuktrC4_e       0         [0, 0.002]
    HC00822_e        0         [0, 0.001]            ctp_e            0         [-0.001, 0.002]
    acmp_e           0         [0, 0.001]            fucfucfucgal...  0         [0, 0.00133]
    adrn_e           0         [0, 0.001]            n5m2masn_g       0         [0, 0.00133]
    allop_e          0         [0, 0.001]            core5_g          0         [0, 0.0012]
    arab_L_e         0         [0, 0.001]            core7_g          0         [0, 0.0012]
    arach_e          0         [0, 0.001]            core8_g          0         [0, 0.0012]
    crn_e            0         [0, 0.001]            3bcrn_e          0         [0, 0.001]
    crvnc_e          0         [0, 0.001]            3ddcrn_e         0         [0, 0.001]
    dca_e            0         [0, 0.001]            3deccrn_e        0         [0, 0.001]
    dd2coa_c         0         [0, 0.001]            3hdececrn_e      0         [0, 0.001]
    doco13ac_e       0         [0, 0.001]            3hexdcrn_e       0         [0, 0.001]
    dopa_e           0         [0, 0.001]            3hpvstet_e       0         [0, 0.001]
    etoh_e           0         [0, 0.001]            3ivcrn_e         0         [0, 0.001]
    fmn_e            0         [0, 0.001]            3mlda_e          0         [0, 0.001]
    gam_e            0         [0, 0.001]            3octdec2crn_e    0         [0, 0.001]
    gluala_e         0         [0, 0.001]            3octdeccrn_e     0         [0, 0.001]
    gtp_e            0         [0, 0.001]            3octdece1crn_e   0         [0, 0.001]
    h2o2_e           0         [0, 0.001]            3tdcrn_e         0         [0, 0.001]
    hista_e          0         [0, 0.001]            3tetd7ecoacrn_e  0         [0, 0.001]
    hpdca_e          0         [0, 0.001]            3thexddcoacrn_e  0         [0, 0.001]
    lgnc_e           0         [0, 0.001]            3ttetddcoacrn_e  0         [0, 0.001]
    lnlc_e           0         [0, 0.001]            abt_e            0         [0, 0.001]
    n2m2nmasn_e      0         [0, 0.001]            aprgstrn_e       0         [0, 0.001]
    ocdca_e          0         [0, 0.001]            c4crn_e          0         [0, 0.001]
    octa_e           0         [0, 0.001]            c51crn_e         0         [0, 0.001]
    ppa_e            0         [0, 0.001]            c5dc_e           0         [0, 0.001]
    prgstrn_e        0         [0, 0.001]            c6crn_e          0         [0, 0.001]
    ptdca_e          0         [0, 0.001]            c8crn_e          0         [0, 0.001]
    ptvst_e          0         [0, 0.001]            clpnd_e          0         [0, 0.001]
    rsv_e            0         [0, 0.001]            ddece1crn_e      0         [0, 0.001]
    tacr_e           0         [0, 0.001]            epoxtac_e        0         [0, 0.001]
    tmndnc_e         0         [0, 0.001]            fol_e            0         [0, 0.001]
    ttdca_e          0         [0, 0.001]            gthox_e          0         [0, 0.001]
    whhdca_e         0         [0, 0.001]            ivcrn_e          0         [0, 0.001]
    T_antigen_g      0         [-0.0002, 0.001]      ndersv_e         0         [0, 0.001]
    dcsptn1_e        0         [-0.001, 0.001]       oxyp_e           0         [0, 0.001]
    nrpphr_e         0         [-0.001, 0.001]       prostgi2_e       0         [0, 0.001]
    pre_prot_r       0         [0, 0.000329]         ptvstm3_e        0         [0, 0.001]
    fald_e           0         [0.001, -0.003]       ribflv_e         0         [0, 0.001]
                                                     rsvlac_e         0         [0, 0.001]
                                                     sulpacmp_e       0         [0, 0.001]
                                                     tststeroneglc_e  0         [0, 0.001]
                                                     ttdcrn_e         0         [0, 0.001]
                                                     Asn_X_Ser_Thr_l  0         [-0.000333, 0.001]
                                                     citr_L_c         0         [-0.001, 0.001]
                                                     citr_L_e         0         [-0.001, 0.001]
                                                     gdp_e            0         [-0.001, 0.001]
                                                     idp_e            0         [-0.001, 0.001]
                                                     itp_e            0         [-0.001, 0.001]
                                                     leuktrA4_e       0         [-0.001, 0.001]
                                                     dem2emgacpai...  0         [0, 0.000329]
                                                     gpi_sig_r        0         [0, 0.000329]



```python
model_HeLaMean1.summary(fva=True)
```

    IN FLUXES                                        OUT FLUXES                                    OBJECTIVES
    -----------------------------------------------  --------------------------------------------  -----------------------
    id                   Flux  Range                 id                   Flux  Range              biomass_reac...  0.0259
    ---------------  --------  --------------------  ---------------  --------  -----------------
    glc_D_e          0.301     [0.251, 4.5]          h_e              0.472     [-0.731, 10]
    fe2_e            0.0676    [1, -1]               lac_L_e          0.439     [-1, 9.63]
    o2_e             0.0672    [0.00261, 10]         pyr_e            0.105     [-0.11, 10]
    arg_L_e          0.0424    [0, 0.084]            fe3_e            0.0676    [1, -1]
    ser_L_e          0.042     [0.042, -0.893]       co2_e            0.0595    [-0.001, 6.62]
    lys_L_e          0.0154    [0.0154, 0.0164]      glyc_e           0.0548    [-0.001, 0.937]
    leu_L_e          0.0141    [0.0141, 0.0151]      h2o_e            0.0448    [-5.19, 10]
    gly_e            0.0134    [0.03, -1.7]          ala_L_e          0.0407    [-0.01, 0.996]
    gln_L_e          0.013     [0.0043, 0.584]       ptrc_e           0.034     [-0.001, 0.184]
    glu_L_e          0.01      [-0.001, 0.01]        urea_e           0.033     [0, 0.111]
    pro_L_e          0.01      [0.01, -0.0972]       so4_e            0.0106    [0.0934, -0.164]
    cys_L_e          0.00979   [-0.0288, 0.0626]     ac_e             0.00859   [-0.001, 0.0674]
    val_L_e          0.00951   [0.00915, 0.094]      ura_e            0.00675   [-0.001, 0.00736]
    asp_L_e          0.00915   [0, 0.01]             gua_e            0.00481   [-0.001, 0.328]
    thr_L_e          0.00811   [0.00811, 0.00811]    ade_e            0.0046    [-0.001, 0.247]
    ile_L_e          0.00742   [0.00742, 0.105]      pi_e             0.00214   [0.00669, -1]
    asn_L_e          0.00725   [0.01, -0.572]        vacc_e           0.002     [-0.001, 0.506]
    phe_L_e          0.00673   [0.00673, 0.066]      hdcea_e          0.001     [-0.001, 0.571]
    tyr_L_e          0.00414   [-0.0551, 0.104]      lnlc_e           0.001     [-0.001, 0.509]
    his_L_e          0.00328   [0.00328, 0.00328]    fol_c            0.001     [-0.001, 0.052]
    met_L_e          0.00297   [0.00197, 0.03]       thym_e           0.001     [-0.001, 0.00969]
    chol_e           0.0029    [0.0029, 0.05]        ribflv_e         0.001     [0, 0.002]
    biomass_other_c  0.0014    [0.0014, 0.0014]      3bcrn_e          0.001     [0, 0.001]
    inost_e          0.00111   [0.00111, -2.83]      lac_D_e          0.001     [0, 0.001]
    dag_hs_e         0.001     [0.001, 0.001]        leuktrD4_e       0.000704  [-0.001, 0.002]
    lpchol_hs_e      0.001     [0.001, 0.001]        dttp_m           0.00066   [0, 0.011]
    pchol_hs_e       0.001     [0.001, 0.001]        Rtotal_e         0.000547  [-0.001, 0.172]
    pe_hs_e          0.001     [0.001, 0.001]        for_e            0.000529  [-0.001, 0.254]
    pe_hs_r          0.001     [0.001, 0.001]        dem2emgacpai...  0.0005    [0.0005, 0.0005]
    ps_hs_e          0.001     [0.001, 0.001]        gpi_sig_r        0.0005    [0.0005, 0.0005]
    tag_hs_e         0.001     [0.001, 0.001]        h2o2_e           0         [-0.001, 10]
    5mthf_e          0.001     [0, 0.001]            gal_e            0         [-0.001, 2.13]
    Lcystin_e        0.001     [0, 0.001]            abt_e            0         [0, 1.7]
    crn_e            0.001     [0, 0.001]            xylt_e           0         [-0.001, 1.7]
    dcsptn1_e        0.001     [0, 0.001]            oxa_e            0         [0, 1.57]
    fad_e            0.001     [0, 0.001]            glyc_S_e         0         [0, 0.917]
    gtp_e            0.001     [0, 0.001]            HC01444_e        0         [-0.001, 0.85]
    hdca_e           0.001     [0, 0.001]            5oxpro_e         0         [0, 0.581]
    mag_hs_e         0.001     [0, 0.001]            ppi_e            0         [0, 0.503]
    ocdca_e          0.001     [0, 0.001]            nrvnc_e          0         [0, 0.368]
    thymd_e          0.001     [0, 0.001]            urate_e          0         [0, 0.328]
    gdp_e            0.001     [0.001, -0.001]       chsterol_e       0         [0, 0.226]
    gmp_e            0.001     [0.001, -0.002]       C02528_e         0         [0, 0.218]
    cytd_e           0.001     [0.001, -0.00736]     cholate_e        0         [0, 0.218]
    dcmp_e           0.001     [0.001, -0.00736]     fuc_L_e          0         [0, 0.189]
    dctp_n           0.001     [0.001, -0.00736]     xolest2_hs_e     0         [-0.001, 0.17]
    dcyt_e           0.001     [0.001, -0.00736]     4hphac_e         0         [0, 0.159]
    duri_e           0.001     [0.001, -0.00736]     tymsf_e          0         [0, 0.159]
    udp_e            0.001     [0.001, -0.00736]     fucfucgalacg...  0         [0, 0.152]
    ump_e            0.001     [0.001, -0.00736]     galgalfucfuc...  0         [0, 0.0962]
    uri_e            0.001     [0.001, -0.00736]     gthox_e          0         [0, 0.0457]
    utp_e            0.001     [0.001, -0.00736]     citr_L_c         0         [-0.001, 0.0357]
    adn_e            0.001     [0.001, -0.00927]     2hb_e            0         [0, 0.028]
    ahcys_e          0.001     [0.001, -0.00927]     q10h2_e          0         [0, 0.0268]
    atp_e            0.001     [0.001, -0.00927]     din_e            0         [-0.001, 0.0178]
    gsn_e            0.001     [0.001, -0.00927]     anth_c           0         [0, 0.0157]
    xmp_e            0.001     [0.001, -0.00927]     bildglcur_e      0         [0, 0.0154]
    dttp_n           0.001     [0.001, -0.01]        bilglcur_e       0         [-0.001, 0.0144]
    dad_2_e          0.001     [0.001, -0.0178]      co_e             0         [0, 0.0134]
    datp_n           0.001     [0.001, -0.0178]      pheme_e          0         [0, 0.0134]
    dgtp_n           0.001     [0.001, -0.0178]      octdececoa_c     0         [-0.0005, 0.0113]
    orn_e            0.001     [0.001, -0.107]       35cgmp_e         0         [0, 0.0103]
    ocdcea_e         0.001     [0.001, -0.49]        camp_e           0         [0, 0.0103]
    man_e            0.001     [0.001, -1.22]        5mta_e           0         [-0.001, 0.0103]
    hco3_e           0.001     [0.001, -3.98]        spmd_e           0         [-0.001, 0.0103]
    leuktrC4_e       0.000704  [0.001, -0.002]       ins_e            0         [-0.001, 0.00927]
    pre_prot_r       0.0005    [0.0005, 0.0005]      ala_B_e          0         [0, 0.00836]
    arachd_e         0.0005    [0.001, -0.002]       sprm_c           0         [0, 0.00564]
    sph1p_e          0.000453  [0.001, -0.469]       Ser_Thr_l        0         [-0.001, 0.005]
    orot_e           0.000398  [0, 0.001]            adrn_e           0         [-0.001, 0.003]
    pglyc_hs_e       0.000378  [0.000378, 0.000378]  1mncam_e         0         [0, 0.002]
    trp_L_e          0.000345  [0.000345, 0.016]     5adtststeron...  0         [0, 0.002]
    datp_m           0.00033   [0, 0.001]            andrstrnglc_e    0         [0, 0.002]
    ha_e             0.00025   [0.001, -0.418]       estroneglc_e     0         [0, 0.002]
    o2s_e            0.000186  [0, 0.001]            imp_e            0         [0, 0.002]
    arab_L_e         7.4e-05   [0, 0.001]            n5m2masn_g       0         [0, 0.002]
    fol_e            0         [-0.003, 0.05]        nac_e            0         [0, 0.002]
    pnto_R_e         0         [0, 0.0103]           core8_g          0         [0, 0.0012]
    3hpvs_e          0         [0, 0.001]            3hdececrn_e      0         [0, 0.001]
    4hpro_LT_m       0         [0, 0.001]            3hpvstet_e       0         [0, 0.001]
    4nph_e           0         [0, 0.001]            3ivcrn_e         0         [0, 0.001]
    5fthf_e          0         [0, 0.001]            3mlda_e          0         [0, 0.001]
    C02470_e         0         [0, 0.001]            3octdec2crn_e    0         [0, 0.001]
    HC00822_e        0         [0, 0.001]            3octdeccrn_e     0         [0, 0.001]
    acac_e           0         [0, 0.001]            3octdece1crn_e   0         [0, 0.001]
    acmp_e           0         [0, 0.001]            3tetd7ecoacrn_e  0         [0, 0.001]
    arach_e          0         [0, 0.001]            3thexddcoacrn_e  0         [0, 0.001]
    bilirub_e        0         [0, 0.001]            3ttetddcoacrn_e  0         [0, 0.001]
    crvnc_e          0         [0, 0.001]            4nphsf_e         0         [0, 0.001]
    dca_e            0         [0, 0.001]            aprgstrn_e       0         [0, 0.001]
    dheas_e          0         [0, 0.001]            c4crn_e          0         [0, 0.001]
    doco13ac_e       0         [0, 0.001]            c5dc_e           0         [0, 0.001]
    dopa_e           0         [0, 0.001]            c6crn_e          0         [0, 0.001]
    etoh_e           0         [0, 0.001]            clpnd_e          0         [0, 0.001]
    fmn_e            0         [0, 0.001]            dmhptcrn_e       0         [0, 0.001]
    gam_e            0         [0, 0.001]            dopasf_e         0         [0, 0.001]
    gluala_e         0         [0, 0.001]            epoxtac_e        0         [0, 0.001]
    hista_e          0         [0, 0.001]            ivcrn_e          0         [0, 0.001]
    hpdca_e          0         [0, 0.001]            leuktrB4_e       0         [0, 0.001]
    lnlnca_e         0         [0, 0.001]            prostgf2_e       0         [0, 0.001]
    lnlncg_e         0         [0, 0.001]            sulpacmp_e       0         [0, 0.001]
    n2m2nmasn_e      0         [0, 0.001]            triodthysuf_e    0         [0, 0.001]
    octa_e           0         [0, 0.001]            tststeroneglc_e  0         [0, 0.001]
    phyt_e           0         [0, 0.001]            5adtststerone_e  0         [-0.001, 0.001]
    prgstrn_e        0         [0, 0.001]            Asn_X_Ser_Thr_l  0         [-0.001, 0.001]
    prostge2_e       0         [0, 0.001]            dlnlcg_e         0         [-0.001, 0.001]
    ptdca_e          0         [0, 0.001]            estradiol_e      0         [-0.001, 0.001]
    strdnc_e         0         [0, 0.001]            estrones_e       0         [-0.001, 0.001]
    tacr_e           0         [0, 0.001]            fald_e           0         [-0.001, 0.001]
    tmndnc_e         0         [0, 0.001]            idp_e            0         [-0.001, 0.001]
    triodthy_e       0         [0, 0.001]            itp_e            0         [-0.001, 0.001]
    tststerone_e     0         [0, 0.001]            meoh_e           0         [-0.001, 0.001]
    ttdca_e          0         [0, 0.001]            leuktrA4_e       0         [0, -0.001]
    whhdca_e         0         [0, 0.001]
    T_antigen_g      0         [-0.0002, 0.001]
    pmtcoa_r         0         [-0.0005, 0.001]
    andrstrn_e       0         [-0.001, 0.001]
    ncam_c           0         [-0.001, 0.001]



```python
model_KeratinocytesMean2.summary(fva=True)
```

    IN FLUXES                                        OUT FLUXES                                     OBJECTIVES
    -----------------------------------------------  ---------------------------------------------  -----------------------
    id                   Flux  Range                 id                   Flux  Range               biomass_reac...  0.0242
    ---------------  --------  --------------------  ---------------  --------  ------------------
    glc_D_e          0.266     [0.199, 4.5]          h_e              0.499     [-0.63, 10]
    o2_e             0.0589    [0.00246, 8.49]       lac_L_e          0.467     [-1, 9.73]
    val_L_e          0.0198    [0.00855, 0.094]      co2_e            0.0549    [-0.001, 10]
    gln_L_e          0.0146    [-0.16, 0.584]        h2o_e            0.0508    [-1.6, 10]
    lys_L_e          0.0144    [0.0144, 0.0154]      pyr_e            0.0364    [-0.11, 9.86]
    ser_L_e          0.0134    [0.042, -1.12]        ala_B_e          0.00953   [-0.001, 1.3]
    leu_L_e          0.0132    [0.0132, 0.0142]      orn_e            0.00855   [0, 0.0865]
    gly_e            0.0121    [0.03, -1.14]         glyc_e           0.00696   [-0.001, 2.73]
    pro_L_e          0.01      [0.01, 0.01]          pi_e             0.00297   [0.0155, -0.893]
    glu_L_e          0.01      [0.01, -0.734]        Rtotal_e         0.002     [-0.001, 0.249]
    arg_L_e          0.00871   [0.00671, 0.00871]    dad_2_e          0.00145   [-0.001, 0.0188]
    thr_L_e          0.00758   [0.00758, 0.00758]    ala_L_e          0.00137   [-0.01, 1.16]
    ile_L_e          0.00694   [0.00694, 0.105]      gua_e            0.00112   [-0.001, 0.342]
    asn_L_e          0.00677   [0.00677, -0.737]     elaid_e          0.001     [-0.001, 0.467]
    phe_L_e          0.00629   [0.00629, 0.066]      prostgh2_e       0.001     [0, 0.002]
    tyr_L_e          0.00387   [-0.0558, 0.104]      3ivcrn_e         0.001     [0, 0.001]
    met_L_e          0.00371   [0.00371, 0.03]       lac_D_e          0.001     [0, 0.001]
    chol_e           0.00217   [0.05, -0.128]        dcsptn1_e        0.001     [0.001, -0.001]
    his_L_e          0.00206   [0.00206, 0.00307]    dgsn_e           0.00076   [-0.001, 0.0188]
    biomass_other_c  0.00131   [0.00131, 0.00131]    for_e            0.000495  [-0.001, 0.395]
    cys_L_e          0.00106   [-0.0262, 0.0626]     5mta_e           0.000495  [0, 0.001]
    pe_hs_e          0.001     [0.000342, 0.001]     sprm_c           0.000495  [0, 0.001]
    adrn_e           0.001     [0, 0.001]            2hb_e            7.3e-05   [0, 0.0273]
    arach_e          0.001     [0, 0.001]            glyc_S_e         0         [0, 2.73]
    carn_e           0.001     [0, 0.001]            gal_e            0         [0, 2.15]
    crn_e            0.001     [0, 0.001]            xylt_e           0         [-0.001, 1.72]
    o2s_e            0.001     [0, 0.001]            acmana_e         0         [0, 1.32]
    orot_e           0.001     [0, 0.001]            taur_e           0         [-0.001, 1.09]
    pchol_hs_e       0.001     [0, 0.001]            fe2_e            0         [-1, 1]
    prostge2_e       0.001     [0, 0.001]            ha_pre1_e        0         [-0.001, 0.849]
    tag_hs_e         0.001     [0, 0.001]            ethamp_r         0         [0, 0.846]
    dag_hs_e         0.001     [-0.000705, 0.001]    5oxpro_e         0         [0, 0.745]
    glyc3p_e         0.001     [-0.000705, 0.001]    vacc_e           0         [-0.001, 0.467]
    lpchol_hs_e      0.001     [-0.000705, 0.001]    urate_e          0         [0, 0.429]
    ps_hs_e          0.001     [-0.000705, 0.001]    hxan_e           0         [-0.001, 0.428]
    cdp_e            0.001     [0.001, -0.002]       ha_e             0         [-0.001, 0.424]
    cmp_e            0.001     [0.001, -0.002]       chsterol_e       0         [0, 0.391]
    ctp_e            0.001     [0.001, -0.002]       cholate_e        0         [0, 0.378]
    cytd_e           0.001     [0.001, -0.00753]     sph1p_e          0         [-0.001, 0.37]
    dctp_n           0.001     [0.001, -0.00753]     eicostet_e       0         [-0.001, 0.355]
    udp_e            0.001     [0.001, -0.00753]     tdchola_e        0         [0, 0.354]
    ump_e            0.001     [0.001, -0.00753]     ade_e            0         [-0.001, 0.347]
    ura_e            0.001     [0.001, -0.00753]     nrvnc_e          0         [0, 0.333]
    uri_e            0.001     [0.001, -0.00753]     xolest2_hs_e     0         [-0.001, 0.248]
    utp_e            0.001     [0.001, -0.00753]     glyb_e           0         [0, 0.177]
    adn_e            0.001     [0.001, -0.0136]      4hphac_e         0         [0, 0.16]
    aicar_e          0.001     [0.001, -0.0136]      acac_e           0         [-0.001, 0.16]
    atp_e            0.001     [0.001, -0.0136]      spc_hs_e         0         [0, 0.14]
    xmp_e            0.001     [0.001, -0.0136]      q10h2_e          0         [0, 0.0918]
    datp_n           0.001     [0.001, -0.0143]      q10_e            0         [-0.001, 0.0903]
    dgtp_n           0.001     [0.001, -0.0188]      cgly_e           0         [-0.001, 0.0898]
    strdnc_e         0.001     [0.001, -0.388]       gthrd_e          0         [-0.001, 0.0898]
    cit_e            0.001     [0.001, -0.743]       pheacgln_e       0         [0, 0.0597]
    hdca_e           0.000892  [0.001, -0.515]       din_e            0         [-0.001, 0.0188]
    ahcys_e          0.000568  [0.001, -0.0136]      35cgmp_e         0         [0, 0.0146]
    inost_e          0.000565  [0.000565, 0.000894]  camp_e           0         [0, 0.0146]
    spmd_e           0.000495  [0, 0.001]            gsn_e            0         [-0.001, 0.0136]
    sphs1p_e         0.000424  [0.001, -0.37]        ins_e            0         [-0.001, 0.0136]
    pglyc_hs_e       0.000353  [0.000353, 0.000353]  bilglcur_e       0         [0, 0.0108]
    pe_hs_r          0.000342  [0.000342, 0.001]     co_e             0         [0, 0.0108]
    trp_L_e          0.000323  [0.000323, 0.000323]  pheme_e          0         [0, 0.0108]
    dttp_n           0.000317  [0.001, -0.00821]     dttp_m           0         [0, 0.00921]
    HC00342_e        0.000233  [0.001, -0.743]       Ser_Thr_l        0         [-0.001, 0.005]
    so4_e            0         [-0.0898, 1]          fuc14galacgl...  0         [0, 0.004]
    fe3_e            0         [-1, 1]               fuc_L_e          0         [0, 0.004]
    asp_L_e          0         [0, 0.01]             galfuc12gal1...  0         [0, 0.004]
    3hpvs_e          0         [0, 0.001]            galfucgalacg...  0         [0, 0.004]
    4hpro_LT_m       0         [0, 0.001]            fald_e           0         [-0.001, 0.003]
    5HPET_c          0         [0, 0.001]            man_e            0         [-0.001, 0.003]
    5fthf_e          0         [0, 0.001]            meoh_e           0         [-0.001, 0.003]
    HC00822_e        0         [0, 0.001]            34dhoxpeg_e      0         [0, 0.002]
    acmp_e           0         [0, 0.001]            5adtststeron...  0         [0, 0.002]
    allop_e          0         [0, 0.001]            andrstrnglc_e    0         [0, 0.002]
    amp_e            0         [0, 0.001]            imp_e            0         [0, 0.002]
    arab_L_e         0         [0, 0.001]            leuktrB4_e       0         [0, 0.002]
    c81coa_c         0         [0, 0.001]            leuktrC4_e       0         [0, 0.002]
    crvnc_e          0         [0, 0.001]            gmp_e            0         [-0.001, 0.002]
    dca_e            0         [0, 0.001]            octdececoa_c     0         [-0.001, 0.002]
    dd2coa_c         0         [0, 0.001]            pmtcoa_r         0         [-0.001, 0.002]
    decdicoa_c       0         [0, 0.001]            fucfucfucgal...  0         [0, 0.00133]
    dheas_e          0         [0, 0.001]            n5m2masn_g       0         [0, 0.00133]
    doco13ac_e       0         [0, 0.001]            core5_g          0         [0, 0.0012]
    dopa_e           0         [0, 0.001]            core7_g          0         [0, 0.0012]
    etoh_e           0         [0, 0.001]            core8_g          0         [0, 0.0012]
    fmn_e            0         [0, 0.001]            3bcrn_e          0         [0, 0.001]
    gam_e            0         [0, 0.001]            3ddcrn_e         0         [0, 0.001]
    gluala_e         0         [0, 0.001]            3deccrn_e        0         [0, 0.001]
    gtp_e            0         [0, 0.001]            3hdececrn_e      0         [0, 0.001]
    h2o2_e           0         [0, 0.001]            3hexdcrn_e       0         [0, 0.001]
    hista_e          0         [0, 0.001]            3hpvstet_e       0         [0, 0.001]
    hpdca_e          0         [0, 0.001]            3mlda_e          0         [0, 0.001]
    lgnc_e           0         [0, 0.001]            3octdec2crn_e    0         [0, 0.001]
    lnlc_e           0         [0, 0.001]            3octdeccrn_e     0         [0, 0.001]
    n2m2nmasn_e      0         [0, 0.001]            3octdece1crn_e   0         [0, 0.001]
    ocdca_e          0         [0, 0.001]            3tdcrn_e         0         [0, 0.001]
    octa_e           0         [0, 0.001]            3tetd7ecoacrn_e  0         [0, 0.001]
    ppa_e            0         [0, 0.001]            3thexddcoacrn_e  0         [0, 0.001]
    prgstrn_e        0         [0, 0.001]            3ttetddcoacrn_e  0         [0, 0.001]
    ptdca_e          0         [0, 0.001]            abt_e            0         [0, 0.001]
    ptvst_e          0         [0, 0.001]            aprgstrn_e       0         [0, 0.001]
    retinol_e        0         [0, 0.001]            c4crn_e          0         [0, 0.001]
    rsv_e            0         [0, 0.001]            c51crn_e         0         [0, 0.001]
    tacr_e           0         [0, 0.001]            c5dc_e           0         [0, 0.001]
    tetdece1coa_c    0         [0, 0.001]            c6crn_e          0         [0, 0.001]
    tmndnc_e         0         [0, 0.001]            c81crn_e         0         [0, 0.001]
    tststerone_e     0         [0, 0.001]            c8crn_e          0         [0, 0.001]
    ttdca_e          0         [0, 0.001]            clpnd_e          0         [0, 0.001]
    whhdca_e         0         [0, 0.001]            ddece1crn_e      0         [0, 0.001]
    T_antigen_g      0         [-0.0002, 0.001]      decdicrn_e       0         [0, 0.001]
    andrstrn_e       0         [-0.001, 0.001]       epoxtac_e        0         [0, 0.001]
    nrpphr_e         0         [-0.001, 0.001]       fol_e            0         [0, 0.001]
    pre_prot_r       0         [0, 0.000329]         gthox_e          0         [0, 0.001]
    hdcea_e          0         [0.001, -0.529]       ivcrn_e          0         [0, 0.001]
    taur_c           0         [0.001, -1.09]        ndersv_e         0         [0, 0.001]
    oxa_e            0         [0, -2.17]            oxyp_e           0         [0, 0.001]
                                                     prostgi2_e       0         [0, 0.001]
                                                     ptvstm3_e        0         [0, 0.001]
                                                     retn_e           0         [0, 0.001]
                                                     ribflv_e         0         [0, 0.001]
                                                     rsvlac_e         0         [0, 0.001]
                                                     sulpacmp_e       0         [0, 0.001]
                                                     tetdece1crn_e    0         [0, 0.001]
                                                     tststeroneglc_e  0         [0, 0.001]
                                                     ttdcrn_e         0         [0, 0.001]
                                                     Asn_X_Ser_Thr_l  0         [-0.000333, 0.001]
                                                     5adtststerone_e  0         [-0.001, 0.001]
                                                     citr_L_c         0         [-0.001, 0.001]
                                                     citr_L_e         0         [-0.001, 0.001]
                                                     gdp_e            0         [-0.001, 0.001]
                                                     idp_e            0         [-0.001, 0.001]
                                                     itp_e            0         [-0.001, 0.001]
                                                     leuktrA4_e       0         [-0.001, 0.001]
                                                     dem2emgacpai...  0         [0, 0.000329]
                                                     gpi_sig_r        0         [0, 0.000329]
                                                     mem2emgacpai...  0         [0, 0.000329]
                                                     dgpi_prot_hs_r   0         [0, 0.000219]



```python
model_HeLaMean2.summary(fva=True)
```

    IN FLUXES                                        OUT FLUXES                                       OBJECTIVES
    -----------------------------------------------  -----------------------------------------------  -----------------------
    id                   Flux  Range                 id                   Flux  Range                 biomass_reac...  0.0266
    ---------------  --------  --------------------  ---------------  --------  --------------------
    glc_D_e          0.304     [0.255, 4.5]          h_e              0.572     [-0.716, 10]
    o2_e             0.0722    [0.00272, 10]         lac_L_e          0.499     [-1, 9.63]
    ser_L_e          0.042     [0.042, -1.22]        pyr_e            0.059     [-0.11, 10]
    arg_L_e          0.0311    [0, 0.084]            co2_e            0.055     [-0.001, 6.61]
    val_L_e          0.0243    [0.00938, 0.094]      glyc_e           0.0346    [-0.001, 1.26]
    phe_L_e          0.018     [0.0069, 0.066]       ala_L_e          0.0285    [-0.01, 1.35]
    lys_L_e          0.0166    [0.0158, 0.0168]      h2o_e            0.0266    [-5.28, 10]
    leu_L_e          0.0145    [0.0145, 0.0155]      urea_e           0.0216    [0, 0.111]
    gln_L_e          0.0133    [-0.279, 0.584]       ptrc_e           0.0206    [-0.001, 0.183]
    pro_L_e          0.01      [0.01, -0.0965]       ac_e             0.00842   [-0.001, 0.0679]
    gly_e            0.0097    [0.03, -1.7]          tyr_L_e          0.00683   [0.0549, -0.104]
    thr_L_e          0.00832   [0.00832, 0.00832]    gua_e            0.00578   [-0.001, 0.357]
    ile_L_e          0.00761   [0.00761, 0.105]      so4_e            0.0056    [0.0933, -0.164]
    asn_L_e          0.00743   [0.01, -0.853]        adn_e            0.00424   [-0.001, 0.00915]
    cys_L_e          0.00683   [-0.0287, 0.0626]     ala_B_e          0.003     [0, 0.00929]
    met_L_e          0.00407   [0.00207, 0.03]       lnlc_e           0.001     [-0.001, 0.509]
    pi_e             0.00323   [-0.00605, 1]         vacc_e           0.001     [-0.001, 0.506]
    chol_e           0.00305   [0.00305, 0.05]       cholate_e        0.001     [0, 0.219]
    his_L_e          0.00236   [0.00236, 0.00336]    5mta_e           0.001     [-0.001, 0.0102]
    biomass_other_c  0.00144   [0.00144, 0.00144]    spmd_e           0.001     [-0.001, 0.0102]
    dag_hs_e         0.001     [0.001, 0.001]        thym_e           0.001     [-0.001, 0.00962]
    lpchol_hs_e      0.001     [0.001, 0.001]        ribflv_e         0.001     [0, 0.002]
    pchol_hs_e       0.001     [0.001, 0.001]        leuktrD4_e       0.001     [-0.001, 0.002]
    pe_hs_e          0.001     [0.001, 0.001]        lac_D_e          0.001     [0, 0.001]
    pe_hs_r          0.001     [0.001, 0.001]        c5dc_e           0.000803  [0, 0.001]
    ps_hs_e          0.001     [0.001, 0.001]        dcmp_e           0.000749  [-0.001, 0.00729]
    tag_hs_e         0.001     [0.001, 0.001]        dttp_m           0.000652  [0, 0.0109]
    4hpro_LT_m       0.001     [0, 0.001]            for_e            0.000543  [-0.001, 0.253]
    carn_e           0.001     [0, 0.001]            uri_e            0.00054   [-0.001, 0.00729]
    crn_e            0.001     [0, 0.001]            1mncam_e         0.000419  [0, 0.002]
    dcsptn1_e        0.001     [0, 0.001]            abt_e            0.000333  [0, 1.7]
    fad_e            0.001     [0, 0.001]            dgpi_prot_hs_r   0.000333  [0.000333, 0.000333]
    gtp_e            0.001     [0, 0.001]            gpi_sig_r        0.000333  [0.000333, 0.000333]
    o2s_e            0.001     [0, 0.001]            hxan_e           0.000311  [-0.001, 0.412]
    orot_e           0.001     [0, 0.001]            3ivcrn_e         0.000197  [0, 0.001]
    strdnc_e         0.001     [0, 0.001]            gal_e            0         [0, 2.12]
    thymd_e          0.001     [0, 0.001]            xylt_e           0         [-0.001, 1.7]
    dlnlcg_e         0.001     [0.001, -0.001]       oxa_e            0         [0, 1.58]
    gdp_e            0.001     [0.001, -0.001]       man_e            0         [-0.001, 1.21]
    gmp_e            0.001     [0.001, -0.002]       glyc_S_e         0         [0, 1.2]
    leuktrC4_e       0.001     [0.001, -0.002]       fe2_e            0         [-1, 1]
    cytd_e           0.001     [0.001, -0.00729]     fe3_e            0         [-1, 1]
    dctp_n           0.001     [0.001, -0.00729]     5oxpro_e         0         [0, 0.938]
    dcyt_e           0.001     [0.001, -0.00729]     ppi_e            0         [0, 0.503]
    duri_e           0.001     [0.001, -0.00729]     ocdcea_e         0         [-0.001, 0.489]
    utp_e            0.001     [0.001, -0.00729]     urate_e          0         [0, 0.415]
    ahcys_e          0.001     [0.001, -0.00915]     nrvnc_e          0         [0, 0.367]
    atp_e            0.001     [0.001, -0.00915]     chsterol_e       0         [0, 0.225]
    gsn_e            0.001     [0.001, -0.00915]     C02528_e         0         [0, 0.217]
    ins_e            0.001     [0.001, -0.00915]     fuc_L_e          0         [0, 0.189]
    xmp_e            0.001     [0.001, -0.00915]     Rtotal_e         0         [-0.001, 0.172]
    dttp_n           0.001     [0.001, -0.00994]     xolest2_hs_e     0         [-0.001, 0.17]
    dad_2_e          0.001     [0.001, -0.0186]      4hphac_e         0         [0, 0.159]
    datp_n           0.001     [0.001, -0.0186]      tymsf_e          0         [0, 0.159]
    dgsn_e           0.001     [0.001, -0.0186]      fucfucgalacg...  0         [0, 0.152]
    dgtp_n           0.001     [0.001, -0.0186]      orn_e            0         [-0.001, 0.107]
    din_e            0.001     [0.001, -0.0186]      galgalfucfuc...  0         [0, 0.0961]
    gchola_e         0.001     [0.001, -0.216]       fol_c            0         [-0.001, 0.052]
    h2o2_e           0.001     [0.001, -10]          glyb_e           0         [0, 0.047]
    inost_e          0.000953  [0.000953, -2.83]     gthox_e          0         [0, 0.0456]
    hco3_e           0.000607  [0.001, -3.97]        citr_L_c         0         [-0.001, 0.0355]
    hdca_e           0.000465  [0, 0.001]            2hb_e            0         [0, 0.0279]
    sph1p_e          0.000465  [0.001, -0.498]       q10h2_e          0         [0, 0.0266]
    ncam_c           0.000419  [-0.001, 0.001]       anth_c           0         [0, 0.0156]
    arachd_e         0.000393  [0.001, -0.002]       bildglcur_e      0         [0, 0.0153]
    pglyc_hs_e       0.000388  [0.000388, 0.000388]  bilglcur_e       0         [-0.001, 0.0143]
    trp_L_e          0.000354  [0.000354, 0.016]     co_e             0         [0, 0.0133]
    pre_prot_r       0.000333  [0.000333, 0.000333]  pheme_e          0         [0, 0.0133]
    datp_m           0.000326  [0, 0.001]            35cgmp_e         0         [0, 0.0102]
    ha_e             0.000167  [0.001, -0.418]       camp_e           0         [0, 0.0102]
    fol_e            0         [-0.003, 0.05]        udp_e            0         [-0.001, 0.00729]
    pnto_R_e         0         [0, 0.0102]           ump_e            0         [-0.001, 0.00729]
    asp_L_e          0         [0, 0.01]             ura_e            0         [-0.001, 0.00729]
    glu_L_e          0         [-0.001, 0.01]        sprm_c           0         [0, 0.00558]
    3hpvs_e          0         [0, 0.001]            Ser_Thr_l        0         [-0.001, 0.005]
    4nph_e           0         [0, 0.001]            adrn_e           0         [-0.001, 0.003]
    5fthf_e          0         [0, 0.001]            5adtststeron...  0         [0, 0.002]
    5mthf_e          0         [0, 0.001]            andrstrnglc_e    0         [0, 0.002]
    C02470_e         0         [0, 0.001]            estroneglc_e     0         [0, 0.002]
    HC00822_e        0         [0, 0.001]            imp_e            0         [0, 0.002]
    Lcystin_e        0         [0, 0.001]            n5m2masn_g       0         [0, 0.002]
    acac_e           0         [0, 0.001]            nac_e            0         [0, 0.002]
    acmp_e           0         [0, 0.001]            13_cis_retn_n    0         [-0.001, 0.002]
    arab_L_e         0         [0, 0.001]            retn_e           0         [-0.001, 0.002]
    arach_e          0         [0, 0.001]            core8_g          0         [0, 0.0012]
    bilirub_e        0         [0, 0.001]            3bcrn_e          0         [0, 0.001]
    crvnc_e          0         [0, 0.001]            3hdececrn_e      0         [0, 0.001]
    dca_e            0         [0, 0.001]            3hpvstet_e       0         [0, 0.001]
    dheas_e          0         [0, 0.001]            3mlda_e          0         [0, 0.001]
    doco13ac_e       0         [0, 0.001]            3octdec2crn_e    0         [0, 0.001]
    dopa_e           0         [0, 0.001]            3octdeccrn_e     0         [0, 0.001]
    etoh_e           0         [0, 0.001]            3octdece1crn_e   0         [0, 0.001]
    fmn_e            0         [0, 0.001]            3tetd7ecoacrn_e  0         [0, 0.001]
    gam_e            0         [0, 0.001]            3thexddcoacrn_e  0         [0, 0.001]
    gluala_e         0         [0, 0.001]            3ttetddcoacrn_e  0         [0, 0.001]
    hista_e          0         [0, 0.001]            4nphsf_e         0         [0, 0.001]
    hpdca_e          0         [0, 0.001]            aprgstrn_e       0         [0, 0.001]
    leuktrA4_e       0         [0, 0.001]            c4crn_e          0         [0, 0.001]
    lnlnca_e         0         [0, 0.001]            c6crn_e          0         [0, 0.001]
    lnlncg_e         0         [0, 0.001]            clpnd_e          0         [0, 0.001]
    mag_hs_e         0         [0, 0.001]            dmhptcrn_e       0         [0, 0.001]
    n2m2nmasn_e      0         [0, 0.001]            dopasf_e         0         [0, 0.001]
    ocdca_e          0         [0, 0.001]            epoxtac_e        0         [0, 0.001]
    octa_e           0         [0, 0.001]            ivcrn_e          0         [0, 0.001]
    phyt_e           0         [0, 0.001]            leuktrB4_e       0         [0, 0.001]
    prgstrn_e        0         [0, 0.001]            prostgf2_e       0         [0, 0.001]
    prostge2_e       0         [0, 0.001]            sulpacmp_e       0         [0, 0.001]
    ptdca_e          0         [0, 0.001]            triodthysuf_e    0         [0, 0.001]
    retinol_e        0         [0, 0.001]            tststeroneglc_e  0         [0, 0.001]
    tacr_e           0         [0, 0.001]            Asn_X_Ser_Thr_l  0         [-0.001, 0.001]
    tmndnc_e         0         [0, 0.001]            estradiol_e      0         [-0.001, 0.001]
    triodthy_e       0         [0, 0.001]            estrones_e       0         [-0.001, 0.001]
    tststerone_e     0         [0, 0.001]            fald_e           0         [-0.001, 0.001]
    ttdca_e          0         [0, 0.001]            idp_e            0         [-0.001, 0.001]
    whhdca_e         0         [0, 0.001]            itp_e            0         [-0.001, 0.001]
    T_antigen_g      0         [-0.0002, 0.001]      meoh_e           0         [-0.001, 0.001]
    pmtcoa_r         0         [-0.000333, 0.001]    taur_c           0         [-0.001, 0.001]
    5adtststerone_e  0         [-0.001, 0.001]       tchola_e         0         [-0.001, 0.001]
    andrstrn_e       0         [-0.001, 0.001]
    eicostet_e       0         [-0.001, 0.001]
    octdececoa_c     0         [0.000333, -0.0112]
    hdcea_e          0         [0.001, -0.57]



```python
model_KeratinocytesMean3.summary(fva=True)
```

    IN FLUXES                                        OUT FLUXES                                     OBJECTIVES
    -----------------------------------------------  ---------------------------------------------  -----------------------
    id                   Flux  Range                 id                   Flux  Range               biomass_reac...  0.0242
    ---------------  --------  --------------------  ---------------  --------  ------------------
    glc_D_e          0.274     [0.203, 4.5]          h_e              0.524     [-0.619, 10]
    o2_e             0.0374    [0.00246, 8.68]       lac_L_e          0.478     [-1, 9.73]
    ser_L_e          0.0281    [0.042, -1.12]        pyr_e            0.0439    [-0.11, 9.86]
    val_L_e          0.0146    [0.00855, 0.094]      h2o_e            0.0406    [-1.54, 10]
    lys_L_e          0.0144    [0.0144, 0.0154]      co2_e            0.027     [-0.001, 10]
    leu_L_e          0.0132    [0.0132, 0.0142]      glyc_e           0.021     [-0.001, 1.6]
    gly_e            0.0111    [0.03, -1.14]         ura_e            0.00434   [-0.001, 0.00753]
    pro_L_e          0.01      [0.01, 0.01]          dgsn_e           0.00355   [-0.001, 0.0188]
    arg_L_e          0.00871   [0.00671, 0.00871]    elaid_e          0.003     [-0.001, 0.461]
    gln_L_e          0.0079    [-0.251, 0.584]       Rtotal_e         0.002     [-0.001, 0.247]
    thr_L_e          0.00758   [0.00758, 0.00758]    gua_e            0.00137   [-0.001, 0.342]
    glu_L_e          0.00736   [0.01, -1.08]         ahcys_e          0.00101   [-0.001, 0.0136]
    ile_L_e          0.00694   [0.00694, 0.105]      ala_B_e          0.001     [-0.001, 1.31]
    asn_L_e          0.00677   [0.00677, -0.828]     leuktrB4_e       0.001     [0, 0.002]
    phe_L_e          0.00629   [0.00629, 0.066]      prostgh2_e       0.001     [0, 0.002]
    met_L_e          0.00472   [0.00371, 0.03]       epoxtac_e        0.001     [0, 0.001]
    tyr_L_e          0.00387   [-0.0558, 0.104]      lac_D_e          0.001     [0, 0.001]
    asp_L_e          0.0025    [0, 0.01]             ndersv_e         0.001     [0, 0.001]
    his_L_e          0.00206   [0.00206, 0.00307]    for_e            0.000495  [-0.001, 0.395]
    biomass_other_c  0.00131   [0.00131, 0.00131]    hdcea_e          0.000423  [-0.001, 0.522]
    cys_L_e          0.00113   [-0.0262, 0.0626]     dad_2_e          0.000228  [-0.001, 0.0188]
    4hpro_LT_m       0.001     [0, 0.001]            oxa_e            0         [0, 2.16]
    adrn_e           0.001     [0, 0.001]            gal_e            0         [0, 2.15]
    arach_e          0.001     [0, 0.001]            xylt_e           0         [-0.001, 1.72]
    carn_e           0.001     [0, 0.001]            glyc_S_e         0         [0, 1.6]
    o2s_e            0.001     [0, 0.001]            acmana_e         0         [0, 1.32]
    pchol_hs_e       0.001     [0, 0.001]            ala_L_e          0         [-0.01, 1.16]
    prostge2_e       0.001     [0, 0.001]            taur_c           0         [-0.001, 1.09]
    rsv_e            0.001     [0, 0.001]            taur_e           0         [-0.001, 1.09]
    tacr_e           0.001     [0, 0.001]            5oxpro_e         0         [0, 1.09]
    tag_hs_e         0.001     [0, 0.001]            fe2_e            0         [-1, 1]
    dag_hs_e         0.001     [-0.000924, 0.001]    ethamp_r         0         [0, 0.852]
    glyc3p_e         0.001     [-0.000924, 0.001]    ha_pre1_e        0         [-0.001, 0.848]
    lpchol_hs_e      0.001     [-0.000924, 0.001]    hdca_e           0         [-0.001, 0.511]
    pe_hs_e          0.001     [-0.000924, 0.001]    vacc_e           0         [-0.001, 0.461]
    pe_hs_r          0.001     [-0.000924, 0.001]    urate_e          0         [0, 0.429]
    ps_hs_e          0.001     [-0.000924, 0.001]    ha_e             0         [-0.001, 0.424]
    dcsptn1_e        0.001     [-0.001, 0.001]       chsterol_e       0         [0, 0.391]
    leuktrA4_e       0.001     [0.001, -0.001]       cholate_e        0         [0, 0.378]
    gmp_e            0.001     [0.001, -0.002]       sph1p_e          0         [-0.001, 0.367]
    cytd_e           0.001     [0.001, -0.00753]     tdchola_e        0         [0, 0.354]
    dctp_n           0.001     [0.001, -0.00753]     eicostet_e       0         [-0.001, 0.352]
    udp_e            0.001     [0.001, -0.00753]     ade_e            0         [-0.001, 0.347]
    ump_e            0.001     [0.001, -0.00753]     nrvnc_e          0         [0, 0.329]
    uri_e            0.001     [0.001, -0.00753]     xolest2_hs_e     0         [-0.001, 0.245]
    utp_e            0.001     [0.001, -0.00753]     4hphac_e         0         [0, 0.16]
    aicar_e          0.001     [0.001, -0.0136]      acac_e           0         [-0.001, 0.16]
    atp_e            0.001     [0.001, -0.0136]      glyb_e           0         [0, 0.128]
    ins_e            0.001     [0.001, -0.0136]      spc_hs_e         0         [0, 0.102]
    xmp_e            0.001     [0.001, -0.0136]      q10h2_e          0         [0, 0.0904]
    datp_n           0.001     [0.001, -0.0143]      cgly_e           0         [-0.001, 0.0898]
    dgtp_n           0.001     [0.001, -0.0188]      gthrd_e          0         [-0.001, 0.0898]
    din_e            0.001     [0.001, -0.0188]      q10_e            0         [-0.001, 0.089]
    strdnc_e         0.001     [0.001, -0.385]       orn_e            0         [0, 0.0865]
    HC00342_e        0.001     [0.001, -1.09]        pheacgln_e       0         [0, 0.0597]
    cit_e            0.001     [0.001, -1.09]        2hb_e            0         [0, 0.0273]
    amp_e            0.000861  [0, 0.001]            35cgmp_e         0         [0, 0.0146]
    inost_e          0.000565  [0.000565, 0.00156]   camp_e           0         [0, 0.0146]
    cmp_e            0.000511  [0.001, -0.002]       adn_e            0         [-0.001, 0.0136]
    sphs1p_e         0.000424  [0.001, -0.367]       gsn_e            0         [-0.001, 0.0136]
    pglyc_hs_e       0.000353  [0.000353, 0.000353]  bilglcur_e       0         [0, 0.0108]
    trp_L_e          0.000323  [0.000323, 0.000323]  co_e             0         [0, 0.0108]
    dttp_n           0.000317  [0.001, -0.00821]     pheme_e          0         [0, 0.0108]
    orot_e           0.000296  [0, 0.001]            dttp_m           0         [0, 0.00921]
    hxan_e           3.5e-05   [0.001, -0.428]       Ser_Thr_l        0         [-0.001, 0.005]
    so4_e            0         [-0.0898, 1]          fuc14galacgl...  0         [0, 0.004]
    fe3_e            0         [-1, 1]               fuc_L_e          0         [0, 0.004]
    pi_e             0         [-0.0162, 0.896]      galfuc12gal1...  0         [0, 0.004]
    3hpvs_e          0         [0, 0.001]            galfucgalacg...  0         [0, 0.004]
    5HPET_c          0         [0, 0.001]            man_e            0         [-0.001, 0.003]
    5fthf_e          0         [0, 0.001]            meoh_e           0         [-0.001, 0.003]
    HC00822_e        0         [0, 0.001]            34dhoxpeg_e      0         [0, 0.002]
    acmp_e           0         [0, 0.001]            5adtststeron...  0         [0, 0.002]
    allop_e          0         [0, 0.001]            andrstrnglc_e    0         [0, 0.002]
    arab_L_e         0         [0, 0.001]            imp_e            0         [0, 0.002]
    c81coa_c         0         [0, 0.001]            leuktrC4_e       0         [0, 0.002]
    crn_e            0         [0, 0.001]            cdp_e            0         [-0.001, 0.002]
    crvnc_e          0         [0, 0.001]            ctp_e            0         [-0.001, 0.002]
    dca_e            0         [0, 0.001]            pmtcoa_r         0         [-0.001, 0.002]
    dd2coa_c         0         [0, 0.001]            fucfucfucgal...  0         [0, 0.00133]
    decdicoa_c       0         [0, 0.001]            n5m2masn_g       0         [0, 0.00133]
    dheas_e          0         [0, 0.001]            core5_g          0         [0, 0.0012]
    doco13ac_e       0         [0, 0.001]            core7_g          0         [0, 0.0012]
    dopa_e           0         [0, 0.001]            core8_g          0         [0, 0.0012]
    etoh_e           0         [0, 0.001]            3bcrn_e          0         [0, 0.001]
    fmn_e            0         [0, 0.001]            3ddcrn_e         0         [0, 0.001]
    gam_e            0         [0, 0.001]            3deccrn_e        0         [0, 0.001]
    gluala_e         0         [0, 0.001]            3hdececrn_e      0         [0, 0.001]
    gtp_e            0         [0, 0.001]            3hexdcrn_e       0         [0, 0.001]
    h2o2_e           0         [0, 0.001]            3hpvstet_e       0         [0, 0.001]
    hista_e          0         [0, 0.001]            3ivcrn_e         0         [0, 0.001]
    hpdca_e          0         [0, 0.001]            3mlda_e          0         [0, 0.001]
    lgnc_e           0         [0, 0.001]            3octdec2crn_e    0         [0, 0.001]
    lnlc_e           0         [0, 0.001]            3octdeccrn_e     0         [0, 0.001]
    n2m2nmasn_e      0         [0, 0.001]            3octdece1crn_e   0         [0, 0.001]
    ocdca_e          0         [0, 0.001]            3tdcrn_e         0         [0, 0.001]
    octa_e           0         [0, 0.001]            3tetd7ecoacrn_e  0         [0, 0.001]
    ppa_e            0         [0, 0.001]            3thexddcoacrn_e  0         [0, 0.001]
    pre_prot_r       0         [0, 0.001]            3ttetddcoacrn_e  0         [0, 0.001]
    prgstrn_e        0         [0, 0.001]            5mta_e           0         [0, 0.001]
    ptdca_e          0         [0, 0.001]            abt_e            0         [0, 0.001]
    ptvst_e          0         [0, 0.001]            aprgstrn_e       0         [0, 0.001]
    retinol_e        0         [0, 0.001]            c4crn_e          0         [0, 0.001]
    spmd_e           0         [0, 0.001]            c51crn_e         0         [0, 0.001]
    tetdece1coa_c    0         [0, 0.001]            c5dc_e           0         [0, 0.001]
    tmndnc_e         0         [0, 0.001]            c6crn_e          0         [0, 0.001]
    tststerone_e     0         [0, 0.001]            c81crn_e         0         [0, 0.001]
    ttdca_e          0         [0, 0.001]            c8crn_e          0         [0, 0.001]
    whhdca_e         0         [0, 0.001]            clpnd_e          0         [0, 0.001]
    T_antigen_g      0         [-0.0002, 0.001]      ddece1crn_e      0         [0, 0.001]
    andrstrn_e       0         [-0.001, 0.001]       decdicrn_e       0         [0, 0.001]
    nrpphr_e         0         [-0.001, 0.001]       dem2emgacpai...  0         [0, 0.001]
    octdececoa_c     0         [0.001, -0.002]       dgpi_prot_hs_r   0         [0, 0.001]
    fald_e           0         [0.001, -0.003]       fol_e            0         [0, 0.001]
                                                     gpi_sig_r        0         [0, 0.001]
                                                     gthox_e          0         [0, 0.001]
                                                     ivcrn_e          0         [0, 0.001]
                                                     mem2emgacpai...  0         [0, 0.001]
                                                     oxyp_e           0         [0, 0.001]
                                                     prostgi2_e       0         [0, 0.001]
                                                     ptvstm3_e        0         [0, 0.001]
                                                     retn_e           0         [0, 0.001]
                                                     ribflv_e         0         [0, 0.001]
                                                     rsvlac_e         0         [0, 0.001]
                                                     sprm_c           0         [0, 0.001]
                                                     sulpacmp_e       0         [0, 0.001]
                                                     tetdece1crn_e    0         [0, 0.001]
                                                     tststeroneglc_e  0         [0, 0.001]
                                                     ttdcrn_e         0         [0, 0.001]
                                                     Asn_X_Ser_Thr_l  0         [-0.000333, 0.001]
                                                     5adtststerone_e  0         [-0.001, 0.001]
                                                     citr_L_c         0         [-0.001, 0.001]
                                                     citr_L_e         0         [-0.001, 0.001]
                                                     gdp_e            0         [-0.001, 0.001]
                                                     idp_e            0         [-0.001, 0.001]
                                                     itp_e            0         [-0.001, 0.001]



```python
model_HeLaMean3.summary(fva=True)
```

    IN FLUXES                                        OUT FLUXES                                       OBJECTIVES
    -----------------------------------------------  -----------------------------------------------  -----------------------
    id                   Flux  Range                 id                   Flux  Range                 biomass_reac...  0.0266
    ---------------  --------  --------------------  ---------------  --------  --------------------
    glc_D_e          0.3       [0.26, 4.5]           h_e              0.536     [-0.715, 10]
    o2_e             0.0575    [0.00272, 10]         lac_L_e          0.536     [-1, 9.62]
    arg_L_e          0.049     [-0.0263, 0.084]      pyr_e            0.0476    [-0.11, 10]
    fe2_e            0.0431    [1, -1]               fe3_e            0.0431    [-1, 1]
    ser_L_e          0.0401    [0.042, -1.22]        ptrc_e           0.0402    [-0.001, 0.181]
    phe_L_e          0.0195    [0.0069, 0.066]       urea_e           0.0394    [0, 0.11]
    lys_L_e          0.0158    [0.0158, 0.0168]      co2_e            0.0358    [-0.001, 6.6]
    leu_L_e          0.0145    [0.0145, 0.0155]      hco3_e           0.0294    [-0.001, 3.96]
    val_L_e          0.0101    [0.00938, 0.094]      glyc_e           0.0138    [-0.001, 1.26]
    glu_L_e          0.01      [-0.001, 0.01]        ac_e             0.00956   [-0.001, 0.0635]
    pro_L_e          0.01      [0.01, -0.0945]       tyr_L_e          0.00832   [0.0549, -0.104]
    asp_L_e          0.00919   [0, 0.01]             gly_e            0.0072    [-0.03, 1.69]
    gln_L_e          0.00867   [-0.277, 0.584]       ura_e            0.00429   [-0.001, 0.00629]
    thr_L_e          0.00832   [0.00832, 0.00832]    gua_e            0.00204   [-0.001, 0.354]
    ile_L_e          0.00761   [0.00761, 0.105]      vacc_e           0.002     [-0.001, 0.505]
    asn_L_e          0.00743   [0.01, -0.851]        ribflv_e         0.002     [0, 0.002]
    met_L_e          0.00407   [0.00207, 0.03]       ade_e            0.00157   [-0.001, 0.32]
    chol_e           0.00305   [0.00305, 0.05]       thym_e           0.00133   [-0.001, 0.00861]
    pi_e             0.00264   [-0.000892, 1]        hdcea_e          0.001     [-0.001, 0.569]
    his_L_e          0.00236   [0.00236, 0.00336]    lnlc_e           0.001     [-0.001, 0.508]
    biomass_other_c  0.00144   [0.00144, 0.00144]    cholate_e        0.001     [0, 0.219]
    dag_hs_e         0.001     [0.001, 0.001]        5mta_e           0.001     [-0.001, 0.00915]
    lpchol_hs_e      0.001     [0.001, 0.001]        ala_B_e          0.001     [0, 0.00829]
    pchol_hs_e       0.001     [0.001, 0.001]        ins_e            0.001     [-0.001, 0.00815]
    pe_hs_e          0.001     [0.001, 0.001]        leuktrD4_e       0.001     [-0.001, 0.002]
    pe_hs_r          0.001     [0.001, 0.001]        c4crn_e          0.001     [0, 0.001]
    ps_hs_e          0.001     [0.001, 0.001]        lac_D_e          0.001     [0, 0.001]
    tag_hs_e         0.001     [0.001, 0.001]        leuktrB4_e       0.001     [0, 0.001]
    4hpro_LT_m       0.001     [0, 0.001]            cys_L_e          0.000761  [0.0287, -0.0626]
    Lcystin_e        0.001     [0, 0.001]            for_e            0.000543  [-0.001, 0.252]
    arab_L_e         0.001     [0, 0.001]            sprm_c           0.0005    [0, 0.00508]
    carn_e           0.001     [0, 0.001]            dgpi_prot_hs_r   0.000333  [0.000333, 0.000333]
    crn_e            0.001     [0, 0.001]            gpi_sig_r        0.000333  [0.000333, 0.000333]
    dcsptn1_e        0.001     [0, 0.001]            dttp_m           0.000326  [0, 0.00994]
    fad_e            0.001     [0, 0.001]            h2o_e            0         [-5.26, 10]
    fmn_e            0.001     [0, 0.001]            h2o2_e           0         [-0.001, 10]
    gtp_e            0.001     [0, 0.001]            gal_e            0         [0, 2.12]
    hdca_e           0.001     [0, 0.001]            abt_e            0         [0, 1.7]
    leuktrA4_e       0.001     [0, 0.001]            oxa_e            0         [0, 1.58]
    ocdca_e          0.001     [0, 0.001]            ala_L_e          0         [-0.01, 1.35]
    strdnc_e         0.001     [0, 0.001]            glyc_S_e         0         [0, 1.2]
    thymd_e          0.001     [0, 0.001]            5oxpro_e         0         [0, 0.936]
    gdp_e            0.001     [0.001, -0.001]       ppi_e            0         [0, 0.5]
    gmp_e            0.001     [0.001, -0.002]       urate_e          0         [0, 0.411]
    leuktrC4_e       0.001     [0.001, -0.002]       nrvnc_e          0         [0, 0.367]
    cytd_e           0.001     [0.001, -0.00629]     chsterol_e       0         [0, 0.225]
    dcmp_e           0.001     [0.001, -0.00629]     C02528_e         0         [0, 0.217]
    dctp_n           0.001     [0.001, -0.00629]     fuc_L_e          0         [0, 0.188]
    dcyt_e           0.001     [0.001, -0.00629]     Rtotal_e         0         [-0.001, 0.172]
    udp_e            0.001     [0.001, -0.00629]     xolest2_hs_e     0         [-0.001, 0.169]
    ump_e            0.001     [0.001, -0.00629]     4hphac_e         0         [0, 0.159]
    utp_e            0.001     [0.001, -0.00629]     tymsf_e          0         [0, 0.159]
    ahcys_e          0.001     [0.001, -0.00815]     fucfucgalacg...  0         [0, 0.152]
    atp_e            0.001     [0.001, -0.00815]     galgalfucfuc...  0         [0, 0.096]
    dttp_n           0.001     [0.001, -0.00894]     fol_c            0         [-0.001, 0.052]
    dad_2_e          0.001     [0.001, -0.0142]      glyb_e           0         [0, 0.047]
    din_e            0.001     [0.001, -0.0142]      gthox_e          0         [0, 0.0456]
    orn_e            0.001     [0.001, -0.105]       citr_L_c         0         [-0.001, 0.0348]
    gchola_e         0.001     [0.001, -0.216]       2hb_e            0         [0, 0.0279]
    ocdcea_e         0.001     [0.001, -0.488]       q10h2_e          0         [0, 0.0261]
    man_e            0.001     [0.001, -1.21]        anth_c           0         [0, 0.0156]
    xylt_e           0.001     [0.001, -1.7]         bildglcur_e      0         [0, 0.0151]
    inost_e          0.000953  [0.000953, -2.83]     bilglcur_e       0         [-0.001, 0.0141]
    mag_hs_e         0.000465  [0, 0.001]            co_e             0         [0, 0.0131]
    sph1p_e          0.000465  [0.001, -0.498]       pheme_e          0         [0, 0.0131]
    pglyc_hs_e       0.000388  [0.000388, 0.000388]  octdececoa_c     0         [-0.000333, 0.0102]
    trp_L_e          0.000354  [0.000354, 0.016]     35cgmp_e         0         [0, 0.00915]
    datp_n           0.000351  [0.000351, 0.000351]  camp_e           0         [0, 0.00915]
    pre_prot_r       0.000333  [0.000333, 0.000333]  spmd_e           0         [-0.001, 0.00915]
    dgtp_n           0.000263  [0.000263, 0.000263]  adn_e            0         [-0.001, 0.00815]
    ha_e             0.000167  [0.001, -0.417]       gsn_e            0         [-0.001, 0.00815]
    dlnlcg_e         0.000128  [0.001, -0.001]       uri_e            0         [-0.001, 0.00629]
    so4_e            0         [-0.0933, 0.164]      Ser_Thr_l        0         [-0.001, 0.005]
    fol_e            0         [-0.003, 0.05]        adrn_e           0         [-0.001, 0.003]
    pnto_R_e         0         [0, 0.00915]          1mncam_e         0         [0, 0.002]
    3hpvs_e          0         [0, 0.001]            5adtststeron...  0         [0, 0.002]
    4nph_e           0         [0, 0.001]            andrstrnglc_e    0         [0, 0.002]
    5fthf_e          0         [0, 0.001]            estroneglc_e     0         [0, 0.002]
    5mthf_e          0         [0, 0.001]            imp_e            0         [0, 0.002]
    C02470_e         0         [0, 0.001]            n5m2masn_g       0         [0, 0.002]
    HC00822_e        0         [0, 0.001]            nac_e            0         [0, 0.002]
    acac_e           0         [0, 0.001]            arachd_e         0         [-0.001, 0.002]
    acmp_e           0         [0, 0.001]            core8_g          0         [0, 0.0012]
    arach_e          0         [0, 0.001]            3bcrn_e          0         [0, 0.001]
    bilirub_e        0         [0, 0.001]            3hdececrn_e      0         [0, 0.001]
    crvnc_e          0         [0, 0.001]            3hpvstet_e       0         [0, 0.001]
    datp_m           0         [0, 0.001]            3ivcrn_e         0         [0, 0.001]
    dca_e            0         [0, 0.001]            3mlda_e          0         [0, 0.001]
    dheas_e          0         [0, 0.001]            3octdec2crn_e    0         [0, 0.001]
    doco13ac_e       0         [0, 0.001]            3octdeccrn_e     0         [0, 0.001]
    dopa_e           0         [0, 0.001]            3octdece1crn_e   0         [0, 0.001]
    etoh_e           0         [0, 0.001]            3tetd7ecoacrn_e  0         [0, 0.001]
    gam_e            0         [0, 0.001]            3thexddcoacrn_e  0         [0, 0.001]
    gluala_e         0         [0, 0.001]            3ttetddcoacrn_e  0         [0, 0.001]
    hista_e          0         [0, 0.001]            4nphsf_e         0         [0, 0.001]
    hpdca_e          0         [0, 0.001]            aprgstrn_e       0         [0, 0.001]
    lnlnca_e         0         [0, 0.001]            c5dc_e           0         [0, 0.001]
    lnlncg_e         0         [0, 0.001]            clpnd_e          0         [0, 0.001]
    n2m2nmasn_e      0         [0, 0.001]            dmhptcrn_e       0         [0, 0.001]
    o2s_e            0         [0, 0.001]            dopasf_e         0         [0, 0.001]
    octa_e           0         [0, 0.001]            epoxtac_e        0         [0, 0.001]
    orot_e           0         [0, 0.001]            ivcrn_e          0         [0, 0.001]
    phyt_e           0         [0, 0.001]            prostgf2_e       0         [0, 0.001]
    prgstrn_e        0         [0, 0.001]            retn_e           0         [0, 0.001]
    prostge2_e       0         [0, 0.001]            sulpacmp_e       0         [0, 0.001]
    ptdca_e          0         [0, 0.001]            triodthysuf_e    0         [0, 0.001]
    retinol_e        0         [0, 0.001]            tststeroneglc_e  0         [0, 0.001]
    tacr_e           0         [0, 0.001]            5adtststerone_e  0         [-0.001, 0.001]
    tmndnc_e         0         [0, 0.001]            Asn_X_Ser_Thr_l  0         [-0.001, 0.001]
    triodthy_e       0         [0, 0.001]            estradiol_e      0         [-0.001, 0.001]
    tststerone_e     0         [0, 0.001]            estrones_e       0         [-0.001, 0.001]
    ttdca_e          0         [0, 0.001]            fald_e           0         [-0.001, 0.001]
    whhdca_e         0         [0, 0.001]            idp_e            0         [-0.001, 0.001]
    T_antigen_g      0         [-0.0002, 0.001]      itp_e            0         [-0.001, 0.001]
    andrstrn_e       0         [-0.001, 0.001]       meoh_e           0         [-0.001, 0.001]
    ncam_c           0         [-0.001, 0.001]       taur_c           0         [-0.001, 0.001]
    c6crn_e          0         [0, -0.001]           tchola_e         0         [-0.001, 0.001]
                                                     pmtcoa_r         0         [0.000333, -0.001]



```python
model_KeratinocytesMean4.summary(fva=True)
```

    IN FLUXES                                        OUT FLUXES                                     OBJECTIVES
    -----------------------------------------------  ---------------------------------------------  -----------------------
    id                   Flux  Range                 id                   Flux  Range               biomass_reac...  0.0242
    ---------------  --------  --------------------  ---------------  --------  ------------------
    glc_D_e          0.267     [0.203, 4.5]          h_e              0.521     [-0.619, 10]
    o2_e             0.0483    [0.00246, 8.68]       lac_L_e          0.471     [-1, 9.73]
    ser_L_e          0.0246    [0.042, -1.12]        pyr_e            0.0488    [-0.11, 9.86]
    lys_L_e          0.0144    [0.0144, 0.0154]      h2o_e            0.0388    [-1.54, 10]
    leu_L_e          0.0132    [0.0132, 0.0142]      co2_e            0.0385    [-0.001, 10]
    val_L_e          0.0115    [0.00855, 0.094]      glyc_e           0.0133    [-0.001, 1.6]
    gly_e            0.0111    [0.03, -1.14]         ala_B_e          0.00683   [-0.001, 1.31]
    pro_L_e          0.01      [0.01, 0.01]          hdcea_e          0.00633   [-0.001, 0.522]
    arg_L_e          0.00871   [0.00671, 0.00871]    pi_e             0.00563   [0.0162, -0.896]
    glu_L_e          0.00815   [0.01, -1.08]         gua_e            0.00535   [-0.001, 0.342]
    gln_L_e          0.0079    [-0.251, 0.584]       dgsn_e           0.00253   [-0.001, 0.0188]
    thr_L_e          0.00758   [0.00758, 0.00758]    ade_e            0.00184   [-0.001, 0.347]
    ile_L_e          0.00694   [0.00694, 0.105]      prostgh2_e       0.001     [0, 0.002]
    asn_L_e          0.00677   [0.00677, -0.828]     3ivcrn_e         0.001     [0, 0.001]
    phe_L_e          0.00629   [0.00629, 0.066]      epoxtac_e        0.001     [0, 0.001]
    asp_L_e          0.00555   [0, 0.01]             lac_D_e          0.001     [0, 0.001]
    met_L_e          0.00472   [0.00371, 0.03]       ndersv_e         0.001     [0, 0.001]
    tyr_L_e          0.00387   [-0.0558, 0.104]      ahcys_e          0.000862  [-0.001, 0.0136]
    his_L_e          0.00206   [0.00206, 0.00307]    orn_e            0.000799  [0, 0.0865]
    biomass_other_c  0.00131   [0.00131, 0.00131]    dad_2_e          0.00068   [-0.001, 0.0188]
    adrn_e           0.001     [0, 0.001]            for_e            0.000495  [-0.001, 0.395]
    amp_e            0.001     [0, 0.001]            2hb_e            0.000149  [0, 0.0273]
    carn_e           0.001     [0, 0.001]            oxa_e            0         [0, 2.16]
    crn_e            0.001     [0, 0.001]            gal_e            0         [0, 2.15]
    gtp_e            0.001     [0, 0.001]            glyc_S_e         0         [0, 1.6]
    lgnc_e           0.001     [0, 0.001]            acmana_e         0         [0, 1.32]
    o2s_e            0.001     [0, 0.001]            ala_L_e          0         [-0.01, 1.16]
    pchol_hs_e       0.001     [0, 0.001]            taur_c           0         [-0.001, 1.09]
    prostge2_e       0.001     [0, 0.001]            taur_e           0         [-0.001, 1.09]
    rsv_e            0.001     [0, 0.001]            5oxpro_e         0         [0, 1.09]
    tacr_e           0.001     [0, 0.001]            fe2_e            0         [-1, 1]
    tag_hs_e         0.001     [0, 0.001]            fe3_e            0         [-1, 1]
    dag_hs_e         0.001     [-0.000924, 0.001]    ethamp_r         0         [0, 0.852]
    glyc3p_e         0.001     [-0.000924, 0.001]    ha_pre1_e        0         [-0.001, 0.848]
    lpchol_hs_e      0.001     [-0.000924, 0.001]    hdca_e           0         [-0.001, 0.511]
    pe_hs_e          0.001     [-0.000924, 0.001]    vacc_e           0         [-0.001, 0.461]
    pe_hs_r          0.001     [-0.000924, 0.001]    urate_e          0         [0, 0.429]
    ps_hs_e          0.001     [-0.000924, 0.001]    ha_e             0         [-0.001, 0.424]
    dcsptn1_e        0.001     [-0.001, 0.001]       chsterol_e       0         [0, 0.391]
    gdp_e            0.001     [0.001, -0.001]       cholate_e        0         [0, 0.378]
    cmp_e            0.001     [0.001, -0.002]       sphs1p_e         0         [-0.001, 0.367]
    gmp_e            0.001     [0.001, -0.002]       tdchola_e        0         [0, 0.354]
    cytd_e           0.001     [0.001, -0.00753]     eicostet_e       0         [-0.001, 0.352]
    dctp_n           0.001     [0.001, -0.00753]     nrvnc_e          0         [0, 0.329]
    udp_e            0.001     [0.001, -0.00753]     Rtotal_e         0         [-0.001, 0.247]
    ump_e            0.001     [0.001, -0.00753]     xolest2_hs_e     0         [-0.001, 0.245]
    ura_e            0.001     [0.001, -0.00753]     digalsgalsid...  0         [0, 0.199]
    uri_e            0.001     [0.001, -0.00753]     4hphac_e         0         [0, 0.16]
    utp_e            0.001     [0.001, -0.00753]     acac_e           0         [-0.001, 0.16]
    adn_e            0.001     [0.001, -0.0136]      glyb_e           0         [0, 0.128]
    aicar_e          0.001     [0.001, -0.0136]      spc_hs_e         0         [0, 0.102]
    atp_e            0.001     [0.001, -0.0136]      q10h2_e          0         [0, 0.0904]
    gsn_e            0.001     [0.001, -0.0136]      cgly_e           0         [-0.001, 0.0898]
    ins_e            0.001     [0.001, -0.0136]      gthrd_e          0         [-0.001, 0.0898]
    xmp_e            0.001     [0.001, -0.0136]      q10_e            0         [-0.001, 0.089]
    datp_n           0.001     [0.001, -0.0143]      pheacgln_e       0         [0, 0.0597]
    dgtp_n           0.001     [0.001, -0.0188]      35cgmp_e         0         [0, 0.0146]
    din_e            0.001     [0.001, -0.0188]      camp_e           0         [0, 0.0146]
    strdnc_e         0.001     [0.001, -0.385]       bilglcur_e       0         [0, 0.0108]
    hxan_e           0.001     [0.001, -0.428]       co_e             0         [0, 0.0108]
    HC00342_e        0.001     [0.001, -1.09]        pheme_e          0         [0, 0.0108]
    cit_e            0.001     [0.001, -1.09]        dttp_m           0         [0, 0.00921]
    xylt_e           0.001     [0.001, -1.72]        Ser_Thr_l        0         [-0.001, 0.005]
    cys_L_e          0.00098   [-0.0262, 0.0626]     fuc14galacgl...  0         [0, 0.004]
    elaid_e          0.000903  [0.001, -0.461]       fuc_L_e          0         [0, 0.004]
    inost_e          0.000565  [0.000565, 0.00156]   galfuc12gal1...  0         [0, 0.004]
    sph1p_e          0.000424  [0.001, -0.367]       galfucgalacg...  0         [0, 0.004]
    pglyc_hs_e       0.000353  [0.000353, 0.000353]  fald_e           0         [-0.001, 0.003]
    trp_L_e          0.000323  [0.000323, 0.000323]  man_e            0         [-0.001, 0.003]
    dttp_n           0.000317  [0.001, -0.00821]     meoh_e           0         [-0.001, 0.003]
    orot_e           0.000296  [0, 0.001]            34dhoxpeg_e      0         [0, 0.002]
    so4_e            0         [-0.0898, 1]          5adtststeron...  0         [0, 0.002]
    3hpvs_e          0         [0, 0.001]            andrstrnglc_e    0         [0, 0.002]
    4hpro_LT_m       0         [0, 0.001]            imp_e            0         [0, 0.002]
    5HPET_c          0         [0, 0.001]            leuktrB4_e       0         [0, 0.002]
    5fthf_e          0         [0, 0.001]            leuktrC4_e       0         [0, 0.002]
    HC00822_e        0         [0, 0.001]            cdp_e            0         [-0.001, 0.002]
    acmp_e           0         [0, 0.001]            ctp_e            0         [-0.001, 0.002]
    allop_e          0         [0, 0.001]            octdececoa_c     0         [-0.001, 0.002]
    arab_L_e         0         [0, 0.001]            pmtcoa_r         0         [-0.001, 0.002]
    arach_e          0         [0, 0.001]            fucfucfucgal...  0         [0, 0.00133]
    c81coa_c         0         [0, 0.001]            n5m2masn_g       0         [0, 0.00133]
    crvnc_e          0         [0, 0.001]            core5_g          0         [0, 0.0012]
    dca_e            0         [0, 0.001]            core7_g          0         [0, 0.0012]
    dd2coa_c         0         [0, 0.001]            core8_g          0         [0, 0.0012]
    decdicoa_c       0         [0, 0.001]            3bcrn_e          0         [0, 0.001]
    dheas_e          0         [0, 0.001]            3ddcrn_e         0         [0, 0.001]
    doco13ac_e       0         [0, 0.001]            3deccrn_e        0         [0, 0.001]
    dopa_e           0         [0, 0.001]            3hdececrn_e      0         [0, 0.001]
    etoh_e           0         [0, 0.001]            3hexdcrn_e       0         [0, 0.001]
    fmn_e            0         [0, 0.001]            3hpvstet_e       0         [0, 0.001]
    gam_e            0         [0, 0.001]            3mlda_e          0         [0, 0.001]
    gluala_e         0         [0, 0.001]            3octdec2crn_e    0         [0, 0.001]
    h2o2_e           0         [0, 0.001]            3octdeccrn_e     0         [0, 0.001]
    hista_e          0         [0, 0.001]            3octdece1crn_e   0         [0, 0.001]
    hpdca_e          0         [0, 0.001]            3tdcrn_e         0         [0, 0.001]
    lnlc_e           0         [0, 0.001]            3tetd7ecoacrn_e  0         [0, 0.001]
    n2m2nmasn_e      0         [0, 0.001]            3thexddcoacrn_e  0         [0, 0.001]
    ocdca_e          0         [0, 0.001]            3ttetddcoacrn_e  0         [0, 0.001]
    octa_e           0         [0, 0.001]            5mta_e           0         [0, 0.001]
    ppa_e            0         [0, 0.001]            abt_e            0         [0, 0.001]
    pre_prot_r       0         [0, 0.001]            aprgstrn_e       0         [0, 0.001]
    prgstrn_e        0         [0, 0.001]            c4crn_e          0         [0, 0.001]
    ptvst_e          0         [0, 0.001]            c51crn_e         0         [0, 0.001]
    retinol_e        0         [0, 0.001]            c5dc_e           0         [0, 0.001]
    spmd_e           0         [0, 0.001]            c6crn_e          0         [0, 0.001]
    tetdece1coa_c    0         [0, 0.001]            c81crn_e         0         [0, 0.001]
    tmndnc_e         0         [0, 0.001]            c8crn_e          0         [0, 0.001]
    tststerone_e     0         [0, 0.001]            clpnd_e          0         [0, 0.001]
    ttdca_e          0         [0, 0.001]            ddece1crn_e      0         [0, 0.001]
    whhdca_e         0         [0, 0.001]            decdicrn_e       0         [0, 0.001]
    T_antigen_g      0         [-0.0002, 0.001]      dem2emgacpai...  0         [0, 0.001]
    andrstrn_e       0         [-0.001, 0.001]       dgpi_prot_hs_r   0         [0, 0.001]
    nrpphr_e         0         [-0.001, 0.001]       fol_e            0         [0, 0.001]
                                                     gpi_sig_r        0         [0, 0.001]
                                                     gthox_e          0         [0, 0.001]
                                                     ivcrn_e          0         [0, 0.001]
                                                     mem2emgacpai...  0         [0, 0.001]
                                                     oxyp_e           0         [0, 0.001]
                                                     prostgi2_e       0         [0, 0.001]
                                                     ptvstm3_e        0         [0, 0.001]
                                                     retn_e           0         [0, 0.001]
                                                     ribflv_e         0         [0, 0.001]
                                                     rsvlac_e         0         [0, 0.001]
                                                     sprm_c           0         [0, 0.001]
                                                     sulpacmp_e       0         [0, 0.001]
                                                     tetdece1crn_e    0         [0, 0.001]
                                                     tststeroneglc_e  0         [0, 0.001]
                                                     ttdcrn_e         0         [0, 0.001]
                                                     Asn_X_Ser_Thr_l  0         [-0.000333, 0.001]
                                                     5adtststerone_e  0         [-0.001, 0.001]
                                                     citr_L_c         0         [-0.001, 0.001]
                                                     citr_L_e         0         [-0.001, 0.001]
                                                     idp_e            0         [-0.001, 0.001]
                                                     itp_e            0         [-0.001, 0.001]
                                                     leuktrA4_e       0         [-0.001, 0.001]
                                                     ptdca_e          0         [0, -0.001]



```python

model_HeLaMean4.summary(fva=True)
```

    IN FLUXES                                        OUT FLUXES                                       OBJECTIVES
    -----------------------------------------------  -----------------------------------------------  -----------------------
    id                   Flux  Range                 id                   Flux  Range                 biomass_reac...  0.0266
    ---------------  --------  --------------------  ---------------  --------  --------------------
    glc_D_e          0.289     [0.255, 4.5]          h_e              0.552     [-0.715, 10]
    o2_e             0.131     [0.00272, 10]         lac_L_e          0.483     [-1, 9.62]
    fol_e            0.05      [-0.003, 0.05]        h2o_e            0.101     [-5.04, 10]
    arg_L_e          0.0429    [-0.027, 0.084]       co2_e            0.0903    [-0.001, 6.61]
    val_L_e          0.0245    [0.00938, 0.094]      pyr_e            0.0865    [-0.11, 10]
    ser_L_e          0.0227    [0.042, -1.22]        fol_c            0.05      [-0.001, 0.052]
    lys_L_e          0.0158    [0.0158, 0.0168]      ptrc_e           0.0333    [-0.001, 0.183]
    leu_L_e          0.0145    [0.0145, 0.0155]      urea_e           0.0333    [0, 0.111]
    pro_L_e          0.01      [0.01, -0.0965]       ac_e             0.00879   [-0.001, 0.0669]
    gln_L_e          0.00907   [-0.279, 0.584]       ura_e            0.00553   [-0.001, 0.00729]
    glu_L_e          0.00883   [-0.001, 0.01]        urate_e          0.005     [0, 0.415]
    thr_L_e          0.00832   [0.00832, 0.00832]    hxan_e           0.00257   [-0.001, 0.412]
    ile_L_e          0.00761   [0.00761, 0.105]      ade_e            0.00165   [-0.001, 0.321]
    asn_L_e          0.00743   [0.01, -0.853]        thym_e           0.00133   [-0.001, 0.00962]
    phe_L_e          0.0069    [0.0069, 0.066]       lnlc_e           0.00111   [-0.001, 0.508]
    tyr_L_e          0.00425   [-0.0549, 0.104]      glyc_e           0.001     [-0.001, 1.26]
    met_L_e          0.00407   [0.00207, 0.03]       cholate_e        0.001     [0, 0.219]
    chol_e           0.00305   [0.00305, 0.05]       5mta_e           0.001     [-0.001, 0.0102]
    his_L_e          0.00236   [0.00236, 0.00336]    spmd_e           0.001     [-0.001, 0.0102]
    biomass_other_c  0.00144   [0.00144, 0.00144]    ala_B_e          0.001     [0, 0.00929]
    cys_L_e          0.00124   [-0.0287, 0.0626]     ribflv_e         0.001     [0, 0.002]
    dag_hs_e         0.001     [0.001, 0.001]        leuktrD4_e       0.001     [-0.001, 0.002]
    lpchol_hs_e      0.001     [0.001, 0.001]        c6crn_e          0.001     [0, 0.001]
    pchol_hs_e       0.001     [0.001, 0.001]        lac_D_e          0.001     [0, 0.001]
    pe_hs_e          0.001     [0.001, 0.001]        5adtststerone_e  0.001     [-0.001, 0.001]
    pe_hs_r          0.001     [0.001, 0.001]        gua_e            0.000776  [-0.001, 0.357]
    ps_hs_e          0.001     [0.001, 0.001]        for_e            0.000543  [-0.001, 0.253]
    tag_hs_e         0.001     [0.001, 0.001]        hdcea_e          0.000535  [-0.001, 0.57]
    carn_e           0.001     [0, 0.001]            dgpi_prot_hs_r   0.000333  [0.000333, 0.000333]
    crn_e            0.001     [0, 0.001]            gpi_sig_r        0.000333  [0.000333, 0.000333]
    dcsptn1_e        0.001     [0, 0.001]            dttp_m           0.000326  [0, 0.0109]
    fad_e            0.001     [0, 0.001]            Asn_X_Ser_Thr_l  6.7e-05   [-0.001, 0.001]
    gtp_e            0.001     [0, 0.001]            h2o2_e           0         [-0.001, 10]
    mag_hs_e         0.001     [0, 0.001]            gal_e            0         [0, 2.12]
    o2s_e            0.001     [0, 0.001]            abt_e            0         [0, 1.7]
    orot_e           0.001     [0, 0.001]            xylt_e           0         [-0.001, 1.7]
    thymd_e          0.001     [0, 0.001]            oxa_e            0         [0, 1.58]
    tststerone_e     0.001     [0, 0.001]            ala_L_e          0         [-0.01, 1.35]
    gdp_e            0.001     [0.001, -0.001]       glyc_S_e         0         [0, 1.2]
    gmp_e            0.001     [0.001, -0.002]       fe2_e            0         [-1, 1]
    leuktrC4_e       0.001     [0.001, -0.002]       fe3_e            0         [-1, 1]
    cytd_e           0.001     [0.001, -0.00729]     5oxpro_e         0         [0, 0.938]
    dcmp_e           0.001     [0.001, -0.00729]     vacc_e           0         [-0.001, 0.505]
    dctp_n           0.001     [0.001, -0.00729]     ppi_e            0         [0, 0.503]
    dcyt_e           0.001     [0.001, -0.00729]     ha_e             0         [-0.001, 0.418]
    duri_e           0.001     [0.001, -0.00729]     nrvnc_e          0         [0, 0.367]
    ump_e            0.001     [0.001, -0.00729]     chsterol_e       0         [0, 0.225]
    utp_e            0.001     [0.001, -0.00729]     C02528_e         0         [0, 0.217]
    adn_e            0.001     [0.001, -0.00915]     fuc_L_e          0         [0, 0.189]
    ahcys_e          0.001     [0.001, -0.00915]     Rtotal_e         0         [-0.001, 0.172]
    atp_e            0.001     [0.001, -0.00915]     xolest2_hs_e     0         [-0.001, 0.17]
    gsn_e            0.001     [0.001, -0.00915]     4hphac_e         0         [0, 0.159]
    ins_e            0.001     [0.001, -0.00915]     tymsf_e          0         [0, 0.159]
    xmp_e            0.001     [0.001, -0.00915]     fucfucgalacg...  0         [0, 0.152]
    dttp_n           0.001     [0.001, -0.00994]     galgalfucfuc...  0         [0, 0.0961]
    dad_2_e          0.001     [0.001, -0.0176]      glyb_e           0         [0, 0.047]
    datp_n           0.001     [0.001, -0.0176]      gthox_e          0         [0, 0.0456]
    dgsn_e           0.001     [0.001, -0.0176]      citr_L_c         0         [-0.001, 0.0355]
    dgtp_n           0.001     [0.001, -0.0176]      2hb_e            0         [0, 0.0279]
    orn_e            0.001     [0.001, -0.107]       q10h2_e          0         [0, 0.0266]
    gchola_e         0.001     [0.001, -0.216]       anth_c           0         [0, 0.0156]
    ocdcea_e         0.001     [0.001, -0.489]       bildglcur_e      0         [0, 0.0153]
    inost_e          0.000953  [0.000953, -2.83]     bilglcur_e       0         [-0.001, 0.0143]
    man_e            0.0008    [0.001, -1.21]        co_e             0         [0, 0.0133]
    sph1p_e          0.000465  [0.001, -0.498]       pheme_e          0         [0, 0.0133]
    whhdca_e         0.000461  [0, 0.001]            octdececoa_c     0         [-0.000333, 0.0112]
    pglyc_hs_e       0.000388  [0.000388, 0.000388]  35cgmp_e         0         [0, 0.0102]
    trp_L_e          0.000354  [0.000354, 0.016]     camp_e           0         [0, 0.0102]
    pre_prot_r       0.000333  [0.000333, 0.000333]  uri_e            0         [-0.001, 0.00729]
    gly_e            0.000248  [0.03, -1.7]          sprm_c           0         [0, 0.00558]
    udp_e            0.000242  [0.001, -0.00729]     Ser_Thr_l        0         [-0.001, 0.005]
    strdnc_e         0.000114  [0, 0.001]            adrn_e           0         [-0.001, 0.003]
    n2m2nmasn_e      6.7e-05   [0, 0.001]            1mncam_e         0         [0, 0.002]
    pi_e             0         [-0.00605, 1]         5adtststeron...  0         [0, 0.002]
    so4_e            0         [-0.0933, 0.164]      andrstrnglc_e    0         [0, 0.002]
    pnto_R_e         0         [0, 0.0102]           estroneglc_e     0         [0, 0.002]
    asp_L_e          0         [0, 0.01]             imp_e            0         [0, 0.002]
    3hpvs_e          0         [0, 0.001]            n5m2masn_g       0         [0, 0.002]
    4hpro_LT_m       0         [0, 0.001]            nac_e            0         [0, 0.002]
    4nph_e           0         [0, 0.001]            arachd_e         0         [-0.001, 0.002]
    5fthf_e          0         [0, 0.001]            core8_g          0         [0, 0.0012]
    5mthf_e          0         [0, 0.001]            3bcrn_e          0         [0, 0.001]
    C02470_e         0         [0, 0.001]            3hdececrn_e      0         [0, 0.001]
    HC00822_e        0         [0, 0.001]            3hpvstet_e       0         [0, 0.001]
    Lcystin_e        0         [0, 0.001]            3ivcrn_e         0         [0, 0.001]
    acac_e           0         [0, 0.001]            3mlda_e          0         [0, 0.001]
    acmp_e           0         [0, 0.001]            3octdec2crn_e    0         [0, 0.001]
    arab_L_e         0         [0, 0.001]            3octdeccrn_e     0         [0, 0.001]
    arach_e          0         [0, 0.001]            3octdece1crn_e   0         [0, 0.001]
    bilirub_e        0         [0, 0.001]            3tetd7ecoacrn_e  0         [0, 0.001]
    crvnc_e          0         [0, 0.001]            3thexddcoacrn_e  0         [0, 0.001]
    datp_m           0         [0, 0.001]            3ttetddcoacrn_e  0         [0, 0.001]
    dca_e            0         [0, 0.001]            4nphsf_e         0         [0, 0.001]
    dheas_e          0         [0, 0.001]            aprgstrn_e       0         [0, 0.001]
    doco13ac_e       0         [0, 0.001]            c4crn_e          0         [0, 0.001]
    dopa_e           0         [0, 0.001]            c5dc_e           0         [0, 0.001]
    etoh_e           0         [0, 0.001]            clpnd_e          0         [0, 0.001]
    fmn_e            0         [0, 0.001]            dmhptcrn_e       0         [0, 0.001]
    gam_e            0         [0, 0.001]            dopasf_e         0         [0, 0.001]
    gluala_e         0         [0, 0.001]            epoxtac_e        0         [0, 0.001]
    hdca_e           0         [0, 0.001]            ivcrn_e          0         [0, 0.001]
    hista_e          0         [0, 0.001]            leuktrB4_e       0         [0, 0.001]
    hpdca_e          0         [0, 0.001]            prostgf2_e       0         [0, 0.001]
    leuktrA4_e       0         [0, 0.001]            retn_e           0         [0, 0.001]
    lnlnca_e         0         [0, 0.001]            sulpacmp_e       0         [0, 0.001]
    lnlncg_e         0         [0, 0.001]            triodthysuf_e    0         [0, 0.001]
    ocdca_e          0         [0, 0.001]            tststeroneglc_e  0         [0, 0.001]
    octa_e           0         [0, 0.001]            dlnlcg_e         0         [-0.001, 0.001]
    phyt_e           0         [0, 0.001]            eicostet_e       0         [-0.001, 0.001]
    prgstrn_e        0         [0, 0.001]            estradiol_e      0         [-0.001, 0.001]
    prostge2_e       0         [0, 0.001]            fald_e           0         [-0.001, 0.001]
    ptdca_e          0         [0, 0.001]            idp_e            0         [-0.001, 0.001]
    retinol_e        0         [0, 0.001]            itp_e            0         [-0.001, 0.001]
    tacr_e           0         [0, 0.001]            taur_c           0         [-0.001, 0.001]
    tmndnc_e         0         [0, 0.001]            tchola_e         0         [-0.001, 0.001]
    triodthy_e       0         [0, 0.001]
    ttdca_e          0         [0, 0.001]
    T_antigen_g      0         [-0.0002, 0.001]
    pmtcoa_r         0         [-0.000333, 0.001]
    andrstrn_e       0         [-0.001, 0.001]
    estrones_e       0         [-0.001, 0.001]
    meoh_e           0         [-0.001, 0.001]
    ncam_c           0         [-0.001, 0.001]



```python
model_KeratinocytesMean5.summary(fva=True)
```

    IN FLUXES                                        OUT FLUXES                                    OBJECTIVES
    -----------------------------------------------  --------------------------------------------  -----------------------
    id                   Flux  Range                 id                   Flux  Range              biomass_reac...  0.0242
    ---------------  --------  --------------------  ---------------  --------  -----------------
    glc_D_e          0.276     [0.203, 4.5]          h_e              0.534     [-0.62, 10]
    o2_e             0.0259    [0.00246, 8.68]       lac_L_e          0.51      [-1, 9.73]
    ser_L_e          0.0244    [0.042, -1.12]        pyr_e            0.0245    [-0.11, 9.86]
    phe_L_e          0.0158    [0.00629, 0.066]      h2o_e            0.0236    [-1.54, 10]
    lys_L_e          0.0144    [0.0144, 0.0154]      co2_e            0.0174    [-0.001, 10]
    leu_L_e          0.0132    [0.0132, 0.0142]      glyc_e           0.0123    [-0.001, 1.6]
    gly_e            0.0111    [0.03, -1.14]         tyr_L_e          0.00566   [0.0558, -0.104]
    pro_L_e          0.01      [0.01, 0.01]          ura_e            0.00382   [-0.001, 0.00853]
    glu_L_e          0.00982   [0.01, -1.08]         urate_e          0.00202   [0, 0.429]
    arg_L_e          0.00871   [0.00671, 0.00871]    Rtotal_e         0.002     [-0.001, 0.247]
    val_L_e          0.00855   [0.00855, 0.094]      dad_2_e          0.00145   [-0.001, 0.0198]
    asp_L_e          0.00855   [0, 0.01]             2hb_e            0.00101   [0, 0.0273]
    thr_L_e          0.00758   [0.00758, 0.00758]    ala_B_e          0.001     [-0.001, 1.31]
    gln_L_e          0.00744   [-0.252, 0.584]       gthrd_e          0.001     [-0.001, 0.0898]
    ile_L_e          0.00694   [0.00694, 0.105]      leuktrB4_e       0.001     [0, 0.002]
    asn_L_e          0.00677   [0.00677, -0.829]     epoxtac_e        0.001     [0, 0.001]
    met_L_e          0.00472   [0.00371, 0.03]       lac_D_e          0.001     [0, 0.001]
    his_L_e          0.00206   [0.00206, 0.00307]    ndersv_e         0.001     [0, 0.001]
    pi_e             0.00137   [-0.0162, 0.899]      oxyp_e           0.001     [0, 0.001]
    biomass_other_c  0.00131   [0.00131, 0.00131]    dgsn_e           0.00076   [-0.001, 0.0198]
    4hpro_LT_m       0.001     [0, 0.001]            ethamp_r         0.000576  [0, 0.852]
    allop_e          0.001     [0, 0.001]            hdca_e           0.000576  [-0.001, 0.511]
    amp_e            0.001     [0, 0.001]            so4_e            0.0005    [0.0898, -1]
    carn_e           0.001     [0, 0.001]            for_e            0.000495  [-0.001, 0.395]
    o2s_e            0.001     [0, 0.001]            oxa_e            0         [0, 2.16]
    pchol_hs_e       0.001     [0, 0.001]            gal_e            0         [0, 2.15]
    rsv_e            0.001     [0, 0.001]            xylt_e           0         [-0.001, 1.72]
    tacr_e           0.001     [0, 0.001]            glyc_S_e         0         [0, 1.6]
    tag_hs_e         0.001     [0, 0.001]            acmana_e         0         [0, 1.32]
    dag_hs_e         0.001     [-0.000924, 0.001]    ala_L_e          0         [-0.01, 1.16]
    glyc3p_e         0.001     [-0.000924, 0.001]    taur_c           0         [-0.001, 1.09]
    lpchol_hs_e      0.001     [-0.000924, 0.001]    taur_e           0         [-0.001, 1.09]
    pe_hs_e          0.001     [-0.000924, 0.001]    5oxpro_e         0         [0, 1.09]
    pe_hs_r          0.001     [-0.000924, 0.001]    HC00342_e        0         [-0.001, 1.09]
    ps_hs_e          0.001     [-0.000924, 0.001]    fe2_e            0         [-1, 1]
    leuktrA4_e       0.001     [0.001, -0.001]       ha_pre1_e        0         [-0.001, 0.849]
    cytd_e           0.001     [0.001, -0.00853]     hdcea_e          0         [-0.001, 0.522]
    dctp_n           0.001     [0.001, -0.00853]     vacc_e           0         [-0.001, 0.461]
    udp_e            0.001     [0.001, -0.00853]     elaid_e          0         [-0.001, 0.461]
    ump_e            0.001     [0.001, -0.00853]     hxan_e           0         [-0.001, 0.428]
    uri_e            0.001     [0.001, -0.00853]     ha_e             0         [-0.001, 0.424]
    utp_e            0.001     [0.001, -0.00853]     chsterol_e       0         [0, 0.391]
    aicar_e          0.001     [0.001, -0.0136]      sph1p_e          0         [-0.001, 0.367]
    atp_e            0.001     [0.001, -0.0136]      tdchola_e        0         [0, 0.354]
    xmp_e            0.001     [0.001, -0.0136]      eicostet_e       0         [-0.001, 0.352]
    datp_n           0.001     [0.001, -0.0143]      ade_e            0         [-0.001, 0.348]
    dgtp_n           0.001     [0.001, -0.0198]      gua_e            0         [-0.001, 0.342]
    cgly_e           0.001     [0.001, -0.0898]      nrvnc_e          0         [0, 0.329]
    sphs1p_e         0.001     [0.001, -0.367]       xolest2_hs_e     0         [-0.001, 0.245]
    cit_e            0.001     [0.001, -1.09]        digalsgalsid...  0         [0, 0.199]
    ins_e            0.000899  [0.001, -0.0136]      4hphac_e         0         [0, 0.16]
    cys_L_e          0.000618  [-0.0262, 0.0626]     acac_e           0         [-0.001, 0.16]
    adrn_e           0.0006    [0, 0.001]            glyb_e           0         [0, 0.128]
    lgnc_e           0.000577  [0, 0.001]            spc_hs_e         0         [0, 0.102]
    inost_e          0.000565  [0.000565, 0.00156]   q10h2_e          0         [0, 0.0904]
    pglyc_hs_e       0.000353  [0.000353, 0.000353]  q10_e            0         [-0.001, 0.089]
    trp_L_e          0.000323  [0.000323, 0.000323]  orn_e            0         [0, 0.0865]
    dttp_n           0.000317  [0.001, -0.00921]     pheacgln_e       0         [0, 0.0597]
    orot_e           0.000296  [0, 0.001]            din_e            0         [-0.001, 0.0198]
    adn_e            7.3e-05   [0.001, -0.0136]      35cgmp_e         0         [0, 0.0146]
    fe3_e            0         [-1, 1]               camp_e           0         [0, 0.0146]
    3hpvs_e          0         [0, 0.001]            ahcys_e          0         [-0.001, 0.0136]
    5HPET_c          0         [0, 0.001]            gsn_e            0         [-0.001, 0.0136]
    5fthf_e          0         [0, 0.001]            bilglcur_e       0         [0, 0.0108]
    HC00822_e        0         [0, 0.001]            co_e             0         [0, 0.0108]
    acmp_e           0         [0, 0.001]            pheme_e          0         [0, 0.0108]
    arab_L_e         0         [0, 0.001]            dttp_m           0         [0, 0.0102]
    arach_e          0         [0, 0.001]            duri_e           0         [-0.001, 0.00853]
    c81coa_c         0         [0, 0.001]            Ser_Thr_l        0         [-0.001, 0.005]
    crn_e            0         [0, 0.001]            fuc14galacgl...  0         [0, 0.004]
    crvnc_e          0         [0, 0.001]            galfuc12gal1...  0         [0, 0.004]
    dca_e            0         [0, 0.001]            galfucgalacg...  0         [0, 0.004]
    dd2coa_c         0         [0, 0.001]            man_e            0         [-0.001, 0.003]
    decdicoa_c       0         [0, 0.001]            meoh_e           0         [-0.001, 0.003]
    dheas_e          0         [0, 0.001]            34dhoxpeg_e      0         [0, 0.002]
    doco13ac_e       0         [0, 0.001]            5adtststeron...  0         [0, 0.002]
    dopa_e           0         [0, 0.001]            andrstrnglc_e    0         [0, 0.002]
    etoh_e           0         [0, 0.001]            imp_e            0         [0, 0.002]
    fmn_e            0         [0, 0.001]            leuktrC4_e       0         [0, 0.002]
    gam_e            0         [0, 0.001]            prostgh2_e       0         [0, 0.002]
    gluala_e         0         [0, 0.001]            cdp_e            0         [-0.001, 0.002]
    gtp_e            0         [0, 0.001]            cmp_e            0         [-0.001, 0.002]
    h2o2_e           0         [0, 0.001]            ctp_e            0         [-0.001, 0.002]
    hista_e          0         [0, 0.001]            gmp_e            0         [-0.001, 0.002]
    hpdca_e          0         [0, 0.001]            octdececoa_c     0         [-0.001, 0.002]
    lnlc_e           0         [0, 0.001]            pmtcoa_r         0         [-0.001, 0.002]
    n2m2nmasn_e      0         [0, 0.001]            fucfucfucgal...  0         [0, 0.00133]
    ocdca_e          0         [0, 0.001]            n5m2masn_g       0         [0, 0.00133]
    octa_e           0         [0, 0.001]            core5_g          0         [0, 0.0012]
    ppa_e            0         [0, 0.001]            core7_g          0         [0, 0.0012]
    pre_prot_r       0         [0, 0.001]            core8_g          0         [0, 0.0012]
    prgstrn_e        0         [0, 0.001]            13_cis_retn_n    0         [0, 0.001]
    prostge2_e       0         [0, 0.001]            3bcrn_e          0         [0, 0.001]
    ptdca_e          0         [0, 0.001]            3ddcrn_e         0         [0, 0.001]
    ptvst_e          0         [0, 0.001]            3deccrn_e        0         [0, 0.001]
    retinol_e        0         [0, 0.001]            3hdececrn_e      0         [0, 0.001]
    spmd_e           0         [0, 0.001]            3hexdcrn_e       0         [0, 0.001]
    tetdece1coa_c    0         [0, 0.001]            3hpvstet_e       0         [0, 0.001]
    tmndnc_e         0         [0, 0.001]            3ivcrn_e         0         [0, 0.001]
    tststerone_e     0         [0, 0.001]            3mlda_e          0         [0, 0.001]
    ttdca_e          0         [0, 0.001]            3octdec2crn_e    0         [0, 0.001]
    whhdca_e         0         [0, 0.001]            3octdeccrn_e     0         [0, 0.001]
    T_antigen_g      0         [-0.0002, 0.001]      3octdece1crn_e   0         [0, 0.001]
    andrstrn_e       0         [-0.001, 0.001]       3tdcrn_e         0         [0, 0.001]
    nrpphr_e         0         [-0.001, 0.001]       3tetd7ecoacrn_e  0         [0, 0.001]
    Asn_X_Ser_Thr_l  0         [0.000333, -0.001]    3thexddcoacrn_e  0         [0, 0.001]
    fol_e            0         [0, -0.001]           3ttetddcoacrn_e  0         [0, 0.001]
    fald_e           0         [0.001, -0.003]       5mta_e           0         [0, 0.001]
    fuc_L_e          0         [0, -0.004]           abt_e            0         [0, 0.001]
    cholate_e        0         [0, -0.378]           aprgstrn_e       0         [0, 0.001]
    strdnc_e         0         [0.001, -0.386]       c4crn_e          0         [0, 0.001]
                                                     c51crn_e         0         [0, 0.001]
                                                     c5dc_e           0         [0, 0.001]
                                                     c6crn_e          0         [0, 0.001]
                                                     c81crn_e         0         [0, 0.001]
                                                     c8crn_e          0         [0, 0.001]
                                                     clpnd_e          0         [0, 0.001]
                                                     ddece1crn_e      0         [0, 0.001]
                                                     decdicrn_e       0         [0, 0.001]
                                                     dem2emgacpai...  0         [0, 0.001]
                                                     dgpi_prot_hs_r   0         [0, 0.001]
                                                     gpi_sig_r        0         [0, 0.001]
                                                     gthox_e          0         [0, 0.001]
                                                     ivcrn_e          0         [0, 0.001]
                                                     mem2emgacpai...  0         [0, 0.001]
                                                     prostgi2_e       0         [0, 0.001]
                                                     ptvstm3_e        0         [0, 0.001]
                                                     ribflv_e         0         [0, 0.001]
                                                     rsvlac_e         0         [0, 0.001]
                                                     sprm_c           0         [0, 0.001]
                                                     sulpacmp_e       0         [0, 0.001]
                                                     tetdece1crn_e    0         [0, 0.001]
                                                     tststeroneglc_e  0         [0, 0.001]
                                                     ttdcrn_e         0         [0, 0.001]
                                                     5adtststerone_e  0         [-0.001, 0.001]
                                                     citr_L_c         0         [-0.001, 0.001]
                                                     citr_L_e         0         [-0.001, 0.001]
                                                     dcsptn1_e        0         [-0.001, 0.001]
                                                     gdp_e            0         [-0.001, 0.001]
                                                     idp_e            0         [-0.001, 0.001]
                                                     itp_e            0         [-0.001, 0.001]



```python
model_HeLaMean5.summary(fva=True)

```

    IN FLUXES                                        OUT FLUXES                                       OBJECTIVES
    -----------------------------------------------  -----------------------------------------------  -----------------------
    id                   Flux  Range                 id                   Flux  Range                 biomass_reac...  0.0266
    ---------------  --------  --------------------  ---------------  --------  --------------------
    glc_D_e          0.293     [0.26, 4.5]           h_e              0.56      [-0.712, 10]
    o2_e             0.07      [0.00272, 10]         lac_L_e          0.505     [-1, 9.62]
    fol_e            0.05      [-0.003, 0.05]        co2_e            0.0661    [-0.001, 6.6]
    ser_L_e          0.042     [0.042, -1.22]        h2o_e            0.0625    [-5.26, 10]
    arg_L_e          0.0252    [-0.0263, 0.084]      pyr_e            0.0503    [-0.11, 10]
    gln_L_e          0.0209    [-0.277, 0.584]       fol_c            0.05      [-0.001, 0.052]
    lys_L_e          0.0158    [0.0158, 0.0168]      ala_L_e          0.0221    [-0.01, 1.35]
    leu_L_e          0.0145    [0.0145, 0.0155]      ptrc_e           0.0167    [-0.001, 0.181]
    val_L_e          0.013     [0.00938, 0.094]      urea_e           0.0157    [0, 0.11]
    phe_L_e          0.0115    [0.0069, 0.066]       glyc_e           0.00983   [-0.001, 1.26]
    asp_L_e          0.01      [0, 0.01]             ura_e            0.00629   [-0.001, 0.00629]
    pro_L_e          0.01      [0.01, -0.0945]       ac_e             0.00606   [-0.001, 0.0632]
    thr_L_e          0.00832   [0.00832, 0.00832]    glyc_S_e         0.00487   [0, 1.2]
    ile_L_e          0.00761   [0.00761, 0.105]      gua_e            0.00461   [-0.001, 0.353]
    asn_L_e          0.00743   [0.01, -0.851]        ade_e            0.00364   [-0.001, 0.319]
    met_L_e          0.00407   [0.00207, 0.03]       gly_e            0.00271   [-0.03, 1.69]
    his_L_e          0.00336   [0.00236, 0.00336]    lnlc_e           0.001     [-0.001, 0.507]
    chol_e           0.00305   [0.00305, 0.05]       chsterol_e       0.001     [0, 0.225]
    pi_e             0.00169   [-0.00184, 1]         thym_e           0.001     [-0.001, 0.00861]
    biomass_other_c  0.00144   [0.00144, 0.00144]    nac_e            0.001     [0, 0.002]
    dag_hs_e         0.001     [0.001, 0.001]        leuktrC4_e       0.001     [-0.001, 0.002]
    lpchol_hs_e      0.001     [0.001, 0.001]        3mlda_e          0.001     [0, 0.001]
    pchol_hs_e       0.001     [0.001, 0.001]        c4crn_e          0.001     [0, 0.001]
    pe_hs_e          0.001     [0.001, 0.001]        lac_D_e          0.001     [0, 0.001]
    pe_hs_r          0.001     [0.001, 0.001]        ribflv_e         0.001     [0, 0.001]
    ps_hs_e          0.001     [0.001, 0.001]        2hb_e            0.000668  [0, 0.0279]
    tag_hs_e         0.001     [0.001, 0.001]        dttp_m           0.000652  [0, 0.00994]
    crn_e            0.001     [0, 0.001]            for_e            0.000543  [-0.001, 0.252]
    dcsptn1_e        0.001     [0, 0.001]            hdcea_e          0.000535  [-0.001, 0.569]
    fmn_e            0.001     [0, 0.001]            tyr_L_e          0.000354  [0.0549, -0.104]
    gam_e            0.001     [0, 0.001]            dgpi_prot_hs_r   0.000333  [0.000333, 0.000333]
    gtp_e            0.001     [0, 0.001]            gpi_sig_r        0.000333  [0.000333, 0.000333]
    hista_e          0.001     [0, 0.001]            Asn_X_Ser_Thr_l  6.7e-05   [-0.001, 0.001]
    o2s_e            0.001     [0, 0.001]            gal_e            0         [-0.001, 2.12]
    orot_e           0.001     [0, 0.001]            abt_e            0         [0, 1.7]
    strdnc_e         0.001     [0, 0.001]            xylt_e           0         [-0.001, 1.7]
    thymd_e          0.001     [0, 0.001]            oxa_e            0         [0, 1.58]
    ncam_c           0.001     [-0.001, 0.001]       fe2_e            0         [-1, 1]
    gdp_e            0.001     [0.001, -0.001]       fe3_e            0         [-1, 1]
    gmp_e            0.001     [0.001, -0.002]       5oxpro_e         0         [0, 0.936]
    leuktrD4_e       0.001     [0.001, -0.002]       vacc_e           0         [-0.001, 0.505]
    cytd_e           0.001     [0.001, -0.00629]     ppi_e            0         [0, 0.501]
    dcmp_e           0.001     [0.001, -0.00629]     ocdcea_e         0         [0, 0.463]
    dctp_n           0.001     [0.001, -0.00629]     ha_e             0         [-0.001, 0.417]
    dcyt_e           0.001     [0.001, -0.00629]     urate_e          0         [0, 0.411]
    udp_e            0.001     [0.001, -0.00629]     nrvnc_e          0         [0, 0.366]
    ump_e            0.001     [0.001, -0.00629]     cholate_e        0         [0, 0.219]
    uri_e            0.001     [0.001, -0.00629]     C02528_e         0         [0, 0.217]
    utp_e            0.001     [0.001, -0.00629]     gchola_e         0         [-0.001, 0.215]
    adn_e            0.001     [0.001, -0.00815]     Rtotal_e         0         [-0.001, 0.171]
    atp_e            0.001     [0.001, -0.00815]     4hphac_e         0         [0, 0.159]
    gsn_e            0.001     [0.001, -0.00815]     tymsf_e          0         [0, 0.159]
    ins_e            0.001     [0.001, -0.00815]     fucfucgalacg...  0         [0, 0.151]
    xmp_e            0.001     [0.001, -0.00815]     galgalfucfuc...  0         [0, 0.0959]
    dttp_n           0.001     [0.001, -0.00894]     glyb_e           0         [0, 0.047]
    datp_n           0.001     [0.001, -0.0098]      gthox_e          0         [0, 0.0456]
    dad_2_e          0.001     [0.001, -0.0139]      citr_L_c         0         [-0.001, 0.0348]
    orn_e            0.001     [0.001, -0.105]       q10h2_e          0         [0, 0.0261]
    xolest2_hs_e     0.001     [0.001, -0.169]       anth_c           0         [0, 0.0156]
    hco3_e           0.001     [0.001, -3.96]        bildglcur_e      0         [0, 0.0151]
    h2o2_e           0.001     [0.001, -10]          bilglcur_e       0         [-0.001, 0.0141]
    inost_e          0.000953  [0.000953, -2.83]     co_e             0         [0, 0.0131]
    man_e            0.0008    [0.001, -1.21]        pheme_e          0         [0, 0.0131]
    ahcys_e          0.000668  [0.001, -0.00815]     octdececoa_c     0         [-0.000333, 0.0102]
    sph1p_e          0.000465  [0.001, -0.498]       35cgmp_e         0         [0, 0.00915]
    pglyc_hs_e       0.000388  [0.000388, 0.000388]  camp_e           0         [0, 0.00915]
    acac_e           0.000384  [0, 0.001]            5mta_e           0         [-0.001, 0.00915]
    trp_L_e          0.000354  [0.000354, 0.016]     spmd_e           0         [-0.001, 0.00915]
    pre_prot_r       0.000333  [0.000333, 0.000333]  ala_B_e          0         [0, 0.00829]
    datp_m           0.000326  [0, 0.001]            sprm_c           0         [0, 0.00508]
    Lcystin_e        0.000285  [0, 0.001]            Ser_Thr_l        0         [-0.001, 0.005]
    dgtp_n           0.000263  [0.000263, 0.000263]  adrn_e           0         [-0.001, 0.003]
    n2m2nmasn_e      6.7e-05   [0, 0.001]            1mncam_e         0         [0, 0.002]
    so4_e            0         [-0.0933, 0.164]      5adtststeron...  0         [0, 0.002]
    cys_L_e          0         [-0.0287, 0.0626]     andrstrnglc_e    0         [0, 0.002]
    glu_L_e          0         [-0.001, 0.01]        estroneglc_e     0         [0, 0.002]
    pnto_R_e         0         [0, 0.00915]          imp_e            0         [0, 0.002]
    3hpvs_e          0         [0, 0.001]            n5m2masn_g       0         [0, 0.002]
    4hpro_LT_m       0         [0, 0.001]            arachd_e         0         [-0.001, 0.002]
    4nph_e           0         [0, 0.001]            core8_g          0         [0, 0.0012]
    5fthf_e          0         [0, 0.001]            3bcrn_e          0         [0, 0.001]
    5mthf_e          0         [0, 0.001]            3hdececrn_e      0         [0, 0.001]
    C02470_e         0         [0, 0.001]            3hpvstet_e       0         [0, 0.001]
    HC00822_e        0         [0, 0.001]            3ivcrn_e         0         [0, 0.001]
    acmp_e           0         [0, 0.001]            3octdec2crn_e    0         [0, 0.001]
    arab_L_e         0         [0, 0.001]            3octdeccrn_e     0         [0, 0.001]
    arach_e          0         [0, 0.001]            3octdece1crn_e   0         [0, 0.001]
    bilirub_e        0         [0, 0.001]            3tetd7ecoacrn_e  0         [0, 0.001]
    carn_e           0         [0, 0.001]            3thexddcoacrn_e  0         [0, 0.001]
    crvnc_e          0         [0, 0.001]            3ttetddcoacrn_e  0         [0, 0.001]
    dca_e            0         [0, 0.001]            4nphsf_e         0         [0, 0.001]
    dheas_e          0         [0, 0.001]            aprgstrn_e       0         [0, 0.001]
    doco13ac_e       0         [0, 0.001]            c5dc_e           0         [0, 0.001]
    dopa_e           0         [0, 0.001]            c6crn_e          0         [0, 0.001]
    etoh_e           0         [0, 0.001]            clpnd_e          0         [0, 0.001]
    gluala_e         0         [0, 0.001]            dmhptcrn_e       0         [0, 0.001]
    hdca_e           0         [0, 0.001]            dopasf_e         0         [0, 0.001]
    hpdca_e          0         [0, 0.001]            epoxtac_e        0         [0, 0.001]
    leuktrA4_e       0         [0, 0.001]            ivcrn_e          0         [0, 0.001]
    lnlnca_e         0         [0, 0.001]            leuktrB4_e       0         [0, 0.001]
    lnlncg_e         0         [0, 0.001]            prostgf2_e       0         [0, 0.001]
    mag_hs_e         0         [0, 0.001]            retn_e           0         [0, 0.001]
    ocdca_e          0         [0, 0.001]            sulpacmp_e       0         [0, 0.001]
    octa_e           0         [0, 0.001]            triodthysuf_e    0         [0, 0.001]
    phyt_e           0         [0, 0.001]            tststeroneglc_e  0         [0, 0.001]
    prgstrn_e        0         [0, 0.001]            dlnlcg_e         0         [-0.001, 0.001]
    prostge2_e       0         [0, 0.001]            estradiol_e      0         [-0.001, 0.001]
    ptdca_e          0         [0, 0.001]            fald_e           0         [-0.001, 0.001]
    retinol_e        0         [0, 0.001]            idp_e            0         [-0.001, 0.001]
    tacr_e           0         [0, 0.001]            itp_e            0         [-0.001, 0.001]
    tmndnc_e         0         [0, 0.001]            meoh_e           0         [-0.001, 0.001]
    triodthy_e       0         [0, 0.001]            taur_c           0         [-0.001, 0.001]
    tststerone_e     0         [0, 0.001]            tchola_e         0         [-0.001, 0.001]
    ttdca_e          0         [0, 0.001]
    whhdca_e         0         [0, 0.001]
    T_antigen_g      0         [-0.0002, 0.001]
    pmtcoa_r         0         [-0.000333, 0.001]
    5adtststerone_e  0         [-0.001, 0.001]
    andrstrn_e       0         [-0.001, 0.001]
    eicostet_e       0         [-0.001, 0.001]
    estrones_e       0         [-0.001, 0.001]
    fuc_L_e          0         [0, -0.188]
    HC01444_e        0         [0.001, -1.03]



```python
FBA1=model_HeLaMean4.optimize()
FBA2=model_KeratinocytesMean4.optimize() 

result = pd.concat([FBA1.fluxes, FBA2.fluxes], axis=1)
result.columns = ["HeLa","Keratinocytes"]
result.to_csv("CellLinesFBA4.csv")

FBA1.fluxes.to_csv("HeLaFBA4.csv")
FBA2.fluxes.to_csv("KeratinocytesFBA4.csv")
```


```python
FBA1=model_HeLaMean5.optimize()
FBA2=model_KeratinocytesMean5.optimize() 

result = pd.concat([FBA1.fluxes, FBA2.fluxes], axis=1)
result.columns = ["HeLa","Keratinocytes"]
result.to_csv("CellLinesFBA5.csv")

FBA1.fluxes.to_csv("HeLaFBA5.csv")
FBA2.fluxes.to_csv("KeratinocytesFBA5.csv")
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

output = mp.Queue(maxsize=12)

# define a example function
def multiple_randomized_analysis(modelC, modelN, Recon, pos, metabolites, nrx=5, penaltyFac=1000):
    cConf_scores=randomize_model(modelC)
    nConf_scores=randomize_model(modelN)

    nConf_scores["EX_arg_L_LPAREN_e_RPAREN_"]=3
    nConf_scores["EX_cys_L_LPAREN_e_RPAREN_"]=3
    nConf_scores["EX_gln_L_LPAREN_e_RPAREN_"]=3
    nConf_scores["EX_gly_LPAREN_e_RPAREN_"]=3
    nConf_scores["EX_his_L_LPAREN_e_RPAREN_"]=3
    nConf_scores["EX_ile_L_LPAREN_e_RPAREN_"]=3
    nConf_scores["EX_leu_L_LPAREN_e_RPAREN_"]=3
    nConf_scores["EX_lys_L_LPAREN_e_RPAREN_"]=3
    nConf_scores["EX_met_L_LPAREN_e_RPAREN_"]=3
    nConf_scores["EX_phe_L_LPAREN_e_RPAREN_"]=3
    nConf_scores["EX_ser_L_LPAREN_e_RPAREN_"]=3
    nConf_scores["EX_thr_L_LPAREN_e_RPAREN_"]=3
    nConf_scores["EX_trp_L_LPAREN_e_RPAREN_"]=3
    nConf_scores["EX_tyr_L_LPAREN_e_RPAREN_"]=3
    nConf_scores["EX_thm_LPAREN_e_RPAREN_"]=3
    nConf_scores["EX_val_L_LPAREN_e_RPAREN_"]=3
    nConf_scores["EX_glc_LPAREN_e_RPAREN_"]=3
    nConf_scores["EX_pyr_LPAREN_e_RPAREN_"]=3
    nConf_scores["EX_ca2_LPAREN_e_RPAREN_"]=3
    nConf_scores["EX_cl_LPAREN_e_RPAREN_"]=3
    nConf_scores["EX_k_LPAREN_e_RPAREN_"]=3
    nConf_scores["EX_na1_LPAREN_e_RPAREN_"]=3
    nConf_scores["EX_fe2_LPAREN_e_RPAREN_"]=3
    nConf_scores["EX_fe3_LPAREN_e_RPAREN_"]=3
    nConf_scores["EX_so4_LPAREN_e_RPAREN_"]=3
    nConf_scores["EX_pi_LPAREN_e_RPAREN_"]=3
    nConf_scores["EX_h_LPAREN_e_RPAREN_"]=3
    nConf_scores["EX_h2o_LPAREN_e_RPAREN_"]=3
    nConf_scores["EX_o2_LPAREN_e_RPAREN_"]=3
    nConf_scores["EX_co2_LPAREN_e_RPAREN_"]=3
    nConf_scores["EX_lac_L_LPAREN_e_RPAREN_"]=3
    nConf_scores["EX_lac_D_LPAREN_e_RPAREN_"]=-1
    nConf_scores["EX_pro_L_LPAREN_e_RPAREN_"]=3
    nConf_scores["EX_ala_L_LPAREN_e_RPAREN_"]=3
    nConf_scores["EX_asn_L_LPAREN_e_RPAREN_"]=3
    nConf_scores["EX_asp_L_LPAREN_e_RPAREN_"]=3
    nConf_scores["EX_glu_L_LPAREN_e_RPAREN_"]=3
    nConf_scores["EX_chol_LPAREN_e_RPAREN_"]=3
    nConf_scores["EX_fol_LPAREN_e_RPAREN_"]=3
    nConf_scores["EX_inost_LPAREN_e_RPAREN_"]=3
    nConf_scores["EX_ncam_LPAREN_e_RPAREN_"]=3
    nConf_scores["EX_pnto_R_LPAREN_e_RPAREN_"]=3
    nConf_scores["EX_pydxn_LPAREN_e_RPAREN_"]=3
    nConf_scores["EX_ribflv_LPAREN_e_RPAREN_"]=3
    
    cConf_scores["EX_arg_L_LPAREN_e_RPAREN_"]=3
    cConf_scores["EX_cys_L_LPAREN_e_RPAREN_"]=3
    cConf_scores["EX_gln_L_LPAREN_e_RPAREN_"]=3
    cConf_scores["EX_gly_LPAREN_e_RPAREN_"]=3
    cConf_scores["EX_his_L_LPAREN_e_RPAREN_"]=3
    cConf_scores["EX_ile_L_LPAREN_e_RPAREN_"]=3
    cConf_scores["EX_leu_L_LPAREN_e_RPAREN_"]=3
    cConf_scores["EX_lys_L_LPAREN_e_RPAREN_"]=3
    cConf_scores["EX_met_L_LPAREN_e_RPAREN_"]=3
    cConf_scores["EX_phe_L_LPAREN_e_RPAREN_"]=3
    cConf_scores["EX_ser_L_LPAREN_e_RPAREN_"]=3
    cConf_scores["EX_thr_L_LPAREN_e_RPAREN_"]=3
    cConf_scores["EX_trp_L_LPAREN_e_RPAREN_"]=3
    cConf_scores["EX_tyr_L_LPAREN_e_RPAREN_"]=3
    cConf_scores["EX_thm_LPAREN_e_RPAREN_"]=3
    cConf_scores["EX_val_L_LPAREN_e_RPAREN_"]=3
    cConf_scores["EX_glc_LPAREN_e_RPAREN_"]=3
    cConf_scores["EX_pyr_LPAREN_e_RPAREN_"]=3
    cConf_scores["EX_ca2_LPAREN_e_RPAREN_"]=3
    cConf_scores["EX_cl_LPAREN_e_RPAREN_"]=3
    cConf_scores["EX_k_LPAREN_e_RPAREN_"]=3
    cConf_scores["EX_na1_LPAREN_e_RPAREN_"]=3
    cConf_scores["EX_fe2_LPAREN_e_RPAREN_"]=3
    cConf_scores["EX_fe3_LPAREN_e_RPAREN_"]=3
    cConf_scores["EX_so4_LPAREN_e_RPAREN_"]=3
    cConf_scores["EX_pi_LPAREN_e_RPAREN_"]=3
    cConf_scores["EX_h_LPAREN_e_RPAREN_"]=3
    cConf_scores["EX_h2o_LPAREN_e_RPAREN_"]=3
    cConf_scores["EX_o2_LPAREN_e_RPAREN_"]=3
    cConf_scores["EX_co2_LPAREN_e_RPAREN_"]=3
    cConf_scores["EX_lac_L_LPAREN_e_RPAREN_"]=3
    cConf_scores["EX_lac_D_LPAREN_e_RPAREN_"]=-1
    cConf_scores["EX_pro_L_LPAREN_e_RPAREN_"]=3
    cConf_scores["EX_ala_L_LPAREN_e_RPAREN_"]=3
    cConf_scores["EX_asn_L_LPAREN_e_RPAREN_"]=3
    cConf_scores["EX_asp_L_LPAREN_e_RPAREN_"]=3
    cConf_scores["EX_glu_L_LPAREN_e_RPAREN_"]=3
    cConf_scores["EX_chol_LPAREN_e_RPAREN_"]=3
    cConf_scores["EX_fol_LPAREN_e_RPAREN_"]=3
    cConf_scores["EX_inost_LPAREN_e_RPAREN_"]=3
    cConf_scores["EX_ncam_LPAREN_e_RPAREN_"]=3
    cConf_scores["EX_pnto_R_LPAREN_e_RPAREN_"]=3
    cConf_scores["EX_pydxn_LPAREN_e_RPAREN_"]=3
    cConf_scores["EX_ribflv_LPAREN_e_RPAREN_"]=3
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
 
    opt_c= CORDA(model=Recon, confidence=cConf_scores, n=nrx, met_prod=metabolites,  penalty_factor=penaltyFac ) 
    opt_n= CORDA(model=Recon, confidence=nConf_scores, n=nrx, met_prod=metabolites,  penalty_factor=penaltyFac ) 
    opt_c.build()
    opt_n.build()
    model_c=opt_c.cobra_model(name="C")
    model_n=opt_n.cobra_model(name="N")
     
    targets=get_drugable_targets(model_n, model_c, "biopsys" )
    f="Targets_pf_1000_Clines_" + str(pos) + ".tsv"                                                                  
    targets.to_csv(f)

# Setup a list of processes that we want to run
pool = mp.Pool(processes=12)

for i in range(1,1001):
    print(i)
    pool.apply_async(multiple_randomized_analysis, args=(conf_HeLaMean,conf_keratinocytesMean,Recon2,i, metas))

pool.close()
pool.join()

```

    1
    2
    3
    4
    5
    6
    7
    8
    9
    10
    11
    12
    13
    14
    15
    16
    17
    18
    19
    20
    21
    22
    23
    24
    25
    26
    27
    28
    29
    30
    31
    32
    33
    34
    35
    36
    37
    38
    39
    40
    41
    42
    43
    44
    45
    46
    47
    48
    49
    50
    51
    52
    53
    54
    55
    56
    57
    58
    59
    60
    61
    62
    63
    64
    65
    66
    67
    68
    69
    70
    71
    72
    73
    74
    75
    76
    77
    78
    79
    80
    81
    82
    83
    84
    85
    86
    87
    88
    89
    90
    91
    92
    93
    94
    95
    96
    97
    98
    99
    100
    101
    102
    103
    104
    105
    106
    107
    108
    109
    110
    111
    112
    113
    114
    115
    116
    117
    118
    119
    120
    121
    122
    123
    124
    125
    126
    127
    128
    129
    130
    131
    132
    133
    134
    135
    136
    137
    138
    139
    140
    141
    142
    143
    144
    145
    146
    147
    148
    149
    150
    151
    152
    153
    154
    155
    156
    157
    158
    159
    160
    161
    162
    163
    164
    165
    166
    167
    168
    169
    170
    171
    172
    173
    174
    175
    176
    177
    178
    179
    180
    181
    182
    183
    184
    185
    186
    187
    188
    189
    190
    191
    192
    193
    194
    195
    196
    197
    198
    199
    200
    201
    202
    203
    204
    205
    206
    207
    208
    209
    210
    211
    212
    213
    214
    215
    216
    217
    218
    219
    220
    221
    222
    223
    224
    225
    226
    227
    228
    229
    230
    231
    232
    233
    234
    235
    236
    237
    238
    239
    240
    241
    242
    243
    244
    245
    246
    247
    248
    249
    250
    251
    252
    253
    254
    255
    256
    257
    258
    259
    260
    261
    262
    263
    264
    265
    266
    267
    268
    269
    270
    271
    272
    273
    274
    275
    276
    277
    278
    279
    280
    281
    282
    283
    284
    285
    286
    287
    288
    289
    290
    291
    292
    293
    294
    295
    296
    297
    298
    299
    300
    301
    302
    303
    304
    305
    306
    307
    308
    309
    310
    311
    312
    313
    314
    315
    316
    317
    318
    319
    320
    321
    322
    323
    324
    325
    326
    327
    328
    329
    330
    331
    332
    333
    334
    335
    336
    337
    338
    339
    340
    341
    342
    343
    344
    345
    346
    347
    348
    349
    350
    351
    352
    353
    354
    355
    356
    357
    358
    359
    360
    361
    362
    363
    364
    365
    366
    367
    368
    369
    370
    371
    372
    373
    374
    375
    376
    377
    378
    379
    380
    381
    382
    383
    384
    385
    386
    387
    388
    389
    390
    391
    392
    393
    394
    395
    396
    397
    398
    399
    400
    401
    402
    403
    404
    405
    406
    407
    408
    409
    410
    411
    412
    413
    414
    415
    416
    417
    418
    419
    420
    421
    422
    423
    424
    425
    426
    427
    428
    429
    430
    431
    432
    433
    434
    435
    436
    437
    438
    439
    440
    441
    442
    443
    444
    445
    446
    447
    448
    449
    450
    451
    452
    453
    454
    455
    456
    457
    458
    459
    460
    461
    462
    463
    464
    465
    466
    467
    468
    469
    470
    471
    472
    473
    474
    475
    476
    477
    478
    479
    480
    481
    482
    483
    484
    485
    486
    487
    488
    489
    490
    491
    492
    493
    494
    495
    496
    497
    498
    499
    500
    501
    502
    503
    504
    505
    506
    507
    508
    509
    510
    511
    512
    513
    514
    515
    516
    517
    518
    519
    520
    521
    522
    523
    524
    525
    526
    527
    528
    529
    530
    531
    532
    533
    534
    535
    536
    537
    538
    539
    540
    541
    542
    543
    544
    545
    546
    547
    548
    549
    550
    551
    552
    553
    554
    555
    556
    557
    558
    559
    560
    561
    562
    563
    564
    565
    566
    567
    568
    569
    570
    571
    572
    573
    574
    575
    576
    577
    578
    579
    580
    581
    582
    583
    584
    585
    586
    587
    588
    589
    590
    591
    592
    593
    594
    595
    596
    597
    598
    599
    600
    601
    602
    603
    604
    605
    606
    607
    608
    609
    610
    611
    612
    613
    614
    615
    616
    617
    618
    619
    620
    621
    622
    623
    624
    625
    626
    627
    628
    629
    630
    631
    632
    633
    634
    635
    636
    637
    638
    639
    640
    641
    642
    643
    644
    645
    646
    647
    648
    649
    650
    651
    652
    653
    654
    655
    656
    657
    658
    659
    660
    661
    662
    663
    664
    665
    666
    667
    668
    669
    670
    671
    672
    673
    674
    675
    676
    677
    678
    679
    680
    681
    682
    683
    684
    685
    686
    687
    688
    689
    690
    691
    692
    693
    694
    695
    696
    697
    698
    699
    700
    701
    702
    703
    704
    705
    706
    707
    708
    709
    710
    711
    712
    713
    714
    715
    716
    717
    718
    719
    720
    721
    722
    723
    724
    725
    726
    727
    728
    729
    730
    731
    732
    733
    734
    735
    736
    737
    738
    739
    740
    741
    742
    743
    744
    745
    746
    747
    748
    749
    750
    751
    752
    753
    754
    755
    756
    757
    758
    759
    760
    761
    762
    763
    764
    765
    766
    767
    768
    769
    770
    771
    772
    773
    774
    775
    776
    777
    778
    779
    780
    781
    782
    783
    784
    785
    786
    787
    788
    789
    790
    791
    792
    793
    794
    795
    796
    797
    798
    799
    800
    801
    802
    803
    804
    805
    806
    807
    808
    809
    810
    811
    812
    813
    814
    815
    816
    817
    818
    819
    820
    821
    822
    823
    824
    825
    826
    827
    828
    829
    830
    831
    832
    833
    834
    835
    836
    837
    838
    839



```python
cobra.io.write_sbml_model(model_HeLaMean1, "Full_HeLa_model_0219_DEMEM6429_n10_pf100_v1.sbml")
cobra.io.write_sbml_model(model_KeratinocytesMean1, "Full_Kerat_model_0219_DEMEM6429_n10_pf5000_v1.sbml")

cobra.io.write_sbml_model(model_HeLaMean2, "Full_HeLa_model_0219_DEMEM6429_n10_pf500_v2.sbml")
cobra.io.write_sbml_model(model_KeratinocytesMean2, "Full_Kerat_model_0219_DEMEM6429_n10_pf5000_v2.sbml")

cobra.io.write_sbml_model(model_HeLaMean3, "Full_HeLa_model_0219_DEMEM6429_n10_pf1000_v3.sbml")
cobra.io.write_sbml_model(model_KeratinocytesMean3, "Full_Kerat_model_0219_DEMEM6429_n10_pf5000_v3.sbml")

cobra.io.write_sbml_model(model_HeLaMean4, "Full_HeLa_model_0219_DEMEM6429_n10_pf2500_v4.sbml")
cobra.io.write_sbml_model(model_KeratinocytesMean4, "Full_Kerat_model_0219_DEMEM6429_n10_pf5000_v4.sbml")

cobra.io.write_sbml_model(model_HeLaMean5, "Full_HeLa_model_0219_DEMEM6429_n10_pf5000_v5.sbml")
cobra.io.write_sbml_model(model_KeratinocytesMean5, "Full_Kerat_model_0219_DEMEM6429_n10_pf5000_v5.sbml")
```



```python
print(model_HeLaMean1.reactions.index)
```


```python
print(model_HeLaMean2.reactions.index)
```


```python
from operator import attrgetter

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
        FBA = model.optimize(error_value=None)
        solution = FBA

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

