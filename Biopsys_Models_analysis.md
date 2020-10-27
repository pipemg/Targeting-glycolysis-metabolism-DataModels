

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
    rxex.bounds=(0.0,10)

Recon2.reactions.get_by_id("DM_4hrpo").lower_bound=-0.012
Recon2.reactions.get_by_id("DM_datp_n_").lower_bound=-0.012
Recon2.reactions.get_by_id("DM_dctp_n_").lower_bound=-0.012
Recon2.reactions.get_by_id("DM_dgtp_n_").lower_bound=-0.012 
Recon2.reactions.get_by_id("DM_dttp_n_").lower_bound=-0.012 
Recon2.reactions.get_by_id("DM_Lcystin").lower_bound=-0.024
Recon2.reactions.get_by_id("DM_pe_hs_LPAREN_r_RPAREN_").lower_bound=-0.012
Recon2.reactions.get_by_id("EX_2hb_LPAREN_e_RPAREN_").lower_bound=-0.012
Recon2.reactions.get_by_id("EX_34hpp_").lower_bound=-0.012
Recon2.reactions.get_by_id("EX_3hpvs_LPAREN_e_RPAREN_").lower_bound=-0.012
Recon2.reactions.get_by_id("EX_3mob_LPAREN_e_RPAREN_").lower_bound=-0.012
Recon2.reactions.get_by_id("EX_4mop_LPAREN_e_RPAREN_").lower_bound=-0.012
Recon2.reactions.get_by_id("EX_ac_LPAREN_e_RPAREN_").lower_bound=-0.012
Recon2.reactions.get_by_id("EX_acac_LPAREN_e_RPAREN_").lower_bound=-0.012
Recon2.reactions.get_by_id("EX_acald_LPAREN_e_RPAREN_").lower_bound=-0.012
Recon2.reactions.get_by_id("EX_acetone_LPAREN_e_RPAREN_").lower_bound=-0.012
Recon2.reactions.get_by_id("EX_acgam_LPAREN_e_RPAREN_").lower_bound=-0.012 
Recon2.reactions.get_by_id("EX_acmana_LPAREN_e_RPAREN_").lower_bound=-0.012
Recon2.reactions.get_by_id("EX_ade_LPAREN_e_RPAREN_").lower_bound=-0.012
Recon2.reactions.get_by_id("EX_adn_LPAREN_e_RPAREN_").lower_bound=-0.012
Recon2.reactions.get_by_id("EX_adpcbl_LPAREN_e_RPAREN_").lower_bound=-0.012
Recon2.reactions.get_by_id("EX_akg_LPAREN_e_RPAREN_").lower_bound=-0.012
Recon2.reactions.get_by_id("EX_ala_B_LPAREN_e_RPAREN_").lower_bound=-0.012
Recon2.reactions.get_by_id("EX_ala_D_LPAREN_e_RPAREN_").lower_bound=-0.04
Recon2.reactions.get_by_id("EX_ala_L_LPAREN_e_RPAREN_").lower_bound=-0.04
Recon2.reactions.get_by_id("EX_am9csa_LPAREN_e_RPAREN_").lower_bound=-0.012
Recon2.reactions.get_by_id("EX_amp_LPAREN_e_RPAREN_").lower_bound=-0.012
Recon2.reactions.get_by_id("EX_arab_L_LPAREN_e_RPAREN_").lower_bound=-0.012
Recon2.reactions.get_by_id("EX_arachd_LPAREN_e_RPAREN_").lower_bound=-0.012 
Recon2.reactions.get_by_id("EX_arg_L_LPAREN_e_RPAREN_").lower_bound=-0.012
Recon2.reactions.get_by_id("EX_asn_L_LPAREN_e_RPAREN_").lower_bound=-0.012
Recon2.reactions.get_by_id("EX_asp_L_LPAREN_e_RPAREN_").lower_bound=-0.012
Recon2.reactions.get_by_id("EX_atp_LPAREN_e_RPAREN_").lower_bound=-0.8
Recon2.reactions.get_by_id("EX_bhb_LPAREN_e_RPAREN_").lower_bound=-0.012
Recon2.reactions.get_by_id("EX_btn_LPAREN_e_RPAREN_").lower_bound=-0.2
Recon2.reactions.get_by_id("EX_ca2_LPAREN_e_RPAREN_").lower_bound=-0.2
Recon2.reactions.get_by_id("EX_carn_LPAREN_e_RPAREN_").lower_bound=-0.012
Recon2.reactions.get_by_id("EX_caro_LPAREN_e_RPAREN_").lower_bound=-0.012
Recon2.reactions.get_by_id("EX_cgly_LPAREN_e_RPAREN_").lower_bound=-0.012
Recon2.reactions.get_by_id("EX_chol_LPAREN_e_RPAREN_").lower_bound=-0.012
Recon2.reactions.get_by_id("EX_chsterol_LPAREN_e_RPAREN_").lower_bound=-2
Recon2.reactions.get_by_id("EX_cit_LPAREN_e_RPAREN_").lower_bound=-0.012
Recon2.reactions.get_by_id("EX_cl_LPAREN_e_RPAREN_").lower_bound=-1
Recon2.reactions.get_by_id("EX_CLPND_LPAREN_e_RPAREN_").lower_bound=-0.012
Recon2.reactions.get_by_id("EX_co2_LPAREN_e_RPAREN_").lower_bound=0
Recon2.reactions.get_by_id("EX_creat_LPAREN_e_RPAREN_").lower_bound=-0.012
Recon2.reactions.get_by_id("EX_crn_LPAREN_e_RPAREN_").lower_bound=-0.012
Recon2.reactions.get_by_id("EX_crvnc_LPAREN_e_RPAREN_").lower_bound=-0.012
Recon2.reactions.get_by_id("EX_csa_LPAREN_e_RPAREN_").lower_bound=-0.012
Recon2.reactions.get_by_id("EX_csn_LPAREN_e_RPAREN_").lower_bound=-0.012
Recon2.reactions.get_by_id("EX_cys_L_LPAREN_e_RPAREN_").lower_bound=-0.012
Recon2.reactions.get_by_id("EX_cytd_LPAREN_e_RPAREN_").lower_bound=-0.012
Recon2.reactions.get_by_id("EX_dad_2_LPAREN_e_RPAREN_").lower_bound=-0.012
Recon2.reactions.get_by_id("EX_dag_hs_LPAREN_e_RPAREN_").lower_bound=-0.012
Recon2.reactions.get_by_id("EX_dcmp_LPAREN_e_RPAREN_").lower_bound=-0.012
Recon2.reactions.get_by_id("EX_dcyt_LPAREN_e_RPAREN_").lower_bound=-0.012
Recon2.reactions.get_by_id("EX_ddca_LPAREN_e_RPAREN_").lower_bound=-0.012
Recon2.reactions.get_by_id("EX_dgsn_LPAREN_e_RPAREN_").lower_bound=-0.012
Recon2.reactions.get_by_id("EX_dhdascb_LPAREN_e_RPAREN_").lower_bound=-0.012
Recon2.reactions.get_by_id("EX_din_LPAREN_e_RPAREN_").lower_bound=-0.012
Recon2.reactions.get_by_id("EX_dopa_LPAREN_e_RPAREN_").lower_bound=-0.012
Recon2.reactions.get_by_id("EX_drib_LPAREN_e_RPAREN_").lower_bound=-0.012
Recon2.reactions.get_by_id("EX_etoh_LPAREN_e_RPAREN_").lower_bound=-0.012
Recon2.reactions.get_by_id("EX_fald_LPAREN_e_RPAREN_").lower_bound=-0.012
Recon2.reactions.get_by_id("EX_fe2_LPAREN_e_RPAREN_").lower_bound=-0.012
Recon2.reactions.get_by_id("EX_fe2_LPAREN_e_RPAREN_").lower_bound=-0.545
Recon2.reactions.get_by_id("EX_fe3_LPAREN_e_RPAREN_").lower_bound=-0.012
Recon2.reactions.get_by_id("EX_fmn_LPAREN_e_RPAREN_").lower_bound=-0.012
Recon2.reactions.get_by_id("EX_fol_LPAREN_e_RPAREN_").lower_bound=-0.012
Recon2.reactions.get_by_id("EX_for_LPAREN_e_RPAREN_").lower_bound=-0.012
Recon2.reactions.get_by_id("EX_fru_LPAREN_e_RPAREN_").lower_bound=-0.012
Recon2.reactions.get_by_id("EX_fuc_L_LPAREN_e_RPAREN_").lower_bound=-0.012
Recon2.reactions.get_by_id("EX_fum_LPAREN_e_RPAREN_").lower_bound=-0.012
Recon2.reactions.get_by_id("EX_gal_LPAREN_e_RPAREN_").lower_bound=-0.012
Recon2.reactions.get_by_id("EX_gam_LPAREN_e_RPAREN_").lower_bound=-0.012
Recon2.reactions.get_by_id("EX_gchola_LPAREN_e_RPAREN_").lower_bound=-0.012
Recon2.reactions.get_by_id("EX_glc_LPAREN_e_RPAREN_").lower_bound=-1
Recon2.reactions.get_by_id("EX_glcur_LPAREN_e_RPAREN_").lower_bound=-0.032
Recon2.reactions.get_by_id("EX_gln_L_LPAREN_e_RPAREN_").lower_bound=-0.09
Recon2.reactions.get_by_id("EX_glu_L_LPAREN_e_RPAREN_").lower_bound=-0.012
Recon2.reactions.get_by_id("EX_gluala_LPAREN_e_RPAREN_").lower_bound=-0.012
Recon2.reactions.get_by_id("EX_gly_LPAREN_e_RPAREN_").lower_bound=-0.03
Recon2.reactions.get_by_id("EX_glyb_LPAREN_e_RPAREN_").lower_bound=-0.012
Recon2.reactions.get_by_id("EX_glyc_LPAREN_e_RPAREN_").lower_bound=-0.21
Recon2.reactions.get_by_id("EX_glyc3p_LPAREN_e_RPAREN_").lower_bound=-0.012
Recon2.reactions.get_by_id("EX_glygn2_LPAREN_e_RPAREN_").lower_bound=-0.012
Recon2.reactions.get_by_id("EX_gsn_LPAREN_e_RPAREN_").lower_bound=-0.012
Recon2.reactions.get_by_id("EX_gthox_LPAREN_e_RPAREN_").lower_bound=-0.012
Recon2.reactions.get_by_id("EX_gthrd_LPAREN_e_RPAREN_").lower_bound=-0.09
Recon2.reactions.get_by_id("EX_gua_LPAREN_e_RPAREN_").lower_bound=-0.012
Recon2.reactions.get_by_id("EX_h_LPAREN_e_RPAREN_").lower_bound=-1
Recon2.reactions.get_by_id("EX_h2o_LPAREN_e_RPAREN_").lower_bound=-10
Recon2.reactions.get_by_id("EX_ha_LPAREN_e_RPAREN_").lower_bound=-0.012
Recon2.reactions.get_by_id("EX_HC00250_LPAREN_e_RPAREN_").lower_bound=-0.012
Recon2.reactions.get_by_id("EX_HC01609_LPAREN_e_RPAREN_").lower_bound=-0.012
Recon2.reactions.get_by_id("EX_HC01610_LPAREN_e_RPAREN_").lower_bound=-0.012
Recon2.reactions.get_by_id("EX_hdca_LPAREN_e_RPAREN_").lower_bound=-0.3
Recon2.reactions.get_by_id("EX_hdcea_LPAREN_e_RPAREN_").lower_bound=-0.012
Recon2.reactions.get_by_id("EX_his_L_LPAREN_e_RPAREN_").lower_bound=-0.012
Recon2.reactions.get_by_id("EX_hxan_LPAREN_e_RPAREN_").lower_bound=-0.012
Recon2.reactions.get_by_id("EX_ile_L_LPAREN_e_RPAREN_").lower_bound=-0.012
Recon2.reactions.get_by_id("EX_inost_LPAREN_e_RPAREN_").lower_bound=-0.012
Recon2.reactions.get_by_id("EX_k_LPAREN_e_RPAREN_").lower_bound=-0.2
Recon2.reactions.get_by_id("EX_lac_L_LPAREN_e_RPAREN_").lower_bound=-0.09
Recon2.reactions.get_by_id("EX_lcts_LPAREN_e_RPAREN_").lower_bound=-0.5
Recon2.reactions.get_by_id("EX_leu_L_LPAREN_e_RPAREN_").lower_bound=-0.012
Recon2.reactions.get_by_id("EX_leuktrA4_LPAREN_e_RPAREN_").lower_bound=-0.012
Recon2.reactions.get_by_id("EX_leuktrD4_LPAREN_e_RPAREN_").lower_bound=-0.012
Recon2.reactions.get_by_id("EX_leuktrE4_LPAREN_e_RPAREN_").lower_bound=-0.012
Recon2.reactions.get_by_id("EX_lnlc_LPAREN_e_RPAREN_").lower_bound=-0.06
Recon2.reactions.get_by_id("EX_lnlnca_LPAREN_e_RPAREN_").lower_bound=-0.012
Recon2.reactions.get_by_id("EX_lnlncg_LPAREN_e_RPAREN_").lower_bound=-0.012
Recon2.reactions.get_by_id("EX_lpchol_hs_LPAREN_e_RPAREN_").lower_bound=-0.012
Recon2.reactions.get_by_id("EX_lys_L_LPAREN_e_RPAREN_").lower_bound=-0.03 
Recon2.reactions.get_by_id("EX_mag_hs_LPAREN_e_RPAREN_").lower_bound=-0.012
Recon2.reactions.get_by_id("EX_mal_L_LPAREN_e_RPAREN_").lower_bound=-0.012
Recon2.reactions.get_by_id("EX_malt_LPAREN_e_RPAREN_").lower_bound=-0.012 
Recon2.reactions.get_by_id("EX_man_LPAREN_e_RPAREN_").lower_bound=-0.012 
Recon2.reactions.get_by_id("EX_meoh_LPAREN_e_RPAREN_").lower_bound=-0.012
Recon2.reactions.get_by_id("EX_met_L_LPAREN_e_RPAREN_").lower_bound=-0.012
Recon2.reactions.get_by_id("EX_n2m2nmasn_LPAREN_e_RPAREN_").lower_bound=-0.012
Recon2.reactions.get_by_id("EX_na1_LPAREN_e_RPAREN_").lower_bound=-1
Recon2.reactions.get_by_id("EX_nac_LPAREN_e_RPAREN_").lower_bound=-0.012
Recon2.reactions.get_by_id("EX_ncam_LPAREN_e_RPAREN_").lower_bound=-0.012
Recon2.reactions.get_by_id("EX_nh4_LPAREN_e_RPAREN_").lower_bound=-0.012
Recon2.reactions.get_by_id("EX_no2_LPAREN_e_RPAREN_").lower_bound=-0.012
Recon2.reactions.get_by_id("EX_o2_LPAREN_e_RPAREN_").lower_bound=-1
Recon2.reactions.get_by_id("EX_o2s_LPAREN_e_RPAREN_").lower_bound=-0.012
Recon2.reactions.get_by_id("EX_ocdca_LPAREN_e_RPAREN_").lower_bound=-0.07
Recon2.reactions.get_by_id("EX_ocdcea_LPAREN_e_RPAREN_").lower_bound=-0.012 ### ??
Recon2.reactions.get_by_id("EX_octa_LPAREN_e_RPAREN_").lower_bound=-0.012
Recon2.reactions.get_by_id("EX_orn_LPAREN_e_RPAREN_").lower_bound=-0.012
Recon2.reactions.get_by_id("EX_oxa_LPAREN_e_RPAREN_").lower_bound=-0.012
Recon2.reactions.get_by_id("EX_pe_hs_LPAREN_e_RPAREN_").lower_bound=-0.012
Recon2.reactions.get_by_id("EX_pglyc_hs_LPAREN_e_RPAREN_").lower_bound=-0.012
Recon2.reactions.get_by_id("EX_phe_L_LPAREN_e_RPAREN_").lower_bound=-0.012
Recon2.reactions.get_by_id("EX_pheme_LPAREN_e_RPAREN_").lower_bound=-0.012
Recon2.reactions.get_by_id("EX_pi_LPAREN_e_RPAREN_").lower_bound=-0.04
Recon2.reactions.get_by_id("EX_pnto_R_LPAREN_e_RPAREN_").lower_bound=-0.012
Recon2.reactions.get_by_id("EX_ppa_LPAREN_e_RPAREN_").lower_bound=-0.012
Recon2.reactions.get_by_id("EX_pro_L_LPAREN_e_RPAREN_").lower_bound=-0.012
Recon2.reactions.get_by_id("EX_prostgh2_LPAREN_e_RPAREN_").lower_bound=-0.012
Recon2.reactions.get_by_id("EX_ps_hs_LPAREN_e_RPAREN_").lower_bound=-0.012
Recon2.reactions.get_by_id("EX_ptrc_LPAREN_e_RPAREN_").lower_bound=-0.012
Recon2.reactions.get_by_id("EX_ptvstlac_LPAREN_e_RPAREN_").lower_bound=-0.012
Recon2.reactions.get_by_id("EX_pydam_LPAREN_e_RPAREN_").lower_bound=-1
Recon2.reactions.get_by_id("EX_pydx_LPAREN_e_RPAREN_").lower_bound=-0.012
Recon2.reactions.get_by_id("EX_pydx5p_LPAREN_e_RPAREN_").lower_bound=-0.012
Recon2.reactions.get_by_id("EX_pydxn_LPAREN_e_RPAREN_").lower_bound=-0.012
Recon2.reactions.get_by_id("EX_pyr_LPAREN_e_RPAREN_").lower_bound=-0.012
Recon2.reactions.get_by_id("EX_q10h2_LPAREN_e_RPAREN_").lower_bound=-0.012
Recon2.reactions.get_by_id("EX_retfa_LPAREN_e_RPAREN_").lower_bound=-0.012
Recon2.reactions.get_by_id("EX_retinol_LPAREN_e_RPAREN_").lower_bound=-0.012
Recon2.reactions.get_by_id("EX_retn_LPAREN_e_RPAREN_").lower_bound=-0.012
Recon2.reactions.get_by_id("EX_rib_D_LPAREN_e_RPAREN_").lower_bound=-0.012
Recon2.reactions.get_by_id("EX_ribflv_LPAREN_e_RPAREN_").lower_bound=-0.012
Recon2.reactions.get_by_id("EX_sbt_DASH_d_LPAREN_e_RPAREN_").lower_bound=-0.012
Recon2.reactions.get_by_id("EX_sel_LPAREN_e_RPAREN_").lower_bound=-0.012
Recon2.reactions.get_by_id("EX_ser_L_LPAREN_e_RPAREN_").lower_bound=-0.012
Recon2.reactions.get_by_id("EX_so4_LPAREN_e_RPAREN_").lower_bound=-0.05
Recon2.reactions.get_by_id("EX_sph1p_LPAREN_e_RPAREN_").lower_bound=-0.012
Recon2.reactions.get_by_id("EX_spmd_LPAREN_e_RPAREN_").lower_bound=-0.012
Recon2.reactions.get_by_id("EX_strch1_LPAREN_e_RPAREN_").lower_bound=-0.012
Recon2.reactions.get_by_id("EX_strch2_LPAREN_e_RPAREN_").lower_bound=-0.012
Recon2.reactions.get_by_id("EX_sucr_LPAREN_e_RPAREN_").lower_bound=-0.012
Recon2.reactions.get_by_id("EX_tag_hs_LPAREN_e_RPAREN_").lower_bound=-0.012
Recon2.reactions.get_by_id("EX_thm_LPAREN_e_RPAREN_").lower_bound=-0.012
Recon2.reactions.get_by_id("EX_thr_L_LPAREN_e_RPAREN_").lower_bound=-0.012
Recon2.reactions.get_by_id("EX_thymd_LPAREN_e_RPAREN_").lower_bound=-0.012
Recon2.reactions.get_by_id("EX_tre_LPAREN_e_RPAREN_").lower_bound=-0.012
Recon2.reactions.get_by_id("EX_trp_L_LPAREN_e_RPAREN_").lower_bound=-0.012
Recon2.reactions.get_by_id("EX_ttdca_LPAREN_e_RPAREN_").lower_bound=-0.012
Recon2.reactions.get_by_id("EX_tyr_L_LPAREN_e_RPAREN_").lower_bound=-0.012
Recon2.reactions.get_by_id("EX_ura_LPAREN_e_RPAREN_").lower_bound=-0.012
Recon2.reactions.get_by_id("EX_urea_LPAREN_e_RPAREN_").lower_bound=-0.4
Recon2.reactions.get_by_id("EX_uri_LPAREN_e_RPAREN_").lower_bound=-0.012
Recon2.reactions.get_by_id("EX_utp_LPAREN_e_RPAREN_").lower_bound=-0.02 
Recon2.reactions.get_by_id("EX_val_L_LPAREN_e_RPAREN_").lower_bound=-0.03
Recon2.reactions.get_by_id("EX_xmp_LPAREN_e_RPAREN_").lower_bound=-0.012
Recon2.reactions.get_by_id("EX_xyl_D_LPAREN_e_RPAREN_").lower_bound=-0.3
Recon2.reactions.get_by_id("EX_xylt_LPAREN_e_RPAREN_").lower_bound=-0.012


Recon2.reactions.biomass_reaction.bounds=(0.008, 0.05)
Recon2.reactions.EX_lac_D_LPAREN_e_RPAREN_.bounds=(-0.001,0.001)
```


```python
Recon2.reactions.get_by_id("biomass_DNA").reaction
```




    '0.941642857142857 datp_n + 0.674428571428572 dctp_n + 0.707 dgtp_n + 0.935071428571429 dttp_n --> biomass_DNA_c'




```python
import pandas as pd
ub_scores_matrix =pd.read_csv("Binary_Byposys_RECON2_2.tsv",sep=",",index_col=0,header=0)
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


conf_CancerBiopsy = {}
conf_NormalBiopsy = {}

for r in Recon2.reactions:
    if(r.gene_reaction_rule!=''):
        conf_CancerBiopsy[r.id]=reaction_confidence(r.gene_reaction_rule,confidence_scores_matrix["MaxCancerBiopsy"])
        conf_NormalBiopsy[r.id]=reaction_confidence(r.gene_reaction_rule,confidence_scores_matrix["MaxNormalBiopsy"])
    else:
        conf_CancerBiopsy[r.id]=1 
        conf_NormalBiopsy[r.id]=1
```


```python
%matplotlib inline
import pandas as pd

df=pd.DataFrame({'Cancer': conf_CancerBiopsy})
df.Cancer.value_counts()

```




     1.0    3450
     3.0    1441
     0.0    1381
    -1.0    1275
     2.0     238
    Name: Cancer, dtype: int64




```python
conf_CancerBiopsy["DM_4hrpo"]=3
conf_CancerBiopsy["DM_datp_n_"]=3
conf_CancerBiopsy["DM_dctp_n_"]=3
conf_CancerBiopsy["DM_dgtp_n_"]=3
conf_CancerBiopsy["DM_dttp_n_"]=3
conf_CancerBiopsy["DM_Lcystin"]=3
conf_CancerBiopsy["DM_pe_hs_LPAREN_r_RPAREN_"]=3
conf_CancerBiopsy["EX_2hb_LPAREN_e_RPAREN_"]=3
conf_CancerBiopsy["EX_34hpp_"]=3
conf_CancerBiopsy["EX_3hpvs_LPAREN_e_RPAREN_"]=3
conf_CancerBiopsy["EX_3mob_LPAREN_e_RPAREN_"]=3
conf_CancerBiopsy["EX_4mop_LPAREN_e_RPAREN_"]=3
conf_CancerBiopsy["EX_ac_LPAREN_e_RPAREN_"]=3
conf_CancerBiopsy["EX_acac_LPAREN_e_RPAREN_"]=3
conf_CancerBiopsy["EX_acald_LPAREN_e_RPAREN_"]=3
conf_CancerBiopsy["EX_acetone_LPAREN_e_RPAREN_"]=3
conf_CancerBiopsy["EX_acgam_LPAREN_e_RPAREN_"]=3
conf_CancerBiopsy["EX_acmana_LPAREN_e_RPAREN_"]=3
conf_CancerBiopsy["EX_ade_LPAREN_e_RPAREN_"]=3
conf_CancerBiopsy["EX_adn_LPAREN_e_RPAREN_"]=3
conf_CancerBiopsy["EX_adpcbl_LPAREN_e_RPAREN_"]=3
conf_CancerBiopsy["EX_akg_LPAREN_e_RPAREN_"]=3
conf_CancerBiopsy["EX_ala_B_LPAREN_e_RPAREN_"]=3
conf_CancerBiopsy["EX_ala_D_LPAREN_e_RPAREN_"]=3
conf_CancerBiopsy["EX_ala_L_LPAREN_e_RPAREN_"]=3
conf_CancerBiopsy["EX_am9csa_LPAREN_e_RPAREN_"]=3
conf_CancerBiopsy["EX_amp_LPAREN_e_RPAREN_"]=3
conf_CancerBiopsy["EX_arab_L_LPAREN_e_RPAREN_"]=3
conf_CancerBiopsy["EX_arachd_LPAREN_e_RPAREN_"]=3
conf_CancerBiopsy["EX_arg_L_LPAREN_e_RPAREN_"]=3
conf_CancerBiopsy["EX_asn_L_LPAREN_e_RPAREN_"]=3
conf_CancerBiopsy["EX_asp_L_LPAREN_e_RPAREN_"]=3
conf_CancerBiopsy["EX_atp_LPAREN_e_RPAREN_"]=3
conf_CancerBiopsy["EX_bhb_LPAREN_e_RPAREN_"]=3
conf_CancerBiopsy["EX_btn_LPAREN_e_RPAREN_"]=3
conf_CancerBiopsy["EX_ca2_LPAREN_e_RPAREN_"]=3
conf_CancerBiopsy["EX_carn_LPAREN_e_RPAREN_"]=3
conf_CancerBiopsy["EX_caro_LPAREN_e_RPAREN_"]=3
conf_CancerBiopsy["EX_cgly_LPAREN_e_RPAREN_"]=3
conf_CancerBiopsy["EX_chol_LPAREN_e_RPAREN_"]=3
conf_CancerBiopsy["EX_chsterol_LPAREN_e_RPAREN_"]=3
conf_CancerBiopsy["EX_cit_LPAREN_e_RPAREN_"]=3
conf_CancerBiopsy["EX_cl_LPAREN_e_RPAREN_"]=3
conf_CancerBiopsy["EX_CLPND_LPAREN_e_RPAREN_"]=3
conf_CancerBiopsy["EX_co2_LPAREN_e_RPAREN_"]=3
conf_CancerBiopsy["EX_creat_LPAREN_e_RPAREN_"]=3
conf_CancerBiopsy["EX_crn_LPAREN_e_RPAREN_"]=3
conf_CancerBiopsy["EX_crvnc_LPAREN_e_RPAREN_"]=3
conf_CancerBiopsy["EX_csa_LPAREN_e_RPAREN_"]=3
conf_CancerBiopsy["EX_csn_LPAREN_e_RPAREN_"]=3
conf_CancerBiopsy["EX_cys_L_LPAREN_e_RPAREN_"]=3
conf_CancerBiopsy["EX_cytd_LPAREN_e_RPAREN_"]=3
conf_CancerBiopsy["EX_dad_2_LPAREN_e_RPAREN_"]=3
conf_CancerBiopsy["EX_dag_hs_LPAREN_e_RPAREN_"]=3
conf_CancerBiopsy["EX_dcmp_LPAREN_e_RPAREN_"]=3
conf_CancerBiopsy["EX_dcyt_LPAREN_e_RPAREN_"]=3
conf_CancerBiopsy["EX_ddca_LPAREN_e_RPAREN_"]=3
conf_CancerBiopsy["EX_dgsn_LPAREN_e_RPAREN_"]=3
conf_CancerBiopsy["EX_dhdascb_LPAREN_e_RPAREN_"]=3
conf_CancerBiopsy["EX_din_LPAREN_e_RPAREN_"]=3
conf_CancerBiopsy["EX_dopa_LPAREN_e_RPAREN_"]=3
conf_CancerBiopsy["EX_drib_LPAREN_e_RPAREN_"]=3
conf_CancerBiopsy["EX_etoh_LPAREN_e_RPAREN_"]=3
conf_CancerBiopsy["EX_fald_LPAREN_e_RPAREN_"]=3
conf_CancerBiopsy["EX_fe2_LPAREN_e_RPAREN_"]=3
conf_CancerBiopsy["EX_fe2_LPAREN_e_RPAREN_"]=3
conf_CancerBiopsy["EX_fe3_LPAREN_e_RPAREN_"]=3
conf_CancerBiopsy["EX_fmn_LPAREN_e_RPAREN_"]=3
conf_CancerBiopsy["EX_fol_LPAREN_e_RPAREN_"]=3
conf_CancerBiopsy["EX_for_LPAREN_e_RPAREN_"]=3
conf_CancerBiopsy["EX_fru_LPAREN_e_RPAREN_"]=3
conf_CancerBiopsy["EX_fuc_L_LPAREN_e_RPAREN_"]=3
conf_CancerBiopsy["EX_fum_LPAREN_e_RPAREN_"]=3
conf_CancerBiopsy["EX_gal_LPAREN_e_RPAREN_"]=3
conf_CancerBiopsy["EX_gam_LPAREN_e_RPAREN_"]=3
conf_CancerBiopsy["EX_gchola_LPAREN_e_RPAREN_"]=3
conf_CancerBiopsy["EX_glc_LPAREN_e_RPAREN_"]=3
conf_CancerBiopsy["EX_glcur_LPAREN_e_RPAREN_"]=3
conf_CancerBiopsy["EX_gln_L_LPAREN_e_RPAREN_"]=3
conf_CancerBiopsy["EX_glu_L_LPAREN_e_RPAREN_"]=3
conf_CancerBiopsy["EX_gluala_LPAREN_e_RPAREN_"]=3
conf_CancerBiopsy["EX_gly_LPAREN_e_RPAREN_"]=3
conf_CancerBiopsy["EX_glyb_LPAREN_e_RPAREN_"]=3
conf_CancerBiopsy["EX_glyc_LPAREN_e_RPAREN_"]=3
conf_CancerBiopsy["EX_glyc3p_LPAREN_e_RPAREN_"]=3
conf_CancerBiopsy["EX_glygn2_LPAREN_e_RPAREN_"]=3
conf_CancerBiopsy["EX_gsn_LPAREN_e_RPAREN_"]=3
conf_CancerBiopsy["EX_gthox_LPAREN_e_RPAREN_"]=3
conf_CancerBiopsy["EX_gthrd_LPAREN_e_RPAREN_"]=3
conf_CancerBiopsy["EX_gua_LPAREN_e_RPAREN_"]=3
conf_CancerBiopsy["EX_h_LPAREN_e_RPAREN_"]=3
conf_CancerBiopsy["EX_h2o_LPAREN_e_RPAREN_"]=3
conf_CancerBiopsy["EX_ha_LPAREN_e_RPAREN_"]=3
conf_CancerBiopsy["EX_HC00250_LPAREN_e_RPAREN_"]=3
conf_CancerBiopsy["EX_HC01609_LPAREN_e_RPAREN_"]=3
conf_CancerBiopsy["EX_HC01610_LPAREN_e_RPAREN_"]=3
conf_CancerBiopsy["EX_hdca_LPAREN_e_RPAREN_"]=3
conf_CancerBiopsy["EX_hdcea_LPAREN_e_RPAREN_"]=3
conf_CancerBiopsy["EX_his_L_LPAREN_e_RPAREN_"]=3
conf_CancerBiopsy["EX_hxan_LPAREN_e_RPAREN_"]=3
conf_CancerBiopsy["EX_ile_L_LPAREN_e_RPAREN_"]=3
conf_CancerBiopsy["EX_inost_LPAREN_e_RPAREN_"]=3
conf_CancerBiopsy["EX_k_LPAREN_e_RPAREN_"]=3
conf_CancerBiopsy["EX_lac_L_LPAREN_e_RPAREN_"]=3
conf_CancerBiopsy["EX_lac_D_LPAREN_e_RPAREN_"]=-1
conf_CancerBiopsy["EX_lcts_LPAREN_e_RPAREN_"]=3
conf_CancerBiopsy["EX_leu_L_LPAREN_e_RPAREN_"]=3
conf_CancerBiopsy["EX_leuktrA4_LPAREN_e_RPAREN_"]=3
conf_CancerBiopsy["EX_leuktrD4_LPAREN_e_RPAREN_"]=3
conf_CancerBiopsy["EX_leuktrE4_LPAREN_e_RPAREN_"]=3
conf_CancerBiopsy["EX_lnlc_LPAREN_e_RPAREN_"]=3
conf_CancerBiopsy["EX_lnlnca_LPAREN_e_RPAREN_"]=3
conf_CancerBiopsy["EX_lnlncg_LPAREN_e_RPAREN_"]=3
conf_CancerBiopsy["EX_lpchol_hs_LPAREN_e_RPAREN_"]=3
conf_CancerBiopsy["EX_lys_L_LPAREN_e_RPAREN_"]=3
conf_CancerBiopsy["EX_mag_hs_LPAREN_e_RPAREN_"]=3
conf_CancerBiopsy["EX_mal_L_LPAREN_e_RPAREN_"]=3
conf_CancerBiopsy["EX_malt_LPAREN_e_RPAREN_"]=3
conf_CancerBiopsy["EX_man_LPAREN_e_RPAREN_"]=3
conf_CancerBiopsy["EX_meoh_LPAREN_e_RPAREN_"]=3
conf_CancerBiopsy["EX_met_L_LPAREN_e_RPAREN_"]=3
conf_CancerBiopsy["EX_n2m2nmasn_LPAREN_e_RPAREN_"]=3
conf_CancerBiopsy["EX_na1_LPAREN_e_RPAREN_"]=3
conf_CancerBiopsy["EX_nac_LPAREN_e_RPAREN_"]=3
conf_CancerBiopsy["EX_ncam_LPAREN_e_RPAREN_"]=3
conf_CancerBiopsy["EX_nh4_LPAREN_e_RPAREN_"]=3
conf_CancerBiopsy["EX_no2_LPAREN_e_RPAREN_"]=3
conf_CancerBiopsy["EX_o2_LPAREN_e_RPAREN_"]=3
conf_CancerBiopsy["EX_o2s_LPAREN_e_RPAREN_"]=3
conf_CancerBiopsy["EX_ocdca_LPAREN_e_RPAREN_"]=3
conf_CancerBiopsy["EX_ocdcea_LPAREN_e_RPAREN_"]=3
conf_CancerBiopsy["EX_octa_LPAREN_e_RPAREN_"]=3
conf_CancerBiopsy["EX_orn_LPAREN_e_RPAREN_"]=3
conf_CancerBiopsy["EX_oxa_LPAREN_e_RPAREN_"]=3
conf_CancerBiopsy["EX_pe_hs_LPAREN_e_RPAREN_"]=3
conf_CancerBiopsy["EX_pglyc_hs_LPAREN_e_RPAREN_"]=3
conf_CancerBiopsy["EX_phe_L_LPAREN_e_RPAREN_"]=3
conf_CancerBiopsy["EX_pheme_LPAREN_e_RPAREN_"]=3
conf_CancerBiopsy["EX_pi_LPAREN_e_RPAREN_"]=3
conf_CancerBiopsy["EX_pnto_R_LPAREN_e_RPAREN_"]=3
conf_CancerBiopsy["EX_ppa_LPAREN_e_RPAREN_"]=3
conf_CancerBiopsy["EX_pro_L_LPAREN_e_RPAREN_"]=3
conf_CancerBiopsy["EX_prostgh2_LPAREN_e_RPAREN_"]=3
conf_CancerBiopsy["EX_ps_hs_LPAREN_e_RPAREN_"]=3
conf_CancerBiopsy["EX_ptrc_LPAREN_e_RPAREN_"]=3
conf_CancerBiopsy["EX_ptvstlac_LPAREN_e_RPAREN_"]=3
conf_CancerBiopsy["EX_pydam_LPAREN_e_RPAREN_"]=3
conf_CancerBiopsy["EX_pydx_LPAREN_e_RPAREN_"]=3
conf_CancerBiopsy["EX_pydx5p_LPAREN_e_RPAREN_"]=3
conf_CancerBiopsy["EX_pydxn_LPAREN_e_RPAREN_"]=3
conf_CancerBiopsy["EX_pyr_LPAREN_e_RPAREN_"]=3
conf_CancerBiopsy["EX_q10h2_LPAREN_e_RPAREN_"]=3
conf_CancerBiopsy["EX_retfa_LPAREN_e_RPAREN_"]=3
conf_CancerBiopsy["EX_retinol_LPAREN_e_RPAREN_"]=3
conf_CancerBiopsy["EX_retn_LPAREN_e_RPAREN_"]=3
conf_CancerBiopsy["EX_rib_D_LPAREN_e_RPAREN_"]=3
conf_CancerBiopsy["EX_ribflv_LPAREN_e_RPAREN_"]=3
conf_CancerBiopsy["EX_sbt_DASH_d_LPAREN_e_RPAREN_"]=3
conf_CancerBiopsy["EX_sel_LPAREN_e_RPAREN_"]=3
conf_CancerBiopsy["EX_ser_L_LPAREN_e_RPAREN_"]=3
conf_CancerBiopsy["EX_so4_LPAREN_e_RPAREN_"]=3
conf_CancerBiopsy["EX_sph1p_LPAREN_e_RPAREN_"]=3
conf_CancerBiopsy["EX_spmd_LPAREN_e_RPAREN_"]=3
conf_CancerBiopsy["EX_strch1_LPAREN_e_RPAREN_"]=3
conf_CancerBiopsy["EX_strch2_LPAREN_e_RPAREN_"]=3
conf_CancerBiopsy["EX_sucr_LPAREN_e_RPAREN_"]=3
conf_CancerBiopsy["EX_tag_hs_LPAREN_e_RPAREN_"]=3
conf_CancerBiopsy["EX_thm_LPAREN_e_RPAREN_"]=3
conf_CancerBiopsy["EX_thr_L_LPAREN_e_RPAREN_"]=3
conf_CancerBiopsy["EX_thymd_LPAREN_e_RPAREN_"]=3
conf_CancerBiopsy["EX_tre_LPAREN_e_RPAREN_"]=3
conf_CancerBiopsy["EX_trp_L_LPAREN_e_RPAREN_"]=3
conf_CancerBiopsy["EX_ttdca_LPAREN_e_RPAREN_"]=3
conf_CancerBiopsy["EX_tyr_L_LPAREN_e_RPAREN_"]=3
conf_CancerBiopsy["EX_ura_LPAREN_e_RPAREN_"]=3
conf_CancerBiopsy["EX_urea_LPAREN_e_RPAREN_"]=3
conf_CancerBiopsy["EX_uri_LPAREN_e_RPAREN_"]=3
conf_CancerBiopsy["EX_utp_LPAREN_e_RPAREN_"]=3
conf_CancerBiopsy["EX_val_L_LPAREN_e_RPAREN_"]=3
conf_CancerBiopsy["EX_xmp_LPAREN_e_RPAREN_"]=3
conf_CancerBiopsy["EX_xyl_D_LPAREN_e_RPAREN_"]=3
conf_CancerBiopsy["EX_xylt_LPAREN_e_RPAREN_"]=3


```


```python
conf_NormalBiopsy["DM_4hrpo"]=3
conf_NormalBiopsy["DM_datp_n_"]=3
conf_NormalBiopsy["DM_dctp_n_"]=3
conf_NormalBiopsy["DM_dgtp_n_"]=3
conf_NormalBiopsy["DM_dttp_n_"]=3
conf_NormalBiopsy["DM_Lcystin"]=3
conf_NormalBiopsy["DM_pe_hs_LPAREN_r_RPAREN_"]=3
conf_NormalBiopsy["EX_2hb_LPAREN_e_RPAREN_"]=3
conf_NormalBiopsy["EX_34hpp_"]=3
conf_NormalBiopsy["EX_3hpvs_LPAREN_e_RPAREN_"]=3
conf_NormalBiopsy["EX_3mob_LPAREN_e_RPAREN_"]=3
conf_NormalBiopsy["EX_4mop_LPAREN_e_RPAREN_"]=3
conf_NormalBiopsy["EX_ac_LPAREN_e_RPAREN_"]=3
conf_NormalBiopsy["EX_acac_LPAREN_e_RPAREN_"]=3
conf_NormalBiopsy["EX_acald_LPAREN_e_RPAREN_"]=3
conf_NormalBiopsy["EX_acetone_LPAREN_e_RPAREN_"]=3
conf_NormalBiopsy["EX_acgam_LPAREN_e_RPAREN_"]=3
conf_NormalBiopsy["EX_acmana_LPAREN_e_RPAREN_"]=3
conf_NormalBiopsy["EX_ade_LPAREN_e_RPAREN_"]=3
conf_NormalBiopsy["EX_adn_LPAREN_e_RPAREN_"]=3
conf_NormalBiopsy["EX_adpcbl_LPAREN_e_RPAREN_"]=3
conf_NormalBiopsy["EX_akg_LPAREN_e_RPAREN_"]=3
conf_NormalBiopsy["EX_ala_B_LPAREN_e_RPAREN_"]=3
conf_NormalBiopsy["EX_ala_D_LPAREN_e_RPAREN_"]=3
conf_NormalBiopsy["EX_ala_L_LPAREN_e_RPAREN_"]=3
conf_NormalBiopsy["EX_am9csa_LPAREN_e_RPAREN_"]=3
conf_NormalBiopsy["EX_amp_LPAREN_e_RPAREN_"]=3
conf_NormalBiopsy["EX_arab_L_LPAREN_e_RPAREN_"]=3
conf_NormalBiopsy["EX_arachd_LPAREN_e_RPAREN_"]=3
conf_NormalBiopsy["EX_arg_L_LPAREN_e_RPAREN_"]=3
conf_NormalBiopsy["EX_asn_L_LPAREN_e_RPAREN_"]=3
conf_NormalBiopsy["EX_asp_L_LPAREN_e_RPAREN_"]=3
conf_NormalBiopsy["EX_atp_LPAREN_e_RPAREN_"]=3
conf_NormalBiopsy["EX_bhb_LPAREN_e_RPAREN_"]=3
conf_NormalBiopsy["EX_btn_LPAREN_e_RPAREN_"]=3
conf_NormalBiopsy["EX_ca2_LPAREN_e_RPAREN_"]=3
conf_NormalBiopsy["EX_carn_LPAREN_e_RPAREN_"]=3
conf_NormalBiopsy["EX_caro_LPAREN_e_RPAREN_"]=3
conf_NormalBiopsy["EX_cgly_LPAREN_e_RPAREN_"]=3
conf_NormalBiopsy["EX_chol_LPAREN_e_RPAREN_"]=3
conf_NormalBiopsy["EX_chsterol_LPAREN_e_RPAREN_"]=3
conf_NormalBiopsy["EX_cit_LPAREN_e_RPAREN_"]=3
conf_NormalBiopsy["EX_cl_LPAREN_e_RPAREN_"]=3
conf_NormalBiopsy["EX_CLPND_LPAREN_e_RPAREN_"]=3
conf_NormalBiopsy["EX_co2_LPAREN_e_RPAREN_"]=3
conf_NormalBiopsy["EX_creat_LPAREN_e_RPAREN_"]=3
conf_NormalBiopsy["EX_crn_LPAREN_e_RPAREN_"]=3
conf_NormalBiopsy["EX_crvnc_LPAREN_e_RPAREN_"]=3
conf_NormalBiopsy["EX_csa_LPAREN_e_RPAREN_"]=3
conf_NormalBiopsy["EX_csn_LPAREN_e_RPAREN_"]=3
conf_NormalBiopsy["EX_cys_L_LPAREN_e_RPAREN_"]=3
conf_NormalBiopsy["EX_cytd_LPAREN_e_RPAREN_"]=3
conf_NormalBiopsy["EX_dad_2_LPAREN_e_RPAREN_"]=3
conf_NormalBiopsy["EX_dag_hs_LPAREN_e_RPAREN_"]=3
conf_NormalBiopsy["EX_dcmp_LPAREN_e_RPAREN_"]=3
conf_NormalBiopsy["EX_dcyt_LPAREN_e_RPAREN_"]=3
conf_NormalBiopsy["EX_ddca_LPAREN_e_RPAREN_"]=3
conf_NormalBiopsy["EX_dgsn_LPAREN_e_RPAREN_"]=3
conf_NormalBiopsy["EX_dhdascb_LPAREN_e_RPAREN_"]=3
conf_NormalBiopsy["EX_din_LPAREN_e_RPAREN_"]=3
conf_NormalBiopsy["EX_dopa_LPAREN_e_RPAREN_"]=3
conf_NormalBiopsy["EX_drib_LPAREN_e_RPAREN_"]=3
conf_NormalBiopsy["EX_etoh_LPAREN_e_RPAREN_"]=3
conf_NormalBiopsy["EX_fald_LPAREN_e_RPAREN_"]=3
conf_NormalBiopsy["EX_fe2_LPAREN_e_RPAREN_"]=3
conf_NormalBiopsy["EX_fe2_LPAREN_e_RPAREN_"]=3
conf_NormalBiopsy["EX_fe3_LPAREN_e_RPAREN_"]=3
conf_NormalBiopsy["EX_fmn_LPAREN_e_RPAREN_"]=3
conf_NormalBiopsy["EX_fol_LPAREN_e_RPAREN_"]=3
conf_NormalBiopsy["EX_for_LPAREN_e_RPAREN_"]=3
conf_NormalBiopsy["EX_fru_LPAREN_e_RPAREN_"]=3
conf_NormalBiopsy["EX_fuc_L_LPAREN_e_RPAREN_"]=3
conf_NormalBiopsy["EX_fum_LPAREN_e_RPAREN_"]=3
conf_NormalBiopsy["EX_gal_LPAREN_e_RPAREN_"]=3
conf_NormalBiopsy["EX_gam_LPAREN_e_RPAREN_"]=3
conf_NormalBiopsy["EX_gchola_LPAREN_e_RPAREN_"]=3
conf_NormalBiopsy["EX_glc_LPAREN_e_RPAREN_"]=3
conf_NormalBiopsy["EX_glcur_LPAREN_e_RPAREN_"]=3
conf_NormalBiopsy["EX_gln_L_LPAREN_e_RPAREN_"]=3
conf_NormalBiopsy["EX_glu_L_LPAREN_e_RPAREN_"]=3
conf_NormalBiopsy["EX_gluala_LPAREN_e_RPAREN_"]=3
conf_NormalBiopsy["EX_gly_LPAREN_e_RPAREN_"]=3
conf_NormalBiopsy["EX_glyb_LPAREN_e_RPAREN_"]=3
conf_NormalBiopsy["EX_glyc_LPAREN_e_RPAREN_"]=3
conf_NormalBiopsy["EX_glyc3p_LPAREN_e_RPAREN_"]=3
conf_NormalBiopsy["EX_glygn2_LPAREN_e_RPAREN_"]=3
conf_NormalBiopsy["EX_gsn_LPAREN_e_RPAREN_"]=3
conf_NormalBiopsy["EX_gthox_LPAREN_e_RPAREN_"]=3
conf_NormalBiopsy["EX_gthrd_LPAREN_e_RPAREN_"]=3
conf_NormalBiopsy["EX_gua_LPAREN_e_RPAREN_"]=3
conf_NormalBiopsy["EX_h_LPAREN_e_RPAREN_"]=3
conf_NormalBiopsy["EX_h2o_LPAREN_e_RPAREN_"]=3
conf_NormalBiopsy["EX_ha_LPAREN_e_RPAREN_"]=3
conf_NormalBiopsy["EX_HC00250_LPAREN_e_RPAREN_"]=3
conf_NormalBiopsy["EX_HC01609_LPAREN_e_RPAREN_"]=3
conf_NormalBiopsy["EX_HC01610_LPAREN_e_RPAREN_"]=3
conf_NormalBiopsy["EX_hdca_LPAREN_e_RPAREN_"]=3
conf_NormalBiopsy["EX_hdcea_LPAREN_e_RPAREN_"]=3
conf_NormalBiopsy["EX_his_L_LPAREN_e_RPAREN_"]=3
conf_NormalBiopsy["EX_hxan_LPAREN_e_RPAREN_"]=3
conf_NormalBiopsy["EX_ile_L_LPAREN_e_RPAREN_"]=3
conf_NormalBiopsy["EX_inost_LPAREN_e_RPAREN_"]=3
conf_NormalBiopsy["EX_k_LPAREN_e_RPAREN_"]=3
conf_NormalBiopsy["EX_lac_L_LPAREN_e_RPAREN_"]=3
conf_NormalBiopsy["EX_lac_D_LPAREN_e_RPAREN_"]=3
conf_NormalBiopsy["EX_lcts_LPAREN_e_RPAREN_"]=3
conf_NormalBiopsy["EX_leu_L_LPAREN_e_RPAREN_"]=3
conf_NormalBiopsy["EX_leuktrA4_LPAREN_e_RPAREN_"]=3
conf_NormalBiopsy["EX_leuktrD4_LPAREN_e_RPAREN_"]=3
conf_NormalBiopsy["EX_leuktrE4_LPAREN_e_RPAREN_"]=3
conf_NormalBiopsy["EX_lnlc_LPAREN_e_RPAREN_"]=3
conf_NormalBiopsy["EX_lnlnca_LPAREN_e_RPAREN_"]=3
conf_NormalBiopsy["EX_lnlncg_LPAREN_e_RPAREN_"]=3
conf_NormalBiopsy["EX_lpchol_hs_LPAREN_e_RPAREN_"]=3
conf_NormalBiopsy["EX_lys_L_LPAREN_e_RPAREN_"]=3
conf_NormalBiopsy["EX_mag_hs_LPAREN_e_RPAREN_"]=3
conf_NormalBiopsy["EX_mal_L_LPAREN_e_RPAREN_"]=3
conf_NormalBiopsy["EX_malt_LPAREN_e_RPAREN_"]=3
conf_NormalBiopsy["EX_man_LPAREN_e_RPAREN_"]=3
conf_NormalBiopsy["EX_meoh_LPAREN_e_RPAREN_"]=3
conf_NormalBiopsy["EX_met_L_LPAREN_e_RPAREN_"]=3
conf_NormalBiopsy["EX_n2m2nmasn_LPAREN_e_RPAREN_"]=3
conf_NormalBiopsy["EX_na1_LPAREN_e_RPAREN_"]=3
conf_NormalBiopsy["EX_nac_LPAREN_e_RPAREN_"]=3
conf_NormalBiopsy["EX_ncam_LPAREN_e_RPAREN_"]=3
conf_NormalBiopsy["EX_nh4_LPAREN_e_RPAREN_"]=3
conf_NormalBiopsy["EX_no2_LPAREN_e_RPAREN_"]=3
conf_NormalBiopsy["EX_o2_LPAREN_e_RPAREN_"]=3
conf_NormalBiopsy["EX_o2s_LPAREN_e_RPAREN_"]=3
conf_NormalBiopsy["EX_ocdca_LPAREN_e_RPAREN_"]=3
conf_NormalBiopsy["EX_ocdcea_LPAREN_e_RPAREN_"]=3
conf_NormalBiopsy["EX_octa_LPAREN_e_RPAREN_"]=3
conf_NormalBiopsy["EX_orn_LPAREN_e_RPAREN_"]=3
conf_NormalBiopsy["EX_oxa_LPAREN_e_RPAREN_"]=3
conf_NormalBiopsy["EX_pe_hs_LPAREN_e_RPAREN_"]=3
conf_NormalBiopsy["EX_pglyc_hs_LPAREN_e_RPAREN_"]=3
conf_NormalBiopsy["EX_phe_L_LPAREN_e_RPAREN_"]=3
conf_NormalBiopsy["EX_pheme_LPAREN_e_RPAREN_"]=3
conf_NormalBiopsy["EX_pi_LPAREN_e_RPAREN_"]=3
conf_NormalBiopsy["EX_pnto_R_LPAREN_e_RPAREN_"]=3
conf_NormalBiopsy["EX_ppa_LPAREN_e_RPAREN_"]=3
conf_NormalBiopsy["EX_pro_L_LPAREN_e_RPAREN_"]=3
conf_NormalBiopsy["EX_prostgh2_LPAREN_e_RPAREN_"]=3
conf_NormalBiopsy["EX_ps_hs_LPAREN_e_RPAREN_"]=3
conf_NormalBiopsy["EX_ptrc_LPAREN_e_RPAREN_"]=3
conf_NormalBiopsy["EX_ptvstlac_LPAREN_e_RPAREN_"]=3
conf_NormalBiopsy["EX_pydam_LPAREN_e_RPAREN_"]=3
conf_NormalBiopsy["EX_pydx_LPAREN_e_RPAREN_"]=3
conf_NormalBiopsy["EX_pydx5p_LPAREN_e_RPAREN_"]=3
conf_NormalBiopsy["EX_pydxn_LPAREN_e_RPAREN_"]=3
conf_NormalBiopsy["EX_pyr_LPAREN_e_RPAREN_"]=3
conf_NormalBiopsy["EX_q10h2_LPAREN_e_RPAREN_"]=3
conf_NormalBiopsy["EX_retfa_LPAREN_e_RPAREN_"]=3
conf_NormalBiopsy["EX_retinol_LPAREN_e_RPAREN_"]=3
conf_NormalBiopsy["EX_retn_LPAREN_e_RPAREN_"]=3
conf_NormalBiopsy["EX_rib_D_LPAREN_e_RPAREN_"]=3
conf_NormalBiopsy["EX_ribflv_LPAREN_e_RPAREN_"]=3
conf_NormalBiopsy["EX_sbt_DASH_d_LPAREN_e_RPAREN_"]=3
conf_NormalBiopsy["EX_sel_LPAREN_e_RPAREN_"]=3
conf_NormalBiopsy["EX_ser_L_LPAREN_e_RPAREN_"]=3
conf_NormalBiopsy["EX_so4_LPAREN_e_RPAREN_"]=3
conf_NormalBiopsy["EX_sph1p_LPAREN_e_RPAREN_"]=3
conf_NormalBiopsy["EX_spmd_LPAREN_e_RPAREN_"]=3
conf_NormalBiopsy["EX_strch1_LPAREN_e_RPAREN_"]=3
conf_NormalBiopsy["EX_strch2_LPAREN_e_RPAREN_"]=3
conf_NormalBiopsy["EX_sucr_LPAREN_e_RPAREN_"]=3
conf_NormalBiopsy["EX_tag_hs_LPAREN_e_RPAREN_"]=3
conf_NormalBiopsy["EX_thm_LPAREN_e_RPAREN_"]=3
conf_NormalBiopsy["EX_thr_L_LPAREN_e_RPAREN_"]=3
conf_NormalBiopsy["EX_thymd_LPAREN_e_RPAREN_"]=3
conf_NormalBiopsy["EX_tre_LPAREN_e_RPAREN_"]=3
conf_NormalBiopsy["EX_trp_L_LPAREN_e_RPAREN_"]=3
conf_NormalBiopsy["EX_ttdca_LPAREN_e_RPAREN_"]=3
conf_NormalBiopsy["EX_tyr_L_LPAREN_e_RPAREN_"]=3
conf_NormalBiopsy["EX_ura_LPAREN_e_RPAREN_"]=3
conf_NormalBiopsy["EX_urea_LPAREN_e_RPAREN_"]=3
conf_NormalBiopsy["EX_uri_LPAREN_e_RPAREN_"]=3
conf_NormalBiopsy["EX_utp_LPAREN_e_RPAREN_"]=3
conf_NormalBiopsy["EX_val_L_LPAREN_e_RPAREN_"]=3
conf_NormalBiopsy["EX_xmp_LPAREN_e_RPAREN_"]=3
conf_NormalBiopsy["EX_xyl_D_LPAREN_e_RPAREN_"]=3
conf_NormalBiopsy["EX_xylt_LPAREN_e_RPAREN_"]=3
```


```python

conf_CancerBiopsy["biomass_reaction"]=3
conf_CancerBiopsy["biomass_DNA"]=3
conf_CancerBiopsy["biomass_RNA"]=3
conf_CancerBiopsy["biomass_carbohydrate"]=3
conf_CancerBiopsy["biomass_lipid"]=3
conf_CancerBiopsy["biomass_other"]=3
conf_CancerBiopsy["biomass_protein"]=3
conf_CancerBiopsy["DM_atp_c_"]=3


conf_NormalBiopsy["biomass_reaction"]=3
conf_NormalBiopsy["biomass_DNA"]=3
conf_NormalBiopsy["biomass_RNA"]=3
conf_NormalBiopsy["biomass_carbohydrate"]=3
conf_NormalBiopsy["biomass_lipid"]=3
conf_NormalBiopsy["biomass_other"]=3
conf_NormalBiopsy["biomass_protein"]=3
conf_NormalBiopsy["DM_atp_c_"]=3
```


```python
metas = ['0.014 biomass_DNA_c + 0.058 biomass_RNA_c + 0.071 biomass_carbohydrate_c + 0.097 biomass_lipid_c + 0.054 biomass_other_c + 0.706 biomass_protein_c --> ', '3pg_c', '4abut_c', '4hpro_LT_c', 'accoa_m', 'accoa_m --> coa_m', 'ade_c', 'adn_c', 'adp_c', 'akg_m', 'ala_B_c', 'ala_L_c', 'amet_c', 'amp_c', 'arg_L_c', 'asn_L_c', 'asp_D_c', 'asp_L_c', 'atp_c', 'bhb_c', 'cdp_c', 'CE1936_c', 'chol_c', 'chsterol_c', 'cit_c', 'citr_L_c', 'cmp_c', 'creat_c', 'crm_hs_c', 'crtn_c', 'ctp_c', 'cys_L_c', 'dag_hs_c', 'dhap_c', 'e4p_c', 'f6p_c', 'fdp_c', 'fru_c', 'fum_c', 'g1p_c', 'g3p_c', 'g6p_c', 'gdp_c', 'glc_D_c', 'glc_D_e', 'glc_D_e --> glc_D_c', 'gln_L_m', 'gln_L_c', 'gln_L_e --> gln_L_c', 'glu_L_c', 'glyb_c', 'gly_c', 'gmp_c', 'gthox_c', 'gthrd_c', 'gua_c', 'HC00342_c', 'his_L_c', 'hxan_c', 'icit_c', 'ile_L_c', 'lac_L_c', 'leu_L_c', 'leu_L_c',  'lys_L_c', 'mag_hs_c', 'mal_L_c', 'met_L_c', 'nad_c', 'nadh_c', 'nadh_m', 'nad_m', 'nadp_c', 'nadph_c', 'nadph_m', 'nadph_m', 'nadp_m', 'oaa_m', 'orn_c', 'pa_hs_c', 'pe_hs_c', 'pep_c', 'phe_L_c', 'pmtcoa_c --> coa_c', 'pro_D_c', 'pro_L_c', 'ps_hs_c', 'ptrc_c', 'pyr_c', 'pyr_m', 'r5p_c', 'ru5p_D_c', 's7p_c', 'ser_L_c', 'spmd_c', 'succ_c', 'succoa_m --> coa_m', 'tag_hs_c', 'thr_L_c', 'trp_L_c', 'tym_c', 'tyr_L_c', 'udp_c', 'ump_c', 'utp_c', 'val_L_c']
```


```python
%%time

from corda import CORDA


opt_NormalBiopsy = CORDA(model=Recon2, confidence=conf_NormalBiopsy, n=5, met_prod=metas,  penalty_factor=1000) 
opt_NormalBiopsy.build()
print(opt_NormalBiopsy)
model_NormalBiopsy=opt_NormalBiopsy.cobra_model(name="NormalBiopsy")
print(model_NormalBiopsy.optimize())
```

    build status: reconstruction complete
    Inc. reactions: 2351/7891
     - unclear: 182/640
     - exclude: 250/2368
     - low and medium: 721/3507
     - high: 1198/1376
    
    <Solution 0.022 at 0x7f80e0007438>
    CPU times: user 22min 4s, sys: 644 ms, total: 22min 4s
    Wall time: 22min 4s



```python
%%time

from corda import CORDA


opt_CancerBiopsy = CORDA(model=Recon2, confidence=conf_CancerBiopsy, n=5, met_prod=metas,  penalty_factor=1000) 
opt_CancerBiopsy.build()
print(opt_CancerBiopsy)
model_CancerBiopsy=opt_CancerBiopsy.cobra_model(name="CancerBiopsy")
print(model_CancerBiopsy.optimize())
```

    build status: reconstruction complete
    Inc. reactions: 2720/7891
     - unclear: 283/1381
     - exclude: 155/1276
     - low and medium: 734/3499
     - high: 1548/1735
    
    <Solution 0.033 at 0x7f80dfe476d8>
    CPU times: user 23min 32s, sys: 644 ms, total: 23min 32s
    Wall time: 23min 32s



```python
cp=model_NormalBiopsy.copy()
cp.optimize()
cp.summary(fva=True)
```

    IN FLUXES                                        OUT FLUXES                                      OBJECTIVES
    -----------------------------------------------  ----------------------------------------------  ----------------------
    id                   Flux  Range                 id                   Flux  Range                biomass_reac...  0.022
    ---------------  --------  --------------------  ---------------  --------  -------------------
    glc_D_e          0.205     [0.00605, 1]          h_e              0.431     [-0.977, 5.4]
    o2_e             0.0294    [-0.0149, 1]          lac_L_e          0.359     [-0.09, 2.81]
    lys_L_e          0.013     [0.013, 0.013]        pyr_e            0.0503    [-0.012, 2.75]
    atp_e            0.0124    [-0.0325, 0.8]        pi_e             0.0202    [0, 1.76]
    gthrd_e          0.012     [-0.023, 0.09]        glyc_S_e         0.0152    [0, 1.33]
    leu_L_e          0.012     [0.012, 0.012]        nh4_e            0.0147    [-0.012, 1.3]
    ps_hs_e          0.012     [0.00391, 0.012]      leuktrD4_e       0.012     [0, 0.012]
    leuktrA4_e       0.012     [0, 0.012]            xmp_e            0.0104    [-0.012, 0.832]
    o2s_e            0.012     [0, 0.012]            dag_hs_e         0.00785   [-0.000385, 0.0235]
    ser_L_e          0.012     [0.012, -0.387]       akg_e            0.00776   [-0.012, 0.371]
    man_e            0.012     [0.012, -0.691]       ala_L_e          0.00411   [-0.04, 0.359]
    gly_e            0.0119    [0.03, -0.369]        Rtotal_e         0.00325   [0, 0.395]
    gln_L_e          0.0114    [0.09, -0.288]        dgtp_m           0.00323   [0, 0.856]
    pro_L_e          0.00907   [0.00907, 0.00907]    co2_e            0.00314   [0, 1.79]
    arg_L_e          0.0079    [0.0079, 0.0079]      lac_D_e          0.001     [0, 0.001]
    val_L_e          0.00776   [0.00776, 0.00776]    for_e            0.000449  [-0.012, 0.04]
    fum_e            0.00776   [0.012, -0.407]       c8crn_e          0.000449  [0, 0.012]
    thr_L_e          0.00688   [0.00688, 0.00688]    ac_e             0         [-0.012, 2.03]
    ile_L_e          0.00629   [0.00629, 0.012]      dgtp_n           0         [-0.012, 0.844]
    asn_L_e          0.00615   [0.00615, 0.00615]    datp_m           0         [0, 0.833]
    phe_L_e          0.00571   [0.00571, 0.012]      amp_e            0         [-0.012, 0.833]
    datp_n           0.00374   [0.012, -0.821]       abt_e            0         [0, 0.774]
    lpchol_hs_e      0.00363   [0, 0.0119]           xylt_e           0         [-0.012, 0.762]
    tyr_L_e          0.00351   [0.00351, 0.012]      oxa_e            0         [0, 0.623]
    met_L_e          0.00337   [0.00337, 0.012]      aicar_e          0         [0, 0.614]
    his_L_e          0.00278   [0, 0.012]            ppi_e            0         [0, 0.586]
    acac_e           0.00269   [-0.00931, 0.012]     fe3_e            0         [-0.012, 0.545]
    pe_hs_e          0.00137   [0.012, -0.0189]      elaid_e          0         [0, 0.428]
    biomass_other_c  0.00119   [0.00119, 0.00119]    hdcea_e          0         [0, 0.428]
    uri_e            0.00118   [0.012, -0.0538]      ocdcea_e         0         [0, 0.428]
    cys_L_e          0.00102   [0.00102, 0.012]      vacc_e           0         [0, 0.428]
    cytd_e           0.000859  [0.012, -0.0538]      mal_L_e          0         [-0.012, 0.407]
    inost_e          0.000513  [0.000513, 0.000513]  5oxpro_e         0         [0, 0.398]
    arachd_e         0.000449  [0, 0.012]            cit_e            0         [-0.012, 0.393]
    crn_e            0.000449  [0, 0.012]            gal_e            0         [0, 0.363]
    sph1p_e          0.000385  [0.000385, 0.012]     ade_e            0         [0, 0.356]
    pglyc_hs_e       0.000321  [0.000321, 0.012]     prpp_e           0         [0, 0.356]
    trp_L_e          0.000293  [0.000293, 0.000293]  nrvnc_e          0         [0, 0.272]
    dttp_n           0.000288  [0.000288, -0.0117]   cgly_e           0         [-0.012, 0.101]
    dctp_n           0.000208  [0.012, -0.0538]      glu_L_e          0         [-0.012, 0.101]
    h2o_e            0         [-3.27, 3.59]         clpnd_e          0         [0, 0.096]
    fe2_e            0         [-0.012, 0.545]       tmndnc_e         0         [0, 0.096]
    hdca_e           0         [-0.0116, 0.3]        eicostet_e       0         [0, 0.084]
    ocdca_e          0         [0, 0.07]             dlnlcg_e         0         [0, 0.072]
    lnlc_e           0         [-0.012, 0.06]        lnlncg_e         0         [0, 0.072]
    3hpvs_e          0         [0, 0.012]            lnlnca_e         0         [-0.012, 0.072]
    3mob_e           0         [0, 0.012]            gthox_e          0         [0, 0.0565]
    4hpro_LT_m       0         [0, 0.012]            utp_e            0         [-0.02, 0.0458]
    4mop_e           0         [0, 0.012]            pheme_e          0         [0, 0.0429]
    arab_L_e         0         [0, 0.012]            bilglcur_e       0         [0, 0.041]
    asp_L_e          0         [0, 0.012]            co_e             0         [0, 0.041]
    carn_e           0         [0, 0.012]            bildglcur_e      0         [0, 0.04]
    crvnc_e          0         [0, 0.012]            fald_e           0         [-0.012, 0.024]
    csa_e            0         [0, 0.012]            meoh_e           0         [-0.012, 0.024]
    dopa_e           0         [0, 0.012]            pe_hs_r          0         [-0.012, 0.0189]
    etoh_e           0         [0, 0.012]            cholate_e        0         [0, 0.0156]
    fol_e            0         [0, 0.012]            3mlda_e          0         [0, 0.0126]
    gam_e            0         [0, 0.012]            34dhoxpeg_e      0         [0, 0.012]
    gluala_e         0         [0, 0.012]            3bcrn_e          0         [0, 0.012]
    ha_e             0         [0, 0.012]            3ddcrn_e         0         [0, 0.012]
    mag_hs_e         0         [0, 0.012]            3deccrn_e        0         [0, 0.012]
    n2m2nmasn_e      0         [0, 0.012]            3hdececrn_e      0         [0, 0.012]
    octa_e           0         [0, 0.012]            3hexdcrn_e       0         [0, 0.012]
    ptvstlac_e       0         [0, 0.012]            3hpvstet_e       0         [0, 0.012]
    retinol_e        0         [0, 0.012]            3ivcrn_e         0         [0, 0.012]
    thymd_e          0         [0, 0.012]            3octdec2crn_e    0         [0, 0.012]
    gchola_e         0         [-0.00355, 0.012]     3octdeccrn_e     0         [0, 0.012]
    ppa_e            0         [-0.00926, 0.012]     3octdece1crn_e   0         [0, 0.012]
    bhb_e            0         [-0.00931, 0.012]     3tdcrn_e         0         [0, 0.012]
    glyc3p_e         0         [-0.0119, 0.012]      3tetd7ecoacrn_e  0         [0, 0.012]
    HC01609_e        0         [-0.012, 0.012]       3thexddcoacrn_e  0         [0, 0.012]
    HC01610_e        0         [-0.012, 0.012]       3ttetddcoacrn_e  0         [0, 0.012]
    pnto_R_e         0         [0, 0.011]            adrn_e           0         [0, 0.012]
    dcmp_e           0         [0.012, -0.0538]      am19cs_e         0         [0, 0.012]
                                                     am1csa_e         0         [0, 0.012]
                                                     am9csa_e         0         [0, 0.012]
                                                     c6crn_e          0         [0, 0.012]
                                                     fol_c            0         [0, 0.012]
                                                     ivcrn_e          0         [0, 0.012]
                                                     leuktrB4_e       0         [0, 0.012]
                                                     n5m2masn_g       0         [0, 0.012]
                                                     ptvstm3_e        0         [0, 0.012]
                                                     retn_e           0         [0, 0.012]
                                                     1glyc_hs_e       0         [0, 0.0117]
                                                     ethamp_r         0         [0, 0.0116]
                                                     fuc14galacgl...  0         [0, 0.0116]
                                                     octdececoa_c     0         [0, 0.011]
                                                     so4_e            0         [0, 0.011]
                                                     taur_c           0         [0, 0.011]
                                                     taur_e           0         [0, 0.011]
                                                     ahcys_e          0         [0, 0.00863]
                                                     4hphac_e         0         [0, 0.00849]
                                                     pchol_hs_e       0         [0, 0.00809]
                                                     pheacgln_e       0         [0, 0.00629]
                                                     HC00250_e        0         [0, 0.00549]
                                                     4mptnl_e         0         [0, 0.00355]
                                                     5adtststerone_e  0         [0, 0.00355]
                                                     andrstrn_e       0         [0, 0.00355]
                                                     andrstrnglc_e    0         [0, 0.00355]
                                                     aprgstrn_e       0         [0, 0.00355]
                                                     chsterol_e       0         [0, 0.00355]
                                                     glyc_e           0         [0.012, -0.21]



```python
cp=model_CancerBiopsy.copy()
cp.optimize()
cp.summary(fva=True)
```

    IN FLUXES                                        OUT FLUXES                                   OBJECTIVES
    -----------------------------------------------  -------------------------------------------  -----------------------
    id                   Flux  Range                 id                  Flux  Range              biomass_reac...  0.0334
    ---------------  --------  --------------------  ---------------  -------  -----------------
    glc_D_e          0.357     [0.197, 1]            h_e              0.793    [-0.411, 3.29]
    o2_e             0.2       [-0.0058, 1]          lac_L_e          0.456    [-0.09, 2.7]
    gln_L_e          0.0886    [0.09, -0.153]        pyr_e            0.174    [-0.012, 1.44]
    lys_L_e          0.0198    [0.0198, 0.03]        nh4_e            0.102    [-0.012, 0.475]
    Lcystin_c        0.0175    [-0.00522, 0.024]     cit_e            0.0771   [-0.012, 0.348]
    arg_L_e          0.012     [0.012, 0.012]        ala_L_e          0.0716   [-0.0169, 0.294]
    asn_L_e          0.012     [0.00933, 0.012]      h2o_e            0.0349   [-0.371, 1.69]
    4mop_e           0.012     [0.00622, 0.012]      co2_e            0.0242   [0, 1.73]
    4hpro_LT_m       0.012     [0, 0.012]            ac_e             0.0237   [-0.012, 0.504]
    asp_L_e          0.012     [0, 0.012]            urate_e          0.0228   [0, 0.0847]
    dgsn_e           0.012     [0, 0.012]            so4_e            0.0221   [0, 0.0584]
    o2s_e            0.012     [0, 0.012]            pi_e             0.0142   [-0.04, 0.119]
    acac_e           0.012     [-0.012, 0.012]       bhb_e            0.012    [-0.012, 0.012]
    gsn_e            0.012     [0.012, -0.0228]      glyc_e           0.0107   [-0.21, 0.299]
    dctp_n           0.012     [0.012, -0.0409]      uri_e            0.00708  [-0.012, 0.0634]
    ps_hs_e          0.012     [0.012, -0.147]       dag_hs_e         0.00665  [-0.012, 0.183]
    ser_L_e          0.012     [0.012, -0.287]       cys_L_e          0.00622  [-0.012, 0.0464]
    akg_e            0.012     [0.012, -0.331]       tchola_e         0.00511  [0, 0.0546]
    val_L_e          0.0118    [-0.000222, 0.03]     thymd_e          0.00268  [-0.012, 0.0994]
    thr_L_e          0.0104    [0.0104, 0.012]       3deccrn_e        0.00204  [0, 0.012]
    ile_L_e          0.00956   [0.00956, 0.012]      lac_D_e          0.001    [-0.001, 0.001]
    phe_L_e          0.00867   [0.00867, 0.012]      fe3_e            0        [-0.012, 0.545]
    3mob_e           0.0074    [0.012, -0.0182]      oxa_e            0        [0, 0.452]
    leu_L_e          0.00622   [0.00622, 0.012]      vacc_e           0        [0, 0.452]
    chol_e           0.00552   [-0.0119, 0.012]      elaid_e          0        [0, 0.449]
    tyr_L_e          0.00533   [0.00533, 0.012]      no_e             0        [0, 0.429]
    met_L_e          0.00511   [0.00511, 0.012]      gal_e            0        [0, 0.401]
    gchola_e         0.00511   [0.012, -0.0465]      ocdcea_e         0        [0, 0.349]
    ppa_e            0.00489   [0.012, -0.0625]      xylt_e           0        [-0.012, 0.328]
    his_L_e          0.00422   [0.00422, 0.012]      hdcea_e          0        [0, 0.328]
    pe_hs_e          0.00302   [0.012, -0.147]       glyc_S_e         0        [0, 0.299]
    arachd_e         0.00204   [0, 0.012]            acmana_e         0        [0, 0.294]
    crn_e            0.00204   [0, 0.012]            pro_L_e          0        [-0.012, 0.287]
    for_e            0.002     [0.012, -0.0999]      ha_pre1_e        0        [0, 0.193]
    biomass_other_c  0.0018    [0.0018, 0.0018]      nrvnc_e          0        [0, 0.175]
    dad_2_e          0.00179   [0, 0.012]            ethamp_r         0        [0, 0.159]
    cytd_e           0.00117   [0.012, -0.0634]      pe_hs_r          0        [-0.012, 0.147]
    gly_e            0.000893  [0.03, -0.269]        ha_e             0        [-0.012, 0.0898]
    inost_e          0.000779  [0.000779, 0.000779]  clpnd_e          0        [0, 0.06]
    hdca_e           0.000584  [-0.0294, 0.3]        dlnlcg_e         0        [0, 0.06]
    sph1p_e          0.000584  [0.012, -0.147]       eicostet_e       0        [0, 0.06]
    pglyc_hs_e       0.000487  [0.000487, 0.000487]  lnlnca_e         0        [0, 0.06]
    trp_L_e          0.000444  [0.000444, 0.012]     lnlncg_e         0        [0, 0.06]
    datp_n           0.00044   [0.00044, 0.00044]    taur_c           0        [0, 0.0584]
    dttp_n           0.000437  [0.000437, -0.0524]   taur_e           0        [0, 0.0584]
    dgtp_n           0.000331  [0.000331, 0.000331]  chsterol_e       0        [0, 0.0524]
    fe2_e            0         [-0.012, 0.545]       4mptnl_e         0        [0, 0.0511]
    ocdca_e          0         [0, 0.07]             aprgstrn_e       0        [0, 0.0511]
    lnlc_e           0         [0, 0.06]             5adtststerone_e  0        [0, 0.0495]
    utp_e            0         [0, 0.02]             andrstrn_e       0        [0, 0.0495]
    3hpvs_e          0         [0, 0.012]            fuc13galacgl...  0        [0, 0.048]
    amp_e            0         [0, 0.012]            fuc14galacgl...  0        [0, 0.048]
    arab_L_e         0         [0, 0.012]            galfucgalacg...  0        [0, 0.0456]
    crvnc_e          0         [0, 0.012]            adn_e            0        [-0.00868, 0.045]
    dcyt_e           0         [0, 0.012]            5adtststeron...  0        [0, 0.0441]
    din_e            0         [0, 0.012]            andrstrnglc_e    0        [0, 0.0441]
    dopa_e           0         [0, 0.012]            tststeroneglc_e  0        [0, 0.044]
    etoh_e           0         [0, 0.012]            man_e            0        [-0.012, 0.036]
    fmn_e            0         [0, 0.012]            fuc_L_e          0        [0, 0.0347]
    gam_e            0         [0, 0.012]            13_cis_retng...  0        [0, 0.024]
    gluala_e         0         [0, 0.012]            gthox_e          0        [0, 0.024]
    leuktrA4_e       0         [0, 0.012]            retnglc_e        0        [0, 0.024]
    n2m2nmasn_e      0         [0, 0.012]            fald_e           0        [-0.012, 0.024]
    orn_e            0         [0, 0.012]            meoh_e           0        [-0.012, 0.024]
    ptrc_e           0         [0, 0.012]            glyb_e           0        [0, 0.0238]
    ptvstlac_e       0         [0, 0.012]            1mncam_e         0        [0, 0.0236]
    retinol_e        0         [0, 0.012]            Rtotal_e         0        [0, 0.0235]
    spmd_e           0         [0, 0.012]            lpchol_hs_e      0        [0, 0.0235]
    2hb_e            0         [-0.00156, 0.012]     pchol_hs_e       0        [0, 0.0235]
    HC01609_e        0         [-0.012, 0.012]       ppi_e            0        [0, 0.02]
    HC01610_e        0         [-0.012, 0.012]       ump_e            0        [0, 0.02]
                                                     fucfucfucgal...  0        [0, 0.016]
                                                     34dhoxpeg_e      0        [0, 0.012]
                                                     3bcrn_e          0        [0, 0.012]
                                                     3ddcrn_e         0        [0, 0.012]
                                                     3hdececrn_e      0        [0, 0.012]
                                                     3hexdcrn_e       0        [0, 0.012]
                                                     3hpvstet_e       0        [0, 0.012]
                                                     3ivcrn_e         0        [0, 0.012]
                                                     3octdec2crn_e    0        [0, 0.012]
                                                     3octdeccrn_e     0        [0, 0.012]
                                                     3octdece1crn_e   0        [0, 0.012]
                                                     3tdcrn_e         0        [0, 0.012]
                                                     3tetd7ecoacrn_e  0        [0, 0.012]
                                                     3thexddcoacrn_e  0        [0, 0.012]
                                                     3ttetddcoacrn_e  0        [0, 0.012]
                                                     4abutn_e         0        [0, 0.012]
                                                     Asn_X_Ser_Thr_l  0        [0, 0.012]
                                                     CE1940_e         0        [0, 0.012]
                                                     abt_e            0        [0, 0.012]
                                                     ade_e            0        [0, 0.012]
                                                     adrn_e           0        [0, 0.012]
                                                     c6crn_e          0        [0, 0.012]
                                                     leuktrB4_e       0        [0, 0.012]
                                                     leuktrD4_e       0        [0, 0.012]
                                                     n5m2masn_g       0        [0, 0.012]
                                                     prpp_e           0        [0, 0.012]
                                                     ptvst_e          0        [0, 0.012]
                                                     ribflv_e         0        [0, 0.012]
                                                     retn_e           0        [-0.012, 0.012]
                                                     anth_c           0        [0, 0.0116]
                                                     c5dc_e           0        [0, 0.0102]
                                                     3mlda_e          0        [0, 0.00778]
                                                     5mta_e           0        [0, 0.00689]
                                                     sprm_c           0        [0, 0.00689]
                                                     sprm_e           0        [0, 0.00689]
                                                     4hphac_e         0        [0, 0.00667]
                                                     q10h2_e          0        [0, 0.00667]
                                                     ivcrn_e          0        [0, 0.00578]
                                                     bilglcur_e       0        [0, 0.00528]
                                                     co_e             0        [0, 0.00528]
                                                     pheme_e          0        [0, 0.00528]
                                                     pheacgln_e       0        [0, 0.00333]
                                                     nac_e            0        [0.0116, -0.012]
                                                     gthrd_e          0        [0.0584, -0.06]



```python
model_CancerBiopsy.optimize()
model_CancerBiopsy.metabolites.atp_c.summary()

print("ATP production")
print(0.681+0.681+0.0178+0.0117)

print("ATP consumption (not in biomass)")
print(0.324+0.321+0.012+0.012+0.0107)
```

    PRODUCING REACTIONS -- ATP(4-) (atp_c)
    --------------------------------------
    %      FLUX  RXN ID      REACTION
    ---  ------  ----------  --------------------------------------------------
    50%  0.686   PGK         3pg_c + atp_c <=> 13dpg_c + adp_c
    50%  0.686   PYK         adp_c + h_c + pep_c --> atp_c + pyr_c
    1%   0.0115  CYTK2       atp_c + dcmp_c <=> adp_c + dcdp_c
    
    CONSUMING REACTIONS -- ATP(4-) (atp_c)
    --------------------------------------
    %      FLUX  RXN ID      REACTION
    ---  ------  ----------  --------------------------------------------------
    50%  0.69    biomass...  0.716189801699717 ala_L_c + 0.508866855524079 a...
    24%  0.325   PFK         atp_c + f6p_c --> adp_c + fdp_c + h_c
    23%  0.322   HEX1        atp_c + glc_D_c --> adp_c + g6p_c + h_c
    1%   0.0142  NICRNS      atp_c + nicrns_c --> adp_c + h_c + nicrnt_c
    1%   0.012   HEX10       atp_c + gam_c --> adp_c + gam6p_c + h_c
    1%   0.0107  GLYK        atp_c + glyc_c --> adp_c + glyc3p_c + h_c
    ATP production
    1.3915000000000002
    ATP consumption (not in biomass)
    0.6797000000000001



```python
model_NormalBiopsy.optimize()
model_NormalBiopsy.metabolites.atp_c.summary()
model_NormalBiopsy.metabolites.atp_m.summary()

print("ATP production")
print(0.387+0.385+0.0763+0.0858)

print("ATP consumption (not in biomass)")
print(0.186+0.162+0.012+0.012+0.0763+0.0104)
```

    PRODUCING REACTIONS -- ATP(4-) (atp_c)
    --------------------------------------
    %      FLUX  RXN ID      REACTION
    ---  ------  ----------  --------------------------------------------------
    50%   0.416  PGK         3pg_c + atp_c <=> 13dpg_c + adp_c
    50%   0.416  PYK         adp_c + h_c + pep_c --> atp_c + pyr_c
    
    CONSUMING REACTIONS -- ATP(4-) (atp_c)
    --------------------------------------
    %      FLUX  RXN ID      REACTION
    ---  ------  ----------  --------------------------------------------------
    55%   0.454  biomass...  0.716189801699717 ala_L_c + 0.508866855524079 a...
    25%   0.206  PFK         atp_c + f6p_c --> adp_c + fdp_c + h_c
    17%   0.138  HEX7        atp_c + fru_c --> adp_c + f6p_c + h_c
    1%    0.012  HEX4        atp_c + man_c --> adp_c + h_c + man6p_c
    PRODUCING REACTIONS -- ATP(4-) (atp_m)
    --------------------------------------
    %       FLUX  RXN ID    REACTION
    ----  ------  --------  --------------------------------------------------
    100%  0.0579  ATPS4m    adp_m + 4.0 h_i + pi_m --> atp_m + h2o_m + 3.0 h_m
    
    CONSUMING REACTIONS -- ATP(4-) (atp_m)
    --------------------------------------
    %       FLUX  RXN ID    REACTION
    ----  ------  --------  --------------------------------------------------
    112%  0.0651  ADK1m     amp_m + atp_m <=> 2.0 adp_m
    ATP production
    0.9341
    ATP consumption (not in biomass)
    0.45870000000000005



```python

# Production of ATP from glucose in anaerobic conditions for Cancer Biopsies
# test for max ATP hydrolysis flux from only glucose

closedModel=model_CancerBiopsy.copy()
closedModel.objective="DM_atp_c_"

for rx in closedModel.exchanges:
    rx.upper_bound= 10
    rx.lower_bound= -0.01
    


######################################################################
## Glucose aerobic
######################################################################

closedModel.reactions.EX_h2o_LPAREN_e_RPAREN_.bounds=(-10,10)
closedModel.reactions.EX_h_LPAREN_e_RPAREN_.bounds=(-10,10)
closedModel.reactions.EX_glc_LPAREN_e_RPAREN_.bounds=(-1,-1)
closedModel.reactions.EX_gln_L_LPAREN_e_RPAREN_.bounds=(0,0)
closedModel.reactions.EX_o2_LPAREN_e_RPAREN_.bounds=(-1,-1)
closedModel.reactions.EX_co2_LPAREN_e_RPAREN_.bounds=(0,10)
closedModel.reactions.EX_pi_LPAREN_e_RPAREN_.bounds=(-10,10)


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


######################################################################
## Glucose anaerobic
######################################################################

closedModel.reactions.EX_h2o_LPAREN_e_RPAREN_.bounds=(-10,10)
closedModel.reactions.EX_h_LPAREN_e_RPAREN_.bounds=(-10,10)
closedModel.reactions.EX_glc_LPAREN_e_RPAREN_.bounds=(-1,-1)
closedModel.reactions.EX_gln_L_LPAREN_e_RPAREN_.bounds=(0,0)
closedModel.reactions.EX_o2_LPAREN_e_RPAREN_.bounds=(0,0)
closedModel.reactions.EX_co2_LPAREN_e_RPAREN_.bounds=(0,10)
closedModel.reactions.EX_pi_LPAREN_e_RPAREN_.bounds=(-10,10)


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

    ===========================
    Glucose aerobic
    ===========================
    Oxigen use o2_e <--  -1.0
    H2O h2o_e <=>  0.624498876
    CO2 production co2_e -->  0.629498876
    Glucose consumption glc_D_e <--  -1.0
    Glutamine consumption gln_L_e -->  0.0
    ATP production 1.95814078667
    IN FLUXES                                        OUT FLUXES                                   OBJECTIVES
    -----------------------------------------------  -------------------------------------------  ---------------
    id                   Flux  Range                 id               Flux  Range                 DM_atp_c_  1.96
    ---------------  --------  --------------------  -----------  --------  --------------------
    glc_D_e          1         [1, 1]                h_e          2.11      [1.63, 2.36]
    o2_e             1         [1, 1]                lac_D_e      1.7       [-0.01, 1.9]
    3mob_e           0.01      [0.01, 0.01]          co2_e        0.629     [0.498, 1.05]
    dctp_n           0.01      [0.01, 0.01]          h2o_e        0.624     [0.468, 1]
    dgsn_e           0.01      [0.01, 0.01]          pyr_e        0.297     [-0.01, 0.461]
    din_e            0.01      [0.01, 0.01]          nh4_e        0.194     [-0.00548, 0.223]
    gam_e            0.01      [0.01, 0.01]          ac_e         0.0592    [0.0297, 0.0592]
    glyc_e           0.01      [0.01, 0.01]          pi_e         0.0512    [0.0317, 0.101]
    gsn_e            0.01      [0.01, 0.01]          urate_e      0.0297    [0.0297, 0.0297]
    ps_hs_e          0.01      [0.01, 0.01]          so4_e        0.0296    [0, 0.0296]
    ser_L_e          0.01      [0.01, 0.01]          dag_hs_e     0.028     [0.00849, 0.028]
    val_L_e          0.01      [0.01, 0.01]          oxa_e        0.0257    [0, 0.0257]
    xylt_e           0.01      [0.01, 0.01]          ala_L_e      0.0227    [0.0128, 0.123]
    ile_L_e          0.01      [0.00229, 0.01]       for_e        0.0198    [-5e-05, 0.0298]
    asn_L_e          0.01      [0.00224, 0.01]       uri_e        0.019     [0.00903, 0.039]
    tyr_L_e          0.01      [0.00128, 0.01]       ppa_e        0.01      [-0.01, 0.0252]
    trp_L_e          0.01      [0.000106, 0.01]      34dhoxpeg_e  0.01      [0, 0.01]
    cytd_e           0.01      [5.4e-05, 0.01]       leuktrB4_e   0.01      [0, 0.01]
    4hpro_LT_m       0.01      [0, 0.01]             retn_e       0.01      [0, 0.01]
    crn_e            0.01      [0, 0.01]             nac_e        0.00989   [0, 0.00989]
    dopa_e           0.01      [0, 0.01]             4hphac_e     0.00872   [0, 0.00872]
    etoh_e           0.01      [0, 0.01]             3ivcrn_e     0.00634   [0, 0.01]
    gthrd_e          0.01      [0, 0.01]             CE1940_e     0.005     [0, 0.005]
    leuktrA4_e       0.01      [0, 0.01]             gthox_e      0.005     [0, 0.005]
    orn_e            0.01      [0, 0.01]             c5dc_e       0.00366   [0, 0.00526]
    retinol_e        0.01      [0, 0.01]             ptvst_e      0.000898  [0, 0.01]
    Lcystin_c        0.01      [-0.00481, 0.01]      thymd_e      0.000213  [0.000213, 0.000213]
    2hb_e            0.01      [-0.0075, 0.01]       lac_L_e      0         [-0.01, 1.9]
    pe_hs_e          0.01      [-0.0095, 0.01]       cit_e        0         [-0.01, 0.104]
    pe_hs_r          0.01      [-0.0095, 0.01]       ppi_e        0         [-0.01, 0.02]
    meoh_e           0.01      [0.01, -0.01]         3bcrn_e      0         [0, 0.01]
    gly_e            0.01      [0.01, -0.0157]       abt_e        0         [0, 0.01]
    cys_L_e          0.01      [0.01, -0.0196]       ivcrn_e      0         [0, 0.01]
    pro_L_e          0.01      [0.01, -0.0906]       leuktrD4_e   0         [0, 0.01]
    lys_L_e          0.0084    [0.00474, 0.01]       lpchol_hs_e  0         [0, 0.01]
    spmd_e           0.005     [0, 0.005]            ribflv_e     0         [0, 0.01]
    leu_L_e          0.00436   [-0.00564, 0.01]      vacc_e       0         [0, 0.01]
    arachd_e         0.00366   [0, 0.00549]          HC01609_e    0         [-0.01, 0.01]
    arg_L_e          0.00287   [0.00287, 0.00287]    HC01610_e    0         [-0.01, 0.01]
    asp_L_e          0.00282   [0, 0.01]             acac_e       0         [-0.01, 0.01]
    thr_L_e          0.0025    [0.0025, 0.01]        ade_e        0         [-0.01, 0.01]
    phe_L_e          0.00208   [0.00208, 0.00208]    amp_e        0         [-0.01, 0.01]
    akg_e            0.00203   [0.01, -0.104]        bhb_e        0         [-0.01, 0.01]
    chol_e           0.00132   [0.00132, 0.01]       fald_e       0         [-0.01, 0.01]
    met_L_e          0.00122   [0.00122, 0.00122]    fe2_e        0         [-0.01, 0.01]
    his_L_e          0.00101   [0.00101, 0.00101]    fe3_e        0         [-0.01, 0.01]
    ptvstlac_e       0.000898  [0, 0.01]             prpp_e       0         [-0.01, 0.01]
    biomass_other_c  0.000432  [0.000432, 0.000432]  sprm_c       0         [-0.01, 0.01]
    adn_e            0.00043   [0.00043, 0.00043]    sprm_e       0         [-0.01, 0.01]
    inost_e          0.000187  [0.000187, 0.000187]  taur_e       0         [-0.01, 0.01]
    sph1p_e          0.00014   [0.00014, 0.00014]    ump_e        0         [-0.01, 0.01]
    Rtotal_e         0.00014   [0.00014, -0.01]      anth_c       0         [0, 0.00989]
    pglyc_hs_e       0.000117  [0.000117, 0.000117]  glyb_e       0         [0, 0.00868]
    datp_n           0.000105  [0.000105, 0.000105]  4abutn_e     0         [0, 0.005]
    dttp_n           0.000105  [0.000105, 0.000105]  3deccrn_e    0         [0, 0.00366]
    dgtp_n           7.9e-05   [7.9e-05, 7.9e-05]    c6crn_e      0         [0, 0.00366]
    arab_L_e         0         [0, 0.01]
    fmn_e            0         [0, 0.01]
    o2s_e            0         [0, 0.01]
    ocdca_e          0         [0, 0.01]
    pchol_hs_e       0         [0, 0.01]
    utp_e            0         [0, 0.01]
    4mop_e           0         [-0.00564, 0.01]
    hdca_e           0         [0, 0.00563]
    ptrc_e           0         [0, 0.005]
    taur_c           0         [0.01, -0.01]
    ===========================
    Glucose anaerobic
    ===========================
    Oxigen use o2_e -->  0.0
    H2O h2o_e <=>  -0.00965014711111
    CO2 production co2_e -->  0.0
    Glucose consumption glc_D_e <--  -1.0
    Glutamine consumption gln_L_e -->  0.0
    ATP production 1.86504778133
    IN FLUXES                                        OUT FLUXES                                 OBJECTIVES
    -----------------------------------------------  -----------------------------------------  ---------------
    id                   Flux  Range                 id              Flux  Range                DM_atp_c_  1.87
    ---------------  --------  --------------------  -----------  -------  -------------------
    glc_D_e          1         [1, 1]                h_e          2.02     [1.92, 2.08]
    dctp_n           0.01      [0.01, 0.01]          lac_D_e      1.99     [-0.01, 2.07]
    gam_e            0.01      [0.01, 0.01]          pyr_e        0.0332   [-0.01, 0.129]
    o2s_e            0.01      [0.01, 0.01]          pi_e         0.0317   [0.0317, 0.0919]
    ps_hs_e          0.01      [0.01, 0.01]          nh4_e        0.0173   [0.00721, 0.049]
    leu_L_e          0.01      [-0.00461, 0.01]      thymd_e      0.00992  [0.00593, 0.00992]
    for_e            0.00976   [0.00976, -0.0117]    dag_hs_e     0.00849  [0.00849, 0.0187]
    h2o_e            0.00965   [-0.0325, 0.0762]     4mop_e       0.00564  [0.00564, -0.01]
    glyc_e           0.00965   [-0.00754, 0.01]      ser_L_e      0.00277  [0.00681, -0.01]
    akg_e            0.00757   [-0.00851, 0.01]      3mob_e       0.00194  [0.00718, -0.00743]
    val_L_e          0.00476   [-0.00461, 0.01]      lac_L_e      0        [-0.01, 2.07]
    lys_L_e          0.00474   [0.00474, 0.00474]    co2_e        0        [0, 0.0825]
    gly_e            0.00431   [0.00431, -0.0157]    ppi_e        0        [-0.01, 0.02]
    arg_L_e          0.00287   [0.00287, 0.00287]    cys_L_e      0        [-0.01, 0.0196]
    pro_L_e          0.00257   [0.0033, -0.0289]     glyc_S_e     0        [0, 0.0168]
    thr_L_e          0.0025    [0.0025, 0.01]        ppa_e        0        [0, 0.0153]
    ile_L_e          0.00229   [0.00229, 0.00229]    ala_L_e      0        [-0.00404, 0.0128]
    asn_L_e          0.00224   [0.00224, 0.01]       cit_e        0        [-0.01, 0.0118]
    phe_L_e          0.00208   [0.00208, 0.00208]    oxa_e        0        [0, 0.0105]
    asp_L_e          0.00168   [0, 0.01]             3bcrn_e      0        [0, 0.01]
    chol_e           0.00132   [0.00132, 0.00132]    abt_e        0        [0, 0.01]
    tyr_L_e          0.00128   [0.00128, 0.00128]    fe2_e        0        [0, 0.01]
    met_L_e          0.00122   [0.00122, 0.00122]    ivcrn_e      0        [0, 0.01]
    his_L_e          0.00101   [0.00101, 0.00101]    leuktrB4_e   0        [0, 0.01]
    acac_e           0.000898  [-0.01, 0.01]         leuktrD4_e   0        [0, 0.01]
    pe_hs_e          0.000497  [-0.0095, 0.01]       lpchol_hs_e  0        [0, 0.01]
    biomass_other_c  0.000432  [0.000432, 0.000432]  ptvst_e      0        [0, 0.01]
    adn_e            0.00043   [0.00043, 0.00043]    retn_e       0        [0, 0.01]
    uri_e            0.000428  [0.000428, -0.032]    ribflv_e     0        [0, 0.01]
    cytd_e           0.000312  [0, 0.01]             HC01609_e    0        [-0.01, 0.01]
    gsn_e            0.000289  [0.000289, 0.00428]   ade_e        0        [-0.01, 0.01]
    inost_e          0.000187  [0.000187, 0.000187]  amp_e        0        [-0.01, 0.01]
    Lcystin_c        0.000186  [-0.00481, 0.01]      bhb_e        0        [-0.01, 0.01]
    sph1p_e          0.00014   [0.00014, 0.00014]    fald_e       0        [-0.01, 0.01]
    Rtotal_e         0.00014   [0.00014, -0.01]      sprm_e       0        [-0.01, 0.01]
    pglyc_hs_e       0.000117  [0.000117, 0.000117]  taur_c       0        [-0.01, 0.01]
    trp_L_e          0.000106  [0.000106, 0.000106]  taur_e       0        [-0.01, 0.01]
    datp_n           0.000105  [0.000105, 0.000105]  ump_e        0        [-0.01, 0.01]
    dttp_n           0.000105  [0.000105, 0.000105]  3deccrn_e    0        [0, 0.00049]
    dgtp_n           7.9e-05   [7.9e-05, 7.9e-05]
    4hpro_LT_m       0         [0, 0.01]
    arab_L_e         0         [0, 0.01]
    crn_e            0         [0, 0.01]
    etoh_e           0         [0, 0.01]
    fe3_e            0         [0, 0.01]
    fmn_e            0         [0, 0.01]
    gthrd_e          0         [0, 0.01]
    leuktrA4_e       0         [0, 0.01]
    pchol_hs_e       0         [0, 0.01]
    ptvstlac_e       0         [0, 0.01]
    retinol_e        0         [0, 0.01]
    utp_e            0         [0, 0.01]
    xylt_e           0         [0, 0.01]
    2hb_e            0         [-0.0075, 0.01]
    pe_hs_r          0         [-0.0095, 0.01]
    HC01610_e        0         [-0.01, 0.01]
    meoh_e           0         [-0.01, 0.01]
    prpp_e           0         [-0.01, 0.01]
    din_e            0         [0, 0.00528]
    arachd_e         0         [0, 0.00049]
    hdca_e           0         [0, 0.00014]
    urate_e          0         [0, -0.00528]
    sprm_c           0         [0.01, -0.01]
    ac_e             0         [0, -0.0168]



```python

# Production of ATP from glucose in anaerobic conditions for Cancer Biopsies
# test for max ATP hydrolysis flux from only glucose

closedModel=model_NormalBiopsy.copy()
closedModel.objective="DM_atp_c_"

for rx in closedModel.exchanges:
    rx.upper_bound= 10
    rx.lower_bound= -0.01
    


######################################################################
## Glucose aerobic
######################################################################

closedModel.reactions.EX_h2o_LPAREN_e_RPAREN_.bounds=(-10,10)
closedModel.reactions.EX_h_LPAREN_e_RPAREN_.bounds=(-10,10)
closedModel.reactions.EX_glc_LPAREN_e_RPAREN_.bounds=(-1,-1)
closedModel.reactions.EX_gln_L_LPAREN_e_RPAREN_.bounds=(0,0)
closedModel.reactions.EX_o2_LPAREN_e_RPAREN_.bounds=(-1,-1)
closedModel.reactions.EX_co2_LPAREN_e_RPAREN_.bounds=(0,10)
closedModel.reactions.EX_pi_LPAREN_e_RPAREN_.bounds=(-10,10)


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


######################################################################
## Glucose anaerobic
######################################################################

closedModel.reactions.EX_h2o_LPAREN_e_RPAREN_.bounds=(-10,10)
closedModel.reactions.EX_h_LPAREN_e_RPAREN_.bounds=(-10,10)
closedModel.reactions.EX_glc_LPAREN_e_RPAREN_.bounds=(-1,-1)
closedModel.reactions.EX_gln_L_LPAREN_e_RPAREN_.bounds=(0,0)
closedModel.reactions.EX_o2_LPAREN_e_RPAREN_.bounds=(0,0)
closedModel.reactions.EX_co2_LPAREN_e_RPAREN_.bounds=(0,10)
closedModel.reactions.EX_pi_LPAREN_e_RPAREN_.bounds=(-10,10)


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

    ===========================
    Glucose aerobic
    ===========================
    Oxigen use o2_e <--  -1.0
    H2O h2o_e <=>  1.90562345067
    CO2 production co2_e -->  0.008321576
    Glucose consumption glc_D_e <--  -1.0
    Glutamine consumption gln_L_e -->  0.0
    ATP production 3.91671747733
    IN FLUXES                                        OUT FLUXES                                 OBJECTIVES
    -----------------------------------------------  -----------------------------------------  ---------------
    id                   Flux  Range                 id              Flux  Range                DM_atp_c_  3.92
    ---------------  --------  --------------------  ----------  --------  -------------------
    glc_D_e          1         [1, 1]                h_e         2.1       [2.04, 2.41]
    o2_e             1         [1, 1]                pyr_e       2.02      [1.77, 2.09]
    akg_e            0.01      [0.01, 0.01]          h2o_e       1.91      [1.58, 1.96]
    gam_e            0.01      [0.01, 0.01]          lac_L_e     0.0743    [-0.01, 0.375]
    glyc3p_e         0.01      [0.01, 0.01]          nh4_e       0.0318    [0.0287, 0.0531]
    ha_e             0.01      [0.01, 0.01]          ac_e        0.0241    [0.02, 0.0693]
    man_e            0.01      [0.01, 0.01]          abt_e       0.02      [0, 0.04]
    o2s_e            0.01      [0.01, 0.01]          cit_e       0.0124    [-0.01, 0.0163]
    glu_L_e          0.00852   [-0.00148, 0.01]      co2_e       0.00832   [0, 0.0483]
    acac_e           0.0057    [0.00258, 0.01]       oxa_e       0.00472   [0, 0.049]
    lys_L_e          0.00474   [0.00474, 0.00474]    bhb_e       0.00472   [0.0016, 0.00902]
    4hpro_LT_m       0.00472   [0.0016, 0.00902]     5oxpro_e    0.00404   [0, 0.01]
    leu_L_e          0.00436   [0.00436, 0.00436]    dctp_n      0.00309   [-0.01, 0.015]
    gly_e            0.00431   [0.00431, -0.0199]    Rtotal_e    0.00118   [0.00118, 0.0397]
    gluala_e         0.00404   [0, 0.01]             tmndnc_e    0.000979  [0.000979, 0.01]
    etoh_e           0.00376   [0, 0.01]             for_e       0.000163  [0.000163, 0.0202]
    pro_L_e          0.0033    [0.0033, 0.0033]      lac_D_e     0         [0, 0.385]
    dcmp_e           0.00316   [0.01, -0.0196]       pi_e        0         [0, 0.0838]
    arg_L_e          0.00287   [0.00287, 0.00287]    ppi_e       0         [-0.01, 0.0379]
    val_L_e          0.00282   [0.00282, 0.00282]    mal_L_e     0         [-0.01, 0.0374]
    thr_L_e          0.0025    [0.0025, 0.0025]      ala_L_e     0         [-0.01, 0.0302]
    fum_e            0.00245   [0.01, -0.0374]       xylt_e      0         [-0.01, 0.03]
    ile_L_e          0.00229   [0.00229, 0.00229]    glyc_S_e    0         [0, 0.0242]
    asn_L_e          0.00224   [0.00224, 0.00224]    ade_e       0         [-0.01, 0.0196]
    phe_L_e          0.00208   [0.00208, 0.00208]    amp_e       0         [-0.01, 0.0196]
    ser_L_e          0.00168   [0.01, -0.0142]       prpp_e      0         [-0.01, 0.0196]
    ps_hs_e          0.00151   [0.00142, 0.01]       fe2_e       0         [0, 0.01]
    lpchol_hs_e      0.00132   [0.00132, 0.0099]     fol_c       0         [0, 0.01]
    tyr_L_e          0.00128   [0.00128, 0.00128]    leuktrB4_e  0         [0, 0.01]
    met_L_e          0.00122   [0.00122, 0.00122]    leuktrD4_e  0         [0, 0.01]
    his_L_e          0.00101   [0.00101, 0.00101]    retn_e      0         [0, 0.01]
    crvnc_e          0.000979  [0.000979, 0.01]      HC01610_e   0         [-0.00208, 0.01]
    pe_hs_e          0.000497  [-0.0095, 0.01]       cgly_e      0         [-0.01, 0.01]
    biomass_other_c  0.000432  [0.000432, 0.000432]  fald_e      0         [-0.01, 0.01]
    atp_e            0.00043   [0.00043, 0.01]       glyc_e      0         [-0.01, 0.01]
    uri_e            0.000428  [-0.00957, 0.01]      taur_c      0         [-0.01, 0.01]
    cys_L_e          0.000373  [0.000373, 0.000373]  taur_e      0         [-0.01, 0.01]
    cytd_e           0.000312  [0.000312, 0.01]      dgtp_m      0         [0, 0.00992]
    xmp_e            0.000289  [0.000289, -0.00928]  1glyc_hs_e  0         [0, 0.00988]
    inost_e          0.000187  [0.000187, 0.000187]  pchol_hs_e  0         [0, 0.00858]
    sph1p_e          0.00014   [0.00014, 0.00014]    dag_hs_e    0         [-8.5e-05, 0.00849]
    pglyc_hs_e       0.000117  [0.000117, 0.01]
    trp_L_e          0.000106  [0.000106, 0.000106]
    datp_n           0.000105  [0.000105, 0.01]
    dttp_n           0.000105  [0.000105, 0.000105]
    dgtp_n           7.9e-05   [7.9e-05, 0.01]
    arab_L_e         0         [0, 0.01]
    fe3_e            0         [0, 0.01]
    fol_e            0         [0, 0.01]
    leuktrA4_e       0         [0, 0.01]
    mag_hs_e         0         [0, 0.01]
    retinol_e        0         [0, 0.01]
    gthrd_e          0         [-0.00148, 0.01]
    HC01609_e        0         [-0.00208, 0.01]
    pe_hs_r          0         [-0.0095, 0.01]
    utp_e            0         [-0.00957, 0.01]
    meoh_e           0         [-0.01, 0.01]
    hdca_e           0         [0, 0.00014]
    datp_m           0         [0, -0.00989]
    ===========================
    Glucose anaerobic
    ===========================
    Oxigen use o2_e -->  0.0
    H2O h2o_e <=>  -0.0855259146667
    CO2 production co2_e -->  0.008321576
    Glucose consumption glc_D_e <--  -1.0
    Glutamine consumption gln_L_e -->  0.0
    ATP production 1.91671747733
    IN FLUXES                                        OUT FLUXES                                OBJECTIVES
    -----------------------------------------------  ----------------------------------------  ---------------
    id                   Flux  Range                 id              Flux  Range               DM_atp_c_  1.92
    ---------------  --------  --------------------  ----------  --------  ------------------
    glc_D_e          1         [1, 1]                h_e         2.1       [2.04, 2.27]
    h2o_e            0.0855    [0.0428, 0.324]       lac_L_e     2.08      [-0.01, 2.12]
    akg_e            0.01      [0.01, 0.01]          ac_e        0.0317    [0.02, 0.0693]
    gam_e            0.01      [0.01, 0.01]          nh4_e       0.0287    [0.0287, 0.0531]
    glyc3p_e         0.01      [0.01, 0.01]          abt_e       0.02      [0, 0.04]
    ha_e             0.01      [0.01, 0.01]          pyr_e       0.0125    [-0.01, 0.0938]
    man_e            0.01      [0.01, 0.01]          cit_e       0.0111    [-0.01, 0.0163]
    o2s_e            0.01      [0.01, 0.01]          co2_e       0.00832   [0, 0.0483]
    etoh_e           0.01      [0, 0.01]             5oxpro_e    0.00404   [0, 0.01]
    glu_L_e          0.00852   [-0.00148, 0.01]      utp_e       0.002     [0.00957, -0.01]
    lys_L_e          0.00474   [0.00474, 0.00474]    oxa_e       0.0016    [0, 0.0101]
    leu_L_e          0.00436   [0.00436, 0.00436]    bhb_e       0.0016    [0.0016, 0.00902]
    gly_e            0.00431   [0.00431, -0.00578]   Rtotal_e    0.00118   [0.00118, 0.0322]
    gluala_e         0.00404   [0, 0.01]             tmndnc_e    0.000979  [0.000979, 0.01]
    pro_L_e          0.0033    [0.0033, 0.0033]      for_e       0.000163  [0.000163, 0.0202]
    arg_L_e          0.00287   [0.00287, 0.00287]    pi_e        0         [0, 0.0838]
    val_L_e          0.00282   [0.00282, 0.00282]    ppi_e       0         [-0.01, 0.0379]
    acac_e           0.00258   [0.00258, 0.01]       fum_e       0         [-0.01, 0.0374]
    thr_L_e          0.0025    [0.0025, 0.0025]      ala_L_e     0         [-0.00909, 0.0302]
    uri_e            0.00243   [-0.00957, 0.01]      xylt_e      0         [-0.01, 0.03]
    ile_L_e          0.00229   [0.00229, 0.00229]    glyc_S_e    0         [0, 0.0242]
    asn_L_e          0.00224   [0.00224, 0.00224]    dcmp_e      0         [-0.01, 0.0196]
    phe_L_e          0.00208   [0.00208, 0.00208]    ade_e       0         [-0.01, 0.0196]
    ser_L_e          0.00177   [0.01, -0.0142]       amp_e       0         [-0.01, 0.0196]
    4hpro_LT_m       0.0016    [0.0016, 0.00902]     prpp_e      0         [-0.01, 0.0196]
    ps_hs_e          0.00142   [0.00142, 0.01]       dctp_n      0         [-0.01, 0.015]
    lpchol_hs_e      0.00132   [0.00132, 0.00239]    fe2_e       0         [0, 0.01]
    tyr_L_e          0.00128   [0.00128, 0.00128]    fol_c       0         [0, 0.01]
    met_L_e          0.00122   [0.00122, 0.00122]    leuktrB4_e  0         [0, 0.01]
    mal_L_e          0.0011    [0.01, -0.0374]       leuktrD4_e  0         [0, 0.01]
    his_L_e          0.00101   [0.00101, 0.00101]    HC01610_e   0         [-0.00208, 0.01]
    crvnc_e          0.000979  [0.000979, 0.01]      cgly_e      0         [-0.01, 0.01]
    pe_hs_e          0.000497  [-0.0095, 0.01]       fald_e      0         [-0.01, 0.01]
    biomass_other_c  0.000432  [0.000432, 0.000432]  glyc_e      0         [-0.01, 0.01]
    atp_e            0.00043   [0.00043, 0.01]       taur_c      0         [-0.01, 0.01]
    cytd_e           0.000388  [0.000312, 0.01]      taur_e      0         [-0.01, 0.01]
    cys_L_e          0.000373  [0.000373, 0.000373]  dgtp_m      0         [0, 0.00992]
    xmp_e            0.000289  [0.000289, -0.00928]  datp_m      0         [0, 0.00989]
    inost_e          0.000187  [0.000187, 0.000187]  1glyc_hs_e  0         [0, 0.00988]
    sph1p_e          0.00014   [0.00014, 0.00014]    pchol_hs_e  0         [0, 0.00107]
    pglyc_hs_e       0.000117  [0.000117, 0.01]      retinol_e   0         [0, -0.01]
    trp_L_e          0.000106  [0.000106, 0.000106]
    datp_n           0.000105  [0.000105, 0.01]
    dttp_n           0.000105  [0.000105, 0.000105]
    dag_hs_e         8.5e-05   [8.5e-05, -0.00849]
    dgtp_n           7.9e-05   [7.9e-05, 0.01]
    arab_L_e         0         [0, 0.01]
    fe3_e            0         [0, 0.01]
    fol_e            0         [0, 0.01]
    leuktrA4_e       0         [0, 0.01]
    mag_hs_e         0         [0, 0.01]
    gthrd_e          0         [-0.00148, 0.01]
    HC01609_e        0         [-0.00208, 0.01]
    pe_hs_r          0         [-0.0095, 0.01]
    meoh_e           0         [-0.01, 0.01]
    hdca_e           0         [0, 0.00014]
    retn_e           0         [0, -0.01]
    lac_D_e          0         [0, -2.13]



```python

# Production of ATP from glucose in anaerobic conditions for Cancer Biopsies
# test for max ATP hydrolysis flux from only glucose

closedModel=model_CancerBiopsy.copy()
closedModel.objective="DM_atp_c_"

for rx in closedModel.exchanges:
    rx.upper_bound= 10
    rx.lower_bound= -0.01
    

######################################################################
## Glutamine aerobic
######################################################################

for rx in closedModel.exchanges:
    rx.upper_bound= 10
    rx.lower_bound= -0.1
    
closedModel.reactions.EX_h2o_LPAREN_e_RPAREN_.bounds=(-10,10)
closedModel.reactions.EX_h_LPAREN_e_RPAREN_.bounds=(-10,10)
closedModel.reactions.EX_glc_LPAREN_e_RPAREN_.bounds=(0,0)
closedModel.reactions.EX_gln_L_LPAREN_e_RPAREN_.bounds=(-1,-1)
closedModel.reactions.EX_o2_LPAREN_e_RPAREN_.bounds=(-1,-1)
closedModel.reactions.EX_co2_LPAREN_e_RPAREN_.bounds=(0,10)
closedModel.reactions.EX_pi_LPAREN_e_RPAREN_.bounds=(-10,10)

FBA = cobra.flux_analysis.pfba(closedModel)
print("===========================")
print("Glutamine aerobic")
print("===========================")
print("Oxigen use",closedModel.reactions.EX_o2_LPAREN_e_RPAREN_.reaction,FBA["EX_o2_LPAREN_e_RPAREN_"])
print("H2O",closedModel.reactions.EX_h2o_LPAREN_e_RPAREN_.reaction,FBA["EX_h2o_LPAREN_e_RPAREN_"])
print("CO2 production",closedModel.reactions.EX_co2_LPAREN_e_RPAREN_.reaction,FBA["EX_co2_LPAREN_e_RPAREN_"])
print("Glucose consumption",closedModel.reactions.EX_glc_LPAREN_e_RPAREN_.reaction,FBA["EX_glc_LPAREN_e_RPAREN_"])
print("Glutamine consumption",closedModel.reactions.EX_gln_L_LPAREN_e_RPAREN_.reaction,FBA["EX_gln_L_LPAREN_e_RPAREN_"])
print("ATP production",FBA["DM_atp_c_"])

closedModel.summary(fva=True)

######################################################################
## Glutamine anaerobic
######################################################################

closedModel.reactions.EX_h2o_LPAREN_e_RPAREN_.bounds=(-10,10)
closedModel.reactions.EX_h_LPAREN_e_RPAREN_.bounds=(-10,10)
closedModel.reactions.EX_glc_LPAREN_e_RPAREN_.bounds=(0,0)
closedModel.reactions.EX_gln_L_LPAREN_e_RPAREN_.bounds=(-1,-1)
closedModel.reactions.EX_o2_LPAREN_e_RPAREN_.bounds=(0,0)
closedModel.reactions.EX_co2_LPAREN_e_RPAREN_.bounds=(0,10)
closedModel.reactions.EX_pi_LPAREN_e_RPAREN_.bounds=(-10,10)

FBA = cobra.flux_analysis.pfba(closedModel)

print("===========================")
print("Glutamine anaerobic")
print("===========================")
print("Oxigen use",closedModel.reactions.EX_o2_LPAREN_e_RPAREN_.reaction,FBA["EX_o2_LPAREN_e_RPAREN_"])
print("H2O",closedModel.reactions.EX_h2o_LPAREN_e_RPAREN_.reaction,FBA["EX_h2o_LPAREN_e_RPAREN_"])
print("CO2 production",closedModel.reactions.EX_co2_LPAREN_e_RPAREN_.reaction,FBA["EX_co2_LPAREN_e_RPAREN_"])
print("Glucose consumption",closedModel.reactions.EX_glc_LPAREN_e_RPAREN_.reaction,FBA["EX_glc_LPAREN_e_RPAREN_"])
print("Glutamine consumption",closedModel.reactions.EX_gln_L_LPAREN_e_RPAREN_.reaction,FBA["EX_gln_L_LPAREN_e_RPAREN_"])
print("ATP production",FBA["DM_atp_c_"])

closedModel.summary(fva=True)

```

    ===========================
    Glutamine aerobic
    ===========================
    Oxigen use o2_e <--  -1.0
    H2O h2o_e <=>  -0.1
    CO2 production co2_e -->  0.0
    Glucose consumption glc_D_e -->  0.0
    Glutamine consumption gln_L_e <--  -1.0
    ATP production 1.1342295279
    IN FLUXES                                        OUT FLUXES                                   OBJECTIVES
    -----------------------------------------------  -------------------------------------------  ---------------
    id                   Flux  Range                 id               Flux  Range                 DM_atp_c_  1.13
    ---------------  --------  --------------------  -----------  --------  --------------------
    gln_L_e          1         [1, 1]                nh4_e        1.2       [1.2, 1.55]
    o2_e             1         [1, 1]                h_e          0.97      [-0.988, 1.4]
    3mob_e           0.1       [0.1, 0.1]            cit_e        0.688     [-0.1, 0.951]
    bhb_e            0.1       [0.1, 0.1]            ala_L_e      0.654     [0.193, 0.768]
    chol_e           0.1       [0.1, 0.1]            pro_L_e      0.566     [0.435, 1.32]
    cys_L_e          0.1       [0.1, 0.1]            pi_e         0.315     [0.315, 0.815]
    din_e            0.1       [0.1, 0.1]            ac_e         0.174     [0.174, 0.274]
    fe2_e            0.1       [0.1, 0.1]            urate_e      0.174     [0.174, 0.174]
    gam_e            0.1       [0.1, 0.1]            acac_e       0.1       [0.1, 0.1]
    glyc_e           0.1       [0.1, 0.1]            fald_e       0.1       [0.1, 0.1]
    lac_D_e          0.1       [0.1, 0.1]            fe3_e        0.1       [0.1, 0.1]
    lac_L_e          0.1       [0.1, 0.1]            vacc_e       0.1       [0.1, 0.1]
    meoh_e           0.1       [0.1, 0.1]            taur_c       0.1       [-0.1, 0.1]
    o2s_e            0.1       [0.1, 0.1]            glyb_e       0.0987    [0.0987, 0.0987]
    ocdca_e          0.1       [0.1, 0.1]            dag_hs_e     0.0985    [0.0985, 0.0985]
    ppa_e            0.1       [0.1, 0.1]            gly_e        0.0957    [0.0957, 0.0957]
    ps_hs_e          0.1       [0.1, 0.1]            uri_e        0.0736    [0.0736, 0.274]
    ser_L_e          0.1       [0.1, 0.1]            Lcystin_c    0.0498    [0.0498, 0.0498]
    val_L_e          0.1       [0.1, 0.1]            akg_e        0.0315    [-0.1, 0.425]
    xylt_e           0.1       [0.1, 0.1]            pyr_e        0.0138    [-0.1, 0.475]
    taur_e           0.1       [0.1, -0.1]           retn_e       0.0104    [0, 0.1]
    h2o_e            0.1       [0.963, -1]           for_e        0.000163  [0.000163, 0.000163]
    gsn_e            0.0746    [0.0746, 0.0746]      co2_e        0         [0, 1.5]
    dctp_n           0.0744    [0.0744, 0.0744]      ppi_e        0         [-0.1, 0.2]
    retinol_e        0.0104    [0, 0.1]              3bcrn_e      0         [0, 0.1]
    lys_L_e          0.00474   [0.00474, 0.00474]    ivcrn_e      0         [0, 0.1]
    leu_L_e          0.00436   [-0.0956, 0.1]        lpchol_hs_e  0         [0, 0.1]
    arg_L_e          0.00287   [0.00287, 0.00287]    ribflv_e     0         [0, 0.1]
    thr_L_e          0.0025    [0.0025, 0.0025]      HC01609_e    0         [-0.1, 0.1]
    ile_L_e          0.00229   [0.00229, 0.1]        ade_e        0         [-0.1, 0.1]
    asn_L_e          0.00224   [0.00224, 0.00224]    amp_e        0         [-0.1, 0.1]
    phe_L_e          0.00208   [0.00208, 0.00208]    prpp_e       0         [-0.1, 0.1]
    tyr_L_e          0.00128   [0.00128, 0.00128]    sprm_c       0         [-0.1, 0.1]
    met_L_e          0.00122   [0.00122, 0.00122]    sprm_e       0         [-0.1, 0.1]
    his_L_e          0.00101   [0.00101, 0.00101]    ump_e        0         [-0.1, 0.1]
    pe_hs_e          0.000497  [-0.0995, 0.1]        crn_e        0         [0, -0.1]
    hdca_e           0.00049   [0.00049, 0.00063]
    biomass_other_c  0.000432  [0.000432, 0.000432]
    adn_e            0.00043   [0.00043, 0.00043]
    inost_e          0.000187  [0.000187, 0.000187]
    sph1p_e          0.00014   [0.00014, 0.00014]
    Rtotal_e         0.00014   [0.00014, -0.1]
    pglyc_hs_e       0.000117  [0.000117, 0.000117]
    trp_L_e          0.000106  [0.000106, 0.000106]
    datp_n           0.000105  [0.000105, 0.000105]
    dttp_n           0.000105  [0.000105, 0.000105]
    dgtp_n           7.9e-05   [7.9e-05, 7.9e-05]
    cytd_e           5.4e-05   [5.4e-05, 5.4e-05]
    2hb_e            0         [0, 0.1]
    asp_L_e          0         [0, 0.1]
    etoh_e           0         [0, 0.1]
    fmn_e            0         [0, 0.1]
    pchol_hs_e       0         [0, 0.1]
    utp_e            0         [0, 0.1]
    4mop_e           0         [-0.0956, 0.1]
    pe_hs_r          0         [-0.0995, 0.1]
    HC01610_e        0         [-0.1, 0.1]
    ===========================
    Glutamine anaerobic
    ===========================
    Oxigen use o2_e -->  0.0
    H2O h2o_e <=>  0.153133899429
    CO2 production co2_e -->  0.0
    Glucose consumption glc_D_e -->  0.0
    Glutamine consumption gln_L_e <--  -1.0
    ATP production 0.133201932762
    IN FLUXES                                        OUT FLUXES                                   OBJECTIVES
    -----------------------------------------------  -------------------------------------------  ----------------
    id                   Flux  Range                 id               Flux  Range                 DM_atp_c_  0.133
    ---------------  --------  --------------------  -----------  --------  --------------------
    gln_L_e          1         [1, 1]                nh4_e        1.1       [1.06, 1.53]
    bhb_e            0.1       [0.1, 0.1]            pro_L_e      0.61      [0.575, 0.918]
    cys_L_e          0.1       [0.1, 0.1]            ala_L_e      0.565     [0.193, 0.665]
    gam_e            0.1       [0.1, 0.1]            akg_e        0.457     [0.0783, 0.562]
    glyc_e           0.1       [0.1, 0.1]            h2o_e        0.153     [0.446, -0.59]
    lac_D_e          0.1       [0.1, 0.1]            acac_e       0.1       [0.1, 0.1]
    lac_L_e          0.1       [0.1, 0.1]            fald_e       0.1       [0.1, 0.1]
    meoh_e           0.1       [0.1, 0.1]            dag_hs_e     0.0985    [0.0985, 0.0985]
    o2s_e            0.1       [0.1, 0.1]            gly_e        0.0957    [0.0957, 0.0957]
    ps_hs_e          0.1       [0.1, 0.1]            cit_e        0.0689    [-0.1, 0.258]
    ser_L_e          0.1       [0.1, 0.1]            HC01609_e    0.0633    [-0.1, 0.1]
    xylt_e           0.1       [0.1, 0.1]            Lcystin_c    0.0498    [0.0498, 0.0498]
    asp_L_e          0.1       [0.0381, 0.1]         3bcrn_e      0.0381    [0.0381, 0.1]
    HC01610_e        0.0633    [0.1, -0.1]           dctp_n       0.0273    [0.0273, 0.0273]
    ppa_e            0.0474    [0.0474, 0.0474]      pi_e         0.00978   [0.00978, 0.51]
    crn_e            0.0381    [0.0381, 0.1]         for_e        0.000163  [0.000163, 0.000163]
    h_e              0.031     [-0.283, 1.22]        co2_e        0         [0, 0.946]
    cytd_e           0.0277    [0.0277, 0.0277]      pyr_e        0         [-0.1, 0.372]
    lys_L_e          0.00474   [0.00474, 0.00474]    ppi_e        0         [-0.1, 0.2]
    4mop_e           0.00436   [-0.0956, 0.1]        lpchol_hs_e  0         [0, 0.1]
    arg_L_e          0.00287   [0.00287, 0.00287]    ribflv_e     0         [0, 0.1]
    3mob_e           0.00282   [-0.0972, 0.1]        ade_e        0         [-0.1, 0.1]
    thr_L_e          0.0025    [0.0025, 0.0025]      amp_e        0         [-0.1, 0.1]
    ile_L_e          0.00229   [0.00229, 0.00229]    prpp_e       0         [-0.1, 0.1]
    asn_L_e          0.00224   [0.00224, 0.00224]    sprm_c       0         [-0.1, 0.1]
    phe_L_e          0.00208   [0.00208, 0.00208]    taur_c       0         [-0.1, 0.1]
    chol_e           0.00132   [0.00132, 0.00132]    taur_e       0         [-0.1, 0.1]
    tyr_L_e          0.00128   [0.00128, 0.00128]    ump_e        0         [-0.1, 0.1]
    met_L_e          0.00122   [0.00122, 0.00122]    ac_e         0         [0, 0.0619]
    his_L_e          0.00101   [0.00101, 0.00101]    ivcrn_e      0         [0, 0.0619]
    pe_hs_e          0.000497  [-0.0995, 0.1]        etoh_e       0         [0, -0.0619]
    biomass_other_c  0.000432  [0.000432, 0.000432]
    adn_e            0.00043   [0.00043, 0.00043]
    uri_e            0.000428  [0.000428, -0.2]
    gsn_e            0.000289  [0.000289, 0.000289]
    inost_e          0.000187  [0.000187, 0.000187]
    sph1p_e          0.00014   [0.00014, 0.00014]
    Rtotal_e         0.00014   [0.00014, -0.1]
    pglyc_hs_e       0.000117  [0.000117, 0.000117]
    trp_L_e          0.000106  [0.000106, 0.000106]
    datp_n           0.000105  [0.000105, 0.000105]
    dttp_n           0.000105  [0.000105, 0.000105]
    dgtp_n           7.9e-05   [7.9e-05, 7.9e-05]
    fmn_e            0         [0, 0.1]
    pchol_hs_e       0         [0, 0.1]
    utp_e            0         [0, 0.1]
    leu_L_e          0         [-0.0956, 0.1]
    val_L_e          0         [-0.0972, 0.1]
    pe_hs_r          0         [-0.0995, 0.1]
    sprm_e           0         [-0.1, 0.1]
    retinol_e        0         [0, 0.0619]
    hdca_e           0         [0, 0.00014]
    retn_e           0         [0, -0.0619]


# Production of ATP from glucose in anaerobic conditions for Normal Biopsies
# test for max ATP hydrolysis flux from only glucose

closedModel=model_NormalBiopsy.copy()
closedModel.objective="DM_atp_c_"

for rx in closedModel.exchanges:
    rx.upper_bound= 10
    rx.lower_bound= -0.01
    

######################################################################
## Glutamine aerobic
######################################################################

for rx in closedModel.exchanges:
    rx.upper_bound= 10
    rx.lower_bound= -0.1
    
closedModel.reactions.EX_h2o_LPAREN_e_RPAREN_.bounds=(-10,10)
closedModel.reactions.EX_h_LPAREN_e_RPAREN_.bounds=(-10,10)
closedModel.reactions.EX_glc_LPAREN_e_RPAREN_.bounds=(0,0)
closedModel.reactions.EX_gln_L_LPAREN_e_RPAREN_.bounds=(-1,-1)
closedModel.reactions.EX_o2_LPAREN_e_RPAREN_.bounds=(-1,-1)
closedModel.reactions.EX_co2_LPAREN_e_RPAREN_.bounds=(0,10)
closedModel.reactions.EX_pi_LPAREN_e_RPAREN_.bounds=(-10,10)

FBA = cobra.flux_analysis.pfba(closedModel)
print("===========================")
print("Glutamine aerobic")
print("===========================")
print("Oxigen use",closedModel.reactions.EX_o2_LPAREN_e_RPAREN_.reaction,FBA["EX_o2_LPAREN_e_RPAREN_"])
print("H2O",closedModel.reactions.EX_h2o_LPAREN_e_RPAREN_.reaction,FBA["EX_h2o_LPAREN_e_RPAREN_"])
print("CO2 production",closedModel.reactions.EX_co2_LPAREN_e_RPAREN_.reaction,FBA["EX_co2_LPAREN_e_RPAREN_"])
print("Glucose consumption",closedModel.reactions.EX_glc_LPAREN_e_RPAREN_.reaction,FBA["EX_glc_LPAREN_e_RPAREN_"])
print("Glutamine consumption",closedModel.reactions.EX_gln_L_LPAREN_e_RPAREN_.reaction,FBA["EX_gln_L_LPAREN_e_RPAREN_"])
print("ATP production",FBA["DM_atp_c_"])

closedModel.summary(fva=True)

######################################################################
## Glutamine anaerobic
######################################################################

closedModel.reactions.EX_h2o_LPAREN_e_RPAREN_.bounds=(-10,10)
closedModel.reactions.EX_h_LPAREN_e_RPAREN_.bounds=(-10,10)
closedModel.reactions.EX_glc_LPAREN_e_RPAREN_.bounds=(0,0)
closedModel.reactions.EX_gln_L_LPAREN_e_RPAREN_.bounds=(-1,-1)
closedModel.reactions.EX_o2_LPAREN_e_RPAREN_.bounds=(0,0)
closedModel.reactions.EX_co2_LPAREN_e_RPAREN_.bounds=(0,10)
closedModel.reactions.EX_pi_LPAREN_e_RPAREN_.bounds=(-10,10)

FBA = cobra.flux_analysis.pfba(closedModel)

print("===========================")
print("Glutamine anaerobic")
print("===========================")
print("Oxigen use",closedModel.reactions.EX_o2_LPAREN_e_RPAREN_.reaction,FBA["EX_o2_LPAREN_e_RPAREN_"])
print("H2O",closedModel.reactions.EX_h2o_LPAREN_e_RPAREN_.reaction,FBA["EX_h2o_LPAREN_e_RPAREN_"])
print("CO2 production",closedModel.reactions.EX_co2_LPAREN_e_RPAREN_.reaction,FBA["EX_co2_LPAREN_e_RPAREN_"])
print("Glucose consumption",closedModel.reactions.EX_glc_LPAREN_e_RPAREN_.reaction,FBA["EX_glc_LPAREN_e_RPAREN_"])
print("Glutamine consumption",closedModel.reactions.EX_gln_L_LPAREN_e_RPAREN_.reaction,FBA["EX_gln_L_LPAREN_e_RPAREN_"])
print("ATP production",FBA["DM_atp_c_"])

closedModel.summary(fva=True)


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

res1=pool.apply_async(model_create, [Recon2, conf_CancerBiopsy, 5, metas ,100] )
res2=pool.apply_async(model_create, [Recon2, conf_CancerBiopsy, 5, metas ,500] )
res3=pool.apply_async(model_create, [Recon2, conf_CancerBiopsy, 5, metas ,1000] )
res4=pool.apply_async(model_create, [Recon2, conf_CancerBiopsy, 5, metas ,2500] )
res5=pool.apply_async(model_create, [Recon2, conf_CancerBiopsy, 5, metas ,5000] )

res6=pool.apply_async(model_create, [Recon2, conf_NormalBiopsy, 5, metas ,100] )
res7=pool.apply_async(model_create, [Recon2, conf_NormalBiopsy, 5, metas ,500] )
res8=pool.apply_async(model_create, [Recon2, conf_NormalBiopsy, 5, metas ,1000] )
res9=pool.apply_async(model_create, [Recon2, conf_NormalBiopsy, 5, metas ,2500] )
res0=pool.apply_async(model_create, [Recon2, conf_NormalBiopsy, 5, metas ,5000] )

opt_CancerBiopsy1=res1.get() 
opt_CancerBiopsy2=res2.get() 
opt_CancerBiopsy3=res3.get() 
opt_CancerBiopsy4=res4.get() 
opt_CancerBiopsy5=res5.get() 


opt_NormalBiopsy1=res6.get() 
opt_NormalBiopsy2=res7.get() 
opt_NormalBiopsy3=res8.get() 
opt_NormalBiopsy4=res9.get() 
opt_NormalBiopsy5=res0.get() 
```

    CPU times: user 1min 2s, sys: 9.83 s, total: 1min 12s
    Wall time: 29min 32s



```python
model_CancerBiopsy1=opt_CancerBiopsy1.cobra_model(name="Cancer1")
print(model_CancerBiopsy1.optimize())
model_CancerBiopsy2=opt_CancerBiopsy2.cobra_model(name="Cancer2")
print(model_CancerBiopsy2.optimize())
model_CancerBiopsy3=opt_CancerBiopsy3.cobra_model(name="Cancer3")
print(model_CancerBiopsy3.optimize())
model_CancerBiopsy4=opt_CancerBiopsy4.cobra_model(name="Cancer4")
print(model_CancerBiopsy4.optimize())
model_CancerBiopsy5=opt_CancerBiopsy5.cobra_model(name="Cancer5")
print(model_CancerBiopsy5.optimize())
```

    <Solution 0.029 at 0x7f8066914c50>
    <Solution 0.033 at 0x7f8051d085c0>
    <Solution 0.033 at 0x7f8064580208>
    <Solution 0.033 at 0x7f8052a9bf28>
    <Solution 0.033 at 0x7f805d8331d0>



```python
model_NormalBiopsy1=opt_NormalBiopsy1.cobra_model(name="Norm1")
print(model_NormalBiopsy1.optimize())
model_NormalBiopsy2=opt_NormalBiopsy2.cobra_model(name="Norm2")
print(model_NormalBiopsy2.optimize())
model_NormalBiopsy3=opt_NormalBiopsy3.cobra_model(name="Norm3")
print(model_NormalBiopsy3.optimize())
model_NormalBiopsy4=opt_NormalBiopsy4.cobra_model(name="Norm4")
print(model_NormalBiopsy4.optimize())
model_NormalBiopsy5=opt_NormalBiopsy5.cobra_model(name="Norm5")
print(model_NormalBiopsy5.optimize())

```

    <Solution 0.022 at 0x7f8094c2eba8>
    <Solution 0.022 at 0x7f8078dece48>
    <Solution 0.022 at 0x7f8072053978>
    <Solution 0.022 at 0x7f8085ef2be0>
    <Solution 0.022 at 0x7f807bb0f390>



```python
model_CancerBiopsy1.summary(fva=True)
```

    IN FLUXES                                        OUT FLUXES                                   OBJECTIVES
    -----------------------------------------------  -------------------------------------------  -----------------------
    id                   Flux  Range                 id                   Flux  Range             biomass_reac...  0.0291
    ---------------  --------  --------------------  ---------------  --------  ----------------
    glc_D_e          0.289     [0.159, 1]            h_e              0.6       [-0.0997, 3.29]
    o2_e             0.0594    [-0.0065, 1]          lac_L_e          0.503     [-0.09, 2.71]
    gln_L_e          0.037     [0.09, -0.167]        pyr_e            0.049     [-0.012, 1.42]
    Lcystin_c        0.0236    [-0.00532, 0.024]     cys_L_e          0.0459    [-0.012, 0.0466]
    val_L_e          0.0229    [-0.00174, 0.03]      nh4_e            0.0396    [-0.012, 0.502]
    lys_L_e          0.0172    [0.0172, 0.0292]      ala_L_e          0.0252    [-0.0147, 0.294]
    gly_e            0.0157    [0.03, -0.267]        cit_e            0.0187    [-0.012, 0.361]
    pro_L_e          0.012     [0.012, 0.012]        h2o_e            0.015     [-0.371, 1.42]
    leu_L_e          0.012     [0.00387, 0.012]      co2_e            0.015     [0, 1.4]
    crn_e            0.012     [0, 0.012]            dcmp_e           0.0108    [-0.012, 0.076]
    gam_e            0.012     [0, 0.012]            uri_e            0.0102    [-0.012, 0.076]
    o2s_e            0.012     [0, 0.012]            pi_e             0.00699   [-0.04, 0.135]
    dctp_n           0.012     [0.012, -0.0465]      3ivcrn_e         0.00681   [0, 0.012]
    cytd_e           0.012     [0.012, -0.076]       dag_hs_e         0.00652   [-0.012, 0.184]
    ps_hs_e          0.012     [0.012, -0.156]       3deccrn_e        0.00519   [0, 0.012]
    ser_L_e          0.012     [0.012, -0.285]       ac_e             0.000938  [-0.012, 0.516]
    arg_L_e          0.0105    [0.0105, 0.0105]      urate_e          0.000938  [0, 0.0739]
    thr_L_e          0.0091    [0.0091, 0.012]       for_e            0.000594  [-0.012, 0.102]
    ile_L_e          0.00832   [0.00832, 0.012]      fe3_e            0         [-0.012, 0.545]
    asn_L_e          0.00813   [0.00813, 0.012]      oxa_e            0         [0, 0.465]
    phe_L_e          0.00755   [0.00755, 0.012]      vacc_e           0         [0, 0.453]
    arachd_e         0.00519   [0, 0.012]            elaid_e          0         [0, 0.449]
    chol_e           0.0048    [0.012, -0.0134]      no_e             0         [0, 0.448]
    tyr_L_e          0.00464   [0.00464, 0.012]      gal_e            0         [0, 0.42]
    met_L_e          0.00445   [0.00445, 0.012]      xylt_e           0         [-0.012, 0.347]
    4mop_e           0.00387   [0.00387, 0.012]      akg_e            0         [-0.012, 0.344]
    his_L_e          0.00368   [0.00368, 0.012]      ocdcea_e         0         [0, 0.333]
    gsn_e            0.00199   [-0.0109, 0.012]      hdcea_e          0         [0, 0.328]
    pe_hs_e          0.00181   [0.012, -0.159]       acmana_e         0         [0, 0.305]
    biomass_other_c  0.00157   [0.00157, 0.00157]    glyc_S_e         0         [0, 0.297]
    adn_e            0.00156   [0.00911, -0.0214]    glyc_e           0         [-0.21, 0.297]
    inost_e          0.000678  [0.000678, 0.000678]  ha_pre1_e        0         [0, 0.195]
    hdca_e           0.000509  [-0.031, 0.3]         nrvnc_e          0         [0, 0.178]
    sph1p_e          0.000509  [0.012, -0.163]       ethamp_r         0         [0, 0.17]
    pglyc_hs_e       0.000424  [0.000424, 0.000424]  pe_hs_r          0         [-0.012, 0.159]
    trp_L_e          0.000387  [0.000387, 0.012]     thymd_e          0         [-0.012, 0.1]
    datp_n           0.000384  [0.000384, 0.000384]  ha_e             0         [-0.012, 0.0908]
    dttp_n           0.000381  [0.000381, -0.0581]   ppa_e            0         [-0.012, 0.067]
    dgtp_n           0.000288  [0.000288, 0.000288]  clpnd_e          0         [0, 0.06]
    fe2_e            0         [-0.012, 0.545]       dlnlcg_e         0         [0, 0.06]
    ocdca_e          0         [0, 0.07]             eicostet_e       0         [0, 0.06]
    lnlc_e           0         [0, 0.06]             lnlnca_e         0         [0, 0.06]
    gthrd_e          0         [-0.0586, 0.06]       lnlncg_e         0         [0, 0.06]
    utp_e            0         [0, 0.02]             so4_e            0         [0, 0.0586]
    3hpvs_e          0         [0, 0.012]            taur_c           0         [0, 0.0586]
    4hpro_LT_m       0         [0, 0.012]            taur_e           0         [0, 0.0586]
    amp_e            0         [0, 0.012]            tchola_e         0         [0, 0.0564]
    arab_L_e         0         [0, 0.012]            chsterol_e       0         [0, 0.0545]
    asp_L_e          0         [0, 0.012]            4mptnl_e         0         [0, 0.0532]
    crvnc_e          0         [0, 0.012]            aprgstrn_e       0         [0, 0.0532]
    dcyt_e           0         [0, 0.012]            5adtststerone_e  0         [0, 0.0516]
    din_e            0         [0, 0.012]            andrstrn_e       0         [0, 0.0516]
    dopa_e           0         [0, 0.012]            gchola_e         0         [-0.012, 0.0485]
    etoh_e           0         [0, 0.012]            fuc13galacgl...  0         [0, 0.048]
    fmn_e            0         [0, 0.012]            fuc14galacgl...  0         [0, 0.048]
    gluala_e         0         [0, 0.012]            5adtststeron...  0         [0, 0.0459]
    leuktrA4_e       0         [0, 0.012]            andrstrnglc_e    0         [0, 0.0459]
    n2m2nmasn_e      0         [0, 0.012]            tststeroneglc_e  0         [0, 0.0458]
    orn_e            0         [0, 0.012]            fuc_L_e          0         [0, 0.0364]
    ptrc_e           0         [0, 0.012]            man_e            0         [-0.012, 0.036]
    ptvstlac_e       0         [0, 0.012]            Rtotal_e         0         [0, 0.025]
    retn_e           0         [0, 0.012]            lpchol_hs_e      0         [0, 0.025]
    spmd_e           0         [0, 0.012]            pchol_hs_e       0         [0, 0.025]
    nac_e            0         [-0.0116, 0.012]      gthox_e          0         [0, 0.024]
                                                     fald_e           0         [-0.012, 0.024]
                                                     meoh_e           0         [-0.012, 0.024]
                                                     1mncam_e         0         [0, 0.0236]
                                                     ppi_e            0         [0, 0.02]
                                                     ump_e            0         [0, 0.02]
                                                     3mob_e           0         [-0.012, 0.0197]
                                                     fucfucfucgal...  0         [0, 0.016]
                                                     13_cis_retng...  0         [0, 0.012]
                                                     34dhoxpeg_e      0         [0, 0.012]
                                                     3bcrn_e          0         [0, 0.012]
                                                     3ddcrn_e         0         [0, 0.012]
                                                     3hdececrn_e      0         [0, 0.012]
                                                     3hexdcrn_e       0         [0, 0.012]
                                                     3hpvstet_e       0         [0, 0.012]
                                                     3octdec2crn_e    0         [0, 0.012]
                                                     3octdeccrn_e     0         [0, 0.012]
                                                     3octdece1crn_e   0         [0, 0.012]
                                                     3tdcrn_e         0         [0, 0.012]
                                                     3tetd7ecoacrn_e  0         [0, 0.012]
                                                     3thexddcoacrn_e  0         [0, 0.012]
                                                     3ttetddcoacrn_e  0         [0, 0.012]
                                                     4abutn_e         0         [0, 0.012]
                                                     Asn_X_Ser_Thr_l  0         [0, 0.012]
                                                     CE1940_e         0         [0, 0.012]
                                                     abt_e            0         [0, 0.012]
                                                     ade_e            0         [0, 0.012]
                                                     adrn_e           0         [0, 0.012]
                                                     c5dc_e           0         [0, 0.012]
                                                     c6crn_e          0         [0, 0.012]
                                                     leuktrB4_e       0         [0, 0.012]
                                                     leuktrD4_e       0         [0, 0.012]
                                                     n5m2masn_g       0         [0, 0.012]
                                                     prpp_e           0         [0, 0.012]
                                                     ptvst_e          0         [0, 0.012]
                                                     retnglc_e        0         [0, 0.012]
                                                     ribflv_e         0         [0, 0.012]
                                                     HC01609_e        0         [-0.012, 0.012]
                                                     HC01610_e        0         [-0.012, 0.012]
                                                     acac_e           0         [-0.012, 0.012]
                                                     bhb_e            0         [-0.012, 0.012]
                                                     anth_c           0         [0, 0.0116]
                                                     3mlda_e          0         [0, 0.00832]
                                                     ivcrn_e          0         [0, 0.00813]
                                                     5mta_e           0         [0, 0.00755]
                                                     sprm_c           0         [0, 0.00755]
                                                     sprm_e           0         [0, 0.00755]
                                                     4hphac_e         0         [0, 0.00736]
                                                     q10h2_e          0         [0, 0.00736]
                                                     bilglcur_e       0         [0, 0.00547]
                                                     co_e             0         [0, 0.00547]
                                                     pheme_e          0         [0, 0.00547]
                                                     pheacgln_e       0         [0, 0.00445]
                                                     2hb_e            0         [0.0029, -0.012]



```python
model_CancerBiopsy2.summary(fva=True)
```

    IN FLUXES                                        OUT FLUXES                                   OBJECTIVES
    -----------------------------------------------  -------------------------------------------  -----------------------
    id                   Flux  Range                 id                   Flux  Range             biomass_reac...  0.0334
    ---------------  --------  --------------------  ---------------  --------  ----------------
    glc_D_e          0.362     [0.22, 1]             h_e              0.698     [-0.0486, 3.2]
    o2_e             0.0845    [-0.0058, 1]          lac_L_e          0.627     [-0.09, 2.67]
    Lcystin_c        0.024     [-0.00522, 0.024]     co2_e            0.0503    [0, 1.34]
    pi_e             0.0208    [0.04, -0.0836]       h2o_e            0.0503    [-0.37, 1.29]
    lys_L_e          0.0198    [0.0198, 0.03]        pyr_e            0.0474    [-0.012, 1.26]
    gly_e            0.018     [0.03, -0.269]        cys_L_e          0.0339    [-0.012, 0.0464]
    gthrd_e          0.012     [-0.0584, 0.06]       glyc_S_e         0.016     [0, 0.299]
    arg_L_e          0.012     [0.012, 0.012]        ac_e             0.0129    [-0.012, 0.473]
    4mop_e           0.012     [0.00622, 0.012]      so4_e            0.0126    [0, 0.0584]
    asp_L_e          0.012     [0, 0.012]            urate_e          0.012     [0, 0.072]
    crn_e            0.012     [0, 0.012]            34dhoxpeg_e      0.012     [0, 0.012]
    din_e            0.012     [0, 0.012]            leuktrD4_e       0.012     [0, 0.012]
    dopa_e           0.012     [0, 0.012]            3bcrn_e          0.00996   [0, 0.012]
    leuktrA4_e       0.012     [0, 0.012]            uri_e            0.00891   [-0.012, 0.0517]
    o2s_e            0.012     [0, 0.012]            dag_hs_e         0.00665   [-0.012, 0.183]
    3mob_e           0.012     [0.012, -0.0182]      ptvst_e          0.00309   [0, 0.012]
    cytd_e           0.012     [0.012, -0.0517]      3deccrn_e        0.00204   [0, 0.012]
    ps_hs_e          0.012     [0.012, -0.112]       lac_D_e          0.001     [-0.001, 0.001]
    ser_L_e          0.012     [0.012, -0.287]       for_e            0.000681  [-0.012, 0.099]
    gln_L_e          0.0118    [0.09, -0.13]         cit_e            0.000222  [-0.012, 0.327]
    thr_L_e          0.0104    [0.0104, 0.012]       fe3_e            0         [-0.012, 0.545]
    ile_L_e          0.00956   [0.00956, 0.012]      vacc_e           0         [0, 0.449]
    asn_L_e          0.00933   [0.00933, 0.012]      oxa_e            0         [0, 0.448]
    ppa_e            0.00888   [0.012, -0.0616]      elaid_e          0         [0, 0.43]
    phe_L_e          0.00867   [0.00867, 0.012]      nh4_e            0         [-0.012, 0.427]
    leu_L_e          0.00622   [0.00622, 0.012]      gal_e            0         [0, 0.39]
    chol_e           0.00552   [-0.0116, 0.012]      no_e             0         [0, 0.385]
    tyr_L_e          0.00533   [0.00533, 0.012]      hdcea_e          0         [0, 0.326]
    met_L_e          0.00511   [0.00511, 0.012]      ocdcea_e         0         [0, 0.31]
    val_L_e          0.00458   [-0.000222, 0.03]     glyc_e           0         [-0.21, 0.299]
    his_L_e          0.00422   [0.00422, 0.012]      ala_L_e          0         [-0.0169, 0.294]
    ptvstlac_e       0.00309   [0, 0.012]            pro_L_e          0         [-0.012, 0.287]
    pe_hs_e          0.00302   [0.012, -0.112]       acmana_e         0         [0, 0.287]
    arachd_e         0.00204   [0, 0.012]            ha_pre1_e        0         [0, 0.184]
    gam_e            0.00184   [0, 0.012]            nrvnc_e          0         [0, 0.165]
    biomass_other_c  0.0018    [0.0018, 0.0018]      ethamp_r         0         [0, 0.124]
    adn_e            0.00179   [0.00868, -0.021]     pe_hs_r          0         [-0.012, 0.112]
    gsn_e            0.00121   [-0.0108, 0.012]      ha_e             0         [-0.012, 0.0852]
    inost_e          0.000779  [0.000779, 0.000779]  thymd_e          0         [-0.012, 0.0757]
    hdca_e           0.000584  [-0.0291, 0.3]        clpnd_e          0         [0, 0.06]
    sph1p_e          0.000584  [0.012, -0.112]       dlnlcg_e         0         [0, 0.06]
    pglyc_hs_e       0.000487  [0.000487, 0.000487]  eicostet_e       0         [0, 0.06]
    trp_L_e          0.000444  [0.000444, 0.012]     lnlnca_e         0         [0, 0.06]
    datp_n           0.00044   [0.00044, 0.00044]    lnlncg_e         0         [0, 0.06]
    dttp_n           0.000437  [0.000437, -0.0408]   taur_c           0         [0, 0.0584]
    dgtp_n           0.000331  [0.000331, 0.000331]  taur_e           0         [0, 0.0584]
    dctp_n           0.000315  [0.000315, 0.000315]  tchola_e         0         [0, 0.0536]
    fe2_e            0         [-0.012, 0.545]       chsterol_e       0         [0, 0.0513]
    ocdca_e          0         [0, 0.07]             4mptnl_e         0         [0, 0.0501]
    lnlc_e           0         [0, 0.06]             aprgstrn_e       0         [0, 0.0501]
    utp_e            0         [0, 0.02]             5adtststerone_e  0         [0, 0.0486]
    3hpvs_e          0         [0, 0.012]            andrstrn_e       0         [0, 0.0486]
    4hpro_LT_m       0         [0, 0.012]            fuc13galacgl...  0         [0, 0.048]
    amp_e            0         [0, 0.012]            fuc14galacgl...  0         [0, 0.048]
    arab_L_e         0         [0, 0.012]            gchola_e         0         [-0.012, 0.0456]
    crvnc_e          0         [0, 0.012]            5adtststeron...  0         [0, 0.0431]
    dcyt_e           0         [0, 0.012]            andrstrnglc_e    0         [0, 0.0431]
    etoh_e           0         [0, 0.012]            tststeroneglc_e  0         [0, 0.0429]
    fmn_e            0         [0, 0.012]            man_e            0         [-0.012, 0.036]
    gluala_e         0         [0, 0.012]            fuc_L_e          0         [0, 0.0338]
    n2m2nmasn_e      0         [0, 0.012]            13_cis_retng...  0         [0, 0.024]
    ptrc_e           0         [0, 0.012]            gthox_e          0         [0, 0.024]
    retinol_e        0         [0, 0.012]            retnglc_e        0         [0, 0.024]
    spmd_e           0         [0, 0.012]            fald_e           0         [-0.012, 0.024]
    2hb_e            0         [-0.00156, 0.012]     meoh_e           0         [-0.012, 0.024]
    nac_e            0         [-0.0116, 0.012]      1mncam_e         0         [0, 0.0236]
    xylt_e           0         [0.012, -0.319]       glyb_e           0         [0, 0.0235]
                                                     Rtotal_e         0         [0, 0.0231]
                                                     lpchol_hs_e      0         [0, 0.0231]
                                                     pchol_hs_e       0         [0, 0.0231]
                                                     ppi_e            0         [0, 0.02]
                                                     ump_e            0         [0, 0.02]
                                                     fucfucfucgal...  0         [0, 0.016]
                                                     3ddcrn_e         0         [0, 0.012]
                                                     3hdececrn_e      0         [0, 0.012]
                                                     3hexdcrn_e       0         [0, 0.012]
                                                     3hpvstet_e       0         [0, 0.012]
                                                     3ivcrn_e         0         [0, 0.012]
                                                     3octdec2crn_e    0         [0, 0.012]
                                                     3octdeccrn_e     0         [0, 0.012]
                                                     3octdece1crn_e   0         [0, 0.012]
                                                     3tdcrn_e         0         [0, 0.012]
                                                     3tetd7ecoacrn_e  0         [0, 0.012]
                                                     3thexddcoacrn_e  0         [0, 0.012]
                                                     3ttetddcoacrn_e  0         [0, 0.012]
                                                     4abutn_e         0         [0, 0.012]
                                                     Asn_X_Ser_Thr_l  0         [0, 0.012]
                                                     CE1940_e         0         [0, 0.012]
                                                     abt_e            0         [0, 0.012]
                                                     ade_e            0         [0, 0.012]
                                                     adrn_e           0         [0, 0.012]
                                                     c6crn_e          0         [0, 0.012]
                                                     leuktrB4_e       0         [0, 0.012]
                                                     n5m2masn_g       0         [0, 0.012]
                                                     prpp_e           0         [0, 0.012]
                                                     ribflv_e         0         [0, 0.012]
                                                     HC01609_e        0         [-0.012, 0.012]
                                                     HC01610_e        0         [-0.012, 0.012]
                                                     acac_e           0         [-0.012, 0.012]
                                                     bhb_e            0         [-0.012, 0.012]
                                                     retn_e           0         [-0.012, 0.012]
                                                     anth_c           0         [0, 0.0116]
                                                     c5dc_e           0         [0, 0.0102]
                                                     3mlda_e          0         [0, 0.00778]
                                                     5mta_e           0         [0, 0.00689]
                                                     sprm_c           0         [0, 0.00689]
                                                     sprm_e           0         [0, 0.00689]
                                                     4hphac_e         0         [0, 0.00667]
                                                     q10h2_e          0         [0, 0.00667]
                                                     ivcrn_e          0         [0, 0.00578]
                                                     bilglcur_e       0         [0, 0.00528]
                                                     co_e             0         [0, 0.00528]
                                                     pheme_e          0         [0, 0.00528]
                                                     pheacgln_e       0         [0, 0.00333]



```python
model_CancerBiopsy3.summary(fva=True)
```

    IN FLUXES                                        OUT FLUXES                                    OBJECTIVES
    -----------------------------------------------  --------------------------------------------  -----------------------
    id                   Flux  Range                 id                   Flux  Range              biomass_reac...  0.0334
    ---------------  --------  --------------------  ---------------  --------  -----------------
    glc_D_e          0.329     [0.197, 1]            h_e              0.716     [-0.411, 3.29]
    o2_e             0.137     [-0.0058, 1]          lac_L_e          0.578     [-0.09, 2.7]
    gln_L_e          0.0462    [0.09, -0.153]        cit_e            0.0521    [-0.012, 0.348]
    lys_L_e          0.0198    [0.0198, 0.03]        h2o_e            0.0312    [-0.371, 1.69]
    val_L_e          0.019     [-0.000222, 0.03]     ala_L_e          0.0291    [-0.0169, 0.294]
    arg_L_e          0.012     [0.012, 0.012]        co2_e            0.0288    [0, 1.73]
    4mop_e           0.012     [0.00622, 0.012]      ac_e             0.0281    [-0.012, 0.504]
    crn_e            0.012     [0, 0.012]            urate_e          0.0281    [0, 0.0847]
    din_e            0.012     [0, 0.012]            nh4_e            0.0245    [-0.012, 0.475]
    o2s_e            0.012     [0, 0.012]            uri_e            0.0197    [-0.012, 0.0634]
    acac_e           0.012     [0.012, -0.012]       pi_e             0.0133    [-0.04, 0.119]
    3mob_e           0.012     [0.012, -0.0182]      pyr_e            0.0122    [-0.012, 1.44]
    gsn_e            0.012     [0.012, -0.0228]      bhb_e            0.012     [0.012, -0.012]
    dctp_n           0.012     [0.012, -0.0409]      3bcrn_e          0.00968   [0, 0.012]
    cys_L_e          0.012     [0.012, -0.0464]      dag_hs_e         0.00571   [-0.012, 0.183]
    ppa_e            0.012     [0.012, -0.0625]      Lcystin_c        0.00522   [0.00522, -0.024]
    cytd_e           0.012     [0.012, -0.0634]      glyc_e           0.00238   [-0.21, 0.299]
    ps_hs_e          0.012     [0.012, -0.147]       3deccrn_e        0.00214   [0, 0.012]
    pro_L_e          0.012     [0.012, -0.287]       lac_D_e          0.001     [-0.001, 0.001]
    ser_L_e          0.012     [0.012, -0.287]       thymd_e          0.000891  [-0.012, 0.0994]
    akg_e            0.012     [0.012, -0.331]       3ivcrn_e         0.000187  [0, 0.012]
    thr_L_e          0.0104    [0.0104, 0.012]       fe3_e            0         [-0.012, 0.545]
    xylt_e           0.00969   [0.012, -0.328]       oxa_e            0         [0, 0.452]
    ile_L_e          0.00956   [0.00956, 0.012]      vacc_e           0         [0, 0.452]
    asn_L_e          0.00933   [0.00933, 0.012]      elaid_e          0         [0, 0.449]
    phe_L_e          0.00867   [0.00867, 0.012]      no_e             0         [0, 0.429]
    leu_L_e          0.00622   [0.00622, 0.012]      gal_e            0         [0, 0.401]
    chol_e           0.00552   [-0.0119, 0.012]      ocdcea_e         0         [0, 0.349]
    tyr_L_e          0.00533   [0.00533, 0.012]      hdcea_e          0         [0, 0.328]
    dgsn_e           0.00526   [0, 0.012]            glyc_S_e         0         [0, 0.299]
    met_L_e          0.00511   [0.00511, 0.012]      acmana_e         0         [0, 0.294]
    his_L_e          0.00422   [0.00422, 0.012]      gly_e            0         [-0.03, 0.269]
    arachd_e         0.00214   [0, 0.012]            ha_pre1_e        0         [0, 0.193]
    pe_hs_e          0.00208   [0.012, -0.147]       nrvnc_e          0         [0, 0.175]
    biomass_other_c  0.0018    [0.0018, 0.0018]      ethamp_r         0         [0, 0.159]
    adn_e            0.00179   [0.00868, -0.045]     pe_hs_r          0         [-0.012, 0.147]
    inost_e          0.000779  [0.000779, 0.000779]  ha_e             0         [-0.012, 0.0898]
    hdca_e           0.000584  [-0.0294, 0.3]        clpnd_e          0         [0, 0.06]
    sph1p_e          0.000584  [0.012, -0.147]       dlnlcg_e         0         [0, 0.06]
    pglyc_hs_e       0.000487  [0.000487, 0.000487]  eicostet_e       0         [0, 0.06]
    trp_L_e          0.000444  [0.000444, 0.012]     lnlnca_e         0         [0, 0.06]
    datp_n           0.00044   [0.00044, 0.00044]    lnlncg_e         0         [0, 0.06]
    dttp_n           0.000437  [0.000437, -0.0524]   so4_e            0         [0, 0.0584]
    dgtp_n           0.000331  [0.000331, 0.000331]  taur_c           0         [0, 0.0584]
    for_e            0.00021   [0.012, -0.0999]      taur_e           0         [0, 0.0584]
    fe2_e            0         [-0.012, 0.545]       tchola_e         0         [0, 0.0546]
    ocdca_e          0         [0, 0.07]             chsterol_e       0         [0, 0.0524]
    gthrd_e          0         [-0.0584, 0.06]       4mptnl_e         0         [0, 0.0511]
    utp_e            0         [0, 0.02]             aprgstrn_e       0         [0, 0.0511]
    3hpvs_e          0         [0, 0.012]            5adtststerone_e  0         [0, 0.0495]
    4hpro_LT_m       0         [0, 0.012]            andrstrn_e       0         [0, 0.0495]
    amp_e            0         [0, 0.012]            fuc13galacgl...  0         [0, 0.048]
    arab_L_e         0         [0, 0.012]            fuc14galacgl...  0         [0, 0.048]
    asp_L_e          0         [0, 0.012]            gchola_e         0         [-0.012, 0.0465]
    crvnc_e          0         [0, 0.012]            galfucgalacg...  0         [0, 0.0456]
    dad_2_e          0         [0, 0.012]            5adtststeron...  0         [0, 0.0441]
    dcyt_e           0         [0, 0.012]            andrstrnglc_e    0         [0, 0.0441]
    dopa_e           0         [0, 0.012]            tststeroneglc_e  0         [0, 0.044]
    etoh_e           0         [0, 0.012]            man_e            0         [-0.012, 0.036]
    fmn_e            0         [0, 0.012]            fuc_L_e          0         [0, 0.0347]
    gam_e            0         [0, 0.012]            13_cis_retng...  0         [0, 0.024]
    gluala_e         0         [0, 0.012]            gthox_e          0         [0, 0.024]
    leuktrA4_e       0         [0, 0.012]            retnglc_e        0         [0, 0.024]
    n2m2nmasn_e      0         [0, 0.012]            fald_e           0         [-0.012, 0.024]
    orn_e            0         [0, 0.012]            meoh_e           0         [-0.012, 0.024]
    ptrc_e           0         [0, 0.012]            glyb_e           0         [0, 0.0238]
    ptvstlac_e       0         [0, 0.012]            1mncam_e         0         [0, 0.0236]
    retinol_e        0         [0, 0.012]            Rtotal_e         0         [0, 0.0235]
    spmd_e           0         [0, 0.012]            lpchol_hs_e      0         [0, 0.0235]
    2hb_e            0         [-0.00156, 0.012]     pchol_hs_e       0         [0, 0.0235]
    retn_e           0         [-0.012, 0.012]       ppi_e            0         [0, 0.02]
                                                     ump_e            0         [0, 0.02]
                                                     fucfucfucgal...  0         [0, 0.016]
                                                     34dhoxpeg_e      0         [0, 0.012]
                                                     3ddcrn_e         0         [0, 0.012]
                                                     3hdececrn_e      0         [0, 0.012]
                                                     3hexdcrn_e       0         [0, 0.012]
                                                     3hpvstet_e       0         [0, 0.012]
                                                     3octdec2crn_e    0         [0, 0.012]
                                                     3octdeccrn_e     0         [0, 0.012]
                                                     3octdece1crn_e   0         [0, 0.012]
                                                     3tdcrn_e         0         [0, 0.012]
                                                     3tetd7ecoacrn_e  0         [0, 0.012]
                                                     3thexddcoacrn_e  0         [0, 0.012]
                                                     3ttetddcoacrn_e  0         [0, 0.012]
                                                     4abutn_e         0         [0, 0.012]
                                                     Asn_X_Ser_Thr_l  0         [0, 0.012]
                                                     CE1940_e         0         [0, 0.012]
                                                     abt_e            0         [0, 0.012]
                                                     ade_e            0         [0, 0.012]
                                                     adrn_e           0         [0, 0.012]
                                                     c6crn_e          0         [0, 0.012]
                                                     leuktrB4_e       0         [0, 0.012]
                                                     leuktrD4_e       0         [0, 0.012]
                                                     n5m2masn_g       0         [0, 0.012]
                                                     prpp_e           0         [0, 0.012]
                                                     ptvst_e          0         [0, 0.012]
                                                     ribflv_e         0         [0, 0.012]
                                                     HC01609_e        0         [-0.012, 0.012]
                                                     HC01610_e        0         [-0.012, 0.012]
                                                     anth_c           0         [0, 0.0116]
                                                     c5dc_e           0         [0, 0.0102]
                                                     3mlda_e          0         [0, 0.00778]
                                                     5mta_e           0         [0, 0.00689]
                                                     sprm_c           0         [0, 0.00689]
                                                     sprm_e           0         [0, 0.00689]
                                                     4hphac_e         0         [0, 0.00667]
                                                     q10h2_e          0         [0, 0.00667]
                                                     ivcrn_e          0         [0, 0.00578]
                                                     bilglcur_e       0         [0, 0.00528]
                                                     co_e             0         [0, 0.00528]
                                                     pheme_e          0         [0, 0.00528]
                                                     pheacgln_e       0         [0, 0.00333]
                                                     nac_e            0         [0.0116, -0.012]
                                                     lnlc_e           0         [0, -0.06]



```python
model_CancerBiopsy4.summary(fva=True)
```

    IN FLUXES                                        OUT FLUXES                                   OBJECTIVES
    -----------------------------------------------  -------------------------------------------  -----------------------
    id                   Flux  Range                 id                   Flux  Range             biomass_reac...  0.0334
    ---------------  --------  --------------------  ---------------  --------  ----------------
    glc_D_e          0.352     [0.22, 1]             h_e              0.7       [-0.357, 3.21]
    o2_e             0.0656    [-0.0058, 1]          lac_L_e          0.653     [-0.09, 2.68]
    val_L_e          0.03      [-0.000222, 0.03]     cys_L_e          0.0464    [-0.012, 0.0464]
    Lcystin_c        0.024     [-0.00522, 0.024]     pyr_e            0.0361    [-0.012, 1.33]
    pi_e             0.0218    [0.04, -0.0836]       co2_e            0.0315    [0, 1.71]
    lys_L_e          0.0198    [0.0198, 0.03]        h2o_e            0.0315    [-0.37, 1.68]
    gly_e            0.018     [0.03, -0.269]        ac_e             0.012     [-0.012, 0.477]
    gln_L_e          0.0158    [0.09, -0.142]        urate_e          0.012     [0, 0.0669]
    arg_L_e          0.012     [0.012, 0.012]        bhb_e            0.012     [0.012, -0.012]
    leu_L_e          0.012     [0.00622, 0.012]      glyc_S_e         0.0107    [0, 0.299]
    din_e            0.012     [0, 0.012]            leuktrD4_e       0.00979   [0, 0.012]
    leuktrA4_e       0.012     [0, 0.012]            dag_hs_e         0.00571   [-0.012, 0.183]
    o2s_e            0.012     [0, 0.012]            leuktrB4_e       0.00221   [0, 0.012]
    acac_e           0.012     [-0.012, 0.012]       nh4_e            0.00204   [-0.012, 0.451]
    ps_hs_e          0.012     [0.012, -0.112]       3deccrn_e        0.00204   [0, 0.012]
    pro_L_e          0.012     [0.012, -0.287]       lac_D_e          0.001     [-0.001, 0.001]
    ser_L_e          0.012     [0.012, -0.287]       for_e            0.000681  [-0.012, 0.0871]
    thr_L_e          0.0104    [0.0104, 0.012]       fe3_e            0         [-0.012, 0.545]
    gthrd_e          0.00979   [-0.0584, 0.06]       elaid_e          0         [0, 0.449]
    ile_L_e          0.00956   [0.00956, 0.012]      oxa_e            0         [0, 0.446]
    asn_L_e          0.00933   [0.00933, 0.012]      no_e             0         [0, 0.406]
    phe_L_e          0.00867   [0.00867, 0.012]      gal_e            0         [0, 0.39]
    4mop_e           0.00622   [0.00622, 0.012]      ocdcea_e         0         [0, 0.34]
    chol_e           0.00552   [-0.0116, 0.012]      cit_e            0         [-0.012, 0.335]
    tyr_L_e          0.00533   [0.00533, 0.012]      hdcea_e          0         [0, 0.327]
    met_L_e          0.00511   [0.00511, 0.012]      xylt_e           0         [-0.012, 0.319]
    his_L_e          0.00422   [0.00422, 0.012]      glyc_e           0         [-0.21, 0.299]
    pe_hs_e          0.00208   [0.012, -0.112]       acmana_e         0         [0, 0.287]
    arachd_e         0.00204   [0, 0.012]            ha_pre1_e        0         [0, 0.185]
    crn_e            0.00204   [0, 0.012]            nrvnc_e          0         [0, 0.17]
    gam_e            0.00182   [0, 0.012]            ethamp_r         0         [0, 0.124]
    biomass_other_c  0.0018    [0.0018, 0.0018]      pe_hs_r          0         [-0.012, 0.112]
    adn_e            0.00179   [0.00868, -0.021]     ha_e             0         [-0.012, 0.086]
    uri_e            0.00178   [0.012, -0.0517]      thymd_e          0         [-0.012, 0.0757]
    ala_L_e          0.00131   [0.0169, -0.294]      ppa_e            0         [-0.012, 0.0616]
    cytd_e           0.0013    [0.012, -0.0517]      clpnd_e          0         [0, 0.06]
    gsn_e            0.00121   [-0.0108, 0.012]      dlnlcg_e         0         [0, 0.06]
    inost_e          0.000779  [0.000779, 0.000779]  eicostet_e       0         [0, 0.06]
    hdca_e           0.000584  [-0.0291, 0.3]        lnlnca_e         0         [0, 0.06]
    sph1p_e          0.000584  [0.012, -0.112]       lnlncg_e         0         [0, 0.06]
    pglyc_hs_e       0.000487  [0.000487, 0.000487]  so4_e            0         [0, 0.0584]
    trp_L_e          0.000444  [0.000444, 0.012]     taur_c           0         [0, 0.0584]
    datp_n           0.00044   [0.00044, 0.00044]    taur_e           0         [0, 0.0584]
    dttp_n           0.000437  [0.000437, -0.0408]   tchola_e         0         [0, 0.0535]
    dgtp_n           0.000331  [0.000331, 0.000331]  chsterol_e       0         [0, 0.0512]
    dctp_n           0.000315  [0.000315, 0.000315]  4mptnl_e         0         [0, 0.0501]
    fe2_e            0         [-0.012, 0.545]       aprgstrn_e       0         [0, 0.0501]
    ocdca_e          0         [0, 0.07]             5adtststerone_e  0         [0, 0.0486]
    lnlc_e           0         [0, 0.06]             andrstrn_e       0         [0, 0.0486]
    utp_e            0         [0, 0.02]             fuc13galacgl...  0         [0, 0.048]
    3hpvs_e          0         [0, 0.012]            fuc14galacgl...  0         [0, 0.048]
    4hpro_LT_m       0         [0, 0.012]            gchola_e         0         [-0.012, 0.0456]
    amp_e            0         [0, 0.012]            galfucgalacg...  0         [0, 0.0443]
    arab_L_e         0         [0, 0.012]            5adtststeron...  0         [0, 0.0431]
    asp_L_e          0         [0, 0.012]            andrstrnglc_e    0         [0, 0.0431]
    crvnc_e          0         [0, 0.012]            tststeroneglc_e  0         [0, 0.0429]
    dcyt_e           0         [0, 0.012]            man_e            0         [-0.012, 0.036]
    dopa_e           0         [0, 0.012]            fuc_L_e          0         [0, 0.0338]
    etoh_e           0         [0, 0.012]            13_cis_retng...  0         [0, 0.024]
    fmn_e            0         [0, 0.012]            gthox_e          0         [0, 0.024]
    gluala_e         0         [0, 0.012]            retnglc_e        0         [0, 0.024]
    n2m2nmasn_e      0         [0, 0.012]            1mncam_e         0         [0, 0.0236]
    orn_e            0         [0, 0.012]            glyb_e           0         [0, 0.0234]
    ptrc_e           0         [0, 0.012]            Rtotal_e         0         [0, 0.0231]
    ptvstlac_e       0         [0, 0.012]            lpchol_hs_e      0         [0, 0.0231]
    retinol_e        0         [0, 0.012]            pchol_hs_e       0         [0, 0.0231]
    spmd_e           0         [0, 0.012]            ppi_e            0         [0, 0.02]
    2hb_e            0         [-0.00156, 0.012]     ump_e            0         [0, 0.02]
    vacc_e           0         [0, -0.452]           3mob_e           0         [-0.012, 0.0182]
                                                     fucfucfucgal...  0         [0, 0.016]
                                                     34dhoxpeg_e      0         [0, 0.012]
                                                     3bcrn_e          0         [0, 0.012]
                                                     3ddcrn_e         0         [0, 0.012]
                                                     3hdececrn_e      0         [0, 0.012]
                                                     3hexdcrn_e       0         [0, 0.012]
                                                     3hpvstet_e       0         [0, 0.012]
                                                     3ivcrn_e         0         [0, 0.012]
                                                     3octdec2crn_e    0         [0, 0.012]
                                                     3octdeccrn_e     0         [0, 0.012]
                                                     3octdece1crn_e   0         [0, 0.012]
                                                     3tdcrn_e         0         [0, 0.012]
                                                     3tetd7ecoacrn_e  0         [0, 0.012]
                                                     3thexddcoacrn_e  0         [0, 0.012]
                                                     3ttetddcoacrn_e  0         [0, 0.012]
                                                     4abutn_e         0         [0, 0.012]
                                                     Asn_X_Ser_Thr_l  0         [0, 0.012]
                                                     CE1940_e         0         [0, 0.012]
                                                     abt_e            0         [0, 0.012]
                                                     ade_e            0         [0, 0.012]
                                                     adrn_e           0         [0, 0.012]
                                                     c6crn_e          0         [0, 0.012]
                                                     n5m2masn_g       0         [0, 0.012]
                                                     prpp_e           0         [0, 0.012]
                                                     ptvst_e          0         [0, 0.012]
                                                     ribflv_e         0         [0, 0.012]
                                                     HC01609_e        0         [-0.012, 0.012]
                                                     HC01610_e        0         [-0.012, 0.012]
                                                     fald_e           0         [-0.012, 0.012]
                                                     retn_e           0         [-0.012, 0.012]
                                                     anth_c           0         [0, 0.0116]
                                                     c5dc_e           0         [0, 0.0102]
                                                     3mlda_e          0         [0, 0.00778]
                                                     5mta_e           0         [0, 0.00689]
                                                     sprm_c           0         [0, 0.00689]
                                                     sprm_e           0         [0, 0.00689]
                                                     4hphac_e         0         [0, 0.00667]
                                                     q10h2_e          0         [0, 0.00667]
                                                     ivcrn_e          0         [0, 0.00578]
                                                     bilglcur_e       0         [0, 0.00528]
                                                     co_e             0         [0, 0.00528]
                                                     pheme_e          0         [0, 0.00528]
                                                     pheacgln_e       0         [0, 0.00333]
                                                     nac_e            0         [0.0116, -0.012]



```python
model_CancerBiopsy5.summary(fva=True)
```

    IN FLUXES                                        OUT FLUXES                                   OBJECTIVES
    -----------------------------------------------  -------------------------------------------  -----------------------
    id                   Flux  Range                 id                   Flux  Range             biomass_reac...  0.0334
    ---------------  --------  --------------------  ---------------  --------  ----------------
    glc_D_e          0.365     [0.22, 1]             h_e              0.725     [-0.0545, 3.21]
    gln_L_e          0.09      [0.09, -0.13]         lac_L_e          0.524     [-0.09, 2.68]
    o2_e             0.0694    [-0.0058, 1]          cit_e            0.0986    [-0.012, 0.339]
    Lcystin_c        0.024     [-0.00522, 0.024]     pyr_e            0.0872    [-0.012, 1.27]
    pi_e             0.0218    [0.04, -0.0836]       nh4_e            0.087     [-0.012, 0.427]
    lys_L_e          0.0198    [0.0198, 0.03]        ala_L_e          0.06      [-0.0169, 0.294]
    gly_e            0.018     [0.03, -0.269]        cys_L_e          0.0464    [-0.012, 0.0464]
    gthrd_e          0.012     [-0.0584, 0.06]       glyc_e           0.0315    [-0.21, 0.299]
    arg_L_e          0.012     [0.012, 0.012]        pro_L_e          0.012     [-0.012, 0.287]
    asn_L_e          0.012     [0.00933, 0.012]      leuktrD4_e       0.012     [0, 0.012]
    leu_L_e          0.012     [0.00622, 0.012]      ac_e             0.00983   [-0.012, 0.473]
    asp_L_e          0.012     [0, 0.012]            urate_e          0.00983   [0, 0.0721]
    leuktrA4_e       0.012     [0, 0.012]            HC01610_e        0.00788   [-0.012, 0.012]
    o2s_e            0.012     [0, 0.012]            dag_hs_e         0.00571   [-0.012, 0.183]
    3mob_e           0.012     [0.012, -0.0182]      3bcrn_e          0.00312   [0, 0.012]
    ppa_e            0.012     [0.012, -0.0617]      3deccrn_e        0.00204   [0, 0.012]
    ps_hs_e          0.012     [0.012, -0.112]       lac_D_e          0.001     [-0.001, 0.001]
    ser_L_e          0.012     [0.012, -0.287]       for_e            0.000681  [-0.012, 0.0992]
    akg_e            0.012     [0.012, -0.324]       co2_e            0         [0, 1.34]
    val_L_e          0.0118    [-0.000222, 0.03]     h2o_e            0         [-0.37, 1.33]
    thr_L_e          0.0104    [0.0104, 0.012]       fe3_e            0         [-0.012, 0.545]
    din_e            0.00983   [0, 0.012]            vacc_e           0         [0, 0.45]
    ile_L_e          0.00956   [0.00956, 0.012]      oxa_e            0         [0, 0.449]
    phe_L_e          0.00867   [0.00867, 0.012]      elaid_e          0         [0, 0.44]
    gam_e            0.00794   [0, 0.012]            no_e             0         [0, 0.401]
    HC01609_e        0.00788   [0.012, -0.012]       gal_e            0         [0, 0.39]
    gluala_e         0.00717   [0, 0.012]            hdcea_e          0         [0, 0.326]
    4mop_e           0.00622   [0.00622, 0.012]      ocdcea_e         0         [0, 0.314]
    chol_e           0.00552   [-0.0116, 0.012]      glyc_S_e         0         [0, 0.299]
    tyr_L_e          0.00533   [0.00533, 0.012]      acmana_e         0         [0, 0.287]
    crn_e            0.00517   [0, 0.012]            ha_pre1_e        0         [0, 0.184]
    met_L_e          0.00511   [0.00511, 0.012]      nrvnc_e          0         [0, 0.166]
    his_L_e          0.00422   [0.00422, 0.012]      ethamp_r         0         [0, 0.124]
    pe_hs_e          0.00208   [0.012, -0.112]       pe_hs_r          0         [-0.012, 0.112]
    arachd_e         0.00204   [0, 0.012]            ha_e             0         [-0.012, 0.0856]
    biomass_other_c  0.0018    [0.0018, 0.0018]      thymd_e          0         [-0.012, 0.0757]
    adn_e            0.00179   [0.00868, -0.021]     clpnd_e          0         [0, 0.06]
    uri_e            0.00178   [0.012, -0.0517]      dlnlcg_e         0         [0, 0.06]
    cytd_e           0.0013    [0.012, -0.0517]      lnlnca_e         0         [0, 0.06]
    gsn_e            0.00121   [-0.0108, 0.012]      lnlncg_e         0         [0, 0.06]
    inost_e          0.000779  [0.000779, 0.000779]  so4_e            0         [0, 0.0584]
    hdca_e           0.000584  [-0.0291, 0.3]        taur_c           0         [0, 0.0584]
    sph1p_e          0.000584  [0.012, -0.112]       taur_e           0         [0, 0.0584]
    pglyc_hs_e       0.000487  [0.000487, 0.000487]  tchola_e         0         [0, 0.0537]
    trp_L_e          0.000444  [0.000444, 0.012]     chsterol_e       0         [0, 0.0514]
    datp_n           0.00044   [0.00044, 0.00044]    4mptnl_e         0         [0, 0.0501]
    dttp_n           0.000437  [0.000437, -0.0408]   aprgstrn_e       0         [0, 0.0501]
    dgtp_n           0.000331  [0.000331, 0.000331]  5adtststerone_e  0         [0, 0.0487]
    dctp_n           0.000315  [0.000315, 0.000315]  andrstrn_e       0         [0, 0.0487]
    ocdca_e          0         [0, 0.07]             fuc13galacgl...  0         [0, 0.048]
    utp_e            0         [0, 0.02]             fuc14galacgl...  0         [0, 0.048]
    3hpvs_e          0         [0, 0.012]            gchola_e         0         [-0.012, 0.0457]
    4hpro_LT_m       0         [0, 0.012]            5adtststeron...  0         [0, 0.0431]
    amp_e            0         [0, 0.012]            andrstrnglc_e    0         [0, 0.0431]
    arab_L_e         0         [0, 0.012]            tststeroneglc_e  0         [0, 0.043]
    crvnc_e          0         [0, 0.012]            galfucgalacg...  0         [0, 0.0424]
    dcyt_e           0         [0, 0.012]            man_e            0         [-0.012, 0.036]
    dopa_e           0         [0, 0.012]            fuc_L_e          0         [0, 0.0338]
    etoh_e           0         [0, 0.012]            13_cis_retng...  0         [0, 0.024]
    fmn_e            0         [0, 0.012]            gthox_e          0         [0, 0.024]
    n2m2nmasn_e      0         [0, 0.012]            retnglc_e        0         [0, 0.024]
    ptrc_e           0         [0, 0.012]            fald_e           0         [-0.012, 0.024]
    ptvstlac_e       0         [0, 0.012]            meoh_e           0         [-0.012, 0.024]
    retinol_e        0         [0, 0.012]            1mncam_e         0         [0, 0.0236]
    spmd_e           0         [0, 0.012]            glyb_e           0         [0, 0.0235]
    2hb_e            0         [-0.00156, 0.012]     Rtotal_e         0         [0, 0.0231]
    nac_e            0         [-0.0116, 0.012]      lpchol_hs_e      0         [0, 0.0231]
    pheme_e          0         [0, -0.00528]         pchol_hs_e       0         [0, 0.0231]
    retn_e           0         [0.012, -0.012]       ppi_e            0         [0, 0.02]
    eicostet_e       0         [0, -0.06]            udp_e            0         [0, 0.02]
    xylt_e           0         [0.012, -0.319]       ump_e            0         [0, 0.02]
                                                     fucfucfucgal...  0         [0, 0.016]
                                                     34dhoxpeg_e      0         [0, 0.012]
                                                     3ddcrn_e         0         [0, 0.012]
                                                     3hdececrn_e      0         [0, 0.012]
                                                     3hexdcrn_e       0         [0, 0.012]
                                                     3hpvstet_e       0         [0, 0.012]
                                                     3ivcrn_e         0         [0, 0.012]
                                                     3octdec2crn_e    0         [0, 0.012]
                                                     3octdeccrn_e     0         [0, 0.012]
                                                     3octdece1crn_e   0         [0, 0.012]
                                                     3tdcrn_e         0         [0, 0.012]
                                                     3tetd7ecoacrn_e  0         [0, 0.012]
                                                     3thexddcoacrn_e  0         [0, 0.012]
                                                     3ttetddcoacrn_e  0         [0, 0.012]
                                                     4abutn_e         0         [0, 0.012]
                                                     Asn_X_Ser_Thr_l  0         [0, 0.012]
                                                     CE1940_e         0         [0, 0.012]
                                                     abt_e            0         [0, 0.012]
                                                     ade_e            0         [0, 0.012]
                                                     adrn_e           0         [0, 0.012]
                                                     c6crn_e          0         [0, 0.012]
                                                     leuktrB4_e       0         [0, 0.012]
                                                     n5m2masn_g       0         [0, 0.012]
                                                     prpp_e           0         [0, 0.012]
                                                     ptvst_e          0         [0, 0.012]
                                                     ribflv_e         0         [0, 0.012]
                                                     acac_e           0         [-0.012, 0.012]
                                                     bhb_e            0         [-0.012, 0.012]
                                                     anth_c           0         [0, 0.0116]
                                                     c5dc_e           0         [0, 0.0102]
                                                     3mlda_e          0         [0, 0.00778]
                                                     5mta_e           0         [0, 0.00689]
                                                     sprm_c           0         [0, 0.00689]
                                                     sprm_e           0         [0, 0.00689]
                                                     4hphac_e         0         [0, 0.00667]
                                                     q10h2_e          0         [0, 0.00667]
                                                     ivcrn_e          0         [0, 0.00578]
                                                     bilglcur_e       0         [0, 0.00528]
                                                     co_e             0         [0, 0.00528]
                                                     pheacgln_e       0         [0, 0.00333]
                                                     lnlc_e           0         [0, -0.06]
                                                     fe2_e            0         [0.012, -0.545]



```python
model_NormalBiopsy1.summary(fva=True)
```

    IN FLUXES                                        OUT FLUXES                                      OBJECTIVES
    -----------------------------------------------  ----------------------------------------------  ----------------------
    id                   Flux  Range                 id                   Flux  Range                biomass_reac...  0.022
    ---------------  --------  --------------------  ---------------  --------  -------------------
    o2_e             0.227     [0.193, 1]            h_e              0.372     [-0.954, 3.14]
    glc_D_e          0.116     [0.00605, 0.874]      glyc_S_e         0.18      [0, 0.786]
    lac_L_e          0.09      [-0.072, 0.09]        ala_L_e          0.168     [-0.04, 0.359]
    gthrd_e          0.0838    [-0.023, 0.09]        mal_L_e          0.113     [-0.012, 0.289]
    glyc_e           0.0427    [-0.012, 0.21]        cgly_e           0.0802    [-0.012, 0.101]
    atp_e            0.0407    [-0.0325, 0.8]        ac_e             0.0667    [-0.012, 0.15]
    lys_L_e          0.013     [0.013, 0.013]        xmp_e            0.0388    [-0.000794, 0.832]
    leu_L_e          0.012     [0.012, 0.012]        pyr_e            0.0376    [-0.012, 0.15]
    ps_hs_e          0.012     [0.00391, 0.012]      oxa_e            0.012     [0, 0.488]
    3mob_e           0.012     [0, 0.012]            3bcrn_e          0.012     [0, 0.012]
    4hpro_LT_m       0.012     [0, 0.012]            dctp_n           0.0118    [-0.012, 0.0538]
    asp_L_e          0.012     [0, 0.012]            utp_e            0.0108    [-0.02, 0.0458]
    crn_e            0.012     [0, 0.012]            pi_e             0.0108    [0, 1.76]
    gam_e            0.012     [0, 0.012]            ppi_e            0.0108    [0, 0.586]
    glyc3p_e         0.012     [-0.0119, 0.012]      dag_hs_e         0.00785   [-0.000385, 0.0235]
    dcmp_e           0.012     [0.012, -0.0538]      acac_e           0.00751   [0.00931, -0.012]
    uri_e            0.012     [0.012, -0.0538]      Rtotal_e         0.00363   [0, 0.217]
    glu_L_e          0.012     [0.012, -0.101]       h2o_e            0.00192   [1.73, -3]
    akg_e            0.012     [0.012, -0.15]        gthox_e          0.00179   [0, 0.0565]
    ser_L_e          0.012     [0.012, -0.387]       lac_D_e          0.001     [0, 0.001]
    nh4_e            0.012     [0.012, -1.21]        co2_e            0.000449  [0, 1.35]
    gly_e            0.0119    [0.03, -0.312]        for_e            0.000449  [-0.012, 0.04]
    bhb_e            0.0102    [-0.00931, 0.012]     amp_e            0         [-0.012, 0.833]
    pro_L_e          0.00907   [0.00907, 0.00907]    abt_e            0         [0, 0.552]
    gln_L_e          0.00796   [0.09, -0.202]        fe3_e            0         [-0.012, 0.545]
    arg_L_e          0.0079    [0.0079, 0.0079]      xylt_e           0         [-0.012, 0.54]
    val_L_e          0.00776   [0.00776, 0.00776]    hdcea_e          0         [0, 0.427]
    thr_L_e          0.00688   [0.00688, 0.00688]    man_e            0         [-0.012, 0.422]
    ile_L_e          0.00629   [0.00629, 0.012]      vacc_e           0         [0, 0.373]
    asn_L_e          0.00615   [0.00615, 0.00615]    ade_e            0         [0, 0.356]
    phe_L_e          0.00571   [0.00571, 0.012]      prpp_e           0         [0, 0.356]
    lpchol_hs_e      0.00363   [0, 0.0119]           aicar_e          0         [0, 0.342]
    tyr_L_e          0.00351   [0.00351, 0.012]      5oxpro_e         0         [0, 0.309]
    met_L_e          0.00337   [0.00337, 0.012]      elaid_e          0         [0, 0.304]
    his_L_e          0.00278   [0, 0.012]            fum_e            0         [-0.012, 0.289]
    pe_hs_e          0.00137   [0.012, -0.0189]      ocdcea_e         0         [0, 0.219]
    biomass_other_c  0.00119   [0.00119, 0.00119]    gal_e            0         [0, 0.171]
    cys_L_e          0.00102   [0.00102, 0.012]      nrvnc_e          0         [0, 0.111]
    cytd_e           0.000859  [0.012, -0.0538]      dlnlcg_e         0         [0, 0.072]
    inost_e          0.000513  [0.000513, 0.000513]  lnlncg_e         0         [0, 0.072]
    arachd_e         0.000449  [0, 0.012]            pheme_e          0         [0, 0.0307]
    hdca_e           0.000385  [-0.0116, 0.3]        bilglcur_e       0         [0, 0.0291]
    sph1p_e          0.000385  [0.000385, 0.012]     co_e             0         [0, 0.0291]
    pglyc_hs_e       0.000321  [0.000321, 0.012]     bildglcur_e      0         [0, 0.0283]
    trp_L_e          0.000293  [0.000293, 0.000293]  fald_e           0         [-0.012, 0.024]
    datp_n           0.00029   [0.012, -0.372]       meoh_e           0         [-0.012, 0.024]
    dttp_n           0.000288  [0.000288, -0.0117]   pe_hs_r          0         [-0.012, 0.0189]
    dgtp_n           0.000218  [0.000218, 0.000218]  cholate_e        0         [0, 0.0156]
    fe2_e            0         [-0.012, 0.545]       3mlda_e          0         [0, 0.0126]
    lnlc_e           0         [-0.012, 0.06]        34dhoxpeg_e      0         [0, 0.012]
    3hpvs_e          0         [0, 0.012]            3ddcrn_e         0         [0, 0.012]
    4mop_e           0         [0, 0.012]            3deccrn_e        0         [0, 0.012]
    arab_L_e         0         [0, 0.012]            3hdececrn_e      0         [0, 0.012]
    carn_e           0         [0, 0.012]            3hexdcrn_e       0         [0, 0.012]
    cit_e            0         [0, 0.012]            3hpvstet_e       0         [0, 0.012]
    crvnc_e          0         [0, 0.012]            3ivcrn_e         0         [0, 0.012]
    csa_e            0         [0, 0.012]            3octdec2crn_e    0         [0, 0.012]
    dopa_e           0         [0, 0.012]            3octdeccrn_e     0         [0, 0.012]
    etoh_e           0         [0, 0.012]            3octdece1crn_e   0         [0, 0.012]
    fol_e            0         [0, 0.012]            3tdcrn_e         0         [0, 0.012]
    gluala_e         0         [0, 0.012]            3tetd7ecoacrn_e  0         [0, 0.012]
    ha_e             0         [0, 0.012]            3thexddcoacrn_e  0         [0, 0.012]
    leuktrA4_e       0         [0, 0.012]            3ttetddcoacrn_e  0         [0, 0.012]
    lnlnca_e         0         [0, 0.012]            am19cs_e         0         [0, 0.012]
    mag_hs_e         0         [0, 0.012]            am1csa_e         0         [0, 0.012]
    n2m2nmasn_e      0         [0, 0.012]            am9csa_e         0         [0, 0.012]
    o2s_e            0         [0, 0.012]            c6crn_e          0         [0, 0.012]
    octa_e           0         [0, 0.012]            c8crn_e          0         [0, 0.012]
    ptvstlac_e       0         [0, 0.012]            eicostet_e       0         [0, 0.012]
    thymd_e          0         [0, 0.012]            fol_c            0         [0, 0.012]
    gchola_e         0         [-0.00355, 0.012]     ivcrn_e          0         [0, 0.012]
    pnto_R_e         0         [0, 0.011]            leuktrB4_e       0         [0, 0.012]
    ahcys_e          0         [0, -0.00863]         leuktrD4_e       0         [0, 0.012]
    adrn_e           0         [0, -0.012]           n5m2masn_g       0         [0, 0.012]
    clpnd_e          0         [0, -0.024]           ptvstm3_e        0         [0, 0.012]
                                                     HC01609_e        0         [-0.012, 0.012]
                                                     HC01610_e        0         [-0.012, 0.012]
                                                     1glyc_hs_e       0         [0, 0.0117]
                                                     ethamp_r         0         [0, 0.0116]
                                                     fuc14galacgl...  0         [0, 0.0116]
                                                     octdececoa_c     0         [0, 0.011]
                                                     so4_e            0         [0, 0.011]
                                                     taur_c           0         [0, 0.011]
                                                     taur_e           0         [0, 0.011]
                                                     4hphac_e         0         [0, 0.00849]
                                                     pchol_hs_e       0         [0, 0.00809]
                                                     pheacgln_e       0         [0, 0.00629]
                                                     HC00250_e        0         [0, 0.00549]
                                                     4mptnl_e         0         [0, 0.00355]
                                                     5adtststerone_e  0         [0, 0.00355]
                                                     andrstrn_e       0         [0, 0.00355]
                                                     andrstrnglc_e    0         [0, 0.00355]
                                                     aprgstrn_e       0         [0, 0.00355]
                                                     chsterol_e       0         [0, 0.00355]
                                                     ocdca_e          0         [0, -0.07]



```python
model_NormalBiopsy2.summary(fva= True)
```

    IN FLUXES                                        OUT FLUXES                                    OBJECTIVES
    -----------------------------------------------  --------------------------------------------  ----------------------
    id                   Flux  Range                 id                   Flux  Range              biomass_reac...  0.022
    ---------------  --------  --------------------  ---------------  --------  -----------------
    glc_D_e          0.139     [0.00605, 1]          h_e              0.34      [-0.872, 5.4]
    o2_e             0.0726    [-0.0107, 1]          lac_L_e          0.161     [-0.09, 2.81]
    lys_L_e          0.013     [0.013, 0.013]        pyr_e            0.15      [-0.012, 2.74]
    leu_L_e          0.012     [0.012, 0.012]        h2o_e            0.106     [3.23, -3.59]
    gam_e            0.012     [0, 0.012]            nh4_e            0.0205    [-0.012, 1.3]
    o2s_e            0.012     [0, 0.012]            cit_e            0.012     [-0.012, 0.385]
    acac_e           0.012     [-0.00931, 0.012]     ac_e             0.011     [0, 2.03]
    akg_e            0.012     [0.012, -0.364]       glyc_S_e         0.00715   [0, 1.33]
    ser_L_e          0.012     [0.012, -0.387]       bhb_e            0.00662   [0.00931, -0.012]
    man_e            0.012     [0.012, -0.691]       oxa_e            0.00645   [0, 0.62]
    gly_e            0.0119    [0.03, -0.369]        tmndnc_e         0.00539   [0, 0.024]
    glu_L_e          0.0109    [0.012, -0.101]       leuktrD4_e       0.00494   [0, 0.012]
    pro_L_e          0.00907   [0.00907, 0.00907]    3bcrn_e          0.00491   [0, 0.012]
    arg_L_e          0.0079    [0.0079, 0.0079]      5oxpro_e         0.00397   [0, 0.392]
    val_L_e          0.00776   [0.00776, 0.00776]    dcmp_e           0.00341   [-0.012, 0.0538]
    crn_e            0.0076    [0, 0.012]            Rtotal_e         0.00325   [0, 0.395]
    gln_L_e          0.00717   [0.09, -0.282]        abt_e            0.00284   [0, 0.774]
    thr_L_e          0.00688   [0.00688, 0.00688]    3ivcrn_e         0.00269   [0, 0.012]
    4hpro_LT_m       0.00645   [0, 0.012]            lac_D_e          0.001     [0, 0.001]
    ile_L_e          0.00629   [0.00629, 0.012]      for_e            0.000449  [-0.012, 0.04]
    asn_L_e          0.00615   [0.00615, 0.00615]    co2_e            0         [0, 1.79]
    phe_L_e          0.00571   [0.00571, 0.012]      pi_e             0         [0, 1.76]
    crvnc_e          0.00539   [0, 0.012]            dgtp_m           0         [0, 0.856]
    gthrd_e          0.00494   [-0.023, 0.09]        datp_m           0         [0, 0.833]
    leuktrA4_e       0.00494   [0, 0.012]            amp_e            0         [-0.012, 0.833]
    gluala_e         0.00397   [0, 0.012]            xylt_e           0         [-0.012, 0.762]
    ps_hs_e          0.00391   [0.00391, 0.012]      aicar_e          0         [0, 0.614]
    lpchol_hs_e      0.00363   [0, 0.0119]           ppi_e            0         [0, 0.586]
    dctp_n           0.00362   [0.012, -0.0538]      fe3_e            0         [-0.012, 0.545]
    tyr_L_e          0.00351   [0.00351, 0.012]      elaid_e          0         [0, 0.428]
    met_L_e          0.00337   [0.00337, 0.012]      hdcea_e          0         [0, 0.428]
    his_L_e          0.00278   [0, 0.012]            ocdcea_e         0         [0, 0.428]
    ha_e             0.00142   [0, 0.012]            fum_e            0         [-0.012, 0.399]
    pe_hs_e          0.00137   [0.012, -0.0189]      mal_L_e          0         [-0.012, 0.399]
    biomass_other_c  0.00119   [0.00119, 0.00119]    gal_e            0         [0, 0.363]
    atp_e            0.00118   [-0.0325, 0.8]        ala_L_e          0         [-0.04, 0.359]
    utp_e            0.00118   [0.02, -0.0458]       ade_e            0         [0, 0.356]
    cys_L_e          0.00102   [0.00102, 0.012]      prpp_e           0         [0, 0.356]
    cytd_e           0.000859  [0.012, -0.0538]      nrvnc_e          0         [0, 0.272]
    xmp_e            0.000794  [0.012, -0.832]       cgly_e           0         [-0.012, 0.101]
    inost_e          0.000513  [0.000513, 0.000513]  dlnlcg_e         0         [0, 0.072]
    asp_L_e          0.000413  [0, 0.012]            lnlncg_e         0         [0, 0.072]
    sph1p_e          0.000385  [0.000385, 0.012]     gthox_e          0         [0, 0.0565]
    pglyc_hs_e       0.000321  [0.000321, 0.012]     uri_e            0         [-0.012, 0.0538]
    trp_L_e          0.000293  [0.000293, 0.000293]  pheme_e          0         [0, 0.038]
    datp_n           0.00029   [0.012, -0.821]       bilglcur_e       0         [0, 0.0366]
    dttp_n           0.000288  [0.000288, -0.0117]   co_e             0         [0, 0.0366]
    dag_hs_e         0.000235  [0.000385, -0.0235]   bildglcur_e      0         [0, 0.0359]
    dgtp_n           0.000218  [0.012, -0.844]       clpnd_e          0         [0, 0.024]
    fe2_e            0         [-0.012, 0.545]       fald_e           0         [-0.012, 0.024]
    hdca_e           0         [-0.0116, 0.3]        meoh_e           0         [-0.012, 0.024]
    glyc_e           0         [-0.012, 0.21]        pe_hs_r          0         [-0.012, 0.0189]
    lnlc_e           0         [-0.012, 0.06]        cholate_e        0         [0, 0.0156]
    3mob_e           0         [0, 0.012]            3mlda_e          0         [0, 0.0126]
    4mop_e           0         [0, 0.012]            34dhoxpeg_e      0         [0, 0.012]
    arab_L_e         0         [0, 0.012]            3ddcrn_e         0         [0, 0.012]
    arachd_e         0         [0, 0.012]            3deccrn_e        0         [0, 0.012]
    carn_e           0         [0, 0.012]            3hdececrn_e      0         [0, 0.012]
    csa_e            0         [0, 0.012]            3hexdcrn_e       0         [0, 0.012]
    dopa_e           0         [0, 0.012]            3hpvstet_e       0         [0, 0.012]
    etoh_e           0         [0, 0.012]            3octdec2crn_e    0         [0, 0.012]
    fol_e            0         [0, 0.012]            3octdeccrn_e     0         [0, 0.012]
    lnlnca_e         0         [0, 0.012]            3octdece1crn_e   0         [0, 0.012]
    mag_hs_e         0         [0, 0.012]            3tdcrn_e         0         [0, 0.012]
    n2m2nmasn_e      0         [0, 0.012]            3tetd7ecoacrn_e  0         [0, 0.012]
    octa_e           0         [0, 0.012]            3thexddcoacrn_e  0         [0, 0.012]
    ptvstlac_e       0         [0, 0.012]            3ttetddcoacrn_e  0         [0, 0.012]
    retinol_e        0         [0, 0.012]            adrn_e           0         [0, 0.012]
    thymd_e          0         [0, 0.012]            am19cs_e         0         [0, 0.012]
    gchola_e         0         [-0.00355, 0.012]     am1csa_e         0         [0, 0.012]
    glyc3p_e         0         [-0.0119, 0.012]      am9csa_e         0         [0, 0.012]
    pnto_R_e         0         [0, 0.011]            c6crn_e          0         [0, 0.012]
    octdececoa_c     0         [0, -0.011]           c8crn_e          0         [0, 0.012]
    vacc_e           0         [0, -0.428]           eicostet_e       0         [0, 0.012]
                                                     fol_c            0         [0, 0.012]
                                                     ivcrn_e          0         [0, 0.012]
                                                     leuktrB4_e       0         [0, 0.012]
                                                     n5m2masn_g       0         [0, 0.012]
                                                     ptvstm3_e        0         [0, 0.012]
                                                     retn_e           0         [0, 0.012]
                                                     HC01609_e        0         [-0.012, 0.012]
                                                     HC01610_e        0         [-0.012, 0.012]
                                                     1glyc_hs_e       0         [0, 0.0117]
                                                     ethamp_r         0         [0, 0.0116]
                                                     fuc14galacgl...  0         [0, 0.0116]
                                                     so4_e            0         [0, 0.011]
                                                     taur_c           0         [0, 0.011]
                                                     taur_e           0         [0, 0.011]
                                                     ahcys_e          0         [0, 0.00863]
                                                     4hphac_e         0         [0, 0.00849]
                                                     pchol_hs_e       0         [0, 0.00809]
                                                     pheacgln_e       0         [0, 0.00629]
                                                     HC00250_e        0         [0, 0.00549]
                                                     4mptnl_e         0         [0, 0.00355]
                                                     5adtststerone_e  0         [0, 0.00355]
                                                     andrstrn_e       0         [0, 0.00355]
                                                     andrstrnglc_e    0         [0, 0.00355]
                                                     aprgstrn_e       0         [0, 0.00355]
                                                     chsterol_e       0         [0, 0.00355]
                                                     3hpvs_e          0         [0, -0.012]
                                                     ocdca_e          0         [0, -0.07]



```python
model_NormalBiopsy3.summary(fva=True)
```

    IN FLUXES                                        OUT FLUXES                                      OBJECTIVES
    -----------------------------------------------  ----------------------------------------------  ----------------------
    id                   Flux  Range                 id                   Flux  Range                biomass_reac...  0.022
    ---------------  --------  --------------------  ---------------  --------  -------------------
    glc_D_e          0.153     [0.00605, 1]          h_e              0.379     [-0.977, 5.4]
    o2_e             0.102     [-0.0149, 1]          lac_L_e          0.272     [-0.09, 2.81]
    lys_L_e          0.013     [0.013, 0.013]        h2o_e            0.0766    [3.27, -3.59]
    leu_L_e          0.012     [0.012, 0.012]        pyr_e            0.0491    [-0.012, 2.75]
    gam_e            0.012     [0, 0.012]            cit_e            0.0475    [-0.012, 0.393]
    o2s_e            0.012     [0, 0.012]            nh4_e            0.0112    [-0.012, 1.3]
    glyc3p_e         0.012     [-0.0119, 0.012]      3bcrn_e          0.0112    [0, 0.012]
    glu_L_e          0.012     [0.012, -0.101]       glyc_S_e         0.0111    [0, 1.33]
    akg_e            0.012     [0.012, -0.371]       dcmp_e           0.00915   [-0.012, 0.0538]
    ser_L_e          0.012     [0.012, -0.387]       cgly_e           0.00424   [-0.012, 0.101]
    fum_e            0.012     [0.012, -0.407]       dag_hs_e         0.00374   [-0.000385, 0.0235]
    mal_L_e          0.012     [0.012, -0.407]       Rtotal_e         0.00363   [0, 0.395]
    gly_e            0.0119    [0.03, -0.369]        tmndnc_e         0.00269   [0, 0.096]
    3mob_e           0.0115    [0, 0.012]            lac_D_e          0.001     [0, 0.001]
    crn_e            0.0112    [0, 0.012]            for_e            0.000449  [-0.012, 0.04]
    cytd_e           0.01      [0.012, -0.0538]      ac_e             0         [-0.012, 2.03]
    pro_L_e          0.00907   [0.00907, 0.00907]    co2_e            0         [0, 1.79]
    arg_L_e          0.0079    [0.0079, 0.0079]      pi_e             0         [0, 1.76]
    ps_hs_e          0.00788   [0.00391, 0.012]      dgtp_m           0         [0, 0.856]
    val_L_e          0.00776   [0.00776, 0.00776]    datp_m           0         [0, 0.833]
    gln_L_e          0.00717   [0.09, -0.288]        amp_e            0         [-0.012, 0.833]
    thr_L_e          0.00688   [0.00688, 0.00688]    abt_e            0         [0, 0.774]
    ile_L_e          0.00629   [0.00629, 0.012]      xylt_e           0         [-0.012, 0.762]
    asn_L_e          0.00615   [0.00615, 0.00615]    man_e            0         [-0.012, 0.691]
    ocdca_e          0.00586   [0, 0.07]             oxa_e            0         [0, 0.623]
    phe_L_e          0.00571   [0.00571, 0.012]      aicar_e          0         [0, 0.614]
    gthrd_e          0.00424   [-0.023, 0.09]        ppi_e            0         [0, 0.586]
    lpchol_hs_e      0.00363   [0, 0.0119]           fe3_e            0         [-0.012, 0.545]
    tyr_L_e          0.00351   [0.00351, 0.012]      elaid_e          0         [0, 0.428]
    met_L_e          0.00337   [0.00337, 0.012]      hdcea_e          0         [0, 0.428]
    his_L_e          0.00278   [0, 0.012]            ocdcea_e         0         [0, 0.428]
    crvnc_e          0.00269   [0, 0.012]            vacc_e           0         [0, 0.428]
    acac_e           0.00269   [-0.00931, 0.012]     5oxpro_e         0         [0, 0.398]
    pe_hs_e          0.00137   [0.012, -0.0189]      gal_e            0         [0, 0.363]
    biomass_other_c  0.00119   [0.00119, 0.00119]    ala_L_e          0         [-0.04, 0.359]
    atp_e            0.00118   [-0.0325, 0.8]        ade_e            0         [0, 0.356]
    utp_e            0.00118   [0.02, -0.0458]       prpp_e           0         [0, 0.356]
    cys_L_e          0.00102   [0.00102, 0.012]      nrvnc_e          0         [0, 0.272]
    xmp_e            0.000794  [0.012, -0.832]       clpnd_e          0         [0, 0.096]
    inost_e          0.000513  [0.000513, 0.000513]  eicostet_e       0         [0, 0.084]
    hdca_e           0.000385  [-0.0116, 0.3]        dlnlcg_e         0         [0, 0.072]
    sph1p_e          0.000385  [0.000385, 0.012]     lnlncg_e         0         [0, 0.072]
    pglyc_hs_e       0.000321  [0.000321, 0.012]     lnlnca_e         0         [-0.012, 0.072]
    trp_L_e          0.000293  [0.000293, 0.000293]  gthox_e          0         [0, 0.0565]
    datp_n           0.00029   [0.012, -0.821]       uri_e            0         [-0.012, 0.0538]
    dttp_n           0.000288  [0.000288, -0.0117]   pheme_e          0         [0, 0.0429]
    dgtp_n           0.000218  [0.012, -0.844]       bilglcur_e       0         [0, 0.041]
    dctp_n           0.000208  [0.012, -0.0538]      co_e             0         [0, 0.041]
    fe2_e            0         [-0.012, 0.545]       bildglcur_e      0         [0, 0.04]
    glyc_e           0         [-0.012, 0.21]        fald_e           0         [-0.012, 0.024]
    lnlc_e           0         [-0.012, 0.06]        meoh_e           0         [-0.012, 0.024]
    3hpvs_e          0         [0, 0.012]            pe_hs_r          0         [-0.012, 0.0189]
    4hpro_LT_m       0         [0, 0.012]            cholate_e        0         [0, 0.0156]
    4mop_e           0         [0, 0.012]            3mlda_e          0         [0, 0.0126]
    arab_L_e         0         [0, 0.012]            34dhoxpeg_e      0         [0, 0.012]
    arachd_e         0         [0, 0.012]            3ddcrn_e         0         [0, 0.012]
    asp_L_e          0         [0, 0.012]            3deccrn_e        0         [0, 0.012]
    carn_e           0         [0, 0.012]            3hdececrn_e      0         [0, 0.012]
    csa_e            0         [0, 0.012]            3hexdcrn_e       0         [0, 0.012]
    dopa_e           0         [0, 0.012]            3hpvstet_e       0         [0, 0.012]
    etoh_e           0         [0, 0.012]            3ivcrn_e         0         [0, 0.012]
    fol_e            0         [0, 0.012]            3octdec2crn_e    0         [0, 0.012]
    gluala_e         0         [0, 0.012]            3octdeccrn_e     0         [0, 0.012]
    ha_e             0         [0, 0.012]            3octdece1crn_e   0         [0, 0.012]
    leuktrA4_e       0         [0, 0.012]            3tdcrn_e         0         [0, 0.012]
    mag_hs_e         0         [0, 0.012]            3tetd7ecoacrn_e  0         [0, 0.012]
    n2m2nmasn_e      0         [0, 0.012]            3thexddcoacrn_e  0         [0, 0.012]
    octa_e           0         [0, 0.012]            3ttetddcoacrn_e  0         [0, 0.012]
    ptvstlac_e       0         [0, 0.012]            adrn_e           0         [0, 0.012]
    retinol_e        0         [0, 0.012]            am19cs_e         0         [0, 0.012]
    thymd_e          0         [0, 0.012]            am1csa_e         0         [0, 0.012]
    gchola_e         0         [-0.00355, 0.012]     am9csa_e         0         [0, 0.012]
    ppa_e            0         [-0.00926, 0.012]     c6crn_e          0         [0, 0.012]
    bhb_e            0         [-0.00931, 0.012]     c8crn_e          0         [0, 0.012]
    pnto_R_e         0         [0, 0.011]            fol_c            0         [0, 0.012]
    octdececoa_c     0         [0, -0.011]           ivcrn_e          0         [0, 0.012]
                                                     leuktrB4_e       0         [0, 0.012]
                                                     leuktrD4_e       0         [0, 0.012]
                                                     n5m2masn_g       0         [0, 0.012]
                                                     ptvstm3_e        0         [0, 0.012]
                                                     retn_e           0         [0, 0.012]
                                                     HC01609_e        0         [-0.012, 0.012]
                                                     HC01610_e        0         [-0.012, 0.012]
                                                     1glyc_hs_e       0         [0, 0.0117]
                                                     ethamp_r         0         [0, 0.0116]
                                                     fuc14galacgl...  0         [0, 0.0116]
                                                     so4_e            0         [0, 0.011]
                                                     taur_c           0         [0, 0.011]
                                                     taur_e           0         [0, 0.011]
                                                     ahcys_e          0         [0, 0.00863]
                                                     4hphac_e         0         [0, 0.00849]
                                                     pchol_hs_e       0         [0, 0.00809]
                                                     pheacgln_e       0         [0, 0.00629]
                                                     HC00250_e        0         [0, 0.00549]
                                                     4mptnl_e         0         [0, 0.00355]
                                                     5adtststerone_e  0         [0, 0.00355]
                                                     andrstrn_e       0         [0, 0.00355]
                                                     andrstrnglc_e    0         [0, 0.00355]
                                                     aprgstrn_e       0         [0, 0.00355]
                                                     chsterol_e       0         [0, 0.00355]



```python

model_NormalBiopsy4.summary(fva=True)
```

    IN FLUXES                                        OUT FLUXES                                      OBJECTIVES
    -----------------------------------------------  ----------------------------------------------  ----------------------
    id                   Flux  Range                 id                   Flux  Range                biomass_reac...  0.022
    ---------------  --------  --------------------  ---------------  --------  -------------------
    o2_e             0.175     [0.16, 1]             h_e              0.364     [-0.844, 3.47]
    glc_D_e          0.123     [0.00605, 0.998]      cit_e            0.135     [-0.012, 0.297]
    gthrd_e          0.09      [-0.023, 0.09]        glyc_S_e         0.114     [0, 0.866]
    lac_L_e          0.0618    [-0.072, 0.09]        ala_L_e          0.103     [-0.04, 0.359]
    lys_L_e          0.013     [0.013, 0.013]        cgly_e           0.09      [-0.012, 0.101]
    gly_e            0.013     [0.03, -0.369]        ac_e             0.0478    [-0.012, 0.15]
    leu_L_e          0.012     [0.012, 0.012]        h2o_e            0.0452    [1.69, -3.13]
    ps_hs_e          0.012     [0.00391, 0.012]      pyr_e            0.0283    [-0.012, 0.15]
    asp_L_e          0.012     [0, 0.012]            nh4_e            0.012     [-0.012, 1.3]
    crn_e            0.012     [0, 0.012]            3bcrn_e          0.012     [0, 0.012]
    gam_e            0.012     [0, 0.012]            Rtotal_e         0.0116    [0, 0.243]
    lpchol_hs_e      0.012     [0, 0.012]            bhb_e            0.00931   [0.00931, -0.012]
    o2s_e            0.012     [0, 0.012]            dctp_m           0.00885   [0, 0.0658]
    acac_e           0.012     [-0.00931, 0.012]     glyb_e           0.00837   [0, 0.0124]
    glyc3p_e         0.012     [0.012, -0.012]       dag_hs_e         0.00785   [-0.000385, 0.0359]
    glu_L_e          0.012     [0.012, -0.101]       tmndnc_e         0.00269   [0, 0.096]
    akg_e            0.012     [0.012, -0.15]        lac_D_e          0.001     [0, 0.001]
    fum_e            0.012     [0.012, -0.255]       for_e            0.000449  [-0.012, 0.04]
    mal_L_e          0.012     [0.012, -0.255]       pheme_e          0.000142  [0, 0.032]
    ser_L_e          0.012     [0.012, -0.387]       pi_e             0         [0, 1.77]
    pro_L_e          0.00907   [0.00907, 0.00907]    co2_e            0         [0, 1.23]
    dcmp_e           0.00906   [0.012, -0.0538]      amp_e            0         [-0.012, 0.833]
    gln_L_e          0.00796   [0.09, -0.192]        abt_e            0         [0, 0.634]
    arg_L_e          0.0079    [0.0079, 0.0079]      xylt_e           0         [-0.012, 0.622]
    val_L_e          0.00776   [0.00776, 0.00776]    ppi_e            0         [0, 0.59]
    thr_L_e          0.00688   [0.00688, 0.00688]    fe3_e            0         [-0.012, 0.545]
    ile_L_e          0.00629   [0.00629, 0.012]      oxa_e            0         [0, 0.511]
    asn_L_e          0.00615   [0.00615, 0.00615]    man_e            0         [-0.012, 0.485]
    phe_L_e          0.00571   [0.00571, 0.012]      hdcea_e          0         [0, 0.429]
    tyr_L_e          0.00351   [0.00351, 0.012]      dgtp_m           0         [0, 0.42]
    met_L_e          0.00337   [0.00337, 0.012]      vacc_e           0         [0, 0.417]
    his_L_e          0.00278   [0, 0.012]            ade_e            0         [0, 0.359]
    crvnc_e          0.00269   [0, 0.012]            prpp_e           0         [0, 0.359]
    pe_hs_e          0.00137   [0.012, -0.0189]      elaid_e          0         [0, 0.346]
    biomass_other_c  0.00119   [0.00119, 0.00119]    aicar_e          0         [0, 0.345]
    atp_e            0.00118   [-0.0325, 0.8]        5oxpro_e         0         [0, 0.299]
    uri_e            0.00118   [0.012, -0.0538]      ocdcea_e         0         [0, 0.26]
    3mob_e           0.00114   [0, 0.012]            gal_e            0         [0, 0.196]
    cys_L_e          0.00102   [0.00102, 0.012]      nrvnc_e          0         [0, 0.138]
    cytd_e           0.000859  [0.012, -0.0538]      clpnd_e          0         [0, 0.096]
    xmp_e            0.000794  [0.012, -0.832]       eicostet_e       0         [0, 0.084]
    inost_e          0.000513  [0.000513, 0.000513]  dlnlcg_e         0         [0, 0.072]
    sph1p_e          0.000385  [0.000385, 0.012]     lnlncg_e         0         [0, 0.072]
    pglyc_hs_e       0.000321  [0.000321, 0.012]     lnlnca_e         0         [-0.012, 0.072]
    trp_L_e          0.000293  [0.000293, 0.000293]  gthox_e          0         [0, 0.0565]
    datp_n           0.00029   [0.012, -0.41]        dctp_n           0         [-0.012, 0.0538]
    dttp_n           0.000288  [0.000288, -0.0117]   utp_e            0         [-0.02, 0.0458]
    dgtp_n           0.000218  [0.012, -0.477]       bilglcur_e       0         [0, 0.0309]
    fe2_e            0.000142  [-0.012, 0.545]       co_e             0         [0, 0.0309]
    hdca_e           0         [-0.0116, 0.3]        bildglcur_e      0         [0, 0.0305]
    glyc_e           0         [-0.012, 0.21]        fald_e           0         [-0.012, 0.024]
    ocdca_e          0         [0, 0.07]             meoh_e           0         [-0.012, 0.024]
    lnlc_e           0         [-0.012, 0.06]        pe_hs_r          0         [-0.012, 0.0189]
    3hpvs_e          0         [0, 0.012]            cholate_e        0         [0, 0.0156]
    4hpro_LT_m       0         [0, 0.012]            3mlda_e          0         [0, 0.0126]
    4mop_e           0         [0, 0.012]            34dhoxpeg_e      0         [0, 0.012]
    arab_L_e         0         [0, 0.012]            3ddcrn_e         0         [0, 0.012]
    arachd_e         0         [0, 0.012]            3deccrn_e        0         [0, 0.012]
    carn_e           0         [0, 0.012]            3hdececrn_e      0         [0, 0.012]
    csa_e            0         [0, 0.012]            3hexdcrn_e       0         [0, 0.012]
    dopa_e           0         [0, 0.012]            3hpvstet_e       0         [0, 0.012]
    etoh_e           0         [0, 0.012]            3ivcrn_e         0         [0, 0.012]
    fol_e            0         [0, 0.012]            3octdec2crn_e    0         [0, 0.012]
    gluala_e         0         [0, 0.012]            3octdeccrn_e     0         [0, 0.012]
    ha_e             0         [0, 0.012]            3octdece1crn_e   0         [0, 0.012]
    leuktrA4_e       0         [0, 0.012]            3tdcrn_e         0         [0, 0.012]
    mag_hs_e         0         [0, 0.012]            3tetd7ecoacrn_e  0         [0, 0.012]
    n2m2nmasn_e      0         [0, 0.012]            3thexddcoacrn_e  0         [0, 0.012]
    octa_e           0         [0, 0.012]            3ttetddcoacrn_e  0         [0, 0.012]
    ptvstlac_e       0         [0, 0.012]            adrn_e           0         [0, 0.012]
    retinol_e        0         [0, 0.012]            am19cs_e         0         [0, 0.012]
    thymd_e          0         [0, 0.012]            am1csa_e         0         [0, 0.012]
    gchola_e         0         [-0.00355, 0.012]     am9csa_e         0         [0, 0.012]
    pnto_R_e         0         [0, 0.011]            c6crn_e          0         [0, 0.012]
    ahcys_e          0         [0, -0.00863]         c8crn_e          0         [0, 0.012]
                                                     fol_c            0         [0, 0.012]
                                                     ivcrn_e          0         [0, 0.012]
                                                     leuktrB4_e       0         [0, 0.012]
                                                     leuktrD4_e       0         [0, 0.012]
                                                     n5m2masn_g       0         [0, 0.012]
                                                     ptvstm3_e        0         [0, 0.012]
                                                     retn_e           0         [0, 0.012]
                                                     HC01609_e        0         [-0.012, 0.012]
                                                     HC01610_e        0         [-0.012, 0.012]
                                                     1glyc_hs_e       0         [0, 0.0117]
                                                     ethamp_r         0         [0, 0.0116]
                                                     fuc14galacgl...  0         [0, 0.0116]
                                                     octdececoa_c     0         [0, 0.011]
                                                     so4_e            0         [0, 0.011]
                                                     taur_c           0         [0, 0.011]
                                                     taur_e           0         [0, 0.011]
                                                     4hphac_e         0         [0, 0.00849]
                                                     pchol_hs_e       0         [0, 0.00809]
                                                     pheacgln_e       0         [0, 0.00629]
                                                     HC00250_e        0         [0, 0.00549]
                                                     4mptnl_e         0         [0, 0.00355]
                                                     5adtststerone_e  0         [0, 0.00355]
                                                     andrstrn_e       0         [0, 0.00355]
                                                     andrstrnglc_e    0         [0, 0.00355]
                                                     aprgstrn_e       0         [0, 0.00355]
                                                     chsterol_e       0         [0, 0.00355]



```python

model_NormalBiopsy5.summary(fva=True)
```

    IN FLUXES                                        OUT FLUXES                                      OBJECTIVES
    -----------------------------------------------  ----------------------------------------------  ----------------------
    id                   Flux  Range                 id                   Flux  Range                biomass_reac...  0.022
    ---------------  --------  --------------------  ---------------  --------  -------------------
    glc_D_e          0.131     [0.00605, 1]          h_e              0.317     [-0.977, 5.41]
    o2_e             0.0844    [-0.0149, 1]          lac_L_e          0.162     [-0.09, 2.81]
    ala_L_e          0.0158    [0.04, -0.359]        pyr_e            0.16      [-0.012, 2.75]
    lys_L_e          0.013     [0.013, 0.013]        h2o_e            0.148     [3.26, -3.59]
    leu_L_e          0.012     [0.012, 0.012]        bhb_e            0.00931   [0.00931, -0.012]
    gam_e            0.012     [0, 0.012]            cit_e            0.00748   [-0.012, 0.393]
    o2s_e            0.012     [0, 0.012]            leuktrD4_e       0.00397   [0, 0.012]
    acac_e           0.012     [-0.00931, 0.012]     Rtotal_e         0.00363   [0, 0.396]
    glyc3p_e         0.012     [0.012, -0.012]       dctp_m           0.00307   [0, 0.0658]
    akg_e            0.012     [0.012, -0.371]       tmndnc_e         0.00269   [0, 0.096]
    man_e            0.012     [0.012, -0.691]       dag_hs_e         0.000729  [-0.000385, 0.0359]
    pro_L_e          0.00907   [0.00907, 0.00907]    for_e            0.000449  [-0.012, 0.04]
    arg_L_e          0.0079    [0.0079, 0.0079]      co2_e            0         [0, 1.79]
    val_L_e          0.00776   [0.00776, 0.00776]    pi_e             0         [0, 1.77]
    gly_e            0.0072    [0.03, -0.369]        glyc_S_e         0         [0, 1.33]
    gln_L_e          0.00717   [0.09, -0.288]        dgtp_m           0         [0, 0.856]
    thr_L_e          0.00688   [0.00688, 0.00688]    amp_e            0         [-0.012, 0.833]
    ile_L_e          0.00629   [0.00629, 0.012]      abt_e            0         [0, 0.774]
    asn_L_e          0.00615   [0.00615, 0.00615]    xylt_e           0         [-0.012, 0.762]
    phe_L_e          0.00571   [0.00571, 0.012]      oxa_e            0         [0, 0.62]
    ps_hs_e          0.00487   [0.00391, 0.012]      aicar_e          0         [0, 0.614]
    4hpro_LT_m       0.00465   [0, 0.012]            ppi_e            0         [0, 0.59]
    3mob_e           0.00434   [0, 0.012]            fe3_e            0         [-0.012, 0.545]
    gthrd_e          0.00397   [-0.023, 0.09]        elaid_e          0         [0, 0.429]
    leuktrA4_e       0.00397   [0, 0.012]            hdcea_e          0         [0, 0.429]
    lpchol_hs_e      0.00363   [0, 0.012]            ocdcea_e         0         [0, 0.429]
    tyr_L_e          0.00351   [0.00351, 0.012]      vacc_e           0         [0, 0.429]
    fum_e            0.00341   [0.012, -0.407]       5oxpro_e         0         [0, 0.398]
    met_L_e          0.00337   [0.00337, 0.012]      gal_e            0         [0, 0.363]
    dcmp_e           0.00307   [0.012, -0.0538]      ade_e            0         [0, 0.359]
    glyc_e           0.00283   [-0.012, 0.21]        prpp_e           0         [0, 0.359]
    his_L_e          0.00278   [0, 0.012]            nrvnc_e          0         [0, 0.272]
    crvnc_e          0.00269   [0, 0.012]            cgly_e           0         [-0.012, 0.101]
    pe_hs_e          0.00137   [0.012, -0.0189]      glu_L_e          0         [-0.012, 0.101]
    biomass_other_c  0.00119   [0.00119, 0.00119]    clpnd_e          0         [0, 0.096]
    atp_e            0.00118   [-0.0325, 0.8]        eicostet_e       0         [0, 0.084]
    utp_e            0.00118   [0.02, -0.0458]       dlnlcg_e         0         [0, 0.072]
    cys_L_e          0.00102   [0.00102, 0.012]      lnlncg_e         0         [0, 0.072]
    cytd_e           0.000859  [0.012, -0.0538]      lnlnca_e         0         [-0.012, 0.072]
    xmp_e            0.000794  [0.012, -0.832]       gthox_e          0         [0, 0.0565]
    inost_e          0.000513  [0.000513, 0.000513]  pheme_e          0         [0, 0.043]
    hdca_e           0.000385  [-0.0116, 0.3]        bilglcur_e       0         [0, 0.041]
    sph1p_e          0.000385  [0.000385, 0.012]     co_e             0         [0, 0.041]
    pglyc_hs_e       0.000321  [0.000321, 0.012]     bildglcur_e      0         [0, 0.0401]
    ser_L_e          0.000301  [0.012, -0.387]       fald_e           0         [-0.012, 0.024]
    trp_L_e          0.000293  [0.000293, 0.000293]  pe_hs_r          0         [-0.012, 0.0189]
    datp_n           0.00029   [0.012, -0.821]       cholate_e        0         [0, 0.0156]
    dttp_n           0.000288  [0.000288, -0.0117]   3mlda_e          0         [0, 0.0126]
    dgtp_n           0.000218  [0.012, -0.844]       glyb_e           0         [0, 0.0124]
    dctp_n           0.000208  [0.012, -0.0538]      34dhoxpeg_e      0         [0, 0.012]
    fe2_e            0         [-0.012, 0.545]       3bcrn_e          0         [0, 0.012]
    ocdca_e          0         [0, 0.07]             3ddcrn_e         0         [0, 0.012]
    lnlc_e           0         [-0.012, 0.06]        3deccrn_e        0         [0, 0.012]
    3hpvs_e          0         [0, 0.012]            3hdececrn_e      0         [0, 0.012]
    4mop_e           0         [0, 0.012]            3hexdcrn_e       0         [0, 0.012]
    arab_L_e         0         [0, 0.012]            3hpvstet_e       0         [0, 0.012]
    asp_L_e          0         [0, 0.012]            3ivcrn_e         0         [0, 0.012]
    carn_e           0         [0, 0.012]            3octdec2crn_e    0         [0, 0.012]
    csa_e            0         [0, 0.012]            3octdeccrn_e     0         [0, 0.012]
    etoh_e           0         [0, 0.012]            3octdece1crn_e   0         [0, 0.012]
    fol_e            0         [0, 0.012]            3tdcrn_e         0         [0, 0.012]
    gluala_e         0         [0, 0.012]            3tetd7ecoacrn_e  0         [0, 0.012]
    ha_e             0         [0, 0.012]            3thexddcoacrn_e  0         [0, 0.012]
    mag_hs_e         0         [0, 0.012]            3ttetddcoacrn_e  0         [0, 0.012]
    n2m2nmasn_e      0         [0, 0.012]            am19cs_e         0         [0, 0.012]
    octa_e           0         [0, 0.012]            am1csa_e         0         [0, 0.012]
    ptvstlac_e       0         [0, 0.012]            am9csa_e         0         [0, 0.012]
    retinol_e        0         [0, 0.012]            c6crn_e          0         [0, 0.012]
    thymd_e          0         [0, 0.012]            c8crn_e          0         [0, 0.012]
    gchola_e         0         [-0.00355, 0.012]     fol_c            0         [0, 0.012]
    ppa_e            0         [-0.00926, 0.012]     ivcrn_e          0         [0, 0.012]
    pnto_R_e         0         [0, 0.011]            leuktrB4_e       0         [0, 0.012]
    octdececoa_c     0         [0, -0.011]           n5m2masn_g       0         [0, 0.012]
    adrn_e           0         [0, -0.012]           ptvstm3_e        0         [0, 0.012]
    meoh_e           0         [0.012, -0.024]       retn_e           0         [0, 0.012]
    uri_e            0         [0.012, -0.0538]      HC01609_e        0         [-0.012, 0.012]
    mal_L_e          0         [0.012, -0.407]       HC01610_e        0         [-0.012, 0.012]
    nh4_e            0         [0.012, -1.3]         1glyc_hs_e       0         [0, 0.0117]
    ac_e             0         [0, -2.03]            ethamp_r         0         [0, 0.0116]
                                                     fuc14galacgl...  0         [0, 0.0116]
                                                     so4_e            0         [0, 0.011]
                                                     taur_c           0         [0, 0.011]
                                                     taur_e           0         [0, 0.011]
                                                     ahcys_e          0         [0, 0.00863]
                                                     4hphac_e         0         [0, 0.00849]
                                                     pchol_hs_e       0         [0, 0.00809]
                                                     pheacgln_e       0         [0, 0.00629]
                                                     HC00250_e        0         [0, 0.00549]
                                                     4mptnl_e         0         [0, 0.00355]
                                                     5adtststerone_e  0         [0, 0.00355]
                                                     andrstrn_e       0         [0, 0.00355]
                                                     andrstrnglc_e    0         [0, 0.00355]
                                                     aprgstrn_e       0         [0, 0.00355]
                                                     chsterol_e       0         [0, 0.00355]
                                                     lac_D_e          0         [0, 0.001]
                                                     arachd_e         0         [0, -0.012]
                                                     crn_e            0         [0, -0.012]
                                                     dopa_e           0         [0, -0.012]



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
from corda import CORDA
import string
import multiprocessing as mp


# define a randomized analysis function
def multiple_randomized_analysis(modelC, modelN, Recon, pos, metabolites, nrx=5, penaltyFac=1000):
    f="ByopsiesPFbootstraping/Targets_pf_"+str(penaltyFac)+"_Biopsies_n_"+str(nrx)+"_"+str(pos)+".csv" 
    print(f)    
    cConf_scores=randomize_model(modelC)
    nConf_scores=randomize_model(modelN)
    
    cConf_scores["DM_4hrpo"]=3
    cConf_scores["DM_datp_n_"]=3
    cConf_scores["DM_dctp_n_"]=3
    cConf_scores["DM_dgtp_n_"]=3
    cConf_scores["DM_dttp_n_"]=3
    cConf_scores["DM_Lcystin"]=3
    cConf_scores["DM_pe_hs_LPAREN_r_RPAREN_"]=3
    cConf_scores["EX_2hb_LPAREN_e_RPAREN_"]=3
    cConf_scores["EX_34hpp_"]=3
    cConf_scores["EX_3hpvs_LPAREN_e_RPAREN_"]=3
    cConf_scores["EX_3mob_LPAREN_e_RPAREN_"]=3
    cConf_scores["EX_4mop_LPAREN_e_RPAREN_"]=3
    cConf_scores["EX_ac_LPAREN_e_RPAREN_"]=3
    cConf_scores["EX_acac_LPAREN_e_RPAREN_"]=3
    cConf_scores["EX_acald_LPAREN_e_RPAREN_"]=3
    cConf_scores["EX_acetone_LPAREN_e_RPAREN_"]=3
    cConf_scores["EX_acgam_LPAREN_e_RPAREN_"]=3
    cConf_scores["EX_acmana_LPAREN_e_RPAREN_"]=3
    cConf_scores["EX_ade_LPAREN_e_RPAREN_"]=3
    cConf_scores["EX_adn_LPAREN_e_RPAREN_"]=3
    cConf_scores["EX_adpcbl_LPAREN_e_RPAREN_"]=3
    cConf_scores["EX_akg_LPAREN_e_RPAREN_"]=3
    cConf_scores["EX_ala_B_LPAREN_e_RPAREN_"]=3
    cConf_scores["EX_ala_D_LPAREN_e_RPAREN_"]=3
    cConf_scores["EX_ala_L_LPAREN_e_RPAREN_"]=3
    cConf_scores["EX_am9csa_LPAREN_e_RPAREN_"]=3
    cConf_scores["EX_amp_LPAREN_e_RPAREN_"]=3
    cConf_scores["EX_arab_L_LPAREN_e_RPAREN_"]=3
    cConf_scores["EX_arachd_LPAREN_e_RPAREN_"]=3
    cConf_scores["EX_arg_L_LPAREN_e_RPAREN_"]=3
    cConf_scores["EX_asn_L_LPAREN_e_RPAREN_"]=3
    cConf_scores["EX_asp_L_LPAREN_e_RPAREN_"]=3
    cConf_scores["EX_atp_LPAREN_e_RPAREN_"]=3
    cConf_scores["EX_bhb_LPAREN_e_RPAREN_"]=3
    cConf_scores["EX_btn_LPAREN_e_RPAREN_"]=3
    cConf_scores["EX_ca2_LPAREN_e_RPAREN_"]=3
    cConf_scores["EX_carn_LPAREN_e_RPAREN_"]=3
    cConf_scores["EX_caro_LPAREN_e_RPAREN_"]=3
    cConf_scores["EX_cgly_LPAREN_e_RPAREN_"]=3
    cConf_scores["EX_chol_LPAREN_e_RPAREN_"]=3
    cConf_scores["EX_chsterol_LPAREN_e_RPAREN_"]=3
    cConf_scores["EX_cit_LPAREN_e_RPAREN_"]=3
    cConf_scores["EX_cl_LPAREN_e_RPAREN_"]=3
    cConf_scores["EX_CLPND_LPAREN_e_RPAREN_"]=3
    cConf_scores["EX_co2_LPAREN_e_RPAREN_"]=3
    cConf_scores["EX_creat_LPAREN_e_RPAREN_"]=3
    cConf_scores["EX_crn_LPAREN_e_RPAREN_"]=3
    cConf_scores["EX_crvnc_LPAREN_e_RPAREN_"]=3
    cConf_scores["EX_csa_LPAREN_e_RPAREN_"]=3
    cConf_scores["EX_csn_LPAREN_e_RPAREN_"]=3
    cConf_scores["EX_cys_L_LPAREN_e_RPAREN_"]=3
    cConf_scores["EX_cytd_LPAREN_e_RPAREN_"]=3
    cConf_scores["EX_dad_2_LPAREN_e_RPAREN_"]=3
    cConf_scores["EX_dag_hs_LPAREN_e_RPAREN_"]=3
    cConf_scores["EX_dcmp_LPAREN_e_RPAREN_"]=3
    cConf_scores["EX_dcyt_LPAREN_e_RPAREN_"]=3
    cConf_scores["EX_ddca_LPAREN_e_RPAREN_"]=3
    cConf_scores["EX_dgsn_LPAREN_e_RPAREN_"]=3
    cConf_scores["EX_dhdascb_LPAREN_e_RPAREN_"]=3
    cConf_scores["EX_din_LPAREN_e_RPAREN_"]=3
    cConf_scores["EX_dopa_LPAREN_e_RPAREN_"]=3
    cConf_scores["EX_drib_LPAREN_e_RPAREN_"]=3
    cConf_scores["EX_etoh_LPAREN_e_RPAREN_"]=3
    cConf_scores["EX_fald_LPAREN_e_RPAREN_"]=3
    cConf_scores["EX_fe2_LPAREN_e_RPAREN_"]=3
    cConf_scores["EX_fe2_LPAREN_e_RPAREN_"]=3
    cConf_scores["EX_fe3_LPAREN_e_RPAREN_"]=3
    cConf_scores["EX_fmn_LPAREN_e_RPAREN_"]=3
    cConf_scores["EX_fol_LPAREN_e_RPAREN_"]=3
    cConf_scores["EX_for_LPAREN_e_RPAREN_"]=3
    cConf_scores["EX_fru_LPAREN_e_RPAREN_"]=3
    cConf_scores["EX_fuc_L_LPAREN_e_RPAREN_"]=3
    cConf_scores["EX_fum_LPAREN_e_RPAREN_"]=3
    cConf_scores["EX_gal_LPAREN_e_RPAREN_"]=3
    cConf_scores["EX_gam_LPAREN_e_RPAREN_"]=3
    cConf_scores["EX_gchola_LPAREN_e_RPAREN_"]=3
    cConf_scores["EX_glc_LPAREN_e_RPAREN_"]=3
    cConf_scores["EX_glcur_LPAREN_e_RPAREN_"]=3
    cConf_scores["EX_gln_L_LPAREN_e_RPAREN_"]=3
    cConf_scores["EX_glu_L_LPAREN_e_RPAREN_"]=3
    cConf_scores["EX_gluala_LPAREN_e_RPAREN_"]=3
    cConf_scores["EX_gly_LPAREN_e_RPAREN_"]=3
    cConf_scores["EX_glyb_LPAREN_e_RPAREN_"]=3
    cConf_scores["EX_glyc_LPAREN_e_RPAREN_"]=3
    cConf_scores["EX_glyc3p_LPAREN_e_RPAREN_"]=3
    cConf_scores["EX_glygn2_LPAREN_e_RPAREN_"]=3
    cConf_scores["EX_gsn_LPAREN_e_RPAREN_"]=3
    cConf_scores["EX_gthox_LPAREN_e_RPAREN_"]=3
    cConf_scores["EX_gthrd_LPAREN_e_RPAREN_"]=3
    cConf_scores["EX_gua_LPAREN_e_RPAREN_"]=3
    cConf_scores["EX_h_LPAREN_e_RPAREN_"]=3
    cConf_scores["EX_h2o_LPAREN_e_RPAREN_"]=3
    cConf_scores["EX_ha_LPAREN_e_RPAREN_"]=3
    cConf_scores["EX_HC00250_LPAREN_e_RPAREN_"]=3
    cConf_scores["EX_HC01609_LPAREN_e_RPAREN_"]=3
    cConf_scores["EX_HC01610_LPAREN_e_RPAREN_"]=3
    cConf_scores["EX_hdca_LPAREN_e_RPAREN_"]=3
    cConf_scores["EX_hdcea_LPAREN_e_RPAREN_"]=3
    cConf_scores["EX_his_L_LPAREN_e_RPAREN_"]=3
    cConf_scores["EX_hxan_LPAREN_e_RPAREN_"]=3
    cConf_scores["EX_ile_L_LPAREN_e_RPAREN_"]=3
    cConf_scores["EX_inost_LPAREN_e_RPAREN_"]=3
    cConf_scores["EX_k_LPAREN_e_RPAREN_"]=3
    cConf_scores["EX_lac_L_LPAREN_e_RPAREN_"]=3
    cConf_scores["EX_lac_D_LPAREN_e_RPAREN_"]=-1
    cConf_scores["EX_lcts_LPAREN_e_RPAREN_"]=3
    cConf_scores["EX_leu_L_LPAREN_e_RPAREN_"]=3
    cConf_scores["EX_leuktrA4_LPAREN_e_RPAREN_"]=3
    cConf_scores["EX_leuktrD4_LPAREN_e_RPAREN_"]=3
    cConf_scores["EX_leuktrE4_LPAREN_e_RPAREN_"]=3
    cConf_scores["EX_lnlc_LPAREN_e_RPAREN_"]=3
    cConf_scores["EX_lnlnca_LPAREN_e_RPAREN_"]=3
    cConf_scores["EX_lnlncg_LPAREN_e_RPAREN_"]=3
    cConf_scores["EX_lpchol_hs_LPAREN_e_RPAREN_"]=3
    cConf_scores["EX_lys_L_LPAREN_e_RPAREN_"]=3
    cConf_scores["EX_mag_hs_LPAREN_e_RPAREN_"]=3
    cConf_scores["EX_mal_L_LPAREN_e_RPAREN_"]=3
    cConf_scores["EX_malt_LPAREN_e_RPAREN_"]=3
    cConf_scores["EX_man_LPAREN_e_RPAREN_"]=3
    cConf_scores["EX_meoh_LPAREN_e_RPAREN_"]=3
    cConf_scores["EX_met_L_LPAREN_e_RPAREN_"]=3
    cConf_scores["EX_n2m2nmasn_LPAREN_e_RPAREN_"]=3
    cConf_scores["EX_na1_LPAREN_e_RPAREN_"]=3
    cConf_scores["EX_nac_LPAREN_e_RPAREN_"]=3
    cConf_scores["EX_ncam_LPAREN_e_RPAREN_"]=3
    cConf_scores["EX_nh4_LPAREN_e_RPAREN_"]=3
    cConf_scores["EX_no2_LPAREN_e_RPAREN_"]=3
    cConf_scores["EX_o2_LPAREN_e_RPAREN_"]=3
    cConf_scores["EX_o2s_LPAREN_e_RPAREN_"]=3
    cConf_scores["EX_ocdca_LPAREN_e_RPAREN_"]=3
    cConf_scores["EX_ocdcea_LPAREN_e_RPAREN_"]=3
    cConf_scores["EX_octa_LPAREN_e_RPAREN_"]=3
    cConf_scores["EX_orn_LPAREN_e_RPAREN_"]=3
    cConf_scores["EX_oxa_LPAREN_e_RPAREN_"]=3
    cConf_scores["EX_pe_hs_LPAREN_e_RPAREN_"]=3
    cConf_scores["EX_pglyc_hs_LPAREN_e_RPAREN_"]=3
    cConf_scores["EX_phe_L_LPAREN_e_RPAREN_"]=3
    cConf_scores["EX_pheme_LPAREN_e_RPAREN_"]=3
    cConf_scores["EX_pi_LPAREN_e_RPAREN_"]=3
    cConf_scores["EX_pnto_R_LPAREN_e_RPAREN_"]=3
    cConf_scores["EX_ppa_LPAREN_e_RPAREN_"]=3
    cConf_scores["EX_pro_L_LPAREN_e_RPAREN_"]=3
    cConf_scores["EX_prostgh2_LPAREN_e_RPAREN_"]=3
    cConf_scores["EX_ps_hs_LPAREN_e_RPAREN_"]=3
    cConf_scores["EX_ptrc_LPAREN_e_RPAREN_"]=3
    cConf_scores["EX_ptvstlac_LPAREN_e_RPAREN_"]=3
    cConf_scores["EX_pydam_LPAREN_e_RPAREN_"]=3
    cConf_scores["EX_pydx_LPAREN_e_RPAREN_"]=3
    cConf_scores["EX_pydx5p_LPAREN_e_RPAREN_"]=3
    cConf_scores["EX_pydxn_LPAREN_e_RPAREN_"]=3
    cConf_scores["EX_pyr_LPAREN_e_RPAREN_"]=3
    cConf_scores["EX_q10h2_LPAREN_e_RPAREN_"]=3
    cConf_scores["EX_retfa_LPAREN_e_RPAREN_"]=3
    cConf_scores["EX_retinol_LPAREN_e_RPAREN_"]=3
    cConf_scores["EX_retn_LPAREN_e_RPAREN_"]=3
    cConf_scores["EX_rib_D_LPAREN_e_RPAREN_"]=3
    cConf_scores["EX_ribflv_LPAREN_e_RPAREN_"]=3
    cConf_scores["EX_sbt_DASH_d_LPAREN_e_RPAREN_"]=3
    cConf_scores["EX_sel_LPAREN_e_RPAREN_"]=3
    cConf_scores["EX_ser_L_LPAREN_e_RPAREN_"]=3
    cConf_scores["EX_so4_LPAREN_e_RPAREN_"]=3
    cConf_scores["EX_sph1p_LPAREN_e_RPAREN_"]=3
    cConf_scores["EX_spmd_LPAREN_e_RPAREN_"]=3
    cConf_scores["EX_strch1_LPAREN_e_RPAREN_"]=3
    cConf_scores["EX_strch2_LPAREN_e_RPAREN_"]=3
    cConf_scores["EX_sucr_LPAREN_e_RPAREN_"]=3
    cConf_scores["EX_tag_hs_LPAREN_e_RPAREN_"]=3
    cConf_scores["EX_thm_LPAREN_e_RPAREN_"]=3
    cConf_scores["EX_thr_L_LPAREN_e_RPAREN_"]=3
    cConf_scores["EX_thymd_LPAREN_e_RPAREN_"]=3
    cConf_scores["EX_tre_LPAREN_e_RPAREN_"]=3
    cConf_scores["EX_trp_L_LPAREN_e_RPAREN_"]=3
    cConf_scores["EX_ttdca_LPAREN_e_RPAREN_"]=3
    cConf_scores["EX_tyr_L_LPAREN_e_RPAREN_"]=3
    cConf_scores["EX_ura_LPAREN_e_RPAREN_"]=3
    cConf_scores["EX_urea_LPAREN_e_RPAREN_"]=3
    cConf_scores["EX_uri_LPAREN_e_RPAREN_"]=3
    cConf_scores["EX_utp_LPAREN_e_RPAREN_"]=3
    cConf_scores["EX_val_L_LPAREN_e_RPAREN_"]=3
    cConf_scores["EX_xmp_LPAREN_e_RPAREN_"]=3
    cConf_scores["EX_xyl_D_LPAREN_e_RPAREN_"]=3
    cConf_scores["EX_xylt_LPAREN_e_RPAREN_"]=3
    cConf_scores["biomass_reaction"]=3
    cConf_scores["biomass_DNA"]=3
    cConf_scores["biomass_RNA"]=3
    cConf_scores["biomass_carbohydrate"]=3
    cConf_scores["biomass_lipid"]=3
    cConf_scores["biomass_other"]=3
    cConf_scores["biomass_protein"]=3
    cConf_scores["DM_atp_c_"]=3

    nConf_scores["DM_4hrpo"]=3
    nConf_scores["DM_datp_n_"]=3
    nConf_scores["DM_dctp_n_"]=3
    nConf_scores["DM_dgtp_n_"]=3
    nConf_scores["DM_dttp_n_"]=3
    nConf_scores["DM_Lcystin"]=3
    nConf_scores["DM_pe_hs_LPAREN_r_RPAREN_"]=3
    nConf_scores["EX_2hb_LPAREN_e_RPAREN_"]=3
    nConf_scores["EX_34hpp_"]=3
    nConf_scores["EX_3hpvs_LPAREN_e_RPAREN_"]=3
    nConf_scores["EX_3mob_LPAREN_e_RPAREN_"]=3
    nConf_scores["EX_4mop_LPAREN_e_RPAREN_"]=3
    nConf_scores["EX_ac_LPAREN_e_RPAREN_"]=3
    nConf_scores["EX_acac_LPAREN_e_RPAREN_"]=3
    nConf_scores["EX_acald_LPAREN_e_RPAREN_"]=3
    nConf_scores["EX_acetone_LPAREN_e_RPAREN_"]=3
    nConf_scores["EX_acgam_LPAREN_e_RPAREN_"]=3
    nConf_scores["EX_acmana_LPAREN_e_RPAREN_"]=3
    nConf_scores["EX_ade_LPAREN_e_RPAREN_"]=3
    nConf_scores["EX_adn_LPAREN_e_RPAREN_"]=3
    nConf_scores["EX_adpcbl_LPAREN_e_RPAREN_"]=3
    nConf_scores["EX_akg_LPAREN_e_RPAREN_"]=3
    nConf_scores["EX_ala_B_LPAREN_e_RPAREN_"]=3
    nConf_scores["EX_ala_D_LPAREN_e_RPAREN_"]=3
    nConf_scores["EX_ala_L_LPAREN_e_RPAREN_"]=3
    nConf_scores["EX_am9csa_LPAREN_e_RPAREN_"]=3
    nConf_scores["EX_amp_LPAREN_e_RPAREN_"]=3
    nConf_scores["EX_arab_L_LPAREN_e_RPAREN_"]=3
    nConf_scores["EX_arachd_LPAREN_e_RPAREN_"]=3
    nConf_scores["EX_arg_L_LPAREN_e_RPAREN_"]=3
    nConf_scores["EX_asn_L_LPAREN_e_RPAREN_"]=3
    nConf_scores["EX_asp_L_LPAREN_e_RPAREN_"]=3
    nConf_scores["EX_atp_LPAREN_e_RPAREN_"]=3
    nConf_scores["EX_bhb_LPAREN_e_RPAREN_"]=3
    nConf_scores["EX_btn_LPAREN_e_RPAREN_"]=3
    nConf_scores["EX_ca2_LPAREN_e_RPAREN_"]=3
    nConf_scores["EX_carn_LPAREN_e_RPAREN_"]=3
    nConf_scores["EX_caro_LPAREN_e_RPAREN_"]=3
    nConf_scores["EX_cgly_LPAREN_e_RPAREN_"]=3
    nConf_scores["EX_chol_LPAREN_e_RPAREN_"]=3
    nConf_scores["EX_chsterol_LPAREN_e_RPAREN_"]=3
    nConf_scores["EX_cit_LPAREN_e_RPAREN_"]=3
    nConf_scores["EX_cl_LPAREN_e_RPAREN_"]=3
    nConf_scores["EX_CLPND_LPAREN_e_RPAREN_"]=3
    nConf_scores["EX_co2_LPAREN_e_RPAREN_"]=3
    nConf_scores["EX_creat_LPAREN_e_RPAREN_"]=3
    nConf_scores["EX_crn_LPAREN_e_RPAREN_"]=3
    nConf_scores["EX_crvnc_LPAREN_e_RPAREN_"]=3
    nConf_scores["EX_csa_LPAREN_e_RPAREN_"]=3
    nConf_scores["EX_csn_LPAREN_e_RPAREN_"]=3
    nConf_scores["EX_cys_L_LPAREN_e_RPAREN_"]=3
    nConf_scores["EX_cytd_LPAREN_e_RPAREN_"]=3
    nConf_scores["EX_dad_2_LPAREN_e_RPAREN_"]=3
    nConf_scores["EX_dag_hs_LPAREN_e_RPAREN_"]=3
    nConf_scores["EX_dcmp_LPAREN_e_RPAREN_"]=3
    nConf_scores["EX_dcyt_LPAREN_e_RPAREN_"]=3
    nConf_scores["EX_ddca_LPAREN_e_RPAREN_"]=3
    nConf_scores["EX_dgsn_LPAREN_e_RPAREN_"]=3
    nConf_scores["EX_dhdascb_LPAREN_e_RPAREN_"]=3
    nConf_scores["EX_din_LPAREN_e_RPAREN_"]=3
    nConf_scores["EX_dopa_LPAREN_e_RPAREN_"]=3
    nConf_scores["EX_drib_LPAREN_e_RPAREN_"]=3
    nConf_scores["EX_etoh_LPAREN_e_RPAREN_"]=3
    nConf_scores["EX_fald_LPAREN_e_RPAREN_"]=3
    nConf_scores["EX_fe2_LPAREN_e_RPAREN_"]=3
    nConf_scores["EX_fe2_LPAREN_e_RPAREN_"]=3
    nConf_scores["EX_fe3_LPAREN_e_RPAREN_"]=3
    nConf_scores["EX_fmn_LPAREN_e_RPAREN_"]=3
    nConf_scores["EX_fol_LPAREN_e_RPAREN_"]=3
    nConf_scores["EX_for_LPAREN_e_RPAREN_"]=3
    nConf_scores["EX_fru_LPAREN_e_RPAREN_"]=3
    nConf_scores["EX_fuc_L_LPAREN_e_RPAREN_"]=3
    nConf_scores["EX_fum_LPAREN_e_RPAREN_"]=3
    nConf_scores["EX_gal_LPAREN_e_RPAREN_"]=3
    nConf_scores["EX_gam_LPAREN_e_RPAREN_"]=3
    nConf_scores["EX_gchola_LPAREN_e_RPAREN_"]=3
    nConf_scores["EX_glc_LPAREN_e_RPAREN_"]=3
    nConf_scores["EX_glcur_LPAREN_e_RPAREN_"]=3
    nConf_scores["EX_gln_L_LPAREN_e_RPAREN_"]=3
    nConf_scores["EX_glu_L_LPAREN_e_RPAREN_"]=3
    nConf_scores["EX_gluala_LPAREN_e_RPAREN_"]=3
    nConf_scores["EX_gly_LPAREN_e_RPAREN_"]=3
    nConf_scores["EX_glyb_LPAREN_e_RPAREN_"]=3
    nConf_scores["EX_glyc_LPAREN_e_RPAREN_"]=3
    nConf_scores["EX_glyc3p_LPAREN_e_RPAREN_"]=3
    nConf_scores["EX_glygn2_LPAREN_e_RPAREN_"]=3
    nConf_scores["EX_gsn_LPAREN_e_RPAREN_"]=3
    nConf_scores["EX_gthox_LPAREN_e_RPAREN_"]=3
    nConf_scores["EX_gthrd_LPAREN_e_RPAREN_"]=3
    nConf_scores["EX_gua_LPAREN_e_RPAREN_"]=3
    nConf_scores["EX_h_LPAREN_e_RPAREN_"]=3
    nConf_scores["EX_h2o_LPAREN_e_RPAREN_"]=3
    nConf_scores["EX_ha_LPAREN_e_RPAREN_"]=3
    nConf_scores["EX_HC00250_LPAREN_e_RPAREN_"]=3
    nConf_scores["EX_HC01609_LPAREN_e_RPAREN_"]=3
    nConf_scores["EX_HC01610_LPAREN_e_RPAREN_"]=3
    nConf_scores["EX_hdca_LPAREN_e_RPAREN_"]=3
    nConf_scores["EX_hdcea_LPAREN_e_RPAREN_"]=3
    nConf_scores["EX_his_L_LPAREN_e_RPAREN_"]=3
    nConf_scores["EX_hxan_LPAREN_e_RPAREN_"]=3
    nConf_scores["EX_ile_L_LPAREN_e_RPAREN_"]=3
    nConf_scores["EX_inost_LPAREN_e_RPAREN_"]=3
    nConf_scores["EX_k_LPAREN_e_RPAREN_"]=3
    nConf_scores["EX_lac_L_LPAREN_e_RPAREN_"]=3
    nConf_scores["EX_lac_D_LPAREN_e_RPAREN_"]=-1
    nConf_scores["EX_lcts_LPAREN_e_RPAREN_"]=3
    nConf_scores["EX_leu_L_LPAREN_e_RPAREN_"]=3
    nConf_scores["EX_leuktrA4_LPAREN_e_RPAREN_"]=3
    nConf_scores["EX_leuktrD4_LPAREN_e_RPAREN_"]=3
    nConf_scores["EX_leuktrE4_LPAREN_e_RPAREN_"]=3
    nConf_scores["EX_lnlc_LPAREN_e_RPAREN_"]=3
    nConf_scores["EX_lnlnca_LPAREN_e_RPAREN_"]=3
    nConf_scores["EX_lnlncg_LPAREN_e_RPAREN_"]=3
    nConf_scores["EX_lpchol_hs_LPAREN_e_RPAREN_"]=3
    nConf_scores["EX_lys_L_LPAREN_e_RPAREN_"]=3
    nConf_scores["EX_mag_hs_LPAREN_e_RPAREN_"]=3
    nConf_scores["EX_mal_L_LPAREN_e_RPAREN_"]=3
    nConf_scores["EX_malt_LPAREN_e_RPAREN_"]=3
    nConf_scores["EX_man_LPAREN_e_RPAREN_"]=3
    nConf_scores["EX_meoh_LPAREN_e_RPAREN_"]=3
    nConf_scores["EX_met_L_LPAREN_e_RPAREN_"]=3
    nConf_scores["EX_n2m2nmasn_LPAREN_e_RPAREN_"]=3
    nConf_scores["EX_na1_LPAREN_e_RPAREN_"]=3
    nConf_scores["EX_nac_LPAREN_e_RPAREN_"]=3
    nConf_scores["EX_ncam_LPAREN_e_RPAREN_"]=3
    nConf_scores["EX_nh4_LPAREN_e_RPAREN_"]=3
    nConf_scores["EX_no2_LPAREN_e_RPAREN_"]=3
    nConf_scores["EX_o2_LPAREN_e_RPAREN_"]=3
    nConf_scores["EX_o2s_LPAREN_e_RPAREN_"]=3
    nConf_scores["EX_ocdca_LPAREN_e_RPAREN_"]=3
    nConf_scores["EX_ocdcea_LPAREN_e_RPAREN_"]=3
    nConf_scores["EX_octa_LPAREN_e_RPAREN_"]=3
    nConf_scores["EX_orn_LPAREN_e_RPAREN_"]=3
    nConf_scores["EX_oxa_LPAREN_e_RPAREN_"]=3
    nConf_scores["EX_pe_hs_LPAREN_e_RPAREN_"]=3
    nConf_scores["EX_pglyc_hs_LPAREN_e_RPAREN_"]=3
    nConf_scores["EX_phe_L_LPAREN_e_RPAREN_"]=3
    nConf_scores["EX_pheme_LPAREN_e_RPAREN_"]=3
    nConf_scores["EX_pi_LPAREN_e_RPAREN_"]=3
    nConf_scores["EX_pnto_R_LPAREN_e_RPAREN_"]=3
    nConf_scores["EX_ppa_LPAREN_e_RPAREN_"]=3
    nConf_scores["EX_pro_L_LPAREN_e_RPAREN_"]=3
    nConf_scores["EX_prostgh2_LPAREN_e_RPAREN_"]=3
    nConf_scores["EX_ps_hs_LPAREN_e_RPAREN_"]=3
    nConf_scores["EX_ptrc_LPAREN_e_RPAREN_"]=3
    nConf_scores["EX_ptvstlac_LPAREN_e_RPAREN_"]=3
    nConf_scores["EX_pydam_LPAREN_e_RPAREN_"]=3
    nConf_scores["EX_pydx_LPAREN_e_RPAREN_"]=3
    nConf_scores["EX_pydx5p_LPAREN_e_RPAREN_"]=3
    nConf_scores["EX_pydxn_LPAREN_e_RPAREN_"]=3
    nConf_scores["EX_pyr_LPAREN_e_RPAREN_"]=3
    nConf_scores["EX_q10h2_LPAREN_e_RPAREN_"]=3
    nConf_scores["EX_retfa_LPAREN_e_RPAREN_"]=3
    nConf_scores["EX_retinol_LPAREN_e_RPAREN_"]=3
    nConf_scores["EX_retn_LPAREN_e_RPAREN_"]=3
    nConf_scores["EX_rib_D_LPAREN_e_RPAREN_"]=3
    nConf_scores["EX_ribflv_LPAREN_e_RPAREN_"]=3
    nConf_scores["EX_sbt_DASH_d_LPAREN_e_RPAREN_"]=3
    nConf_scores["EX_sel_LPAREN_e_RPAREN_"]=3
    nConf_scores["EX_ser_L_LPAREN_e_RPAREN_"]=3
    nConf_scores["EX_so4_LPAREN_e_RPAREN_"]=3
    nConf_scores["EX_sph1p_LPAREN_e_RPAREN_"]=3
    nConf_scores["EX_spmd_LPAREN_e_RPAREN_"]=3
    nConf_scores["EX_strch1_LPAREN_e_RPAREN_"]=3
    nConf_scores["EX_strch2_LPAREN_e_RPAREN_"]=3
    nConf_scores["EX_sucr_LPAREN_e_RPAREN_"]=3
    nConf_scores["EX_tag_hs_LPAREN_e_RPAREN_"]=3
    nConf_scores["EX_thm_LPAREN_e_RPAREN_"]=3
    nConf_scores["EX_thr_L_LPAREN_e_RPAREN_"]=3
    nConf_scores["EX_thymd_LPAREN_e_RPAREN_"]=3
    nConf_scores["EX_tre_LPAREN_e_RPAREN_"]=3
    nConf_scores["EX_trp_L_LPAREN_e_RPAREN_"]=3
    nConf_scores["EX_ttdca_LPAREN_e_RPAREN_"]=3
    nConf_scores["EX_tyr_L_LPAREN_e_RPAREN_"]=3
    nConf_scores["EX_ura_LPAREN_e_RPAREN_"]=3
    nConf_scores["EX_urea_LPAREN_e_RPAREN_"]=3
    nConf_scores["EX_uri_LPAREN_e_RPAREN_"]=3
    nConf_scores["EX_utp_LPAREN_e_RPAREN_"]=3
    nConf_scores["EX_val_L_LPAREN_e_RPAREN_"]=3
    nConf_scores["EX_xmp_LPAREN_e_RPAREN_"]=3
    nConf_scores["EX_xyl_D_LPAREN_e_RPAREN_"]=3
    nConf_scores["EX_xylt_LPAREN_e_RPAREN_"]=3
    nConf_scores["biomass_reaction"]=3
    nConf_scores["biomass_DNA"]=3
    nConf_scores["biomass_RNA"]=3
    nConf_scores["biomass_carbohydrate"]=3
    nConf_scores["biomass_lipid"]=3
    nConf_scores["biomass_other"]=3
    nConf_scores["biomass_protein"]=3
    nConf_scores["DM_atp_c_"]=3
   
    opt_c= CORDA(model=Recon, confidence=cConf_scores, n=5, met_prod=metabolites,  penalty_factor=1000 ) 
    opt_n= CORDA(model=Recon, confidence=nConf_scores, n=5, met_prod=metabolites,  penalty_factor=1000 )
    opt_c.build()
    opt_n.build()
    model_c=opt_c.cobra_model(name="C")
    model_n=opt_n.cobra_model(name="N")
    cobra.io.write_sbml_model(model_c, "c_model_randomized_"+str(i)+".sbml")
    cobra.io.write_sbml_model(model_n, "n_model_randomized_"+str(i)+".sbml")

    
    targets=get_drugable_targets(model_n, model_c, "biopsys" )
    targets.to_csv(f)
    print("finished"+str(i))
    
```


```python
output = mp.Queue(maxsize=12)
pool = mp.Pool(processes=12)
    
for i in range(656,1001):
    pool.apply_async(multiple_randomized_analysis, args=(conf_CancerBiopsy,conf_NormalBiopsy,Recon2,i, metas))
pool.close()
pool.join()
```


```python

```
