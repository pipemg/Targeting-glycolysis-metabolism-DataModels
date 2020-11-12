# Generating the biopsis models 


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
#Recon204 = cobra.io.load_matlab_model("Models/Recon2.v04.mat")
Recon2.reactions.OIVD3m.gene_reaction_rule='HGNC:2698 and HGNC:987 and HGNC:986 and HGNC:2898 or HGNC:2698 and HGNC:2898 and HGNC:987 and HGNC:986'
Recon2.reactions.OIVD1m.gene_reaction_rule='HGNC:2698 and HGNC:987 and HGNC:986 and HGNC:2898 or HGNC:2698 and HGNC:2898 and HGNC:987 and HGNC:986'
Recon2.reactions.OIVD2m.gene_reaction_rule='HGNC:2698 and HGNC:987 and HGNC:986 and HGNC:2898 or HGNC:2698 and HGNC:2898 and HGNC:987 and HGNC:986'
```


Structure of metabolic reactions

```python
print("Reactions:", len(Recon2.reactions))
print("Metabolites:", len(Recon2.metabolites))
```

    Reactions: 7785
    Metabolites: 5324



```python
bm = Recon2.reactions.get_by_id("biomass_reaction")
print(bm.build_reaction_string())
```

    0.014 biomass_DNA_c + 0.058 biomass_RNA_c + 0.071 biomass_carbohydrate_c + 0.097 biomass_lipid_c + 0.054 biomass_other_c + 0.706 biomass_protein_c --> 



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
cp=Recon2.copy()
cp.optimize()
cp.summary()
```

    IN FLUXES                  OUT FLUXES           OBJECTIVES
    -------------------------  -------------------  -----------------------
    o2_e             0.41      pi_e       0.515     biomass_reac...  0.0384
    fe2_e            0.406     fe3_e      0.406
    atp_e            0.167     pyr_e      0.278
    so4_e            0.05      co2_e      0.216
    gly_e            0.03      ade_e      0.0822
    Lcystin_c        0.024     adn_e      0.0627
    val_L_e          0.0237    slfcys_e   0.05
    lys_L_e          0.0227    ura_e      0.0405
    utp_e            0.02      nh4_e      0.0337
    gln_L_e          0.0172    elaid_e    0.0188
    arg_L_e          0.012     glyc_e     0.0185
    asn_L_e          0.012     aicar_e    0.0154
    asp_L_e          0.012     gua_e      0.0153
    cytd_e           0.012     3mob_e     0.0101
    ddca_e           0.012     strdnc_e   0.0101
    for_e            0.012     dcsptn1_e  0.00775
    fru_e            0.012     tmndnc_e   0.00775
    glu_L_e          0.012     tchola_e   0.00428
    gsn_e            0.012     3deccrn_e  0.00187
    ile_L_e          0.012     acac_e     0.00187
    inost_e          0.012     lac_L_e    0.00187
    leu_L_e          0.012     3mop_e     0.00102
    lnlnca_e         0.012     dag_hs_e   0.000587
    malt_e           0.012     ac_e       0.000255
    man_e            0.012     tag_hs_e   8.4e-05
    o2s_e            0.012
    ocdcea_e         0.012
    phe_L_e          0.012
    pro_L_e          0.012
    rib_D_e          0.012
    ser_L_e          0.012
    thr_L_e          0.012
    uri_e            0.012
    4mop_e           0.0108
    cys_L_e          0.00806
    cit_e            0.008
    crvnc_e          0.00775
    lnlc_e           0.00775
    ps_hs_e          0.00731
    met_L_e          0.00587
    chol_e           0.00574
    his_L_e          0.00485
    gchola_e         0.00428
    tyr_L_e          0.00409
    orn_e            0.00272
    pe_hs_e          0.00213
    biomass_other_c  0.00207
    crn_e            0.00187
    octa_e           0.00187
    creat_e          0.00179
    sbt_D_e          0.00107
    pe_hs_r          0.000857
    sph1p_e          0.000671
    pglyc_hs_e       0.000559
    trp_L_e          0.000511
    datp_n           0.000506
    dttp_n           0.000502
    dgtp_n           0.00038
    dctp_n           0.000362



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
%matplotlib inline
confidence_scores_matrix["MeanCancerBiopsy"].hist(bins=20)
```




    <matplotlib.axes._subplots.AxesSubplot at 0x7f93ab250ba8>




![png](output_7_1.png)



```python
%matplotlib inline
confidence_scores_matrix["MeanNormalBiopsy"].hist(bins=10)
```




    <matplotlib.axes._subplots.AxesSubplot at 0x7f93ab139240>




![png](output_8_1.png)



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
%matplotlib inline
confidence_scores_matrix["MeanNormalBiopsy"].hist(bins=10)
```




    <matplotlib.axes._subplots.AxesSubplot at 0x7f93ab029780>




![png](output_10_1.png)



```python
%matplotlib inline
confidence_scores_matrix["MeanCancerBiopsy"].hist(bins=10)
```

![png](output_11_1.png)



```python
from corda import reaction_confidence


conf_MeanCancerBiopsy = {}
conf_MeanNormalBiopsy = {}

for r in Recon2.reactions:
    if(r.gene_reaction_rule!=''):
        conf_MeanCancerBiopsy[r.id]=reaction_confidence(r.gene_reaction_rule,confidence_scores_matrix["MaxCancerBiopsy"])
        conf_MeanNormalBiopsy[r.id]=reaction_confidence(r.gene_reaction_rule,confidence_scores_matrix["MaxNormalBiopsy"])
    else:
        conf_MeanCancerBiopsy[r.id]=1 
        conf_MeanNormalBiopsy[r.id]=1
```


```python
%matplotlib inline
import pandas as pd

df=pd.DataFrame({'Cancer': conf_MeanCancerBiopsy})
df.Cancer.value_counts()

```




     1.0    3450
     3.0    1441
     0.0    1381
    -1.0    1275
     2.0     238
    Name: Cancer, dtype: int64




```python
conf_MeanCancerBiopsy["DM_4hrpo"]=3
conf_MeanCancerBiopsy["DM_datp_n_"]=3
conf_MeanCancerBiopsy["DM_dctp_n_"]=3
conf_MeanCancerBiopsy["DM_dgtp_n_"]=3
conf_MeanCancerBiopsy["DM_dttp_n_"]=3
conf_MeanCancerBiopsy["DM_Lcystin"]=3
conf_MeanCancerBiopsy["DM_pe_hs_LPAREN_r_RPAREN_"]=3
conf_MeanCancerBiopsy["EX_2hb_LPAREN_e_RPAREN_"]=3
conf_MeanCancerBiopsy["EX_34hpp_"]=3
conf_MeanCancerBiopsy["EX_3hpvs_LPAREN_e_RPAREN_"]=3
conf_MeanCancerBiopsy["EX_3mob_LPAREN_e_RPAREN_"]=3
conf_MeanCancerBiopsy["EX_4mop_LPAREN_e_RPAREN_"]=3
conf_MeanCancerBiopsy["EX_ac_LPAREN_e_RPAREN_"]=3
conf_MeanCancerBiopsy["EX_acac_LPAREN_e_RPAREN_"]=3
conf_MeanCancerBiopsy["EX_acald_LPAREN_e_RPAREN_"]=3
conf_MeanCancerBiopsy["EX_acetone_LPAREN_e_RPAREN_"]=3
conf_MeanCancerBiopsy["EX_acgam_LPAREN_e_RPAREN_"]=3
conf_MeanCancerBiopsy["EX_acmana_LPAREN_e_RPAREN_"]=3
conf_MeanCancerBiopsy["EX_ade_LPAREN_e_RPAREN_"]=3
conf_MeanCancerBiopsy["EX_adn_LPAREN_e_RPAREN_"]=3
conf_MeanCancerBiopsy["EX_adpcbl_LPAREN_e_RPAREN_"]=3
conf_MeanCancerBiopsy["EX_akg_LPAREN_e_RPAREN_"]=3
conf_MeanCancerBiopsy["EX_ala_B_LPAREN_e_RPAREN_"]=3
conf_MeanCancerBiopsy["EX_ala_D_LPAREN_e_RPAREN_"]=3
conf_MeanCancerBiopsy["EX_ala_L_LPAREN_e_RPAREN_"]=3
conf_MeanCancerBiopsy["EX_am9csa_LPAREN_e_RPAREN_"]=3
conf_MeanCancerBiopsy["EX_amp_LPAREN_e_RPAREN_"]=3
conf_MeanCancerBiopsy["EX_arab_L_LPAREN_e_RPAREN_"]=3
conf_MeanCancerBiopsy["EX_arachd_LPAREN_e_RPAREN_"]=3
conf_MeanCancerBiopsy["EX_arg_L_LPAREN_e_RPAREN_"]=3
conf_MeanCancerBiopsy["EX_asn_L_LPAREN_e_RPAREN_"]=3
conf_MeanCancerBiopsy["EX_asp_L_LPAREN_e_RPAREN_"]=3
conf_MeanCancerBiopsy["EX_atp_LPAREN_e_RPAREN_"]=3
conf_MeanCancerBiopsy["EX_bhb_LPAREN_e_RPAREN_"]=3
conf_MeanCancerBiopsy["EX_btn_LPAREN_e_RPAREN_"]=3
conf_MeanCancerBiopsy["EX_ca2_LPAREN_e_RPAREN_"]=3
conf_MeanCancerBiopsy["EX_carn_LPAREN_e_RPAREN_"]=3
conf_MeanCancerBiopsy["EX_caro_LPAREN_e_RPAREN_"]=3
conf_MeanCancerBiopsy["EX_cgly_LPAREN_e_RPAREN_"]=3
conf_MeanCancerBiopsy["EX_chol_LPAREN_e_RPAREN_"]=3
conf_MeanCancerBiopsy["EX_chsterol_LPAREN_e_RPAREN_"]=3
conf_MeanCancerBiopsy["EX_cit_LPAREN_e_RPAREN_"]=3
conf_MeanCancerBiopsy["EX_cl_LPAREN_e_RPAREN_"]=3
conf_MeanCancerBiopsy["EX_CLPND_LPAREN_e_RPAREN_"]=3
conf_MeanCancerBiopsy["EX_co2_LPAREN_e_RPAREN_"]=3
conf_MeanCancerBiopsy["EX_creat_LPAREN_e_RPAREN_"]=3
conf_MeanCancerBiopsy["EX_crn_LPAREN_e_RPAREN_"]=3
conf_MeanCancerBiopsy["EX_crvnc_LPAREN_e_RPAREN_"]=3
conf_MeanCancerBiopsy["EX_csa_LPAREN_e_RPAREN_"]=3
conf_MeanCancerBiopsy["EX_csn_LPAREN_e_RPAREN_"]=3
conf_MeanCancerBiopsy["EX_cys_L_LPAREN_e_RPAREN_"]=3
conf_MeanCancerBiopsy["EX_cytd_LPAREN_e_RPAREN_"]=3
conf_MeanCancerBiopsy["EX_dad_2_LPAREN_e_RPAREN_"]=3
conf_MeanCancerBiopsy["EX_dag_hs_LPAREN_e_RPAREN_"]=3
conf_MeanCancerBiopsy["EX_dcmp_LPAREN_e_RPAREN_"]=3
conf_MeanCancerBiopsy["EX_dcyt_LPAREN_e_RPAREN_"]=3
conf_MeanCancerBiopsy["EX_ddca_LPAREN_e_RPAREN_"]=3
conf_MeanCancerBiopsy["EX_dgsn_LPAREN_e_RPAREN_"]=3
conf_MeanCancerBiopsy["EX_dhdascb_LPAREN_e_RPAREN_"]=3
conf_MeanCancerBiopsy["EX_din_LPAREN_e_RPAREN_"]=3
conf_MeanCancerBiopsy["EX_dopa_LPAREN_e_RPAREN_"]=3
conf_MeanCancerBiopsy["EX_drib_LPAREN_e_RPAREN_"]=3
conf_MeanCancerBiopsy["EX_etoh_LPAREN_e_RPAREN_"]=3
conf_MeanCancerBiopsy["EX_fald_LPAREN_e_RPAREN_"]=3
conf_MeanCancerBiopsy["EX_fe2_LPAREN_e_RPAREN_"]=3
conf_MeanCancerBiopsy["EX_fe2_LPAREN_e_RPAREN_"]=3
conf_MeanCancerBiopsy["EX_fe3_LPAREN_e_RPAREN_"]=3
conf_MeanCancerBiopsy["EX_fmn_LPAREN_e_RPAREN_"]=3
conf_MeanCancerBiopsy["EX_fol_LPAREN_e_RPAREN_"]=3
conf_MeanCancerBiopsy["EX_for_LPAREN_e_RPAREN_"]=3
conf_MeanCancerBiopsy["EX_fru_LPAREN_e_RPAREN_"]=3
conf_MeanCancerBiopsy["EX_fuc_L_LPAREN_e_RPAREN_"]=3
conf_MeanCancerBiopsy["EX_fum_LPAREN_e_RPAREN_"]=3
conf_MeanCancerBiopsy["EX_gal_LPAREN_e_RPAREN_"]=3
conf_MeanCancerBiopsy["EX_gam_LPAREN_e_RPAREN_"]=3
conf_MeanCancerBiopsy["EX_gchola_LPAREN_e_RPAREN_"]=3
conf_MeanCancerBiopsy["EX_glc_LPAREN_e_RPAREN_"]=3
conf_MeanCancerBiopsy["EX_glcur_LPAREN_e_RPAREN_"]=3
conf_MeanCancerBiopsy["EX_gln_L_LPAREN_e_RPAREN_"]=3
conf_MeanCancerBiopsy["EX_glu_L_LPAREN_e_RPAREN_"]=3
conf_MeanCancerBiopsy["EX_gluala_LPAREN_e_RPAREN_"]=3
conf_MeanCancerBiopsy["EX_gly_LPAREN_e_RPAREN_"]=3
conf_MeanCancerBiopsy["EX_glyb_LPAREN_e_RPAREN_"]=3
conf_MeanCancerBiopsy["EX_glyc_LPAREN_e_RPAREN_"]=3
conf_MeanCancerBiopsy["EX_glyc3p_LPAREN_e_RPAREN_"]=3
conf_MeanCancerBiopsy["EX_glygn2_LPAREN_e_RPAREN_"]=3
conf_MeanCancerBiopsy["EX_gsn_LPAREN_e_RPAREN_"]=3
conf_MeanCancerBiopsy["EX_gthox_LPAREN_e_RPAREN_"]=3
conf_MeanCancerBiopsy["EX_gthrd_LPAREN_e_RPAREN_"]=3
conf_MeanCancerBiopsy["EX_gua_LPAREN_e_RPAREN_"]=3
conf_MeanCancerBiopsy["EX_h_LPAREN_e_RPAREN_"]=3
conf_MeanCancerBiopsy["EX_h2o_LPAREN_e_RPAREN_"]=3
conf_MeanCancerBiopsy["EX_ha_LPAREN_e_RPAREN_"]=3
conf_MeanCancerBiopsy["EX_HC00250_LPAREN_e_RPAREN_"]=3
conf_MeanCancerBiopsy["EX_HC01609_LPAREN_e_RPAREN_"]=3
conf_MeanCancerBiopsy["EX_HC01610_LPAREN_e_RPAREN_"]=3
conf_MeanCancerBiopsy["EX_hdca_LPAREN_e_RPAREN_"]=3
conf_MeanCancerBiopsy["EX_hdcea_LPAREN_e_RPAREN_"]=3
conf_MeanCancerBiopsy["EX_his_L_LPAREN_e_RPAREN_"]=3
conf_MeanCancerBiopsy["EX_hxan_LPAREN_e_RPAREN_"]=3
conf_MeanCancerBiopsy["EX_ile_L_LPAREN_e_RPAREN_"]=3
conf_MeanCancerBiopsy["EX_inost_LPAREN_e_RPAREN_"]=3
conf_MeanCancerBiopsy["EX_k_LPAREN_e_RPAREN_"]=3
conf_MeanCancerBiopsy["EX_lac_L_LPAREN_e_RPAREN_"]=3
conf_MeanCancerBiopsy["EX_lcts_LPAREN_e_RPAREN_"]=3
conf_MeanCancerBiopsy["EX_leu_L_LPAREN_e_RPAREN_"]=3
conf_MeanCancerBiopsy["EX_leuktrA4_LPAREN_e_RPAREN_"]=3
conf_MeanCancerBiopsy["EX_leuktrD4_LPAREN_e_RPAREN_"]=3
conf_MeanCancerBiopsy["EX_leuktrE4_LPAREN_e_RPAREN_"]=3
conf_MeanCancerBiopsy["EX_lnlc_LPAREN_e_RPAREN_"]=3
conf_MeanCancerBiopsy["EX_lnlnca_LPAREN_e_RPAREN_"]=3
conf_MeanCancerBiopsy["EX_lnlncg_LPAREN_e_RPAREN_"]=3
conf_MeanCancerBiopsy["EX_lpchol_hs_LPAREN_e_RPAREN_"]=3
conf_MeanCancerBiopsy["EX_lys_L_LPAREN_e_RPAREN_"]=3
conf_MeanCancerBiopsy["EX_mag_hs_LPAREN_e_RPAREN_"]=3
conf_MeanCancerBiopsy["EX_mal_L_LPAREN_e_RPAREN_"]=3
conf_MeanCancerBiopsy["EX_malt_LPAREN_e_RPAREN_"]=3
conf_MeanCancerBiopsy["EX_man_LPAREN_e_RPAREN_"]=3
conf_MeanCancerBiopsy["EX_meoh_LPAREN_e_RPAREN_"]=3
conf_MeanCancerBiopsy["EX_met_L_LPAREN_e_RPAREN_"]=3
conf_MeanCancerBiopsy["EX_n2m2nmasn_LPAREN_e_RPAREN_"]=3
conf_MeanCancerBiopsy["EX_na1_LPAREN_e_RPAREN_"]=3
conf_MeanCancerBiopsy["EX_nac_LPAREN_e_RPAREN_"]=3
conf_MeanCancerBiopsy["EX_ncam_LPAREN_e_RPAREN_"]=3
conf_MeanCancerBiopsy["EX_nh4_LPAREN_e_RPAREN_"]=3
conf_MeanCancerBiopsy["EX_no2_LPAREN_e_RPAREN_"]=3
conf_MeanCancerBiopsy["EX_o2_LPAREN_e_RPAREN_"]=3
conf_MeanCancerBiopsy["EX_o2s_LPAREN_e_RPAREN_"]=3
conf_MeanCancerBiopsy["EX_ocdca_LPAREN_e_RPAREN_"]=3
conf_MeanCancerBiopsy["EX_ocdcea_LPAREN_e_RPAREN_"]=3
conf_MeanCancerBiopsy["EX_octa_LPAREN_e_RPAREN_"]=3
conf_MeanCancerBiopsy["EX_orn_LPAREN_e_RPAREN_"]=3
conf_MeanCancerBiopsy["EX_oxa_LPAREN_e_RPAREN_"]=3
conf_MeanCancerBiopsy["EX_pe_hs_LPAREN_e_RPAREN_"]=3
conf_MeanCancerBiopsy["EX_pglyc_hs_LPAREN_e_RPAREN_"]=3
conf_MeanCancerBiopsy["EX_phe_L_LPAREN_e_RPAREN_"]=3
conf_MeanCancerBiopsy["EX_pheme_LPAREN_e_RPAREN_"]=3
conf_MeanCancerBiopsy["EX_pi_LPAREN_e_RPAREN_"]=3
conf_MeanCancerBiopsy["EX_pnto_R_LPAREN_e_RPAREN_"]=3
conf_MeanCancerBiopsy["EX_ppa_LPAREN_e_RPAREN_"]=3
conf_MeanCancerBiopsy["EX_pro_L_LPAREN_e_RPAREN_"]=3
conf_MeanCancerBiopsy["EX_prostgh2_LPAREN_e_RPAREN_"]=3
conf_MeanCancerBiopsy["EX_ps_hs_LPAREN_e_RPAREN_"]=3
conf_MeanCancerBiopsy["EX_ptrc_LPAREN_e_RPAREN_"]=3
conf_MeanCancerBiopsy["EX_ptvstlac_LPAREN_e_RPAREN_"]=3
conf_MeanCancerBiopsy["EX_pydam_LPAREN_e_RPAREN_"]=3
conf_MeanCancerBiopsy["EX_pydx_LPAREN_e_RPAREN_"]=3
conf_MeanCancerBiopsy["EX_pydx5p_LPAREN_e_RPAREN_"]=3
conf_MeanCancerBiopsy["EX_pydxn_LPAREN_e_RPAREN_"]=3
conf_MeanCancerBiopsy["EX_pyr_LPAREN_e_RPAREN_"]=3
conf_MeanCancerBiopsy["EX_q10h2_LPAREN_e_RPAREN_"]=3
conf_MeanCancerBiopsy["EX_retfa_LPAREN_e_RPAREN_"]=3
conf_MeanCancerBiopsy["EX_retinol_LPAREN_e_RPAREN_"]=3
conf_MeanCancerBiopsy["EX_retn_LPAREN_e_RPAREN_"]=3
conf_MeanCancerBiopsy["EX_rib_D_LPAREN_e_RPAREN_"]=3
conf_MeanCancerBiopsy["EX_ribflv_LPAREN_e_RPAREN_"]=3
conf_MeanCancerBiopsy["EX_sbt_DASH_d_LPAREN_e_RPAREN_"]=3
conf_MeanCancerBiopsy["EX_sel_LPAREN_e_RPAREN_"]=3
conf_MeanCancerBiopsy["EX_ser_L_LPAREN_e_RPAREN_"]=3
conf_MeanCancerBiopsy["EX_so4_LPAREN_e_RPAREN_"]=3
conf_MeanCancerBiopsy["EX_sph1p_LPAREN_e_RPAREN_"]=3
conf_MeanCancerBiopsy["EX_spmd_LPAREN_e_RPAREN_"]=3
conf_MeanCancerBiopsy["EX_strch1_LPAREN_e_RPAREN_"]=3
conf_MeanCancerBiopsy["EX_strch2_LPAREN_e_RPAREN_"]=3
conf_MeanCancerBiopsy["EX_sucr_LPAREN_e_RPAREN_"]=3
conf_MeanCancerBiopsy["EX_tag_hs_LPAREN_e_RPAREN_"]=3
conf_MeanCancerBiopsy["EX_thm_LPAREN_e_RPAREN_"]=3
conf_MeanCancerBiopsy["EX_thr_L_LPAREN_e_RPAREN_"]=3
conf_MeanCancerBiopsy["EX_thymd_LPAREN_e_RPAREN_"]=3
conf_MeanCancerBiopsy["EX_tre_LPAREN_e_RPAREN_"]=3
conf_MeanCancerBiopsy["EX_trp_L_LPAREN_e_RPAREN_"]=3
conf_MeanCancerBiopsy["EX_ttdca_LPAREN_e_RPAREN_"]=3
conf_MeanCancerBiopsy["EX_tyr_L_LPAREN_e_RPAREN_"]=3
conf_MeanCancerBiopsy["EX_ura_LPAREN_e_RPAREN_"]=3
conf_MeanCancerBiopsy["EX_urea_LPAREN_e_RPAREN_"]=3
conf_MeanCancerBiopsy["EX_uri_LPAREN_e_RPAREN_"]=3
conf_MeanCancerBiopsy["EX_utp_LPAREN_e_RPAREN_"]=3
conf_MeanCancerBiopsy["EX_val_L_LPAREN_e_RPAREN_"]=3
conf_MeanCancerBiopsy["EX_xmp_LPAREN_e_RPAREN_"]=3
conf_MeanCancerBiopsy["EX_xyl_D_LPAREN_e_RPAREN_"]=3
conf_MeanCancerBiopsy["EX_xylt_LPAREN_e_RPAREN_"]=3


```


```python
conf_MeanNormalBiopsy["DM_4hrpo"]=3
conf_MeanNormalBiopsy["DM_datp_n_"]=3
conf_MeanNormalBiopsy["DM_dctp_n_"]=3
conf_MeanNormalBiopsy["DM_dgtp_n_"]=3
conf_MeanNormalBiopsy["DM_dttp_n_"]=3
conf_MeanNormalBiopsy["DM_Lcystin"]=3
conf_MeanNormalBiopsy["DM_pe_hs_LPAREN_r_RPAREN_"]=3
conf_MeanNormalBiopsy["EX_2hb_LPAREN_e_RPAREN_"]=3
conf_MeanNormalBiopsy["EX_34hpp_"]=3
conf_MeanNormalBiopsy["EX_3hpvs_LPAREN_e_RPAREN_"]=3
conf_MeanNormalBiopsy["EX_3mob_LPAREN_e_RPAREN_"]=3
conf_MeanNormalBiopsy["EX_4mop_LPAREN_e_RPAREN_"]=3
conf_MeanNormalBiopsy["EX_ac_LPAREN_e_RPAREN_"]=3
conf_MeanNormalBiopsy["EX_acac_LPAREN_e_RPAREN_"]=3
conf_MeanNormalBiopsy["EX_acald_LPAREN_e_RPAREN_"]=3
conf_MeanNormalBiopsy["EX_acetone_LPAREN_e_RPAREN_"]=3
conf_MeanNormalBiopsy["EX_acgam_LPAREN_e_RPAREN_"]=3
conf_MeanNormalBiopsy["EX_acmana_LPAREN_e_RPAREN_"]=3
conf_MeanNormalBiopsy["EX_ade_LPAREN_e_RPAREN_"]=3
conf_MeanNormalBiopsy["EX_adn_LPAREN_e_RPAREN_"]=3
conf_MeanNormalBiopsy["EX_adpcbl_LPAREN_e_RPAREN_"]=3
conf_MeanNormalBiopsy["EX_akg_LPAREN_e_RPAREN_"]=3
conf_MeanNormalBiopsy["EX_ala_B_LPAREN_e_RPAREN_"]=3
conf_MeanNormalBiopsy["EX_ala_D_LPAREN_e_RPAREN_"]=3
conf_MeanNormalBiopsy["EX_ala_L_LPAREN_e_RPAREN_"]=3
conf_MeanNormalBiopsy["EX_am9csa_LPAREN_e_RPAREN_"]=3
conf_MeanNormalBiopsy["EX_amp_LPAREN_e_RPAREN_"]=3
conf_MeanNormalBiopsy["EX_arab_L_LPAREN_e_RPAREN_"]=3
conf_MeanNormalBiopsy["EX_arachd_LPAREN_e_RPAREN_"]=3
conf_MeanNormalBiopsy["EX_arg_L_LPAREN_e_RPAREN_"]=3
conf_MeanNormalBiopsy["EX_asn_L_LPAREN_e_RPAREN_"]=3
conf_MeanNormalBiopsy["EX_asp_L_LPAREN_e_RPAREN_"]=3
conf_MeanNormalBiopsy["EX_atp_LPAREN_e_RPAREN_"]=3
conf_MeanNormalBiopsy["EX_bhb_LPAREN_e_RPAREN_"]=3
conf_MeanNormalBiopsy["EX_btn_LPAREN_e_RPAREN_"]=3
conf_MeanNormalBiopsy["EX_ca2_LPAREN_e_RPAREN_"]=3
conf_MeanNormalBiopsy["EX_carn_LPAREN_e_RPAREN_"]=3
conf_MeanNormalBiopsy["EX_caro_LPAREN_e_RPAREN_"]=3
conf_MeanNormalBiopsy["EX_cgly_LPAREN_e_RPAREN_"]=3
conf_MeanNormalBiopsy["EX_chol_LPAREN_e_RPAREN_"]=3
conf_MeanNormalBiopsy["EX_chsterol_LPAREN_e_RPAREN_"]=3
conf_MeanNormalBiopsy["EX_cit_LPAREN_e_RPAREN_"]=3
conf_MeanNormalBiopsy["EX_cl_LPAREN_e_RPAREN_"]=3
conf_MeanNormalBiopsy["EX_CLPND_LPAREN_e_RPAREN_"]=3
conf_MeanNormalBiopsy["EX_co2_LPAREN_e_RPAREN_"]=3
conf_MeanNormalBiopsy["EX_creat_LPAREN_e_RPAREN_"]=3
conf_MeanNormalBiopsy["EX_crn_LPAREN_e_RPAREN_"]=3
conf_MeanNormalBiopsy["EX_crvnc_LPAREN_e_RPAREN_"]=3
conf_MeanNormalBiopsy["EX_csa_LPAREN_e_RPAREN_"]=3
conf_MeanNormalBiopsy["EX_csn_LPAREN_e_RPAREN_"]=3
conf_MeanNormalBiopsy["EX_cys_L_LPAREN_e_RPAREN_"]=3
conf_MeanNormalBiopsy["EX_cytd_LPAREN_e_RPAREN_"]=3
conf_MeanNormalBiopsy["EX_dad_2_LPAREN_e_RPAREN_"]=3
conf_MeanNormalBiopsy["EX_dag_hs_LPAREN_e_RPAREN_"]=3
conf_MeanNormalBiopsy["EX_dcmp_LPAREN_e_RPAREN_"]=3
conf_MeanNormalBiopsy["EX_dcyt_LPAREN_e_RPAREN_"]=3
conf_MeanNormalBiopsy["EX_ddca_LPAREN_e_RPAREN_"]=3
conf_MeanNormalBiopsy["EX_dgsn_LPAREN_e_RPAREN_"]=3
conf_MeanNormalBiopsy["EX_dhdascb_LPAREN_e_RPAREN_"]=3
conf_MeanNormalBiopsy["EX_din_LPAREN_e_RPAREN_"]=3
conf_MeanNormalBiopsy["EX_dopa_LPAREN_e_RPAREN_"]=3
conf_MeanNormalBiopsy["EX_drib_LPAREN_e_RPAREN_"]=3
conf_MeanNormalBiopsy["EX_etoh_LPAREN_e_RPAREN_"]=3
conf_MeanNormalBiopsy["EX_fald_LPAREN_e_RPAREN_"]=3
conf_MeanNormalBiopsy["EX_fe2_LPAREN_e_RPAREN_"]=3
conf_MeanNormalBiopsy["EX_fe2_LPAREN_e_RPAREN_"]=3
conf_MeanNormalBiopsy["EX_fe3_LPAREN_e_RPAREN_"]=3
conf_MeanNormalBiopsy["EX_fmn_LPAREN_e_RPAREN_"]=3
conf_MeanNormalBiopsy["EX_fol_LPAREN_e_RPAREN_"]=3
conf_MeanNormalBiopsy["EX_for_LPAREN_e_RPAREN_"]=3
conf_MeanNormalBiopsy["EX_fru_LPAREN_e_RPAREN_"]=3
conf_MeanNormalBiopsy["EX_fuc_L_LPAREN_e_RPAREN_"]=3
conf_MeanNormalBiopsy["EX_fum_LPAREN_e_RPAREN_"]=3
conf_MeanNormalBiopsy["EX_gal_LPAREN_e_RPAREN_"]=3
conf_MeanNormalBiopsy["EX_gam_LPAREN_e_RPAREN_"]=3
conf_MeanNormalBiopsy["EX_gchola_LPAREN_e_RPAREN_"]=3
conf_MeanNormalBiopsy["EX_glc_LPAREN_e_RPAREN_"]=3
conf_MeanNormalBiopsy["EX_glcur_LPAREN_e_RPAREN_"]=3
conf_MeanNormalBiopsy["EX_gln_L_LPAREN_e_RPAREN_"]=3
conf_MeanNormalBiopsy["EX_glu_L_LPAREN_e_RPAREN_"]=3
conf_MeanNormalBiopsy["EX_gluala_LPAREN_e_RPAREN_"]=3
conf_MeanNormalBiopsy["EX_gly_LPAREN_e_RPAREN_"]=3
conf_MeanNormalBiopsy["EX_glyb_LPAREN_e_RPAREN_"]=3
conf_MeanNormalBiopsy["EX_glyc_LPAREN_e_RPAREN_"]=3
conf_MeanNormalBiopsy["EX_glyc3p_LPAREN_e_RPAREN_"]=3
conf_MeanNormalBiopsy["EX_glygn2_LPAREN_e_RPAREN_"]=3
conf_MeanNormalBiopsy["EX_gsn_LPAREN_e_RPAREN_"]=3
conf_MeanNormalBiopsy["EX_gthox_LPAREN_e_RPAREN_"]=3
conf_MeanNormalBiopsy["EX_gthrd_LPAREN_e_RPAREN_"]=3
conf_MeanNormalBiopsy["EX_gua_LPAREN_e_RPAREN_"]=3
conf_MeanNormalBiopsy["EX_h_LPAREN_e_RPAREN_"]=3
conf_MeanNormalBiopsy["EX_h2o_LPAREN_e_RPAREN_"]=3
conf_MeanNormalBiopsy["EX_ha_LPAREN_e_RPAREN_"]=3
conf_MeanNormalBiopsy["EX_HC00250_LPAREN_e_RPAREN_"]=3
conf_MeanNormalBiopsy["EX_HC01609_LPAREN_e_RPAREN_"]=3
conf_MeanNormalBiopsy["EX_HC01610_LPAREN_e_RPAREN_"]=3
conf_MeanNormalBiopsy["EX_hdca_LPAREN_e_RPAREN_"]=3
conf_MeanNormalBiopsy["EX_hdcea_LPAREN_e_RPAREN_"]=3
conf_MeanNormalBiopsy["EX_his_L_LPAREN_e_RPAREN_"]=3
conf_MeanNormalBiopsy["EX_hxan_LPAREN_e_RPAREN_"]=3
conf_MeanNormalBiopsy["EX_ile_L_LPAREN_e_RPAREN_"]=3
conf_MeanNormalBiopsy["EX_inost_LPAREN_e_RPAREN_"]=3
conf_MeanNormalBiopsy["EX_k_LPAREN_e_RPAREN_"]=3
conf_MeanNormalBiopsy["EX_lac_L_LPAREN_e_RPAREN_"]=3
conf_MeanNormalBiopsy["EX_lcts_LPAREN_e_RPAREN_"]=3
conf_MeanNormalBiopsy["EX_leu_L_LPAREN_e_RPAREN_"]=3
conf_MeanNormalBiopsy["EX_leuktrA4_LPAREN_e_RPAREN_"]=3
conf_MeanNormalBiopsy["EX_leuktrD4_LPAREN_e_RPAREN_"]=3
conf_MeanNormalBiopsy["EX_leuktrE4_LPAREN_e_RPAREN_"]=3
conf_MeanNormalBiopsy["EX_lnlc_LPAREN_e_RPAREN_"]=3
conf_MeanNormalBiopsy["EX_lnlnca_LPAREN_e_RPAREN_"]=3
conf_MeanNormalBiopsy["EX_lnlncg_LPAREN_e_RPAREN_"]=3
conf_MeanNormalBiopsy["EX_lpchol_hs_LPAREN_e_RPAREN_"]=3
conf_MeanNormalBiopsy["EX_lys_L_LPAREN_e_RPAREN_"]=3
conf_MeanNormalBiopsy["EX_mag_hs_LPAREN_e_RPAREN_"]=3
conf_MeanNormalBiopsy["EX_mal_L_LPAREN_e_RPAREN_"]=3
conf_MeanNormalBiopsy["EX_malt_LPAREN_e_RPAREN_"]=3
conf_MeanNormalBiopsy["EX_man_LPAREN_e_RPAREN_"]=3
conf_MeanNormalBiopsy["EX_meoh_LPAREN_e_RPAREN_"]=3
conf_MeanNormalBiopsy["EX_met_L_LPAREN_e_RPAREN_"]=3
conf_MeanNormalBiopsy["EX_n2m2nmasn_LPAREN_e_RPAREN_"]=3
conf_MeanNormalBiopsy["EX_na1_LPAREN_e_RPAREN_"]=3
conf_MeanNormalBiopsy["EX_nac_LPAREN_e_RPAREN_"]=3
conf_MeanNormalBiopsy["EX_ncam_LPAREN_e_RPAREN_"]=3
conf_MeanNormalBiopsy["EX_nh4_LPAREN_e_RPAREN_"]=3
conf_MeanNormalBiopsy["EX_no2_LPAREN_e_RPAREN_"]=3
conf_MeanNormalBiopsy["EX_o2_LPAREN_e_RPAREN_"]=3
conf_MeanNormalBiopsy["EX_o2s_LPAREN_e_RPAREN_"]=3
conf_MeanNormalBiopsy["EX_ocdca_LPAREN_e_RPAREN_"]=3
conf_MeanNormalBiopsy["EX_ocdcea_LPAREN_e_RPAREN_"]=3
conf_MeanNormalBiopsy["EX_octa_LPAREN_e_RPAREN_"]=3
conf_MeanNormalBiopsy["EX_orn_LPAREN_e_RPAREN_"]=3
conf_MeanNormalBiopsy["EX_oxa_LPAREN_e_RPAREN_"]=3
conf_MeanNormalBiopsy["EX_pe_hs_LPAREN_e_RPAREN_"]=3
conf_MeanNormalBiopsy["EX_pglyc_hs_LPAREN_e_RPAREN_"]=3
conf_MeanNormalBiopsy["EX_phe_L_LPAREN_e_RPAREN_"]=3
conf_MeanNormalBiopsy["EX_pheme_LPAREN_e_RPAREN_"]=3
conf_MeanNormalBiopsy["EX_pi_LPAREN_e_RPAREN_"]=3
conf_MeanNormalBiopsy["EX_pnto_R_LPAREN_e_RPAREN_"]=3
conf_MeanNormalBiopsy["EX_ppa_LPAREN_e_RPAREN_"]=3
conf_MeanNormalBiopsy["EX_pro_L_LPAREN_e_RPAREN_"]=3
conf_MeanNormalBiopsy["EX_prostgh2_LPAREN_e_RPAREN_"]=3
conf_MeanNormalBiopsy["EX_ps_hs_LPAREN_e_RPAREN_"]=3
conf_MeanNormalBiopsy["EX_ptrc_LPAREN_e_RPAREN_"]=3
conf_MeanNormalBiopsy["EX_ptvstlac_LPAREN_e_RPAREN_"]=3
conf_MeanNormalBiopsy["EX_pydam_LPAREN_e_RPAREN_"]=3
conf_MeanNormalBiopsy["EX_pydx_LPAREN_e_RPAREN_"]=3
conf_MeanNormalBiopsy["EX_pydx5p_LPAREN_e_RPAREN_"]=3
conf_MeanNormalBiopsy["EX_pydxn_LPAREN_e_RPAREN_"]=3
conf_MeanNormalBiopsy["EX_pyr_LPAREN_e_RPAREN_"]=3
conf_MeanNormalBiopsy["EX_q10h2_LPAREN_e_RPAREN_"]=3
conf_MeanNormalBiopsy["EX_retfa_LPAREN_e_RPAREN_"]=3
conf_MeanNormalBiopsy["EX_retinol_LPAREN_e_RPAREN_"]=3
conf_MeanNormalBiopsy["EX_retn_LPAREN_e_RPAREN_"]=3
conf_MeanNormalBiopsy["EX_rib_D_LPAREN_e_RPAREN_"]=3
conf_MeanNormalBiopsy["EX_ribflv_LPAREN_e_RPAREN_"]=3
conf_MeanNormalBiopsy["EX_sbt_DASH_d_LPAREN_e_RPAREN_"]=3
conf_MeanNormalBiopsy["EX_sel_LPAREN_e_RPAREN_"]=3
conf_MeanNormalBiopsy["EX_ser_L_LPAREN_e_RPAREN_"]=3
conf_MeanNormalBiopsy["EX_so4_LPAREN_e_RPAREN_"]=3
conf_MeanNormalBiopsy["EX_sph1p_LPAREN_e_RPAREN_"]=3
conf_MeanNormalBiopsy["EX_spmd_LPAREN_e_RPAREN_"]=3
conf_MeanNormalBiopsy["EX_strch1_LPAREN_e_RPAREN_"]=3
conf_MeanNormalBiopsy["EX_strch2_LPAREN_e_RPAREN_"]=3
conf_MeanNormalBiopsy["EX_sucr_LPAREN_e_RPAREN_"]=3
conf_MeanNormalBiopsy["EX_tag_hs_LPAREN_e_RPAREN_"]=3
conf_MeanNormalBiopsy["EX_thm_LPAREN_e_RPAREN_"]=3
conf_MeanNormalBiopsy["EX_thr_L_LPAREN_e_RPAREN_"]=3
conf_MeanNormalBiopsy["EX_thymd_LPAREN_e_RPAREN_"]=3
conf_MeanNormalBiopsy["EX_tre_LPAREN_e_RPAREN_"]=3
conf_MeanNormalBiopsy["EX_trp_L_LPAREN_e_RPAREN_"]=3
conf_MeanNormalBiopsy["EX_ttdca_LPAREN_e_RPAREN_"]=3
conf_MeanNormalBiopsy["EX_tyr_L_LPAREN_e_RPAREN_"]=3
conf_MeanNormalBiopsy["EX_ura_LPAREN_e_RPAREN_"]=3
conf_MeanNormalBiopsy["EX_urea_LPAREN_e_RPAREN_"]=3
conf_MeanNormalBiopsy["EX_uri_LPAREN_e_RPAREN_"]=3
conf_MeanNormalBiopsy["EX_utp_LPAREN_e_RPAREN_"]=3
conf_MeanNormalBiopsy["EX_val_L_LPAREN_e_RPAREN_"]=3
conf_MeanNormalBiopsy["EX_xmp_LPAREN_e_RPAREN_"]=3
conf_MeanNormalBiopsy["EX_xyl_D_LPAREN_e_RPAREN_"]=3
conf_MeanNormalBiopsy["EX_xylt_LPAREN_e_RPAREN_"]=3
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
```


```python
metas = ['0.014 biomass_DNA_c + 0.058 biomass_RNA_c + 0.071 biomass_carbohydrate_c + 0.097 biomass_lipid_c + 0.054 biomass_other_c + 0.706 biomass_protein_c --> ', '3pg_c', '4abut_c', '4hpro_LT_c', 'accoa_m', 'accoa_m --> coa_m', 'ade_c', 'adn_c', 'adp_c', 'akg_m', 'ala_B_c', 'ala_L_c', 'amet_c', 'amp_c', 'arg_L_c', 'asn_L_c', 'asp_D_c', 'asp_L_c', 'atp_c', 'bhb_c', 'cdp_c', 'CE1936_c', 'chol_c', 'chsterol_c', 'cit_c', 'citr_L_c', 'cmp_c', 'creat_c', 'crm_hs_c', 'crtn_c', 'ctp_c', 'cys_L_c', 'dag_hs_c', 'dhap_c', 'e4p_c', 'f6p_c', 'fdp_c', 'fru_c', 'fum_c', 'g1p_c', 'g3p_c', 'g6p_c', 'gdp_c', 'glc_D_c', 'glc_D_e', 'glc_D_e --> glc_D_c', 'gln_L_m', 'gln_L_c', 'gln_L_e --> gln_L_c', 'glu_L_c', 'glyb_c', 'gly_c', 'gmp_c', 'gthox_c', 'gthrd_c', 'gua_c', 'HC00342_c', 'his_L_c', 'hxan_c', 'icit_c', 'ile_L_c', 'lac_L_c', 'leu_L_c', 'leu_L_c',  'lys_L_c', 'mag_hs_c', 'mal_L_c', 'met_L_c', 'nad_c', 'nadh_c', 'nadh_m', 'nad_m', 'nadp_c', 'nadph_c', 'nadph_m', 'nadph_m', 'nadp_m', 'oaa_m', 'orn_c', 'pa_hs_c', 'pe_hs_c', 'pep_c', 'phe_L_c', 'pmtcoa_c --> coa_c', 'pro_D_c', 'pro_L_c', 'ps_hs_c', 'ptrc_c', 'pyr_c', 'pyr_m', 'r5p_c', 'ru5p_D_c', 's7p_c', 'ser_L_c', 'spmd_c', 'succ_c', 'succoa_m --> coa_m', 'tag_hs_c', 'thr_L_c', 'trp_L_c', 'tym_c', 'tyr_L_c', 'udp_c', 'ump_c', 'utp_c', 'val_L_c']
```


```python
for item in conf_MeanNormalBiopsy.items():
    print(item[0], "\t",item[1])
```

    r1516 	 3.0
    RE2594C 	 1
    BBHOX 	 3.0
    LVSTACOXD6MEhep 	 3.0
    PMTCOAFABP1tc 	 -1.0
    EX_ptvstm3_LPAREN_e_RPAREN_ 	 1
    RE1829M 	 -1.0
    r1866 	 -1.0
    G14T7g 	 1.0
    RE3144M 	 -1.0
    NACHEX26ly 	 1.0
    r0853 	 1
    r1592 	 3.0
    FA120ACPH 	 -1.0
    MM7B1g 	 3.0
    RE1050C 	 3.0
    ST3GAL21g 	 -1.0
    PPDOx 	 3.0
    HSD3A1r 	 -1.0
    r2316 	 3.0
    r0694 	 -1.0
    NDPK3m 	 -1.0
    2HBt2 	 1.0
    TMDOATtev 	 1
    EX_sl_L_LPAREN_e_RPAREN_ 	 1
    RE0702C 	 3.0
    EX_trp_L_LPAREN_e_RPAREN_ 	 3
    S6T11g 	 -1.0
    r2060 	 -1.0
    RE2972R 	 0.0
    RE1956R 	 0.0
    r2081 	 0.0
    ST3GAL23g 	 -1.0
    r2226 	 -1.0
    MCLACCYSR 	 1
    RE2988X 	 3.0
    r0695 	 -1.0
    r0591 	 0.0
    G14T4g 	 1.0
    EBASTINEtr 	 1
    C120CPT1 	 0
    SEBCOAPET 	 -1.0
    r1816 	 -1.0
    UDPGD 	 3.0
    RE2991X 	 3.0
    r0246 	 3.0
    r0962 	 1
    RE1448R 	 1
    r2179 	 -1.0
    GCNTg 	 -1.0
    FACOAL206 	 -1.0
    RE2349C 	 -1.0
    TMDM1hr 	 -1.0
    PI4PP 	 1
    S2TASE3ly 	 0.0
    RE3103R 	 3.0
    LVSTitr 	 1
    P4504B1r 	 -1.0
    LALDO2x 	 3.0
    r1964 	 -1.0
    FAOXC221C201x 	 -1.0
    NDPK5m 	 -1.0
    RNDR2 	 2.0
    EX_35dsmv_LPAREN_e_RPAREN_ 	 1
    r2002 	 -1.0
    1HIBUPGLUC_Sthv 	 1
    MM8Cg 	 3.0
    C6CRNtcx 	 -1.0
    EX_fvstet_LPAREN_e_RPAREN_ 	 1
    KHK 	 -1.0
    r1546 	 3.0
    RE1816R 	 3.0
    r0590 	 -1.0
    MALTly 	 -1.0
    EX_glygn2_LPAREN_e_RPAREN_ 	 3
    ESTSULT 	 1.0
    CHOLtn 	 1
    AM19CSteb 	 1
    GLNALANaEx 	 -1.0
    MI34PP 	 -1.0
    RE3534M 	 1.0
    CREATtmdiffir 	 1
    GluForTx 	 -1.0
    r2214 	 -1.0
    PROSTGE2t3 	 -1.0
    TLACFVShc 	 1
    7DHFtm 	 1
    ACMPGLUTtep 	 -1.0
    P4502C8 	 -1.0
    DSPVStev 	 1
    LNLNCGCPT2 	 -1.0
    5HXKYNOXDA 	 3.0
    VALATB0tc 	 1.0
    RE1919C 	 1
    ALAALACNc 	 3.0
    DHCR72r 	 0.0
    FAOXC121_3Em 	 1.0
    EX_psyl_LPAREN_e_RPAREN_ 	 1
    GACMTRc 	 -1.0
    r1576 	 3.0
    RE0565C 	 3.0
    RE3535R 	 1
    CYTK2n 	 3.0
    RE1818C 	 1
    AICART 	 1.0
    ADRNt 	 1
    MAGt 	 1
    LAPCOAe 	 1
    EX_sph1p_LPAREN_e_RPAREN_ 	 3
    r0666 	 -1.0
    GALNACT3g 	 0.0
    EX_nh4_LPAREN_e_RPAREN_ 	 3
    FAOXC170150m 	 3.0
    r1477 	 -1.0
    RE3038R 	 0.0
    RE3236C 	 1
    CRVSM23hc 	 -1.0
    AG13T12g 	 1.0
    RE1651C 	 1
    EX_HC02161_LPAREN_e_RPAREN_ 	 1
    DOCOSCOAtxc 	 -1.0
    SULPACMPtev 	 -1.0
    EX_appnn_LPAREN_e_RPAREN_ 	 1
    BILGLCURte 	 1.0
    EX_dmhptcrn_LPAREN_e_RPAREN_ 	 1
    r2521 	 1
    r0085 	 1.0
    O2tm 	 1
    ENMAN3g 	 1
    ACACT6p 	 3.0
    ADNCNT3tc 	 0.0
    RE3384M 	 3.0
    r2401 	 -1.0
    RE0689C 	 0.0
    RE3014R 	 -1.0
    RE1096M 	 -1.0
    RE2651R 	 1
    r1027 	 1
    RE2919M 	 3.0
    BTND1 	 -1.0
    GGH_7THFl 	 1.0
    EX_34dhoxpeg_LPAREN_e_RPAREN_ 	 1
    r0940 	 1
    EX_oxyp_LPAREN_e_RPAREN_ 	 1
    EHGLAT2m 	 3.0
    GP1CALPHAte 	 1
    RE2384C 	 1
    GALASE18ly 	 0
    EX_cyan_LPAREN_e_RPAREN_ 	 1
    RE3243C 	 1
    CRMPte 	 1
    BALAtmr 	 1
    r1401 	 1
    RE3286R 	 -1.0
    SO4OXAtex2 	 -1.0
    EX_dttp_LPAREN_e_RPAREN_ 	 1
    SPRMTDe 	 -1.0
    r0380 	 0
    EX_i_LPAREN_e_RPAREN_ 	 1
    RIBFLVte 	 -1.0
    EX_4nphsf_LPAREN_e_RPAREN_ 	 1
    EX_h2o2_LPAREN_e_RPAREN_ 	 1
    CH25H 	 -1.0
    FACOAL191 	 -1.0
    NADPtru 	 1
    RE2920X 	 -1.0
    AGDC 	 0
    RE3430C 	 1.0
    CRTSLtr 	 1
    r0082 	 2.0
    KYN 	 1.0
    GALSIDEtg 	 1
    r2069 	 -1.0
    RE1254C 	 1
    RE3525C 	 3.0
    EX_alaala_LPAREN_e_RPAREN_ 	 1
    r0537 	 3.0
    EX_56dhpvs_LPAREN_e_RPAREN_ 	 1
    EX_crvsm23_LPAREN_e_RPAREN_ 	 1
    EX_so4_LPAREN_e_RPAREN_ 	 3
    DESAT20_1 	 0
    r1292 	 1
    GLZitr 	 1
    PTRCAT1 	 3.0
    SEBACACT 	 1
    CGLYt3_2_ 	 -1.0
    CYTK1n 	 3.0
    ACP1_FMN_ 	 2.0
    RE3238C 	 1
    NTD3 	 3.0
    RE3153R 	 1
    RE3513N 	 2.0
    ACCOAgt 	 1.0
    HESTRATRIOLte 	 1
    NACHEXA8ly 	 -1.0
    DECDPtm 	 1
    RE3009C 	 0
    3HPVSTETtev 	 1
    r1752 	 -1.0
    fumt 	 1
    EX_ahdt_LPAREN_e_RPAREN_ 	 1
    FAOXC24C22x 	 -1.0
    HS3ly 	 1.0
    RE2974N 	 1
    PUNP2 	 2.0
    r0242 	 -1.0
    r2177 	 -1.0
    MCDm 	 -1.0
    FA182ACPH 	 -1.0
    RE3018C 	 1
    GALASE12ly 	 0
    56DHPVSitr 	 1
    5ADTSTSTERONESte 	 1
    C180CPT1 	 0
    DOLMANP_Lter 	 1
    GLZtd 	 1
    r1728 	 -1.0
    ASPPROASCT1 	 3.0
    r1012 	 1
    r1799 	 -1.0
    ST8SIA53g 	 -1.0
    TTDCPT1 	 0
    AMACR2p 	 -1.0
    SPODMn 	 3.0
    FUMAC 	 -1.0
    THRATB0tc 	 1.0
    AM1C9CSitr 	 1
    ORETNF2 	 1
    r1687 	 -1.0
    APNNOXte 	 1
    biomass_carbohydrate 	 3
    EX_crvsm1_LPAREN_e_RPAREN_ 	 1
    EX_ins_LPAREN_e_RPAREN_ 	 1
    r2015 	 -1.0
    r2066 	 -1.0
    INSKm 	 1
    CYTDt 	 -1.0
    MALSO4tm 	 -1.0
    10FTHF5GLUtl 	 1
    r1613 	 -1.0
    r1487 	 2.0
    RTOTALt 	 1
    RE3261C 	 -1.0
    DIGALSIDEtl 	 1
    AMACRr 	 -1.0
    FATP8t 	 -1.0
    r1614 	 -1.0
    DESAT18_7 	 -1.0
    r0604 	 3.0
    RE3112R 	 1
    CRTNsyn 	 1
    RE3264R 	 0
    r1450 	 3.0
    TAURt4_2_r 	 0.0
    G14T10g 	 1.0
    XANDp 	 -1.0
    RE3464R 	 -1.0
    RE3110R 	 3.0
    TSTSTERONEt 	 1
    r1661 	 3.0
    RE3251M 	 -1.0
    C2M26DCOAHLm 	 3.0
    RE1925C 	 1
    RE3245C 	 1
    r1666 	 -1.0
    EX_13dmt_LPAREN_e_RPAREN_ 	 1
    CYSTGL 	 -1.0
    SMVACIDitr 	 1
    GALGLUSIDEtl 	 1
    FALDtly 	 1
    PACCOAL 	 1
    C4x 	 3.0
    EX_acac_LPAREN_e_RPAREN_ 	 3
    RE0918G 	 -1.0
    r0791 	 -1.0
    GLY3Pt 	 1
    DATPtm 	 1
    RE0937E 	 0.0
    CBASPte 	 1
    r0915 	 -1.0
    EX_urea_LPAREN_e_RPAREN_ 	 3
    CYTK6 	 3.0
    r1530 	 0.0
    PGM 	 3.0
    LPS3 	 3.0
    NAGLCAly 	 1
    EX_pheme_LPAREN_e_RPAREN_ 	 3
    r1726 	 -1.0
    6BHGLZitr 	 1
    EX_prostgf2_LPAREN_e_RPAREN_ 	 1
    FBA4 	 3.0
    BALAVECSEC 	 -1.0
    RE1099G 	 -1.0
    HSD17B2r 	 -1.0
    r1568 	 3.0
    EICOSTETCRNt 	 -1.0
    NACHEXA13ly 	 -1.0
    r2344 	 0
    VALt5m 	 1
    r2243 	 -1.0
    r1451 	 3.0
    r1880 	 -1.0
    ACGAGBSIDEtl 	 1
    RE0864C 	 1
    DMNONCRNCPT2 	 1
    r2032 	 -1.0
    EX_no_LPAREN_e_RPAREN_ 	 1
    r0765 	 -1.0
    FUCASEe 	 -1.0
    KSII_CORE4tly 	 1
    DCSPTN1COAtxc 	 -1.0
    PGL 	 3.0
    DHEAStr 	 1
    sink_citr_LPAREN_c_RPAREN_ 	 1
    LSTNM2hr 	 3.0
    RE1904C 	 1
    DUDPtn 	 1
    OCTDECE1CRNe 	 1
    FAOXC143C123m 	 3.0
    r1020 	 1
    RE2248C 	 1
    ADNK1m 	 1
    TDCRNe 	 1
    RE3460C 	 1
    r2219 	 -1.0
    HIVCRNe 	 1
    GPAMm_hs 	 2.0
    RE1630R 	 1
    RE1915C 	 1
    RE3537C 	 1
    r0652 	 -1.0
    LACZly 	 0
    RE2151R 	 0
    NACHEX7ly 	 1.0
    NTD2 	 3.0
    FAOXC163GC142m 	 3.0
    RE2541E 	 -1.0
    ABO4g 	 -1.0
    r1625 	 -1.0
    r1316 	 -1.0
    PAIL4P_HStn 	 1
    AM4N9CSitr 	 1
    ELAIDCPT2 	 -1.0
    THRFVShc 	 1
    r0669 	 -1.0
    HSD17B3r 	 -1.0
    RE0580L 	 1
    r0753 	 2.0
    NTD7 	 3.0
    RBK_D 	 -1.0
    r2157 	 -1.0
    FUCGALGBSIDEtg 	 1
    ILETAm 	 -1.0
    2HIBUPGLUC_Sthv 	 1
    ARTPLM3 	 1
    r2172 	 -1.0
    EX_crtsl_LPAREN_e_RPAREN_ 	 1
    ETF 	 1.0
    FAOXC163Gm 	 1.0
    FAOXC16C16OHm 	 3.0
    r0830 	 -1.0
    RE3230R 	 1
    CBLtle 	 1
    FAOXMC10OHMC10r 	 -1.0
    OCTDECCPT2 	 -1.0
    COKECBESr 	 1.0
    EX_uri_LPAREN_e_RPAREN_ 	 3
    RE2154C 	 1
    LSTNM7itr 	 1
    LS3 	 1
    ACNACNGALGBSIDEte 	 1
    r1536 	 -1.0
    RE0908R 	 -1.0
    HMGCOAtm 	 1
    RE2377C 	 1
    r1970 	 -1.0
    GLCNACASE1ly 	 -1.0
    NDPK7m 	 -1.0
    NDPK9m 	 -1.0
    MI3PP 	 3.0
    NDPK6 	 0
    NDP6 	 3.0
    CRVSM23hr 	 3.0
    CBPPer 	 1.0
    DUTPDP 	 1
    LEUPHELAT2tc 	 1
    AM1C4N9CSteb 	 1
    OCTDEC2ACBP 	 3.0
    r1629 	 3.0
    RSVGLUChc 	 1
    EX_npthl_LPAREN_e_RPAREN_ 	 1
    ALDD20x 	 3.0
    HOXG 	 0.0
    AGLPR 	 1
    RE2677E 	 1
    ACGALFUCGALACGALFUC12GAL14ACGLCGALGLUSIDEtg 	 1
    4MTOLBUTAMIDEte 	 1
    r1876 	 -1.0
    SELCYSLY2 	 -1.0
    RE1804M 	 -1.0
    EX_g1p_LPAREN_e_RPAREN_ 	 1
    r1743 	 -1.0
    RE3141X 	 3.0
    EX_thyox_L_LPAREN_e_RPAREN_ 	 1
    r1844 	 -1.0
    S6TASE4ly 	 -1.0
    CMPACNAtn 	 1
    RE1815M 	 1.0
    PECTCHLe 	 1
    r1431 	 -1.0
    r2045 	 -1.0
    GTPCI 	 2.0
    EX_sarcs_LPAREN_e_RPAREN_ 	 1
    FAOXC204184m2 	 3.0
    CSDPASEly 	 1
    r2012 	 -1.0
    P4508B11r 	 -1.0
    NMNS 	 2.0
    EX_lst4exp_LPAREN_e_RPAREN_ 	 1
    GLYLEUPEPT1tc 	 -1.0
    r0584 	 1
    ACACT1r 	 -1.0
    SPHMYLNtl 	 1
    RE1441R 	 1
    r0921 	 1
    FUCACNGAL14ACGLCGALGLUSIDEtg 	 1
    CARIBUP_Sthv 	 1
    GDPFUCtg 	 -1.0
    RE2862C 	 1
    ECOAH9m 	 3.0
    FAOXC15ATPx 	 -1.0
    CYSAMPtev 	 1
    GMPtn 	 1
    TDPm 	 1
    r0283 	 -1.0
    r0510 	 3.0
    GLCNACT1g 	 0.0
    CYSSERNaEx 	 3.0
    CRGLZtev 	 1
    HEXCCPT2 	 -1.0
    EX_3hpvstet_LPAREN_e_RPAREN_ 	 1
    C50CPT1 	 0
    RE3146R 	 2.0
    NDP7g 	 3.0
    r1291 	 1
    r1894 	 -1.0
    RE2920M 	 3.0
    RE3160C 	 1
    EX_glgchlo_LPAREN_e_RPAREN_ 	 1
    NADPtxu 	 1
    RE3301G 	 3.0
    G14T19g 	 1.0
    A_MANASE 	 -1.0
    r2491 	 -1.0
    MECOALm 	 1
    r1700 	 -1.0
    r1741 	 -1.0
    RE3179M 	 3.0
    S6TASE16ly 	 3.0
    NDP7er 	 3.0
    FATP1t 	 1.0
    DOCOSDIACTD 	 1
    PROSTGI2tr 	 1
    EX_nrpphr_LPAREN_e_RPAREN_ 	 1
    EX_bildglcur_LPAREN_e_RPAREN_ 	 1
    OXAtp 	 1
    RE2909C 	 1
    RE1532X 	 3.0
    TETPENT6t 	 1
    RE3397M 	 1
    FUC14GALACGLCGALGLUSIDEte 	 1
    RE1835M 	 1.0
    r1997 	 -1.0
    r1766 	 -1.0
    GLC3MEACPtev 	 1
    CHSTEROLSULT 	 3.0
    UGALNACter 	 0.0
    FACOAL150 	 2.0
    C3DCe 	 1
    r0838 	 1
    XOLTRI27tc 	 1
    RE0512C 	 1
    r0221 	 0.0
    RE3120R 	 3.0
    r0365 	 1
    LYStip 	 1
    EX_4hpro_LPAREN_e_RPAREN_ 	 1
    r2112 	 0.0
    r2103 	 1
    RE3562X 	 2.0
    GLCNACT_U 	 1
    DGTPtn 	 1
    RE2865C 	 1
    EX_31dmt_LPAREN_e_RPAREN_ 	 1
    PVSATPtu 	 0.0
    EX_peplys_LPAREN_e_RPAREN_ 	 1
    r2399 	 -1.0
    RE3342M 	 3.0
    r1751 	 -1.0
    EX_octa_LPAREN_e_RPAREN_ 	 3
    OCDCAFATPtc 	 1.0
    r1500 	 -1.0
    RE3153C 	 1
    r2020 	 -1.0
    EX_dtmp_LPAREN_e_RPAREN_ 	 1
    PTVSTitr 	 1
    GLCNACASE4ly 	 -1.0
    FUT91g 	 -1.0
    NaKt 	 3.0
    RE3079X 	 -1.0
    r1901 	 -1.0
    IPDPtx 	 1
    ITCOAL1m 	 1.0
    MTAP 	 -1.0
    VITD2Hm 	 -1.0
    NTD1 	 3.0
    r1863 	 -1.0
    RE0383C 	 3.0
    DHCR71r 	 0.0
    r1798 	 -1.0
    ORETNtn2 	 1
    S6TASE10ly 	 -1.0
    PTRCOX1 	 0.0
    1331TAALThr 	 1
    S6TASE1ly 	 3.0
    FAOXC163_4Z_7Z_10Zm 	 1.0
    r1260 	 2.0
    r0055 	 1
    r0268 	 1
    XSERtg 	 1
    CEPTE 	 -1.0
    EX_1531tacr_LPAREN_e_RPAREN_ 	 1
    NTPP10 	 3.0
    NDPK1m 	 -1.0
    CYSSNAT5tc 	 -1.0
    GASNASE3ly 	 1.0
    RE3108C 	 1
    PI34P4Pn 	 1
    PFK 	 3.0
    CRVSATPthc 	 1
    DUMPtn 	 1
    FRDPtr 	 1
    UGT1A6r 	 0
    DNDPt9m 	 -1.0
    ATVLAChc 	 1
    KYNAKGAT 	 -1.0
    DM_4hrpo 	 3
    biomass_other 	 3
    ST8SIA56g 	 -1.0
    r0936 	 1
    ACONT 	 0.0
    TREH 	 -1.0
    1513TACRtep 	 1
    BACCLm 	 -1.0
    r0767 	 -1.0
    C10DCCACT 	 -1.0
    NOS2 	 -1.0
    PA_HStn 	 1
    KYN3OX 	 -1.0
    r0715 	 3.0
    RAHY 	 1
    FAOXTC142TC122m 	 3.0
    GLCNACT4g 	 -1.0
    RE3493C 	 1
    C182OHc 	 0
    FAOXC183806x 	 0.0
    r0763 	 -1.0
    RE0925C 	 -1.0
    6EPVShc 	 1
    5THFtl 	 1
    DCATDr 	 1
    r2488 	 -1.0
    DOLPGT2_Uer 	 1
    TETTET6t 	 1
    RE3436C 	 1
    NACHEXA10ly 	 -1.0
    RE3003M 	 1.0
    RSVhc 	 1
    RE1906C 	 1
    r2224 	 -1.0
    CSND 	 1
    RE2151C 	 -1.0
    TDPDRR 	 1
    PI4PLCn 	 -1.0
    DMHPTCRNCPT2 	 1
    HMGCOASim 	 -1.0
    AG13T4g 	 1.0
    r0781 	 0
    GASNASEly 	 1.0
    EX_maltttr_LPAREN_e_RPAREN_ 	 1
    MI1345PKn 	 -1.0
    RE2513N 	 -1.0
    r2404 	 -1.0
    r1992 	 -1.0
    r2043 	 -1.0
    r2322 	 3.0
    1331TACRhr 	 1
    RE1635X 	 1.0
    RE2853C 	 1
    RN0032R 	 -1.0
    r1375 	 1
    ME1m 	 -1.0
    GPIMTer_U 	 1
    r0149 	 3.0
    SO4CLtex2 	 -1.0
    HSAT1ly 	 1
    PI34P5K 	 -1.0
    C6COAt 	 3.0
    EX_drib_LPAREN_e_RPAREN_ 	 3
    r1986 	 -1.0
    SUBERICACT 	 1
    24_25VITD3Hm 	 -1.0
    CTPtn 	 1
    CRVSATPtu 	 0.0
    EX_limnen_LPAREN_e_RPAREN_ 	 1
    OXYPthc 	 -1.0
    PCHOL_HStg 	 1
    FAOXC150130m 	 3.0
    RE1587C 	 1
    RE3521M 	 1.0
    RE2319M 	 1.0
    EX_gdp_LPAREN_e_RPAREN_ 	 1
    ACMPGLUTthc 	 -1.0
    EX_am1cglc_LPAREN_e_RPAREN_ 	 1
    ACGALFUCGALACGALFUC12GAL14ACGLCGALGLUSIDEte 	 1
    FAOXTC122m 	 1.0
    L_LACt2r 	 1.0
    NTD8l 	 2.0
    RE3502X 	 0
    r2005 	 -1.0
    r1813 	 -1.0
    RE3596M 	 0
    EX_gam_LPAREN_e_RPAREN_ 	 3
    r0741 	 1
    r0720 	 -1.0
    SIAASEly 	 -1.0
    RE3470X 	 0
    MM6ag 	 0.0
    RE1519X 	 1
    RAI1 	 1
    r1328 	 1
    COQ5m 	 1
    RE3596X 	 0
    RE2251C 	 1
    RE1711M 	 -1.0
    RE3161R 	 3.0
    PTDCACRNt 	 -1.0
    ACGALtlg 	 1
    S6T18g 	 1.0
    FUCACNGALACGLCGALGLUSIDEte 	 1
    r1758 	 -1.0
    ACGAGBSIDEtg 	 1
    r1380 	 3.0
    RE3520M 	 1
    GLZABCteb 	 1
    CYSTLEUrBATtc 	 -1.0
    LSTO2r 	 2.0
    BILDGLCURtr 	 1
    r2116 	 0.0
    RE3559M 	 3.0
    ARGB0AT3tc 	 -1.0
    ACGSm 	 -1.0
    NAt3_1 	 0
    RE3587N 	 3.0
    r1660 	 3.0
    PPCOACm 	 -1.0
    NDERSVteb 	 0.0
    r1763 	 -1.0
    NTD6 	 3.0
    P5CRm 	 -1.0
    M4MPDOL_Uter 	 1
    RE2917X 	 3.0
    r0450 	 -1.0
    FATP9t 	 -1.0
    RE2759X 	 -1.0
    OCBTm 	 -1.0
    4HATVLACtep 	 1
    TIGCRNe 	 1
    RE3550X 	 1
    r0527 	 1
    EX_h_LPAREN_e_RPAREN_ 	 3
    r2385 	 -1.0
    r0722 	 3.0
    GF6PTA 	 0.0
    BILGLCURt 	 -1.0
    5ADTSTSTERONEtr 	 1
    DSPVSitr 	 1
    CDS 	 3.0
    GUMGCHLe 	 1
    G1M7MASNBterg 	 1
    r2315 	 3.0
    GLCAASE5ly 	 3.0
    TMDM3hr 	 1
    r0474 	 -1.0
    MHGLZhr 	 3.0
    DGCHOLtx 	 1
    GMPtg 	 1
    UGT1A5r2 	 0
    AM1ALCSitr 	 1
    6HTSTSTERONEtr 	 1
    FAOXC122_3E_6Em 	 1.0
    EX_imp_LPAREN_e_RPAREN_ 	 1
    RE3506R 	 -1.0
    PGLYCP 	 -1.0
    C9BRxtc 	 -1.0
    RE0935C 	 1
    ACN13ACNGALGBSIDEtg 	 1
    ST6GALNAC25 	 3.0
    EICOSTETCPT2 	 -1.0
    4HDEBRISOQUINEte 	 1
    r2130 	 1
    H8MTer_L 	 -1.0
    UMPK 	 3.0
    6HLVSTACIDitr 	 1
    RE1941R 	 0.0
    r0383 	 3.0
    ADNCYC 	 1.0
    GACPAILter 	 1
    6HSMVhep 	 3.0
    RE0827E 	 -1.0
    r2328 	 3.0
    LNLNCAt 	 -1.0
    r0571 	 1
    r2029 	 -1.0
    RE3189M 	 3.0
    r0480 	 1
    RE0574C 	 1
    RE1522M 	 3.0
    GSNKm 	 1
    r2175 	 -1.0
    PAFHe 	 -1.0
    ALA_DASH_DTDe 	 1
    FAOXC164_4Z_7Z_10Z_13Zm 	 1.0
    MALtm 	 -1.0
    HXt 	 1
    4ABUTtm 	 1
    FACOAL60i 	 1
    RE3308C 	 1
    r1683 	 -1.0
    34HPLFM 	 1
    G14T9g 	 1.0
    RE3554R 	 3.0
    r1947 	 -1.0
    KHK2 	 -1.0
    G14T12g 	 1.0
    KHK3 	 -1.0
    EX_gal_LPAREN_e_RPAREN_ 	 3
    METyLATthc 	 -1.0
    GLCAT7g 	 -1.0
    r1183 	 1.0
    r0693 	 -1.0
    4HPRO_LTte 	 1
    M4ATAer 	 -1.0
    BTNtn 	 1
    C4CRNCPT2 	 -1.0
    TMDOATPtsc 	 1
    r2021 	 -1.0
    RE3016R 	 1
    r1711 	 -1.0
    NACHEXA1ly 	 -1.0
    PCm 	 -1.0
    GLRASE 	 1
    r1019 	 1
    CYSTALArBATtc 	 -1.0
    r0648 	 3.0
    GALFUCGALACGLCGAL14ACGLCGALGLUSIDEtg 	 1
    RE2636R 	 -1.0
    EX_decdicrn_ 	 1
    EX_c101crn_ 	 1
    2OXOADPTm 	 -1.0
    NTD11 	 3.0
    PMI12346PHn 	 1.0
    FPGS8m 	 -1.0
    EX_ddecrn_ 	 1
    r1071 	 1
    FATP3t 	 -1.0
    CHSTEROLt3 	 1
    GLUPROASCT1 	 3.0
    FAOXC142_5Z_8Zm 	 3.0
    r0963 	 3.0
    r1936 	 -1.0
    RE2514L 	 3.0
    3SALATAim 	 3.0
    GLCtg 	 3.0
    GGH_5THFe 	 1.0
    GTHRDtr 	 1
    BMTer_U 	 1
    EX_3mop_LPAREN_e_RPAREN_ 	 1
    AM9CSAitr 	 1
    RE3346M 	 3.0
    NDP3ex 	 3.0
    G14Tg 	 1.0
    MGSA2 	 1
    THFt2 	 -1.0
    ATVETHGLUChc 	 0
    AMCOXO 	 1
    EX_lpchol_hs_LPAREN_e_RPAREN_ 	 3
    r1561 	 3.0
    PROt4_2_r 	 -1.0
    CYSTGLUex 	 -1.0
    TAGAT_Dt 	 1
    RE1520X 	 -1.0
    MERACMPtep 	 -1.0
    TKT2 	 1.0
    RE3237C 	 1
    RE3557R 	 3.0
    r0926 	 1
    SFGTH 	 3.0
    RE3367X 	 -1.0
    RE2649X 	 -1.0
    RE2523X 	 1.0
    r1836 	 -1.0
    ACALDt 	 1
    UGALNACtg 	 -1.0
    r2149 	 -1.0
    RE1635R 	 3.0
    LEUKTRD4t 	 1
    EX_5adtststerone_LPAREN_e_RPAREN_ 	 1
    S6T3g 	 1.0
    r2152 	 -1.0
    r2024 	 -1.0
    H2Oter 	 1
    SMVGLUChep 	 0
    RE2813C 	 1
    RE3526C 	 3.0
    CLPNDCPT2 	 -1.0
    ACCOAtr 	 1.0
    RE3181C 	 1
    r2490 	 -1.0
    42A12BOOX 	 3.0
    PI4PLC 	 3.0
    EX_am19cs_LPAREN_e_RPAREN_ 	 1
    r2400 	 -1.0
    B_MANNASEly 	 -1.0
    LYStn 	 1
    RE2127C 	 1
    DOCOSACTDr 	 1
    TYRCBOX 	 -1.0
    NCAMDe 	 1
    r2091 	 0.0
    r1029 	 1
    SERGLNexR 	 -1.0
    FAOXC9BRC7BRm 	 1
    RE1905R 	 -1.0
    sink_octdececoa_LPAREN_c_RPAREN_ 	 1
    RE3339C 	 0
    INSt5le 	 0.0
    CORE3GTg 	 -1.0
    RE2427M 	 -1.0
    r0424 	 2.0
    EX_s3meacmp_LPAREN_e_RPAREN_ 	 1
    DHCR242r 	 3.0
    HEXDICOAACBPx 	 3.0
    r2391 	 -1.0
    RE1236C 	 1
    G14T11g 	 1.0
    r1911 	 -1.0
    ST8SIA52g 	 -1.0
    S6T1g 	 -1.0
    DOLICHOL_Uter 	 1
    r2188 	 -1.0
    ADCim 	 1
    1513DTALThr 	 1
    G14T18g 	 1.0
    RE3176R 	 1
    PTPAT 	 -1.0
    RE3038X 	 0.0
    6MELVSTthep 	 1
    ITPtn 	 1
    r1522 	 3.0
    DCIm 	 1.0
    GSNt5le 	 0.0
    RE2622R 	 0
    S3TASE3ly 	 1
    KDNH 	 1
    Htx 	 1
    VALB0AT3tc 	 -1.0
    RN0023C 	 0.0
    GALU 	 3.0
    ADK3 	 -1.0
    SGPL12r 	 3.0
    r0907 	 1
    r0859 	 1
    DNDPt13m 	 -1.0
    PIt2m 	 3.0
    CYSGLNNaEx 	 -1.0
    PPCOAOm 	 2.0
    CYANtm 	 1
    34DHXMANDACOX 	 1.0
    r1996 	 -1.0
    PRDX 	 3.0
    GGH_10FTHF7GLUl 	 1.0
    HMBS 	 -1.0
    EX_ddca_LPAREN_e_RPAREN_ 	 3
    EX_isolvstacid_LPAREN_e_RPAREN_ 	 1
    NACHEXA22ly 	 -1.0
    r0691 	 -1.0
    FAOXC2051843x 	 0.0
    EX_cbasp_LPAREN_e_RPAREN_ 	 1
    GAMYe 	 -1.0
    r1732 	 -1.0
    GLCAASE6ly 	 3.0
    EX_amp_LPAREN_e_RPAREN_ 	 3
    RE1631C 	 1
    S3T2g 	 -1.0
    RE1818X 	 1.0
    NTD5m 	 -1.0
    r2129 	 0.0
    r0047 	 3.0
    FAOXC163C143m 	 3.0
    EBP1r 	 2.0
    RE3533C 	 1
    r1860 	 -1.0
    RE3560X 	 -1.0
    NRPPHRtu 	 2.0
    TPI 	 3.0
    RE3435R 	 2.0
    r1984 	 -1.0
    PAFt 	 1
    RE2514N 	 -1.0
    ADEt 	 -1.0
    EX_oxy1rb_LPAREN_e_RPAREN_ 	 1
    XOL7AONEtr 	 1
    RE3452C 	 1
    r2384 	 -1.0
    FAOXC184_6Z_9Z_12Z_15Zx 	 0.0
    FALDH 	 3.0
    r1991 	 -1.0
    TYRASE 	 -1.0
    NDPK1n 	 0
    RE1816M 	 1.0
    EX_gchola_LPAREN_e_RPAREN_ 	 3
    r2194 	 -1.0
    r0021 	 3.0
    r0434 	 1
    ADMDC 	 2.0
    NICRNS 	 1
    3SPYRSP 	 1
    r0186 	 -1.0
    EX_pvsgluc_LPAREN_e_RPAREN_ 	 1
    r2191 	 -1.0
    DOLPGT3_Uer 	 1
    EX_xoltri27_LPAREN_e_RPAREN_ 	 1
    CRVSM22hc 	 1
    DGK1 	 1.0
    ACOAHi 	 -1.0
    GLCNACASE3ly 	 -1.0
    EX_leuktrD4_LPAREN_e_RPAREN_ 	 3
    r2413 	 -1.0
    RDH1a 	 3.0
    FAOXC226m 	 1.0
    RN0031C 	 0.0
    EX_malttr_LPAREN_e_RPAREN_ 	 1
    r2447 	 0.0
    r1835 	 -1.0
    NACHEX6ly 	 1.0
    3SALATAi 	 1.0
    r1610 	 3.0
    ELAIDCRNt 	 -1.0
    ATVLACGLCURhc 	 0
    URIDK2m 	 -1.0
    RE2525C 	 1
    RE0958E 	 0
    r1723 	 -1.0
    RE2908X 	 3.0
    r2526 	 -1.0
    NAtx 	 1
    EX_adrnl_LPAREN_e_RPAREN_ 	 1
    r2407 	 -1.0
    DIGALSGALSIDEte 	 1
    RE3444X 	 -1.0
    r1912 	 -1.0
    r1557 	 3.0
    r1856 	 -1.0
    DNDPt46m 	 -1.0
    NAGAly 	 -1.0
    G14T14g 	 1.0
    ANDRSTRNGLCtr 	 1
    3HKYNAKGAT 	 -1.0
    FAOXC101C102x 	 0.0
    r2506 	 1
    r2517 	 2.0
    APOC_LYS_BTNP 	 1
    ASPDxt 	 1
    r1818 	 -1.0
    CYOOm3 	 -1.0
    MHGLZABCt 	 1
    4HATVACIDteb 	 -1.0
    RE3229R 	 1
    NACHEX17ly 	 1.0
    DNDPt43m 	 -1.0
    ATVLAChr 	 1
    EBASTINEOHtr 	 1
    RE3000M 	 2.0
    TETDEC2CRNe 	 1
    r1146 	 1
    RE2666G 	 -1.0
    TAURIBUP_Sthv 	 1
    LNLNCGCRNt 	 -1.0
    STACMPtev 	 1
    RE3534C 	 1
    r0430 	 -1.0
    H7_TAer 	 -1.0
    HISt4 	 3.0
    RE3454C 	 1
    DEDOLR_U 	 1
    r0620 	 -1.0
    EX_strdnc_LPAREN_e_RPAREN_ 	 1
    FVSTETitr 	 1
    EX_HC02216_LPAREN_e_RPAREN_ 	 1
    FAOXC130110m 	 3.0
    EX_arachcoa_LPAREN_e_RPAREN_ 	 1
    GALGT4 	 -1.0
    NACHEX23ly 	 1.0
    RE3286C 	 1
    B3GNT312g 	 -1.0
    r1045 	 0.0
    FAOXCPRIST3x 	 -1.0
    NACHEXA4ly 	 -1.0
    DPHMBDCm 	 1
    THRPHELAT2tc 	 1
    GLCNACASE5ly 	 -1.0
    EX_13_cis_retnglc_LPAREN_e_RPAREN_ 	 1
    OMHDEACIDTD 	 1
    r1006 	 1
    EX_tetdece1crn_ 	 1
    C60CPT1 	 0
    RE0577X 	 -1.0
    r1013 	 1
    OMPDC 	 1.0
    OCD11CRNCACT 	 1
    RE1573M 	 1.0
    EX_ser_L_LPAREN_e_RPAREN_ 	 3
    EX_pvs_LPAREN_e_RPAREN_ 	 1
    DM_hretn_n_ 	 1
    RETH2e 	 1
    XAOx 	 -1.0
    DM_core5_g_ 	 1
    MAN1_7Ber 	 -1.0
    EX_crn_LPAREN_e_RPAREN_ 	 3
    r2156 	 -1.0
    RE3342X 	 3.0
    RE2146C 	 -1.0
    G14T17g 	 1.0
    EX_25hvitd2_LPAREN_e_RPAREN_ 	 1
    EX_mag_hs_LPAREN_e_RPAREN_ 	 3
    r0747 	 -1.0
    RE3562R 	 2.0
    S3TASE1ly 	 1
    RE3287R 	 -1.0
    FACOAL161 	 2.0
    FVStep 	 1
    DMATTx 	 2.0
    WHTTDCAte 	 1
    PRISTCOAtx 	 1
    sink_5hpet_LPAREN_c_RPAREN_ 	 1
    O2Stm 	 1
    ACACT5p 	 3.0
    TCYNTtm 	 1
    BTND1n 	 -1.0
    D_LACtm 	 0.0
    DNDPt37m 	 -1.0
    r2232 	 -1.0
    r2378 	 -1.0
    DCK2n 	 -1.0
    GQ1BALPHAte 	 1
    H2CO3D 	 3.0
    DITPtn 	 1
    FUT14g 	 -1.0
    NAt 	 3.0
    DNDPt60m 	 -1.0
    EX_tststeroneglc_LPAREN_e_RPAREN_ 	 1
    PYAM5Ptm 	 1
    RE3496C 	 1
    RE3159X 	 -1.0
    6OHFVSGLUitr 	 1
    EX_q10h2_LPAREN_e_RPAREN_ 	 3
    PRPPS 	 1.0
    C16DCe 	 1
    r0546 	 2.0
    AKR1C42 	 -1.0
    OCCOAtm 	 1
    EX_crvsm24_LPAREN_e_RPAREN_ 	 1
    r1784 	 -1.0
    RE3089X 	 3.0
    r2160 	 -1.0
    r1492 	 -1.0
    CSPG_Dtly 	 1
    GLYt7_311_r 	 -1.0
    FACOAL1832 	 2.0
    5ADTSTSTERONEte 	 1
    NTD12 	 -1.0
    RE1916X 	 1.0
    RE1835X 	 -1.0
    CDIPTr 	 3.0
    RE3093X 	 3.0
    PUNP1 	 2.0
    r0193 	 -1.0
    r1621 	 -1.0
    UDPGLCAter 	 0.0
    SARCStm 	 1
    r1921 	 -1.0
    4HATVLACthc 	 -1.0
    r2217 	 -1.0
    SMVGLUCLAChep 	 1
    ASNATB0tc 	 1.0
    EX_prostge1_LPAREN_e_RPAREN_ 	 1
    r2331 	 3.0
    r1973 	 -1.0
    MG2er 	 1.0
    S6T24g 	 1.0
    r0497 	 3.0
    FAOXC121_5Em 	 3.0
    MCOATAm 	 -1.0
    RE0124C 	 1
    ABTti 	 1
    r0987 	 1
    r0696 	 -1.0
    r1995 	 -1.0
    METS 	 -1.0
    3HIBUPGLUC_Sthv 	 1
    4NPHSFte 	 1
    ACNAMlt 	 -1.0
    FBA2 	 3.0
    OCCOAtx 	 1
    KSII_CORE4t 	 1
    r1324 	 1
    EX_xoltri25_LPAREN_e_RPAREN_ 	 1
    EX_aldstrn_LPAREN_e_RPAREN_ 	 1
    7HPVSitr 	 1
    r0610 	 -1.0
    PI3P4K 	 0.0
    CYTK12 	 3.0
    r0081 	 -1.0
    THMTPt 	 -1.0
    3HPVSTETCOAitm 	 1
    EX_xylt_LPAREN_e_RPAREN_ 	 3
    RE1834C 	 0
    FAOXC182C182OHm 	 3.0
    r1749 	 -1.0
    r1656 	 3.0
    r1418 	 3.0
    DM_dctp_m_ 	 1
    DELACCRVSM23hc 	 1
    LVSTACOXD6Hhep 	 3.0
    PIK3 	 3.0
    r1615 	 -1.0
    EICOSTETt 	 1
    4HATVLACteb 	 -1.0
    r0756 	 2.0
    r2255 	 -1.0
    RE3597X 	 0
    r0786 	 3.0
    BGLUTCHLe 	 1
    r2393 	 -1.0
    ARGSL 	 -1.0
    FAOXTC102C101m 	 3.0
    UPP3S 	 -1.0
    r1421 	 1
    r0034 	 1
    r2503 	 -1.0
    3HSMVhep 	 3.0
    NACHEX27ly 	 1.0
    FUT16g 	 -1.0
    r1881 	 -1.0
    FAOXCPRIST1x 	 -1.0
    r2041 	 -1.0
    DOPAVESSEC 	 0.0
    RE1582L 	 1
    r1169 	 -1.0
    DNDPt34m 	 1
    EX_oh1_LPAREN_e_RPAREN_ 	 1
    RE3172R 	 1
    EX_HC02204_LPAREN_e_RPAREN_ 	 1
    RE0924C 	 -1.0
    DHEASABCCte 	 1.0
    EX_itp_LPAREN_e_RPAREN_ 	 1
    TSTSTERONEGLCte 	 1.0
    LNS14DM 	 -1.0
    SLCBK1 	 0.0
    B3GALT44g 	 -1.0
    34HPPOR 	 -1.0
    C120CPT2 	 -1.0
    RE3250X 	 3.0
    EX_cysacmp_LPAREN_e_RPAREN_ 	 1
    FUCGAL14ACGLCGALGLUSIDEtg 	 1
    r2330 	 3.0
    RE3233L 	 -1.0
    XOLDIOLONEtm 	 1
    PGLer 	 -1.0
    RE1860E 	 3.0
    FADtru 	 1
    H6ET3er 	 -1.0
    RE2636C 	 1
    3NTD7l 	 2.0
    FACOAL246_1 	 2.0
    r2198 	 -1.0
    GLACter 	 1
    CO2ter 	 1
    PPPItn 	 1
    FE3R2e 	 3.0
    EX_3octdec2crn_ 	 1
    RE2954C 	 3.0
    RE3201C 	 1
    r1605 	 3.0
    DLNLCGCRNt 	 -1.0
    P45021A2r 	 0
    RNDR4 	 2.0
    C4STMO2Pr 	 3.0
    RETNGLCtr 	 1
    RN0027C 	 1
    THRGLNNaEx 	 -1.0
    CHSTEROLt 	 3.0
    r1593 	 3.0
    DOLPMT1_Uer 	 1
    TETTET6CPT1 	 0
    r2099 	 1
    r0762 	 -1.0
    7BHGLZGLCtev 	 1
    PROSTGH2t 	 1
    RE3443X 	 3.0
    MMEm 	 -1.0
    AIRCr_PRASCS 	 3.0
    r1021 	 1
    RE3010C 	 1.0
    FAOXC122x 	 1.0
    1a_24_25VITD3Hm 	 1
    HIStiDF 	 -1.0
    6EPVStep 	 1
    HDECAACBP 	 3.0
    FAOXC225m 	 1.0
    r0886 	 1
    HEXDCRNe 	 1
    CMPSASn 	 -1.0
    S6T19g 	 1.0
    r0911 	 1
    RE3244R 	 1
    CYTDK2m 	 1
    HEXDIACTD 	 1
    RE1533M 	 3.0
    EX_3tdcrn_ 	 1
    r2136 	 -1.0
    r0627 	 2.0
    C101CRNe 	 1
    r0013 	 -1.0
    FA160ACPH 	 -1.0
    r2018 	 -1.0
    FOLOAT2tc 	 -1.0
    TRDR3 	 0.0
    RE3469C 	 1
    r1890 	 -1.0
    C142ACBP 	 3.0
    r0578 	 3.0
    SPH1Ptr 	 1
    RE2642C 	 1
    1HMDGLUCitr 	 1
    GLYSARCNc 	 3.0
    LEUGLYPEPT1tc 	 -1.0
    ACACt2m 	 0.0
    r1952 	 -1.0
    r0097 	 0.0
    DDECCRNe 	 1
    DM_mem2emgacpail_prot_hs_r_ 	 1
    RTOTALCRNCPT1 	 0
    EX_galside_hs_LPAREN_e_RPAREN_ 	 1
    NODe 	 1
    FAOXC205C184x 	 -1.0
    C180CRNt 	 -1.0
    RE2296C 	 1
    EX_CE1950_LPAREN_e_RPAREN_ 	 1
    GUMTCHOLe 	 1
    GLUt6 	 1.0
    r2210 	 -1.0
    GLCAASE7ly 	 3.0
    GALASE1ly 	 0
    r1915 	 -1.0
    RE1806C 	 1
    DCYTDn 	 -1.0
    RE1922C 	 1
    r0384 	 -1.0
    RE1521X 	 -1.0
    RE2051C 	 3.0
    PI345P3Pn 	 3.0
    RE2122C 	 1
    PCHOLPg_hs 	 3.0
    DM_Ser_Gly_Ala_X_Gly_ly_ 	 1
    r1563 	 3.0
    r2008 	 -1.0
    r0679 	 3.0
    RE2445C 	 0.0
    UAGDP 	 2.0
    r1584 	 3.0
    PVSGLUCitr 	 1
    IPDDI 	 3.0
    r1501 	 -1.0
    FAOXC225_4Z_7Z_10Z_13Z_16Zm 	 1.0
    r2155 	 -1.0
    EX_dcmp_LPAREN_e_RPAREN_ 	 3
    RE3232C 	 1
    EX_2425dhvitd2_LPAREN_e_RPAREN_ 	 1
    GAL3ST11 	 -1.0
    EX_lac_D_LPAREN_e_RPAREN_ 	 1
    r1488 	 1
    MI14P4P 	 1
    NABTNOm 	 3.0
    PHYQt 	 1
    DNDPt55m 	 -1.0
    RE3400M 	 3.0
    CHSTEROLt1 	 -1.0
    TETPENT6COAtx 	 1
    HIBDm 	 -1.0
    C6CRNe 	 1
    ARGLYSex 	 -1.0
    MM7Cbg 	 3.0
    GAO2 	 1.0
    UMPK6n 	 3.0
    RE3170R 	 3.0
    RBTt 	 1
    NRVNCCOAtx 	 1
    RE2541C 	 -1.0
    RN0001R 	 2.0
    DOCOSACT 	 0.0
    EX_fucfucfucgalacglc13galacglcgal14acglcgalgluside_hs_LPAREN_e_RPAREN_ 	 1
    AGTix 	 -1.0
    ACNAMtn 	 1
    ABTD 	 1
    CBR2 	 1.0
    r1587 	 3.0
    EX_n2m2nmasn_LPAREN_e_RPAREN_ 	 3
    EX_Lcystin_LPAREN_e_RPAREN_ 	 1
    CBL2tm 	 1.0
    ESTRADIOLt 	 1
    DNDPt31m 	 -1.0
    r1822 	 -1.0
    r2520 	 -1.0
    RE1923C 	 1
    RE2705C 	 1
    FAOXC10DCC8DCx 	 -1.0
    r1746 	 -1.0
    ALLOPtepvb 	 -1.0
    RE2150C 	 -1.0
    RE1534M 	 3.0
    NTD4 	 3.0
    S2TASE2ly 	 0.0
    ECOAH12m 	 3.0
    r2436 	 -1.0
    RE3338M 	 3.0
    r1025 	 1
    RE0944C 	 1
    ALAALAPEPT1tc 	 -1.0
    r0726 	 3.0
    EX_ocdcea_LPAREN_e_RPAREN_ 	 3
    RN0031R 	 0.0
    IDOAASE4ly 	 -1.0
    TSACMSULhc 	 1
    S6TASE11ly 	 3.0
    ARTFR31 	 1
    r0560 	 -1.0
    PTE5x 	 0
    CHOLK 	 -1.0
    EX_sbt_DASH_d_LPAREN_e_RPAREN_ 	 3
    C141ACBP 	 3.0
    r1427 	 -1.0
    BIOCYTtn 	 1
    LSTO1r 	 2.0
    HMGCOAtx 	 1
    RE0344M 	 1.0
    FAOXC200180x 	 0.0
    r1332 	 1
    SUBERCROT 	 3.0
    r2269 	 -1.0
    r0446 	 1
    CAROtr 	 1
    EX_CE5797_LPAREN_e_RPAREN_ 	 1
    AM1CCStev 	 1
    r0022 	 3.0
    RE3164C 	 1
    r1364 	 1
    6BHGLZGLCtev 	 1
    GLU5Km 	 1.0
    RE3308R 	 1
    EX_ebastine_LPAREN_e_RPAREN_ 	 1
    RE3582X 	 1.0
    ADPT 	 3.0
    r2086 	 0.0
    CHLPCTD 	 0.0
    41R1H2MAE12BOOX 	 3.0
    CK 	 0
    P4507B12r 	 -1.0
    RE1938R 	 0.0
    RE1629C 	 1
    ENMAN6g 	 1
    r2238 	 -1.0
    UNK3 	 1
    r0975 	 1
    RE2912M 	 2.0
    FAOXC241C221x 	 0.0
    SARDHm 	 -1.0
    RE3011M 	 -1.0
    MM8Ag 	 3.0
    RE3498R 	 1
    PTDCACRNCPT2 	 -1.0
    DORNOp 	 -1.0
    AMY1e 	 0
    GLNS 	 3.0
    SUCOAS1m 	 1.0
    CHOLESTTDe 	 1
    r2050 	 -1.0
    EX_leuktrC4_LPAREN_e_RPAREN_ 	 1
    ENMAN4g 	 1
    IVCRNe 	 1
    FAOXC184C163x 	 -1.0
    CRVSthc 	 -1.0
    RE3455C 	 1
    CARIBUPGLU_Sitr 	 1
    RE3174R 	 3.0
    RE3381E 	 -1.0
    DUTPDPn 	 3.0
    THMTP 	 -1.0
    DMHPTCRNt 	 1
    ASNALANaEx 	 -1.0
    FPGS3 	 -1.0
    r2502 	 -1.0
    EX_mthgxl_LPAREN_e_RPAREN_ 	 1
    RE3301C 	 -1.0
    TS3 	 1
    PPBNGS 	 2.0
    CHOLtr 	 1
    r0682 	 -1.0
    RE3044C 	 1
    EX_4ohmdz_LPAREN_e_RPAREN_ 	 1
    GALNACT1g 	 0.0
    GLYCTO1p 	 -1.0
    RE3176C 	 1
    AKR1C1 	 3.0
    ACMPGLUTdt 	 1
    RE3503C 	 1
    r1761 	 -1.0
    r0649 	 3.0
    RE2985M 	 0.0
    O2tp 	 1
    S6TASE7ly 	 -1.0
    PYRt2m 	 0.0
    RE3125C 	 1
    URATEt 	 1
    7AHGLZitr 	 1
    1HIBUPGLUitr 	 1
    ADNK1 	 2.0
    CYTK8 	 3.0
    r1902 	 -1.0
    THFtl 	 1
    RE2524C 	 1
    EX_HC02154_LPAREN_e_RPAREN_ 	 1
    FAOXC2452256x 	 0.0
    BETBGTtc 	 -1.0
    GLYLEUHYDROc 	 3.0
    2HATVLAChc 	 1
    EX_na1_LPAREN_e_RPAREN_ 	 3
    PTPATe 	 1
    AM1C4N9CShc 	 1
    MCCCrm 	 1.0
    RE2031M 	 1
    sink_dd2coa_LPAREN_c_RPAREN_ 	 1
    r2133 	 0.0
    RE2032M 	 1
    r2163 	 -1.0
    RE3166C 	 1
    DPCOAK 	 -1.0
    EX_56eppvs_LPAREN_e_RPAREN_ 	 1
    RN0020R 	 3.0
    EX_ndersv_LPAREN_e_RPAREN_ 	 1
    GHMT2r 	 -1.0
    RE3154C 	 1
    EX_dgtp_LPAREN_e_RPAREN_ 	 1
    RE3157X 	 3.0
    ASNSERNaEx 	 -1.0
    ARGN 	 3.0
    OPAHir 	 -1.0
    ANDRSTRNtr 	 1
    EX_6melvst_LPAREN_e_RPAREN_ 	 1
    RE3119R 	 1
    r1503 	 -1.0
    THYOXt 	 1.0
    ORNTArm 	 3.0
    GUR1PP 	 1
    RAtn 	 1
    DOPABMO 	 0.0
    RE3570C 	 1
    r1607 	 3.0
    ADAe 	 -1.0
    GLNASNNaEx 	 -1.0
    MEVK1c 	 3.0
    AM9CSAhr 	 3.0
    MACACI 	 -1.0
    RE3250C 	 1
    CSNt 	 -1.0
    r1549 	 3.0
    RE1526M 	 3.0
    G3PD1 	 -1.0
    EX_idp_LPAREN_e_RPAREN_ 	 1
    EX_6ahglz_LPAREN_e_RPAREN_ 	 1
    RE2866C 	 1
    GSNt2r 	 0.0
    EX_HC02217_LPAREN_e_RPAREN_ 	 1
    DM_Ser_Thr_ly_ 	 1
    EX_aprgstrn_LPAREN_e_RPAREN_ 	 1
    P4503A7r 	 0
    MAL_Lte 	 1
    RE1845C 	 -1.0
    MELATN23DOX 	 -1.0
    CARhPTtc 	 -1.0
    DOLDPP_Uer 	 1
    r1313 	 0.0
    HSD3B2r 	 -1.0
    RE2030M 	 1
    H2O2tp 	 1
    r0813 	 -1.0
    DOPAENT4tc 	 -1.0
    EX_ivcrn_ 	 1
    RE3390M 	 3.0
    r2089 	 0.0
    RE2514C 	 3.0
    r1769 	 -1.0
    DTTPtn 	 1
    HDD2COAtx 	 1
    ACOX22x 	 -1.0
    r1560 	 3.0
    DCYTt 	 -1.0
    RE3514R 	 -1.0
    RE3564X 	 3.0
    C161CPT12 	 0
    EX_ptrc_LPAREN_e_RPAREN_ 	 3
    EX_HC02160_LPAREN_e_RPAREN_ 	 1
    r0626 	 1
    r1608 	 3.0
    RE3414C 	 1
    TMDM3itr 	 1
    RE3431C 	 1
    RE1527M 	 3.0
    ACALDtx 	 1
    2AMACSULT 	 0
    r1180 	 1.0
    DESFVSteb 	 1
    CYSGLTH 	 1
    LNLNCACPT1 	 0
    DLNLCGt 	 -1.0
    MI1346PKn 	 -1.0
    r0596 	 3.0
    CRNCARtp 	 0.0
    r1977 	 -1.0
    A4GNT1g 	 -1.0
    GALASE7ly 	 0
    r1933 	 -1.0
    ACGALK 	 1
    RE1062C 	 1
    ASNt4 	 3.0
    RE2398R 	 -1.0
    FPGS6 	 -1.0
    TYRB0AT3tc 	 -1.0
    AHEXASE2ly 	 1.0
    r0598 	 1
    GMAND 	 -1.0
    EX_oxy7rb_LPAREN_e_RPAREN_ 	 1
    r0908 	 1
    EX_acgalfucgalacgalfuc12gal14acglcgalgluside_hs_LPAREN_e_RPAREN_ 	 1
    r1734 	 -1.0
    G2M8MASNterg 	 1
    CO2tp 	 1
    ME2 	 3.0
    r0670 	 -1.0
    HDCACBP 	 3.0
    RE3495C 	 1
    r1842 	 -1.0
    DM_4abut_LPAREN_n_RPAREN_ 	 1
    RE1527X 	 3.0
    MI134PP 	 3.0
    EX_5ohfvsglu_LPAREN_e_RPAREN_ 	 1
    AQCOBALt 	 1
    RADH4 	 1
    EX_arachd_LPAREN_e_RPAREN_ 	 3
    RE3413C 	 1
    r2166 	 -1.0
    EX_bglc_LPAREN_e_RPAREN_ 	 1
    r1948 	 -1.0
    PTVSTATPtu 	 0.0
    SLDx 	 3.0
    CLFORtex 	 -1.0
    S4T3g 	 -1.0
    NACHEX14ly 	 1.0
    r2184 	 -1.0
    EX_csn_LPAREN_e_RPAREN_ 	 3
    MCPST 	 -1.0
    RE3258R 	 0.0
    RE3511R 	 -1.0
    RE3136C 	 2.0
    RE3095X 	 1
    2HBO 	 3.0
    r2287 	 -1.0
    EX_HC02172_LPAREN_e_RPAREN_ 	 1
    FUT17g 	 -1.0
    EX_3hibupglu_S_LPAREN_e_RPAREN_ 	 1
    DHFtl 	 1
    r1181 	 1.0
    EX_7bhglzglc_LPAREN_e_RPAREN_ 	 1
    FVStu 	 -1.0
    RE2624X 	 -1.0
    r2325 	 3.0
    XOLTRIOLtm 	 1
    r1499 	 -1.0
    FACOAL120i 	 1
    PROSTGF2t 	 -1.0
    RE1096C 	 1
    6HLVSTitr 	 1
    H8MTer_U 	 -1.0
    DNDPt7m 	 -1.0
    AM1ACSitr 	 1
    35DHPVSthc 	 -1.0
    FA141ACPH 	 -1.0
    r1649 	 3.0
    FAOXC5OHc 	 0
    r1962 	 -1.0
    FUT32g 	 1.0
    GQ1Bte 	 1
    r1717 	 -1.0
    r1872 	 -1.0
    PRGSTRNt 	 1
    HC00342te 	 1
    S2L2FN2M2MASNt 	 1
    PLA2_2 	 2.0
    r0437 	 1
    RE0828C 	 0.0
    r1502 	 -1.0
    RE2318R 	 3.0
    GAO1g 	 1.0
    r0641 	 3.0
    RE2533C 	 1
    RE3520C 	 1
    FAOXC184C164m 	 3.0
    RE3076X 	 1
    RE3084X 	 3.0
    RE3440R 	 -1.0
    r1544 	 3.0
    ADPCOAPTE 	 -1.0
    RE2861C 	 1
    RE0936E 	 0.0
    P45019A1r 	 -1.0
    r1955 	 -1.0
    ST6GALNAC24 	 3.0
    r0318 	 0.0
    EX_estradiol_LPAREN_e_RPAREN_ 	 1
    r2343 	 0
    ACMPGLUitr 	 1
    C162OHc 	 0
    34DHOXPEGt 	 1
    RE1835C 	 0
    CRNrtx 	 1
    EX_lgnc_LPAREN_e_RPAREN_ 	 1
    MG1er 	 -1.0
    UAG2EMAi 	 -1.0
    Rtotaltp 	 1
    M4BTAer 	 -1.0
    UGT1A9r 	 -1.0
    r1889 	 -1.0
    MANt1r 	 3.0
    C160CPT2 	 -1.0
    NACHEXA7ly 	 -1.0
    r1322 	 1
    RE1811R 	 -1.0
    MMSAD3m 	 -1.0
    AM1CGLChr 	 1
    RDH2 	 -1.0
    NMPTRCOX 	 0.0
    FMNAT 	 -1.0
    35CGMPtn 	 1
    H6_ET2er 	 1.0
    ACNAM9PL2 	 2.0
    DMHPTCRNCPT1 	 1
    SIAASE3ly 	 -1.0
    EX_mal_L_LPAREN_e_RPAREN_ 	 3
    DNDPt12m 	 -1.0
    FAOXC4C4DCc 	 0
    3ISPVSteb 	 0.0
    EX_7hpvs_LPAREN_e_RPAREN_ 	 1
    ASNS1 	 1.0
    r2209 	 -1.0
    r1061 	 1
    r2411 	 -1.0
    ACACT4p 	 3.0
    PS_HStg 	 1
    PVSitr 	 1
    RE2048N 	 2.0
    PHETHPTOX2 	 0
    PROAKGOX1r 	 2.0
    r0871 	 1
    HIStN1 	 -1.0
    FAOXC11BRC9BRx 	 -1.0
    NACASPAH 	 0.0
    EX_1ohmdz_LPAREN_e_RPAREN_ 	 1
    GALASE6ly 	 0
    RE1266C 	 1
    C18OHc 	 0
    r1770 	 -1.0
    FUT93g 	 -1.0
    PIK5n 	 1
    PI45P5Pn 	 1
    EX_ser_D_LPAREN_e_RPAREN_ 	 1
    14HMDZALThr 	 3.0
    RE3532C 	 1
    DHEASt 	 -1.0
    25VITD3Hm 	 0.0
    C160CRNt 	 -1.0
    CYSTAm 	 3.0
    BAAT5x 	 1
    RADH 	 1
    r1459 	 1
    RE3583X 	 -1.0
    M14NTg 	 -1.0
    NFDLACtep 	 1
    r1975 	 -1.0
    ACN23ACNGALGBSIDEte 	 1
    LNELDCCPT1 	 0
    NACHEXA17ly 	 -1.0
    FAOXC122_3Z_6Zx 	 3.0
    C141OHe 	 1
    UGT1A1r 	 0
    TYRATB0tc 	 1.0
    CSPG_Etly 	 1
    PI5P4Kn 	 1
    2HATVLACthc 	 -1.0
    ACNMLr 	 1
    IDHPOXOX2b 	 -1.0
    RE2871C 	 1
    r2236 	 -1.0
    RE0935E 	 0
    r2057 	 -1.0
    RE3456C 	 1
    ACS 	 1.0
    r2309 	 -1.0
    3HXKYNDCL 	 -1.0
    STCOATxc 	 -1.0
    RE2974G 	 -1.0
    NTD5l 	 2.0
    r1966 	 -1.0
    r2111 	 0.0
    12HTACRtu 	 1
    NTD9l 	 2.0
    RE1700C 	 1
    GLC3MEACPhr 	 0
    3HPVSthc 	 -1.0
    SOAT12 	 -1.0
    LEUTA 	 -1.0
    KHte 	 1
    r1143 	 1
    RE2864C 	 1
    r0122 	 -1.0
    S4TASE3ly 	 -1.0
    10FTHF5GLUtm 	 1
    UMPKn 	 3.0
    NTD2l 	 2.0
    EX_chtn_LPAREN_e_RPAREN_ 	 1
    r1515 	 3.0
    RE2404R 	 0
    r2296 	 -1.0
    r1078 	 1
    P4502C92 	 -1.0
    HEXCt 	 1
    CYSO 	 -1.0
    S6TASE17ly 	 3.0
    LST4EXPthc 	 1
    EX_eaflatoxin_LPAREN_e_RPAREN_ 	 1
    CBPSam 	 -1.0
    r1068 	 1
    PI45P3Kn 	 1
    EX_leugly_LPAREN_e_RPAREN_ 	 1
    r1852 	 -1.0
    LEUKTRE4t 	 0
    RE3195M 	 -1.0
    RE2851C 	 1
    r1677 	 -1.0
    RDH1 	 -1.0
    RE1099C 	 -1.0
    r0385 	 -1.0
    r1432 	 -1.0
    PTVSTLACitr 	 1
    CHOLATEt 	 -1.0
    r2259 	 -1.0
    r1786 	 -1.0
    RE2625M 	 -1.0
    ACNACNGAL14ACGLCGALGLUSIDEtg 	 1
    r1514 	 3.0
    CYOR_u10m 	 0
    CITMCOAHm 	 1
    RE3012C 	 3.0
    NACHEX20ly 	 1.0
    RE3072X 	 -1.0
    RETHe 	 1
    EX_pi_LPAREN_e_RPAREN_ 	 3
    3HPVStep 	 1
    FAOXC164C165m 	 3.0
    RE1516M 	 -1.0
    C226CRNt 	 1
    S4T5g 	 1
    RE3087X 	 -1.0
    EX_aqcobal_LPAREN_e_RPAREN_ 	 1
    5DHFtl 	 1
    r1914 	 -1.0
    BIDGLCURr 	 0
    TAUPAT1c 	 -1.0
    r2484 	 -1.0
    RE3163C 	 1
    DHCHOLESTANATEtm 	 1
    MI14PP 	 3.0
    r0909 	 1
    AG13T6g 	 1.0
    DM_dgpi_prot_hs_r_ 	 1
    PA_HStg 	 1
    RSVLAChv 	 1
    EX_coumarin_LPAREN_e_RPAREN_ 	 1
    C12DCACT 	 1
    EX_val_L_LPAREN_e_RPAREN_ 	 3
    RE1819M 	 0
    EX_3hpvs_LPAREN_e_RPAREN_ 	 3
    CARVEOLte 	 1
    r2501 	 -1.0
    ISOLVSTAChep 	 1
    r0782 	 -1.0
    PI3P3Pn 	 1
    P45017A2r 	 -1.0
    RE3421R 	 1
    PE_HStg 	 1
    RE1809R 	 -1.0
    LSTNM4tev 	 1
    STS1r 	 -1.0
    MTHGXLt 	 1
    r0698 	 3.0
    RE3506C 	 1
    C142CPT1 	 0
    r1654 	 3.0
    r1827 	 -1.0
    r2183 	 -1.0
    FACOAL2251 	 2.0
    O2Stx 	 1
    EX_gmp_LPAREN_e_RPAREN_ 	 1
    CYTDtn 	 1
    CSNAT2x 	 3.0
    CYTK14n 	 3.0
    PROFVShc 	 1
    PAN4PP 	 1
    BTNDm 	 -1.0
    r1604 	 3.0
    r2051 	 -1.0
    RE1517M 	 -1.0
    EX_ptvst_LPAREN_e_RPAREN_ 	 1
    FAS120COA 	 -1.0
    RNDR3 	 2.0
    r1755 	 -1.0
    r2256 	 -1.0
    r1864 	 -1.0
    EX_thym_LPAREN_e_RPAREN_ 	 1
    EX_deoxfvs_LPAREN_e_RPAREN_ 	 1
    MGACONm 	 1
    r1645 	 3.0
    RE3447X 	 -1.0
    TCHOLAt 	 1.0
    DM_pe_hs_LPAREN_r_RPAREN_ 	 3
    4HPROLTASCT1 	 3.0
    r1173 	 1.0
    ACGALFUCGALACGALFUCGALACGLCGAL14ACGLCGALGLUSIDEtg 	 1
    r1434 	 1
    RE1952R 	 0.0
    RE3288C 	 1
    GULLACter 	 1
    RE3422C 	 1
    FAOXC183_6Z_9Z_12Zm 	 3.0
    NTD6l 	 2.0
    PDE4 	 -1.0
    ADK1 	 3.0
    GLYt4 	 3.0
    10FTHF6GLUtm 	 1
    r1671 	 -1.0
    RE1473C 	 1
    RE2973G 	 3.0
    PUNP7 	 2.0
    PRPNCOAHYDm 	 3.0
    FAOXC182806m 	 2.0
    EX_c8crn_ 	 1
    DTDPtn 	 1
    RE1134R 	 -1.0
    r2234 	 -1.0
    EX_atvacid_LPAREN_e_RPAREN_ 	 1
    FAOXC103C102m 	 3.0
    EX_HC02187_LPAREN_e_RPAREN_ 	 1
    EX_xoltri24_LPAREN_e_RPAREN_ 	 1
    RE3502C 	 1.0
    LNLNCACRNt 	 -1.0
    3MOX4HOXPGALDOX 	 2.0
    HSPGt 	 1
    P45027A1m 	 -1.0
    IDHPOXOXb 	 -1.0
    DHORTS 	 0.0
    DCSPTN1CPT2 	 -1.0
    EX_23cump_LPAREN_e_RPAREN_ 	 1
    RE2985X 	 3.0
    CSPG_At 	 1
    NACHEX18ly 	 1.0
    GALGT3 	 -1.0
    r2233 	 -1.0
    PYDXK 	 3.0
    EHGLATm 	 3.0
    PVSHtu 	 1
    5THFtm 	 1
    FACOAL204 	 2.0
    r2525 	 -1.0
    r0426 	 3.0
    B3GNT313g 	 -1.0
    RE1531M 	 3.0
    2HATVACIDOXDhc 	 3.0
    NRPPHRVESSEC 	 0.0
    UMPK4 	 3.0
    FUCFUC12GAL14ACGLCGALGLUSIDEtg 	 1
    ACNGALACGLCGAL14ACGLCGALGLUSIDEtg 	 1
    5HLTDL 	 -1.0
    FAOXC5C3x 	 -1.0
    CSBPASEly 	 1
    3AIB_Dtm 	 1
    SEBACIDTD 	 1
    SMVACIDATPteb 	 -1.0
    EX_strch2_LPAREN_e_RPAREN_ 	 3
    FAOXC16OHC16r 	 -1.0
    r1446 	 3.0
    URIt5le 	 0.0
    DM_oretn_n_ 	 1
    RE3010R 	 1
    3HIBUPGLUC_Sitr 	 1
    6MELVSTitr 	 1
    PGS 	 2.0
    RE3228C 	 1
    r2213 	 -1.0
    r0466 	 1
    DEOXFVShc 	 -1.0
    r1617 	 -1.0
    ARGDCm 	 -1.0
    RE3152R 	 1
    r1002 	 1
    GALASE17ly 	 0
    FAH2 	 -1.0
    r1538 	 -1.0
    EX_HC00822_LPAREN_e_RPAREN_ 	 1
    EX_hexc_LPAREN_e_RPAREN_ 	 1
    ORNALArBATtc 	 -1.0
    DHPR 	 0.0
    EX_udpg_LPAREN_e_RPAREN_ 	 1
    RE3560C 	 0
    RE2373C 	 1
    TXASr 	 -1.0
    r1917 	 -1.0
    RE1699C 	 1
    r2001 	 -1.0
    FAOXC15C13m 	 3.0
    r2273 	 -1.0
    ALAATB0tc 	 1.0
    RE1632C 	 1
    EX_lstn_LPAREN_e_RPAREN_ 	 1
    FAOXC13C11m 	 3.0
    FVSitx 	 1
    EX_leu_L_LPAREN_e_RPAREN_ 	 3
    FACOAL245_2 	 2.0
    RE3019R 	 -1.0
    r0642 	 3.0
    DM_1a25dhvitd3_LPAREN_n_RPAREN_ 	 1
    CRGLZABCt 	 1
    RE3041N 	 2.0
    GLCtly 	 1
    PAIL_HStn 	 1
    r0644 	 -1.0
    r1641 	 3.0
    DCSPTN1COAtx 	 1
    AM1A4NCSteb 	 1
    6CSMVhep 	 3.0
    N2M2NMASNtly 	 1
    15DMTitr 	 1
    RE2909M 	 3.0
    RE0691C 	 1
    r1999 	 -1.0
    r1893 	 -1.0
    AG13T11g 	 1.0
    r0361 	 3.0
    STS2 	 -1.0
    ALDD20xm 	 3.0
    DHFtm 	 1
    RE3410C 	 1
    FAOXC6DCC4DCx 	 -1.0
    RE2319X 	 1.0
    r1832 	 -1.0
    ABUTt4_2_r 	 -1.0
    FAOXC142_5Z_8Zx 	 0.0
    3HXKYNOXDA 	 3.0
    r0983 	 -1.0
    LRAT1 	 1
    CHOLt4 	 -1.0
    RE3194M 	 3.0
    DNDPt45m 	 -1.0
    EX_c51crn_ 	 1
    BTNDe 	 -1.0
    r0997 	 3.0
    BCRNe 	 1
    RE2995M 	 2.0
    CHSTEROLtg 	 3.0
    r1657 	 3.0
    CYSACMPAChc 	 1
    EX_sucr_LPAREN_e_RPAREN_ 	 3
    HSD3B3r 	 -1.0
    r0739 	 3.0
    PA_HSter 	 1
    RE3446C 	 -1.0
    r2146 	 -1.0
    H2O2tm 	 1
    RE3075C 	 1
    PTVSTM3eb 	 0.0
    ADNt5le 	 0.0
    DESAT18_8 	 -1.0
    FAOXC10C10OHm 	 1
    5OHFVSteb 	 1
    OAGD3tg 	 1
    EX_no2_LPAREN_e_RPAREN_ 	 3
    RE1523M 	 3.0
    LPCHOLt 	 1
    RE3449C 	 1
    CLPNDCPT1 	 0
    r2516 	 0.0
    r2268 	 -1.0
    RE0512M 	 3.0
    GSNt 	 -1.0
    HESTRATRIOLtr 	 1
    PCHOLPm_hs 	 3.0
    r0822 	 1
    ACGALFUCGALACGALFUCGALACGLCGAL14ACGLCGALGLUSIDEte 	 1
    FAOXTC122TC101m 	 0.0
    G12MT1_U 	 1
    CYSASNNaEx 	 -1.0
    NTD2m 	 -1.0
    FATP6t 	 -1.0
    C4DCCACT 	 -1.0
    RE2526C 	 1
    CRNt 	 -1.0
    AG13T7g 	 1.0
    PI45P5P 	 0.0
    ALR2 	 3.0
    ORPT 	 1.0
    UDPGLCter 	 1
    r1808 	 -1.0
    LVSTOXD6Hhep 	 3.0
    ADA 	 -1.0
    DOPASULT 	 1.0
    r1940 	 -1.0
    CYTK6n 	 3.0
    r0842 	 1
    r2515 	 1
    CYSB0AT3tc 	 -1.0
    TETHEX3t 	 1
    EX_4mop_LPAREN_e_RPAREN_ 	 3
    RE3341M 	 1
    UGLT 	 -1.0
    ADSL1 	 3.0
    sink_decdicoa_LPAREN_c_RPAREN_ 	 1
    GLACO 	 3.0
    LVSTPGPtu 	 -1.0
    IDOAASE3ly 	 -1.0
    CVM23GLUChc 	 0
    GLNLASEer 	 1
    EX_gua_LPAREN_e_RPAREN_ 	 3
    RE0922C 	 -1.0
    r0514 	 1
    MMTSADm 	 1
    r1481 	 3.0
    r2142 	 -1.0
    r1076 	 1
    r1754 	 -1.0
    MCDp 	 -1.0
    r1903 	 -1.0
    DM_core8_g_ 	 1
    RE3244C 	 1
    EX_pro_D_LPAREN_e_RPAREN_ 	 1
    RE1905C 	 1
    RE1587L 	 1
    RE1099L 	 -1.0
    RE2857C 	 1
    r1573 	 3.0
    L_LACt4r 	 -1.0
    GLNATB0tc 	 1.0
    r2280 	 -1.0
    RE3044N 	 2.0
    DM_dctp_n_ 	 3
    r1883 	 -1.0
    RE2383R 	 -1.0
    COAtn 	 1
    CSm 	 3.0
    RE1818M 	 1.0
    RE3343M 	 1.0
    r1088 	 1
    MERACMPthc 	 -1.0
    EX_acmana_LPAREN_e_RPAREN_ 	 3
    UREAt 	 0.0
    PI5P4K 	 0.0
    AKGtp 	 1
    THRGLNexR 	 -1.0
    RE3339M 	 1.0
    r0355 	 3.0
    ACMPGLUTTRsc 	 1
    RE1050N 	 -1.0
    r1116 	 1
    COAtl 	 1
    MI145PP 	 1.0
    RE2632M 	 -1.0
    RE2660C 	 1
    r1600 	 3.0
    EX_6htststerone_LPAREN_e_RPAREN_ 	 1
    r0728 	 -1.0
    r1523 	 3.0
    SPODMm 	 3.0
    IDOURtly 	 -1.0
    RETFAt2 	 1
    RE3013C 	 1
    RE3575X 	 3.0
    TRPt 	 -1.0
    ALAGLNNaEx 	 -1.0
    EX_tmdm1_LPAREN_e_RPAREN_ 	 1
    RE3123R 	 1
    1a25DHVITD3TRn 	 1
    MELATNOX 	 -1.0
    RE1828C 	 1
    CBL2OR 	 1
    THYMDt1 	 -1.0
    CMPSAS 	 -1.0
    r0399 	 0
    XOLDIOLONEt 	 1
    APRGSTRNte 	 1
    r1820 	 -1.0
    DNDPt26m 	 1
    CRTSLtm 	 1
    PPA 	 3.0
    EX_hista_LPAREN_e_RPAREN_ 	 1
    PROSTGE2t 	 3.0
    r0954 	 0.0
    12HTACRitr 	 1
    FAOXC101C8m 	 0.0
    RE2155C 	 1
    RE3336M 	 3.0
    FADtm 	 1
    RE2958C 	 1
    RE3150C 	 1
    r2485 	 -1.0
    EX_9_cis_retfa_LPAREN_e_RPAREN_ 	 1
    IBUPGLUCtchep 	 1
    56EPPVSteb 	 0.0
    EX_34dhphe_LPAREN_e_RPAREN_ 	 1
    ALAyLATthc 	 -1.0
    EX_o2s_LPAREN_e_RPAREN_ 	 3
    EX_gthrd_LPAREN_e_RPAREN_ 	 3
    r2193 	 -1.0
    r1909 	 -1.0
    RE3630C 	 1
    CERT2gt 	 0.0
    RE2252C 	 1
    RE3230C 	 1
    r2039 	 -1.0
    NTD7l 	 2.0
    r1904 	 -1.0
    P45027A16m 	 -1.0
    CITMCOALm 	 1
    FAOXC204_5Z_8Z_11Z_14Zx 	 0.0
    LEUtec 	 -1.0
    THMt3 	 3.0
    PSERT 	 2.0
    EX_5ohfvs_LPAREN_e_RPAREN_ 	 1
    CKc 	 3.0
    PTVSTtu 	 -1.0
    r0386 	 -1.0
    r1920 	 -1.0
    H3ETer 	 -1.0
    C2M26DCOAHLx 	 2.0
    GLYCtm 	 1
    r0523 	 0
    r1575 	 3.0
    ACN23ACNGALGBSIDEtg 	 1
    r2397 	 -1.0
    TSACGLUCtev 	 1
    AM4NCSteb 	 1
    NRVNCCOAtxc 	 -1.0
    P45011B11m 	 -1.0
    RE2677G 	 1
    B3GNT39g 	 -1.0
    FACOAL2042 	 2.0
    CVM1GLUChc 	 0
    THRD_L 	 -1.0
    RE3624X 	 3.0
    CRVSM24tev 	 1
    ABTArm 	 -1.0
    M16N4Tg 	 1
    Rtotaltl 	 1
    P45019A2r 	 -1.0
    RE3191M 	 3.0
    ASNtN1 	 -1.0
    r1946 	 -1.0
    r2094 	 0.0
    LGNCFATPtc 	 -1.0
    PEAMNO 	 3.0
    UGT1A10r 	 0
    H4ETer 	 -1.0
    B3GNT32g 	 -1.0
    P45027A14m 	 -1.0
    3HAO 	 -1.0
    DADAe 	 -1.0
    PROSTGD2t 	 3.0
    r0731 	 -1.0
    ACSMCT1 	 -1.0
    r0331 	 1
    AMPDA 	 2.0
    DOLPMT4_Uer 	 1
    EX_nadp_LPAREN_e_RPAREN_ 	 1
    S6T12g 	 -1.0
    FUT911g 	 -1.0
    B3GNT12g 	 -1.0
    Htm 	 2.0
    r1176 	 1.0
    RE0581C 	 1
    RE1627C 	 1
    FUT94g 	 -1.0
    DM_dem2emgacpail_prot_hs_r_ 	 1
    FOLABCCte 	 1.0
    RE1818R 	 3.0
    M7MASNBterg 	 1
    RE2235R 	 -1.0
    XOL27OHtm 	 1
    r1147 	 -1.0
    r2382 	 -1.0
    biomass_RNA 	 3
    GLCAT8g 	 -1.0
    HEXCCOAtx 	 1
    PDXPP 	 0
    EPOXTAChr 	 1
    r2150 	 -1.0
    RE3147C 	 1
    S6T8g 	 -1.0
    EX_dhf_LPAREN_e_RPAREN_ 	 1
    r1748 	 -1.0
    r1588 	 3.0
    r1961 	 -1.0
    GLUVESSEC 	 -1.0
    DIGALSIDEtg 	 1
    NADH2_u10m 	 -1.0
    RE3022C 	 1
    BTNt4i 	 -1.0
    RE3020C 	 1
    G1M6MASNB1terg 	 1
    RE3011R 	 1
    r2534 	 -1.0
    RE2910X 	 3.0
    RE3633C 	 1
    RE2633C 	 1
    RE3156X 	 3.0
    C6DCc 	 0
    r1528 	 1
    ALADGLYexR 	 -1.0
    6DHFtl 	 1
    r1762 	 -1.0
    RE0920C 	 -1.0
    GLNB0AT3tc 	 -1.0
    XOLTRI24te 	 1
    r2013 	 -1.0
    RE1834X 	 -1.0
    r2284 	 -1.0
    RE2410N 	 0.0
    5ADTSTSTERONEGLCtr 	 1
    FAOXC260240x 	 0.0
    EX_taur_LPAREN_e_RPAREN_ 	 1
    PItg 	 1
    CLPNDCOAtx 	 1
    RE2638C 	 -1.0
    PDE4g 	 -1.0
    RE3448C 	 3.0
    6HSMVACIDhep 	 -1.0
    FA140ACPH 	 -1.0
    ILEB0AT3tc 	 -1.0
    RE1836M 	 3.0
    34DHXMANDACOX_NADP_ 	 2.0
    DHCRD1 	 3.0
    LYStiDF 	 3.0
    r0390 	 -1.0
    DOLPGT3_Ler 	 -1.0
    GARFT 	 -1.0
    r1531 	 0.0
    EX_3ispvs_LPAREN_e_RPAREN_ 	 1
    ADNK4 	 1
    r0776 	 2.0
    HSD17B9r 	 2.0
    r2059 	 -1.0
    CSAtu 	 -1.0
    EX_25hvitd3_LPAREN_e_RPAREN_ 	 1
    r1519 	 3.0
    ALOX15 	 3.0
    PRPNCOAHYDx 	 -1.0
    ALATHRNaEx 	 3.0
    FAOXC183_9Z_12Z_15Zm 	 3.0
    RE2908M 	 3.0
    EX_nac_LPAREN_e_RPAREN_ 	 3
    PROD2 	 -1.0
    LTA4H 	 3.0
    RE3345M 	 2.0
    r1651 	 3.0
    RE0916G 	 -1.0
    1331TACRtev 	 1
    r2310 	 -1.0
    r2102 	 1
    FUCASE2ly 	 -1.0
    4HOXPACDOX_NADP_ 	 2.0
    5AOPtm 	 1
    r2482 	 -1.0
    TMDK1 	 2.0
    ALACYSNaEx 	 3.0
    BAMPPALDOX 	 3.0
    r1879 	 -1.0
    r2182 	 -1.0
    VLCSp 	 -1.0
    EX_desfvs_LPAREN_e_RPAREN_ 	 1
    EX_tststerone_LPAREN_e_RPAREN_ 	 1
    r2317 	 3.0
    r2079 	 0.0
    CSASULPtev 	 1
    TRIODTHYt 	 1.0
    RE0927R 	 0
    r2253 	 -1.0
    r0672 	 -1.0
    RE2474R 	 3.0
    FOLOAT1tc 	 -1.0
    EX_11_cis_retfa_LPAREN_e_RPAREN_ 	 1
    r2085 	 0.0
    LPS 	 -1.0
    15DMTtu 	 1
    r2508 	 1
    S23T3g 	 3.0
    DOLK_U 	 1
    H6_ETer 	 -1.0
    RE0571C 	 1
    C14STRr 	 0.0
    RE3503N 	 2.0
    xmpt 	 1
    DLNLCGCPT2 	 -1.0
    OXYP1CONJ 	 1
    RE3496N 	 3.0
    AGPex 	 1
    RE3392M 	 3.0
    GAL3ST12 	 -1.0
    r0712 	 -1.0
    r1925 	 -1.0
    HSAT2ly 	 1
    RE1899C 	 1
    r2164 	 -1.0
    r1662 	 3.0
    HEX7 	 3.0
    EX_gly_LPAREN_e_RPAREN_ 	 3
    NACHEXA14ly 	 -1.0
    TETPENT3COAtx 	 1
    EX_dad_2_LPAREN_e_RPAREN_ 	 3
    H2Otm 	 -1.0
    NADtm 	 1
    M8MASNterg 	 1
    LGNCCOAtx 	 1
    RE3564C 	 3.0
    RE2405C 	 -1.0
    CHOLPtl 	 1
    r2003 	 -1.0
    RE2445E 	 0.0
    EX_CE4881_LPAREN_e_RPAREN_ 	 1
    FACOAL1822 	 2.0
    r1301 	 1
    ST8SIA12 	 -1.0
    Uritm 	 -1.0
    NPTHLte 	 1
    DADA 	 -1.0
    r1464 	 1
    CTPS1 	 0.0
    SFCYSc 	 1
    FAOXC14DCC12DCx 	 -1.0
    r0445 	 1
    NS26T2g 	 -1.0
    EX_hyptaur_LPAREN_e_RPAREN_ 	 1
    CARIBUP_Sitr 	 1
    RE2992X 	 1.0
    RE0566C 	 1
    FACOAL1813 	 2.0
    DGAT 	 -1.0
    EX_4nph_LPAREN_e_RPAREN_ 	 1
    RE3250M 	 3.0
    RE3561M 	 1
    EX_5homeprazole_LPAREN_e_RPAREN_ 	 1
    r1299 	 1
    RE3162R 	 1
    C161CPT2 	 -1.0
    r2061 	 -1.0
    PROGLYPEPT1tc 	 -1.0
    MALSO3tm 	 -1.0
    CYTDt2r 	 0.0
    r0825 	 1
    r2205 	 -1.0
    r1826 	 -1.0
    DOCO13ECOAtxc 	 -1.0
    SPR 	 -1.0
    ACODA 	 -1.0
    EX_fvs_LPAREN_e_RPAREN_ 	 1
    FPGS4m 	 -1.0
    PIt7 	 -1.0
    EX_rib_D_LPAREN_e_RPAREN_ 	 3
    r2288 	 -1.0
    ASPGLUm 	 0.0
    XOLTRI27te 	 1
    MEOHtr 	 1
    RTOT_2 	 1
    S4TASE1ly 	 -1.0
    VLCS2p 	 -1.0
    r2195 	 -1.0
    r1005 	 1
    LINKDEG3ly 	 1
    THRS 	 -1.0
    PPMI12346Ptn 	 1
    SBCOAACOTx 	 -1.0
    3ISPVStep 	 1
    P4502A6 	 -1.0
    C60CRNt 	 -1.0
    RN0031X 	 0.0
    SMVACIDhep 	 3.0
    PSDm_hs 	 0
    UMPK4n 	 3.0
    EX_glysar_LPAREN_e_RPAREN_ 	 1
    ADRNCPT2 	 -1.0
    HEXDICOAACBP 	 3.0
    GCHOLAtx 	 1
    EX_fuc13galacglcgal14acglcgalgluside_hs_LPAREN_e_RPAREN_ 	 1
    RE3104R 	 1
    r1884 	 -1.0
    RE2875C 	 1
    OIVD2m 	 -1.0
    sink_c81coa_LPAREN_c_RPAREN_ 	 1
    FUCGALGBSIDEte 	 1
    GAMt1r 	 3.0
    r0809 	 1
    BDHm 	 0.0
    DESAT22_2p 	 1
    PGCD 	 3.0
    DPGase 	 3.0
    r2366 	 0
    APPNNte 	 1
    S6T22g 	 -1.0
    CLPNDCRNt 	 -1.0
    GLCAT5g 	 -1.0
    r0393 	 3.0
    EX_pydxn_LPAREN_e_RPAREN_ 	 3
    r1518 	 3.0
    biomass_lipid 	 3
    r2027 	 -1.0
    RE2814R 	 1
    RE2130C 	 1
    OCTDECCPT1 	 0
    GLCAT4g 	 -1.0
    HPDCAt 	 1
    RE1134M 	 -1.0
    56DHPVShc 	 3.0
    RE0456N 	 -1.0
    RE2410C 	 1
    EX_co2_LPAREN_e_RPAREN_ 	 3
    3HPCOAHYD 	 3.0
    RE0912C 	 0
    GLYKm 	 -1.0
    EX_sphs1p_LPAREN_e_RPAREN_ 	 1
    r0156 	 1
    MCLOR 	 3.0
    3SALACBOXL 	 0.0
    AK2LGCHOLt 	 1
    FAOXC12DCC10DCx 	 -1.0
    GALK 	 -1.0
    DOLPH_Uer 	 1
    4HMDGLUCitr 	 1
    ACNACNGALGBSIDEtg 	 1
    FUCFUCFUCGALACGLCGAL14ACGLCGALGLUSIDEte 	 1
    MDZGLChr 	 1
    EX_nrvnc_LPAREN_e_RPAREN_ 	 1
    LEUATB0tc 	 1.0
    P4508B13r 	 -1.0
    RE3393M 	 3.0
    r2265 	 -1.0
    UDPG1P 	 1
    DNDPt52m 	 -1.0
    LEUKTRC4t 	 -1.0
    RETH1 	 1
    EX_avite2_LPAREN_e_RPAREN_ 	 1
    r1819 	 -1.0
    TACRtu 	 -1.0
    SPRn 	 -1.0
    RE3563X 	 -1.0
    RE2149C 	 -1.0
    PPCDC 	 -1.0
    RE2996X 	 -1.0
    r1720 	 -1.0
    NNATm 	 -1.0
    RE2635C 	 1
    EX_citr_L_LPAREN_e_RPAREN_ 	 1
    6OHFVShc 	 -1.0
    RE0908C 	 1
    EX_12htacr_LPAREN_e_RPAREN_ 	 1
    PDHm 	 -1.0
    r0502 	 -1.0
    GGH_6THFe 	 1.0
    OCOAT1m 	 1.0
    3M4HDXPAC 	 2.0
    P45011B21m 	 -1.0
    ALDD2xm 	 3.0
    DNDPt39m 	 -1.0
    RE1834M 	 1.0
    4HBZCOAFm 	 1
    4HMDGLUCtev 	 1
    EX_tmndnc_LPAREN_e_RPAREN_ 	 1
    MLTG1e 	 -1.0
    AG13T8g 	 1.0
    UDPGALtg 	 -1.0
    RE1702C 	 1
    EX_prgstrn_LPAREN_e_RPAREN_ 	 1
    ESTRIOLtr 	 1
    r2067 	 -1.0
    CSPG_Btly 	 1
    RE3404M 	 3.0
    r1455 	 1
    2HATVLACGLUCitr 	 1
    RE0582N 	 1
    FACOAL40im 	 -1.0
    r1871 	 -1.0
    EX_dopasf_LPAREN_e_RPAREN_ 	 1
    6EPVSthc 	 -1.0
    r0403 	 -1.0
    RE3079C 	 1
    glyc3pte 	 1
    CGMPt 	 1.0
    BZtr 	 1
    RE2221C 	 1
    TMACMPitr 	 1
    ASNPHELAT2tc 	 -1.0
    r0986 	 1
    RE0922R 	 0
    RE3432X 	 0
    EX_progly_LPAREN_e_RPAREN_ 	 1
    FUCFUC12GAL14ACGLCGALGLUSIDEte 	 1
    r1681 	 -1.0
    EX_asn_L_LPAREN_e_RPAREN_ 	 3
    P4504F121r 	 -1.0
    EX_biocyt_LPAREN_e_RPAREN_ 	 1
    EX_HC02201_LPAREN_e_RPAREN_ 	 1
    FRDPtc 	 1
    r2406 	 -1.0
    C141OHc 	 0
    DHCRD2 	 3.0
    SELCYSLY 	 -1.0
    RE1527C 	 1
    GLYt7_211_r 	 -1.0
    AM1ALCStep 	 1
    r0942 	 2.0
    ADSS 	 0.0
    DDCRNe 	 1
    DINt 	 -1.0
    RE3289R 	 -1.0
    GLUB0AT3tc 	 -1.0
    RE0921R 	 0
    FAOXC10080m 	 3.0
    ASNtm 	 1
    r2176 	 -1.0
    PUNP6 	 2.0
    5MTHFt2 	 -1.0
    SPHS1Pte 	 1
    NACASPtm 	 1
    ETHAK 	 1.0
    r1643 	 3.0
    CYSTHRNaEx 	 3.0
    EX_crglz_LPAREN_e_RPAREN_ 	 1
    r2327 	 3.0
    r2362 	 0
    DNDPt58m 	 -1.0
    EX_HC01700_LPAREN_e_RPAREN_ 	 1
    6MELACAChep 	 -1.0
    SCP2x 	 3.0
    r1648 	 3.0
    EX_inost_LPAREN_e_RPAREN_ 	 3
    r2000 	 -1.0
    SELCYSTS 	 -1.0
    ESTRIOLGLCte 	 1.0
    ANTHte 	 1
    AHEXASEly 	 1.0
    FVSCOAhc 	 1
    RE3168R 	 1
    RE3218L 	 -1.0
    CRVSM24teb 	 0.0
    ATVACIDOATPtu 	 -1.0
    IMACTD_m 	 3.0
    Htr 	 1
    OIVD1m 	 -1.0
    RE3235C 	 1
    RE3386M 	 3.0
    EX_fe3_LPAREN_e_RPAREN_ 	 3
    r1630 	 3.0
    NDPK8 	 0
    GGLUCT 	 3.0
    r1779 	 -1.0
    r1916 	 -1.0
    AMETt2m 	 -1.0
    Uritl 	 -1.0
    PGDIr 	 2.0
    FAOXC163C142x 	 -1.0
    FAS80COA_L 	 -1.0
    SMVLAChep 	 1
    r0127 	 -1.0
    NDPK9 	 0
    FAOXC204C184m 	 3.0
    ADNtl 	 -1.0
    RE3161C 	 3.0
    LSTNM1tev 	 1
    1OHMDZtep 	 1
    BDG2HCGHD 	 1
    EX_dtdp_LPAREN_e_RPAREN_ 	 1
    3HSMVACIDhep 	 -1.0
    r2228 	 -1.0
    LNELDCCPT2 	 -1.0
    ST3GAL31g 	 3.0
    LSTNtd 	 1
    RE3343X 	 -1.0
    EX_6hmsmvacid_LPAREN_e_RPAREN_ 	 1
    HACD9m 	 3.0
    FACOAL1812 	 2.0
    EX_q10_LPAREN_e_RPAREN_ 	 1
    r2276 	 -1.0
    EX_tsacmgluc_LPAREN_e_RPAREN_ 	 1
    RE1815C 	 1
    GLNtm 	 1
    TMDPPK 	 1
    FACOAL203 	 2.0
    r0779 	 3.0
    r2420 	 -1.0
    r0870 	 1
    r2120 	 0.0
    RE1134C 	 1
    GALFUCGALACGLCGAL14ACGLCGALGLUSIDEte 	 1
    UDPXYLtg 	 0.0
    PRISTANALtx 	 1
    RE1635C 	 1
    r1564 	 3.0
    r0732 	 3.0
    TYROXDAc 	 3.0
    r2199 	 -1.0
    EX_HC02220_LPAREN_e_RPAREN_ 	 1
    r0656 	 3.0
    ACSOMT 	 -1.0
    r1525 	 3.0
    MAOLNOR 	 3.0
    RE2346C 	 1
    RE2876C 	 1
    B3GNT35g 	 -1.0
    DARGOp 	 -1.0
    EX_carveol_LPAREN_e_RPAREN_ 	 1
    31DMTtu 	 1
    ALATA_L 	 -1.0
    EX_HC02200_LPAREN_e_RPAREN_ 	 1
    HISTASE 	 -1.0
    FAOXC226C205m 	 3.0
    CITRtm 	 -1.0
    r0545 	 2.0
    DESFVSitr 	 1
    OXAHCOtex 	 0.0
    MM6B1ag 	 3.0
    GLCAT3g 	 -1.0
    r1950 	 -1.0
    SIAASE 	 -1.0
    ACGALK2 	 1
    r1165 	 -1.0
    RE3036N 	 1
    EX_tmd_LPAREN_e_RPAREN_ 	 1
    OXYPR1tehv 	 1
    FAOXC2442246x 	 0.0
    SPMS 	 -1.0
    GALASE15ly 	 0
    r2386 	 -1.0
    RE2360N 	 -1.0
    HAtly 	 1
    SERTHRNaEx 	 3.0
    MANt4 	 -1.0
    LPS2 	 -1.0
    AG13T1g 	 1.0
    RE2524X 	 1.0
    H2O2itr 	 1
    PYDAMtr 	 1
    FDH 	 -1.0
    r2304 	 -1.0
    MDZtd 	 1
    RE3307X 	 0
    FAOXC161C161OHm 	 3.0
    RE1531X 	 3.0
    ADEtl 	 -1.0
    C12DCe 	 1
    FPGS2m 	 -1.0
    5OHFVSitr 	 1
    GALt4 	 -1.0
    FACOAL200 	 2.0
    RE0549C 	 1
    EX_HC01361_LPAREN_e_RPAREN_ 	 1
    DDECE1CRNe 	 1
    DM_dgtp_n_ 	 3
    r1800 	 -1.0
    EX_glyphe_LPAREN_e_RPAREN_ 	 1
    DHEAtr 	 1
    PRDXl 	 -1.0
    FAOXC141C121x 	 -1.0
    25HVITD3tm 	 1
    PTVSTM13te 	 1
    S6T21g 	 1.0
    AM1A4NCShc 	 1
    NDP8ex 	 3.0
    BUP2 	 -1.0
    SIAT4Bg 	 -1.0
    2HATVACIDGLUChr 	 0
    EX_allop_LPAREN_e_RPAREN_ 	 1
    EX_fmn_LPAREN_e_RPAREN_ 	 3
    ADPGLC 	 1.0
    FAOXC120100m 	 3.0
    PYDXNK 	 3.0
    DM_dttp_m_ 	 1
    MGCHrm 	 -1.0
    RE3394M 	 3.0
    RE3476M 	 0
    SALMCOM 	 3.0
    r2444 	 1
    THMMPtm 	 1
    ILEt5m 	 1
    r1867 	 -1.0
    2AMACHYD 	 -1.0
    35DSMVteb 	 1
    r1805 	 -1.0
    IDHPOXOX3 	 -1.0
    r2117 	 0.0
    6MSMVhep 	 3.0
    RE2541L 	 3.0
    SPHMYLNtg 	 1
    GGH_10FTHF5GLUl 	 1.0
    RE2635R 	 -1.0
    VITD2t 	 1
    RE2474C 	 -1.0
    LNLCCPT2 	 -1.0
    PPDOy 	 3.0
    RE2454M 	 3.0
    DPMVDx 	 -1.0
    RE3175R 	 1
    G14T6g 	 1.0
    EX_24nph_LPAREN_e_RPAREN_ 	 1
    DAGt 	 1
    EX_HC01441_LPAREN_e_RPAREN_ 	 1
    GLCt4 	 3.0
    EX_tetdec2crn_ 	 1
    GTACMPitr 	 1
    UGT1A4r 	 0
    RE2150R 	 0
    r1821 	 -1.0
    PCLAD 	 -1.0
    RE1826M 	 -1.0
    GTPtn 	 1
    RE2814M 	 1
    GLCNACASE2ly 	 -1.0
    7THFtl 	 1
    RE2880C 	 1
    EX_3octdece1crn_ 	 1
    r1982 	 -1.0
    RE2131C 	 1
    r1117 	 1
    r0084 	 2.0
    RE1817M 	 1.0
    r0573 	 -1.0
    r0431 	 -1.0
    RE1897C 	 1
    IMACTD 	 3.0
    PUNP4 	 2.0
    RE3490C 	 1
    DIDPtn 	 1
    RE1906R 	 -1.0
    FAOXC181C161m 	 3.0
    FAOXC18C18OHm 	 3.0
    MMCD 	 -1.0
    EX_coa_LPAREN_e_RPAREN_ 	 1
    r0913 	 -1.0
    RE2640C 	 -1.0
    B3GNT315g 	 -1.0
    r2162 	 -1.0
    FAOXC102C103x 	 0.0
    RE0916R 	 -1.0
    GQ1BALPHAtg 	 1
    LGNCCRNt 	 -1.0
    AGPRim 	 1.0
    TACRDtsc 	 1
    r1716 	 -1.0
    r1051 	 1
    PGDI 	 2.0
    r0016 	 0
    ATVACIDtdu 	 1
    FAOXC12C12OHm 	 2.0
    HSD3B13r 	 -1.0
    CRTSTRNtm 	 1
    PDE1 	 1.0
    LCYSTAT 	 1.0
    TMDtd 	 1
    EX_3ddcrn_ 	 1
    r1636 	 3.0
    ADK3m 	 3.0
    GLC3MEACPitr 	 1
    S4T4g 	 1
    3HSMVACIDteb 	 1
    H2CO3Dm 	 -1.0
    RE2995X 	 2.0
    MDHm 	 3.0
    RE1525C 	 1
    EX_prostgd2_LPAREN_e_RPAREN_ 	 1
    BACCL 	 -1.0
    r0777 	 2.0
    EX_fucacngal14acglcgalgluside_hs_LPAREN_e_RPAREN_ 	 1
    RE2768M 	 1
    EX_estradiolglc_LPAREN_e_RPAREN_ 	 1
    AM19CSALThr 	 3.0
    EX_lstnm7_LPAREN_e_RPAREN_ 	 1
    EX_mdzglc_LPAREN_e_RPAREN_ 	 1
    ACCOAL 	 0.0
    r2218 	 -1.0
    RE1447N 	 1
    r2313 	 3.0
    DNDPt20m 	 -1.0
    OAGD3te 	 1
    PUNP3 	 2.0
    r1170 	 -1.0
    ALDSTRNtm 	 1
    LSTNM5tev 	 1
    ARACHDCOAtx 	 1
    FAOXC102C81m 	 0.0
    EX_fucacgalfucgalacglcgalgluside_hs_LPAREN_e_RPAREN_ 	 1
    FAOXC122m 	 1.0
    r2023 	 -1.0
    r2035 	 -1.0
    PHEt4 	 3.0
    PCRNtm 	 -1.0
    r0552 	 1.0
    r1865 	 -1.0
    TMABADH 	 3.0
    r2279 	 -1.0
    G14T20g 	 1.0
    GALASE16ly 	 0
    ESTRADIOLtr 	 1
    HACD1m 	 3.0
    r2095 	 0.0
    AMETr 	 1
    P450SCC1m 	 1
    RE1836C 	 1
    r2493 	 -1.0
    DURIt 	 -1.0
    CHSTEROLtrc 	 1
    S6TASE20ly 	 3.0
    GALGALFUCFUCGALACGLCGALACGLCGAL14ACGLCGALGLUSIDEtg 	 1
    r1702 	 -1.0
    RE1135R 	 -1.0
    DMANTIPYRINEte 	 1
    DM_avite2_c_ 	 1
    VALB0AT2tc 	 1
    KCCt 	 3.0
    r1559 	 3.0
    IDHPOXOX4 	 -1.0
    C100CRNt 	 -1.0
    r0826 	 1
    HPDCACRNt 	 -1.0
    EX_HC02210_LPAREN_e_RPAREN_ 	 1
    GALT 	 -1.0
    r0409 	 1
    ASP1DC 	 -1.0
    EX_34hpp_ 	 3
    AHCYStr 	 1
    r1855 	 -1.0
    r1983 	 -1.0
    H2O2syn 	 3.0
    FACOAL160i 	 2.0
    HS1ly 	 1.0
    FUCFUCGALACGLCGALGLUSIDEte 	 1
    RE3534R 	 3.0
    AATAi 	 -1.0
    EX_dhdascb_LPAREN_e_RPAREN_ 	 3
    TETDECE1CRNe 	 1
    ADNtm 	 -1.0
    EX_14hmdz_LPAREN_e_RPAREN_ 	 1
    SEAHCYSHYD 	 2.0
    CMPACNAtg 	 3.0
    TMLYSter 	 1
    EX_3hdececrn_ 	 1
    O2ter 	 1
    COUMARINte 	 1
    r0961 	 -1.0
    CSNATer 	 0.0
    FAOXC11C9m 	 0.0
    ADNt4 	 0.0
    SAMHISTA 	 0.0
    FAOXC122_3Z_6Zm 	 1.0
    r0357 	 3.0
    r2143 	 -1.0
    TSTSTERONEGLCtr 	 1
    SPRMS 	 -1.0
    EX_4hatvacid_LPAREN_e_RPAREN_ 	 1
    RE2874C 	 1
    RE1448N 	 1
    RE1309M 	 -1.0
    PTDCACRNCPT1 	 0
    GQ1Btg 	 1
    34DHPHEt 	 -1.0
    TAXOLte 	 1
    RE3150R 	 1
    FE2tm 	 1
    r2062 	 -1.0
    HDCAter 	 1
    AM1CGLCteb 	 1
    CYTK9n 	 3.0
    r1457 	 1
    S6T9g 	 -1.0
    r1578 	 3.0
    EX_glcur_LPAREN_e_RPAREN_ 	 3
    HCO3_CLt 	 -1.0
    ALKP 	 -1.0
    HEXCCPT1 	 0
    RE1100G 	 -1.0
    MESCOALm 	 1
    AKR1C41 	 -1.0
    ACALDtm 	 1
    7BHGLZABCt 	 1
    LEULEUPEPT1tc 	 -1.0
    r0363 	 3.0
    DADNt4 	 -1.0
    RE1526C 	 1
    OAGT3te 	 1
    r2270 	 -1.0
    GLYitr 	 1
    PIt8 	 1.0
    RE3391M 	 3.0
    r1606 	 3.0
    RE1829C 	 1
    RE3228R 	 1
    RE2870C 	 1
    SEAHCYStn 	 1
    LEUB0AT2tc 	 1
    r2250 	 -1.0
    TDPDRE 	 1
    ORETNF 	 1
    r1834 	 -1.0
    S6T25g 	 -1.0
    MMCDm 	 -1.0
    EX_fuc14galacglcgalgluside_hs_LPAREN_e_RPAREN_ 	 1
    r0729 	 -1.0
    RE3095C 	 3.0
    DASCBH 	 1
    r1184 	 1.0
    EX_4bhglz_LPAREN_e_RPAREN_ 	 1
    24_25VITD2Hm 	 -1.0
    STRDNCCOAtxc 	 -1.0
    FACOAL180i 	 2.0
    PPAer 	 1.0
    r0157 	 1
    ACACT10m 	 3.0
    C8DCc 	 0
    ATPasel 	 1.0
    S4TASE2ly 	 -1.0
    ACOAD10m 	 2.0
    r1526 	 3.0
    r1858 	 -1.0
    EX_tdchola_LPAREN_e_RPAREN_ 	 1
    MDH 	 3.0
    r0787 	 3.0
    r2425 	 -1.0
    UMPK3 	 3.0
    RE1520M 	 3.0
    EX_prostgh2_LPAREN_e_RPAREN_ 	 3
    C14OHc 	 0
    FAOXC182_9Z_12Zm 	 3.0
    RE2914X 	 3.0
    r2105 	 1
    PNTOt5 	 -1.0
    FUCACGALFUCGALACGLCGALGLUSIDEtg 	 1
    PEPCKm 	 1.0
    M13N4Tg 	 3.0
    GUAD 	 -1.0
    XOLTRI24tc 	 1
    EX_idour_LPAREN_e_RPAREN_ 	 1
    r1154 	 -1.0
    6BHGLZGLCitr 	 1
    EX_pnto_R_LPAREN_e_RPAREN_ 	 3
    MAN6PI 	 -1.0
    EX_estrones_LPAREN_e_RPAREN_ 	 1
    RE3233G 	 -1.0
    EPOXTACteb 	 1
    r1637 	 3.0
    35DSMVitr 	 1
    r2291 	 -1.0
    4NPHSULT 	 1.0
    EX_xmp_LPAREN_e_RPAREN_ 	 3
    EX_4hphac_LPAREN_e_RPAREN_ 	 1
    RE2132C 	 1
    r1512 	 1
    FAOXC205_5Z_8Z_11Z_14Z_17Zm 	 3.0
    RE2658C 	 1
    G1M8MASNterg 	 1
    EX_CE1940_LPAREN_e_RPAREN_ 	 1
    SRTNENT4tc 	 -1.0
    RE3383M 	 3.0
    r1383 	 1.0
    PCRNtc 	 1
    LYSt4 	 1.0
    r1602 	 3.0
    EX_stacmp_LPAREN_e_RPAREN_ 	 1
    r1706 	 -1.0
    EPOXTACitr 	 1
    r2379 	 -1.0
    r0994 	 1
    ETFQO 	 0.0
    r1162 	 1
    PRISTCOAtcx 	 2.0
    FAOXC227C226m 	 3.0
    GLYCLTDy 	 2.0
    r0760 	 -1.0
    11DOCRTSLtr 	 1
    RE3171R 	 1
    FAOXC184m 	 1.0
    ESTRIOLGLCtr 	 1
    AM1ACStep 	 1
    RE1539C 	 1
    CSEPASEly 	 1
    EX_psyltchol_LPAREN_e_RPAREN_ 	 1
    SELCYSTGL 	 -1.0
    ARACHCRNt 	 -1.0
    r2483 	 -1.0
    INSTt4_2 	 -1.0
    EX_cys_L_LPAREN_e_RPAREN_ 	 3
    SPHS1Ptr 	 1
    ADKd 	 -1.0
    TMDM5hr 	 -1.0
    RE3010X 	 0
    r0651 	 -1.0
    T2M26DCOAHLm 	 3.0
    P45046A1r 	 -1.0
    RE2917M 	 3.0
    r2314 	 3.0
    RE3563M 	 3.0
    FAOXC184_6Z_9Z_12Z_15Zm 	 3.0
    r1064 	 1
    RE2149R 	 0
    GP1CALPHAtg 	 1
    3OHACMPhr 	 -1.0
    SMVACIDtev 	 1
    HSAT3ly 	 1
    RE2888N 	 -1.0
    3SPYRSPm 	 1
    GULNter 	 1
    r1595 	 3.0
    ESTRONEGLCt 	 -1.0
    ALAt2rL 	 -1.0
    r1787 	 -1.0
    TSACMSULtev 	 1
    EX_Rtotal_LPAREN_e_RPAREN_ 	 1
    TYRt4 	 1.0
    LEUKTRF4t 	 0
    4HATVACIDhc 	 -1.0
    ABO5g 	 -1.0
    r2220 	 -1.0
    41R2A1H12BOOX 	 3.0
    ARGtm 	 -1.0
    EX_fucfuc132galacglcgal14acglcgalgluside_hs_LPAREN_e_RPAREN_ 	 1
    EX_HC02191_LPAREN_e_RPAREN_ 	 1
    AMACRp 	 -1.0
    AM4N9CStev 	 1
    r2082 	 0.0
    RE3597C 	 1.0
    FAOXC225C204m 	 3.0
    THMDt4 	 0.0
    BGLUTDECHOe 	 1
    RE3241R 	 3.0
    ASNTHRNaEx 	 -1.0
    r0432 	 0.0
    PMI1346PHn 	 1.0
    5ADTSTSTERONEGLCte 	 1.0
    r0638 	 -1.0
    RE2404C 	 -1.0
    4PYRDX 	 1
    EX_am4ncs_LPAREN_e_RPAREN_ 	 1
    EX_ttdca_LPAREN_e_RPAREN_ 	 3
    C30CPT1 	 0
    MDRPD 	 1
    ACOAD9m 	 2.0
    r2071 	 -1.0
    PCHOLHSTDe 	 1
    RE2909X 	 3.0
    IDPtn 	 1
    FAOXC163C164x 	 0.0
    PI5P3Ker 	 3.0
    B3GNT34g 	 -1.0
    r2248 	 -1.0
    r2289 	 -1.0
    CBPS 	 0.0
    RE3234R 	 1
    NDPK4 	 0
    r1001 	 1
    RE0567C 	 1
    RE2632C 	 1
    NACHEX19ly 	 1.0
    RE3075X 	 -1.0
    ACONTm 	 3.0
    ARSA 	 -1.0
    RE1096R 	 -1.0
    r0464 	 3.0
    r0553 	 1.0
    r1721 	 -1.0
    r0836 	 -1.0
    FACOAL170 	 2.0
    RE1582C 	 1
    FUMtm 	 -1.0
    RE3536C 	 1
    r0068 	 0.0
    FACOAL241 	 2.0
    EX_HC00342_LPAREN_e_RPAREN_ 	 1
    ADSL2 	 3.0
    r1497 	 -1.0
    RE1938C 	 1
    MM7Cag 	 3.0
    ANTIPYRENEte 	 1
    r1778 	 -1.0
    FAOXC161C141m 	 3.0
    NICRNTtn 	 1
    NACHEX16ly 	 1.0
    FTHFLm 	 2.0
    AM4N9CShc 	 1
    r1495 	 -1.0
    6DHFtm 	 1
    ACETONEt2m 	 0.0
    EX_HC01577_LPAREN_e_RPAREN_ 	 1
    r0990 	 1
    r0472 	 -1.0
    CBPter 	 1
    PI345P3P 	 3.0
    PPItr 	 1
    r0028 	 3.0
    r1556 	 3.0
    FAOXC181_9Em 	 3.0
    EX_adpcbl_LPAREN_e_RPAREN_ 	 3
    5HTRPDOX 	 -1.0
    RE3301R 	 3.0
    RE2319R 	 3.0
    r0208 	 0
    SERB0AT3tc 	 -1.0
    SBPP1er 	 0.0
    DHEASULT 	 0.0
    UDPGP 	 1.0
    FBP26 	 3.0
    r1474 	 -1.0
    CHSTEROLt2 	 1
    RE3121R 	 1
    CYTK14 	 3.0
    r1599 	 3.0
    RE0573N 	 1
    APOCF 	 1
    r2357 	 0
    LEUKTRB4tr 	 1
    INSt 	 -1.0
    GLYGLYCNc 	 3.0
    r0733 	 3.0
    EX_dgchol_LPAREN_e_RPAREN_ 	 1
    r0974 	 1
    GALNACT2g 	 0.0
    r0683 	 1
    AM1CGLCitr 	 1
    ASPDt6 	 1.0
    EX_doco13ac_ 	 1
    FOLt2le 	 1
    ASNCYSNaEx 	 -1.0
    RN0022X 	 0.0
    C101CPT1 	 0
    RE2078M 	 1
    MERCPLACCYSt 	 1
    r2246 	 -1.0
    CYStec 	 -1.0
    LINKDEG1ly 	 1
    RE3231C 	 1
    M4CET3er 	 -1.0
    SLDt 	 1
    FAOXC204_5Z_8Z_11Z_14Zm 	 3.0
    RE3450C 	 1
    r0925 	 1
    DESAT18_6 	 -1.0
    RE2405R 	 0
    S6T5g 	 -1.0
    RE3147R 	 1
    EX_rbt_LPAREN_e_RPAREN_ 	 1
    IZPN 	 -1.0
    UPPN 	 -1.0
    GTHPm 	 3.0
    ARACHDFATPtc 	 -1.0
    1531TALThr 	 1
    HSD17B4x 	 3.0
    DESAT16_2 	 3.0
    r2499 	 -1.0
    FAOXC162C162OHm 	 3.0
    7DHCHSTEROLtr 	 1
    NADS2 	 1.0
    B3GALTg 	 1
    ADRNCPT1 	 0
    FUT92g 	 -1.0
    r2144 	 -1.0
    RE0875C 	 1
    ACMPdt 	 1
    RE3232R 	 1
    DCSPTN1t 	 1
    EX_cynt_LPAREN_e_RPAREN_ 	 1
    r2283 	 -1.0
    PDE1g 	 -1.0
    r1318 	 1
    r1773 	 -1.0
    NH4t3r 	 -1.0
    ACACtx 	 1
    LSTNtu 	 -1.0
    NACHEX21ly 	 1.0
    4MOPte 	 1
    G6Pter 	 -1.0
    GLXtm 	 1
    DNDPt50m 	 -1.0
    FAOXC9C7m 	 0.0
    INSK 	 1
    RDH3a 	 3.0
    CSPG_Atly 	 1
    LINOFATPtc 	 1
    S6TASE24ly 	 3.0
    G1M7MASNCterg 	 1
    ATVACIDitr 	 1
    ESTRONEGLCtr 	 1
    EX_lnlnca_LPAREN_e_RPAREN_ 	 3
    RE0578M 	 1.0
    LSTNM7hr 	 0
    ILEPHELAT2tc 	 1
    r1678 	 -1.0
    RE1628C 	 1
    RE3083X 	 -1.0
    RE3518C 	 1
    RE0568C 	 1
    DOLPMT3_Uer 	 1
    EX_ptvstm13_LPAREN_e_RPAREN_ 	 1
    r0818 	 1
    RE0579C 	 0
    G14T5g 	 1.0
    r2409 	 -1.0
    14HMDZhr 	 3.0
    SUCCCROT 	 3.0
    r1384 	 1
    TXA2te 	 1
    RE3005M 	 3.0
    RE2633R 	 -1.0
    FUT33g 	 1.0
    r0364 	 3.0
    EX_1hmdgluc_LPAREN_e_RPAREN_ 	 1
    OMHPALTD 	 1
    NAHCO3_HCLt 	 -1.0
    RE2111M 	 -1.0
    r2295 	 -1.0
    RE2948C 	 1
    PAPtg 	 1
    4HATVLAChc 	 1
    DPCOAPPe 	 1
    RE1956X 	 0.0
    NRVNCCPT1 	 0
    EX_smv_LPAREN_e_RPAREN_ 	 1
    G3PD2m 	 1.0
    3AIBt 	 1
    RE3251C 	 1
    r1570 	 3.0
    r2539 	 1
    EX_avite1_LPAREN_e_RPAREN_ 	 1
    LVSTACIDtu 	 -1.0
    DHCR241r 	 3.0
    RE3252C 	 1
    FALDtm 	 1
    IDOAASE2ly 	 -1.0
    15DMTtep 	 1
    DAGK_hs 	 3.0
    DHORD9 	 -1.0
    SACCD4m 	 0.0
    r1529 	 3.0
    EX_3deccrn_ 	 1
    r2216 	 -1.0
    RE3446X 	 2.0
    UGT1A2r 	 0
    EX_2425dhvitd3_LPAREN_e_RPAREN_ 	 1
    FAOXC142_5E_8Em 	 3.0
    OCD11COACPT1 	 1
    RE3114R 	 2.0
    r0707 	 -1.0
    UDPXYLter 	 1
    TMDS 	 1.0
    AASAD3m 	 1
    r2145 	 -1.0
    FAOXTC101TC102m 	 2.0
    r0587 	 -1.0
    r1653 	 3.0
    OCDCEAtr 	 -1.0
    ST6GALNAC62 	 -1.0
    r0778 	 2.0
    S2T3g 	 1.0
    NKCCt 	 0.0
    TMDPP 	 1.0
    EX_so3_LPAREN_e_RPAREN_ 	 1
    EX_ppa_LPAREN_e_RPAREN_ 	 3
    RE2912X 	 2.0
    RE2513C 	 3.0
    FVSCOAitx 	 1
    GLYt2r 	 0
    LGTHL 	 3.0
    r0800 	 1
    FAOXC101x 	 1.0
    PPOR 	 -1.0
    FUT99g 	 -1.0
    PIt9 	 1.0
    RE0512X 	 3.0
    RE2026C 	 1
    BPNT 	 -1.0
    FPGS6m 	 -1.0
    LSTNRATt 	 1
    3HPVSTETCOAhcx 	 -1.0
    EX_camp_LPAREN_e_RPAREN_ 	 1
    RE3132R 	 2.0
    ASCBOX 	 1
    HIShPTtc 	 -1.0
    GLYATm 	 -1.0
    BGLUGCHe 	 1
    LNELDCt 	 1
    CARPEPT1tc 	 -1.0
    r2303 	 -1.0
    EX_3ump_LPAREN_e_RPAREN_ 	 1
    RE2562C 	 1
    HSD17B42x 	 3.0
    COAtm 	 3.0
    FPGS5 	 -1.0
    RE0690C 	 0.0
    RAI4 	 1
    CHLtm 	 1
    RE2858C 	 1
    RE3180C 	 1
    NACHEXA6ly 	 -1.0
    GLYAMDTRc 	 3.0
    FUMSO4tm 	 -1.0
    GPIAT 	 -1.0
    MM6B1bg 	 3.0
    r2261 	 -1.0
    DEOXFVStev 	 1
    r2108 	 1
    RE3192M 	 3.0
    r0236 	 1
    RE3513R 	 1
    6HMSMVACIDhep 	 -1.0
    URIt4 	 0.0
    RE0572N 	 1
    EX_fe2_LPAREN_e_RPAREN_ 	 3
    NADtru 	 1
    r1631 	 3.0
    DNDPt62m 	 -1.0
    EX_dspvs_LPAREN_e_RPAREN_ 	 1
    r1730 	 -1.0
    RE2921M 	 3.0
    FUCGAL14ACGLCGALGLUSIDEte 	 1
    FAOXC141_5Zm 	 3.0
    4BHGLZABCt 	 1
    S6TASE13ly 	 3.0
    OIVD3m 	 -1.0
    SALMCOM2 	 1.0
    DIGALSGALSIDEtg 	 1
    r1715 	 -1.0
    FAOXC225C226x 	 0.0
    EX_5adtststeroneglc_LPAREN_e_RPAREN_ 	 1
    RE0344X 	 -1.0
    TMNDNCCRNt 	 -1.0
    ICDHy 	 2.0
    RE2999M 	 1.0
    THD1m 	 -1.0
    CYTK2 	 3.0
    ACMPGLUthc 	 -1.0
    ASCBSVCTtc 	 -1.0
    35DHPVStep 	 1
    RE1532M 	 3.0
    RSVSPONhc 	 -1.0
    SMVtu 	 -1.0
    ESTRONESt 	 2.0
    RE2852C 	 1
    FAOXC184C164x 	 -1.0
    RE1810M 	 1
    r0989 	 1
    EX_1hibupglu_S_LPAREN_e_RPAREN_ 	 1
    EX_asp_D_LPAREN_e_RPAREN_ 	 1
    r0267 	 1
    r1713 	 -1.0
    PSFLIP 	 2.0
    r2342 	 0
    TMDM5OATt 	 1
    ADRNCRNt 	 -1.0
    EX_CE5560_LPAREN_e_RPAREN_ 	 1
    IDOAASE1ly 	 -1.0
    FACOAL224 	 2.0
    DHPM2 	 -1.0
    HTDCACBP 	 3.0
    LEUKTRB4t 	 0
    r2410 	 -1.0
    RE3477C 	 1
    RE2877C 	 1
    RE1630C 	 1
    r2104 	 1
    CRTSTRNt 	 1
    r0165 	 -1.0
    RE2888C 	 3.0
    C180CPT2 	 -1.0
    GBSIDEtl 	 1
    EX_prostgi2_LPAREN_e_RPAREN_ 	 1
    CRVSM1hr 	 3.0
    EX_HC02196_LPAREN_e_RPAREN_ 	 1
    TMNDNCCOAtx 	 1
    FADtx 	 1
    HMGCOARc 	 3.0
    RE3015R 	 -1.0
    RE1817C 	 1
    DECDICRNe 	 1
    RE1943R 	 0.0
    ILETA 	 -1.0
    MM5bg 	 3.0
    LPASE 	 0.0
    r1918 	 -1.0
    H7MTer_U 	 1
    EX_phyQ_LPAREN_e_RPAREN_ 	 1
    DSAT 	 1
    RE2273E 	 0.0
    r2121 	 0.0
    ATVACIDMCTtu 	 0.0
    RE3090X 	 3.0
    RE3557M 	 1.0
    PI3P5K 	 -1.0
    RE1903R 	 1
    G6PDH2rer 	 -1.0
    MTHFR3 	 -1.0
    FAOXC102x 	 1.0
    DOLPGT1_Uer 	 1
    DUTPDPm 	 3.0
    RE2975C 	 1
    r2230 	 -1.0
    r2064 	 -1.0
    DURIPP 	 2.0
    ELAIDt 	 -1.0
    RE3340X 	 3.0
    EX_HC01446_LPAREN_e_RPAREN_ 	 1
    RIBt 	 1
    6BHGLZGLChr 	 1
    BALAPAT1tc 	 -1.0
    r0139 	 3.0
    FAOXC160140x 	 0.0
    r2398 	 -1.0
    r2167 	 -1.0
    r0504 	 -1.0
    RE2270C 	 1
    EX_2hatvacidgluc_LPAREN_e_RPAREN_ 	 1
    RE1573X 	 3.0
    2MCITt 	 1
    RE3521X 	 1.0
    BILGLCURtr 	 1
    GLBRAN 	 3.0
    EX_crvnc_LPAREN_e_RPAREN_ 	 3
    PAFS 	 1
    RE2910M 	 3.0
    DOPAtu 	 2.0
    SQLSr 	 3.0
    SGALSIDEtl 	 1
    CRNtim 	 -1.0
    H2Otn 	 1
    r1757 	 -1.0
    r2073 	 -1.0
    EX_acgalfucgalacgalfucgalacglcgal14acglcgalgluside_hs_LPAREN_e_RPAREN_ 	 1
    RE1812R 	 -1.0
    PPPG9tm 	 1
    NBAHH_ir 	 3.0
    RE0688C 	 0.0
    r2258 	 -1.0
    S4T1g 	 -1.0
    CITt4_4 	 -1.0
    LRAT 	 -1.0
    DHAAt1r 	 3.0
    EX_glyleu_LPAREN_e_RPAREN_ 	 1
    RE2398C 	 1
    FAOXC183806m 	 2.0
    NDERSVhc 	 0.0
    RE3346R 	 3.0
    VACCt 	 1
    r2305 	 -1.0
    RE2658R 	 1
    RE2051R 	 3.0
    MI13456Ptn 	 1
    EX_chol_LPAREN_e_RPAREN_ 	 3
    RE2768C 	 1
    r2044 	 -1.0
    3SALAOX 	 1
    LNLNCACPT2 	 -1.0
    r1768 	 -1.0
    RE2334C 	 1
    PRODt2r 	 -1.0
    EX_perillyl_LPAREN_e_RPAREN_ 	 1
    ATP1ter 	 1
    STRDNCt 	 1
    TRIPVSitr 	 1
    Am19CStev 	 1
    CAT2p 	 3.0
    FAOXC225x 	 1.0
    CITtam 	 -1.0
    RE3104C 	 1
    HPYRDCm 	 1
    r2239 	 -1.0
    r0633 	 -1.0
    GLNt4 	 3.0
    SACCD3m 	 0.0
    OMHDOCOSACTD 	 1
    CYTK7 	 3.0
    GGH_10FTHF6GLUl 	 1.0
    RE3106R 	 1
    ASPCTr 	 0.0
    r1850 	 -1.0
    FAOXC12DCTc 	 -1.0
    FAOXC226C227m 	 3.0
    RE3571C 	 1
    NACHEXA16ly 	 -1.0
    r1686 	 -1.0
    EX_mepi_LPAREN_e_RPAREN_ 	 1
    r1923 	 -1.0
    r1447 	 3.0
    PCFLOPm 	 -1.0
    ACACT8p 	 3.0
    MCITS 	 1
    HSD17B1 	 1.0
    FAOXC61m 	 1.0
    NRPPHRSULT 	 0
    r2380 	 -1.0
    EX_vitd2_LPAREN_e_RPAREN_ 	 1
    4MPTNLte 	 1
    GLCMter 	 1
    GTHOm 	 3.0
    FAOXC22C22DCHYr 	 -1.0
    r2007 	 -1.0
    r1854 	 -1.0
    LST4EXPitr 	 1
    r0178 	 0.0
    PPAP 	 3.0
    OCDCAFAPMtc 	 3.0
    RE3173C 	 1
    RN0023R 	 0.0
    r1440 	 1
    RE0924R 	 0
    r1953 	 -1.0
    AGLPC 	 1
    P4504F122r 	 -1.0
    D3AIBTm 	 -1.0
    r0531 	 1
    CYTK13 	 3.0
    MICITDr 	 1
    EX_HC02195_LPAREN_e_RPAREN_ 	 1
    AMY2e 	 0
    G6PPer 	 1.0
    r1109 	 1
    RE0690X 	 -1.0
    r0860 	 1
    PE_HStm 	 1
    SPC_HSt 	 1
    r2260 	 -1.0
    r2323 	 3.0
    BHBt 	 1
    13HTACRitr 	 1
    GT1Atg 	 1
    HISyLATthc 	 -1.0
    GLUDym 	 3.0
    r2168 	 -1.0
    RE1530M 	 -1.0
    r2244 	 -1.0
    RDH4 	 -1.0
    RE2444C 	 1
    MI134P4P 	 1
    3ISPVShc 	 1
    FAOXC181C181OHm 	 3.0
    6BHGLZtev 	 1
    r2403 	 -1.0
    GALSGLT1le 	 0.0
    r0840 	 1
    EX_taxol_LPAREN_e_RPAREN_ 	 1
    r2318 	 3.0
    ARACHCPT2 	 -1.0
    r0671 	 -1.0
    GGH_5DHFe 	 1.0
    PRO1x 	 -1.0
    r2022 	 -1.0
    r1891 	 -1.0
    B3GALT42g 	 -1.0
    IDOURte 	 1
    CORE5GTg 	 1
    r1722 	 -1.0
    r2159 	 -1.0
    r2396 	 -1.0
    QUILSYN 	 1
    EX_ade_LPAREN_e_RPAREN_ 	 3
    r0451 	 -1.0
    r1030 	 1
    EX_gt1a_hs_LPAREN_e_RPAREN_ 	 1
    GGH_5DHFl 	 1.0
    EX_prostge2_LPAREN_e_RPAREN_ 	 1
    SPHK21c 	 0.0
    EX_glc_LPAREN_e_RPAREN_ 	 3
    P4504F123r 	 0.0
    r1987 	 -1.0
    FUCACNGALACGLCGALGLUSIDEtg 	 1
    r0345 	 3.0
    EX_glyb_LPAREN_e_RPAREN_ 	 3
    BCDO 	 -1.0
    RE2382R 	 -1.0
    1a_25VITD2Hm 	 1
    SERALANaEx 	 3.0
    r2087 	 0.0
    B3GALT41g 	 -1.0
    RE2147C 	 -1.0
    FAOXC181_11Em 	 3.0
    FVSitr 	 1
    RN0029R 	 -1.0
    NTD4e 	 0.0
    COQ3m 	 -1.0
    r0941 	 1
    RE0938E 	 0.0
    RE3474R 	 -1.0
    SELt4_3 	 -1.0
    PHYTt 	 1
    A_MANASEly 	 -1.0
    RE2112R 	 -1.0
    DM_Lcystin 	 3
    CDPDAGtm 	 1
    SPHMDAc 	 1
    ADRNLPVESSEC 	 0.0
    EX_o2_LPAREN_e_RPAREN_ 	 3
    RE3154R 	 1
    FUCFUC132GALACGLCGAL14ACGLCGALGLUSIDEte 	 1
    r2092 	 1
    CRGLZhr 	 1
    NTD1m 	 -1.0
    NDPK8n 	 0
    CORE8GTg 	 1
    RE2973N 	 3.0
    NADPHtru 	 1
    AMPtr 	 1
    THRSERNaEx 	 3.0
    ANDRSTRNte 	 1
    r1548 	 3.0
    EX_ascb_L_LPAREN_e_RPAREN_ 	 1
    TMDM1OATt 	 1
    RE3129N 	 1
    RN0013C 	 1
    RE2152C 	 1
    EX_profvs_LPAREN_e_RPAREN_ 	 1
    r1329 	 1
    BILIRED 	 3.0
    URCN 	 -1.0
    RE2873C 	 1
    CLPNDt 	 1
    NDP10ex 	 3.0
    GLUt7l 	 1
    HPETFABP1tc 	 -1.0
    r2107 	 1
    r2127 	 0.0
    r1567 	 3.0
    r0618 	 0
    r1772 	 -1.0
    r2519 	 1
    FAOXC184C163m 	 3.0
    TMNDNCCPT1 	 0
    SUBEACTD 	 1
    APOCFm 	 1
    SUCD1m 	 -1.0
    1513TACRtu 	 1
    r2374 	 -1.0
    C12DCTD 	 1
    HAS1 	 2.0
    FAOXC120100x 	 0.0
    GLYBt4_2_r 	 -1.0
    r2212 	 -1.0
    VLCS2r 	 -1.0
    3AIBtm 	 1
    4MPTNLtm 	 1
    C100CPT2 	 -1.0
    PVSGLUChc 	 1
    r1934 	 -1.0
    r1968 	 -1.0
    EX_c4crn_ 	 1
    PTVSTLACtev 	 1
    CRVSM1teb 	 0.0
    r2131 	 0.0
    XOLEST2HSTDle 	 1
    EX_lcts_LPAREN_e_RPAREN_ 	 3
    RE3522C 	 1
    r2311 	 -1.0
    LCADi 	 3.0
    CYTK5 	 3.0
    ALCD2if 	 3.0
    RE3040R 	 0.0
    31DMTtep 	 1
    RE3124R 	 3.0
    CYSTS 	 -1.0
    r2058 	 -1.0
    FAOXC164C143m 	 3.0
    AMPTASECG 	 0.0
    r2004 	 -1.0
    CHOLATEt2 	 -1.0
    RE3014C 	 1
    r1790 	 -1.0
    METAT 	 3.0
    RE1830M 	 -1.0
    EX_HC00250_LPAREN_e_RPAREN_ 	 3
    L_LACDcm 	 -1.0
    RE1523X 	 -1.0
    r1077 	 1
    FACOAL1821 	 2.0
    FAOXC2051843m 	 2.0
    RE2521C 	 1
    1531TACRteb 	 1
    LSTN1GLUCtev 	 1
    r0119 	 3.0
    ST3GAL22g 	 -1.0
    FUCFUCGALACGLCGALGLUSIDEtg 	 1
    H2OGLYAQPt 	 3.0
    2OXOADOXm 	 -1.0
    AKGt4_3 	 -1.0
    NCCt 	 -1.0
    r1913 	 -1.0
    r1905 	 -1.0
    CYTK9 	 3.0
    FAOXC142C142OHm 	 3.0
    LNLCCRNt 	 1
    FUT15g 	 -1.0
    CHOLESACATc 	 -1.0
    FUCFUCFUCGALACGLC13GALACGLCGAL14ACGLCGALGLUSIDEtg 	 1
    FADDPle 	 0.0
    6AHGLZtev 	 1
    AGMTm 	 -1.0
    CYTDK1 	 0
    r0191 	 1
    LNLNCGCPT1 	 0
    PI45P4P 	 1
    r0517 	 1.0
    MM5ag 	 3.0
    EX_antipyrene_LPAREN_e_RPAREN_ 	 1
    6BHGLZGLCABCt 	 1
    UGT1A3r 	 0
    FVShc 	 1
    THFtm 	 1
    FAOXC102C81x 	 -1.0
    AHCYStd 	 1
    P4507A1r 	 -1.0
    CHTNASE 	 -1.0
    PPA2m 	 1
    r1900 	 -1.0
    r1738 	 -1.0
    DOPAMT 	 3.0
    r0559 	 -1.0
    SERt4 	 3.0
    r1674 	 -1.0
    TRPHYDRO2 	 -1.0
    HPCLx 	 1
    H2Otly 	 1
    PIK3n 	 1
    MM8Ber 	 -1.0
    r2180 	 -1.0
    FKYNH 	 -1.0
    DCSPTN1CRNt 	 -1.0
    FAOXC121x 	 1.0
    RE2850C 	 1
    EX_pydx_LPAREN_e_RPAREN_ 	 3
    r0381 	 1
    RE1100L 	 -1.0
    OAGT3tg 	 1
    FAOXC22C20x 	 -1.0
    HSD3B11 	 -1.0
    SELMETAT 	 3.0
    GCC2cm 	 -1.0
    r2319 	 3.0
    RE2269E 	 -1.0
    RE2848C 	 1
    DGSNt 	 -1.0
    EX_lvst_LPAREN_e_RPAREN_ 	 1
    r1873 	 -1.0
    AM1C9CStev 	 1
    r1939 	 -1.0
    r2377 	 -1.0
    B3GALT43g 	 -1.0
    EX_meoh_LPAREN_e_RPAREN_ 	 3
    OCD11CRNCPT2 	 1
    O2tn 	 1
    24_25DHVITD2t 	 1
    6OHFVSitr 	 1
    HPYRDC 	 1
    FAOXC102C101m 	 3.0
    r0655 	 -1.0
    r2381 	 -1.0
    RE3367C 	 0.0
    SPMDOX 	 1
    PIK4 	 0.0
    PPItm 	 1
    SMPD3g 	 0.0
    ARGSS 	 1.0
    r2132 	 1
    r1956 	 -1.0
    EX_csa_LPAREN_e_RPAREN_ 	 3
    ATPtn 	 1
    DNDPt3m 	 -1.0
    S6T23g 	 -1.0
    OCTDECCACT 	 -1.0
    PCREATtmdiffir 	 1
    COAtp 	 1
    RE2659R 	 1
    CHOLESTle 	 1.0
    EX_HC02198_LPAREN_e_RPAREN_ 	 1
    DM_pnto_R 	 1
    RE2525X 	 1.0
    12HTACRhr 	 3.0
    ACACT9p 	 3.0
    6AHGLZABCt 	 1
    EX_rsv_LPAREN_e_RPAREN_ 	 1
    r2222 	 -1.0
    RE0453C 	 1
    PROSTGE1t3 	 -1.0
    r2465 	 1.0
    AG13T9g 	 1.0
    ACOAD1fm 	 2.0
    r0027 	 0.0
    RDH2a 	 3.0
    RE2637X 	 -1.0
    RE3445X 	 3.0
    NNATr 	 -1.0
    r1634 	 3.0
    P5CDm 	 -1.0
    RE0702N 	 -1.0
    PI45PLC 	 3.0
    ALOX12 	 3.0
    AGLPH 	 1
    G14T2g 	 1.0
    RE1525M 	 3.0
    MEPIVESSte 	 1
    IVCOAACBP 	 3.0
    PETHCT 	 -1.0
    GFUCS 	 -1.0
    1OHMDZitr 	 1
    EX_pchol_hs_LPAREN_e_RPAREN_ 	 1
    MM6B2g 	 3.0
    ILEATB0tc 	 1.0
    ATVLACtu 	 -1.0
    r1771 	 -1.0
    3ISPVSthc 	 -1.0
    DURIK1 	 2.0
    RE3010M 	 0
    FAOXC61x 	 1.0
    RE3260C 	 -1.0
    r1853 	 -1.0
    RE0575C 	 1
    GLUCYS 	 -1.0
    MI145PKn 	 -1.0
    3DPHBH1 	 1
    ARACHDtr 	 1
    EX_am1csa_LPAREN_e_RPAREN_ 	 1
    3MOBte 	 1
    RE3341X 	 3.0
    GLYC3Ptmc 	 1
    r1760 	 -1.0
    GALFUC12GAL14ACGLCGALGLUSIDEtg 	 1
    DOLPMT2_Uer 	 1
    DESAT24_1 	 -1.0
    HKt 	 -1.0
    RE3247X 	 3.0
    ADK1m 	 3.0
    AG13T3g 	 1.0
    EX_utp_LPAREN_e_RPAREN_ 	 3
    6OHFVSGLUhc 	 1
    TMDM3OATt 	 1
    r0395 	 -1.0
    ICDHyp 	 2.0
    HCO3_NAt 	 0.0
    RE3123C 	 1
    r0841 	 1
    S23T2g 	 -1.0
    r2292 	 -1.0
    r2312 	 -1.0
    EX_7bhglz_LPAREN_e_RPAREN_ 	 1
    RE3013R 	 -1.0
    NACHEX25ly 	 1.0
    RE3446R 	 2.0
    r2367 	 0
    SCPx 	 3.0
    FAOXC164C165x 	 0.0
    EX_4hatvlac_LPAREN_e_RPAREN_ 	 1
    GLYPROPEPT1tc 	 -1.0
    EX_acmp_LPAREN_e_RPAREN_ 	 1
    r2202 	 -1.0
    MEVK1x 	 3.0
    r0899 	 1
    SUCCt2m 	 -1.0
    FAOXC61_3Zm 	 1.0
    r2063 	 -1.0
    HYXNt 	 -1.0
    CBLATm 	 1.0
    ACMPGLUTitr 	 1
    UREAtm 	 -1.0
    FPGS4 	 -1.0
    DNDPt19m 	 -1.0
    RE3012R 	 3.0
    r0750 	 -1.0
    TOLBUTAMIDEte 	 1
    r1693 	 -1.0
    C142OHc 	 0
    BUTt2r 	 0.0
    RE1803C 	 1
    RE2908C 	 1
    r0074 	 -1.0
    r1990 	 -1.0
    EX_malt_LPAREN_e_RPAREN_ 	 3
    r2329 	 3.0
    r1993 	 -1.0
    NACHEXA3ly 	 -1.0
    FAH1 	 -1.0
    6HLVSTthep 	 1
    r1951 	 -1.0
    PTHPSn 	 2.0
    EX_malthx_LPAREN_e_RPAREN_ 	 1
    r1775 	 -1.0
    r0744 	 3.0
    RE3111M 	 -1.0
    GLCAASE4ly 	 3.0
    r1374 	 -1.0
    ENMAN5g 	 1
    r1739 	 -1.0
    r2016 	 -1.0
    FAOXC220200x 	 0.0
    EX_glc3meacp_LPAREN_e_RPAREN_ 	 1
    FAOXC142C122x 	 -1.0
    MMCDp 	 -1.0
    EX_tolbutamide_LPAREN_e_RPAREN_ 	 1
    r1969 	 -1.0
    RE0578X 	 -1.0
    RNMK 	 2.0
    7BHGLZtev 	 1
    r1796 	 -1.0
    r1949 	 -1.0
    RE1441G 	 -1.0
    SCP21cx 	 2.0
    r1603 	 3.0
    CYTK4 	 3.0
    RE0926E 	 0
    5ADTSTSTERONESULT 	 0
    CYSTSERex 	 -1.0
    SGPL11r 	 3.0
    r0051 	 3.0
    RE3122C 	 1
    r1837 	 -1.0
    5FTHFt2 	 -1.0
    24NPHte 	 1
    EX_am1acs_LPAREN_e_RPAREN_ 	 1
    EX_udp_LPAREN_e_RPAREN_ 	 1
    56DHPVStev 	 1
    PHCDm 	 -1.0
    ADHAPtx 	 1
    RE3193M 	 3.0
    UDPGLCtg 	 3.0
    EX_ribflv_LPAREN_e_RPAREN_ 	 3
    EX_tdechola_LPAREN_e_RPAREN_ 	 1
    P4502C9 	 -1.0
    r1298 	 1
    r1809 	 -1.0
    r2126 	 0.0
    TETPENT6CPT1 	 0
    INSt2 	 0.0
    CYSGLUexR 	 -1.0
    RE3443M 	 3.0
    RE3165C 	 3.0
    RE2514E 	 -1.0
    S6T10g 	 -1.0
    DNDPt49m 	 -1.0
    r0413 	 -1.0
    EX_retinol_cis_11_LPAREN_e_RPAREN_ 	 1
    AP4AH1 	 -1.0
    FUT18g 	 -1.0
    r0783 	 3.0
    AM1CSAtep 	 1
    r2201 	 -1.0
    RE2382C 	 1
    AM1ALCShr 	 1
    FAOXC226205m 	 1.0
    RE1653C 	 1
    EX_gtacmp_LPAREN_e_RPAREN_ 	 1
    r0993 	 -1.0
    PVSOATPtu 	 -1.0
    EX_4hdebrisoquine_LPAREN_e_RPAREN_ 	 1
    DHGLZtev 	 1
    RE0688E 	 -1.0
    PI34P5Kn 	 1
    ITCOALm 	 0.0
    P45021A1r 	 0
    EX_thmmp_LPAREN_e_RPAREN_ 	 1
    RN0027R 	 -1.0
    3HBCOARc 	 3.0
    SRTNMTX 	 -1.0
    RE3264C 	 -1.0
    THRASNNaEx 	 -1.0
    LVSTACIDhep 	 -1.0
    FUT11g 	 -1.0
    r2281 	 -1.0
    GCHOLAt2 	 -1.0
    r1937 	 -1.0
    RE3041C 	 3.0
    EX_glu_L_LPAREN_e_RPAREN_ 	 3
    r1885 	 -1.0
    S23Tg 	 -1.0
    G5SDym 	 1.0
    FAOXC141_5Em 	 3.0
    RE3522R 	 -1.0
    EX_andrstrnglc_LPAREN_e_RPAREN_ 	 1
    RE1308M 	 -1.0
    UGT1A5r 	 0
    TREHe 	 -1.0
    ST3GAL61g 	 -1.0
    RE2296X 	 1.0
    GTHPe 	 3.0
    FAOXC161140m 	 3.0
    EX_CE4633_LPAREN_e_RPAREN_ 	 1
    RE3511M 	 -1.0
    RE2910C 	 1
    EX_phe_L_LPAREN_e_RPAREN_ 	 3
    PROtm 	 1
    r2352 	 0
    FA181ACPH 	 -1.0
    r0179 	 1
    r0721 	 -1.0
    DESAT18_5 	 3.0
    NH4tp 	 1
    DGSNtm 	 1
    ACMPitr 	 1
    RE2428M 	 -1.0
    RE2704C 	 1
    S6T6g 	 -1.0
    FAOXC4020m 	 2.0
    r0839 	 1
    EX_lstnm4_LPAREN_e_RPAREN_ 	 1
    EX_co_LPAREN_e_RPAREN_ 	 1
    RE1064C 	 1
    r0834 	 -1.0
    RE2223M 	 1
    GALASE14ly 	 0
    r0754 	 2.0
    FAOXC16DCr 	 -1.0
    EX_whhdca_LPAREN_e_RPAREN_ 	 1
    4HGLSDm 	 -1.0
    EX_6melvacid_LPAREN_e_RPAREN_ 	 1
    EX_triodthysuf_LPAREN_e_RPAREN_ 	 1
    KCC2t 	 3.0
    PROt2rL 	 -1.0
    PI4P3Kn 	 1
    RE3446M 	 2.0
    DCK1n 	 -1.0
    EX_6ohfvsglu_LPAREN_e_RPAREN_ 	 1
    DOGULNO2 	 1
    HYPTROX 	 1
    EX_HC02207_LPAREN_e_RPAREN_ 	 1
    r0400 	 1
    RE0918C 	 1
    RE1927C 	 1
    RE2997X 	 3.0
    DKMPPD 	 1
    HYPTROXe 	 1
    RE3155C 	 1
    RE3111R 	 -1.0
    DNAMTn 	 2.0
    P45039A1r 	 -1.0
    GALt1r 	 3.0
    HISD 	 0.0
    r1753 	 -1.0
    SBPP3er 	 0.0
    EX_vacc_LPAREN_e_RPAREN_ 	 1
    COt 	 1
    FAOXC11090m 	 3.0
    EX_lac_L_LPAREN_e_RPAREN_ 	 3
    GTACMPhr 	 1
    EX_ump_LPAREN_e_RPAREN_ 	 1
    RE0944E 	 0
    S4T6g 	 -1.0
    H5MTer_U 	 1
    P4502D6 	 -1.0
    GALASE10ly 	 0
    MANter 	 1
    r2394 	 -1.0
    GLYOp 	 1
    EX_spc_hs_LPAREN_e_RPAREN_ 	 1
    EX_1513tacr_LPAREN_e_RPAREN_ 	 1
    DMNONCRNt 	 1
    EX_malcoa_LPAREN_e_RPAREN_ 	 1
    GLB1 	 0
    CSNAT2m 	 0.0
    RE3345C 	 -1.0
    r0410 	 3.0
    GALASE11ly 	 0
    SIAT9g 	 -1.0
    GLNCYSNaEx 	 -1.0
    FAOXC165C164x 	 -1.0
    r1043 	 0.0
    MM7B2g 	 3.0
    RE3352C 	 1
    1531TACRitr 	 1
    r0402 	 -1.0
    FUCFUCFUCGALACGLC13GALACGLCGAL14ACGLCGALGLUSIDEte 	 1
    CRVSM23tev 	 1
    AM19CShr 	 1
    r1847 	 -1.0
    FACOAL244_1 	 2.0
    EX_cdp_LPAREN_e_RPAREN_ 	 1
    P4502C93 	 -1.0
    LVSTACIDitr 	 1
    CSCPASEly 	 1
    EX_1331tacr_LPAREN_e_RPAREN_ 	 1
    RE0689X 	 -1.0
    RE2863C 	 1
    TYRTA 	 1.0
    MDZGLCitr 	 1
    RE3144C 	 1
    THYMDtl 	 -1.0
    r2154 	 -1.0
    C181CPT1 	 0
    r0407 	 3.0
    EX_thrfvs_LPAREN_e_RPAREN_ 	 1
    FAOXC14C14OHm 	 3.0
    PSYGCHe 	 1
    ACNGALACGLCGAL14ACGLCGALGLUSIDEte 	 1
    P4504F81r 	 -1.0
    DALAxt 	 1
    PIter 	 -1.0
    EX_pglyc_hs_LPAREN_e_RPAREN_ 	 3
    ESTRONEtr 	 1
    GRTT 	 -1.0
    RE3499C 	 1
    2HATVLACteb 	 -1.0
    RE3033R 	 1
    EX_htaxol_LPAREN_e_RPAREN_ 	 1
    CYTK7n 	 3.0
    RE3336X 	 3.0
    r1823 	 -1.0
    RE2989X 	 3.0
    RE0583C 	 1
    EX_4mptnl_LPAREN_e_RPAREN_ 	 1
    RE3155R 	 1
    6HSMVACIDteb 	 1
    RE0578C 	 0
    r1959 	 -1.0
    RE3475N 	 2.0
    RE3248X 	 -1.0
    SERASNNaEx 	 -1.0
    3HPVSCOAitx 	 1
    NTD5 	 3.0
    r2473 	 -1.0
    DNDPt59m 	 -1.0
    ACACt2 	 1.0
    DM_13_cis_oretn_n_ 	 1
    CSAtd 	 1
    MLTG1ly 	 -1.0
    C161CRNt 	 -1.0
    NACDe 	 1
    FOLOATPtc 	 -1.0
    r1377 	 1
    r0556 	 0
    AM1C9CSteb 	 1
    RE2133C 	 1
    FTHFL 	 2.0
    EX_c5dc_ 	 1
    TCHOLAt3 	 1.0
    FAOXC181C161x 	 -1.0
    C80CPT1 	 0
    OCTAt 	 1
    PRGNLONESULT 	 3.0
    DMNONCOACRNCPT1 	 1
    BUTt2m 	 0.0
    LGNCCOAtcx 	 -1.0
    r1978 	 -1.0
    r2203 	 -1.0
    CYTK1m 	 -1.0
    FAOXC2242046x 	 0.0
    ACtg 	 1
    DOPACHRMISO 	 -1.0
    3MOBt2im 	 1
    RE3525M 	 1
    SERPT 	 -1.0
    r2028 	 -1.0
    RE3494C 	 1
    RE3268R 	 1
    3AIBTm 	 1
    r0706 	 3.0
    r1789 	 -1.0
    EX_6bhglzglc_LPAREN_e_RPAREN_ 	 1
    r1011 	 1
    FACOAL100i 	 1
    FUC13GALACGLCGAL14ACGLCGALGLUSIDEte 	 1
    C4OHc 	 0
    RE3038C 	 0.0
    NMNATm 	 -1.0
    EX_dgmp_LPAREN_e_RPAREN_ 	 1
    RE1632R 	 3.0
    FAOXC123x 	 1.0
    RE3561X 	 3.0
    C100CPT1 	 0
    r2495 	 -1.0
    RE1530C 	 1
    TDCHOLAtx 	 1
    NTD2e 	 0.0
    RE0827X 	 -1.0
    RE1050L 	 3.0
    RE0927C 	 -1.0
    BTNPL 	 -1.0
    r2266 	 -1.0
    r1172 	 1.0
    LEULEULAPc 	 -1.0
    r2290 	 -1.0
    INSTt2r 	 -1.0
    HACD1x 	 3.0
    EX_acald_LPAREN_e_RPAREN_ 	 3
    EX_gluala_LPAREN_e_RPAREN_ 	 3
    RE2856C 	 1
    AGLPED 	 1
    PI345P5Pn 	 1
    r2235 	 -1.0
    TETTET6CRNt 	 -1.0
    DNDPt29m 	 1
    r1174 	 1.0
    FATP2t 	 -1.0
    RE2655R 	 1
    FAOXTC162TC142m 	 3.0
    EX_3bcrn_ 	 1
    PROFVStev 	 1
    RE3038N 	 1
    DESAT20_2 	 0
    XYLK 	 -1.0
    AOBUTDsm 	 1
    FAOXC5C5OHm 	 1
    C3STDH1Pr 	 3.0
    PYLALDOXm 	 3.0
    CSASULPteb 	 1
    PECGONCOATr 	 1
    r1466 	 1
    EX_CLPND_LPAREN_e_RPAREN_ 	 3
    FAOXC101C102m 	 2.0
    DOCOSACTDe 	 1
    FAOXC8DCC6DCx 	 -1.0
    ACACT1x 	 -1.0
    RE0456M 	 1
    FAOXC15NADx 	 -1.0
    RE3113R 	 1
    3OHACMPtev 	 1
    C2tcx 	 -1.0
    EX_ha_pre1_LPAREN_e_RPAREN_ 	 1
    r0525 	 0.0
    r0774 	 3.0
    r0160 	 -1.0
    FPGS7 	 -1.0
    FUT910g 	 -1.0
    HSD3B7P 	 -1.0
    r0010 	 3.0
    RE1808R 	 1
    MCD 	 -1.0
    r0769 	 -1.0
    Ser_Thrtg 	 1
    RE1635M 	 1.0
    r2494 	 -1.0
    r2053 	 -1.0
    EX_am9csa_LPAREN_e_RPAREN_ 	 3
    RE1898C 	 1
    RE0576C 	 1
    EX_srtn_LPAREN_e_RPAREN_ 	 1
    RE0920R 	 0
    NDPK5n 	 0
    HEDCECRNe 	 1
    RE2622C 	 -1.0
    2HATVLACGLUChr 	 0
    URIK1 	 0
    r0643 	 3.0
    RE1816X 	 1.0
    ADPRDP 	 1.0
    r1919 	 -1.0
    EX_tyr_L_LPAREN_e_RPAREN_ 	 3
    r0718 	 3.0
    3HPVSCOAhc 	 1
    r2507 	 1
    EX_4hmdgluc_LPAREN_e_RPAREN_ 	 1
    AKGMALtm 	 -1.0
    CYSATB0tc 	 1.0
    FATP7t 	 -1.0
    LSTN1GLUCitr 	 1
    EX_rsvlac_LPAREN_e_RPAREN_ 	 1
    4HATVACIDitr 	 1
    r2097 	 0.0
    CYSPHELAT2tc 	 1
    GLNtN1 	 -1.0
    RE0577M 	 1.0
    r1638 	 3.0
    MEOHtly 	 1
    SERATB0tc 	 1.0
    CRVSM22itr 	 1
    PHETA1 	 1.0
    EX_maltpt_LPAREN_e_RPAREN_ 	 1
    FAOXC2251836x 	 0.0
    r0795 	 1
    GLCNACDASg 	 -1.0
    THCHOLSTOICtm 	 1
    r1985 	 -1.0
    DNDPt42m 	 -1.0
    r0548 	 2.0
    r0845 	 -1.0
    r0309 	 -1.0
    r2416 	 -1.0
    r0295 	 -1.0
    CITL 	 1
    RETI3 	 1
    r1814 	 -1.0
    FUT35g 	 1.0
    ALAGLNexR 	 -1.0
    KSII_CORE2tly 	 1
    FACOAL184 	 2.0
    EX_5fthf_LPAREN_e_RPAREN_ 	 1
    UDPDOLPT_U 	 1
    RTOTALFATPc 	 1
    2HATVACIDteb 	 -1.0
    NDPK2n 	 0
    ADSK 	 3.0
    DASCBR 	 3.0
    r0358 	 3.0
    CERT1rt 	 0.0
    r0885 	 -1.0
    THMDt2r 	 0.0
    RE2442C 	 1
    S6T16g 	 -1.0
    FAOXC201C181x 	 -1.0
    GLYPHEPEPT1tc 	 -1.0
    HSD11B1r 	 -1.0
    RE3040X 	 0.0
    CYTK1 	 3.0
    NADHtpu 	 1
    LSTN1GLUChr 	 0
    RE3580X 	 3.0
    DGULND 	 1
    RE0702L 	 3.0
    ACNAMPH 	 0.0
    TMDK1m 	 -1.0
    DOPAt4_2_r 	 -1.0
    LCAT1e 	 -1.0
    RE3097X 	 3.0
    ATVACIDtu 	 -1.0
    ADPCOACROT 	 3.0
    RE1956C 	 0.0
    7BHGLZGLCABCt 	 1
    EX_leuleu_LPAREN_e_RPAREN_ 	 1
    11DOCRTSTRNtm 	 1
    GLYCLTDym 	 1
    S6TASE14ly 	 3.0
    LCYSTCBOXL 	 0.0
    FAOXC162_7Z_10Zm 	 3.0
    VLCSr 	 -1.0
    ARGATB0tc 	 1.0
    FUMTSULtm 	 -1.0
    r1062 	 1
    CYTK10 	 3.0
    CRVSM31itr 	 1
    r0113 	 -1.0
    EX_fuc_L_LPAREN_e_RPAREN_ 	 3
    PROPAT4te 	 -1.0
    FAOXC241221x 	 0.0
    r0083 	 1.0
    DNDPt57m 	 -1.0
    DM_avite1_c_ 	 1
    DM_anth 	 1
    MI1345PP 	 1.0
    r0033 	 -1.0
    r0661 	 -1.0
    DOGULND2 	 1
    GCHOLAte 	 1
    EX_caro_LPAREN_e_RPAREN_ 	 3
    RE2027C 	 1
    ACMPGLUtep 	 -1.0
    DESAT18_3 	 3.0
    FUM 	 3.0
    ALR3 	 3.0
    FAOXC101_3Em 	 1.0
    TRIPVShc 	 3.0
    LSTNM5itr 	 1
    RE2383C 	 1
    PCHOLP_hs 	 3.0
    EX_c81crn_ 	 1
    AG13T10g 	 1.0
    CO2tn 	 1
    RE1952X 	 0.0
    EX_fru_LPAREN_e_RPAREN_ 	 3
    LEUKABCtc 	 1.0
    EX_creat_LPAREN_e_RPAREN_ 	 3
    RE2273C 	 1
    RE0958C 	 1
    r2147 	 -1.0
    EX_adprbp_LPAREN_e_RPAREN_ 	 1
    FAOXC15NADPx 	 -1.0
    r2083 	 0.0
    r0758 	 3.0
    UGT1A8r 	 0
    PPAtm 	 1
    ASNNm 	 -1.0
    r1596 	 3.0
    RE0915C 	 1
    AHCYSte 	 1
    FBA5 	 -1.0
    r1874 	 -1.0
    MTHFCm 	 3.0
    EX_vitd3_LPAREN_e_RPAREN_ 	 1
    DESAT18_9 	 0
    PPAm 	 3.0
    RDH3 	 -1.0
    C226CPT1 	 0
    TETPENT3t 	 1
    r0727 	 -1.0
    FA1821ACPH 	 1
    FAOXC205_5Z_8Z_11Z_14Z_17Zx 	 0.0
    CRVS1M24hc 	 3.0
    TETHEX3COAtx 	 1
    r0377 	 -1.0
    r1391 	 1.0
    RE3514C 	 1
    34DHPHELAT1tc 	 1
    PMI12346PH 	 -1.0
    1a_25VITD3Hm 	 1
    r0617 	 0
    EX_ptvstlac_LPAREN_e_RPAREN_ 	 3
    PCLYSOX 	 3.0
    3HPPD 	 1
    RE0580R 	 1
    THMDt5le 	 0.0
    GLCURtly 	 -1.0
    RE1796R 	 -1.0
    AFLATOXINte 	 1
    GGH_7THFe 	 1.0
    RE0828E 	 -1.0
    EX_carn_LPAREN_e_RPAREN_ 	 3
    r1698 	 -1.0
    r0723 	 3.0
    RE3310C 	 1
    DNDPt22m 	 1
    r0789 	 1
    r1888 	 -1.0
    4MOPt2im 	 1
    r2262 	 -1.0
    IBUPGLUCitr 	 1
    RE0453M 	 -1.0
    F1PGT 	 -1.0
    C204CPT2 	 -1.0
    MM6bg 	 3.0
    APRTO2 	 3.0
    r2181 	 -1.0
    r2489 	 -1.0
    PLA2_2e 	 0
    AMACR2r 	 -1.0
    r1658 	 3.0
    EX_4mtolbutamide_LPAREN_e_RPAREN_ 	 1
    RE2913M 	 3.0
    r0392 	 3.0
    EX_ala_B_LPAREN_e_RPAREN_ 	 3
    TCHOLAtx 	 1
    r1073 	 1
    r0311 	 2.0
    r2308 	 -1.0
    UGT1A7r 	 0
    C3STKR2r 	 3.0
    CREATt4_2_r 	 2.0
    DM_Asn_X_Ser_Thr_ly_ 	 1
    CYTK8n 	 3.0
    B3GNT37g 	 -1.0
    RE3103C 	 3.0
    PMEVKc 	 -1.0
    ST6GALNAC31 	 -1.0
    HGNTOR 	 -1.0
    EX_CE2839_LPAREN_e_RPAREN_ 	 1
    RE3521R 	 3.0
    RE3572X 	 3.0
    DGCHOLte 	 1
    CARIBUP_SGLUhep 	 0
    LINKDEG4ly 	 1
    EX_abt_LPAREN_e_RPAREN_ 	 1
    EX_1mncam_LPAREN_e_RPAREN_ 	 1
    r1177 	 1.0
    GTHO 	 3.0
    RETNGLCt2r 	 1
    r0568 	 3.0
    RE2913X 	 3.0
    56EPPVStev 	 1
    EX_chsterol_LPAREN_e_RPAREN_ 	 3
    RE3458C 	 1
    GLYPROPRO1c 	 -1.0
    RE2915M 	 3.0
    r1998 	 -1.0
    r1705 	 -1.0
    PCHOLABCtc 	 3.0
    RE3247M 	 0.0
    r1444 	 3.0
    FUCASE2e 	 -1.0
    RE0581R 	 0
    RE3167R 	 1
    AG13T5g 	 1.0
    FPGS2 	 -1.0
    r1579 	 3.0
    RE2666C 	 1
    DADNK 	 1
    r1690 	 -1.0
    r0465 	 -1.0
    RE3519X 	 0.0
    FVSGLUChc 	 1
    HISTAVESSEC 	 0.0
    ARTFR207 	 1
    GUAt 	 -1.0
    r2170 	 -1.0
    r0697 	 -1.0
    EX_anth_LPAREN_e_RPAREN_ 	 1
    LDH_L 	 3.0
    FAOXC225C204x 	 -1.0
    1OHMDZhr 	 3.0
    3HPVSTETCOAitx 	 1
    ESTROSABCCte 	 1.0
    COQ6m 	 -1.0
    RN0030C 	 1
    EX_fucfuc12gal14acglcgalgluside_hs_LPAREN_e_RPAREN_ 	 1
    ACCOALm 	 -1.0
    r0423 	 1.0
    EX_sel_LPAREN_e_RPAREN_ 	 3
    C4tcx 	 -1.0
    FAOXC2031836m 	 3.0
    r1756 	 -1.0
    LTD4DP 	 1
    DOCOSADIACTD 	 1
    S6T20g 	 1.0
    RE3259R 	 0.0
    r1972 	 -1.0
    C40CPT1 	 0
    ACACT7p 	 3.0
    RSVitr 	 1
    RE3476X 	 0
    CERK 	 -1.0
    DTTPtm 	 1
    RE3524R 	 1
    RADH3 	 1
    10FTHFtl 	 1
    CYSTA 	 1.0
    B3GNT311g 	 -1.0
    ST6GALNAC23 	 3.0
    2AMADPTm 	 -1.0
    r2512 	 1
    SUBERCACT 	 -1.0
    GALGT1 	 -1.0
    EX_3octdeccrn_ 	 1
    RE3345X 	 2.0
    ACCOACm 	 -1.0
    P4501B1r 	 1
    PGI 	 1.0
    RE1709N 	 -1.0
    GGH_6DHFl 	 1.0
    EX_lstnm5_LPAREN_e_RPAREN_ 	 1
    RE3126C 	 1
    r1859 	 -1.0
    RE3171C 	 1
    RE1518M 	 -1.0
    O2St 	 1
    RE3040C 	 0.0
    2HATVLACGLUCteb 	 -1.0
    EX_tag_hs_LPAREN_e_RPAREN_ 	 3
    EBASTINEte 	 1
    TRYPTAOX 	 3.0
    DNDPt32m 	 -1.0
    DHFR 	 -1.0
    CITt4_2 	 0.0
    2HCO3_NAt 	 1.0
    EX_acmpglut_LPAREN_e_RPAREN_ 	 1
    EX_7ahglz_LPAREN_e_RPAREN_ 	 1
    PHEACGLNt 	 1
    PPM 	 3.0
    FRUt4 	 -1.0
    AM4NCShr 	 1
    3HLVSTACtbc 	 1
    CLHCO3tex2 	 -1.0
    FTHFCL 	 0
    COQ7m 	 -1.0
    EX_6csmvacid_LPAREN_e_RPAREN_ 	 1
    24_25DHVITD3t 	 1
    RE2079R 	 1
    EX_ahandrostanglc_LPAREN_e_RPAREN_ 	 1
    RE2986X 	 -1.0
    EX_dcyt_LPAREN_e_RPAREN_ 	 3
    EX_3hexdcrn_ 	 1
    BHBtm 	 1
    RE3551X 	 1
    r1598 	 3.0
    ALAt2r 	 0
    r0992 	 -1.0
    O16G2e 	 -1.0
    S6TASE22ly 	 -1.0
    RE2112C 	 1
    ISOLVSTtbc 	 -1.0
    RE2878C 	 1
    EX_nrpphrsf_LPAREN_e_RPAREN_ 	 1
    LSTNM1itr 	 1
    RE2854C 	 1
    TRPATB0tc 	 1.0
    EX_35cgmp_LPAREN_e_RPAREN_ 	 1
    FUT12g 	 -1.0
    HOCDACBP 	 3.0
    FAOXC121C10x 	 -1.0
    FVSteb 	 1
    EX_ocdca_LPAREN_e_RPAREN_ 	 3
    NDPK4m 	 -1.0
    GABAVESSEC 	 -1.0
    r2161 	 -1.0
    RE1701C 	 1
    r1928 	 -1.0
    25VITD2Hm 	 0.0
    RE3227C 	 1
    r2356 	 0
    r0470 	 -1.0
    MEOHt2 	 1
    RE1817X 	 1.0
    RE0689E 	 -1.0
    EX_orot_LPAREN_e_RPAREN_ 	 1
    r1321 	 1
    VALTAm 	 -1.0
    LYSMTF2n 	 3.0
    FUCFUCFUCGALACGLCGAL14ACGLCGALGLUSIDEtg 	 1
    EX_retinol_LPAREN_e_RPAREN_ 	 3
    TRIODTHYt2 	 -1.0
    EX_yvite_LPAREN_e_RPAREN_ 	 1
    RE3148C 	 3.0
    ENMAN1g 	 1
    NAGAlby 	 1
    RE3186M 	 2.0
    MMMm 	 3.0
    TRPt4 	 1.0
    FAOXC143C123x 	 -1.0
    r1010 	 1
    EX_gumtchol_LPAREN_e_RPAREN_ 	 1
    EX_gumgchol_LPAREN_e_RPAREN_ 	 1
    6CSMVitr 	 1
    EX_apnnox_LPAREN_e_RPAREN_ 	 1
    EX_retnglc_LPAREN_e_RPAREN_ 	 1
    r1679 	 -1.0
    RE3565C 	 1
    r1692 	 -1.0
    SMVtv 	 1
    r1315 	 -1.0
    r0615 	 -1.0
    S6TASE3ly 	 3.0
    ACMPtu 	 -1.0
    CHOLPtg 	 1
    N4Tg 	 1
    LPS4e 	 2.0
    NNMT 	 1.0
    TRIPVStev 	 1
    2DR1PP 	 1
    FAOXC185_3Z_6Z_9Z_12Z_15Zm 	 1.0
    S6T2g 	 -1.0
    DOLASNT_Uer 	 1
    r2114 	 1
    RE3126R 	 1
    GALNTg 	 3.0
    r0946 	 -1.0
    r1810 	 -1.0
    FAOXC16DCC14DCx 	 -1.0
    r1166 	 -1.0
    r2321 	 3.0
    r2274 	 -1.0
    r1135 	 2.0
    RE3581X 	 -1.0
    DNDPt47m 	 -1.0
    EX_dgsn_LPAREN_e_RPAREN_ 	 3
    RE1099R 	 -1.0
    r2307 	 -1.0
    VALPHELAT2tc 	 1
    P4502C19 	 -1.0
    TCHOLAte 	 1
    RE3554C 	 1
    S6T14g 	 -1.0
    GLCAT9g 	 -1.0
    EX_2hatvacid_LPAREN_e_RPAREN_ 	 1
    ACOAD8m 	 -1.0
    r2204 	 -1.0
    RE2272L 	 3.0
    1a_24_25VITD2Hm 	 1
    GALACGLCGALGBSIDEte 	 1
    RE2124C 	 1
    RE3335M 	 2.0
    FORtrn 	 1
    ATVLACtdhc 	 1
    NDPK2m 	 -1.0
    RE1534X 	 3.0
    ORNt3m 	 -1.0
    RE2766C 	 -1.0
    GBAl 	 0
    BHMT 	 -1.0
    RE3519R 	 0.0
    BUTSMCT1 	 -1.0
    r2249 	 -1.0
    1HMDZGLUChc 	 1
    r0249 	 1
    NDPK2 	 0
    FUT31g 	 1.0
    r1008 	 1
    ACGBGBSIDEtg 	 1
    GGT_U 	 1
    EX_CE0074_LPAREN_e_RPAREN_ 	 1
    MAN2_7Cer 	 1
    HOCTDECCRNe 	 1
    SBTR 	 3.0
    r0995 	 -1.0
    DOPASFt 	 1
    RE1308C 	 1
    TYRPHELAT2tc 	 1
    C3STDH1r 	 3.0
    RE2899C 	 1
    RE2270E 	 0.0
    UAGALDP 	 1
    FUCASEly 	 -1.0
    MAOX 	 0.0
    EX_CE2011_LPAREN_e_RPAREN_ 	 1
    7DHFtl 	 1
    ILELAT1tc 	 1
    RE3265R 	 0
    LST4EXPhr 	 3.0
    BILDGLCURte 	 1.0
    LGNCFATtc 	 3.0
    RE2677C 	 1
    C161OHc 	 0
    r2096 	 0.0
    ALCD22_L 	 3.0
    DOLPGT2_Ler 	 0.0
    7HPVStev 	 1
    ECOAH1m 	 3.0
    r1327 	 1
    MM7Ag 	 0.0
    S6T17g 	 -1.0
    TRPB0AT3tc 	 -1.0
    FAOXC61C4m 	 0.0
    EX_3ivcrn_ 	 1
    EX_k_LPAREN_e_RPAREN_ 	 3
    ACNACNGAL14ACGLCGALGLUSIDEte 	 1
    RE1100C 	 -1.0
    RE3451C 	 1
    FAOXC15BRC13BRx 	 -1.0
    LYSATB0tc 	 1.0
    RE3295C 	 1
    r0398 	 0.0
    r0773 	 -1.0
    DM_gncore2_g_ 	 1
    PETOHMm_hs 	 -1.0
    C10DCe 	 1
    FAOXC143_5Z_8Z_11Zx 	 0.0
    r1803 	 -1.0
    RE3378C 	 1
    ESTRADIOLGLCt 	 1.0
    GAO1 	 1.0
    r2185 	 -1.0
    VITKtl 	 1
    PPASMCT1 	 -1.0
    r2354 	 0
    35DHPVShc 	 3.0
    GLYC3PFADm 	 1
    MHGLZtev 	 1
    r2487 	 -1.0
    EX_but_LPAREN_e_RPAREN_ 	 1
    r2272 	 -1.0
    EX_pydam_LPAREN_e_RPAREN_ 	 3
    EX_hdcea_LPAREN_e_RPAREN_ 	 3
    EX_arab_L_LPAREN_e_RPAREN_ 	 3
    r1759 	 -1.0
    PGPPT 	 -1.0
    r1954 	 -1.0
    r0713 	 -1.0
    SO4tl 	 1
    HYPOE 	 0
    EX_ddece1crn_ 	 1
    PLA2 	 1
    RE1311C 	 1
    PTVSThc 	 -1.0
    FAOXC101_4Em 	 1.0
    MI1346Ptn 	 1
    PYDXPP 	 0
    r1106 	 1
    PE_HSter 	 1
    FUT97g 	 -1.0
    EX_fucfucfucgalacglcgal14acglcgalgluside_hs_LPAREN_e_RPAREN_ 	 1
    GGH_10FTHF5GLUe 	 1.0
    P45017A1r 	 -1.0
    HSPASEly 	 1
    EX_C02528_LPAREN_e_RPAREN_ 	 1
    FAOXC8C6x 	 -1.0
    RE2050R 	 2.0
    FVSTETtev 	 1
    UROLACer 	 1
    7AHGLZhr 	 1
    SPHGNtr 	 1
    EX_glyc_LPAREN_e_RPAREN_ 	 3
    r1882 	 -1.0
    RE0921C 	 -1.0
    OBDHc 	 1
    INStm 	 -1.0
    r2405 	 -1.0
    r0281 	 -1.0
    r0494 	 3.0
    34DHPLACOX 	 2.0
    PNTKm 	 3.0
    r1988 	 -1.0
    PRODt2rL 	 -1.0
    RE1943C 	 1
    RE0830N 	 -1.0
    HOCTDEC2CRNe 	 1
    ECOAH1x 	 2.0
    SRTN23OX 	 -1.0
    r1747 	 -1.0
    CYSLYSL 	 1
    r2360 	 0
    DOGULND1 	 1
    TRPO2 	 -1.0
    r0821 	 1
    HPDCACRNCPT2 	 -1.0
    r2257 	 -1.0
    r2282 	 -1.0
    NDPK1 	 0
    PHACCOAGLNAC 	 1
    RE3259C 	 0.0
    RE3269C 	 0.0
    NP1 	 2.0
    EX_am1c9cs_LPAREN_e_RPAREN_ 	 1
    NACHEX2ly 	 1.0
    DNDPt38m 	 -1.0
    RE3258C 	 0.0
    ATP2ter 	 1
    r0321 	 1
    ACCOAtn 	 1
    VITEtl 	 3.0
    MI13PP 	 1
    PVStep 	 1
    GHMT3m 	 3.0
    r0639 	 3.0
    PI3P4Kn 	 1
    G13MT_U 	 1
    PIK3er 	 3.0
    RE3448X 	 3.0
    ABUTt2r 	 -1.0
    CRNCAR3tp 	 1
    EX_cytd_LPAREN_e_RPAREN_ 	 3
    C81CRNe 	 1
    r0770 	 -1.0
    7AHGLZtev 	 1
    OCDEAFABP1tc 	 -1.0
    r1780 	 -1.0
    r0226 	 1
    RE3050R 	 1
    RE2974R 	 -1.0
    MDZGLCtev 	 1
    DNDPt63m 	 -1.0
    NDPK7n 	 0
    6AHGLZitr 	 1
    FAOXC123C102x 	 -1.0
    RE0702E 	 -1.0
    r1472 	 -1.0
    LEUt5m 	 1
    r1943 	 -1.0
    UDPG4E 	 0.0
    MG3er 	 1.0
    FAOXC2251836m 	 2.0
    r1423 	 1
    EX_tlacfvs_LPAREN_e_RPAREN_ 	 1
    GALNACT5g 	 0.0
    FE3MTP1 	 3.0
    EX_leuktrA4_LPAREN_e_RPAREN_ 	 3
    RE0453N 	 2.0
    FCLTm 	 -1.0
    NACHEXA18ly 	 -1.0
    r1597 	 3.0
    C162OHe 	 1
    FAOXC102m 	 1.0
    5HXKYNDCL 	 -1.0
    4BHGLZhr 	 1
    r1675 	 -1.0
    LYSMTF1n 	 3.0
    TYRTAm 	 3.0
    XOL7AH2tm 	 1
    r1941 	 -1.0
    CSAPASEly 	 1
    FAOXC163_4Z_7Z_10Zx 	 -1.0
    DESFVShc 	 -1.0
    H7_ETer 	 -1.0
    THBPT4ACAMDASE 	 -1.0
    EX_retn_LPAREN_e_RPAREN_ 	 3
    r0668 	 -1.0
    RE0066M 	 -1.0
    r2042 	 -1.0
    RE3566C 	 1
    BALABETAtc 	 0.0
    r1898 	 -1.0
    ASCBt 	 1
    r0636 	 1
    RE3440C 	 1
    EX_smvacid_LPAREN_e_RPAREN_ 	 1
    r1632 	 -1.0
    UGCG 	 3.0
    TLACFVSitr 	 1
    FRTT 	 -1.0
    3HPVSTEThc 	 1
    r0931 	 1
    RE3139X 	 3.0
    r1255 	 2.0
    FAOXC123m 	 1.0
    SPMTDe 	 1
    RPE 	 -1.0
    r1848 	 -1.0
    r0308 	 3.0
    SPH1Pte 	 1
    F1Atg 	 1
    EX_leuktrF4_LPAREN_e_RPAREN_ 	 1
    RE0926C 	 1
    HXANtx 	 1
    RE3520E 	 1
    r1680 	 -1.0
    EX_tsacmsul_LPAREN_e_RPAREN_ 	 1
    TYMSULT 	 1.0
    TRIOK 	 1
    DNDPt18m 	 -1.0
    GLUTCOAACBP 	 3.0
    HSD3B12r 	 -1.0
    C4STMO2r 	 3.0
    r2537 	 1
    ENGASE2ly 	 1
    r1580 	 3.0
    PHEyLATthc 	 -1.0
    RE3513C 	 1
    FATP5t 	 -1.0
    DNDPt56m 	 -1.0
    EX_CE2250_LPAREN_e_RPAREN_ 	 1
    FAOXC246226x 	 0.0
    RE1537X 	 -1.0
    FACOAL181i 	 2.0
    TSTSTERONESULT 	 0
    PROFVSCOAitx 	 1
    r0645 	 1
    r0511 	 3.0
    GLYCK2 	 1
    EX_tmdm3_LPAREN_e_RPAREN_ 	 1
    NDPK7 	 0
    PI45PLCn 	 -1.0
    EX_crtstrn_LPAREN_e_RPAREN_ 	 1
    56DHPVSteb 	 0.0
    RE0923R 	 0
    UMPK2 	 3.0
    EX_3hsmvacid_LPAREN_e_RPAREN_ 	 1
    r0673 	 -1.0
    r1942 	 -1.0
    r0647 	 1
    PNTK 	 3.0
    EX_caribup_R_LPAREN_e_RPAREN_ 	 1
    SPODM 	 3.0
    RE2872C 	 1
    MTHFD2 	 3.0
    RE1830C 	 1
    r2037 	 -1.0
    RE3233C 	 1
    EX_4abut_LPAREN_e_RPAREN_ 	 1
    Uritn 	 1
    THMt2m 	 1
    SFCYSe 	 1
    AKR1D2 	 -1.0
    EX_arach_LPAREN_e_RPAREN_ 	 1
    RE3511C 	 1
    EX_akg_LPAREN_e_RPAREN_ 	 3
    LIPOti 	 -1.0
    B3GALT5g 	 -1.0
    SPHINGStr 	 1
    DRIBt 	 1
    ALCD1 	 3.0
    r2496 	 -1.0
    ETHP 	 -1.0
    FAOXTC182TC162m 	 3.0
    5OHFVSGLUtev 	 1
    5MTHFt2le 	 1
    r1007 	 1
    r1910 	 -1.0
    O2t 	 1
    RE3402M 	 3.0
    TRDR2 	 0.0
    6MELVACIDitr 	 1
    EX_HC02179_LPAREN_e_RPAREN_ 	 1
    r2211 	 -1.0
    ADPACTD 	 1
    ACHVESSEC 	 -1.0
    EX_cysam_LPAREN_e_RPAREN_ 	 1
    FAOXC201181x 	 0.0
    GALNACT4g 	 0.0
    ICDHyrm 	 1.0
    UMPK5n 	 3.0
    ARTPLM3m 	 1
    r0927 	 1
    EX_HC02180_LPAREN_e_RPAREN_ 	 1
    RN0021X 	 0.0
    r0205 	 1.0
    FORt2m 	 1
    r2334 	 3.0
    UGALGTg 	 -1.0
    PEPCK 	 -1.0
    CRNtuNa 	 -1.0
    FAOXC22OHC22r 	 -1.0
    RN0028C 	 0.0
    RE1817R 	 3.0
    r0441 	 1
    GNMT 	 -1.0
    GLUTCOADHm 	 -1.0
    NABTNO 	 3.0
    EX_am1ccs_LPAREN_e_RPAREN_ 	 1
    FACOAL140i 	 1.0
    SMVHYDROhep 	 -1.0
    r0391 	 1.0
    RE3488X 	 0.0
    RE1240C 	 1
    ALLOP2tu 	 0.0
    r0170 	 3.0
    ACTNMO 	 -1.0
    FACOAL1831 	 2.0
    THSACMPhr 	 1
    ABO2g 	 -1.0
    RE3073X 	 1
    RE3162C 	 1
    RE1135L 	 -1.0
    r1059 	 -1.0
    r0339 	 -1.0
    RE1944C 	 1
    ST8SIA11 	 -1.0
    RE3401M 	 3.0
    RE1517X 	 1
    EX_4pyrdx_LPAREN_e_RPAREN_ 	 1
    CORE7GTg 	 1
    RE3437C 	 1
    EX_gum_LPAREN_e_RPAREN_ 	 1
    r2055 	 -1.0
    r1594 	 3.0
    GLUt2m 	 0.0
    P45017A3r 	 -1.0
    r0483 	 -1.0
    PAFH 	 -1.0
    HTAXOLte 	 1
    RE3396M 	 3.0
    RE3417C 	 1
    FUT95g 	 -1.0
    EICOSTETCPT1 	 0
    DDPGAm 	 1
    CRVSM31hc 	 1
    FAOXC9070m 	 3.0
    CTPS2 	 0.0
    XYLTt 	 1
    r1833 	 -1.0
    CARIBUP_Rthv 	 1
    RE2919X 	 -1.0
    P45027A15m 	 -1.0
    r1659 	 3.0
    3OHACMPitr 	 1
    DNDPt1m 	 -1.0
    PHETA1m 	 3.0
    r1960 	 -1.0
    SRTNt6_2_r 	 0.0
    TXA2tr 	 1
    HSD3A2r 	 -1.0
    ABO9g 	 -1.0
    XYLTer 	 -1.0
    r1635 	 3.0
    RE3088X 	 3.0
    LACLt 	 1
    r1733 	 -1.0
    EX_6hsmvacid_LPAREN_e_RPAREN_ 	 1
    r1326 	 1
    24_25DHVITD2tm 	 1
    THRFVStev 	 1
    NADHtru 	 1
    THYMDtm 	 -1.0
    RE3170C 	 3.0
    EX_ahcys_LPAREN_e_RPAREN_ 	 1
    NADPHtxu 	 1
    EX_galgalfucfucgalacglcgalacglcgal14acglcgalgluside_hs_LPAREN_e_RPAREN_ 	 1
    XANtx 	 1
    EX_aicar_LPAREN_e_RPAREN_ 	 1
    EBP2r 	 2.0
    WHTSTSTERONEte 	 1
    EX_adprib_LPAREN_e_RPAREN_ 	 1
    RE3436R 	 2.0
    FAOXC226C225x 	 -1.0
    FAOXC101_4Zx 	 -1.0
    SEBCOACROT 	 3.0
    24_25DHVITD3tm 	 1
    RE3002X 	 1.0
    GLCAASE8ly 	 3.0
    FAOXC164C143x 	 -1.0
    CSNATr 	 0.0
    ACt2r 	 1
    APOC_LYS_BTNPm 	 1
    EX_tsul_LPAREN_e_RPAREN_ 	 1
    FAOXC141C121m 	 3.0
    AM1C9CShr 	 1
    HS2ly 	 1.0
    FAOXC13BRC11BRx 	 -1.0
    OXYPtepv 	 1
    RE3562M 	 2.0
    UMPK7n 	 3.0
    METATB0tc 	 1.0
    r2038 	 -1.0
    C226COAtx 	 1
    DHPR2 	 1
    STRDNCCPT1 	 0
    G14T3g 	 1.0
    r0330 	 1
    TAUBETAtc 	 0.0
    FAS180COA 	 -1.0
    ACOX2x 	 -1.0
    GGH_10FTHF6GLUe 	 1.0
    S4TASE5ly 	 -1.0
    56EPPVSitr 	 1
    RE2269C 	 1
    RE2563C 	 1
    RE3415C 	 1
    RE3173R 	 1
    RE0925R 	 0
    GULNDer 	 1
    FAOXC7C5m 	 0.0
    RSVLACteb 	 0.0
    EX_ca2_LPAREN_e_RPAREN_ 	 3
    P4502E1 	 -1.0
    HDCAFAPMtc 	 3.0
    r0730 	 3.0
    CAMPt 	 1.0
    DNDPt2m 	 -1.0
    r1938 	 -1.0
    S23T4g 	 -1.0
    RE1514M 	 2.0
    LEUGLYHYc 	 1
    VALLAT1tc 	 1
    PGPP_hs 	 -1.0
    SERLYSNaex 	 -1.0
    r2118 	 0.0
    WHDDCAte 	 1
    H2Ot 	 3.0
    RE3533M 	 1.0
    FVSTETGLUtev 	 1
    r1812 	 -1.0
    GALASE5ly 	 0
    r2433 	 -1.0
    DSPVSteb 	 0.0
    S6TASE26ly 	 3.0
    r1616 	 -1.0
    r2302 	 -1.0
    HOCTDACBP 	 3.0
    H2O2tly 	 1
    HMGCOASi 	 3.0
    r2080 	 0.0
    UNK2 	 1
    FT 	 1
    RE2292C 	 1
    RE0937C 	 0.0
    r1254 	 2.0
    LEUB0AT3tc 	 -1.0
    EX_caribupglu_S_LPAREN_e_RPAREN_ 	 1
    ACGAM6PSi 	 1.0
    3HCO3_NAt 	 1.0
    ENGASE3ly 	 1
    GLYBtm 	 1
    M4BET2er 	 1.0
    DEDOLP1_U 	 1
    C4tmc 	 -1.0
    r2196 	 -1.0
    SPHINGStl 	 1
    ATVACIDhr 	 1
    FAOXC185C164m 	 3.0
    ALAR 	 1
    GGH_6DHFe 	 1.0
    EX_lstnm2_LPAREN_e_RPAREN_ 	 1
    LIMNENte 	 1
    ACOATA 	 -1.0
    AM1ACCShr 	 1
    MHISOR 	 0.0
    RE2459C 	 1
    MINOHPtn 	 1
    RE2318X 	 1.0
    PHEtec 	 -1.0
    PTE3x 	 0
    HDDACBP 	 3.0
    4HATVACIDthc 	 -1.0
    APAT2rm 	 -1.0
    ADPMAN 	 1.0
    ETOHMO 	 -1.0
    r2124 	 0.0
    r2227 	 -1.0
    MACOXO 	 3.0
    ARACHFATPtc 	 1.0
    r0541 	 -1.0
    4BHGLZitr 	 1
    EX_3mlda_LPAREN_e_RPAREN_ 	 1
    r0819 	 1
    ALAt4 	 3.0
    13DMTtu 	 1
    r2225 	 -1.0
    RE2914M 	 3.0
    PTVSTM3itr 	 1
    RE3347C 	 1
    RE3532M 	 1.0
    r0947 	 -1.0
    UMPK2n 	 3.0
    r0086 	 1.0
    RE2677N 	 1
    r1612 	 -1.0
    IMPC 	 1.0
    r1566 	 3.0
    EX_dha_LPAREN_e_RPAREN_ 	 1
    RE1711C 	 1
    RE3335R 	 2.0
    PYRSMCT1 	 -1.0
    r2294 	 -1.0
    TRIPVSteb 	 0.0
    RE3573X 	 3.0
    FPGS 	 -1.0
    biomass_DNA 	 3
    HCOUMARINte 	 1
    FAOXC183_6Z_9Z_12Zx 	 0.0
    EX_dcsptn1_LPAREN_e_RPAREN_ 	 1
    RE1796M 	 -1.0
    FOLt2 	 -1.0
    r1552 	 3.0
    r1989 	 -1.0
    RE3220L 	 -1.0
    12HTACRtep 	 1
    RE1539X 	 -1.0
    r1024 	 1
    r0512 	 -1.0
    EX_ala_L_LPAREN_e_RPAREN_ 	 3
    CYSGLYexR 	 -1.0
    r1710 	 -1.0
    EX_fucgalgbside_hs_LPAREN_e_RPAREN_ 	 1
    ENGASEly 	 1
    VACCCRNt 	 -1.0
    EX_txa2_LPAREN_e_RPAREN_ 	 1
    MAN1_6B1er 	 -1.0
    STS1 	 -1.0
    RE3624M 	 1
    FA161ACPH 	 -1.0
    RE0928C 	 -1.0
    FAOXC205C185m 	 3.0
    r2368 	 0
    HRETNtn 	 1
    r0245 	 3.0
    r1899 	 -1.0
    RE3174C 	 3.0
    EX_btn_LPAREN_e_RPAREN_ 	 3
    EX_bvite_LPAREN_e_RPAREN_ 	 1
    EX_acgam_LPAREN_e_RPAREN_ 	 3
    RE2998M 	 0.0
    NACHEX8ly 	 1.0
    UPPDC1 	 3.0
    EX_pydx5p_LPAREN_e_RPAREN_ 	 3
    r2297 	 -1.0
    METLEUex 	 -1.0
    EX_ha_LPAREN_e_RPAREN_ 	 3
    ALDD2x 	 3.0
    EX_retfa_LPAREN_e_RPAREN_ 	 3
    RE1587R 	 -1.0
    RE0828X 	 -1.0
    ACSRTNMT 	 -1.0
    PI4P5Kn 	 1
    ILEtec 	 -1.0
    RE1235C 	 1
    PIK4n 	 1
    r2371 	 -1.0
    r1965 	 -1.0
    RAI3 	 3.0
    RE2375C 	 1
    GGH_7DHFe 	 1.0
    r0879 	 -1.0
    RE3111C 	 1
    r0892 	 1
    r0755 	 2.0
    RE3106C 	 1
    r0629 	 -1.0
    RE3525X 	 1
    FATP4t 	 -1.0
    AM1CCSitr 	 1
    RE2319C 	 1
    RE3231R 	 1
    EX_HC01610_LPAREN_e_RPAREN_ 	 3
    EX_15dmt_LPAREN_e_RPAREN_ 	 1
    RE3148R 	 3.0
    EX_bz_LPAREN_e_RPAREN_ 	 1
    DURIK1m 	 -1.0
    RE3236R 	 1
    CERT1gt 	 0.0
    BTNt2m 	 0.0
    RE1050E 	 -1.0
    FAOXC180x 	 0.0
    ASNB0AT3tc 	 -1.0
    r0761 	 -1.0
    r2395 	 -1.0
    RE2718C 	 0
    RE1952C 	 0.0
    PGLYCABCte 	 -1.0
    r2011 	 -1.0
    LDH_D 	 -1.0
    TRDRm 	 0
    LNLCt 	 -1.0
    ENO 	 3.0
    r2293 	 -1.0
    r1685 	 -1.0
    ST6GALNAC21 	 3.0
    RE3033C 	 1
    4OHMDZhr 	 3.0
    13DMTtep 	 1
    IPDDIx 	 3.0
    EX_fucgalfucgalacglcgalgluside_hs_LPAREN_e_RPAREN_ 	 1
    FAOXC225C226m 	 3.0
    PRAGSr 	 -1.0
    RE2318C 	 1
    EPOXTACtev 	 1
    7THFtm 	 1
    EX_fucgal14acglcgalgluside_hs_LPAREN_e_RPAREN_ 	 1
    RE3421M 	 1
    RE1942C 	 1
    r0394 	 -1.0
    D_3AIBt 	 1
    CLS_hs 	 1.0
    r0547 	 2.0
    r0570 	 3.0
    NACHEXA21ly 	 -1.0
    r2437 	 -1.0
    OROTGLUt 	 -1.0
    AG13T14g 	 1.0
    r2353 	 0
    DASPO1p 	 -1.0
    DPROOp 	 -1.0
    EX_adn_LPAREN_e_RPAREN_ 	 3
    TSACMGLUCitr 	 1
    GLXO2p 	 -1.0
    P4507B11r 	 -1.0
    SUCCOAPET 	 -1.0
    RE1233M 	 2.0
    FAOXC123_3Z_6Z_9Zm 	 1.0
    XYLTD_Dr 	 1
    r2361 	 0
    RE1447M 	 1
    LTC4Sr 	 3.0
    PI5P3K 	 3.0
    RETNtr 	 1
    FRUt1r 	 0
    C10CRNe 	 1
    NCNt 	 -1.0
    LVSTOXD3Hhep 	 3.0
    EX_cholate_LPAREN_e_RPAREN_ 	 1
    RE1942R 	 0.0
    RE3435C 	 1
    dcmpt 	 1
    DOCO13EFATP 	 1.0
    FAOXC163C164Gm 	 3.0
    GLUDC 	 -1.0
    BAAT2x 	 -1.0
    r2110 	 0.0
    GLCNACT5g 	 -1.0
    ARTCOAL3 	 1
    DNDPt6m 	 -1.0
    RE3338X 	 3.0
    ACGBGBSIDEtl 	 1
    MI1456PKn 	 -1.0
    RE3552X 	 1
    SR5AR2r 	 3.0
    EX_am1accs_LPAREN_e_RPAREN_ 	 1
    31DMThr 	 1
    PTE4x 	 0
    PROSTGE2t2 	 -1.0
    EX_leuktrE4_LPAREN_e_RPAREN_ 	 3
    C16OHc 	 0
    CSPG_Ct 	 1
    CSNATp 	 0.0
    SLDxm 	 3.0
    SERDGLNexR 	 -1.0
    GHMT2rm 	 3.0
    FAOXC183C163m 	 3.0
    HPACtr 	 1
    S2TASE5ly 	 1
    UMPtr 	 1
    PHYHx 	 -1.0
    r1591 	 3.0
    RE3464C 	 1
    3HPVSCOAitm 	 1
    sink_tetdece1coa_LPAREN_c_RPAREN_ 	 1
    TMLYSOX 	 -1.0
    GALASE8ly 	 0
    RE1907C 	 1
    EX_hpdca_LPAREN_e_RPAREN_ 	 1
    r2178 	 -1.0
    RE3533R 	 3.0
    r1392 	 1.0
    FORtr 	 1
    r0196 	 1.0
    HXPRT 	 0.0
    ATVLACitr 	 1
    RE3334M 	 1
    EX_cdpea_LPAREN_e_RPAREN_ 	 1
    ABO3g 	 -1.0
    RSVtu 	 -1.0
    LPCOXp 	 -1.0
    6HMSMVACIDteb 	 1
    LTC4CP 	 1
    FAOXC81_5Zm 	 2.0
    RE3218R 	 -1.0
    FACOAL245_1 	 2.0
    SEASMETtn 	 1
    r1532 	 0.0
    r1182 	 1.0
    ELAIDCPT1 	 0
    ASPTAm 	 3.0
    RE2994X 	 3.0
    RE1709C 	 1
    ACYP 	 -1.0
    GALOR 	 3.0
    r2090 	 0.0
    EX_aflatoxin_LPAREN_e_RPAREN_ 	 1
    GPIDA2er 	 0.0
    C204CPT1 	 0
    LGNCCPT2 	 -1.0
    RE0915E 	 0
    RE2318M 	 1.0
    RE3498N 	 2.0
    MI1PS 	 -1.0
    INOSTO 	 -1.0
    LNSTLSr 	 3.0
    GP1Ctg 	 1
    VACCCPT1 	 0
    RE3151C 	 3.0
    RE3021C 	 1
    GLYC_St 	 1
    NACHEXA12ly 	 -1.0
    DM_datp_m_ 	 1
    GCC2am 	 -1.0
    PI345P5P 	 1.0
    r1498 	 -1.0
    RE3344M 	 3.0
    RE3459C 	 1
    C160CPT1 	 0
    MI134PK 	 1
    ACGAM2E 	 -1.0
    LSTNM2itr 	 1
    DURItn 	 1
    RE3587C 	 1
    CSPG_Et 	 1
    RETFAt1 	 1
    PVD3 	 1
    NTD7e 	 0.0
    GLPASE1 	 1.0
    RE3225R 	 3.0
    FVSTETGLUitr 	 1
    HEX4 	 3.0
    GCALDD 	 3.0
    r1259 	 2.0
    r0093 	 0.0
    ATVACIDhc 	 -1.0
    CLOHtex2 	 -1.0
    r0319 	 0.0
    FADDP 	 0.0
    r2306 	 -1.0
    r2115 	 0.0
    XOLTRIOLtr 	 1
    r2208 	 -1.0
    CYSAMOe 	 1
    EX_sprm_LPAREN_e_RPAREN_ 	 1
    FAOXC184x 	 1.0
    EX_ach_LPAREN_e_RPAREN_ 	 1
    RN0028X 	 0.0
    EX_hxan_LPAREN_e_RPAREN_ 	 3
    25HVITD3c 	 0.0
    RE1135C 	 -1.0
    r1668 	 -1.0
    EX_c6crn_ 	 1
    7AHGLZABCt 	 1
    FAOXC185m 	 1.0
    RE1807C 	 1
    r0743 	 3.0
    NTD9 	 3.0
    DRBK 	 -1.0
    r0002 	 1
    r1764 	 -1.0
    SARCOXp 	 -1.0
    CSPG_Bt 	 1
    EX_HC02199_LPAREN_e_RPAREN_ 	 1
    CYTK12n 	 3.0
    H2Otg 	 1
    GLUtr 	 1
    EX_HC00229_LPAREN_e_RPAREN_ 	 1
    B3GNT36g 	 -1.0
    RE3525R 	 1
    GLGNS1 	 1.0
    AM4NCStep 	 1
    RE3398M 	 3.0
    S3MEACMPhc 	 1
    r1815 	 -1.0
    LEUKTRD4tr 	 1
    r2140 	 -1.0
    ASPt6 	 1.0
    RE1806R 	 -1.0
    TTDCAtr 	 -1.0
    r1825 	 -1.0
    r2278 	 -1.0
    EX_dlnlcg_LPAREN_e_RPAREN_ 	 1
    r0692 	 -1.0
    RE0916C 	 1
    r2472 	 -1.0
    RE1508C 	 1
    P45011A1m 	 -1.0
    RE0936C 	 0.0
    RE3017R 	 1
    UGLCNACtg 	 3.0
    SPTix 	 -1.0
    S4T2g 	 -1.0
    r2355 	 0
    TCHOLAt2 	 -1.0
    RE2426C 	 1
    RE3444M 	 1.0
    TETPENT3CPT1 	 0
    RE3597M 	 0
    S2T4g 	 -1.0
    r0456 	 1
    RE2477C 	 1
    RE3476C 	 1.0
    T4HCINNMFM 	 1
    GTHDH 	 0.0
    RE0569C 	 3.0
    r2100 	 0.0
    RE0570C 	 1
    EX_dhap_LPAREN_e_RPAREN_ 	 1
    CRVSM1hc 	 3.0
    EX_5htrp_LPAREN_e_RPAREN_ 	 1
    PHCHGSm 	 1
    RE2626M 	 -1.0
    r0998 	 1
    RE3082X 	 3.0
    G14T16g 	 1.0
    P4502C18 	 -1.0
    r2056 	 -1.0
    r1155 	 1.0
    GLCNACPT_U 	 1
    r1767 	 -1.0
    GAO2g 	 1.0
    FAH3 	 -1.0
    PRISTtx 	 1
    ORNLEUrBATtc 	 -1.0
    CYTK4n 	 3.0
    3HLVSTAChep 	 -1.0
    RE3430X 	 1
    CRVS23M24hc 	 3.0
    RE3243R 	 1
    r0120 	 2.0
    EX_hom_L_LPAREN_e_RPAREN_ 	 1
    RE3389M 	 1
    RE2605C 	 1
    r1179 	 1.0
    EX_HC02208_LPAREN_e_RPAREN_ 	 1
    EX_dopa_LPAREN_e_RPAREN_ 	 3
    ALOX52 	 2.0
    MTHFD2m 	 3.0
    BTNt2 	 0.0
    FAOXC164x 	 1.0
    DCT 	 -1.0
    EX_am1c4n9cs_LPAREN_e_RPAREN_ 	 1
    DOLP_Uter 	 1
    GLCGLUT2 	 -1.0
    r2026 	 -1.0
    GALASE9ly 	 0
    r0797 	 -1.0
    r0153 	 -1.0
    FAOXC102_4Z_7Zx 	 -1.0
    r2369 	 0
    NADtx 	 1
    CRNtx 	 1
    RE3112C 	 1
    XOLTRI25te 	 1
    RE3568C 	 1
    URAt 	 -1.0
    r0121 	 2.0
    PTVSTM13hr 	 -1.0
    DESAT22_1p 	 1
    HEXCCRNt 	 -1.0
    RE1582R 	 -1.0
    RE2048R 	 1
    RE3270C 	 1
    RE3273C 	 -1.0
    RE1815X 	 1.0
    GCCbim 	 -1.0
    r0062 	 1
    DHCR243r 	 3.0
    35DSMVhep 	 3.0
    3HPVSitr 	 1
    SUCRe 	 -1.0
    G16MT_U 	 1
    r2240 	 -1.0
    r1323 	 1
    CRVSM1SPhc 	 1
    HBZOPT10m 	 0.0
    PTVSTGLUChc 	 0
    3HPVSteb 	 0.0
    PEFLIPm 	 2.0
    r2326 	 3.0
    3HLYTCL 	 -1.0
    r2010 	 -1.0
    PRO1xm 	 1
    EX_cit_LPAREN_e_RPAREN_ 	 3
    NACHEX4ly 	 1.0
    r2324 	 3.0
    FAEL205 	 3.0
    RE3399M 	 3.0
    PYLALDOX 	 3.0
    S6TASE19ly 	 3.0
    r1026 	 -1.0
    PHEMEABCte 	 0.0
    EX_tcynt_LPAREN_e_RPAREN_ 	 1
    r1449 	 3.0
    FAOXC5C5DCc 	 0
    AMETtd 	 1
    DEBRISOQUINEt 	 -1.0
    CYTK5n 	 3.0
    ST8SIA54g 	 -1.0
    TYRt 	 -1.0
    NTPP9 	 3.0
    MDZtu 	 -1.0
    NDPK9n 	 0
    r1791 	 -1.0
    ACMPthc 	 -1.0
    RETt 	 -1.0
    r2513 	 1
    FCOAH 	 1
    AM1CSAitr 	 1
    DOLGPP_Uer 	 1
    TDPGDH 	 -1.0
    r0239 	 -1.0
    FAOXC141_7Em 	 3.0
    GALT2g 	 0.0
    11DOCRTSTRNtr 	 1
    RE1516X 	 1
    IBUPGLUCtpvb 	 1
    GALASE13ly 	 0
    PGSr 	 2.0
    r1807 	 -1.0
    RETNt 	 1
    EX_5mta_LPAREN_e_RPAREN_ 	 1
    r0408 	 3.0
    VITD3tm 	 1
    AM4N9CShr 	 1
    P4503A43r 	 -1.0
    GPDDA1 	 1
    S6T13g 	 -1.0
    RE3337X 	 -1.0
    RE3491C 	 1
    GLXtp 	 1
    r1862 	 -1.0
    NAt3_1g 	 -1.0
    NORANMT 	 -1.0
    r2390 	 -1.0
    PAIL45P_HStn 	 1
    ARGNm 	 3.0
    TETPENT6CPT2 	 -1.0
    RE1233C 	 -1.0
    r0614 	 -1.0
    S6T15g 	 1.0
    r0835 	 -1.0
    DOPAQNISO1 	 1
    VITD3t 	 1
    C162ACBP 	 3.0
    RE3220R 	 -1.0
    NNDPR 	 -1.0
    G14T21g 	 1.0
    r1682 	 -1.0
    r2389 	 -1.0
    THMATPe 	 1
    ALOX12R 	 3.0
    PECDCHe 	 1
    10FTHF6GLUtl 	 1
    RE2596C 	 1
    EX_glz_LPAREN_e_RPAREN_ 	 1
    FAOXC184_3Z_6Z_9Z_12Zm 	 1.0
    GK1 	 1.0
    It 	 -1.0
    FAOXC164GC163m 	 3.0
    LVSTtu 	 1
    GLYtm 	 1
    RE3346C 	 3.0
    GTPCIn 	 2.0
    EX_6hlvst_LPAREN_e_RPAREN_ 	 1
    3MOPte 	 1
    HMGLm 	 -1.0
    CYTK3n 	 3.0
    r2119 	 1
    EX_meracmp_LPAREN_e_RPAREN_ 	 1
    GPAM_hs 	 2.0
    RE3338C 	 3.0
    TMNDNCt 	 1
    ST6GALNAC27 	 3.0
    r1585 	 3.0
    PI4P5K 	 -1.0
    r2275 	 -1.0
    r2046 	 -1.0
    EX_whddca_LPAREN_e_RPAREN_ 	 1
    RE3081X 	 -1.0
    EX_12HPET_LPAREN_e_RPAREN_ 	 1
    CRVSitr 	 1
    r0388 	 -1.0
    EX_psylchol_LPAREN_e_RPAREN_ 	 1
    TAURCHAe 	 1
    r1974 	 -1.0
    FADH2tx 	 1
    r2040 	 -1.0
    LYSMTF3n 	 3.0
    PROSTGI2t 	 1
    LTDCL 	 -1.0
    RE3036C 	 1
    PPAt 	 1
    TACRitr 	 1
    RE3430M 	 1
    RTOTALCRNCPT2 	 -1.0
    r2267 	 -1.0
    G14T8g 	 1.0
    DNDPt5m 	 -1.0
    AM4NCSitr 	 1
    HMGLx 	 -1.0
    RE3145X 	 3.0
    NDPK5 	 0
    r1745 	 -1.0
    FOLR2 	 -1.0
    RE3381C 	 -1.0
    DECCRNe 	 1
    TKT1 	 1.0
    DURAD2 	 2.0
    AACOAT 	 1.0
    r2231 	 -1.0
    RE3385M 	 1
    PEPLYStn 	 1
    AM1CCSteb 	 1
    RE0579X 	 -1.0
    URIt 	 -1.0
    FAOXC103C102x 	 -1.0
    RE3474C 	 1
    EX_HC02206_LPAREN_e_RPAREN_ 	 1
    7BHGLZGLCitr 	 1
    C10DCc 	 0
    r1765 	 -1.0
    F6Tg 	 -1.0
    EX_galacglcgalgbside_hs_LPAREN_e_RPAREN_ 	 1
    NACHEX12ly 	 1.0
    RE2990X 	 -1.0
    UMPK7 	 3.0
    FAOXC122C101m 	 0.0
    RE2040C 	 1
    CYTDn 	 -1.0
    RE3564M 	 3.0
    r0934 	 0.0
    r1929 	 -1.0
    ATPS4m 	 -1.0
    DM_atp_c_ 	 3
    r1979 	 -1.0
    FAOXC163_7Z_10Z_13Zm 	 3.0
    FUT98g 	 -1.0
    r2215 	 -1.0
    r1251 	 2.0
    RE2626C 	 1
    3MEACMPhc 	 1
    EHGLAT 	 1.0
    LEUyLAThtc 	 -1.0
    r1843 	 -1.0
    THRCYSNaEx 	 3.0
    THP2Ctp 	 1
    PItn 	 1
    FADH2ETC 	 1.0
    Am1CSAteb 	 1
    EX_orn_LPAREN_e_RPAREN_ 	 3
    r1304 	 1
    DCATDc 	 -1.0
    MI13456PK 	 1
    RE2637C 	 -1.0
    r1922 	 -1.0
    LCTStl 	 1
    NDPK8m 	 -1.0
    PGK 	 3.0
    NACHEXA20ly 	 -1.0
    NTMELYStner 	 1
    GLNTHRNaEx 	 -1.0
    r0210 	 -1.0
    FAOXC164_4Z_7Z_10Z_13Zx 	 -1.0
    3HPVSTETCOAhcm 	 3.0
    RE1521M 	 3.0
    3MEACMPitr 	 1
    r2047 	 -1.0
    RE3574X 	 -1.0
    NT5C 	 3.0
    RE3626M 	 1
    UMPK6 	 3.0
    EX_sulpacmp_LPAREN_e_RPAREN_ 	 1
    RE3596C 	 1.0
    r1527 	 1
    r0681 	 -1.0
    EX_tetpent3_LPAREN_e_RPAREN_ 	 1
    RE3408C 	 1
    CARBIBUP_SGLUthv 	 1
    FAOXC226205x 	 -1.0
    DNDPt48m 	 -1.0
    EX_pro_L_LPAREN_e_RPAREN_ 	 3
    CYSALANaEx 	 3.0
    DNAMTSEn 	 2.0
    TMDM1itr 	 1
    LCTStg 	 1
    NACHEX9ly 	 1.0
    r1586 	 3.0
    34DHPLACOX_NADP_ 	 2.0
    ABO8g 	 -1.0
    r2505 	 1.0
    1331TACRteb 	 1
    PERILLYLte 	 1
    RE3340M 	 3.0
    r2197 	 -1.0
    PDE4n 	 -1.0
    RE1978C 	 1
    RN0030R 	 -1.0
    EX_tre_LPAREN_e_RPAREN_ 	 3
    WHHDCAte 	 1
    SPODMx 	 3.0
    r0539 	 1.0
    MECOAS1m 	 1
    EX_pe_hs_LPAREN_e_RPAREN_ 	 3
    EX_spmd_LPAREN_e_RPAREN_ 	 3
    r1958 	 -1.0
    r1589 	 3.0
    RE3273R 	 3.0
    GBA 	 0
    ALASERNaEx 	 3.0
    RE3372C 	 1
    EX_gln_L_LPAREN_e_RPAREN_ 	 3
    RE3092X 	 -1.0
    DM_dttp_n_ 	 3
    DGNSKm 	 1.0
    r0678 	 -1.0
    NADtn 	 1
    ASNGLNNaEx 	 -1.0
    HEX10 	 3.0
    r1368 	 1
    RSVLACitr 	 1
    RE2335C 	 1
    RE3001M 	 3.0
    GTHP 	 3.0
    PETOHMr_hs 	 -1.0
    PRGNLONEtr 	 1
    PIK5 	 -1.0
    RE3453C 	 1
    6CSMVACIDhep 	 -1.0
    r2054 	 -1.0
    r2365 	 0
    r1887 	 -1.0
    RE3184M 	 2.0
    RE1807M 	 -1.0
    HS4ly 	 1.0
    RE0702M 	 1
    RE3261R 	 0
    EX_tymsf_LPAREN_e_RPAREN_ 	 1
    ORETNtn 	 1
    TTDCPT2 	 -1.0
    EX_3ttetddcoacrn_ 	 1
    EX_fol_LPAREN_e_RPAREN_ 	 3
    r1846 	 -1.0
    EX_asp_L_LPAREN_e_RPAREN_ 	 3
    r2123 	 1
    H2O2t 	 1
    RE0908G 	 -1.0
    GALTg 	 -1.0
    PNTORDe 	 1
    r2187 	 -1.0
    CYTK11 	 3.0
    FAEL204 	 3.0
    r0714 	 -1.0
    RE3241C 	 3.0
    r1378 	 3.0
    RE1846C 	 -1.0
    RE1311M 	 -1.0
    25HVITD3tin 	 -1.0
    r2486 	 -1.0
    r2387 	 -1.0
    EX_cl_LPAREN_e_RPAREN_ 	 3
    NACHEX5ly 	 1.0
    VALTA 	 -1.0
    RE3169C 	 1
    EX_for_LPAREN_e_RPAREN_ 	 3
    RN0022R 	 0.0
    DM_taur_LPAREN_c_RPAREN_ 	 1
    GUACYC 	 0.0
    r1740 	 -1.0
    NTPP11 	 3.0
    EX_thm_LPAREN_e_RPAREN_ 	 3
    RE1309C 	 1
    FADH2tru 	 1
    RE3151R 	 3.0
    ADNK3 	 1
    RE3074X 	 2.0
    A4GALTg 	 -1.0
    6BHGLZABCt 	 1
    CYTK3 	 3.0
    RE3119C 	 1
    SRTNACT 	 -1.0
    PS_HSter 	 1
    RE3015C 	 1
    RE1954C 	 1
    6EPSteb 	 0.0
    ECGISOr 	 1
    RE3125R 	 1
    RE3631C 	 1
    DHGLZhc 	 1
    EX_HC01104_LPAREN_e_RPAREN_ 	 1
    RE3018R 	 -1.0
    FAOXC123_3Z_6Z_9Zx 	 3.0
    RE2249C 	 1
    r0973 	 1
    B3GALT3g 	 -1.0
    NCKt 	 3.0
    r2148 	 -1.0
    r2048 	 -1.0
    THYPX 	 -1.0
    r0611 	 1
    MTRI 	 1
    r1520 	 3.0
    XYLt 	 -1.0
    2HIBUPGLUC_Sitr 	 1
    r0595 	 3.0
    CSASULPhc 	 1
    31DMTitr 	 1
    ARACHCPT1 	 0
    NDPK10 	 0
    FAOXC226C225m 	 3.0
    RE1846X 	 -1.0
    ACGAMK 	 2.0
    C51CPT1 	 0
    1531TACRtev 	 1
    EX_tmdm5_LPAREN_e_RPAREN_ 	 1
    P450LTB4r 	 1
    r2509 	 1
    GGH_6THFl 	 1.0
    ADPRDPm 	 1.0
    RETNGLCt 	 1
    FRDPtcr 	 1
    biomass_protein 	 3
    RE3629C 	 1
    r0317 	 3.0
    r0202 	 -1.0
    RE3432M 	 0
    r1646 	 3.0
    RE3086X 	 3.0
    ATPtm 	 3.0
    FAEL183 	 3.0
    EX_ptth_LPAREN_e_RPAREN_ 	 1
    AACTOOR 	 0.0
    BAMPPALDOXm 	 3.0
    A4GALTc 	 -1.0
    RE2443M 	 1.0
    LEUKTRA4t 	 0
    FACOAL260 	 2.0
    UMPK3n 	 3.0
    EX_5mthf_LPAREN_e_RPAREN_ 	 1
    3HBCOAHLm 	 1
    TTDCAFATPtc 	 1.0
    r2174 	 -1.0
    r2301 	 -1.0
    r1164 	 -1.0
    MOGAT 	 -1.0
    LSTNM5hr 	 3.0
    EX_pectintchol_LPAREN_e_RPAREN_ 	 1
    6AHGLZhr 	 1
    MI3456PK 	 1
    CRNtHa 	 -1.0
    RE2029C 	 1
    EX_epoxtac_LPAREN_e_RPAREN_ 	 1
    r2247 	 -1.0
    PMI1346PH 	 -1.0
    CRVS1tev 	 1
    ADNt 	 -1.0
    r1080 	 1
    FUCGALFUCGALACGLCGALGLUSIDEte 	 1
    1513DTACRhr 	 1
    LCADim 	 3.0
    PUNP5 	 2.0
    r0660 	 3.0
    DURAD 	 2.0
    DPPS 	 1
    RE0827C 	 0.0
    ARTFR51 	 1
    EX_fad_LPAREN_e_RPAREN_ 	 1
    FAOXC2452253x 	 0.0
    TYRDOPO 	 -1.0
    BPNT2 	 -1.0
    EX_HC01609_LPAREN_e_RPAREN_ 	 3
    EX_gthox_LPAREN_e_RPAREN_ 	 3
    THSACMPitr 	 1
    r2492 	 -1.0
    MANtg 	 1
    4HATVACIDtep 	 1
    TYRDOPO3 	 -1.0
    r1797 	 -1.0
    r2364 	 0
    RE3326M 	 -1.0
    RE2349M 	 -1.0
    RE1533X 	 3.0
    r1831 	 -1.0
    r1785 	 -1.0
    FAOXC81C61m 	 0.0
    RE1860C 	 1
    r2190 	 -1.0
    BILIRUBt2 	 -1.0
    RE3272N 	 1
    RE3134C 	 0.0
    5MTAte 	 1
    S6TASE23ly 	 3.0
    PCHOLPr_hs 	 3.0
    CERT2rt 	 0.0
    RE1812C 	 1
    S6TASE18ly 	 3.0
    MAL_Ltx 	 1
    RE3475C 	 1
    RE1514X 	 2.0
    r0752 	 2.0
    GALASE4ly 	 0
    BZt 	 1
    5HTRPVESSEC 	 0.0
    10FTHF7GLUtl 	 1
    r2251 	 -1.0
    r0433 	 3.0
    r1574 	 3.0
    CATm 	 3.0
    FACOAL226 	 2.0
    r2014 	 -1.0
    r0557 	 0
    RE3387M 	 3.0
    G1PTT 	 1
    RE2221M 	 -1.0
    PI4P3K 	 3.0
    GLDBRAN 	 1.0
    RE3420C 	 1
    UDPRIBc 	 1
    EX_C04849_LPAREN_e_RPAREN_ 	 1
    AGTim 	 -1.0
    PFK26 	 3.0
    EX_HC01440_LPAREN_e_RPAREN_ 	 1
    RN0021C 	 0.0
    r1521 	 3.0
    r2338 	 -1.0
    THYMt 	 -1.0
    TMDM5itr 	 1
    DEOXFVSitx 	 1
    GALGALFUCFUCGALACGLCGALACGLCGAL14ACGLCGALGLUSIDEte 	 1
    EX_his_L_LPAREN_e_RPAREN_ 	 3
    CHOLtu 	 -1.0
    r0737 	 0
    RE2081C 	 1
    GLACOm 	 3.0
    RE2235C 	 1
    S2T1g 	 1.0
    NRVNCCPT2 	 -1.0
    XOL7AH2tr 	 1
    RE3234C 	 1
    DSREDUCr 	 3.0
    EX_lstn1gluc_LPAREN_e_RPAREN_ 	 1
    C161CPT1 	 0
    ANDRSTRNGLCte 	 1.0
    r2271 	 -1.0
    NMNATn 	 -1.0
    SUCOASm 	 0.0
    r1437 	 1
    CHATn 	 -1.0
    PYDAMK 	 3.0
    DEDOLP2_U 	 1
    CHOLD2m 	 -1.0
    TETTET6COAtx 	 1
    TMDPK 	 -1.0
    EX_digalsgalside_hs_LPAREN_e_RPAREN_ 	 1
    ACACT1rm 	 3.0
    FAOXC123C102m 	 0.0
    PVSGLUCteb 	 0.0
    DRPA 	 3.0
    RE2360C 	 1
    URIt2r 	 0.0
    METtec 	 -1.0
    RETNtr2 	 1
    H8TAer 	 -1.0
    EX_HC02192_LPAREN_e_RPAREN_ 	 1
    AM1ALCSteb 	 1
    RE1827C 	 1
    1PPDCRp 	 1
    EX_hcoumarin_LPAREN_e_RPAREN_ 	 1
    SBTle 	 1
    2HATVACIDhc 	 -1.0
    RE3557C 	 1
    RE2387C 	 1
    DCMPDA 	 1.0
    PHEMEe 	 1
    r2372 	 -1.0
    EX_ak2lgchol_hs_LPAREN_e_RPAREN_ 	 1
    r1063 	 1
    ENMAN2g 	 1
    SCP22x 	 3.0
    r1963 	 -1.0
    DHDPBMTm 	 -1.0
    GK1m 	 1.0
    DMHPTCRNte 	 1
    r1253 	 2.0
    r1411 	 0
    FACOAL205 	 2.0
    ALDD21 	 3.0
    MTHFTe 	 1
    RE3242R 	 1
    C120CRNt 	 -1.0
    CRVSM24itr 	 1
    r0984 	 -1.0
    C12DCACOT 	 -1.0
    DOLGLCP_Lter 	 1
    C8CRNe 	 1
    r2419 	 -1.0
    EX_galfuc12gal14acglcgalgluside_hs_LPAREN_e_RPAREN_ 	 1
    FAOXC101C8x 	 -1.0
    r1669 	 -1.0
    DHGLZABCt 	 1
    RE2897C 	 1
    ASPTA 	 1.0
    ESTRSABCtc 	 1.0
    RE3412C 	 1
    RE3567C 	 1
    SMPD4 	 0.0
    EX_etoh_LPAREN_e_RPAREN_ 	 3
    r2402 	 -1.0
    N2M2NMASNt 	 1
    RETNGLCt2 	 1
    NADK 	 1.0
    NKCC2t 	 0.0
    BILDGLCURt 	 -1.0
    RE3227R 	 1
    r1945 	 -1.0
    EX_HC02194_LPAREN_e_RPAREN_ 	 1
    OXYP2CONJ 	 1
    r1971 	 -1.0
    EX_hx_LPAREN_e_RPAREN_ 	 1
    RN0032C 	 1
    FAOXC10080x 	 0.0
    CRVStu 	 1
    OAADC 	 1
    ADPtx 	 1
    PSYTCHe 	 1
    GMPS2 	 3.0
    RE1538X 	 -1.0
    EX_HC02202_LPAREN_e_RPAREN_ 	 1
    UDPtl 	 1
    r1150 	 1
    r1633 	 3.0
    EX_fucfucgalacglcgalgluside_hs_LPAREN_e_RPAREN_ 	 1
    r1857 	 -1.0
    C81CPT1 	 0
    RE2644C 	 -1.0
    4ABUTtcn 	 1
    FAOXC101m 	 1.0
    DCYTD 	 1.0
    r1981 	 -1.0
    2HATVACIDtep 	 1
    PSP_L 	 -1.0
    C142OHe 	 1
    FUCGALFUCGALACGLCGALGLUSIDEtg 	 1
    DM_dgtp_m_ 	 1
    ASPte 	 1
    STRDNCCRNt 	 -1.0
    EX_CE5798_LPAREN_e_RPAREN_ 	 1
    r0579 	 3.0
    RE2722C 	 0
    r2388 	 -1.0
    SARCStp 	 1
    DNDPt53m 	 -1.0
    G5SADrm 	 1
    GLXO1 	 3.0
    r1691 	 -1.0
    EX_whtststerone_LPAREN_e_RPAREN_ 	 1
    PItx 	 1
    RE2327C 	 1
    ALDD2y 	 3.0
    r2084 	 1
    r2017 	 -1.0
    RE1518X 	 1
    RE3224R 	 3.0
    r1667 	 -1.0
    PSt3 	 1
    RETH2 	 1
    r1571 	 3.0
    C226CPT2 	 -1.0
    r0181 	 3.0
    IMPD 	 3.0
    C161CRN2t 	 -1.0
    NACHORCTL3le 	 -1.0
    FAOXC6040m 	 2.0
    RE0830C 	 1
    r2101 	 1
    RE3500R 	 0.0
    r1303 	 1
    YVITEt 	 1
    r0594 	 1.0
    ACGALtly 	 1
    r1851 	 -1.0
    RE3500X 	 0.0
    EGMESTr 	 0.0
    EX_glygn5_LPAREN_e_RPAREN_ 	 1
    DOLPGT1_Ler 	 -1.0
    5HOXINDACTO2OX 	 -1.0
    r0220 	 0.0
    AVITE2t 	 1
    PYDXDH 	 -1.0
    TRIODTHYSUFt 	 1
    r0792 	 -1.0
    RE3519C 	 0.0
    PROD2m 	 0
    AVITE1t 	 -1.0
    EX_2mcit_LPAREN_e_RPAREN_ 	 1
    EAFLATOXINte 	 1
    RE3445M 	 1
    r0354 	 3.0
    DHAPtc 	 1
    GBGT1 	 -1.0
    GALGT2 	 -1.0
    GNDc 	 1.0
    FK 	 -1.0
    CYTDt4 	 0.0
    DM_yvite_c_ 	 1
    DTMPKm 	 -1.0
    EX_10fthf_LPAREN_e_RPAREN_ 	 1
    r0812 	 1
    NDP8 	 -1.0
    q10h2tc 	 1
    TSTSTERONEtr 	 1
    DNDPt30m 	 -1.0
    MI1PP 	 3.0
    RETH 	 -1.0
    RN0028R 	 0.0
    TETPENT6CRNt 	 -1.0
    ACALDtr 	 1
    RE2034C 	 1
    FAOXC185_3Z_6Z_9Z_12Z_15Zx 	 3.0
    BETALDHxm 	 1
    EX_crmp_hs_LPAREN_e_RPAREN_ 	 1
    ALAB0AT3tc 	 -1.0
    GLPASE2 	 1.0
    GTMLTe 	 -1.0
    PRFGS 	 -1.0
    r1642 	 3.0
    FAOXC7050m 	 2.0
    PROt2r 	 0
    P5CR 	 -1.0
    4HBZFm 	 1
    PYK 	 3.0
    FAOXC162C142m 	 3.0
    MI145PK 	 3.0
    r2151 	 -1.0
    VACCCPT2 	 -1.0
    SMS 	 -1.0
    RE3335C 	 -1.0
    PI4P3Ker 	 3.0
    CPCTDTX 	 0.0
    FAOXC102C103m 	 2.0
    NDPK10m 	 -1.0
    EX_ppi_LPAREN_e_RPAREN_ 	 1
    FAOXC181_9Zm 	 3.0
    ACSm 	 -1.0
    14MDZtev 	 1
    GLYOX 	 -1.0
    RE1234C 	 1
    EX_din_LPAREN_e_RPAREN_ 	 3
    EX_dhglz_LPAREN_e_RPAREN_ 	 1
    NDPK3 	 0
    M16NTg 	 0.0
    ARTFR12 	 1
    S6TASE25ly 	 -1.0
    RE2992M 	 1.0
    RE2660N 	 -1.0
    GTHS 	 3.0
    r1980 	 -1.0
    FAOXC180 	 -1.0
    RE3226R 	 3.0
    EX_ps_hs_LPAREN_e_RPAREN_ 	 3
    FAS160COA 	 -1.0
    r1547 	 3.0
    r0686 	 -1.0
    BAAT4x 	 -1.0
    Kt3g 	 -1.0
    EX_bhb_LPAREN_e_RPAREN_ 	 3
    SERtN1 	 -1.0
    r0970 	 1
    OMEPRAZOLEte 	 1
    ORNtiDF 	 3.0
    Coqe 	 1
    RE0583N 	 1
    r1551 	 3.0
    CRVNCtr 	 1
    r2263 	 -1.0
    r1533 	 0.0
    4OHPROIMINOtc 	 -1.0
    r1367 	 1
    GLNyLATthc 	 -1.0
    RE2975M 	 -1.0
    ACMPGLUChr 	 0
    PYRt2p 	 1.0
    r2093 	 1
    r2518 	 2.0
    RE3345R 	 2.0
    THMPPtm 	 1
    EX_tripvs_LPAREN_e_RPAREN_ 	 1
    RE2520C 	 1
    r1735 	 -1.0
    r2206 	 -1.0
    RE2386C 	 1
    1513TACRitr 	 1
    AG13T13g 	 1.0
    RE3500C 	 0.0
    PPAn 	 1
    GLYPHEHYc 	 1
    LGNCCPT1 	 0
    FPGS3m 	 -1.0
    r0768 	 -1.0
    PYRt2r 	 1.0
    FAOXC204C205x 	 0.0
    HSD3B13 	 -1.0
    ST3GAL62g 	 -1.0
    r2200 	 -1.0
    RE2916M 	 3.0
    r1861 	 -1.0
    r1171 	 -1.0
    SARCStex 	 1
    MGSA 	 1
    r2370 	 -1.0
    PYAM5POr 	 -1.0
    BTNTe 	 1
    FAOXC182_9E_12Em 	 3.0
    HTDCRNe 	 1
    OXYPR7tehv 	 1
    NAPQIhr 	 1
    r1262 	 2.0
    r1829 	 -1.0
    7HPVShc 	 3.0
    EX_docosac_LPAREN_e_RPAREN_ 	 1
    EX_HC02214_LPAREN_e_RPAREN_ 	 1
    r1994 	 -1.0
    DESAT18_4 	 3.0
    PTVSTtep 	 1
    r1018 	 1
    4MPTNLtr 	 1
    RE3457C 	 1
    ACHEe 	 -1.0
    FE2t 	 1
    AMETtn 	 1
    DALAt2r 	 -1.0
    FPGS5m 	 -1.0
    RE2718G 	 -1.0
    1MNCAMti 	 1
    DTMPK 	 -1.0
    ASAH1 	 3.0
    EX_3mob_LPAREN_e_RPAREN_ 	 3
    RE3242C 	 1
    GPIDAer 	 0.0
    GLYK 	 -1.0
    RE3488C 	 1
    RE3004M 	 3.0
    r2165 	 -1.0
    FBA 	 3.0
    TRIODTHYSULT 	 1.0
    r1907 	 -1.0
    RE1828M 	 -1.0
    RETFA 	 1
    RN0020C 	 1
    TLACFVStev 	 1
    RE0919C 	 -1.0
    DNDPt33m 	 1
    PROIMINOtc 	 -1.0
    RE3308M 	 1
    RE3515C 	 1
    RE3124C 	 3.0
    EX_mdz_LPAREN_e_RPAREN_ 	 1
    SADT 	 3.0
    RE3229C 	 1
    3SALAASPm 	 0.0
    6HMSMVitr 	 1
    EX_pan4p_LPAREN_e_RPAREN_ 	 1
    MDZitr 	 1
    FPGS9m 	 -1.0
    EX_2hibupglu_S_LPAREN_e_RPAREN_ 	 1
    CATp 	 3.0
    RE1100R 	 -1.0
    RE0066R 	 -1.0
    RE3526M 	 1
    LVSTOXD6METhep 	 3.0
    EX_nfdlac_LPAREN_e_RPAREN_ 	 1
    CYSSNAT4te 	 -1.0
    N3Tg 	 -1.0
    CHAT 	 -1.0
    NTD3l 	 2.0
    RE3152C 	 1
    SRTNtu 	 -1.0
    EX_acmpglu_LPAREN_e_RPAREN_ 	 1
    AG13T17g 	 1.0
    PYDXNtr 	 1
    RE3520N 	 1
    MALTSULtm 	 -1.0
    SUCCACT 	 1
    r1140 	 1
    r2122 	 0.0
    DOGULNO1 	 1
    FUC13GALACGLCGAL14ACGLCGALGLUSIDEtg 	 1
    SUCCt4_2 	 0.0
    GHMT3 	 -1.0
    S4TASE4ly 	 -1.0
    ENGASE2 	 3.0
    S2TASE1ly 	 0.0
    RADH2 	 1
    r0166 	 3.0
    RE3006M 	 3.0
    GSNtl 	 -1.0
    EX_HC01787_LPAREN_e_RPAREN_ 	 1
    Htg 	 1
    METB0AT2tc 	 -1.0
    r0389 	 -1.0
    ST6GALNAC22 	 3.0
    EBASTINEOHte 	 1
    ARAB_Lt 	 1
    r1935 	 -1.0
    GGT6 	 1
    r1554 	 3.0
    FUC14GALACGLCGALGLUSIDEtg 	 1
    r2299 	 -1.0
    RE2916X 	 3.0
    FPGSm 	 -1.0
    GLUDxm 	 3.0
    r0603 	 -1.0
    C10OHc 	 0
    RE0579M 	 1.0
    PROB0AT2tc 	 1
    C141CPT1 	 0
    EX_HC02203_LPAREN_e_RPAREN_ 	 1
    GD1B2te 	 1
    6CSMVACIDteb 	 1
    C12OHc 	 0
    FAOXC164C163x 	 -1.0
    EX_fum_LPAREN_e_RPAREN_ 	 3
    RE2987X 	 3.0
    ESTRONESt2 	 -1.0
    CRTSTRNtr 	 1
    r2363 	 0
    5HOXINOXDA 	 3.0
    EX_gtp_LPAREN_e_RPAREN_ 	 1
    ATVACYLGLUChc 	 0
    EX_35dhpvs_LPAREN_e_RPAREN_ 	 1
    ATPtx 	 1
    TDP 	 1
    r1806 	 -1.0
    NACHEX22ly 	 1.0
    UAG4E 	 0.0
    EX_onpthl_LPAREN_e_RPAREN_ 	 1
    CSNAT3x 	 0.0
    FAOXC161C141x 	 -1.0
    THRt4 	 3.0
    FAOXC16BRx 	 -1.0
    r1776 	 -1.0
    LFORKYNHYD 	 1.0
    r1670 	 -1.0
    GNDer 	 1.0
    ABO7g 	 -1.0
    RE2513L 	 3.0
    NTD9e 	 0.0
    EX_urate_LPAREN_e_RPAREN_ 	 1
    EX_lvstacid_LPAREN_e_RPAREN_ 	 1
    EX_cmp_LPAREN_e_RPAREN_ 	 1
    SR5ARr 	 3.0
    FAOXC184_3Z_6Z_9Z_12Zx 	 3.0
    RE3178M 	 3.0
    TETPENT3CPT2 	 -1.0
    MI145P6Kn 	 -1.0
    RE3225C 	 3.0
    RE2156M 	 1
    3MOX4HOXPGALDOX_NADP_ 	 2.0
    r1144 	 1
    EX_tagat_D_LPAREN_e_RPAREN_ 	 1
    CYANt 	 1
    RE1920C 	 1
    r0310 	 3.0
    r1302 	 1
    SQLEr 	 2.0
    FAOXC161_7Zm 	 3.0
    RE2680C 	 1
    CLI2tex 	 -1.0
    2MB2COAc 	 3.0
    RE3388M 	 3.0
    GLUPRT 	 -1.0
    r0280 	 -1.0
    AHCYStn 	 1
    4OHMDZitr 	 1
    RBK 	 -1.0
    RE3095L 	 3.0
    RE3339X 	 -1.0
    r1676 	 -1.0
    r2392 	 -1.0
    r1828 	 -1.0
    FUMSO3tm 	 -1.0
    RE2675C 	 1
    r0829 	 -1.0
    r0650 	 -1.0
    PSFLIPm 	 2.0
    B3GNT51g 	 3.0
    r1783 	 -1.0
    MLTG1 	 -1.0
    CRTSLt 	 1
    NDPK10n 	 0
    r2171 	 -1.0
    6THFtm 	 1
    DCTPtn 	 1
    KSIt 	 1
    GGT5r 	 -1.0
    CEPTC 	 2.0
    r2358 	 0
    FAOXC141C141OHm 	 3.0
    RE2080C 	 1
    r1695 	 -1.0
    2HATVACIDitr 	 1
    B3GNT11g 	 -1.0
    EX_crm_hs_LPAREN_e_RPAREN_ 	 1
    AM1ACSteb 	 1
    r0709 	 2.0
    PAN4PPe 	 1
    DM_bvite_c_ 	 1
    FAOXC140120x 	 0.0
    SCP21x 	 3.0
    H7ET2er 	 1.0
    r2125 	 0.0
    RE1317C 	 1
    LDH_Lm 	 3.0
    NACSMCTte 	 -1.0
    r2376 	 -1.0
    LRAT2 	 1
    r1014 	 1
    ARACHDt2 	 -1.0
    C60CPT2 	 -1.0
    RE2028C 	 1
    DM_datp_n_ 	 3
    HEX1 	 3.0
    2HATVACIDGLUCitr 	 1
    DM_btn 	 1
    ACTLMO 	 -1.0
    EX_HC01444_LPAREN_e_RPAREN_ 	 1
    EX_3thexddcoacrn_ 	 1
    ALDSTRNte 	 1
    STACMPhc 	 1
    RE1537C 	 1
    r2207 	 -1.0
    LNLNCGt 	 -1.0
    r2449 	 0.0
    EX_glyc3p_LPAREN_e_RPAREN_ 	 3
    LSTNM4hr 	 1
    r1581 	 3.0
    EX_adp_ 	 1
    EX_prpp_LPAREN_e_RPAREN_ 	 1
    EX_tdcrn_ 	 1
    EX_lnlncg_LPAREN_e_RPAREN_ 	 3
    r1701 	 -1.0
    EX_met_L_LPAREN_e_RPAREN_ 	 3
    EX_am1a4ncs_LPAREN_e_RPAREN_ 	 1
    RE3239C 	 1
    RE3432C 	 1.0
    EX_lys_L_LPAREN_e_RPAREN_ 	 3
    r0287 	 3.0
    SELADT 	 3.0
    RE3198C 	 1
    PMANM 	 3.0
    r2034 	 -1.0
    TAURtcx 	 0.0
    MCOATA 	 -1.0
    RE1901R 	 1
    RETH1e 	 1
    r0130 	 3.0
    RE2722G 	 -1.0
    EX_gsn_LPAREN_e_RPAREN_ 	 3
    BTNPLm 	 -1.0
    r2245 	 -1.0
    RE2972M 	 0.0
    r1445 	 3.0
    RE0569E 	 1
    EX_ebastineoh_LPAREN_e_RPAREN_ 	 1
    EX_3ohacmp_LPAREN_e_RPAREN_ 	 1
    25HVITD2tm 	 1
    3MOXTYROX 	 3.0
    CORE4GTg 	 0.0
    DNDPt36m 	 1
    r1156 	 0
    PTVSTM13itr 	 1
    RE3334X 	 3.0
    LPSe 	 -1.0
    EX_gltcho_LPAREN_e_RPAREN_ 	 1
    RE1526X 	 3.0
    DM_n5m2masn_g_ 	 1
    PPA2 	 1
    HXANtl 	 -1.0
    r2223 	 -1.0
    HPDCACRNCPT1 	 0
    FAOXC2252053m 	 2.0
    ACGAMtly 	 1
    S3TASE2ly 	 1
    PRGNLONEtm 	 1
    EX_cgly_LPAREN_e_RPAREN_ 	 3
    15DMThr 	 1
    RE3020R 	 -1.0
    NACHEXA19ly 	 -1.0
    RE2250C 	 1
    TMDOATthc 	 -1.0
    H2MTer_U 	 1
    RE3556C 	 1
    r0764 	 -1.0
    M1316Mg 	 3.0
    HPYRR2x 	 3.0
    RE3134R 	 0.0
    ST8SIA55g 	 -1.0
    S6TASE5ly 	 -1.0
    RTOTALCRNt 	 -1.0
    r2221 	 -1.0
    NACHEXA11ly 	 -1.0
    UDPGALt2g 	 1
    RE3532R 	 3.0
    EX_fvstetglu_LPAREN_e_RPAREN_ 	 1
    GUAPRT 	 0.0
    AM1CSAhr 	 3.0
    EX_2hatvlac_LPAREN_e_RPAREN_ 	 1
    r1774 	 -1.0
    GLYB0AT3tc 	 -1.0
    PECGCHLe 	 1
    RE1819X 	 0
    RE3554M 	 1.0
    RE3447M 	 3.0
    LYStm 	 -1.0
    PTVSTLAChc 	 1
    URATEtx 	 1
    LSTNM7TDhc 	 1
    DM_sprm_c_ 	 1
    r1782 	 -1.0
    EX_am4n9cs_LPAREN_e_RPAREN_ 	 1
    AACTtm 	 1
    EX_ptdca_LPAREN_e_RPAREN_ 	 1
    r2237 	 -1.0
    Q10H2e 	 1
    PROGLYPRO1c 	 -1.0
    5OHFVSGLUhc 	 1
    NACHEX15ly 	 1.0
    AGPAT1 	 2.0
    RETI1 	 1
    GASNASE2ly 	 1.0
    PAPSitr 	 1
    RE3289C 	 1
    EX_lnlc_LPAREN_e_RPAREN_ 	 3
    RN0002R 	 1
    Tyr_ggnt 	 1
    GLCNACT2g 	 -1.0
    DPMVDc 	 -1.0
    GLNSERNaEx 	 -1.0
    r1830 	 -1.0
    ALOX5 	 2.0
    RE1525X 	 3.0
    THRGLYexR 	 -1.0
    GALASE19ly 	 0
    RE2067C 	 1
    UMPK5 	 3.0
    RE3168C 	 1
    COAtr 	 1
    RN0002N 	 1
    r0766 	 -1.0
    M16N6Tg 	 1.0
    RE3033N 	 2.0
    LPS2e 	 -1.0
    DM_T_antigen_g_ 	 1
    NACHEX13ly 	 1.0
    SUCCTD 	 1
    NAGA2ly 	 -1.0
    EX_nad_LPAREN_e_RPAREN_ 	 1
    RE3169R 	 1
    ALAtN1 	 -1.0
    r2065 	 -1.0
    RE2993X 	 -1.0
    r1811 	 -1.0
    r2335 	 3.0
    Asn_X_Ser_Thrtr 	 1
    r0734 	 -1.0
    ESTRGLCABCCte 	 1.0
    AG13T15g 	 1.0
    RE2915X 	 3.0
    r0968 	 1
    RE1917C 	 1
    S6TASE8ly 	 -1.0
    ACETONEt2 	 1.0
    PGISr 	 0.0
    EX_lstnm1_LPAREN_e_RPAREN_ 	 1
    AHANDROSTANGLCtr 	 1
    r2535 	 -1.0
    CYTDtm 	 -1.0
    ALLOP1tu 	 0.0
    FAOXC183C163Gm 	 3.0
    HSAT4ly 	 1
    RE2849C 	 1
    ODECOAtx 	 1
    DESAT18_10 	 -1.0
    GCHOLAt 	 -1.0
    EX_glypro_LPAREN_e_RPAREN_ 	 1
    HSD17B8r 	 -1.0
    r0440 	 1
    MAN1PT2 	 0
    PYDX5Ptm 	 1
    DHAPAx 	 1.0
    BAAT1x 	 -1.0
    RE2306C 	 -1.0
    1531TACRhr 	 1
    MTHFD 	 2.0
    EX_caribup_s_LPAREN_e_RPAREN_ 	 1
    RE2888E 	 -1.0
    FORMCOAtx 	 1
    RE3636C 	 1
    r1569 	 3.0
    PROSTGE1t 	 3.0
    TALA 	 3.0
    HDCEAtr 	 -1.0
    ESTRAABCtc 	 1.0
    FAOXCPRIST2x 	 -1.0
    RE3235R 	 1
    PSSA1_hs 	 3.0
    FAOXC102C101x 	 -1.0
    5OXPROt 	 1
    RE3571R 	 2.0
    PPNCL3 	 -1.0
    r0129 	 3.0
    C16DCc 	 0
    r0360 	 3.0
    RE3175C 	 1
    RE1238X 	 1
    GLCAASE1ly 	 3.0
    TTDCRNt 	 -1.0
    FACOAL2252 	 2.0
    RN0001C 	 1
    q10tm 	 1
    2HATVACIDGLUCteb 	 -1.0
    EX_hco3_LPAREN_e_RPAREN_ 	 1
    r1932 	 -1.0
    ACt2m 	 1
    ACGPID 	 -1.0
    r0555 	 -1.0
    EX_lipoate_LPAREN_e_RPAREN_ 	 1
    RE1819C 	 1.0
    RE3310R 	 1
    EX_ibupgluc_LPAREN_e_RPAREN_ 	 1
    FMNALKPle 	 1.0
    RNDR1 	 2.0
    NRPPHRt4_2_r 	 -1.0
    r2186 	 -1.0
    AM9CSAteb 	 1
    FAOXC2242046m 	 3.0
    PTE2x 	 0
    3HLVSTitr 	 1
    PNP 	 2.0
    RE2868C 	 1
    GAPD 	 3.0
    FUT96g 	 -1.0
    r2359 	 0
    r2511 	 1
    RE3172C 	 1
    6HTSTSTERONEte 	 1
    RE1918C 	 1
    r1877 	 -1.0
    A4GNT2g 	 -1.0
    RE0452M 	 1
    PSYTDECHe 	 1
    GLCNACT3g 	 -1.0
    RE1310C 	 1
    S26Tg 	 0.0
    LYSOXp 	 1
    CRGLZitr 	 1
    r2341 	 0
    4HMDZGLUChr 	 0
    SELNPS 	 3.0
    FAOXC102_4Z_7Zm 	 1.0
    ACITL 	 3.0
    PSHSABCtc 	 3.0
    ATVLACh2r 	 1
    GLAl 	 2.0
    MHGLZitr 	 1
    GMPR 	 -1.0
    SIAASE4ly 	 -1.0
    S3MEACMPtev 	 1
    PHEATB0tc 	 1.0
    G14T15g 	 1.0
    6THFtl 	 1
    MI4PP 	 3.0
    FAOXC140120m 	 3.0
    PDX5PO 	 -1.0
    VD3 	 1
    LSTNM2tev 	 1
    AM9CSAtep 	 1
    DGK2m 	 1
    34DHOXPEGOX 	 3.0
    RE2265C 	 1
    NTD10 	 3.0
    RE3381L 	 3.0
    25HVITD2tin 	 1
    RE3560M 	 1.0
    GALGALGALTHCRMtg 	 1
    10FTHF7GLUtm 	 1
    EX_galgalgalthcrm_hs_LPAREN_e_RPAREN_ 	 1
    13DAMPPOX 	 0.0
    RE3344X 	 -1.0
    GTACMPtev 	 1
    RE3177M 	 3.0
    GCHOLAt3 	 1.0
    13DMThr 	 1
    FAOXC200180m 	 3.0
    6HMSMVhep 	 3.0
    RE3110C 	 3.0
    FVSTETGLUhc 	 1
    M13N2Tg 	 -1.0
    4HATVLACOXDhc 	 3.0
    RN0023X 	 0.0
    RE3260R 	 0
    r1650 	 3.0
    PI45P3K 	 1.0
    FAOXC161_9Em 	 3.0
    RE3166R 	 1
    RE3526X 	 1
    EX_6bhglz_LPAREN_e_RPAREN_ 	 1
    RE2513E 	 -1.0
    FAOXC160140m 	 3.0
    r0924 	 -1.0
    DOLPMT_U 	 1
    r0735 	 -1.0
    SMPD3l 	 0.0
    r1652 	 3.0
    PIACGT 	 -1.0
    RE2859C 	 1
    MTHFDm 	 3.0
    r2277 	 -1.0
    EX_ctp_LPAREN_e_RPAREN_ 	 1
    DM_gpi_sig_er_ 	 1
    GLUNm 	 2.0
    GALASE3ly 	 0
    EX_tacr_LPAREN_e_RPAREN_ 	 1
    PTCRTD 	 1
    GCCam 	 -1.0
    RAI2 	 1
    RE3423C 	 1
    RE2717L 	 -1.0
    ADRNCOAtx 	 1
    RE1933C 	 1
    RE3434C 	 1
    r0988 	 1
    LVACLAChep 	 1
    COUCOAFm 	 1
    GALtly 	 1
    RE2068C 	 1
    L_LACtcm 	 -1.0
    r1750 	 -1.0
    RSVtev 	 1
    FAOXC3DC 	 0
    O2Stn 	 1
    HISTAtu 	 2.0
    RE0928R 	 0
    INSTt4 	 3.0
    S6TASE21ly 	 3.0
    CLOXAtex2 	 -1.0
    NTP3e 	 3.0
    DNDPt51m 	 -1.0
    EX_thf_LPAREN_e_RPAREN_ 	 1
    EX_ura_LPAREN_e_RPAREN_ 	 3
    RE3518R 	 2.0
    EX_glygn4_LPAREN_e_RPAREN_ 	 1
    ACN13ACNGALGBSIDEte 	 1
    HPYRtp 	 1
    GGH_7DHFl 	 1.0
    EX_retinol_9_cis_LPAREN_e_RPAREN_ 	 1
    RE2050C 	 1
    SERHL 	 -1.0
    r2254 	 -1.0
    FAOXC12DCx 	 3.0
    GD1Cte 	 1
    DNDPt21m 	 -1.0
    DALAt2rL 	 -1.0
    FOLTle 	 -1.0
    RE1063C 	 1
    sink_pre_prot_LPAREN_r_RPAREN_ 	 1
    DM_5hpet_LPAREN_r_RPAREN_ 	 1
    CPPPGO 	 -1.0
    AKGDm 	 -1.0
    C121CPT1 	 0
    GLYt2rL 	 -1.0
    LAPCOAl 	 2.0
    ETOHtx 	 1
    PPD2CSPp 	 1
    r0142 	 2.0
    r1443 	 3.0
    DNADtn 	 1
    C102CPT1 	 0
    RETFAt 	 1
    FPGS9 	 -1.0
    RBFK 	 -1.0
    SULFOX 	 2.0
    r1708 	 -1.0
    CSNATm 	 0.0
    NAIt 	 -1.0
    RE3444C 	 0
    FAOXC5030m 	 2.0
    STS2r 	 -1.0
    PMEVKx 	 -1.0
    r1906 	 -1.0
    FUCFUC132GALACGLCGAL14ACGLCGALGLUSIDEtg 	 1
    S6TASE15ly 	 3.0
    UDPGNP 	 1
    EX_estriolglc_LPAREN_e_RPAREN_ 	 1
    r0438 	 1
    r1957 	 -1.0
    RIBt2 	 1
    MAN2_6B1er 	 1
    ACS2 	 1.0
    GLCSGLT1le 	 0.0
    RETI2 	 1
    EX_glygly_LPAREN_e_RPAREN_ 	 1
    RE1816C 	 1
    RE1941C 	 1
    ACOAO7p 	 0.0
    NAt5 	 -1.0
    RE2146R 	 0
    5OHFVSGLUitr 	 1
    DNDPt44m 	 -1.0
    RE2869C 	 1
    7HPVSteb 	 0.0
    r2006 	 -1.0
    MI1P_Dtn 	 1
    DPGM 	 3.0
    r0772 	 -1.0
    RE1062M 	 3.0
    r0224 	 -1.0
    EX_oxa_LPAREN_e_RPAREN_ 	 3
    EX_sfcys_LPAREN_e_RPAREN_ 	 1
    GLCATg 	 -1.0
    RE1691M 	 3.0
    RE3395M 	 1
    r1673 	 -1.0
    r0817 	 1
    PTDCAt 	 1
    MDHx 	 1
    r1792 	 -1.0
    GLCter 	 1
    GLCt2_2 	 0.0
    G14T13g 	 1.0
    AGLPET 	 1
    RE3488N 	 2.0
    RE3559X 	 3.0
    ALCD21_L 	 3.0
    r2438 	 -1.0
    H4ET3er 	 -1.0
    RE3470C 	 1.0
    r0443 	 1
    RE2921X 	 -1.0
    r1454 	 1.0
    ALAGLYexR 	 -1.0
    RE2439C 	 1
    r0616 	 -1.0
    r0775 	 2.0
    H3MTer_U 	 1
    RE3190M 	 3.0
    SO4HCOtex 	 2.0
    r1665 	 -1.0
    RE2476C 	 1
    SO4t4_2 	 -1.0
    r1004 	 1
    EX_omeprazole_LPAREN_e_RPAREN_ 	 1
    EX_glyc_S_LPAREN_e_RPAREN_ 	 1
    PIPLC 	 3.0
    DMGDHm 	 -1.0
    r0719 	 3.0
    GLCAT2g 	 -1.0
    GLCURter 	 1
    4NPHte 	 1
    34DHPHAMT 	 3.0
    RE2522C 	 1
    r1777 	 -1.0
    35DHPVSitr 	 1
    3DPHBH2 	 1
    r1926 	 -1.0
    r0009 	 3.0
    SPODMe 	 -1.0
    RE0919R 	 0
    EX_xolest2_hs_LPAREN_e_RPAREN_ 	 1
    LNLCCPT1 	 0
    EX_udpgal_LPAREN_e_RPAREN_ 	 1
    S6T7g 	 -1.0
    EX_gbside_hs_LPAREN_e_RPAREN_ 	 1
    PGMT 	 3.0
    HEXCOAACBP 	 3.0
    r1590 	 3.0
    RE3488R 	 2.0
    r2068 	 -1.0
    GCCcm 	 -1.0
    GLYtp 	 1
    RETNCOA 	 1
    r0801 	 3.0
    RE3525N 	 1
    VALtec 	 -1.0
    ST8SIA51g 	 -1.0
    RETABCtc 	 -1.0
    r1365 	 1
    PTVSTM3hc 	 3.0
    PHEB0AT3tc 	 -1.0
    UDPGLDCg 	 1.0
    r1742 	 -1.0
    r1052 	 1
    RE2051G 	 3.0
    r2333 	 3.0
    GLYSARPEPT1tc 	 -1.0
    r2109 	 1
    DHPM1 	 0.0
    PMTFATPtc 	 1
    CYSACMPitr 	 1
    LCADi_Dm 	 3.0
    RE3486C 	 1
    r1655 	 3.0
    ARGt4 	 1.0
    FAOXC81C61x 	 -1.0
    B3GNT314g 	 -1.0
    C6DCCACT 	 -1.0
    DM_fol 	 1
    NADN 	 1
    B3GNT31g 	 -1.0
    H2Otp 	 1
    r2192 	 -1.0
    r1736 	 -1.0
    TSTSTERONESte 	 1
    RE3220C 	 -1.0
    r1386 	 1
    r0463 	 3.0
    SIAASE2ly 	 -1.0
    DM_ethamp_r_ 	 1
    r0937 	 1
    PYNP2r 	 3.0
    FTHFDH 	 -1.0
    MALTt1r 	 -1.0
    DCTPtm 	 1
    ATPH1e 	 3.0
    r1067 	 1
    NACHEX1ly 	 1.0
    r1003 	 1
    FAS140COA 	 -1.0
    NS26Tg 	 -1.0
    FA1822ACPH 	 1
    RE3185M 	 2.0
    r0637 	 3.0
    r2019 	 -1.0
    GLYVESSEC 	 -1.0
    r1664 	 -1.0
    AM19CSitr 	 1
    HAS2 	 2.0
    3DSPHR 	 1.0
    TETTET6CPT2 	 -1.0
    ENGASE 	 3.0
    GLCAASE9ly 	 3.0
    GCALDDm 	 3.0
    CYTK13n 	 3.0
    PCHOL_HSter 	 1
    PAPStg 	 0
    r1944 	 -1.0
    TRDR 	 0.0
    NACHEX11ly 	 1.0
    H2O2tn 	 1
    P45027A11m 	 -1.0
    CO2t 	 1
    r2031 	 -1.0
    ASPNATm 	 1
    LINKDEG2ly 	 1
    METB0AT3tc 	 -1.0
    PVSGLUCtev 	 1
    r2169 	 -1.0
    r0724 	 3.0
    r1382 	 -1.0
    r2049 	 -1.0
    T2M26DCOAHLx 	 2.0
    biomass_reaction 	 3
    3MLDAt 	 1
    EX_atp_LPAREN_e_RPAREN_ 	 3
    RE3164R 	 1
    r2113 	 0.0
    RE1804C 	 1
    S6TASE12ly 	 3.0
    GULN3D 	 1
    EX_retinal_LPAREN_e_RPAREN_ 	 1
    NADPN 	 1
    ACNAM9PL 	 2.0
    LSTNM7thc 	 1
    FAOXC225_4Z_7Z_10Z_13Z_16Zx 	 -1.0
    C6DCe 	 1
    SMVitr 	 1
    DNDPt35m 	 1
    FACOAL80i 	 1.0
    r2332 	 3.0
    r2106 	 1
    SOAT12r 	 1
    r1517 	 3.0
    EX_crvs_LPAREN_e_RPAREN_ 	 1
    EX_triodthy_LPAREN_e_RPAREN_ 	 1
    EX_Tyr_ggn_LPAREN_e_RPAREN_ 	 1
    r1609 	 3.0
    EX_estroneglc_LPAREN_e_RPAREN_ 	 1
    S3T1g 	 -1.0
    r2373 	 -1.0
    r0575 	 3.0
    r2098 	 0.0
    r1824 	 -1.0
    CHOLtg 	 1
    ACGAMPM 	 -1.0
    HSD17B7r 	 2.0
    EX_duri_LPAREN_e_RPAREN_ 	 1
    r1159 	 1
    DM_pmtcoa_LPAREN_r_RPAREN_ 	 1
    ALAPAT4te 	 -1.0
    TMNDNCCPT2 	 -1.0
    NADNe 	 -1.0
    3HSMVitr 	 1
    10FTHFtm 	 1
    r1325 	 1
    AHC 	 2.0
    IPDPtr 	 1
    ASPDTDe 	 1
    r2141 	 -1.0
    r1684 	 -1.0
    RIBFLVt3 	 -1.0
    RE3586X 	 -1.0
    S2TASE4ly 	 1
    RE3434R 	 2.0
    34HPPte 	 1
    r0422 	 2.0
    3HBCDm 	 3.0
    EX_2hb_LPAREN_e_RPAREN_ 	 3
    SERGLNNaEx 	 -1.0
    KSItly 	 1
    RE2452C 	 1
    DNDPt8m 	 -1.0
    DLNLCGCPT1 	 0
    P45017A4r 	 -1.0
    FAOXC221201x 	 0.0
    r2375 	 -1.0
    EX_acetone_LPAREN_e_RPAREN_ 	 3
    CYTDt5le 	 0.0
    O16G1e 	 -1.0
    RE0923C 	 -1.0
    r1976 	 -1.0
    6MELVACtbc 	 -1.0
    TMDitr 	 1
    RE3409C 	 1
    r0425 	 0.0
    7BHGLZitr 	 1
    ADSELK 	 3.0
    r1017 	 1
    HSD11B2r 	 1.0
    CO2tm 	 1
    FA180ACPH 	 -1.0
    EX_eicostet_LPAREN_e_RPAREN_ 	 1
    RAtn3 	 1
    RE1826C 	 1
    r1448 	 3.0
    ESTRADIOLGLCtr 	 1
    S6TASE6ly 	 -1.0
    EX_ile_L_LPAREN_e_RPAREN_ 	 3
    ICDHxm 	 0.0
    4OHMDZtev 	 1
    VITD2tm 	 1
    DPCOAtl 	 1
    G6PDH1rer 	 -1.0
    RE3307C 	 1.0
    PIPLCn 	 -1.0
    GLCAE1g 	 1
    ABUTD 	 3.0
    RE3308X 	 1
    HDECEACBP 	 3.0
    CHOLATEt3 	 -1.0
    C8DCe 	 1
    TSULt4_3 	 -1.0
    TYMSFt 	 1
    TDCHOLAte 	 1
    RE3113C 	 1
    r0716 	 -1.0
    r0475 	 -1.0
    EX_ncam_LPAREN_e_RPAREN_ 	 3
    RE1827M 	 -1.0
    GLYOXm 	 -1.0
    RE1522X 	 -1.0
    r1696 	 -1.0
    RE2522X 	 1.0
    RE3343C 	 0
    FAOXC8060m 	 2.0
    r2030 	 -1.0
    1HMDGLUChr 	 0
    GT1Ate 	 1
    ADPACDAc 	 1
    PTHPS 	 2.0
    RE0452N 	 1
    r2158 	 -1.0
    EX_tchola_LPAREN_e_RPAREN_ 	 1
    GRTTx 	 2.0
    RE2625C 	 1
    STRDNCCPT2 	 -1.0
    AG13T18g 	 1.0
    PEFLIP 	 2.0
    RE2443C 	 1
    RE2649M 	 1.0
    COAtg 	 1
    EX_hdca_LPAREN_e_RPAREN_ 	 3
    C204CRNt 	 -1.0
    DHAPA 	 1.0
    6BHGLZhr 	 0.0
    GALACGLCGALGBSIDEtg 	 1
    RE0688X 	 -1.0
    RE1957G 	 -1.0
    ACRNtm 	 -1.0
    GALt2_2 	 0.0
    DMATT 	 -1.0
    ORNt4m 	 -1.0
    r2139 	 -1.0
    CYOOm2 	 0
    RE2407C 	 1
    LCYSTATm 	 3.0
    SERCYSNaEx 	 3.0
    RE2041C 	 1
    GBSIDEte 	 1
    RN0022C 	 0.0
    DM_ncam 	 1
    ABUTt2rL 	 -1.0
    GALFUC12GAL14ACGLCGALGLUSIDEte 	 1
    EX_dchac_LPAREN_e_RPAREN_ 	 1
    r0717 	 -1.0
    r2241 	 -1.0
    AKR1D 	 -1.0
    EX_tethex3_LPAREN_e_RPAREN_ 	 1
    RE2453M 	 3.0
    THMP 	 1
    NOS1 	 -1.0
    RE2155R 	 0
    RE3248M 	 3.0
    RE2387R 	 -1.0
    r1672 	 -1.0
    MANtly 	 1
    r2347 	 0
    CORE6GTg 	 -1.0
    UREAt5 	 0.0
    CBR1 	 1.0
    r1493 	 1
    r1725 	 -1.0
    NTD4l 	 2.0
    LEUKTRA4tr 	 1
    LST4EXPTDhc 	 1
    NADKm 	 3.0
    r2514 	 1
    RE2649C 	 0
    UDPACGALtl 	 1
    PGESr 	 -1.0
    DNDPt4m 	 -1.0
    MTHFC 	 2.0
    PNTEHe 	 1
    TSACMGLUChr 	 1
    NACHEXA2ly 	 -1.0
    PNTEH 	 -1.0
    EX_thmtp_LPAREN_e_RPAREN_ 	 1
    GP1Cte 	 1
    CSPG_Dt 	 1
    r2153 	 -1.0
    EX_C02470_LPAREN_e_RPAREN_ 	 1
    B3GNT310g 	 -1.0
    r0788 	 1
    EPCTX 	 -1.0
    ASCBt4 	 -1.0
    RE3576X 	 3.0
    GDPtg 	 1
    RE3485N 	 2.0
    r0522 	 -1.0
    DM_13_cis_retn_n_ 	 1
    PYDXtr 	 1
    BILIRUBtr 	 1
    EX_xyl_D_LPAREN_e_RPAREN_ 	 3
    RE2799C 	 1
    r1583 	 3.0
    EX_pheacgln_LPAREN_e_RPAREN_ 	 1
    RE3335X 	 2.0
    GGH_5THFl 	 1.0
    ALASm 	 1.0
    OCDCEAFATPtc 	 -1.0
    DNDPt54m 	 -1.0
    NMNATr 	 -1.0
    NDPK6n 	 0
    r0757 	 2.0
    AMANK 	 2.0
    FAOXC182C162m 	 3.0
    THYOXt2 	 -1.0
    RE2069C 	 1
    P45027A12m 	 -1.0
    r1788 	 -1.0
    P4502F1 	 -1.0
    RE2129C 	 1
    EX_HC02205_LPAREN_e_RPAREN_ 	 1
    AGLPT 	 1
    RE2333C 	 1
    PGLYCt 	 1
    r1841 	 -1.0
    P4502C94 	 -1.0
    r1707 	 -1.0
    r2189 	 -1.0
    GGH_10FTHF7GLUe 	 1.0
    BTNt3i 	 -1.0
    GLYCTDle 	 1
    r1639 	 3.0
    MALTe 	 -1.0
    RE3577X 	 2.0
    RE3370C 	 1
    EX_leuktrB4_LPAREN_e_RPAREN_ 	 1
    LNELDCCRNt 	 -1.0
    FAOXC121C101m 	 0.0
    H6MTer_U 	 1
    EX_ala_D_LPAREN_e_RPAREN_ 	 3
    EX_csasulp_LPAREN_e_RPAREN_ 	 1
    EX_bilirub_LPAREN_e_RPAREN_ 	 1
    RE2972G 	 0.0
    EX_am1alcs_LPAREN_e_RPAREN_ 	 1
    ME2m 	 -1.0
    r1892 	 -1.0
    r1628 	 3.0
    ETOHt 	 1
    FAOXC101_4Zm 	 1.0
    RE3012M 	 3.0
    RE3411C 	 1
    RE3370R 	 -1.0
    FACOAL240 	 2.0
    HISDC 	 -1.0
    EX_c10crn_ 	 1
    C5DCe 	 1
    L_LACtm 	 0.0
    r1640 	 3.0
    r1456 	 1
    FAOXC1811601m 	 2.0
    C4CRNe 	 1
    4HATVLACitr 	 1
    RE3122R 	 1
    EX_debrisoquine_LPAREN_e_RPAREN_ 	 1
    PPItx 	 1
    RE3226C 	 3.0
    GUMDCHAe 	 1
    r0163 	 -1.0
    r0301 	 3.0
    PHEMEtm 	 1
    FAS100COA 	 -1.0
    RE1310M 	 -1.0
    DNDPt41m 	 -1.0
    ACCOAC 	 1.0
    RE3160R 	 2.0
    r2128 	 1
    r0950 	 1
    UTPtn 	 1
    CYTK10n 	 3.0
    EX_1glyc_hs_LPAREN_e_RPAREN_ 	 1
    FAOXC12DCc 	 0
    RE1135G 	 -1.0
    G12MT2_U 	 1
    FAOXC161_7Em 	 3.0
    r2070 	 -1.0
    r1737 	 -1.0
    FAOXC143_5Z_8Z_11Zm 	 3.0
    RE3120C 	 3.0
    DNDPt40m 	 -1.0
    CLCFTRte 	 0.0
    EX_4abutn_LPAREN_e_RPAREN_ 	 1
    r2434 	 -1.0
    EX_galfucgalacglcgal14acglcgalgluside_hs_LPAREN_e_RPAREN_ 	 1
    NRVNCt 	 1
    ACHtn 	 1
    r0580 	 1
    r2088 	 0.0
    RE2898C 	 1
    r2252 	 -1.0
    GLYSNAT5tc 	 -1.0
    ONPTHLte 	 1
    6HSMVitr 	 1
    C181CPT2 	 -1.0
    r0932 	 1
    GALC 	 -1.0
    r1178 	 1.0
    P45011B12m 	 -1.0
    5HOMEPRAZOLEte 	 1
    C181OHc 	 0
    RE1809C 	 1
    TETPENT3CRNt 	 -1.0
    RE0918R 	 -1.0
    RE3163R 	 1
    r1781 	 -1.0
    RE3121C 	 1
    CDSm 	 3.0
    r0630 	 -1.0
    GD1Ctg 	 1
    3MOPt2im 	 1
    NACHEX24ly 	 1.0
    S2T2g 	 1.0
    NDP7ex 	 3.0
    r1731 	 -1.0
    RE0066C 	 1
    GALSIDEtl 	 1
    5OHFVShc 	 3.0
    C4DCe 	 1
    r2532 	 -1.0
    RE1811C 	 1
    EX_elaid_LPAREN_e_RPAREN_ 	 1
    EX_thymd_LPAREN_e_RPAREN_ 	 3
    r1744 	 -1.0
    CYTD 	 1.0
    r0024 	 1
    RE3218C 	 -1.0
    RE1957R 	 -1.0
    ATPH2e 	 3.0
    DALAOXx 	 -1.0
    r1028 	 -1.0
    VITD3Hm 	 -1.0
    DCSPTN1CPT1 	 0
    r0680 	 3.0
    SBTD_D2 	 3.0
    FPGS7m 	 -1.0
    LSTNM4itr 	 1
    r1167 	 -1.0
    RE2860C 	 1
    HIVCACBP 	 3.0
    CRNATBtc 	 1.0
    ARABR 	 3.0
    EX_3aib_LPAREN_e_RPAREN_ 	 1
    RE3273G 	 3.0
    G3M8MASNterg 	 1
    CHTNASEe 	 -1.0
    G6PDH2r 	 -1.0
    r0796 	 1
    ALLOPOXDhep 	 -1.0
    DAGKn_hs 	 -1.0
    T4HCINNOX 	 -1.0
    HPYRRy 	 2.0
    PI34P3Pn 	 1
    r0708 	 2.0
    RE2624M 	 -1.0
    CO2tg 	 1
    EX_pectingchol_LPAREN_e_RPAREN_ 	 1
    AM1CCShr 	 1
    RPI 	 -1.0
    LACZe 	 0
    4HATVACIDOXDhc 	 3.0
    FAOXC162_7E_10Em 	 3.0
    FERO 	 0.0
    KYNATESYN 	 1
    RN0029C 	 1
    RE2523C 	 1
    RE3340C 	 3.0
    EX_bilglcur_LPAREN_e_RPAREN_ 	 1
    GGNG 	 1.0
    RE2638X 	 -1.0
    r1175 	 1.0
    UMPKm 	 1
    MI14Ptn 	 1
    4BHGLZtev 	 1
    RE3158X 	 3.0
    THMMPt4 	 -1.0
    r1896 	 -1.0
    FAEL184 	 3.0
    r2285 	 -1.0
    RE1916C 	 1
    EX_ac_LPAREN_e_RPAREN_ 	 3
    AG13T16g 	 1.0
    ARGtiDF 	 3.0
    r1930 	 -1.0
    CAATPS 	 3.0
    r1908 	 -1.0
    DNDPt61m 	 -1.0
    ALADGLNexR 	 -1.0
    NADPNe 	 -1.0
    r2510 	 1
    G6PDA 	 0.0
    r0444 	 1
    11DOCRTSLtm 	 1
    S3T3g 	 2.0
    NACHEXA9ly 	 -1.0
    DNDPt27m 	 1
    HMGCOARr 	 3.0
    TAGHSTDe 	 1
    EX_fald_LPAREN_e_RPAREN_ 	 3
    TCHOLABCtc 	 1.0
    RE3066X 	 -1.0
    RE2304E 	 0.0
    RE2440C 	 1
    DATPtn 	 1
    TYR3MO2 	 -1.0
    PAPitr 	 1
    FAOXC2252053x 	 0.0
    r1400 	 1
    GTHRDt 	 1
    LCADi_D 	 3.0
    RE3265C 	 -1.0
    ABO6g 	 -1.0
    r0549 	 3.0
    MMSAD1m 	 -1.0
    AMPTASECGe 	 0.0
    D_LACt2 	 1.0
    HISyLATtc 	 -1.0
    EX_h2o_LPAREN_e_RPAREN_ 	 3
    r1931 	 -1.0
    XYLUR 	 3.0
    r0960 	 -1.0
    AICARte 	 1
    GABABGTtc 	 -1.0
    DAG_HSter 	 1
    XAO2x 	 -1.0
    GCC2bim 	 -1.0
    GLCAE2g 	 -1.0
    r1074 	 1
    r1967 	 -1.0
    AM1ACCStev 	 1
    r1897 	 -1.0
    SERGLYexR 	 -1.0
    HSD3B11r 	 -1.0
    7BHGLZGLChr 	 1
    RE2677R 	 -1.0
    LALDO 	 3.0
    PHYCBOXL 	 -1.0
    14HMDZitr 	 1
    SERtp 	 1
    r2346 	 0
    r2242 	 -1.0
    r2286 	 -1.0
    r1895 	 -1.0
    EX_2hatvlacgluc_LPAREN_e_RPAREN_ 	 1
    EX_hestratriol_LPAREN_e_RPAREN_ 	 1
    AM1ACCSitr 	 1
    7BHGLZhr 	 0.0
    EX_dca_LPAREN_e_RPAREN_ 	 1
    r2229 	 -1.0
    FUCtly 	 1
    EX_HC02197_LPAREN_e_RPAREN_ 	 1
    CHLP 	 -1.0
    PSSA2_hs 	 -1.0
    C16txc 	 -1.0
    GLYGLYPEPT1tc 	 -1.0
    r1817 	 -1.0
    EX_CE2838_LPAREN_e_RPAREN_ 	 1
    2HATVLACOXDhc 	 3.0
    MM5cg 	 3.0
    r2036 	 -1.0
    DM_core7_g_ 	 1
    EX_HC00955_LPAREN_e_RPAREN_ 	 1
    FUT34g 	 1.0
    NNATn 	 -1.0
    2HATVACIDthc 	 -1.0
    C161CPT22 	 -1.0
    EX_adrn_LPAREN_e_RPAREN_ 	 1
    EX_dheas_LPAREN_e_RPAREN_ 	 1
    AMPtp 	 1
    r1148 	 1
    THRALANaEx 	 3.0
    HSPGtly 	 1
    C4STMO1r 	 3.0
    RE2128C 	 1
    GALGALGALTHCRMte 	 1
    FUCACNGAL14ACGLCGALGLUSIDEte 	 1
    EX_man_LPAREN_e_RPAREN_ 	 3
    EX_dpcoa_LPAREN_e_RPAREN_ 	 1
    RN0021R 	 0.0
    EX_andrstrn_LPAREN_e_RPAREN_ 	 1
    SERDGLYexR 	 -1.0
    r1000 	 1
    DMGtm 	 1
    FAOXOHC16C16DCc 	 3.0
    SMVthep 	 -1.0
    CRVSM23teb 	 0.0
    FAOXC163x 	 1.0
    RE3039C 	 1
    r0688 	 3.0
    LPS3e 	 -1.0
    RE3167C 	 1
    RE2220C 	 1
    r0917 	 -1.0
    XYLtly 	 1
    EX_dag_hs_LPAREN_e_RPAREN_ 	 3
    ARACHt 	 1
    ADPRIBt 	 -1.0
    EX_tauribup_S_LPAREN_e_RPAREN_ 	 1
    LSTNM1hr 	 3.0
    PPPGOm 	 -1.0
    RE3403M 	 3.0
    RE2117M 	 -1.0
    r0653 	 3.0
    MALT 	 -1.0
    AM1ACShr 	 1
    RE3337M 	 3.0
    3HPVShc 	 3.0
    FAOXOHC22C22DCc 	 3.0
    FUCACGALFUCGALACGLCGALGLUSIDEte 	 1
    6HLVSTAChep 	 -1.0
    RE3165R 	 3.0
    EX_6ohfvs_LPAREN_e_RPAREN_ 	 1
    RE3562C 	 -1.0
    FAOXC165C164m 	 3.0
    RE3307M 	 0
    NDPK6m 	 -1.0
    NACHEX10ly 	 1.0
    r2320 	 3.0
    NACHEX3ly 	 1.0
    r2298 	 -1.0
    NACHEXA15ly 	 -1.0
    r0340 	 -1.0
    RE3470M 	 0
    r1611 	 3.0
    RSVATPtu 	 0.0
    r1886 	 -1.0
    RE1077C 	 1
    LEUTAm 	 -1.0
    r1553 	 3.0
    R_group_phosphotase_3 	 1
    EX_phyt_LPAREN_e_RPAREN_ 	 1
    CYTK11n 	 3.0
    r2009 	 -1.0
    LSTNitr 	 1
    NTD8 	 3.0
    CSAitr 	 1
    NRVNCCRNt 	 -1.0
    RE1904R 	 -1.0
    HKYNH 	 1.0
    FPGS8 	 -1.0
    RE3367E 	 -1.0
    INStl 	 -1.0
    ABO1g 	 -1.0
    NMNtn 	 1
    3HPVSTETteb 	 0.0
    RE3485C 	 1
    GALGLUSIDEtg 	 1
    ST6GALNAC28 	 3.0
    RE2429M 	 -1.0
    EX_paf_hs_LPAREN_e_RPAREN_ 	 1
    FAOXC61C4x 	 -1.0
    KSII_CORE2t 	 1
    r0634 	 3.0
    DCK1m 	 1
    RE1538C 	 1
    NDPK3n 	 0
    ORNDC 	 3.0
    XOLTRI25tc 	 1
    FUMm 	 3.0
    EX_atvlac_LPAREN_e_RPAREN_ 	 1
    r1801 	 -1.0
    r1168 	 -1.0
    RE2240C 	 1
    r2408 	 -1.0
    EX_strch1_LPAREN_e_RPAREN_ 	 3
    r1044 	 0.0
    DGTPtm 	 1
    GLCAT6g 	 -1.0
    CAt7r 	 0.0
    EX_6epvs_LPAREN_e_RPAREN_ 	 1
    TCYNTt 	 1
    r0173 	 3.0
    CYTDtl 	 -1.0
    EX_dmantipyrine_LPAREN_e_RPAREN_ 	 1
    AHANDROSTANGLCte 	 1.0
    56EPPVShc 	 3.0
    GLYCLTtp 	 1
    NDERSVitr 	 1
    r2471 	 -1.0
    r2538 	 1
    P45027A13m 	 -1.0
    5HOXINDACTOXm 	 3.0
    RE3407C 	 1
    6OHFVSteb 	 1
    EX_mhglz_LPAREN_e_RPAREN_ 	 1
    HOMt4 	 3.0
    FAOXOHMC10DC10c 	 3.0
    UDPGLCAtg 	 1
    CSPG_Ctly 	 1
    NACHEXA5ly 	 -1.0
    EX_3hlvstacid_LPAREN_e_RPAREN_ 	 1
    AM4NC9CSteb 	 1
    6OHFVSGLUtev 	 1
    NRPPHRSFt 	 1
    2HATVLACtep 	 1
    OCT11EFATP 	 1.0
    S6TASE2ly 	 3.0
    GSNtm 	 -1.0
    r2264 	 -1.0
    RE3224C 	 3.0
    EX_pyr_LPAREN_e_RPAREN_ 	 3
    RE3140X 	 -1.0
    FBP 	 -1.0
    FTCD 	 -1.0
    RE3288R 	 -1.0
    r2052 	 -1.0
    RE3521C 	 1
    ATVLACThc 	 -1.0
    S2L2N2M2MASNtly 	 1
    ALCD2yf 	 3.0
    EX_3aib_D_LPAREN_e_RPAREN_ 	 1
    GALASE20ly 	 0
    r1441 	 1
    PRO_Dtde 	 1
    RE1303C 	 1
    EX_arg_L_LPAREN_e_RPAREN_ 	 3
    ARACHCOAtx 	 1
    B3GNT33g 	 -1.0
    RE3019C 	 1
    r2300 	 -1.0
    FAOXC6C4x 	 -1.0
    PEt 	 1
    GD1B2tg 	 1
    EX_HC02213_LPAREN_e_RPAREN_ 	 1
    SGALSIDEtg 	 1
    EX_malthp_LPAREN_e_RPAREN_ 	 1
    ILEB0AT2tc 	 1
    ACMPShc 	 1.0
    RE3045C 	 1
    RE1815R 	 3.0
    NH4tn 	 1
    RE3448M 	 3.0
    ST6GALNAC26 	 3.0
    S2L2FN2M2MASNtly 	 1
    RE2768R 	 -1.0
    NDPK4n 	 0
    CORE2GTg 	 0.0
    RE2867C 	 1
    BDMT_U 	 1
    EX_thr_L_LPAREN_e_RPAREN_ 	 3
    r2173 	 -1.0
    S6TASE9ly 	 -1.0
    ALAASNNaEx 	 -1.0
    r0784 	 1
    S6T4g 	 -1.0
    RE1796C 	 1
    RE3287C 	 1
    RE2272C 	 1
    FE2DMT1 	 3.0
    P4503A5 	 -1.0
    r1924 	 -1.0
    r1565 	 3.0
    r1562 	 3.0
    EX_succ_LPAREN_e_RPAREN_ 	 1
    AG13T2g 	 1.0
    13DMTitr 	 1
    RE0938C 	 0.0
    FOAXC122C101x 	 -1.0
    RE0690E 	 -1.0
    H2ETer 	 -1.0
    FAOXC121_3Zm 	 1.0
    r0558 	 1
    EX_5oxpro_LPAREN_e_RPAREN_ 	 1
    r1627 	 3.0
    EX_3tetd7ecoacrn_ 	 1
    r1436 	 1
    RE3578X 	 -1.0
    FAOXC164m 	 1.0
    3ISPVSitr 	 1
    r1257 	 2.0
    BVITEt 	 1
    CITtbm 	 -1.0



```python
for item in conf_MeanCancerBiopsy.items():
    print(item[0], "\t",item[1])
```

    r1516 	 3.0
    RE2594C 	 1
    BBHOX 	 1.0
    LVSTACOXD6MEhep 	 1.0
    PMTCOAFABP1tc 	 -1.0
    EX_ptvstm3_LPAREN_e_RPAREN_ 	 1
    RE1829M 	 0.0
    r1866 	 0.0
    G14T7g 	 3.0
    RE3144M 	 1.0
    NACHEX26ly 	 3.0
    r0853 	 1
    r1592 	 3.0
    FA120ACPH 	 -1.0
    MM7B1g 	 3.0
    RE1050C 	 3.0
    ST3GAL21g 	 -1.0
    PPDOx 	 3.0
    HSD3A1r 	 -1.0
    r2316 	 1.0
    r0694 	 0.0
    NDPK3m 	 0.0
    2HBt2 	 2.0
    TMDOATtev 	 1
    EX_sl_L_LPAREN_e_RPAREN_ 	 1
    RE0702C 	 3.0
    EX_trp_L_LPAREN_e_RPAREN_ 	 3
    S6T11g 	 0.0
    r2060 	 -1.0
    RE2972R 	 0.0
    RE1956R 	 2.0
    r2081 	 2.0
    ST3GAL23g 	 -1.0
    r2226 	 -1.0
    MCLACCYSR 	 1
    RE2988X 	 3.0
    r0695 	 0.0
    r0591 	 -1.0
    G14T4g 	 3.0
    EBASTINEtr 	 1
    C120CPT1 	 1.0
    SEBCOAPET 	 -1.0
    r1816 	 0.0
    UDPGD 	 2.0
    RE2991X 	 3.0
    r0246 	 3.0
    r0962 	 1
    RE1448R 	 1
    r2179 	 -1.0
    GCNTg 	 1.0
    FACOAL206 	 0.0
    RE2349C 	 -1.0
    TMDM1hr 	 0.0
    PI4PP 	 1
    S2TASE3ly 	 3.0
    RE3103R 	 3.0
    LVSTitr 	 1
    P4504B1r 	 0.0
    LALDO2x 	 3.0
    r1964 	 0.0
    FAOXC221C201x 	 -1.0
    NDPK5m 	 0.0
    RNDR2 	 3.0
    EX_35dsmv_LPAREN_e_RPAREN_ 	 1
    r2002 	 -1.0
    1HIBUPGLUC_Sthv 	 1
    MM8Cg 	 3.0
    C6CRNtcx 	 -1.0
    EX_fvstet_LPAREN_e_RPAREN_ 	 1
    KHK 	 -1.0
    r1546 	 3.0
    RE1816R 	 3.0
    r0590 	 -1.0
    MALTly 	 0.0
    EX_glygn2_LPAREN_e_RPAREN_ 	 3
    ESTSULT 	 2.0
    CHOLtn 	 1
    AM19CSteb 	 1
    GLNALANaEx 	 0.0
    MI34PP 	 1.0
    RE3534M 	 3.0
    CREATtmdiffir 	 1
    GluForTx 	 -1.0
    r2214 	 -1.0
    PROSTGE2t3 	 1.0
    TLACFVShc 	 1
    7DHFtm 	 1
    ACMPGLUTtep 	 0.0
    P4502C8 	 -1.0
    DSPVStev 	 1
    LNLNCGCPT2 	 0.0
    5HXKYNOXDA 	 3.0
    VALATB0tc 	 3.0
    RE1919C 	 1
    ALAALACNc 	 3.0
    DHCR72r 	 0.0
    FAOXC121_3Em 	 3.0
    EX_psyl_LPAREN_e_RPAREN_ 	 1
    GACMTRc 	 -1.0
    r1576 	 3.0
    RE0565C 	 3.0
    RE3535R 	 1
    CYTK2n 	 3.0
    RE1818C 	 1
    AICART 	 3.0
    ADRNt 	 1
    MAGt 	 1
    LAPCOAe 	 1
    EX_sph1p_LPAREN_e_RPAREN_ 	 3
    r0666 	 3.0
    GALNACT3g 	 1.0
    EX_nh4_LPAREN_e_RPAREN_ 	 3
    FAOXC170150m 	 3.0
    r1477 	 1.0
    RE3038R 	 2.0
    RE3236C 	 1
    CRVSM23hc 	 -1.0
    AG13T12g 	 3.0
    RE1651C 	 1
    EX_HC02161_LPAREN_e_RPAREN_ 	 1
    DOCOSCOAtxc 	 -1.0
    SULPACMPtev 	 1.0
    EX_appnn_LPAREN_e_RPAREN_ 	 1
    BILGLCURte 	 3.0
    EX_dmhptcrn_LPAREN_e_RPAREN_ 	 1
    r2521 	 1
    r0085 	 3.0
    O2tm 	 1
    ENMAN3g 	 1
    ACACT6p 	 3.0
    ADNCNT3tc 	 0.0
    RE3384M 	 3.0
    r2401 	 -1.0
    RE0689C 	 -1.0
    RE3014R 	 0.0
    RE1096M 	 0.0
    RE2651R 	 1
    r1027 	 1
    RE2919M 	 3.0
    BTND1 	 -1.0
    GGH_7THFl 	 3.0
    EX_34dhoxpeg_LPAREN_e_RPAREN_ 	 1
    r0940 	 1
    EX_oxyp_LPAREN_e_RPAREN_ 	 1
    EHGLAT2m 	 3.0
    GP1CALPHAte 	 1
    RE2384C 	 1
    GALASE18ly 	 0
    EX_cyan_LPAREN_e_RPAREN_ 	 1
    RE3243C 	 1
    CRMPte 	 1
    BALAtmr 	 1
    r1401 	 1
    RE3286R 	 0.0
    SO4OXAtex2 	 -1.0
    EX_dttp_LPAREN_e_RPAREN_ 	 1
    SPRMTDe 	 0.0
    r0380 	 0
    EX_i_LPAREN_e_RPAREN_ 	 1
    RIBFLVte 	 0.0
    EX_4nphsf_LPAREN_e_RPAREN_ 	 1
    EX_h2o2_LPAREN_e_RPAREN_ 	 1
    CH25H 	 0.0
    FACOAL191 	 0.0
    NADPtru 	 1
    RE2920X 	 1.0
    AGDC 	 0
    RE3430C 	 2.0
    CRTSLtr 	 1
    r0082 	 3.0
    KYN 	 3.0
    GALSIDEtg 	 1
    r2069 	 -1.0
    RE1254C 	 1
    RE3525C 	 1.0
    EX_alaala_LPAREN_e_RPAREN_ 	 1
    r0537 	 3.0
    EX_56dhpvs_LPAREN_e_RPAREN_ 	 1
    EX_crvsm23_LPAREN_e_RPAREN_ 	 1
    EX_so4_LPAREN_e_RPAREN_ 	 3
    DESAT20_1 	 0
    r1292 	 1
    GLZitr 	 1
    PTRCAT1 	 3.0
    SEBACACT 	 1
    CGLYt3_2_ 	 1.0
    CYTK1n 	 3.0
    ACP1_FMN_ 	 3.0
    RE3238C 	 1
    NTD3 	 3.0
    RE3153R 	 1
    RE3513N 	 2.0
    ACCOAgt 	 3.0
    HESTRATRIOLte 	 1
    NACHEXA8ly 	 3.0
    DECDPtm 	 1
    RE3009C 	 0
    3HPVSTETtev 	 1
    r1752 	 0.0
    fumt 	 1
    EX_ahdt_LPAREN_e_RPAREN_ 	 1
    FAOXC24C22x 	 -1.0
    HS3ly 	 1.0
    RE2974N 	 1
    PUNP2 	 3.0
    r0242 	 -1.0
    r2177 	 -1.0
    MCDm 	 -1.0
    FA182ACPH 	 -1.0
    RE3018C 	 1
    GALASE12ly 	 0
    56DHPVSitr 	 1
    5ADTSTSTERONESte 	 1
    C180CPT1 	 1.0
    DOLMANP_Lter 	 1
    GLZtd 	 1
    r1728 	 0.0
    ASPPROASCT1 	 0.0
    r1012 	 1
    r1799 	 0.0
    ST8SIA53g 	 -1.0
    TTDCPT1 	 1.0
    AMACR2p 	 -1.0
    SPODMn 	 3.0
    FUMAC 	 0.0
    THRATB0tc 	 3.0
    AM1C9CSitr 	 1
    ORETNF2 	 1
    r1687 	 0.0
    APNNOXte 	 1
    biomass_carbohydrate 	 3
    EX_crvsm1_LPAREN_e_RPAREN_ 	 1
    EX_ins_LPAREN_e_RPAREN_ 	 1
    r2015 	 -1.0
    r2066 	 -1.0
    INSKm 	 1
    CYTDt 	 -1.0
    MALSO4tm 	 -1.0
    10FTHF5GLUtl 	 1
    r1613 	 0.0
    r1487 	 3.0
    RTOTALt 	 1
    RE3261C 	 -1.0
    DIGALSIDEtl 	 1
    AMACRr 	 -1.0
    FATP8t 	 0.0
    r1614 	 0.0
    DESAT18_7 	 0.0
    r0604 	 3.0
    RE3112R 	 1
    CRTNsyn 	 1
    RE3264R 	 0
    r1450 	 0.0
    TAURt4_2_r 	 3.0
    G14T10g 	 3.0
    XANDp 	 2.0
    RE3464R 	 0.0
    RE3110R 	 3.0
    TSTSTERONEt 	 1
    r1661 	 3.0
    RE3251M 	 0.0
    C2M26DCOAHLm 	 3.0
    RE1925C 	 1
    RE3245C 	 1
    r1666 	 0.0
    EX_13dmt_LPAREN_e_RPAREN_ 	 1
    CYSTGL 	 0.0
    SMVACIDitr 	 1
    GALGLUSIDEtl 	 1
    FALDtly 	 1
    PACCOAL 	 1
    C4x 	 1.0
    EX_acac_LPAREN_e_RPAREN_ 	 3
    RE0918G 	 0.0
    r0791 	 0.0
    GLY3Pt 	 1
    DATPtm 	 1
    RE0937E 	 0.0
    CBASPte 	 1
    r0915 	 0.0
    EX_urea_LPAREN_e_RPAREN_ 	 3
    CYTK6 	 3.0
    r1530 	 1.0
    PGM 	 3.0
    LPS3 	 2.0
    NAGLCAly 	 1
    EX_pheme_LPAREN_e_RPAREN_ 	 3
    r1726 	 0.0
    6BHGLZitr 	 1
    EX_prostgf2_LPAREN_e_RPAREN_ 	 1
    FBA4 	 3.0
    BALAVECSEC 	 -1.0
    RE1099G 	 0.0
    HSD17B2r 	 1.0
    r1568 	 3.0
    EICOSTETCRNt 	 -1.0
    NACHEXA13ly 	 3.0
    r2344 	 0
    VALt5m 	 1
    r2243 	 -1.0
    r1451 	 0.0
    r1880 	 0.0
    ACGAGBSIDEtl 	 1
    RE0864C 	 1
    DMNONCRNCPT2 	 1
    r2032 	 -1.0
    EX_no_LPAREN_e_RPAREN_ 	 1
    r0765 	 0.0
    FUCASEe 	 2.0
    KSII_CORE4tly 	 1
    DCSPTN1COAtxc 	 -1.0
    PGL 	 3.0
    DHEAStr 	 1
    sink_citr_LPAREN_c_RPAREN_ 	 1
    LSTNM2hr 	 1.0
    RE1904C 	 1
    DUDPtn 	 1
    OCTDECE1CRNe 	 1
    FAOXC143C123m 	 3.0
    r1020 	 1
    RE2248C 	 1
    ADNK1m 	 1
    TDCRNe 	 1
    RE3460C 	 1
    r2219 	 -1.0
    HIVCRNe 	 1
    GPAMm_hs 	 0.0
    RE1630R 	 1
    RE1915C 	 1
    RE3537C 	 1
    r0652 	 -1.0
    LACZly 	 0
    RE2151R 	 0
    NACHEX7ly 	 3.0
    NTD2 	 3.0
    FAOXC163GC142m 	 3.0
    RE2541E 	 0.0
    ABO4g 	 -1.0
    r1625 	 0.0
    r1316 	 0.0
    PAIL4P_HStn 	 1
    AM4N9CSitr 	 1
    ELAIDCPT2 	 0.0
    THRFVShc 	 1
    r0669 	 1.0
    HSD17B3r 	 -1.0
    RE0580L 	 1
    r0753 	 3.0
    NTD7 	 3.0
    RBK_D 	 -1.0
    r2157 	 -1.0
    FUCGALGBSIDEtg 	 1
    ILETAm 	 0.0
    2HIBUPGLUC_Sthv 	 1
    ARTPLM3 	 1
    r2172 	 -1.0
    EX_crtsl_LPAREN_e_RPAREN_ 	 1
    ETF 	 3.0
    FAOXC163Gm 	 3.0
    FAOXC16C16OHm 	 3.0
    r0830 	 -1.0
    RE3230R 	 1
    CBLtle 	 1
    FAOXMC10OHMC10r 	 0.0
    OCTDECCPT2 	 0.0
    COKECBESr 	 0.0
    EX_uri_LPAREN_e_RPAREN_ 	 3
    RE2154C 	 1
    LSTNM7itr 	 1
    LS3 	 1
    ACNACNGALGBSIDEte 	 1
    r1536 	 1.0
    RE0908R 	 0.0
    HMGCOAtm 	 1
    RE2377C 	 1
    r1970 	 0.0
    GLCNACASE1ly 	 -1.0
    NDPK7m 	 0.0
    NDPK9m 	 0.0
    MI3PP 	 3.0
    NDPK6 	 0.0
    NDP6 	 3.0
    CRVSM23hr 	 3.0
    CBPPer 	 1.0
    DUTPDP 	 1
    LEUPHELAT2tc 	 1
    AM1C4N9CSteb 	 1
    OCTDEC2ACBP 	 3.0
    r1629 	 3.0
    RSVGLUChc 	 1
    EX_npthl_LPAREN_e_RPAREN_ 	 1
    ALDD20x 	 3.0
    HOXG 	 1.0
    AGLPR 	 1
    RE2677E 	 1
    ACGALFUCGALACGALFUC12GAL14ACGLCGALGLUSIDEtg 	 1
    4MTOLBUTAMIDEte 	 1
    r1876 	 0.0
    SELCYSLY2 	 -1.0
    RE1804M 	 0.0
    EX_g1p_LPAREN_e_RPAREN_ 	 1
    r1743 	 0.0
    RE3141X 	 3.0
    EX_thyox_L_LPAREN_e_RPAREN_ 	 1
    r1844 	 0.0
    S6TASE4ly 	 -1.0
    CMPACNAtn 	 1
    RE1815M 	 3.0
    PECTCHLe 	 1
    r1431 	 1.0
    r2045 	 -1.0
    GTPCI 	 3.0
    EX_sarcs_LPAREN_e_RPAREN_ 	 1
    FAOXC204184m2 	 3.0
    CSDPASEly 	 1
    r2012 	 -1.0
    P4508B11r 	 -1.0
    NMNS 	 3.0
    EX_lst4exp_LPAREN_e_RPAREN_ 	 1
    GLYLEUPEPT1tc 	 -1.0
    r0584 	 1
    ACACT1r 	 0.0
    SPHMYLNtl 	 1
    RE1441R 	 1
    r0921 	 1
    FUCACNGAL14ACGLCGALGLUSIDEtg 	 1
    CARIBUP_Sthv 	 1
    GDPFUCtg 	 -1.0
    RE2862C 	 1
    ECOAH9m 	 3.0
    FAOXC15ATPx 	 0.0
    CYSAMPtev 	 1
    GMPtn 	 1
    TDPm 	 1
    r0283 	 -1.0
    r0510 	 3.0
    GLCNACT1g 	 1.0
    CYSSERNaEx 	 0.0
    CRGLZtev 	 1
    HEXCCPT2 	 0.0
    EX_3hpvstet_LPAREN_e_RPAREN_ 	 1
    C50CPT1 	 1.0
    RE3146R 	 3.0
    NDP7g 	 3.0
    r1291 	 1
    r1894 	 0.0
    RE2920M 	 3.0
    RE3160C 	 1
    EX_glgchlo_LPAREN_e_RPAREN_ 	 1
    NADPtxu 	 1
    RE3301G 	 2.0
    G14T19g 	 3.0
    A_MANASE 	 -1.0
    r2491 	 -1.0
    MECOALm 	 1
    r1700 	 0.0
    r1741 	 0.0
    RE3179M 	 3.0
    S6TASE16ly 	 3.0
    NDP7er 	 3.0
    FATP1t 	 0.0
    DOCOSDIACTD 	 1
    PROSTGI2tr 	 1
    EX_nrpphr_LPAREN_e_RPAREN_ 	 1
    EX_bildglcur_LPAREN_e_RPAREN_ 	 1
    OXAtp 	 1
    RE2909C 	 1
    RE1532X 	 3.0
    TETPENT6t 	 1
    RE3397M 	 1
    FUC14GALACGLCGALGLUSIDEte 	 1
    RE1835M 	 3.0
    r1997 	 -1.0
    r1766 	 0.0
    GLC3MEACPtev 	 1
    CHSTEROLSULT 	 2.0
    UGALNACter 	 2.0
    FACOAL150 	 3.0
    C3DCe 	 1
    r0838 	 1
    XOLTRI27tc 	 1
    RE0512C 	 1
    r0221 	 0.0
    RE3120R 	 3.0
    r0365 	 1
    LYStip 	 1
    EX_4hpro_LPAREN_e_RPAREN_ 	 1
    r2112 	 2.0
    r2103 	 1
    RE3562X 	 3.0
    GLCNACT_U 	 1
    DGTPtn 	 1
    RE2865C 	 1
    EX_31dmt_LPAREN_e_RPAREN_ 	 1
    PVSATPtu 	 -1.0
    EX_peplys_LPAREN_e_RPAREN_ 	 1
    r2399 	 -1.0
    RE3342M 	 3.0
    r1751 	 0.0
    EX_octa_LPAREN_e_RPAREN_ 	 3
    OCDCAFATPtc 	 0.0
    r1500 	 -1.0
    RE3153C 	 1
    r2020 	 -1.0
    EX_dtmp_LPAREN_e_RPAREN_ 	 1
    PTVSTitr 	 1
    GLCNACASE4ly 	 -1.0
    FUT91g 	 -1.0
    NaKt 	 3.0
    RE3079X 	 0.0
    r1901 	 0.0
    IPDPtx 	 1
    ITCOAL1m 	 3.0
    MTAP 	 1.0
    VITD2Hm 	 0.0
    NTD1 	 3.0
    r1863 	 0.0
    RE0383C 	 3.0
    DHCR71r 	 0.0
    r1798 	 0.0
    ORETNtn2 	 1
    S6TASE10ly 	 -1.0
    PTRCOX1 	 0.0
    1331TAALThr 	 1
    S6TASE1ly 	 3.0
    FAOXC163_4Z_7Z_10Zm 	 3.0
    r1260 	 3.0
    r0055 	 1
    r0268 	 1
    XSERtg 	 1
    CEPTE 	 -1.0
    EX_1531tacr_LPAREN_e_RPAREN_ 	 1
    NTPP10 	 3.0
    NDPK1m 	 0.0
    CYSSNAT5tc 	 -1.0
    GASNASE3ly 	 3.0
    RE3108C 	 1
    PI34P4Pn 	 1
    PFK 	 3.0
    CRVSATPthc 	 1
    DUMPtn 	 1
    FRDPtr 	 1
    UGT1A6r 	 0
    DNDPt9m 	 -1.0
    ATVLAChc 	 1
    KYNAKGAT 	 0.0
    DM_4hrpo 	 3
    biomass_other 	 3
    ST8SIA56g 	 -1.0
    r0936 	 1
    ACONT 	 1.0
    TREH 	 0.0
    1513TACRtep 	 1
    BACCLm 	 -1.0
    r0767 	 0.0
    C10DCCACT 	 -1.0
    NOS2 	 -1.0
    PA_HStn 	 1
    KYN3OX 	 0.0
    r0715 	 3.0
    RAHY 	 1
    FAOXTC142TC122m 	 3.0
    GLCNACT4g 	 2.0
    RE3493C 	 1
    C182OHc 	 1.0
    FAOXC183806x 	 -1.0
    r0763 	 0.0
    RE0925C 	 -1.0
    6EPVShc 	 1
    5THFtl 	 1
    DCATDr 	 1
    r2488 	 -1.0
    DOLPGT2_Uer 	 1
    TETTET6t 	 1
    RE3436C 	 1
    NACHEXA10ly 	 3.0
    RE3003M 	 3.0
    RSVhc 	 1
    RE1906C 	 1
    r2224 	 -1.0
    CSND 	 1
    RE2151C 	 -1.0
    TDPDRR 	 1
    PI4PLCn 	 1.0
    DMHPTCRNCPT2 	 1
    HMGCOASim 	 -1.0
    AG13T4g 	 3.0
    r0781 	 0
    GASNASEly 	 3.0
    EX_maltttr_LPAREN_e_RPAREN_ 	 1
    MI1345PKn 	 0.0
    RE2513N 	 -1.0
    r2404 	 -1.0
    r1992 	 -1.0
    r2043 	 -1.0
    r2322 	 1.0
    1331TACRhr 	 1
    RE1635X 	 3.0
    RE2853C 	 1
    RN0032R 	 0.0
    r1375 	 1
    ME1m 	 3.0
    GPIMTer_U 	 1
    r0149 	 2.0
    SO4CLtex2 	 -1.0
    HSAT1ly 	 1
    PI34P5K 	 0.0
    C6COAt 	 1.0
    EX_drib_LPAREN_e_RPAREN_ 	 3
    r1986 	 0.0
    SUBERICACT 	 1
    24_25VITD3Hm 	 1.0
    CTPtn 	 1
    CRVSATPtu 	 -1.0
    EX_limnen_LPAREN_e_RPAREN_ 	 1
    OXYPthc 	 -1.0
    PCHOL_HStg 	 1
    FAOXC150130m 	 3.0
    RE1587C 	 1
    RE3521M 	 3.0
    RE2319M 	 3.0
    EX_gdp_LPAREN_e_RPAREN_ 	 1
    ACMPGLUTthc 	 0.0
    EX_am1cglc_LPAREN_e_RPAREN_ 	 1
    ACGALFUCGALACGALFUC12GAL14ACGLCGALGLUSIDEte 	 1
    FAOXTC122m 	 3.0
    L_LACt2r 	 2.0
    NTD8l 	 2.0
    RE3502X 	 0
    r2005 	 -1.0
    r1813 	 0.0
    RE3596M 	 0
    EX_gam_LPAREN_e_RPAREN_ 	 3
    r0741 	 1
    r0720 	 1.0
    SIAASEly 	 -1.0
    RE3470X 	 0
    MM6ag 	 0.0
    RE1519X 	 1
    RAI1 	 1
    r1328 	 1
    COQ5m 	 1
    RE3596X 	 0
    RE2251C 	 1
    RE1711M 	 -1.0
    RE3161R 	 3.0
    PTDCACRNt 	 -1.0
    ACGALtlg 	 1
    S6T18g 	 0.0
    FUCACNGALACGLCGALGLUSIDEte 	 1
    r1758 	 0.0
    ACGAGBSIDEtg 	 1
    r1380 	 3.0
    RE3520M 	 1
    GLZABCteb 	 1
    CYSTLEUrBATtc 	 -1.0
    LSTO2r 	 3.0
    BILDGLCURtr 	 1
    r2116 	 2.0
    RE3559M 	 3.0
    ARGB0AT3tc 	 -1.0
    ACGSm 	 -1.0
    NAt3_1 	 0
    RE3587N 	 3.0
    r1660 	 3.0
    PPCOACm 	 0.0
    NDERSVteb 	 -1.0
    r1763 	 0.0
    NTD6 	 3.0
    P5CRm 	 -1.0
    M4MPDOL_Uter 	 1
    RE2917X 	 3.0
    r0450 	 0.0
    FATP9t 	 0.0
    RE2759X 	 -1.0
    OCBTm 	 -1.0
    4HATVLACtep 	 1
    TIGCRNe 	 1
    RE3550X 	 1
    r0527 	 1
    EX_h_LPAREN_e_RPAREN_ 	 3
    r2385 	 0.0
    r0722 	 3.0
    GF6PTA 	 3.0
    BILGLCURt 	 0.0
    5ADTSTSTERONEtr 	 1
    DSPVSitr 	 1
    CDS 	 3.0
    GUMGCHLe 	 1
    G1M7MASNBterg 	 1
    r2315 	 1.0
    GLCAASE5ly 	 3.0
    TMDM3hr 	 1
    r0474 	 1.0
    MHGLZhr 	 1.0
    DGCHOLtx 	 1
    GMPtg 	 1
    UGT1A5r2 	 0
    AM1ALCSitr 	 1
    6HTSTSTERONEtr 	 1
    FAOXC122_3E_6Em 	 3.0
    EX_imp_LPAREN_e_RPAREN_ 	 1
    RE3506R 	 0.0
    PGLYCP 	 1.0
    C9BRxtc 	 -1.0
    RE0935C 	 1
    ACN13ACNGALGBSIDEtg 	 1
    ST6GALNAC25 	 3.0
    EICOSTETCPT2 	 0.0
    4HDEBRISOQUINEte 	 1
    r2130 	 1
    H8MTer_L 	 0.0
    UMPK 	 3.0
    6HLVSTACIDitr 	 1
    RE1941R 	 1.0
    r0383 	 3.0
    ADNCYC 	 2.0
    GACPAILter 	 1
    6HSMVhep 	 1.0
    RE0827E 	 0.0
    r2328 	 1.0
    LNLNCAt 	 -1.0
    r0571 	 1
    r2029 	 -1.0
    RE3189M 	 3.0
    r0480 	 1
    RE0574C 	 1
    RE1522M 	 3.0
    GSNKm 	 1
    r2175 	 -1.0
    PAFHe 	 1.0
    ALA_DASH_DTDe 	 1
    FAOXC164_4Z_7Z_10Z_13Zm 	 3.0
    MALtm 	 -1.0
    HXt 	 1
    4ABUTtm 	 1
    FACOAL60i 	 1
    RE3308C 	 1
    r1683 	 0.0
    34HPLFM 	 1
    G14T9g 	 3.0
    RE3554R 	 3.0
    r1947 	 0.0
    KHK2 	 -1.0
    G14T12g 	 3.0
    KHK3 	 -1.0
    EX_gal_LPAREN_e_RPAREN_ 	 3
    METyLATthc 	 0.0
    GLCAT7g 	 2.0
    r1183 	 3.0
    r0693 	 0.0
    4HPRO_LTte 	 1
    M4ATAer 	 -1.0
    BTNtn 	 1
    C4CRNCPT2 	 0.0
    TMDOATPtsc 	 1
    r2021 	 -1.0
    RE3016R 	 1
    r1711 	 0.0
    NACHEXA1ly 	 3.0
    PCm 	 0.0
    GLRASE 	 1
    r1019 	 1
    CYSTALArBATtc 	 -1.0
    r0648 	 1.0
    GALFUCGALACGLCGAL14ACGLCGALGLUSIDEtg 	 1
    RE2636R 	 -1.0
    EX_decdicrn_ 	 1
    EX_c101crn_ 	 1
    2OXOADPTm 	 -1.0
    NTD11 	 3.0
    PMI12346PHn 	 0.0
    FPGS8m 	 -1.0
    EX_ddecrn_ 	 1
    r1071 	 1
    FATP3t 	 -1.0
    CHSTEROLt3 	 1
    GLUPROASCT1 	 0.0
    FAOXC142_5Z_8Zm 	 3.0
    r0963 	 3.0
    r1936 	 0.0
    RE2514L 	 3.0
    3SALATAim 	 3.0
    GLCtg 	 3.0
    GGH_5THFe 	 3.0
    GTHRDtr 	 1
    BMTer_U 	 1
    EX_3mop_LPAREN_e_RPAREN_ 	 1
    AM9CSAitr 	 1
    RE3346M 	 3.0
    NDP3ex 	 3.0
    G14Tg 	 3.0
    MGSA2 	 1
    THFt2 	 -1.0
    ATVETHGLUChc 	 0
    AMCOXO 	 1
    EX_lpchol_hs_LPAREN_e_RPAREN_ 	 3
    r1561 	 3.0
    PROt4_2_r 	 -1.0
    CYSTGLUex 	 0.0
    TAGAT_Dt 	 1
    RE1520X 	 1.0
    MERACMPtep 	 0.0
    TKT2 	 3.0
    RE3237C 	 1
    RE3557R 	 3.0
    r0926 	 1
    SFGTH 	 3.0
    RE3367X 	 0.0
    RE2649X 	 0.0
    RE2523X 	 3.0
    r1836 	 0.0
    ACALDt 	 1
    UGALNACtg 	 2.0
    r2149 	 -1.0
    RE1635R 	 3.0
    LEUKTRD4t 	 1
    EX_5adtststerone_LPAREN_e_RPAREN_ 	 1
    S6T3g 	 0.0
    r2152 	 -1.0
    r2024 	 -1.0
    H2Oter 	 1
    SMVGLUChep 	 0
    RE2813C 	 1
    RE3526C 	 1.0
    CLPNDCPT2 	 0.0
    ACCOAtr 	 3.0
    RE3181C 	 1
    r2490 	 -1.0
    42A12BOOX 	 3.0
    PI4PLC 	 1.0
    EX_am19cs_LPAREN_e_RPAREN_ 	 1
    r2400 	 -1.0
    B_MANNASEly 	 0.0
    LYStn 	 1
    RE2127C 	 1
    DOCOSACTDr 	 1
    TYRCBOX 	 0.0
    NCAMDe 	 1
    r2091 	 2.0
    r1029 	 1
    SERGLNexR 	 -1.0
    FAOXC9BRC7BRm 	 1
    RE1905R 	 -1.0
    sink_octdececoa_LPAREN_c_RPAREN_ 	 1
    RE3339C 	 0
    INSt5le 	 0.0
    CORE3GTg 	 -1.0
    RE2427M 	 -1.0
    r0424 	 3.0
    EX_s3meacmp_LPAREN_e_RPAREN_ 	 1
    DHCR242r 	 3.0
    HEXDICOAACBPx 	 3.0
    r2391 	 0.0
    RE1236C 	 1
    G14T11g 	 3.0
    r1911 	 0.0
    ST8SIA52g 	 -1.0
    S6T1g 	 0.0
    DOLICHOL_Uter 	 1
    r2188 	 -1.0
    ADCim 	 1
    1513DTALThr 	 1
    G14T18g 	 3.0
    RE3176R 	 1
    PTPAT 	 -1.0
    RE3038X 	 0.0
    6MELVSTthep 	 1
    ITPtn 	 1
    r1522 	 3.0
    DCIm 	 3.0
    GSNt5le 	 0.0
    RE2622R 	 0
    S3TASE3ly 	 1
    KDNH 	 1
    Htx 	 1
    VALB0AT3tc 	 -1.0
    RN0023C 	 0.0
    GALU 	 3.0
    ADK3 	 -1.0
    SGPL12r 	 3.0
    r0907 	 1
    r0859 	 1
    DNDPt13m 	 -1.0
    PIt2m 	 3.0
    CYSGLNNaEx 	 0.0
    PPCOAOm 	 3.0
    CYANtm 	 1
    34DHXMANDACOX 	 3.0
    r1996 	 -1.0
    PRDX 	 3.0
    GGH_10FTHF7GLUl 	 3.0
    HMBS 	 -1.0
    EX_ddca_LPAREN_e_RPAREN_ 	 3
    EX_isolvstacid_LPAREN_e_RPAREN_ 	 1
    NACHEXA22ly 	 3.0
    r0691 	 0.0
    FAOXC2051843x 	 -1.0
    EX_cbasp_LPAREN_e_RPAREN_ 	 1
    GAMYe 	 -1.0
    r1732 	 0.0
    GLCAASE6ly 	 3.0
    EX_amp_LPAREN_e_RPAREN_ 	 3
    RE1631C 	 1
    S3T2g 	 0.0
    RE1818X 	 3.0
    NTD5m 	 -1.0
    r2129 	 2.0
    r0047 	 3.0
    FAOXC163C143m 	 3.0
    EBP1r 	 3.0
    RE3533C 	 1
    r1860 	 0.0
    RE3560X 	 0.0
    NRPPHRtu 	 1.0
    TPI 	 3.0
    RE3435R 	 2.0
    r1984 	 0.0
    PAFt 	 1
    RE2514N 	 -1.0
    ADEt 	 -1.0
    EX_oxy1rb_LPAREN_e_RPAREN_ 	 1
    XOL7AONEtr 	 1
    RE3452C 	 1
    r2384 	 0.0
    FAOXC184_6Z_9Z_12Z_15Zx 	 -1.0
    FALDH 	 3.0
    r1991 	 0.0
    TYRASE 	 0.0
    NDPK1n 	 0
    RE1816M 	 3.0
    EX_gchola_LPAREN_e_RPAREN_ 	 3
    r2194 	 -1.0
    r0021 	 2.0
    r0434 	 1
    ADMDC 	 2.0
    NICRNS 	 1
    3SPYRSP 	 1
    r0186 	 0.0
    EX_pvsgluc_LPAREN_e_RPAREN_ 	 1
    r2191 	 -1.0
    DOLPGT3_Uer 	 1
    EX_xoltri27_LPAREN_e_RPAREN_ 	 1
    CRVSM22hc 	 1
    DGK1 	 3.0
    ACOAHi 	 -1.0
    GLCNACASE3ly 	 -1.0
    EX_leuktrD4_LPAREN_e_RPAREN_ 	 3
    r2413 	 -1.0
    RDH1a 	 3.0
    FAOXC226m 	 3.0
    RN0031C 	 0.0
    EX_malttr_LPAREN_e_RPAREN_ 	 1
    r2447 	 0.0
    r1835 	 0.0
    NACHEX6ly 	 3.0
    3SALATAi 	 3.0
    r1610 	 3.0
    ELAIDCRNt 	 -1.0
    ATVLACGLCURhc 	 0
    URIDK2m 	 0.0
    RE2525C 	 1
    RE0958E 	 0
    r1723 	 0.0
    RE2908X 	 3.0
    r2526 	 -1.0
    NAtx 	 1
    EX_adrnl_LPAREN_e_RPAREN_ 	 1
    r2407 	 -1.0
    DIGALSGALSIDEte 	 1
    RE3444X 	 0.0
    r1912 	 0.0
    r1557 	 3.0
    r1856 	 0.0
    DNDPt46m 	 -1.0
    NAGAly 	 0.0
    G14T14g 	 3.0
    ANDRSTRNGLCtr 	 1
    3HKYNAKGAT 	 0.0
    FAOXC101C102x 	 -1.0
    r2506 	 1
    r2517 	 3.0
    APOC_LYS_BTNP 	 1
    ASPDxt 	 1
    r1818 	 0.0
    CYOOm3 	 -1.0
    MHGLZABCt 	 1
    4HATVACIDteb 	 -1.0
    RE3229R 	 1
    NACHEX17ly 	 3.0
    DNDPt43m 	 -1.0
    ATVLAChr 	 1
    EBASTINEOHtr 	 1
    RE3000M 	 3.0
    TETDEC2CRNe 	 1
    r1146 	 1
    RE2666G 	 0.0
    TAURIBUP_Sthv 	 1
    LNLNCGCRNt 	 -1.0
    STACMPtev 	 1
    RE3534C 	 1
    r0430 	 0.0
    H7_TAer 	 -1.0
    HISt4 	 3.0
    RE3454C 	 1
    DEDOLR_U 	 1
    r0620 	 1.0
    EX_strdnc_LPAREN_e_RPAREN_ 	 1
    FVSTETitr 	 1
    EX_HC02216_LPAREN_e_RPAREN_ 	 1
    FAOXC130110m 	 3.0
    EX_arachcoa_LPAREN_e_RPAREN_ 	 1
    GALGT4 	 0.0
    NACHEX23ly 	 3.0
    RE3286C 	 1
    B3GNT312g 	 0.0
    r1045 	 0.0
    FAOXCPRIST3x 	 -1.0
    NACHEXA4ly 	 3.0
    DPHMBDCm 	 1
    THRPHELAT2tc 	 1
    GLCNACASE5ly 	 -1.0
    EX_13_cis_retnglc_LPAREN_e_RPAREN_ 	 1
    OMHDEACIDTD 	 1
    r1006 	 1
    EX_tetdece1crn_ 	 1
    C60CPT1 	 1.0
    RE0577X 	 0.0
    r1013 	 1
    OMPDC 	 3.0
    OCD11CRNCACT 	 1
    RE1573M 	 3.0
    EX_ser_L_LPAREN_e_RPAREN_ 	 3
    EX_pvs_LPAREN_e_RPAREN_ 	 1
    DM_hretn_n_ 	 1
    RETH2e 	 1
    XAOx 	 2.0
    DM_core5_g_ 	 1
    MAN1_7Ber 	 0.0
    EX_crn_LPAREN_e_RPAREN_ 	 3
    r2156 	 -1.0
    RE3342X 	 3.0
    RE2146C 	 -1.0
    G14T17g 	 3.0
    EX_25hvitd2_LPAREN_e_RPAREN_ 	 1
    EX_mag_hs_LPAREN_e_RPAREN_ 	 3
    r0747 	 -1.0
    RE3562R 	 3.0
    S3TASE1ly 	 1
    RE3287R 	 0.0
    FACOAL161 	 3.0
    FVStep 	 1
    DMATTx 	 3.0
    WHTTDCAte 	 1
    PRISTCOAtx 	 1
    sink_5hpet_LPAREN_c_RPAREN_ 	 1
    O2Stm 	 1
    ACACT5p 	 3.0
    TCYNTtm 	 1
    BTND1n 	 -1.0
    D_LACtm 	 2.0
    DNDPt37m 	 -1.0
    r2232 	 -1.0
    r2378 	 0.0
    DCK2n 	 2.0
    GQ1BALPHAte 	 1
    H2CO3D 	 3.0
    DITPtn 	 1
    FUT14g 	 0.0
    NAt 	 3.0
    DNDPt60m 	 -1.0
    EX_tststeroneglc_LPAREN_e_RPAREN_ 	 1
    PYAM5Ptm 	 1
    RE3496C 	 1
    RE3159X 	 1.0
    6OHFVSGLUitr 	 1
    EX_q10h2_LPAREN_e_RPAREN_ 	 3
    PRPPS 	 3.0
    C16DCe 	 1
    r0546 	 3.0
    AKR1C42 	 -1.0
    OCCOAtm 	 1
    EX_crvsm24_LPAREN_e_RPAREN_ 	 1
    r1784 	 0.0
    RE3089X 	 0.0
    r2160 	 -1.0
    r1492 	 -1.0
    CSPG_Dtly 	 1
    GLYt7_311_r 	 -1.0
    FACOAL1832 	 3.0
    5ADTSTSTERONEte 	 1
    NTD12 	 0.0
    RE1916X 	 3.0
    RE1835X 	 0.0
    CDIPTr 	 3.0
    RE3093X 	 3.0
    PUNP1 	 3.0
    r0193 	 0.0
    r1621 	 0.0
    UDPGLCAter 	 2.0
    SARCStm 	 1
    r1921 	 0.0
    4HATVLACthc 	 0.0
    r2217 	 -1.0
    SMVGLUCLAChep 	 1
    ASNATB0tc 	 3.0
    EX_prostge1_LPAREN_e_RPAREN_ 	 1
    r2331 	 1.0
    r1973 	 0.0
    MG2er 	 3.0
    S6T24g 	 0.0
    r0497 	 2.0
    FAOXC121_5Em 	 3.0
    MCOATAm 	 -1.0
    RE0124C 	 1
    ABTti 	 1
    r0987 	 1
    r0696 	 0.0
    r1995 	 -1.0
    METS 	 0.0
    3HIBUPGLUC_Sthv 	 1
    4NPHSFte 	 1
    ACNAMlt 	 -1.0
    FBA2 	 3.0
    OCCOAtx 	 1
    KSII_CORE4t 	 1
    r1324 	 1
    EX_xoltri25_LPAREN_e_RPAREN_ 	 1
    EX_aldstrn_LPAREN_e_RPAREN_ 	 1
    7HPVSitr 	 1
    r0610 	 -1.0
    PI3P4K 	 0.0
    CYTK12 	 3.0
    r0081 	 0.0
    THMTPt 	 -1.0
    3HPVSTETCOAitm 	 1
    EX_xylt_LPAREN_e_RPAREN_ 	 3
    RE1834C 	 0
    FAOXC182C182OHm 	 3.0
    r1749 	 0.0
    r1656 	 3.0
    r1418 	 3.0
    DM_dctp_m_ 	 1
    DELACCRVSM23hc 	 1
    LVSTACOXD6Hhep 	 1.0
    PIK3 	 3.0
    r1615 	 0.0
    EICOSTETt 	 1
    4HATVLACteb 	 -1.0
    r0756 	 3.0
    r2255 	 -1.0
    RE3597X 	 0
    r0786 	 3.0
    BGLUTCHLe 	 1
    r2393 	 0.0
    ARGSL 	 0.0
    FAOXTC102C101m 	 3.0
    UPP3S 	 -1.0
    r1421 	 1
    r0034 	 1
    r2503 	 -1.0
    3HSMVhep 	 1.0
    NACHEX27ly 	 3.0
    FUT16g 	 0.0
    r1881 	 0.0
    FAOXCPRIST1x 	 -1.0
    r2041 	 -1.0
    DOPAVESSEC 	 0.0
    RE1582L 	 1
    r1169 	 1.0
    DNDPt34m 	 1
    EX_oh1_LPAREN_e_RPAREN_ 	 1
    RE3172R 	 1
    EX_HC02204_LPAREN_e_RPAREN_ 	 1
    RE0924C 	 -1.0
    DHEASABCCte 	 -1.0
    EX_itp_LPAREN_e_RPAREN_ 	 1
    TSTSTERONEGLCte 	 3.0
    LNS14DM 	 0.0
    SLCBK1 	 1.0
    B3GALT44g 	 -1.0
    34HPPOR 	 -1.0
    C120CPT2 	 0.0
    RE3250X 	 3.0
    EX_cysacmp_LPAREN_e_RPAREN_ 	 1
    FUCGAL14ACGLCGALGLUSIDEtg 	 1
    r2330 	 1.0
    RE3233L 	 -1.0
    XOLDIOLONEtm 	 1
    PGLer 	 0.0
    RE1860E 	 3.0
    FADtru 	 1
    H6ET3er 	 0.0
    RE2636C 	 1
    3NTD7l 	 2.0
    FACOAL246_1 	 3.0
    r2198 	 -1.0
    GLACter 	 1
    CO2ter 	 1
    PPPItn 	 1
    FE3R2e 	 3.0
    EX_3octdec2crn_ 	 1
    RE2954C 	 3.0
    RE3201C 	 1
    r1605 	 3.0
    DLNLCGCRNt 	 -1.0
    P45021A2r 	 0
    RNDR4 	 3.0
    C4STMO2Pr 	 2.0
    RETNGLCtr 	 1
    RN0027C 	 1
    THRGLNNaEx 	 0.0
    CHSTEROLt 	 3.0
    r1593 	 3.0
    DOLPMT1_Uer 	 1
    TETTET6CPT1 	 1.0
    r2099 	 1
    r0762 	 0.0
    7BHGLZGLCtev 	 1
    PROSTGH2t 	 1
    RE3443X 	 3.0
    MMEm 	 1.0
    AIRCr_PRASCS 	 3.0
    r1021 	 1
    RE3010C 	 3.0
    FAOXC122x 	 1.0
    1a_24_25VITD3Hm 	 1
    HIStiDF 	 -1.0
    6EPVStep 	 1
    HDECAACBP 	 3.0
    FAOXC225m 	 3.0
    r0886 	 1
    HEXDCRNe 	 1
    CMPSASn 	 0.0
    S6T19g 	 0.0
    r0911 	 1
    RE3244R 	 1
    CYTDK2m 	 1
    HEXDIACTD 	 1
    RE1533M 	 3.0
    EX_3tdcrn_ 	 1
    r2136 	 -1.0
    r0627 	 3.0
    C101CRNe 	 1
    r0013 	 2.0
    FA160ACPH 	 0.0
    r2018 	 -1.0
    FOLOAT2tc 	 -1.0
    TRDR3 	 3.0
    RE3469C 	 1
    r1890 	 0.0
    C142ACBP 	 3.0
    r0578 	 3.0
    SPH1Ptr 	 1
    RE2642C 	 1
    1HMDGLUCitr 	 1
    GLYSARCNc 	 3.0
    LEUGLYPEPT1tc 	 -1.0
    ACACt2m 	 2.0
    r1952 	 0.0
    r0097 	 0.0
    DDECCRNe 	 1
    DM_mem2emgacpail_prot_hs_r_ 	 1
    RTOTALCRNCPT1 	 1.0
    EX_galside_hs_LPAREN_e_RPAREN_ 	 1
    NODe 	 1
    FAOXC205C184x 	 1.0
    C180CRNt 	 -1.0
    RE2296C 	 1
    EX_CE1950_LPAREN_e_RPAREN_ 	 1
    GUMTCHOLe 	 1
    GLUt6 	 1.0
    r2210 	 -1.0
    GLCAASE7ly 	 3.0
    GALASE1ly 	 0
    r1915 	 0.0
    RE1806C 	 1
    DCYTDn 	 -1.0
    RE1922C 	 1
    r0384 	 1.0
    RE1521X 	 1.0
    RE2051C 	 3.0
    PI345P3Pn 	 3.0
    RE2122C 	 1
    PCHOLPg_hs 	 2.0
    DM_Ser_Gly_Ala_X_Gly_ly_ 	 1
    r1563 	 3.0
    r2008 	 -1.0
    r0679 	 3.0
    RE2445C 	 0.0
    UAGDP 	 3.0
    r1584 	 3.0
    PVSGLUCitr 	 1
    IPDDI 	 3.0
    r1501 	 -1.0
    FAOXC225_4Z_7Z_10Z_13Z_16Zm 	 3.0
    r2155 	 -1.0
    EX_dcmp_LPAREN_e_RPAREN_ 	 3
    RE3232C 	 1
    EX_2425dhvitd2_LPAREN_e_RPAREN_ 	 1
    GAL3ST11 	 -1.0
    EX_lac_D_LPAREN_e_RPAREN_ 	 1
    r1488 	 1
    MI14P4P 	 1
    NABTNOm 	 3.0
    PHYQt 	 1
    DNDPt55m 	 -1.0
    RE3400M 	 3.0
    CHSTEROLt1 	 -1.0
    TETPENT6COAtx 	 1
    HIBDm 	 3.0
    C6CRNe 	 1
    ARGLYSex 	 0.0
    MM7Cbg 	 3.0
    GAO2 	 3.0
    UMPK6n 	 3.0
    RE3170R 	 3.0
    RBTt 	 1
    NRVNCCOAtx 	 1
    RE2541C 	 0.0
    RN0001R 	 0.0
    DOCOSACT 	 -1.0
    EX_fucfucfucgalacglc13galacglcgal14acglcgalgluside_hs_LPAREN_e_RPAREN_ 	 1
    AGTix 	 -1.0
    ACNAMtn 	 1
    ABTD 	 1
    CBR2 	 2.0
    r1587 	 3.0
    EX_n2m2nmasn_LPAREN_e_RPAREN_ 	 3
    EX_Lcystin_LPAREN_e_RPAREN_ 	 1
    CBL2tm 	 1.0
    ESTRADIOLt 	 1
    DNDPt31m 	 -1.0
    r1822 	 0.0
    r2520 	 -1.0
    RE1923C 	 1
    RE2705C 	 1
    FAOXC10DCC8DCx 	 -1.0
    r1746 	 0.0
    ALLOPtepvb 	 -1.0
    RE2150C 	 -1.0
    RE1534M 	 3.0
    NTD4 	 3.0
    S2TASE2ly 	 3.0
    ECOAH12m 	 3.0
    r2436 	 -1.0
    RE3338M 	 3.0
    r1025 	 1
    RE0944C 	 1
    ALAALAPEPT1tc 	 -1.0
    r0726 	 3.0
    EX_ocdcea_LPAREN_e_RPAREN_ 	 3
    RN0031R 	 2.0
    IDOAASE4ly 	 -1.0
    TSACMSULhc 	 1
    S6TASE11ly 	 3.0
    ARTFR31 	 1
    r0560 	 -1.0
    PTE5x 	 0
    CHOLK 	 0.0
    EX_sbt_DASH_d_LPAREN_e_RPAREN_ 	 3
    C141ACBP 	 3.0
    r1427 	 -1.0
    BIOCYTtn 	 1
    LSTO1r 	 3.0
    HMGCOAtx 	 1
    RE0344M 	 3.0
    FAOXC200180x 	 -1.0
    r1332 	 1
    SUBERCROT 	 1.0
    r2269 	 -1.0
    r0446 	 1
    CAROtr 	 1
    EX_CE5797_LPAREN_e_RPAREN_ 	 1
    AM1CCStev 	 1
    r0022 	 2.0
    RE3164C 	 1
    r1364 	 1
    6BHGLZGLCtev 	 1
    GLU5Km 	 3.0
    RE3308R 	 1
    EX_ebastine_LPAREN_e_RPAREN_ 	 1
    RE3582X 	 1.0
    ADPT 	 3.0
    r2086 	 2.0
    CHLPCTD 	 1.0
    41R1H2MAE12BOOX 	 3.0
    CK 	 0
    P4507B12r 	 -1.0
    RE1938R 	 1.0
    RE1629C 	 1
    ENMAN6g 	 1
    r2238 	 -1.0
    UNK3 	 1
    r0975 	 1
    RE2912M 	 3.0
    FAOXC241C221x 	 -1.0
    SARDHm 	 -1.0
    RE3011M 	 -1.0
    MM8Ag 	 3.0
    RE3498R 	 1
    PTDCACRNCPT2 	 0.0
    DORNOp 	 -1.0
    AMY1e 	 0
    GLNS 	 3.0
    SUCOAS1m 	 3.0
    CHOLESTTDe 	 1
    r2050 	 -1.0
    EX_leuktrC4_LPAREN_e_RPAREN_ 	 1
    ENMAN4g 	 1
    IVCRNe 	 1
    FAOXC184C163x 	 1.0
    CRVSthc 	 -1.0
    RE3455C 	 1
    CARIBUPGLU_Sitr 	 1
    RE3174R 	 3.0
    RE3381E 	 0.0
    DUTPDPn 	 3.0
    THMTP 	 -1.0
    DMHPTCRNt 	 1
    ASNALANaEx 	 0.0
    FPGS3 	 -1.0
    r2502 	 -1.0
    EX_mthgxl_LPAREN_e_RPAREN_ 	 1
    RE3301C 	 -1.0
    TS3 	 1
    PPBNGS 	 2.0
    CHOLtr 	 1
    r0682 	 0.0
    RE3044C 	 1
    EX_4ohmdz_LPAREN_e_RPAREN_ 	 1
    GALNACT1g 	 1.0
    GLYCTO1p 	 -1.0
    RE3176C 	 1
    AKR1C1 	 3.0
    ACMPGLUTdt 	 1
    RE3503C 	 1
    r1761 	 0.0
    r0649 	 1.0
    RE2985M 	 -1.0
    O2tp 	 1
    S6TASE7ly 	 -1.0
    PYRt2m 	 2.0
    RE3125C 	 1
    URATEt 	 1
    7AHGLZitr 	 1
    1HIBUPGLUitr 	 1
    ADNK1 	 3.0
    CYTK8 	 3.0
    r1902 	 0.0
    THFtl 	 1
    RE2524C 	 1
    EX_HC02154_LPAREN_e_RPAREN_ 	 1
    FAOXC2452256x 	 -1.0
    BETBGTtc 	 -1.0
    GLYLEUHYDROc 	 3.0
    2HATVLAChc 	 1
    EX_na1_LPAREN_e_RPAREN_ 	 3
    PTPATe 	 1
    AM1C4N9CShc 	 1
    MCCCrm 	 2.0
    RE2031M 	 1
    sink_dd2coa_LPAREN_c_RPAREN_ 	 1
    r2133 	 2.0
    RE2032M 	 1
    r2163 	 -1.0
    RE3166C 	 1
    DPCOAK 	 -1.0
    EX_56eppvs_LPAREN_e_RPAREN_ 	 1
    RN0020R 	 2.0
    EX_ndersv_LPAREN_e_RPAREN_ 	 1
    GHMT2r 	 -1.0
    RE3154C 	 1
    EX_dgtp_LPAREN_e_RPAREN_ 	 1
    RE3157X 	 0.0
    ASNSERNaEx 	 0.0
    ARGN 	 0.0
    OPAHir 	 0.0
    ANDRSTRNtr 	 1
    EX_6melvst_LPAREN_e_RPAREN_ 	 1
    RE3119R 	 1
    r1503 	 -1.0
    THYOXt 	 1.0
    ORNTArm 	 3.0
    GUR1PP 	 1
    RAtn 	 1
    DOPABMO 	 0.0
    RE3570C 	 1
    r1607 	 3.0
    ADAe 	 1.0
    GLNASNNaEx 	 0.0
    MEVK1c 	 3.0
    AM9CSAhr 	 1.0
    MACACI 	 -1.0
    RE3250C 	 1
    CSNt 	 -1.0
    r1549 	 3.0
    RE1526M 	 3.0
    G3PD1 	 -1.0
    EX_idp_LPAREN_e_RPAREN_ 	 1
    EX_6ahglz_LPAREN_e_RPAREN_ 	 1
    RE2866C 	 1
    GSNt2r 	 0.0
    EX_HC02217_LPAREN_e_RPAREN_ 	 1
    DM_Ser_Thr_ly_ 	 1
    EX_aprgstrn_LPAREN_e_RPAREN_ 	 1
    P4503A7r 	 0
    MAL_Lte 	 1
    RE1845C 	 0.0
    MELATN23DOX 	 2.0
    CARhPTtc 	 1.0
    DOLDPP_Uer 	 1
    r1313 	 0.0
    HSD3B2r 	 -1.0
    RE2030M 	 1
    H2O2tp 	 1
    r0813 	 0.0
    DOPAENT4tc 	 -1.0
    EX_ivcrn_ 	 1
    RE3390M 	 3.0
    r2089 	 2.0
    RE2514C 	 3.0
    r1769 	 0.0
    DTTPtn 	 1
    HDD2COAtx 	 1
    ACOX22x 	 -1.0
    r1560 	 3.0
    DCYTt 	 -1.0
    RE3514R 	 0.0
    RE3564X 	 3.0
    C161CPT12 	 1.0
    EX_ptrc_LPAREN_e_RPAREN_ 	 3
    EX_HC02160_LPAREN_e_RPAREN_ 	 1
    r0626 	 1
    r1608 	 3.0
    RE3414C 	 1
    TMDM3itr 	 1
    RE3431C 	 1
    RE1527M 	 3.0
    ACALDtx 	 1
    2AMACSULT 	 0
    r1180 	 3.0
    DESFVSteb 	 1
    CYSGLTH 	 1
    LNLNCACPT1 	 1.0
    DLNLCGt 	 -1.0
    MI1346PKn 	 0.0
    r0596 	 3.0
    CRNCARtp 	 -1.0
    r1977 	 0.0
    A4GNT1g 	 -1.0
    GALASE7ly 	 0
    r1933 	 0.0
    ACGALK 	 1
    RE1062C 	 1
    ASNt4 	 3.0
    RE2398R 	 0.0
    FPGS6 	 -1.0
    TYRB0AT3tc 	 -1.0
    AHEXASE2ly 	 3.0
    r0598 	 1
    GMAND 	 0.0
    EX_oxy7rb_LPAREN_e_RPAREN_ 	 1
    r0908 	 1
    EX_acgalfucgalacgalfuc12gal14acglcgalgluside_hs_LPAREN_e_RPAREN_ 	 1
    r1734 	 0.0
    G2M8MASNterg 	 1
    CO2tp 	 1
    ME2 	 2.0
    r0670 	 -1.0
    HDCACBP 	 3.0
    RE3495C 	 1
    r1842 	 0.0
    DM_4abut_LPAREN_n_RPAREN_ 	 1
    RE1527X 	 3.0
    MI134PP 	 3.0
    EX_5ohfvsglu_LPAREN_e_RPAREN_ 	 1
    AQCOBALt 	 1
    RADH4 	 1
    EX_arachd_LPAREN_e_RPAREN_ 	 3
    RE3413C 	 1
    r2166 	 -1.0
    EX_bglc_LPAREN_e_RPAREN_ 	 1
    r1948 	 0.0
    PTVSTATPtu 	 -1.0
    SLDx 	 3.0
    CLFORtex 	 0.0
    S4T3g 	 0.0
    NACHEX14ly 	 3.0
    r2184 	 -1.0
    EX_csn_LPAREN_e_RPAREN_ 	 3
    MCPST 	 -1.0
    RE3258R 	 2.0
    RE3511R 	 -1.0
    RE3136C 	 3.0
    RE3095X 	 1
    2HBO 	 3.0
    r2287 	 -1.0
    EX_HC02172_LPAREN_e_RPAREN_ 	 1
    FUT17g 	 0.0
    EX_3hibupglu_S_LPAREN_e_RPAREN_ 	 1
    DHFtl 	 1
    r1181 	 3.0
    EX_7bhglzglc_LPAREN_e_RPAREN_ 	 1
    FVStu 	 0.0
    RE2624X 	 -1.0
    r2325 	 1.0
    XOLTRIOLtm 	 1
    r1499 	 -1.0
    FACOAL120i 	 1
    PROSTGF2t 	 0.0
    RE1096C 	 1
    6HLVSTitr 	 1
    H8MTer_U 	 0.0
    DNDPt7m 	 -1.0
    AM1ACSitr 	 1
    35DHPVSthc 	 0.0
    FA141ACPH 	 -1.0
    r1649 	 3.0
    FAOXC5OHc 	 1.0
    r1962 	 0.0
    FUT32g 	 1.0
    GQ1Bte 	 1
    r1717 	 0.0
    r1872 	 0.0
    PRGSTRNt 	 1
    HC00342te 	 1
    S2L2FN2M2MASNt 	 1
    PLA2_2 	 1.0
    r0437 	 1
    RE0828C 	 -1.0
    r1502 	 -1.0
    RE2318R 	 3.0
    GAO1g 	 3.0
    r0641 	 1.0
    RE2533C 	 1
    RE3520C 	 1
    FAOXC184C164m 	 3.0
    RE3076X 	 1
    RE3084X 	 0.0
    RE3440R 	 0.0
    r1544 	 3.0
    ADPCOAPTE 	 -1.0
    RE2861C 	 1
    RE0936E 	 0.0
    P45019A1r 	 -1.0
    r1955 	 0.0
    ST6GALNAC24 	 3.0
    r0318 	 0.0
    EX_estradiol_LPAREN_e_RPAREN_ 	 1
    r2343 	 0
    ACMPGLUitr 	 1
    C162OHc 	 1.0
    34DHOXPEGt 	 1
    RE1835C 	 0
    CRNrtx 	 1
    EX_lgnc_LPAREN_e_RPAREN_ 	 1
    MG1er 	 0.0
    UAG2EMAi 	 0.0
    Rtotaltp 	 1
    M4BTAer 	 -1.0
    UGT1A9r 	 -1.0
    r1889 	 0.0
    MANt1r 	 3.0
    C160CPT2 	 0.0
    NACHEXA7ly 	 3.0
    r1322 	 1
    RE1811R 	 0
    MMSAD3m 	 0.0
    AM1CGLChr 	 1
    RDH2 	 -1.0
    NMPTRCOX 	 0.0
    FMNAT 	 -1.0
    35CGMPtn 	 1
    H6_ET2er 	 2.0
    ACNAM9PL2 	 3.0
    DMHPTCRNCPT1 	 1
    SIAASE3ly 	 -1.0
    EX_mal_L_LPAREN_e_RPAREN_ 	 3
    DNDPt12m 	 -1.0
    FAOXC4C4DCc 	 1.0
    3ISPVSteb 	 -1.0
    EX_7hpvs_LPAREN_e_RPAREN_ 	 1
    ASNS1 	 2.0
    r2209 	 -1.0
    r1061 	 1
    r2411 	 -1.0
    ACACT4p 	 3.0
    PS_HStg 	 1
    PVSitr 	 1
    RE2048N 	 2.0
    PHETHPTOX2 	 0
    PROAKGOX1r 	 3.0
    r0871 	 1
    HIStN1 	 -1.0
    FAOXC11BRC9BRx 	 0.0
    NACASPAH 	 -1.0
    EX_1ohmdz_LPAREN_e_RPAREN_ 	 1
    GALASE6ly 	 0
    RE1266C 	 1
    C18OHc 	 1.0
    r1770 	 0.0
    FUT93g 	 -1.0
    PIK5n 	 1
    PI45P5Pn 	 1
    EX_ser_D_LPAREN_e_RPAREN_ 	 1
    14HMDZALThr 	 1.0
    RE3532C 	 1
    DHEASt 	 0.0
    25VITD3Hm 	 1.0
    C160CRNt 	 -1.0
    CYSTAm 	 3.0
    BAAT5x 	 1
    RADH 	 1
    r1459 	 1
    RE3583X 	 1.0
    M14NTg 	 -1.0
    NFDLACtep 	 1
    r1975 	 0.0
    ACN23ACNGALGBSIDEte 	 1
    LNELDCCPT1 	 1.0
    NACHEXA17ly 	 3.0
    FAOXC122_3Z_6Zx 	 3.0
    C141OHe 	 1
    UGT1A1r 	 0
    TYRATB0tc 	 3.0
    CSPG_Etly 	 1
    PI5P4Kn 	 1
    2HATVLACthc 	 0.0
    ACNMLr 	 1
    IDHPOXOX2b 	 -1.0
    RE2871C 	 1
    r2236 	 -1.0
    RE0935E 	 0
    r2057 	 -1.0
    RE3456C 	 1
    ACS 	 3.0
    r2309 	 -1.0
    3HXKYNDCL 	 0.0
    STCOATxc 	 -1.0
    RE2974G 	 0.0
    NTD5l 	 2.0
    r1966 	 0.0
    r2111 	 2.0
    12HTACRtu 	 1
    NTD9l 	 2.0
    RE1700C 	 1
    GLC3MEACPhr 	 0
    3HPVSthc 	 0.0
    SOAT12 	 1.0
    LEUTA 	 1.0
    KHte 	 1
    r1143 	 1
    RE2864C 	 1
    r0122 	 -1.0
    S4TASE3ly 	 -1.0
    10FTHF5GLUtm 	 1
    UMPKn 	 3.0
    NTD2l 	 2.0
    EX_chtn_LPAREN_e_RPAREN_ 	 1
    r1515 	 3.0
    RE2404R 	 0
    r2296 	 -1.0
    r1078 	 1
    P4502C92 	 0.0
    HEXCt 	 1
    CYSO 	 -1.0
    S6TASE17ly 	 3.0
    LST4EXPthc 	 1
    EX_eaflatoxin_LPAREN_e_RPAREN_ 	 1
    CBPSam 	 0.0
    r1068 	 1
    PI45P3Kn 	 1
    EX_leugly_LPAREN_e_RPAREN_ 	 1
    r1852 	 0.0
    LEUKTRE4t 	 0
    RE3195M 	 -1.0
    RE2851C 	 1
    r1677 	 0.0
    RDH1 	 -1.0
    RE1099C 	 0.0
    r0385 	 -1.0
    r1432 	 1.0
    PTVSTLACitr 	 1
    CHOLATEt 	 -1.0
    r2259 	 -1.0
    r1786 	 0.0
    RE2625M 	 0.0
    ACNACNGAL14ACGLCGALGLUSIDEtg 	 1
    r1514 	 3.0
    CYOR_u10m 	 0
    CITMCOAHm 	 1
    RE3012C 	 3.0
    NACHEX20ly 	 3.0
    RE3072X 	 0.0
    RETHe 	 1
    EX_pi_LPAREN_e_RPAREN_ 	 3
    3HPVStep 	 1
    FAOXC164C165m 	 3.0
    RE1516M 	 0.0
    C226CRNt 	 1
    S4T5g 	 1
    RE3087X 	 1.0
    EX_aqcobal_LPAREN_e_RPAREN_ 	 1
    5DHFtl 	 1
    r1914 	 0.0
    BIDGLCURr 	 0
    TAUPAT1c 	 0.0
    r2484 	 -1.0
    RE3163C 	 1
    DHCHOLESTANATEtm 	 1
    MI14PP 	 3.0
    r0909 	 1
    AG13T6g 	 3.0
    DM_dgpi_prot_hs_r_ 	 1
    PA_HStg 	 1
    RSVLAChv 	 1
    EX_coumarin_LPAREN_e_RPAREN_ 	 1
    C12DCACT 	 1
    EX_val_L_LPAREN_e_RPAREN_ 	 3
    RE1819M 	 0
    EX_3hpvs_LPAREN_e_RPAREN_ 	 3
    CARVEOLte 	 1
    r2501 	 -1.0
    ISOLVSTAChep 	 1
    r0782 	 1.0
    PI3P3Pn 	 1
    P45017A2r 	 0.0
    RE3421R 	 1
    PE_HStg 	 1
    RE1809R 	 0
    LSTNM4tev 	 1
    STS1r 	 0.0
    MTHGXLt 	 1
    r0698 	 3.0
    RE3506C 	 1
    C142CPT1 	 1.0
    r1654 	 3.0
    r1827 	 0.0
    r2183 	 -1.0
    FACOAL2251 	 3.0
    O2Stx 	 1
    EX_gmp_LPAREN_e_RPAREN_ 	 1
    CYTDtn 	 1
    CSNAT2x 	 1.0
    CYTK14n 	 3.0
    PROFVShc 	 1
    PAN4PP 	 1
    BTNDm 	 -1.0
    r1604 	 3.0
    r2051 	 -1.0
    RE1517M 	 0.0
    EX_ptvst_LPAREN_e_RPAREN_ 	 1
    FAS120COA 	 0.0
    RNDR3 	 3.0
    r1755 	 0.0
    r2256 	 -1.0
    r1864 	 0.0
    EX_thym_LPAREN_e_RPAREN_ 	 1
    EX_deoxfvs_LPAREN_e_RPAREN_ 	 1
    MGACONm 	 1
    r1645 	 3.0
    RE3447X 	 1.0
    TCHOLAt 	 1.0
    DM_pe_hs_LPAREN_r_RPAREN_ 	 3
    4HPROLTASCT1 	 0.0
    r1173 	 3.0
    ACGALFUCGALACGALFUCGALACGLCGAL14ACGLCGALGLUSIDEtg 	 1
    r1434 	 1
    RE1952R 	 2.0
    RE3288C 	 1
    GULLACter 	 1
    RE3422C 	 1
    FAOXC183_6Z_9Z_12Zm 	 3.0
    NTD6l 	 2.0
    PDE4 	 0.0
    ADK1 	 3.0
    GLYt4 	 3.0
    10FTHF6GLUtm 	 1
    r1671 	 0.0
    RE1473C 	 1
    RE2973G 	 3.0
    PUNP7 	 3.0
    PRPNCOAHYDm 	 3.0
    FAOXC182806m 	 3.0
    EX_c8crn_ 	 1
    DTDPtn 	 1
    RE1134R 	 -1.0
    r2234 	 -1.0
    EX_atvacid_LPAREN_e_RPAREN_ 	 1
    FAOXC103C102m 	 3.0
    EX_HC02187_LPAREN_e_RPAREN_ 	 1
    EX_xoltri24_LPAREN_e_RPAREN_ 	 1
    RE3502C 	 3.0
    LNLNCACRNt 	 -1.0
    3MOX4HOXPGALDOX 	 3.0
    HSPGt 	 1
    P45027A1m 	 0.0
    IDHPOXOXb 	 -1.0
    DHORTS 	 0.0
    DCSPTN1CPT2 	 0.0
    EX_23cump_LPAREN_e_RPAREN_ 	 1
    RE2985X 	 0.0
    CSPG_At 	 1
    NACHEX18ly 	 3.0
    GALGT3 	 0.0
    r2233 	 -1.0
    PYDXK 	 3.0
    EHGLATm 	 3.0
    PVSHtu 	 1
    5THFtm 	 1
    FACOAL204 	 3.0
    r2525 	 -1.0
    r0426 	 3.0
    B3GNT313g 	 0.0
    RE1531M 	 3.0
    2HATVACIDOXDhc 	 1.0
    NRPPHRVESSEC 	 0.0
    UMPK4 	 3.0
    FUCFUC12GAL14ACGLCGALGLUSIDEtg 	 1
    ACNGALACGLCGAL14ACGLCGALGLUSIDEtg 	 1
    5HLTDL 	 0.0
    FAOXC5C3x 	 -1.0
    CSBPASEly 	 1
    3AIB_Dtm 	 1
    SEBACIDTD 	 1
    SMVACIDATPteb 	 -1.0
    EX_strch2_LPAREN_e_RPAREN_ 	 3
    FAOXC16OHC16r 	 0.0
    r1446 	 0.0
    URIt5le 	 0.0
    DM_oretn_n_ 	 1
    RE3010R 	 1
    3HIBUPGLUC_Sitr 	 1
    6MELVSTitr 	 1
    PGS 	 2.0
    RE3228C 	 1
    r2213 	 -1.0
    r0466 	 1
    DEOXFVShc 	 -1.0
    r1617 	 0.0
    ARGDCm 	 -1.0
    RE3152R 	 1
    r1002 	 1
    GALASE17ly 	 0
    FAH2 	 -1.0
    r1538 	 -1.0
    EX_HC00822_LPAREN_e_RPAREN_ 	 1
    EX_hexc_LPAREN_e_RPAREN_ 	 1
    ORNALArBATtc 	 -1.0
    DHPR 	 2.0
    EX_udpg_LPAREN_e_RPAREN_ 	 1
    RE3560C 	 0
    RE2373C 	 1
    TXASr 	 -1.0
    r1917 	 0.0
    RE1699C 	 1
    r2001 	 -1.0
    FAOXC15C13m 	 3.0
    r2273 	 -1.0
    ALAATB0tc 	 3.0
    RE1632C 	 1
    EX_lstn_LPAREN_e_RPAREN_ 	 1
    FAOXC13C11m 	 3.0
    FVSitx 	 1
    EX_leu_L_LPAREN_e_RPAREN_ 	 3
    FACOAL245_2 	 3.0
    RE3019R 	 0.0
    r0642 	 3.0
    DM_1a25dhvitd3_LPAREN_n_RPAREN_ 	 1
    CRGLZABCt 	 1
    RE3041N 	 2.0
    GLCtly 	 1
    PAIL_HStn 	 1
    r0644 	 -1.0
    r1641 	 3.0
    DCSPTN1COAtx 	 1
    AM1A4NCSteb 	 1
    6CSMVhep 	 1.0
    N2M2NMASNtly 	 1
    15DMTitr 	 1
    RE2909M 	 3.0
    RE0691C 	 1
    r1999 	 -1.0
    r1893 	 0.0
    AG13T11g 	 3.0
    r0361 	 3.0
    STS2 	 0.0
    ALDD20xm 	 3.0
    DHFtm 	 1
    RE3410C 	 1
    FAOXC6DCC4DCx 	 -1.0
    RE2319X 	 3.0
    r1832 	 0.0
    ABUTt4_2_r 	 -1.0
    FAOXC142_5Z_8Zx 	 -1.0
    3HXKYNOXDA 	 3.0
    r0983 	 -1.0
    LRAT1 	 1
    CHOLt4 	 -1.0
    RE3194M 	 3.0
    DNDPt45m 	 -1.0
    EX_c51crn_ 	 1
    BTNDe 	 -1.0
    r0997 	 1.0
    BCRNe 	 1
    RE2995M 	 3.0
    CHSTEROLtg 	 3.0
    r1657 	 3.0
    CYSACMPAChc 	 1
    EX_sucr_LPAREN_e_RPAREN_ 	 3
    HSD3B3r 	 -1.0
    r0739 	 3.0
    PA_HSter 	 1
    RE3446C 	 -1.0
    r2146 	 -1.0
    H2O2tm 	 1
    RE3075C 	 1
    PTVSTM3eb 	 -1.0
    ADNt5le 	 0.0
    DESAT18_8 	 0.0
    FAOXC10C10OHm 	 1
    5OHFVSteb 	 1
    OAGD3tg 	 1
    EX_no2_LPAREN_e_RPAREN_ 	 3
    RE1523M 	 3.0
    LPCHOLt 	 1
    RE3449C 	 1
    CLPNDCPT1 	 1.0
    r2516 	 2.0
    r2268 	 -1.0
    RE0512M 	 3.0
    GSNt 	 -1.0
    HESTRATRIOLtr 	 1
    PCHOLPm_hs 	 2.0
    r0822 	 1
    ACGALFUCGALACGALFUCGALACGLCGAL14ACGLCGALGLUSIDEte 	 1
    FAOXTC122TC101m 	 1.0
    G12MT1_U 	 1
    CYSASNNaEx 	 0.0
    NTD2m 	 -1.0
    FATP6t 	 0.0
    C4DCCACT 	 -1.0
    RE2526C 	 1
    CRNt 	 0.0
    AG13T7g 	 3.0
    PI45P5P 	 1.0
    ALR2 	 3.0
    ORPT 	 3.0
    UDPGLCter 	 1
    r1808 	 0.0
    LVSTOXD6Hhep 	 1.0
    ADA 	 1.0
    DOPASULT 	 2.0
    r1940 	 0.0
    CYTK6n 	 3.0
    r0842 	 1
    r2515 	 1
    CYSB0AT3tc 	 -1.0
    TETHEX3t 	 1
    EX_4mop_LPAREN_e_RPAREN_ 	 3
    RE3341M 	 1
    UGLT 	 -1.0
    ADSL1 	 3.0
    sink_decdicoa_LPAREN_c_RPAREN_ 	 1
    GLACO 	 3.0
    LVSTPGPtu 	 -1.0
    IDOAASE3ly 	 -1.0
    CVM23GLUChc 	 0
    GLNLASEer 	 1
    EX_gua_LPAREN_e_RPAREN_ 	 3
    RE0922C 	 -1.0
    r0514 	 1
    MMTSADm 	 1
    r1481 	 3.0
    r2142 	 -1.0
    r1076 	 1
    r1754 	 0.0
    MCDp 	 -1.0
    r1903 	 0.0
    DM_core8_g_ 	 1
    RE3244C 	 1
    EX_pro_D_LPAREN_e_RPAREN_ 	 1
    RE1905C 	 1
    RE1587L 	 1
    RE1099L 	 0.0
    RE2857C 	 1
    r1573 	 3.0
    L_LACt4r 	 -1.0
    GLNATB0tc 	 3.0
    r2280 	 -1.0
    RE3044N 	 2.0
    DM_dctp_n_ 	 3
    r1883 	 0.0
    RE2383R 	 0.0
    COAtn 	 1
    CSm 	 3.0
    RE1818M 	 3.0
    RE3343M 	 3.0
    r1088 	 1
    MERACMPthc 	 0.0
    EX_acmana_LPAREN_e_RPAREN_ 	 3
    UREAt 	 0.0
    PI5P4K 	 0.0
    AKGtp 	 1
    THRGLNexR 	 -1.0
    RE3339M 	 3.0
    r0355 	 3.0
    ACMPGLUTTRsc 	 1
    RE1050N 	 -1.0
    r1116 	 1
    COAtl 	 1
    MI145PP 	 1.0
    RE2632M 	 0.0
    RE2660C 	 1
    r1600 	 3.0
    EX_6htststerone_LPAREN_e_RPAREN_ 	 1
    r0728 	 1.0
    r1523 	 3.0
    SPODMm 	 3.0
    IDOURtly 	 -1.0
    RETFAt2 	 1
    RE3013C 	 1
    RE3575X 	 3.0
    TRPt 	 0.0
    ALAGLNNaEx 	 0.0
    EX_tmdm1_LPAREN_e_RPAREN_ 	 1
    RE3123R 	 1
    1a25DHVITD3TRn 	 1
    MELATNOX 	 0.0
    RE1828C 	 1
    CBL2OR 	 1
    THYMDt1 	 -1.0
    CMPSAS 	 0.0
    r0399 	 0
    XOLDIOLONEt 	 1
    APRGSTRNte 	 1
    r1820 	 0.0
    DNDPt26m 	 1
    CRTSLtm 	 1
    PPA 	 3.0
    EX_hista_LPAREN_e_RPAREN_ 	 1
    PROSTGE2t 	 2.0
    r0954 	 2.0
    12HTACRitr 	 1
    FAOXC101C8m 	 1.0
    RE2155C 	 1
    RE3336M 	 3.0
    FADtm 	 1
    RE2958C 	 1
    RE3150C 	 1
    r2485 	 -1.0
    EX_9_cis_retfa_LPAREN_e_RPAREN_ 	 1
    IBUPGLUCtchep 	 1
    56EPPVSteb 	 -1.0
    EX_34dhphe_LPAREN_e_RPAREN_ 	 1
    ALAyLATthc 	 0.0
    EX_o2s_LPAREN_e_RPAREN_ 	 3
    EX_gthrd_LPAREN_e_RPAREN_ 	 3
    r2193 	 -1.0
    r1909 	 0.0
    RE3630C 	 1
    CERT2gt 	 -1.0
    RE2252C 	 1
    RE3230C 	 1
    r2039 	 -1.0
    NTD7l 	 2.0
    r1904 	 0.0
    P45027A16m 	 0.0
    CITMCOALm 	 1
    FAOXC204_5Z_8Z_11Z_14Zx 	 -1.0
    LEUtec 	 -1.0
    THMt3 	 3.0
    PSERT 	 3.0
    EX_5ohfvs_LPAREN_e_RPAREN_ 	 1
    CKc 	 3.0
    PTVSTtu 	 0.0
    r0386 	 -1.0
    r1920 	 0.0
    H3ETer 	 0.0
    C2M26DCOAHLx 	 3.0
    GLYCtm 	 1
    r0523 	 0
    r1575 	 3.0
    ACN23ACNGALGBSIDEtg 	 1
    r2397 	 -1.0
    TSACGLUCtev 	 1
    AM4NCSteb 	 1
    NRVNCCOAtxc 	 -1.0
    P45011B11m 	 -1.0
    RE2677G 	 1
    B3GNT39g 	 0.0
    FACOAL2042 	 3.0
    CVM1GLUChc 	 0
    THRD_L 	 0.0
    RE3624X 	 0.0
    CRVSM24tev 	 1
    ABTArm 	 1.0
    M16N4Tg 	 1
    Rtotaltl 	 1
    P45019A2r 	 -1.0
    RE3191M 	 3.0
    ASNtN1 	 -1.0
    r1946 	 0.0
    r2094 	 2.0
    LGNCFATPtc 	 0.0
    PEAMNO 	 3.0
    UGT1A10r 	 0
    H4ETer 	 0.0
    B3GNT32g 	 0.0
    P45027A14m 	 0.0
    3HAO 	 -1.0
    DADAe 	 1.0
    PROSTGD2t 	 1.0
    r0731 	 1.0
    ACSMCT1 	 -1.0
    r0331 	 1
    AMPDA 	 2.0
    DOLPMT4_Uer 	 1
    EX_nadp_LPAREN_e_RPAREN_ 	 1
    S6T12g 	 0.0
    FUT911g 	 -1.0
    B3GNT12g 	 3.0
    Htm 	 3.0
    r1176 	 3.0
    RE0581C 	 1
    RE1627C 	 1
    FUT94g 	 -1.0
    DM_dem2emgacpail_prot_hs_r_ 	 1
    FOLABCCte 	 3.0
    RE1818R 	 3.0
    M7MASNBterg 	 1
    RE2235R 	 0.0
    XOL27OHtm 	 1
    r1147 	 0.0
    r2382 	 0.0
    biomass_RNA 	 3
    GLCAT8g 	 2.0
    HEXCCOAtx 	 1
    PDXPP 	 0
    EPOXTAChr 	 1
    r2150 	 -1.0
    RE3147C 	 1
    S6T8g 	 0.0
    EX_dhf_LPAREN_e_RPAREN_ 	 1
    r1748 	 0.0
    r1588 	 3.0
    r1961 	 0.0
    GLUVESSEC 	 -1.0
    DIGALSIDEtg 	 1
    NADH2_u10m 	 0
    RE3022C 	 1
    BTNt4i 	 -1.0
    RE3020C 	 1
    G1M6MASNB1terg 	 1
    RE3011R 	 1
    r2534 	 -1.0
    RE2910X 	 3.0
    RE3633C 	 1
    RE2633C 	 1
    RE3156X 	 3.0
    C6DCc 	 1.0
    r1528 	 1
    ALADGLYexR 	 -1.0
    6DHFtl 	 1
    r1762 	 0.0
    RE0920C 	 -1.0
    GLNB0AT3tc 	 -1.0
    XOLTRI24te 	 1
    r2013 	 -1.0
    RE1834X 	 0.0
    r2284 	 -1.0
    RE2410N 	 0.0
    5ADTSTSTERONEGLCtr 	 1
    FAOXC260240x 	 -1.0
    EX_taur_LPAREN_e_RPAREN_ 	 1
    PItg 	 1
    CLPNDCOAtx 	 1
    RE2638C 	 0.0
    PDE4g 	 -1.0
    RE3448C 	 3.0
    6HSMVACIDhep 	 0.0
    FA140ACPH 	 -1.0
    ILEB0AT3tc 	 -1.0
    RE1836M 	 3.0
    34DHXMANDACOX_NADP_ 	 3.0
    DHCRD1 	 2.0
    LYStiDF 	 3.0
    r0390 	 0.0
    DOLPGT3_Ler 	 0.0
    GARFT 	 3.0
    r1531 	 1.0
    EX_3ispvs_LPAREN_e_RPAREN_ 	 1
    ADNK4 	 1
    r0776 	 3.0
    HSD17B9r 	 2.0
    r2059 	 -1.0
    CSAtu 	 -1.0
    EX_25hvitd3_LPAREN_e_RPAREN_ 	 1
    r1519 	 3.0
    ALOX15 	 1.0
    PRPNCOAHYDx 	 1.0
    ALATHRNaEx 	 0.0
    FAOXC183_9Z_12Z_15Zm 	 3.0
    RE2908M 	 3.0
    EX_nac_LPAREN_e_RPAREN_ 	 3
    PROD2 	 -1.0
    LTA4H 	 3.0
    RE3345M 	 3.0
    r1651 	 3.0
    RE0916G 	 0.0
    1331TACRtev 	 1
    r2310 	 -1.0
    r2102 	 1
    FUCASE2ly 	 3.0
    4HOXPACDOX_NADP_ 	 3.0
    5AOPtm 	 1
    r2482 	 -1.0
    TMDK1 	 3.0
    ALACYSNaEx 	 0.0
    BAMPPALDOX 	 3.0
    r1879 	 0.0
    r2182 	 -1.0
    VLCSp 	 0.0
    EX_desfvs_LPAREN_e_RPAREN_ 	 1
    EX_tststerone_LPAREN_e_RPAREN_ 	 1
    r2317 	 1.0
    r2079 	 2.0
    CSASULPtev 	 1
    TRIODTHYt 	 1.0
    RE0927R 	 0
    r2253 	 -1.0
    r0672 	 0.0
    RE2474R 	 1.0
    FOLOAT1tc 	 -1.0
    EX_11_cis_retfa_LPAREN_e_RPAREN_ 	 1
    r2085 	 2.0
    LPS 	 1.0
    15DMTtu 	 1
    r2508 	 1
    S23T3g 	 0.0
    DOLK_U 	 1
    H6_ETer 	 0.0
    RE0571C 	 1
    C14STRr 	 0.0
    RE3503N 	 2.0
    xmpt 	 1
    DLNLCGCPT2 	 0.0
    OXYP1CONJ 	 1
    RE3496N 	 3.0
    AGPex 	 1
    RE3392M 	 3.0
    GAL3ST12 	 -1.0
    r0712 	 0.0
    r1925 	 0.0
    HSAT2ly 	 1
    RE1899C 	 1
    r2164 	 -1.0
    r1662 	 3.0
    HEX7 	 3.0
    EX_gly_LPAREN_e_RPAREN_ 	 3
    NACHEXA14ly 	 3.0
    TETPENT3COAtx 	 1
    EX_dad_2_LPAREN_e_RPAREN_ 	 3
    H2Otm 	 -1.0
    NADtm 	 1
    M8MASNterg 	 1
    LGNCCOAtx 	 1
    RE3564C 	 3.0
    RE2405C 	 -1.0
    CHOLPtl 	 1
    r2003 	 -1.0
    RE2445E 	 0.0
    EX_CE4881_LPAREN_e_RPAREN_ 	 1
    FACOAL1822 	 3.0
    r1301 	 1
    ST8SIA12 	 -1.0
    Uritm 	 -1.0
    NPTHLte 	 1
    DADA 	 1.0
    r1464 	 1
    CTPS1 	 3.0
    SFCYSc 	 1
    FAOXC14DCC12DCx 	 -1.0
    r0445 	 1
    NS26T2g 	 2.0
    EX_hyptaur_LPAREN_e_RPAREN_ 	 1
    CARIBUP_Sitr 	 1
    RE2992X 	 1.0
    RE0566C 	 1
    FACOAL1813 	 3.0
    DGAT 	 0.0
    EX_4nph_LPAREN_e_RPAREN_ 	 1
    RE3250M 	 3.0
    RE3561M 	 1
    EX_5homeprazole_LPAREN_e_RPAREN_ 	 1
    r1299 	 1
    RE3162R 	 1
    C161CPT2 	 0.0
    r2061 	 -1.0
    PROGLYPEPT1tc 	 -1.0
    MALSO3tm 	 -1.0
    CYTDt2r 	 0.0
    r0825 	 1
    r2205 	 -1.0
    r1826 	 0.0
    DOCO13ECOAtxc 	 -1.0
    SPR 	 0.0
    ACODA 	 -1.0
    EX_fvs_LPAREN_e_RPAREN_ 	 1
    FPGS4m 	 -1.0
    PIt7 	 0.0
    EX_rib_D_LPAREN_e_RPAREN_ 	 3
    r2288 	 -1.0
    ASPGLUm 	 2.0
    XOLTRI27te 	 1
    MEOHtr 	 1
    RTOT_2 	 1
    S4TASE1ly 	 -1.0
    VLCS2p 	 0.0
    r2195 	 -1.0
    r1005 	 1
    LINKDEG3ly 	 1
    THRS 	 0.0
    PPMI12346Ptn 	 1
    SBCOAACOTx 	 -1.0
    3ISPVStep 	 1
    P4502A6 	 -1.0
    C60CRNt 	 -1.0
    RN0031X 	 0.0
    SMVACIDhep 	 3.0
    PSDm_hs 	 0
    UMPK4n 	 3.0
    EX_glysar_LPAREN_e_RPAREN_ 	 1
    ADRNCPT2 	 0.0
    HEXDICOAACBP 	 3.0
    GCHOLAtx 	 1
    EX_fuc13galacglcgal14acglcgalgluside_hs_LPAREN_e_RPAREN_ 	 1
    RE3104R 	 1
    r1884 	 0.0
    RE2875C 	 1
    OIVD2m 	 -1.0
    sink_c81coa_LPAREN_c_RPAREN_ 	 1
    FUCGALGBSIDEte 	 1
    GAMt1r 	 3.0
    r0809 	 1
    BDHm 	 0.0
    DESAT22_2p 	 1
    PGCD 	 2.0
    DPGase 	 3.0
    r2366 	 0
    APPNNte 	 1
    S6T22g 	 0.0
    CLPNDCRNt 	 -1.0
    GLCAT5g 	 -1.0
    r0393 	 3.0
    EX_pydxn_LPAREN_e_RPAREN_ 	 3
    r1518 	 3.0
    biomass_lipid 	 3
    r2027 	 -1.0
    RE2814R 	 1
    RE2130C 	 1
    OCTDECCPT1 	 1.0
    GLCAT4g 	 -1.0
    HPDCAt 	 1
    RE1134M 	 0.0
    56DHPVShc 	 1.0
    RE0456N 	 1.0
    RE2410C 	 1
    EX_co2_LPAREN_e_RPAREN_ 	 3
    3HPCOAHYD 	 3.0
    RE0912C 	 0
    GLYKm 	 0.0
    EX_sphs1p_LPAREN_e_RPAREN_ 	 1
    r0156 	 1
    MCLOR 	 3.0
    3SALACBOXL 	 1.0
    AK2LGCHOLt 	 1
    FAOXC12DCC10DCx 	 -1.0
    GALK 	 0.0
    DOLPH_Uer 	 1
    4HMDGLUCitr 	 1
    ACNACNGALGBSIDEtg 	 1
    FUCFUCFUCGALACGLCGAL14ACGLCGALGLUSIDEte 	 1
    MDZGLChr 	 1
    EX_nrvnc_LPAREN_e_RPAREN_ 	 1
    LEUATB0tc 	 3.0
    P4508B13r 	 -1.0
    RE3393M 	 3.0
    r2265 	 -1.0
    UDPG1P 	 1
    DNDPt52m 	 -1.0
    LEUKTRC4t 	 0.0
    RETH1 	 1
    EX_avite2_LPAREN_e_RPAREN_ 	 1
    r1819 	 0.0
    TACRtu 	 -1.0
    SPRn 	 0.0
    RE3563X 	 1.0
    RE2149C 	 -1.0
    PPCDC 	 -1.0
    RE2996X 	 -1.0
    r1720 	 0.0
    NNATm 	 0.0
    RE2635C 	 1
    EX_citr_L_LPAREN_e_RPAREN_ 	 1
    6OHFVShc 	 0.0
    RE0908C 	 1
    EX_12htacr_LPAREN_e_RPAREN_ 	 1
    PDHm 	 1.0
    r0502 	 2.0
    GGH_6THFe 	 3.0
    OCOAT1m 	 3.0
    3M4HDXPAC 	 3.0
    P45011B21m 	 -1.0
    ALDD2xm 	 3.0
    DNDPt39m 	 -1.0
    RE1834M 	 3.0
    4HBZCOAFm 	 1
    4HMDGLUCtev 	 1
    EX_tmndnc_LPAREN_e_RPAREN_ 	 1
    MLTG1e 	 -1.0
    AG13T8g 	 3.0
    UDPGALtg 	 2.0
    RE1702C 	 1
    EX_prgstrn_LPAREN_e_RPAREN_ 	 1
    ESTRIOLtr 	 1
    r2067 	 -1.0
    CSPG_Btly 	 1
    RE3404M 	 3.0
    r1455 	 1
    2HATVLACGLUCitr 	 1
    RE0582N 	 1
    FACOAL40im 	 -1.0
    r1871 	 0.0
    EX_dopasf_LPAREN_e_RPAREN_ 	 1
    6EPVSthc 	 0.0
    r0403 	 0.0
    RE3079C 	 1
    glyc3pte 	 1
    CGMPt 	 1.0
    BZtr 	 1
    RE2221C 	 1
    TMACMPitr 	 1
    ASNPHELAT2tc 	 0.0
    r0986 	 1
    RE0922R 	 0
    RE3432X 	 0
    EX_progly_LPAREN_e_RPAREN_ 	 1
    FUCFUC12GAL14ACGLCGALGLUSIDEte 	 1
    r1681 	 0.0
    EX_asn_L_LPAREN_e_RPAREN_ 	 3
    P4504F121r 	 0.0
    EX_biocyt_LPAREN_e_RPAREN_ 	 1
    EX_HC02201_LPAREN_e_RPAREN_ 	 1
    FRDPtc 	 1
    r2406 	 -1.0
    C141OHc 	 1.0
    DHCRD2 	 2.0
    SELCYSLY 	 -1.0
    RE1527C 	 1
    GLYt7_211_r 	 -1.0
    AM1ALCStep 	 1
    r0942 	 3.0
    ADSS 	 3.0
    DDCRNe 	 1
    DINt 	 -1.0
    RE3289R 	 0.0
    GLUB0AT3tc 	 -1.0
    RE0921R 	 0
    FAOXC10080m 	 3.0
    ASNtm 	 1
    r2176 	 -1.0
    PUNP6 	 3.0
    5MTHFt2 	 -1.0
    SPHS1Pte 	 1
    NACASPtm 	 1
    ETHAK 	 3.0
    r1643 	 3.0
    CYSTHRNaEx 	 0.0
    EX_crglz_LPAREN_e_RPAREN_ 	 1
    r2327 	 1.0
    r2362 	 0
    DNDPt58m 	 -1.0
    EX_HC01700_LPAREN_e_RPAREN_ 	 1
    6MELACAChep 	 0.0
    SCP2x 	 3.0
    r1648 	 3.0
    EX_inost_LPAREN_e_RPAREN_ 	 3
    r2000 	 -1.0
    SELCYSTS 	 0.0
    ESTRIOLGLCte 	 3.0
    ANTHte 	 1
    AHEXASEly 	 3.0
    FVSCOAhc 	 1
    RE3168R 	 1
    RE3218L 	 0.0
    CRVSM24teb 	 -1.0
    ATVACIDOATPtu 	 0.0
    IMACTD_m 	 3.0
    Htr 	 1
    OIVD1m 	 -1.0
    RE3235C 	 1
    RE3386M 	 3.0
    EX_fe3_LPAREN_e_RPAREN_ 	 3
    r1630 	 3.0
    NDPK8 	 0.0
    GGLUCT 	 3.0
    r1779 	 0.0
    r1916 	 0.0
    AMETt2m 	 -1.0
    Uritl 	 -1.0
    PGDIr 	 0.0
    FAOXC163C142x 	 1.0
    FAS80COA_L 	 0.0
    SMVLAChep 	 1
    r0127 	 -1.0
    NDPK9 	 0.0
    FAOXC204C184m 	 3.0
    ADNtl 	 -1.0
    RE3161C 	 3.0
    LSTNM1tev 	 1
    1OHMDZtep 	 1
    BDG2HCGHD 	 1
    EX_dtdp_LPAREN_e_RPAREN_ 	 1
    3HSMVACIDhep 	 0.0
    r2228 	 -1.0
    LNELDCCPT2 	 0.0
    ST3GAL31g 	 0.0
    LSTNtd 	 1
    RE3343X 	 0.0
    EX_6hmsmvacid_LPAREN_e_RPAREN_ 	 1
    HACD9m 	 3.0
    FACOAL1812 	 3.0
    EX_q10_LPAREN_e_RPAREN_ 	 1
    r2276 	 -1.0
    EX_tsacmgluc_LPAREN_e_RPAREN_ 	 1
    RE1815C 	 1
    GLNtm 	 1
    TMDPPK 	 1
    FACOAL203 	 3.0
    r0779 	 3.0
    r2420 	 -1.0
    r0870 	 1
    r2120 	 2.0
    RE1134C 	 1
    GALFUCGALACGLCGAL14ACGLCGALGLUSIDEte 	 1
    UDPXYLtg 	 0.0
    PRISTANALtx 	 1
    RE1635C 	 1
    r1564 	 3.0
    r0732 	 3.0
    TYROXDAc 	 3.0
    r2199 	 -1.0
    EX_HC02220_LPAREN_e_RPAREN_ 	 1
    r0656 	 3.0
    ACSOMT 	 -1.0
    r1525 	 3.0
    MAOLNOR 	 3.0
    RE2346C 	 1
    RE2876C 	 1
    B3GNT35g 	 0.0
    DARGOp 	 -1.0
    EX_carveol_LPAREN_e_RPAREN_ 	 1
    31DMTtu 	 1
    ALATA_L 	 0.0
    EX_HC02200_LPAREN_e_RPAREN_ 	 1
    HISTASE 	 0.0
    FAOXC226C205m 	 3.0
    CITRtm 	 -1.0
    r0545 	 3.0
    DESFVSitr 	 1
    OXAHCOtex 	 2.0
    MM6B1ag 	 3.0
    GLCAT3g 	 -1.0
    r1950 	 0.0
    SIAASE 	 -1.0
    ACGALK2 	 1
    r1165 	 1.0
    RE3036N 	 1
    EX_tmd_LPAREN_e_RPAREN_ 	 1
    OXYPR1tehv 	 1
    FAOXC2442246x 	 -1.0
    SPMS 	 0.0
    GALASE15ly 	 0
    r2386 	 0.0
    RE2360N 	 0.0
    HAtly 	 1
    SERTHRNaEx 	 0.0
    MANt4 	 -1.0
    LPS2 	 1.0
    AG13T1g 	 3.0
    RE2524X 	 3.0
    H2O2itr 	 1
    PYDAMtr 	 1
    FDH 	 0.0
    r2304 	 -1.0
    MDZtd 	 1
    RE3307X 	 0
    FAOXC161C161OHm 	 3.0
    RE1531X 	 3.0
    ADEtl 	 -1.0
    C12DCe 	 1
    FPGS2m 	 -1.0
    5OHFVSitr 	 1
    GALt4 	 -1.0
    FACOAL200 	 3.0
    RE0549C 	 1
    EX_HC01361_LPAREN_e_RPAREN_ 	 1
    DDECE1CRNe 	 1
    DM_dgtp_n_ 	 3
    r1800 	 0.0
    EX_glyphe_LPAREN_e_RPAREN_ 	 1
    DHEAtr 	 1
    PRDXl 	 -1.0
    FAOXC141C121x 	 -1.0
    25HVITD3tm 	 1
    PTVSTM13te 	 1
    S6T21g 	 0.0
    AM1A4NCShc 	 1
    NDP8ex 	 3.0
    BUP2 	 -1.0
    SIAT4Bg 	 -1.0
    2HATVACIDGLUChr 	 0
    EX_allop_LPAREN_e_RPAREN_ 	 1
    EX_fmn_LPAREN_e_RPAREN_ 	 3
    ADPGLC 	 3.0
    FAOXC120100m 	 3.0
    PYDXNK 	 3.0
    DM_dttp_m_ 	 1
    MGCHrm 	 0.0
    RE3394M 	 3.0
    RE3476M 	 0
    SALMCOM 	 3.0
    r2444 	 1
    THMMPtm 	 1
    ILEt5m 	 1
    r1867 	 0.0
    2AMACHYD 	 0.0
    35DSMVteb 	 1
    r1805 	 0.0
    IDHPOXOX3 	 -1.0
    r2117 	 2.0
    6MSMVhep 	 1.0
    RE2541L 	 3.0
    SPHMYLNtg 	 1
    GGH_10FTHF5GLUl 	 3.0
    RE2635R 	 -1.0
    VITD2t 	 1
    RE2474C 	 0.0
    LNLCCPT2 	 0.0
    PPDOy 	 3.0
    RE2454M 	 3.0
    DPMVDx 	 -1.0
    RE3175R 	 1
    G14T6g 	 3.0
    EX_24nph_LPAREN_e_RPAREN_ 	 1
    DAGt 	 1
    EX_HC01441_LPAREN_e_RPAREN_ 	 1
    GLCt4 	 3.0
    EX_tetdec2crn_ 	 1
    GTACMPitr 	 1
    UGT1A4r 	 0
    RE2150R 	 0
    r1821 	 0.0
    PCLAD 	 -1.0
    RE1826M 	 0.0
    GTPtn 	 1
    RE2814M 	 1
    GLCNACASE2ly 	 -1.0
    7THFtl 	 1
    RE2880C 	 1
    EX_3octdece1crn_ 	 1
    r1982 	 0.0
    RE2131C 	 1
    r1117 	 1
    r0084 	 3.0
    RE1817M 	 3.0
    r0573 	 -1.0
    r0431 	 0.0
    RE1897C 	 1
    IMACTD 	 3.0
    PUNP4 	 3.0
    RE3490C 	 1
    DIDPtn 	 1
    RE1906R 	 -1.0
    FAOXC181C161m 	 3.0
    FAOXC18C18OHm 	 3.0
    MMCD 	 -1.0
    EX_coa_LPAREN_e_RPAREN_ 	 1
    r0913 	 0.0
    RE2640C 	 -1.0
    B3GNT315g 	 0.0
    r2162 	 -1.0
    FAOXC102C103x 	 -1.0
    RE0916R 	 0.0
    GQ1BALPHAtg 	 1
    LGNCCRNt 	 -1.0
    AGPRim 	 3.0
    TACRDtsc 	 1
    r1716 	 0.0
    r1051 	 1
    PGDI 	 0.0
    r0016 	 0
    ATVACIDtdu 	 1
    FAOXC12C12OHm 	 3.0
    HSD3B13r 	 -1.0
    CRTSTRNtm 	 1
    PDE1 	 2.0
    LCYSTAT 	 3.0
    TMDtd 	 1
    EX_3ddcrn_ 	 1
    r1636 	 3.0
    ADK3m 	 3.0
    GLC3MEACPitr 	 1
    S4T4g 	 1
    3HSMVACIDteb 	 1
    H2CO3Dm 	 -1.0
    RE2995X 	 3.0
    MDHm 	 3.0
    RE1525C 	 1
    EX_prostgd2_LPAREN_e_RPAREN_ 	 1
    BACCL 	 -1.0
    r0777 	 3.0
    EX_fucacngal14acglcgalgluside_hs_LPAREN_e_RPAREN_ 	 1
    RE2768M 	 1
    EX_estradiolglc_LPAREN_e_RPAREN_ 	 1
    AM19CSALThr 	 1.0
    EX_lstnm7_LPAREN_e_RPAREN_ 	 1
    EX_mdzglc_LPAREN_e_RPAREN_ 	 1
    ACCOAL 	 0.0
    r2218 	 -1.0
    RE1447N 	 1
    r2313 	 1.0
    DNDPt20m 	 -1.0
    OAGD3te 	 1
    PUNP3 	 3.0
    r1170 	 1.0
    ALDSTRNtm 	 1
    LSTNM5tev 	 1
    ARACHDCOAtx 	 1
    FAOXC102C81m 	 1.0
    EX_fucacgalfucgalacglcgalgluside_hs_LPAREN_e_RPAREN_ 	 1
    FAOXC122m 	 3.0
    r2023 	 -1.0
    r2035 	 -1.0
    PHEt4 	 3.0
    PCRNtm 	 -1.0
    r0552 	 3.0
    r1865 	 0.0
    TMABADH 	 3.0
    r2279 	 -1.0
    G14T20g 	 3.0
    GALASE16ly 	 0
    ESTRADIOLtr 	 1
    HACD1m 	 3.0
    r2095 	 2.0
    AMETr 	 1
    P450SCC1m 	 1
    RE1836C 	 1
    r2493 	 -1.0
    DURIt 	 -1.0
    CHSTEROLtrc 	 1
    S6TASE20ly 	 3.0
    GALGALFUCFUCGALACGLCGALACGLCGAL14ACGLCGALGLUSIDEtg 	 1
    r1702 	 0.0
    RE1135R 	 0.0
    DMANTIPYRINEte 	 1
    DM_avite2_c_ 	 1
    VALB0AT2tc 	 1
    KCCt 	 3.0
    r1559 	 3.0
    IDHPOXOX4 	 -1.0
    C100CRNt 	 -1.0
    r0826 	 1
    HPDCACRNt 	 -1.0
    EX_HC02210_LPAREN_e_RPAREN_ 	 1
    GALT 	 -1.0
    r0409 	 1
    ASP1DC 	 1.0
    EX_34hpp_ 	 3
    AHCYStr 	 1
    r1855 	 0.0
    r1983 	 0.0
    H2O2syn 	 3.0
    FACOAL160i 	 3.0
    HS1ly 	 1.0
    FUCFUCGALACGLCGALGLUSIDEte 	 1
    RE3534R 	 3.0
    AATAi 	 0.0
    EX_dhdascb_LPAREN_e_RPAREN_ 	 3
    TETDECE1CRNe 	 1
    ADNtm 	 -1.0
    EX_14hmdz_LPAREN_e_RPAREN_ 	 1
    SEAHCYSHYD 	 3.0
    CMPACNAtg 	 3.0
    TMLYSter 	 1
    EX_3hdececrn_ 	 1
    O2ter 	 1
    COUMARINte 	 1
    r0961 	 0.0
    CSNATer 	 -1.0
    FAOXC11C9m 	 1.0
    ADNt4 	 0.0
    SAMHISTA 	 1.0
    FAOXC122_3Z_6Zm 	 3.0
    r0357 	 3.0
    r2143 	 -1.0
    TSTSTERONEGLCtr 	 1
    SPRMS 	 3.0
    EX_4hatvacid_LPAREN_e_RPAREN_ 	 1
    RE2874C 	 1
    RE1448N 	 1
    RE1309M 	 1.0
    PTDCACRNCPT1 	 1.0
    GQ1Btg 	 1
    34DHPHEt 	 0.0
    TAXOLte 	 1
    RE3150R 	 1
    FE2tm 	 1
    r2062 	 -1.0
    HDCAter 	 1
    AM1CGLCteb 	 1
    CYTK9n 	 3.0
    r1457 	 1
    S6T9g 	 0.0
    r1578 	 3.0
    EX_glcur_LPAREN_e_RPAREN_ 	 3
    HCO3_CLt 	 -1.0
    ALKP 	 0.0
    HEXCCPT1 	 1.0
    RE1100G 	 0.0
    MESCOALm 	 1
    AKR1C41 	 -1.0
    ACALDtm 	 1
    7BHGLZABCt 	 1
    LEULEUPEPT1tc 	 -1.0
    r0363 	 3.0
    DADNt4 	 -1.0
    RE1526C 	 1
    OAGT3te 	 1
    r2270 	 -1.0
    GLYitr 	 1
    PIt8 	 3.0
    RE3391M 	 3.0
    r1606 	 3.0
    RE1829C 	 1
    RE3228R 	 1
    RE2870C 	 1
    SEAHCYStn 	 1
    LEUB0AT2tc 	 1
    r2250 	 -1.0
    TDPDRE 	 1
    ORETNF 	 1
    r1834 	 0.0
    S6T25g 	 2.0
    MMCDm 	 -1.0
    EX_fuc14galacglcgalgluside_hs_LPAREN_e_RPAREN_ 	 1
    r0729 	 1.0
    RE3095C 	 3.0
    DASCBH 	 1
    r1184 	 3.0
    EX_4bhglz_LPAREN_e_RPAREN_ 	 1
    24_25VITD2Hm 	 1.0
    STRDNCCOAtxc 	 -1.0
    FACOAL180i 	 3.0
    PPAer 	 1.0
    r0157 	 1
    ACACT10m 	 3.0
    C8DCc 	 1.0
    ATPasel 	 2.0
    S4TASE2ly 	 -1.0
    ACOAD10m 	 3.0
    r1526 	 3.0
    r1858 	 0.0
    EX_tdchola_LPAREN_e_RPAREN_ 	 1
    MDH 	 3.0
    r0787 	 3.0
    r2425 	 -1.0
    UMPK3 	 3.0
    RE1520M 	 3.0
    EX_prostgh2_LPAREN_e_RPAREN_ 	 3
    C14OHc 	 1.0
    FAOXC182_9Z_12Zm 	 3.0
    RE2914X 	 0.0
    r2105 	 1
    PNTOt5 	 0.0
    FUCACGALFUCGALACGLCGALGLUSIDEtg 	 1
    PEPCKm 	 2.0
    M13N4Tg 	 3.0
    GUAD 	 1.0
    XOLTRI24tc 	 1
    EX_idour_LPAREN_e_RPAREN_ 	 1
    r1154 	 -1.0
    6BHGLZGLCitr 	 1
    EX_pnto_R_LPAREN_e_RPAREN_ 	 3
    MAN6PI 	 -1.0
    EX_estrones_LPAREN_e_RPAREN_ 	 1
    RE3233G 	 -1.0
    EPOXTACteb 	 1
    r1637 	 3.0
    35DSMVitr 	 1
    r2291 	 -1.0
    4NPHSULT 	 2.0
    EX_xmp_LPAREN_e_RPAREN_ 	 3
    EX_4hphac_LPAREN_e_RPAREN_ 	 1
    RE2132C 	 1
    r1512 	 1
    FAOXC205_5Z_8Z_11Z_14Z_17Zm 	 3.0
    RE2658C 	 1
    G1M8MASNterg 	 1
    EX_CE1940_LPAREN_e_RPAREN_ 	 1
    SRTNENT4tc 	 -1.0
    RE3383M 	 3.0
    r1383 	 3.0
    PCRNtc 	 1
    LYSt4 	 3.0
    r1602 	 3.0
    EX_stacmp_LPAREN_e_RPAREN_ 	 1
    r1706 	 0.0
    EPOXTACitr 	 1
    r2379 	 0.0
    r0994 	 1
    ETFQO 	 0.0
    r1162 	 1
    PRISTCOAtcx 	 3.0
    FAOXC227C226m 	 3.0
    GLYCLTDy 	 3.0
    r0760 	 0.0
    11DOCRTSLtr 	 1
    RE3171R 	 1
    FAOXC184m 	 3.0
    ESTRIOLGLCtr 	 1
    AM1ACStep 	 1
    RE1539C 	 1
    CSEPASEly 	 1
    EX_psyltchol_LPAREN_e_RPAREN_ 	 1
    SELCYSTGL 	 0.0
    ARACHCRNt 	 -1.0
    r2483 	 -1.0
    INSTt4_2 	 -1.0
    EX_cys_L_LPAREN_e_RPAREN_ 	 3
    SPHS1Ptr 	 1
    ADKd 	 -1.0
    TMDM5hr 	 0.0
    RE3010X 	 0
    r0651 	 -1.0
    T2M26DCOAHLm 	 3.0
    P45046A1r 	 -1.0
    RE2917M 	 3.0
    r2314 	 1.0
    RE3563M 	 3.0
    FAOXC184_6Z_9Z_12Z_15Zm 	 3.0
    r1064 	 1
    RE2149R 	 0
    GP1CALPHAtg 	 1
    3OHACMPhr 	 0.0
    SMVACIDtev 	 1
    HSAT3ly 	 1
    RE2888N 	 -1.0
    3SPYRSPm 	 1
    GULNter 	 1
    r1595 	 3.0
    ESTRONEGLCt 	 0.0
    ALAt2rL 	 0.0
    r1787 	 0.0
    TSACMSULtev 	 1
    EX_Rtotal_LPAREN_e_RPAREN_ 	 1
    TYRt4 	 3.0
    LEUKTRF4t 	 0
    4HATVACIDhc 	 0.0
    ABO5g 	 -1.0
    r2220 	 -1.0
    41R2A1H12BOOX 	 3.0
    ARGtm 	 -1.0
    EX_fucfuc132galacglcgal14acglcgalgluside_hs_LPAREN_e_RPAREN_ 	 1
    EX_HC02191_LPAREN_e_RPAREN_ 	 1
    AMACRp 	 -1.0
    AM4N9CStev 	 1
    r2082 	 2.0
    RE3597C 	 3.0
    FAOXC225C204m 	 3.0
    THMDt4 	 0.0
    BGLUTDECHOe 	 1
    RE3241R 	 3.0
    ASNTHRNaEx 	 0.0
    r0432 	 0.0
    PMI1346PHn 	 0.0
    5ADTSTSTERONEGLCte 	 3.0
    r0638 	 -1.0
    RE2404C 	 -1.0
    4PYRDX 	 1
    EX_am4ncs_LPAREN_e_RPAREN_ 	 1
    EX_ttdca_LPAREN_e_RPAREN_ 	 3
    C30CPT1 	 1.0
    MDRPD 	 1
    ACOAD9m 	 3.0
    r2071 	 -1.0
    PCHOLHSTDe 	 1
    RE2909X 	 3.0
    IDPtn 	 1
    FAOXC163C164x 	 -1.0
    PI5P3Ker 	 3.0
    B3GNT34g 	 0.0
    r2248 	 -1.0
    r2289 	 -1.0
    CBPS 	 0.0
    RE3234R 	 1
    NDPK4 	 0.0
    r1001 	 1
    RE0567C 	 1
    RE2632C 	 1
    NACHEX19ly 	 3.0
    RE3075X 	 0.0
    ACONTm 	 3.0
    ARSA 	 -1.0
    RE1096R 	 -1.0
    r0464 	 3.0
    r0553 	 3.0
    r1721 	 0.0
    r0836 	 -1.0
    FACOAL170 	 3.0
    RE1582C 	 1
    FUMtm 	 -1.0
    RE3536C 	 1
    r0068 	 0.0
    FACOAL241 	 3.0
    EX_HC00342_LPAREN_e_RPAREN_ 	 1
    ADSL2 	 3.0
    r1497 	 -1.0
    RE1938C 	 1
    MM7Cag 	 3.0
    ANTIPYRENEte 	 1
    r1778 	 0.0
    FAOXC161C141m 	 3.0
    NICRNTtn 	 1
    NACHEX16ly 	 3.0
    FTHFLm 	 3.0
    AM4N9CShc 	 1
    r1495 	 -1.0
    6DHFtm 	 1
    ACETONEt2m 	 2.0
    EX_HC01577_LPAREN_e_RPAREN_ 	 1
    r0990 	 1
    r0472 	 1.0
    CBPter 	 1
    PI345P3P 	 3.0
    PPItr 	 1
    r0028 	 2.0
    r1556 	 3.0
    FAOXC181_9Em 	 3.0
    EX_adpcbl_LPAREN_e_RPAREN_ 	 3
    5HTRPDOX 	 2.0
    RE3301R 	 2.0
    RE2319R 	 3.0
    r0208 	 0
    SERB0AT3tc 	 -1.0
    SBPP1er 	 1.0
    DHEASULT 	 0.0
    UDPGP 	 2.0
    FBP26 	 3.0
    r1474 	 1.0
    CHSTEROLt2 	 1
    RE3121R 	 1
    CYTK14 	 3.0
    r1599 	 3.0
    RE0573N 	 1
    APOCF 	 1
    r2357 	 0
    LEUKTRB4tr 	 1
    INSt 	 -1.0
    GLYGLYCNc 	 3.0
    r0733 	 3.0
    EX_dgchol_LPAREN_e_RPAREN_ 	 1
    r0974 	 1
    GALNACT2g 	 1.0
    r0683 	 1
    AM1CGLCitr 	 1
    ASPDt6 	 1.0
    EX_doco13ac_ 	 1
    FOLt2le 	 1
    ASNCYSNaEx 	 0.0
    RN0022X 	 0.0
    C101CPT1 	 1.0
    RE2078M 	 1
    MERCPLACCYSt 	 1
    r2246 	 -1.0
    CYStec 	 -1.0
    LINKDEG1ly 	 1
    RE3231C 	 1
    M4CET3er 	 0.0
    SLDt 	 1
    FAOXC204_5Z_8Z_11Z_14Zm 	 3.0
    RE3450C 	 1
    r0925 	 1
    DESAT18_6 	 0.0
    RE2405R 	 0
    S6T5g 	 0.0
    RE3147R 	 1
    EX_rbt_LPAREN_e_RPAREN_ 	 1
    IZPN 	 0.0
    UPPN 	 -1.0
    GTHPm 	 3.0
    ARACHDFATPtc 	 -1.0
    1531TALThr 	 1
    HSD17B4x 	 3.0
    DESAT16_2 	 3.0
    r2499 	 -1.0
    FAOXC162C162OHm 	 3.0
    7DHCHSTEROLtr 	 1
    NADS2 	 1.0
    B3GALTg 	 1
    ADRNCPT1 	 1.0
    FUT92g 	 -1.0
    r2144 	 -1.0
    RE0875C 	 1
    ACMPdt 	 1
    RE3232R 	 1
    DCSPTN1t 	 1
    EX_cynt_LPAREN_e_RPAREN_ 	 1
    r2283 	 -1.0
    PDE1g 	 -1.0
    r1318 	 1
    r1773 	 0.0
    NH4t3r 	 -1.0
    ACACtx 	 1
    LSTNtu 	 -1.0
    NACHEX21ly 	 3.0
    4MOPte 	 1
    G6Pter 	 -1.0
    GLXtm 	 1
    DNDPt50m 	 -1.0
    FAOXC9C7m 	 1.0
    INSK 	 1
    RDH3a 	 3.0
    CSPG_Atly 	 1
    LINOFATPtc 	 1
    S6TASE24ly 	 3.0
    G1M7MASNCterg 	 1
    ATVACIDitr 	 1
    ESTRONEGLCtr 	 1
    EX_lnlnca_LPAREN_e_RPAREN_ 	 3
    RE0578M 	 3.0
    LSTNM7hr 	 0
    ILEPHELAT2tc 	 1
    r1678 	 0.0
    RE1628C 	 1
    RE3083X 	 -1.0
    RE3518C 	 1
    RE0568C 	 1
    DOLPMT3_Uer 	 1
    EX_ptvstm13_LPAREN_e_RPAREN_ 	 1
    r0818 	 1
    RE0579C 	 0
    G14T5g 	 3.0
    r2409 	 -1.0
    14HMDZhr 	 1.0
    SUCCCROT 	 1.0
    r1384 	 1
    TXA2te 	 1
    RE3005M 	 3.0
    RE2633R 	 0
    FUT33g 	 1.0
    r0364 	 3.0
    EX_1hmdgluc_LPAREN_e_RPAREN_ 	 1
    OMHPALTD 	 1
    NAHCO3_HCLt 	 -1.0
    RE2111M 	 -1.0
    r2295 	 -1.0
    RE2948C 	 1
    PAPtg 	 1
    4HATVLAChc 	 1
    DPCOAPPe 	 1
    RE1956X 	 0.0
    NRVNCCPT1 	 1.0
    EX_smv_LPAREN_e_RPAREN_ 	 1
    G3PD2m 	 3.0
    3AIBt 	 1
    RE3251C 	 1
    r1570 	 3.0
    r2539 	 1
    EX_avite1_LPAREN_e_RPAREN_ 	 1
    LVSTACIDtu 	 -1.0
    DHCR241r 	 3.0
    RE3252C 	 1
    FALDtm 	 1
    IDOAASE2ly 	 -1.0
    15DMTtep 	 1
    DAGK_hs 	 2.0
    DHORD9 	 -1.0
    SACCD4m 	 -1.0
    r1529 	 3.0
    EX_3deccrn_ 	 1
    r2216 	 -1.0
    RE3446X 	 3.0
    UGT1A2r 	 0
    EX_2425dhvitd3_LPAREN_e_RPAREN_ 	 1
    FAOXC142_5E_8Em 	 3.0
    OCD11COACPT1 	 1
    RE3114R 	 3.0
    r0707 	 -1.0
    UDPXYLter 	 1
    TMDS 	 3.0
    AASAD3m 	 1
    r2145 	 -1.0
    FAOXTC101TC102m 	 3.0
    r0587 	 -1.0
    r1653 	 3.0
    OCDCEAtr 	 -1.0
    ST6GALNAC62 	 -1.0
    r0778 	 3.0
    S2T3g 	 0.0
    NKCCt 	 1.0
    TMDPP 	 3.0
    EX_so3_LPAREN_e_RPAREN_ 	 1
    EX_ppa_LPAREN_e_RPAREN_ 	 3
    RE2912X 	 3.0
    RE2513C 	 3.0
    FVSCOAitx 	 1
    GLYt2r 	 0.0
    LGTHL 	 3.0
    r0800 	 1
    FAOXC101x 	 1.0
    PPOR 	 -1.0
    FUT99g 	 -1.0
    PIt9 	 0.0
    RE0512X 	 3.0
    RE2026C 	 1
    BPNT 	 0.0
    FPGS6m 	 -1.0
    LSTNRATt 	 1
    3HPVSTETCOAhcx 	 -1.0
    EX_camp_LPAREN_e_RPAREN_ 	 1
    RE3132R 	 3.0
    ASCBOX 	 1
    HIShPTtc 	 1.0
    GLYATm 	 -1.0
    BGLUGCHe 	 1
    LNELDCt 	 1
    CARPEPT1tc 	 -1.0
    r2303 	 -1.0
    EX_3ump_LPAREN_e_RPAREN_ 	 1
    RE2562C 	 1
    HSD17B42x 	 3.0
    COAtm 	 2.0
    FPGS5 	 -1.0
    RE0690C 	 -1.0
    RAI4 	 1
    CHLtm 	 1
    RE2858C 	 1
    RE3180C 	 1
    NACHEXA6ly 	 3.0
    GLYAMDTRc 	 3.0
    FUMSO4tm 	 -1.0
    GPIAT 	 1.0
    MM6B1bg 	 3.0
    r2261 	 -1.0
    DEOXFVStev 	 1
    r2108 	 1
    RE3192M 	 3.0
    r0236 	 1
    RE3513R 	 1
    6HMSMVACIDhep 	 0.0
    URIt4 	 0.0
    RE0572N 	 1
    EX_fe2_LPAREN_e_RPAREN_ 	 3
    NADtru 	 1
    r1631 	 3.0
    DNDPt62m 	 -1.0
    EX_dspvs_LPAREN_e_RPAREN_ 	 1
    r1730 	 0.0
    RE2921M 	 3.0
    FUCGAL14ACGLCGALGLUSIDEte 	 1
    FAOXC141_5Zm 	 3.0
    4BHGLZABCt 	 1
    S6TASE13ly 	 3.0
    OIVD3m 	 -1.0
    SALMCOM2 	 0.0
    DIGALSGALSIDEtg 	 1
    r1715 	 0.0
    FAOXC225C226x 	 -1.0
    EX_5adtststeroneglc_LPAREN_e_RPAREN_ 	 1
    RE0344X 	 0.0
    TMNDNCCRNt 	 -1.0
    ICDHy 	 3.0
    RE2999M 	 3.0
    THD1m 	 1.0
    CYTK2 	 3.0
    ACMPGLUthc 	 0.0
    ASCBSVCTtc 	 0.0
    35DHPVStep 	 1
    RE1532M 	 3.0
    RSVSPONhc 	 0.0
    SMVtu 	 -1.0
    ESTRONESt 	 2.0
    RE2852C 	 1
    FAOXC184C164x 	 -1.0
    RE1810M 	 1
    r0989 	 1
    EX_1hibupglu_S_LPAREN_e_RPAREN_ 	 1
    EX_asp_D_LPAREN_e_RPAREN_ 	 1
    r0267 	 1
    r1713 	 0.0
    PSFLIP 	 1.0
    r2342 	 0
    TMDM5OATt 	 1
    ADRNCRNt 	 -1.0
    EX_CE5560_LPAREN_e_RPAREN_ 	 1
    IDOAASE1ly 	 -1.0
    FACOAL224 	 3.0
    DHPM2 	 -1.0
    HTDCACBP 	 3.0
    LEUKTRB4t 	 0
    r2410 	 -1.0
    RE3477C 	 1
    RE2877C 	 1
    RE1630C 	 1
    r2104 	 1
    CRTSTRNt 	 1
    r0165 	 -1.0
    RE2888C 	 3.0
    C180CPT2 	 0.0
    GBSIDEtl 	 1
    EX_prostgi2_LPAREN_e_RPAREN_ 	 1
    CRVSM1hr 	 3.0
    EX_HC02196_LPAREN_e_RPAREN_ 	 1
    TMNDNCCOAtx 	 1
    FADtx 	 1
    HMGCOARc 	 3.0
    RE3015R 	 0.0
    RE1817C 	 1
    DECDICRNe 	 1
    RE1943R 	 1.0
    ILETA 	 1.0
    MM5bg 	 3.0
    LPASE 	 1.0
    r1918 	 0.0
    H7MTer_U 	 1
    EX_phyQ_LPAREN_e_RPAREN_ 	 1
    DSAT 	 1
    RE2273E 	 0.0
    r2121 	 2.0
    ATVACIDMCTtu 	 2.0
    RE3090X 	 3.0
    RE3557M 	 3.0
    PI3P5K 	 0.0
    RE1903R 	 1
    G6PDH2rer 	 0.0
    MTHFR3 	 -1.0
    FAOXC102x 	 1.0
    DOLPGT1_Uer 	 1
    DUTPDPm 	 3.0
    RE2975C 	 1
    r2230 	 -1.0
    r2064 	 -1.0
    DURIPP 	 3.0
    ELAIDt 	 0.0
    RE3340X 	 3.0
    EX_HC01446_LPAREN_e_RPAREN_ 	 1
    RIBt 	 1
    6BHGLZGLChr 	 1
    BALAPAT1tc 	 0.0
    r0139 	 2.0
    FAOXC160140x 	 -1.0
    r2398 	 -1.0
    r2167 	 -1.0
    r0504 	 2.0
    RE2270C 	 1
    EX_2hatvacidgluc_LPAREN_e_RPAREN_ 	 1
    RE1573X 	 3.0
    2MCITt 	 1
    RE3521X 	 3.0
    BILGLCURtr 	 1
    GLBRAN 	 3.0
    EX_crvnc_LPAREN_e_RPAREN_ 	 3
    PAFS 	 1
    RE2910M 	 3.0
    DOPAtu 	 1.0
    SQLSr 	 3.0
    SGALSIDEtl 	 1
    CRNtim 	 -1.0
    H2Otn 	 1
    r1757 	 0.0
    r2073 	 -1.0
    EX_acgalfucgalacgalfucgalacglcgal14acglcgalgluside_hs_LPAREN_e_RPAREN_ 	 1
    RE1812R 	 0
    PPPG9tm 	 1
    NBAHH_ir 	 3.0
    RE0688C 	 -1.0
    r2258 	 -1.0
    S4T1g 	 0.0
    CITt4_4 	 0.0
    LRAT 	 -1.0
    DHAAt1r 	 3.0
    EX_glyleu_LPAREN_e_RPAREN_ 	 1
    RE2398C 	 1
    FAOXC183806m 	 3.0
    NDERSVhc 	 0.0
    RE3346R 	 3.0
    VACCt 	 1
    r2305 	 -1.0
    RE2658R 	 1
    RE2051R 	 3.0
    MI13456Ptn 	 1
    EX_chol_LPAREN_e_RPAREN_ 	 3
    RE2768C 	 1
    r2044 	 -1.0
    3SALAOX 	 1
    LNLNCACPT2 	 0.0
    r1768 	 0.0
    RE2334C 	 1
    PRODt2r 	 0.0
    EX_perillyl_LPAREN_e_RPAREN_ 	 1
    ATP1ter 	 1
    STRDNCt 	 1
    TRIPVSitr 	 1
    Am19CStev 	 1
    CAT2p 	 3.0
    FAOXC225x 	 1.0
    CITtam 	 0.0
    RE3104C 	 1
    HPYRDCm 	 1
    r2239 	 -1.0
    r0633 	 -1.0
    GLNt4 	 3.0
    SACCD3m 	 -1.0
    OMHDOCOSACTD 	 1
    CYTK7 	 3.0
    GGH_10FTHF6GLUl 	 3.0
    RE3106R 	 1
    ASPCTr 	 0.0
    r1850 	 0.0
    FAOXC12DCTc 	 -1.0
    FAOXC226C227m 	 3.0
    RE3571C 	 1
    NACHEXA16ly 	 3.0
    r1686 	 0.0
    EX_mepi_LPAREN_e_RPAREN_ 	 1
    r1923 	 0.0
    r1447 	 0.0
    PCFLOPm 	 -1.0
    ACACT8p 	 3.0
    MCITS 	 1
    HSD17B1 	 0.0
    FAOXC61m 	 3.0
    NRPPHRSULT 	 0
    r2380 	 0.0
    EX_vitd2_LPAREN_e_RPAREN_ 	 1
    4MPTNLte 	 1
    GLCMter 	 1
    GTHOm 	 2.0
    FAOXC22C22DCHYr 	 0.0
    r2007 	 -1.0
    r1854 	 0.0
    LST4EXPitr 	 1
    r0178 	 2.0
    PPAP 	 3.0
    OCDCAFAPMtc 	 3.0
    RE3173C 	 1
    RN0023R 	 2.0
    r1440 	 1
    RE0924R 	 0
    r1953 	 0.0
    AGLPC 	 1
    P4504F122r 	 0.0
    D3AIBTm 	 -1.0
    r0531 	 1
    CYTK13 	 3.0
    MICITDr 	 1
    EX_HC02195_LPAREN_e_RPAREN_ 	 1
    AMY2e 	 0
    G6PPer 	 1.0
    r1109 	 1
    RE0690X 	 0.0
    r0860 	 1
    PE_HStm 	 1
    SPC_HSt 	 1
    r2260 	 -1.0
    r2323 	 1.0
    BHBt 	 1
    13HTACRitr 	 1
    GT1Atg 	 1
    HISyLATthc 	 0.0
    GLUDym 	 3.0
    r2168 	 -1.0
    RE1530M 	 -1.0
    r2244 	 -1.0
    RDH4 	 -1.0
    RE2444C 	 1
    MI134P4P 	 1
    3ISPVShc 	 1
    FAOXC181C181OHm 	 3.0
    6BHGLZtev 	 1
    r2403 	 -1.0
    GALSGLT1le 	 0.0
    r0840 	 1
    EX_taxol_LPAREN_e_RPAREN_ 	 1
    r2318 	 1.0
    ARACHCPT2 	 0.0
    r0671 	 2.0
    GGH_5DHFe 	 3.0
    PRO1x 	 -1.0
    r2022 	 -1.0
    r1891 	 0.0
    B3GALT42g 	 -1.0
    IDOURte 	 1
    CORE5GTg 	 1
    r1722 	 0.0
    r2159 	 -1.0
    r2396 	 -1.0
    QUILSYN 	 1
    EX_ade_LPAREN_e_RPAREN_ 	 3
    r0451 	 1.0
    r1030 	 1
    EX_gt1a_hs_LPAREN_e_RPAREN_ 	 1
    GGH_5DHFl 	 3.0
    EX_prostge2_LPAREN_e_RPAREN_ 	 1
    SPHK21c 	 1.0
    EX_glc_LPAREN_e_RPAREN_ 	 3
    P4504F123r 	 0.0
    r1987 	 0.0
    FUCACNGALACGLCGALGLUSIDEtg 	 1
    r0345 	 3.0
    EX_glyb_LPAREN_e_RPAREN_ 	 3
    BCDO 	 0.0
    RE2382R 	 0.0
    1a_25VITD2Hm 	 1
    SERALANaEx 	 0.0
    r2087 	 2.0
    B3GALT41g 	 -1.0
    RE2147C 	 -1.0
    FAOXC181_11Em 	 3.0
    FVSitr 	 1
    RN0029R 	 0.0
    NTD4e 	 1.0
    COQ3m 	 -1.0
    r0941 	 1
    RE0938E 	 0.0
    RE3474R 	 0.0
    SELt4_3 	 -1.0
    PHYTt 	 1
    A_MANASEly 	 2.0
    RE2112R 	 -1.0
    DM_Lcystin 	 3
    CDPDAGtm 	 1
    SPHMDAc 	 1
    ADRNLPVESSEC 	 0.0
    EX_o2_LPAREN_e_RPAREN_ 	 3
    RE3154R 	 1
    FUCFUC132GALACGLCGAL14ACGLCGALGLUSIDEte 	 1
    r2092 	 1
    CRGLZhr 	 1
    NTD1m 	 -1.0
    NDPK8n 	 0
    CORE8GTg 	 1
    RE2973N 	 3.0
    NADPHtru 	 1
    AMPtr 	 1
    THRSERNaEx 	 0.0
    ANDRSTRNte 	 1
    r1548 	 3.0
    EX_ascb_L_LPAREN_e_RPAREN_ 	 1
    TMDM1OATt 	 1
    RE3129N 	 1
    RN0013C 	 1
    RE2152C 	 1
    EX_profvs_LPAREN_e_RPAREN_ 	 1
    r1329 	 1
    BILIRED 	 3.0
    URCN 	 -1.0
    RE2873C 	 1
    CLPNDt 	 1
    NDP10ex 	 3.0
    GLUt7l 	 1
    HPETFABP1tc 	 -1.0
    r2107 	 1
    r2127 	 2.0
    r1567 	 3.0
    r0618 	 0
    r1772 	 0.0
    r2519 	 1
    FAOXC184C163m 	 3.0
    TMNDNCCPT1 	 1.0
    SUBEACTD 	 1
    APOCFm 	 1
    SUCD1m 	 2.0
    1513TACRtu 	 1
    r2374 	 0.0
    C12DCTD 	 1
    HAS1 	 3.0
    FAOXC120100x 	 -1.0
    GLYBt4_2_r 	 -1.0
    r2212 	 -1.0
    VLCS2r 	 0.0
    3AIBtm 	 1
    4MPTNLtm 	 1
    C100CPT2 	 0.0
    PVSGLUChc 	 1
    r1934 	 0.0
    r1968 	 0.0
    EX_c4crn_ 	 1
    PTVSTLACtev 	 1
    CRVSM1teb 	 -1.0
    r2131 	 2.0
    XOLEST2HSTDle 	 1
    EX_lcts_LPAREN_e_RPAREN_ 	 3
    RE3522C 	 1
    r2311 	 -1.0
    LCADi 	 3.0
    CYTK5 	 3.0
    ALCD2if 	 3.0
    RE3040R 	 2.0
    31DMTtep 	 1
    RE3124R 	 3.0
    CYSTS 	 0.0
    r2058 	 -1.0
    FAOXC164C143m 	 3.0
    AMPTASECG 	 0.0
    r2004 	 -1.0
    CHOLATEt2 	 -1.0
    RE3014C 	 1
    r1790 	 0.0
    METAT 	 3.0
    RE1830M 	 0.0
    EX_HC00250_LPAREN_e_RPAREN_ 	 3
    L_LACDcm 	 -1.0
    RE1523X 	 1.0
    r1077 	 1
    FACOAL1821 	 3.0
    FAOXC2051843m 	 3.0
    RE2521C 	 1
    1531TACRteb 	 1
    LSTN1GLUCtev 	 1
    r0119 	 3.0
    ST3GAL22g 	 -1.0
    FUCFUCGALACGLCGALGLUSIDEtg 	 1
    H2OGLYAQPt 	 3.0
    2OXOADOXm 	 0
    AKGt4_3 	 -1.0
    NCCt 	 -1.0
    r1913 	 0.0
    r1905 	 0.0
    CYTK9 	 3.0
    FAOXC142C142OHm 	 3.0
    LNLCCRNt 	 1
    FUT15g 	 0.0
    CHOLESACATc 	 1.0
    FUCFUCFUCGALACGLC13GALACGLCGAL14ACGLCGALGLUSIDEtg 	 1
    FADDPle 	 -1.0
    6AHGLZtev 	 1
    AGMTm 	 0.0
    CYTDK1 	 0
    r0191 	 1
    LNLNCGCPT1 	 1.0
    PI45P4P 	 1
    r0517 	 1.0
    MM5ag 	 3.0
    EX_antipyrene_LPAREN_e_RPAREN_ 	 1
    6BHGLZGLCABCt 	 1
    UGT1A3r 	 0
    FVShc 	 1
    THFtm 	 1
    FAOXC102C81x 	 1.0
    AHCYStd 	 1
    P4507A1r 	 -1.0
    CHTNASE 	 0.0
    PPA2m 	 1
    r1900 	 0.0
    r1738 	 0.0
    DOPAMT 	 3.0
    r0559 	 -1.0
    SERt4 	 3.0
    r1674 	 0.0
    TRPHYDRO2 	 0.0
    HPCLx 	 1
    H2Otly 	 1
    PIK3n 	 1
    MM8Ber 	 0.0
    r2180 	 -1.0
    FKYNH 	 -1.0
    DCSPTN1CRNt 	 -1.0
    FAOXC121x 	 1.0
    RE2850C 	 1
    EX_pydx_LPAREN_e_RPAREN_ 	 3
    r0381 	 1
    RE1100L 	 0.0
    OAGT3tg 	 1
    FAOXC22C20x 	 -1.0
    HSD3B11 	 0.0
    SELMETAT 	 3.0
    GCC2cm 	 -1.0
    r2319 	 1.0
    RE2269E 	 -1.0
    RE2848C 	 1
    DGSNt 	 -1.0
    EX_lvst_LPAREN_e_RPAREN_ 	 1
    r1873 	 0.0
    AM1C9CStev 	 1
    r1939 	 0.0
    r2377 	 0.0
    B3GALT43g 	 -1.0
    EX_meoh_LPAREN_e_RPAREN_ 	 3
    OCD11CRNCPT2 	 1
    O2tn 	 1
    24_25DHVITD2t 	 1
    6OHFVSitr 	 1
    HPYRDC 	 1
    FAOXC102C101m 	 3.0
    r0655 	 -1.0
    r2381 	 0.0
    RE3367C 	 -1.0
    SPMDOX 	 1
    PIK4 	 0.0
    PPItm 	 1
    SMPD3g 	 0
    ARGSS 	 3.0
    r2132 	 1
    r1956 	 0.0
    EX_csa_LPAREN_e_RPAREN_ 	 3
    ATPtn 	 1
    DNDPt3m 	 -1.0
    S6T23g 	 0.0
    OCTDECCACT 	 -1.0
    PCREATtmdiffir 	 1
    COAtp 	 1
    RE2659R 	 1
    CHOLESTle 	 3.0
    EX_HC02198_LPAREN_e_RPAREN_ 	 1
    DM_pnto_R 	 1
    RE2525X 	 3.0
    12HTACRhr 	 1.0
    ACACT9p 	 3.0
    6AHGLZABCt 	 1
    EX_rsv_LPAREN_e_RPAREN_ 	 1
    r2222 	 -1.0
    RE0453C 	 1
    PROSTGE1t3 	 1.0
    r2465 	 3.0
    AG13T9g 	 3.0
    ACOAD1fm 	 3.0
    r0027 	 3.0
    RDH2a 	 3.0
    RE2637X 	 0.0
    RE3445X 	 0.0
    NNATr 	 0.0
    r1634 	 3.0
    P5CDm 	 -1.0
    RE0702N 	 -1.0
    PI45PLC 	 1.0
    ALOX12 	 1.0
    AGLPH 	 1
    G14T2g 	 3.0
    RE1525M 	 3.0
    MEPIVESSte 	 1
    IVCOAACBP 	 3.0
    PETHCT 	 -1.0
    GFUCS 	 1.0
    1OHMDZitr 	 1
    EX_pchol_hs_LPAREN_e_RPAREN_ 	 1
    MM6B2g 	 3.0
    ILEATB0tc 	 3.0
    ATVLACtu 	 -1.0
    r1771 	 0.0
    3ISPVSthc 	 0.0
    DURIK1 	 3.0
    RE3010M 	 0
    FAOXC61x 	 1.0
    RE3260C 	 -1.0
    r1853 	 0.0
    RE0575C 	 1
    GLUCYS 	 2.0
    MI145PKn 	 0.0
    3DPHBH1 	 1
    ARACHDtr 	 1
    EX_am1csa_LPAREN_e_RPAREN_ 	 1
    3MOBte 	 1
    RE3341X 	 0.0
    GLYC3Ptmc 	 1
    r1760 	 0.0
    GALFUC12GAL14ACGLCGALGLUSIDEtg 	 1
    DOLPMT2_Uer 	 1
    DESAT24_1 	 0.0
    HKt 	 -1.0
    RE3247X 	 0.0
    ADK1m 	 3.0
    AG13T3g 	 3.0
    EX_utp_LPAREN_e_RPAREN_ 	 3
    6OHFVSGLUhc 	 1
    TMDM3OATt 	 1
    r0395 	 2.0
    ICDHyp 	 3.0
    HCO3_NAt 	 -1.0
    RE3123C 	 1
    r0841 	 1
    S23T2g 	 0.0
    r2292 	 -1.0
    r2312 	 -1.0
    EX_7bhglz_LPAREN_e_RPAREN_ 	 1
    RE3013R 	 0.0
    NACHEX25ly 	 3.0
    RE3446R 	 3.0
    r2367 	 0
    SCPx 	 3.0
    FAOXC164C165x 	 -1.0
    EX_4hatvlac_LPAREN_e_RPAREN_ 	 1
    GLYPROPEPT1tc 	 -1.0
    EX_acmp_LPAREN_e_RPAREN_ 	 1
    r2202 	 -1.0
    MEVK1x 	 3.0
    r0899 	 1
    SUCCt2m 	 -1.0
    FAOXC61_3Zm 	 3.0
    r2063 	 -1.0
    HYXNt 	 -1.0
    CBLATm 	 1.0
    ACMPGLUTitr 	 1
    UREAtm 	 0.0
    FPGS4 	 -1.0
    DNDPt19m 	 -1.0
    RE3012R 	 3.0
    r0750 	 -1.0
    TOLBUTAMIDEte 	 1
    r1693 	 0.0
    C142OHc 	 1.0
    BUTt2r 	 2.0
    RE1803C 	 1
    RE2908C 	 1
    r0074 	 -1.0
    r1990 	 0.0
    EX_malt_LPAREN_e_RPAREN_ 	 3
    r2329 	 1.0
    r1993 	 -1.0
    NACHEXA3ly 	 3.0
    FAH1 	 0.0
    6HLVSTthep 	 1
    r1951 	 0.0
    PTHPSn 	 3.0
    EX_malthx_LPAREN_e_RPAREN_ 	 1
    r1775 	 0.0
    r0744 	 3.0
    RE3111M 	 0.0
    GLCAASE4ly 	 3.0
    r1374 	 0.0
    ENMAN5g 	 1
    r1739 	 0.0
    r2016 	 -1.0
    FAOXC220200x 	 -1.0
    EX_glc3meacp_LPAREN_e_RPAREN_ 	 1
    FAOXC142C122x 	 -1.0
    MMCDp 	 -1.0
    EX_tolbutamide_LPAREN_e_RPAREN_ 	 1
    r1969 	 0.0
    RE0578X 	 0.0
    RNMK 	 2.0
    7BHGLZtev 	 1
    r1796 	 0.0
    r1949 	 0.0
    RE1441G 	 0.0
    SCP21cx 	 3.0
    r1603 	 3.0
    CYTK4 	 3.0
    RE0926E 	 0
    5ADTSTSTERONESULT 	 0
    CYSTSERex 	 -1.0
    SGPL11r 	 3.0
    r0051 	 3.0
    RE3122C 	 1
    r1837 	 0.0
    5FTHFt2 	 -1.0
    24NPHte 	 1
    EX_am1acs_LPAREN_e_RPAREN_ 	 1
    EX_udp_LPAREN_e_RPAREN_ 	 1
    56DHPVStev 	 1
    PHCDm 	 -1.0
    ADHAPtx 	 1
    RE3193M 	 3.0
    UDPGLCtg 	 3.0
    EX_ribflv_LPAREN_e_RPAREN_ 	 3
    EX_tdechola_LPAREN_e_RPAREN_ 	 1
    P4502C9 	 0.0
    r1298 	 1
    r1809 	 0.0
    r2126 	 2.0
    TETPENT6CPT1 	 1.0
    INSt2 	 0.0
    CYSGLUexR 	 0.0
    RE3443M 	 3.0
    RE3165C 	 3.0
    RE2514E 	 -1.0
    S6T10g 	 0.0
    DNDPt49m 	 -1.0
    r0413 	 -1.0
    EX_retinol_cis_11_LPAREN_e_RPAREN_ 	 1
    AP4AH1 	 -1.0
    FUT18g 	 0.0
    r0783 	 3.0
    AM1CSAtep 	 1
    r2201 	 -1.0
    RE2382C 	 1
    AM1ALCShr 	 1
    FAOXC226205m 	 3.0
    RE1653C 	 1
    EX_gtacmp_LPAREN_e_RPAREN_ 	 1
    r0993 	 -1.0
    PVSOATPtu 	 0.0
    EX_4hdebrisoquine_LPAREN_e_RPAREN_ 	 1
    DHGLZtev 	 1
    RE0688E 	 0.0
    PI34P5Kn 	 1
    ITCOALm 	 2.0
    P45021A1r 	 0
    EX_thmmp_LPAREN_e_RPAREN_ 	 1
    RN0027R 	 0.0
    3HBCOARc 	 3.0
    SRTNMTX 	 -1.0
    RE3264C 	 -1.0
    THRASNNaEx 	 0.0
    LVSTACIDhep 	 0.0
    FUT11g 	 0.0
    r2281 	 -1.0
    GCHOLAt2 	 -1.0
    r1937 	 0.0
    RE3041C 	 1.0
    EX_glu_L_LPAREN_e_RPAREN_ 	 3
    r1885 	 0.0
    S23Tg 	 0.0
    G5SDym 	 3.0
    FAOXC141_5Em 	 3.0
    RE3522R 	 0.0
    EX_andrstrnglc_LPAREN_e_RPAREN_ 	 1
    RE1308M 	 1.0
    UGT1A5r 	 0
    TREHe 	 -1.0
    ST3GAL61g 	 0.0
    RE2296X 	 3.0
    GTHPe 	 3.0
    FAOXC161140m 	 3.0
    EX_CE4633_LPAREN_e_RPAREN_ 	 1
    RE3511M 	 -1.0
    RE2910C 	 1
    EX_phe_L_LPAREN_e_RPAREN_ 	 3
    PROtm 	 1
    r2352 	 0
    FA181ACPH 	 -1.0
    r0179 	 1
    r0721 	 1.0
    DESAT18_5 	 3.0
    NH4tp 	 1
    DGSNtm 	 1
    ACMPitr 	 1
    RE2428M 	 -1.0
    RE2704C 	 1
    S6T6g 	 0.0
    FAOXC4020m 	 3.0
    r0839 	 1
    EX_lstnm4_LPAREN_e_RPAREN_ 	 1
    EX_co_LPAREN_e_RPAREN_ 	 1
    RE1064C 	 1
    r0834 	 -1.0
    RE2223M 	 1
    GALASE14ly 	 0
    r0754 	 3.0
    FAOXC16DCr 	 1.0
    EX_whhdca_LPAREN_e_RPAREN_ 	 1
    4HGLSDm 	 -1.0
    EX_6melvacid_LPAREN_e_RPAREN_ 	 1
    EX_triodthysuf_LPAREN_e_RPAREN_ 	 1
    KCC2t 	 3.0
    PROt2rL 	 0.0
    PI4P3Kn 	 1
    RE3446M 	 3.0
    DCK1n 	 2.0
    EX_6ohfvsglu_LPAREN_e_RPAREN_ 	 1
    DOGULNO2 	 1
    HYPTROX 	 1
    EX_HC02207_LPAREN_e_RPAREN_ 	 1
    r0400 	 1
    RE0918C 	 1
    RE1927C 	 1
    RE2997X 	 0.0
    DKMPPD 	 1
    HYPTROXe 	 1
    RE3155C 	 1
    RE3111R 	 0.0
    DNAMTn 	 3.0
    P45039A1r 	 -1.0
    GALt1r 	 3.0
    HISD 	 -1.0
    r1753 	 0.0
    SBPP3er 	 1.0
    EX_vacc_LPAREN_e_RPAREN_ 	 1
    COt 	 1
    FAOXC11090m 	 3.0
    EX_lac_L_LPAREN_e_RPAREN_ 	 3
    GTACMPhr 	 1
    EX_ump_LPAREN_e_RPAREN_ 	 1
    RE0944E 	 0
    S4T6g 	 0.0
    H5MTer_U 	 1
    P4502D6 	 -1.0
    GALASE10ly 	 0
    MANter 	 1
    r2394 	 0.0
    GLYOp 	 1
    EX_spc_hs_LPAREN_e_RPAREN_ 	 1
    EX_1513tacr_LPAREN_e_RPAREN_ 	 1
    DMNONCRNt 	 1
    EX_malcoa_LPAREN_e_RPAREN_ 	 1
    GLB1 	 0
    CSNAT2m 	 -1.0
    RE3345C 	 -1.0
    r0410 	 3.0
    GALASE11ly 	 0
    SIAT9g 	 0.0
    GLNCYSNaEx 	 0.0
    FAOXC165C164x 	 -1.0
    r1043 	 0.0
    MM7B2g 	 3.0
    RE3352C 	 1
    1531TACRitr 	 1
    r0402 	 0.0
    FUCFUCFUCGALACGLC13GALACGLCGAL14ACGLCGALGLUSIDEte 	 1
    CRVSM23tev 	 1
    AM19CShr 	 1
    r1847 	 0.0
    FACOAL244_1 	 3.0
    EX_cdp_LPAREN_e_RPAREN_ 	 1
    P4502C93 	 0.0
    LVSTACIDitr 	 1
    CSCPASEly 	 1
    EX_1331tacr_LPAREN_e_RPAREN_ 	 1
    RE0689X 	 0.0
    RE2863C 	 1
    TYRTA 	 3.0
    MDZGLCitr 	 1
    RE3144C 	 1
    THYMDtl 	 -1.0
    r2154 	 -1.0
    C181CPT1 	 1.0
    r0407 	 3.0
    EX_thrfvs_LPAREN_e_RPAREN_ 	 1
    FAOXC14C14OHm 	 3.0
    PSYGCHe 	 1
    ACNGALACGLCGAL14ACGLCGALGLUSIDEte 	 1
    P4504F81r 	 0.0
    DALAxt 	 1
    PIter 	 -1.0
    EX_pglyc_hs_LPAREN_e_RPAREN_ 	 3
    ESTRONEtr 	 1
    GRTT 	 0.0
    RE3499C 	 1
    2HATVLACteb 	 -1.0
    RE3033R 	 1
    EX_htaxol_LPAREN_e_RPAREN_ 	 1
    CYTK7n 	 3.0
    RE3336X 	 3.0
    r1823 	 0.0
    RE2989X 	 3.0
    RE0583C 	 1
    EX_4mptnl_LPAREN_e_RPAREN_ 	 1
    RE3155R 	 1
    6HSMVACIDteb 	 1
    RE0578C 	 0
    r1959 	 0.0
    RE3475N 	 2.0
    RE3248X 	 1.0
    SERASNNaEx 	 0.0
    3HPVSCOAitx 	 1
    NTD5 	 3.0
    r2473 	 -1.0
    DNDPt59m 	 -1.0
    ACACt2 	 2.0
    DM_13_cis_oretn_n_ 	 1
    CSAtd 	 1
    MLTG1ly 	 0.0
    C161CRNt 	 -1.0
    NACDe 	 1
    FOLOATPtc 	 -1.0
    r1377 	 1
    r0556 	 0
    AM1C9CSteb 	 1
    RE2133C 	 1
    FTHFL 	 3.0
    EX_c5dc_ 	 1
    TCHOLAt3 	 0.0
    FAOXC181C161x 	 -1.0
    C80CPT1 	 1.0
    OCTAt 	 1
    PRGNLONESULT 	 2.0
    DMNONCOACRNCPT1 	 1
    BUTt2m 	 2.0
    LGNCCOAtcx 	 -1.0
    r1978 	 0.0
    r2203 	 -1.0
    CYTK1m 	 3.0
    FAOXC2242046x 	 -1.0
    ACtg 	 1
    DOPACHRMISO 	 -1.0
    3MOBt2im 	 1
    RE3525M 	 1
    SERPT 	 0.0
    r2028 	 -1.0
    RE3494C 	 1
    RE3268R 	 1
    3AIBTm 	 1
    r0706 	 0.0
    r1789 	 0.0
    EX_6bhglzglc_LPAREN_e_RPAREN_ 	 1
    r1011 	 1
    FACOAL100i 	 1
    FUC13GALACGLCGAL14ACGLCGALGLUSIDEte 	 1
    C4OHc 	 1.0
    RE3038C 	 0.0
    NMNATm 	 0.0
    EX_dgmp_LPAREN_e_RPAREN_ 	 1
    RE1632R 	 1.0
    FAOXC123x 	 1.0
    RE3561X 	 0.0
    C100CPT1 	 1.0
    r2495 	 -1.0
    RE1530C 	 1
    TDCHOLAtx 	 1
    NTD2e 	 1.0
    RE0827X 	 0.0
    RE1050L 	 3.0
    RE0927C 	 -1.0
    BTNPL 	 -1.0
    r2266 	 -1.0
    r1172 	 3.0
    LEULEULAPc 	 2.0
    r2290 	 -1.0
    INSTt2r 	 0.0
    HACD1x 	 3.0
    EX_acald_LPAREN_e_RPAREN_ 	 3
    EX_gluala_LPAREN_e_RPAREN_ 	 3
    RE2856C 	 1
    AGLPED 	 1
    PI345P5Pn 	 1
    r2235 	 -1.0
    TETTET6CRNt 	 -1.0
    DNDPt29m 	 1
    r1174 	 3.0
    FATP2t 	 -1.0
    RE2655R 	 1
    FAOXTC162TC142m 	 3.0
    EX_3bcrn_ 	 1
    PROFVStev 	 1
    RE3038N 	 1
    DESAT20_2 	 0
    XYLK 	 -1.0
    AOBUTDsm 	 1
    FAOXC5C5OHm 	 1
    C3STDH1Pr 	 2.0
    PYLALDOXm 	 3.0
    CSASULPteb 	 1
    PECGONCOATr 	 1
    r1466 	 1
    EX_CLPND_LPAREN_e_RPAREN_ 	 3
    FAOXC101C102m 	 3.0
    DOCOSACTDe 	 1
    FAOXC8DCC6DCx 	 -1.0
    ACACT1x 	 0.0
    RE0456M 	 1
    FAOXC15NADx 	 0.0
    RE3113R 	 1
    3OHACMPtev 	 1
    C2tcx 	 -1.0
    EX_ha_pre1_LPAREN_e_RPAREN_ 	 1
    r0525 	 -1.0
    r0774 	 3.0
    r0160 	 -1.0
    FPGS7 	 -1.0
    FUT910g 	 -1.0
    HSD3B7P 	 0.0
    r0010 	 3.0
    RE1808R 	 1
    MCD 	 -1.0
    r0769 	 0.0
    Ser_Thrtg 	 1
    RE1635M 	 3.0
    r2494 	 -1.0
    r2053 	 -1.0
    EX_am9csa_LPAREN_e_RPAREN_ 	 3
    RE1898C 	 1
    RE0576C 	 1
    EX_srtn_LPAREN_e_RPAREN_ 	 1
    RE0920R 	 0
    NDPK5n 	 0
    HEDCECRNe 	 1
    RE2622C 	 -1.0
    2HATVLACGLUChr 	 0
    URIK1 	 0
    r0643 	 3.0
    RE1816X 	 3.0
    ADPRDP 	 3.0
    r1919 	 0.0
    EX_tyr_L_LPAREN_e_RPAREN_ 	 3
    r0718 	 3.0
    3HPVSCOAhc 	 1
    r2507 	 1
    EX_4hmdgluc_LPAREN_e_RPAREN_ 	 1
    AKGMALtm 	 -1.0
    CYSATB0tc 	 3.0
    FATP7t 	 0.0
    LSTN1GLUCitr 	 1
    EX_rsvlac_LPAREN_e_RPAREN_ 	 1
    4HATVACIDitr 	 1
    r2097 	 2.0
    CYSPHELAT2tc 	 1
    GLNtN1 	 -1.0
    RE0577M 	 3.0
    r1638 	 3.0
    MEOHtly 	 1
    SERATB0tc 	 3.0
    CRVSM22itr 	 1
    PHETA1 	 3.0
    EX_maltpt_LPAREN_e_RPAREN_ 	 1
    FAOXC2251836x 	 -1.0
    r0795 	 1
    GLCNACDASg 	 -1.0
    THCHOLSTOICtm 	 1
    r1985 	 0.0
    DNDPt42m 	 -1.0
    r0548 	 3.0
    r0845 	 3.0
    r0309 	 -1.0
    r2416 	 -1.0
    r0295 	 -1.0
    CITL 	 1
    RETI3 	 1
    r1814 	 0.0
    FUT35g 	 1.0
    ALAGLNexR 	 -1.0
    KSII_CORE2tly 	 1
    FACOAL184 	 3.0
    EX_5fthf_LPAREN_e_RPAREN_ 	 1
    UDPDOLPT_U 	 1
    RTOTALFATPc 	 1
    2HATVACIDteb 	 -1.0
    NDPK2n 	 0
    ADSK 	 3.0
    DASCBR 	 3.0
    r0358 	 3.0
    CERT1rt 	 -1.0
    r0885 	 -1.0
    THMDt2r 	 0.0
    RE2442C 	 1
    S6T16g 	 0.0
    FAOXC201C181x 	 -1.0
    GLYPHEPEPT1tc 	 -1.0
    HSD11B1r 	 0.0
    RE3040X 	 0.0
    CYTK1 	 3.0
    NADHtpu 	 1
    LSTN1GLUChr 	 0
    RE3580X 	 0.0
    DGULND 	 1
    RE0702L 	 3.0
    ACNAMPH 	 1.0
    TMDK1m 	 -1.0
    DOPAt4_2_r 	 -1.0
    LCAT1e 	 0.0
    RE3097X 	 3.0
    ATVACIDtu 	 -1.0
    ADPCOACROT 	 1.0
    RE1956C 	 0.0
    7BHGLZGLCABCt 	 1
    EX_leuleu_LPAREN_e_RPAREN_ 	 1
    11DOCRTSTRNtm 	 1
    GLYCLTDym 	 1
    S6TASE14ly 	 3.0
    LCYSTCBOXL 	 1.0
    FAOXC162_7Z_10Zm 	 3.0
    VLCSr 	 0.0
    ARGATB0tc 	 3.0
    FUMTSULtm 	 -1.0
    r1062 	 1
    CYTK10 	 3.0
    CRVSM31itr 	 1
    r0113 	 0.0
    EX_fuc_L_LPAREN_e_RPAREN_ 	 3
    PROPAT4te 	 -1.0
    FAOXC241221x 	 -1.0
    r0083 	 3.0
    DNDPt57m 	 -1.0
    DM_avite1_c_ 	 1
    DM_anth 	 1
    MI1345PP 	 1.0
    r0033 	 -1.0
    r0661 	 1.0
    DOGULND2 	 1
    GCHOLAte 	 1
    EX_caro_LPAREN_e_RPAREN_ 	 3
    RE2027C 	 1
    ACMPGLUtep 	 0.0
    DESAT18_3 	 3.0
    FUM 	 3.0
    ALR3 	 3.0
    FAOXC101_3Em 	 3.0
    TRIPVShc 	 1.0
    LSTNM5itr 	 1
    RE2383C 	 1
    PCHOLP_hs 	 2.0
    EX_c81crn_ 	 1
    AG13T10g 	 3.0
    CO2tn 	 1
    RE1952X 	 0.0
    EX_fru_LPAREN_e_RPAREN_ 	 3
    LEUKABCtc 	 0.0
    EX_creat_LPAREN_e_RPAREN_ 	 3
    RE2273C 	 1
    RE0958C 	 1
    r2147 	 -1.0
    EX_adprbp_LPAREN_e_RPAREN_ 	 1
    FAOXC15NADPx 	 0.0
    r2083 	 2.0
    r0758 	 3.0
    UGT1A8r 	 0
    PPAtm 	 1
    ASNNm 	 0.0
    r1596 	 3.0
    RE0915C 	 1
    AHCYSte 	 1
    FBA5 	 -1.0
    r1874 	 0.0
    MTHFCm 	 3.0
    EX_vitd3_LPAREN_e_RPAREN_ 	 1
    DESAT18_9 	 0
    PPAm 	 3.0
    RDH3 	 -1.0
    C226CPT1 	 1.0
    TETPENT3t 	 1
    r0727 	 1.0
    FA1821ACPH 	 1
    FAOXC205_5Z_8Z_11Z_14Z_17Zx 	 -1.0
    CRVS1M24hc 	 1.0
    TETHEX3COAtx 	 1
    r0377 	 2.0
    r1391 	 3.0
    RE3514C 	 1
    34DHPHELAT1tc 	 1
    PMI12346PH 	 0.0
    1a_25VITD3Hm 	 1
    r0617 	 0
    EX_ptvstlac_LPAREN_e_RPAREN_ 	 3
    PCLYSOX 	 3.0
    3HPPD 	 1
    RE0580R 	 1
    THMDt5le 	 0.0
    GLCURtly 	 -1.0
    RE1796R 	 0.0
    AFLATOXINte 	 1
    GGH_7THFe 	 3.0
    RE0828E 	 0.0
    EX_carn_LPAREN_e_RPAREN_ 	 3
    r1698 	 0.0
    r0723 	 3.0
    RE3310C 	 1
    DNDPt22m 	 1
    r0789 	 1
    r1888 	 0.0
    4MOPt2im 	 1
    r2262 	 -1.0
    IBUPGLUCitr 	 1
    RE0453M 	 0.0
    F1PGT 	 1.0
    C204CPT2 	 0.0
    MM6bg 	 3.0
    APRTO2 	 3.0
    r2181 	 -1.0
    r2489 	 -1.0
    PLA2_2e 	 0.0
    AMACR2r 	 -1.0
    r1658 	 3.0
    EX_4mtolbutamide_LPAREN_e_RPAREN_ 	 1
    RE2913M 	 3.0
    r0392 	 3.0
    EX_ala_B_LPAREN_e_RPAREN_ 	 3
    TCHOLAtx 	 1
    r1073 	 1
    r0311 	 3.0
    r2308 	 -1.0
    UGT1A7r 	 0
    C3STKR2r 	 3.0
    CREATt4_2_r 	 3.0
    DM_Asn_X_Ser_Thr_ly_ 	 1
    CYTK8n 	 3.0
    B3GNT37g 	 0.0
    RE3103C 	 3.0
    PMEVKc 	 0.0
    ST6GALNAC31 	 0.0
    HGNTOR 	 0.0
    EX_CE2839_LPAREN_e_RPAREN_ 	 1
    RE3521R 	 3.0
    RE3572X 	 3.0
    DGCHOLte 	 1
    CARIBUP_SGLUhep 	 0
    LINKDEG4ly 	 1
    EX_abt_LPAREN_e_RPAREN_ 	 1
    EX_1mncam_LPAREN_e_RPAREN_ 	 1
    r1177 	 3.0
    GTHO 	 2.0
    RETNGLCt2r 	 1
    r0568 	 3.0
    RE2913X 	 0.0
    56EPPVStev 	 1
    EX_chsterol_LPAREN_e_RPAREN_ 	 3
    RE3458C 	 1
    GLYPROPRO1c 	 -1.0
    RE2915M 	 3.0
    r1998 	 -1.0
    r1705 	 0.0
    PCHOLABCtc 	 3.0
    RE3247M 	 -1.0
    r1444 	 0.0
    FUCASE2e 	 2.0
    RE0581R 	 0
    RE3167R 	 1
    AG13T5g 	 3.0
    FPGS2 	 -1.0
    r1579 	 3.0
    RE2666C 	 1
    DADNK 	 1
    r1690 	 0.0
    r0465 	 -1.0
    RE3519X 	 0.0
    FVSGLUChc 	 1
    HISTAVESSEC 	 0.0
    ARTFR207 	 1
    GUAt 	 -1.0
    r2170 	 -1.0
    r0697 	 0.0
    EX_anth_LPAREN_e_RPAREN_ 	 1
    LDH_L 	 3.0
    FAOXC225C204x 	 1.0
    1OHMDZhr 	 1.0
    3HPVSTETCOAitx 	 1
    ESTROSABCCte 	 -1.0
    COQ6m 	 -1.0
    RN0030C 	 1
    EX_fucfuc12gal14acglcgalgluside_hs_LPAREN_e_RPAREN_ 	 1
    ACCOALm 	 0.0
    r0423 	 3.0
    EX_sel_LPAREN_e_RPAREN_ 	 3
    C4tcx 	 -1.0
    FAOXC2031836m 	 3.0
    r1756 	 0.0
    LTD4DP 	 1
    DOCOSADIACTD 	 1
    S6T20g 	 0.0
    RE3259R 	 2.0
    r1972 	 0.0
    C40CPT1 	 1.0
    ACACT7p 	 3.0
    RSVitr 	 1
    RE3476X 	 0
    CERK 	 0.0
    DTTPtm 	 1
    RE3524R 	 1
    RADH3 	 1
    10FTHFtl 	 1
    CYSTA 	 3.0
    B3GNT311g 	 0.0
    ST6GALNAC23 	 3.0
    2AMADPTm 	 -1.0
    r2512 	 1
    SUBERCACT 	 -1.0
    GALGT1 	 0.0
    EX_3octdeccrn_ 	 1
    RE3345X 	 3.0
    ACCOACm 	 -1.0
    P4501B1r 	 1
    PGI 	 3.0
    RE1709N 	 1.0
    GGH_6DHFl 	 3.0
    EX_lstnm5_LPAREN_e_RPAREN_ 	 1
    RE3126C 	 1
    r1859 	 0.0
    RE3171C 	 1
    RE1518M 	 0.0
    O2St 	 1
    RE3040C 	 0.0
    2HATVLACGLUCteb 	 -1.0
    EX_tag_hs_LPAREN_e_RPAREN_ 	 3
    EBASTINEte 	 1
    TRYPTAOX 	 3.0
    DNDPt32m 	 -1.0
    DHFR 	 1.0
    CITt4_2 	 -1.0
    2HCO3_NAt 	 -1.0
    EX_acmpglut_LPAREN_e_RPAREN_ 	 1
    EX_7ahglz_LPAREN_e_RPAREN_ 	 1
    PHEACGLNt 	 1
    PPM 	 3.0
    FRUt4 	 -1.0
    AM4NCShr 	 1
    3HLVSTACtbc 	 1
    CLHCO3tex2 	 -1.0
    FTHFCL 	 0
    COQ7m 	 -1.0
    EX_6csmvacid_LPAREN_e_RPAREN_ 	 1
    24_25DHVITD3t 	 1
    RE2079R 	 1
    EX_ahandrostanglc_LPAREN_e_RPAREN_ 	 1
    RE2986X 	 1.0
    EX_dcyt_LPAREN_e_RPAREN_ 	 3
    EX_3hexdcrn_ 	 1
    BHBtm 	 1
    RE3551X 	 1
    r1598 	 3.0
    ALAt2r 	 0.0
    r0992 	 -1.0
    O16G2e 	 0.0
    S6TASE22ly 	 -1.0
    RE2112C 	 1
    ISOLVSTtbc 	 -1.0
    RE2878C 	 1
    EX_nrpphrsf_LPAREN_e_RPAREN_ 	 1
    LSTNM1itr 	 1
    RE2854C 	 1
    TRPATB0tc 	 3.0
    EX_35cgmp_LPAREN_e_RPAREN_ 	 1
    FUT12g 	 0.0
    HOCDACBP 	 3.0
    FAOXC121C10x 	 1.0
    FVSteb 	 1
    EX_ocdca_LPAREN_e_RPAREN_ 	 3
    NDPK4m 	 0.0
    GABAVESSEC 	 -1.0
    r2161 	 -1.0
    RE1701C 	 1
    r1928 	 0.0
    25VITD2Hm 	 1.0
    RE3227C 	 1
    r2356 	 0
    r0470 	 1.0
    MEOHt2 	 1
    RE1817X 	 3.0
    RE0689E 	 0.0
    EX_orot_LPAREN_e_RPAREN_ 	 1
    r1321 	 1
    VALTAm 	 0.0
    LYSMTF2n 	 3.0
    FUCFUCFUCGALACGLCGAL14ACGLCGALGLUSIDEtg 	 1
    EX_retinol_LPAREN_e_RPAREN_ 	 3
    TRIODTHYt2 	 0.0
    EX_yvite_LPAREN_e_RPAREN_ 	 1
    RE3148C 	 3.0
    ENMAN1g 	 1
    NAGAlby 	 1
    RE3186M 	 3.0
    MMMm 	 3.0
    TRPt4 	 3.0
    FAOXC143C123x 	 -1.0
    r1010 	 1
    EX_gumtchol_LPAREN_e_RPAREN_ 	 1
    EX_gumgchol_LPAREN_e_RPAREN_ 	 1
    6CSMVitr 	 1
    EX_apnnox_LPAREN_e_RPAREN_ 	 1
    EX_retnglc_LPAREN_e_RPAREN_ 	 1
    r1679 	 0.0
    RE3565C 	 1
    r1692 	 0.0
    SMVtv 	 1
    r1315 	 0.0
    r0615 	 -1.0
    S6TASE3ly 	 3.0
    ACMPtu 	 -1.0
    CHOLPtg 	 1
    N4Tg 	 1
    LPS4e 	 0.0
    NNMT 	 3.0
    TRIPVStev 	 1
    2DR1PP 	 1
    FAOXC185_3Z_6Z_9Z_12Z_15Zm 	 3.0
    S6T2g 	 0.0
    DOLASNT_Uer 	 1
    r2114 	 1
    RE3126R 	 1
    GALNTg 	 3.0
    r0946 	 -1.0
    r1810 	 0.0
    FAOXC16DCC14DCx 	 -1.0
    r1166 	 1.0
    r2321 	 1.0
    r2274 	 -1.0
    r1135 	 2.0
    RE3581X 	 -1.0
    DNDPt47m 	 -1.0
    EX_dgsn_LPAREN_e_RPAREN_ 	 3
    RE1099R 	 0.0
    r2307 	 -1.0
    VALPHELAT2tc 	 1
    P4502C19 	 0.0
    TCHOLAte 	 1
    RE3554C 	 1
    S6T14g 	 0.0
    GLCAT9g 	 2.0
    EX_2hatvacid_LPAREN_e_RPAREN_ 	 1
    ACOAD8m 	 -1.0
    r2204 	 -1.0
    RE2272L 	 3.0
    1a_24_25VITD2Hm 	 1
    GALACGLCGALGBSIDEte 	 1
    RE2124C 	 1
    RE3335M 	 3.0
    FORtrn 	 1
    ATVLACtdhc 	 1
    NDPK2m 	 0.0
    RE1534X 	 3.0
    ORNt3m 	 -1.0
    RE2766C 	 -1.0
    GBAl 	 0
    BHMT 	 -1.0
    RE3519R 	 2.0
    BUTSMCT1 	 -1.0
    r2249 	 -1.0
    1HMDZGLUChc 	 1
    r0249 	 1
    NDPK2 	 0.0
    FUT31g 	 1.0
    r1008 	 1
    ACGBGBSIDEtg 	 1
    GGT_U 	 1
    EX_CE0074_LPAREN_e_RPAREN_ 	 1
    MAN2_7Cer 	 1
    HOCTDECCRNe 	 1
    SBTR 	 3.0
    r0995 	 -1.0
    DOPASFt 	 1
    RE1308C 	 1
    TYRPHELAT2tc 	 1
    C3STDH1r 	 2.0
    RE2899C 	 1
    RE2270E 	 0.0
    UAGALDP 	 1
    FUCASEly 	 3.0
    MAOX 	 0.0
    EX_CE2011_LPAREN_e_RPAREN_ 	 1
    7DHFtl 	 1
    ILELAT1tc 	 1
    RE3265R 	 0
    LST4EXPhr 	 1.0
    BILDGLCURte 	 3.0
    LGNCFATtc 	 3.0
    RE2677C 	 1
    C161OHc 	 1.0
    r2096 	 2.0
    ALCD22_L 	 3.0
    DOLPGT2_Ler 	 3.0
    7HPVStev 	 1
    ECOAH1m 	 3.0
    r1327 	 1
    MM7Ag 	 0.0
    S6T17g 	 0.0
    TRPB0AT3tc 	 -1.0
    FAOXC61C4m 	 1.0
    EX_3ivcrn_ 	 1
    EX_k_LPAREN_e_RPAREN_ 	 3
    ACNACNGAL14ACGLCGALGLUSIDEte 	 1
    RE1100C 	 0.0
    RE3451C 	 1
    FAOXC15BRC13BRx 	 0.0
    LYSATB0tc 	 3.0
    RE3295C 	 1
    r0398 	 2.0
    r0773 	 0.0
    DM_gncore2_g_ 	 1
    PETOHMm_hs 	 -1.0
    C10DCe 	 1
    FAOXC143_5Z_8Z_11Zx 	 -1.0
    r1803 	 0.0
    RE3378C 	 1
    ESTRADIOLGLCt 	 1.0
    GAO1 	 3.0
    r2185 	 -1.0
    VITKtl 	 1
    PPASMCT1 	 -1.0
    r2354 	 0
    35DHPVShc 	 1.0
    GLYC3PFADm 	 1
    MHGLZtev 	 1
    r2487 	 -1.0
    EX_but_LPAREN_e_RPAREN_ 	 1
    r2272 	 -1.0
    EX_pydam_LPAREN_e_RPAREN_ 	 3
    EX_hdcea_LPAREN_e_RPAREN_ 	 3
    EX_arab_L_LPAREN_e_RPAREN_ 	 3
    r1759 	 0.0
    PGPPT 	 -1.0
    r1954 	 0.0
    r0713 	 0.0
    SO4tl 	 1
    HYPOE 	 0
    EX_ddece1crn_ 	 1
    PLA2 	 1
    RE1311C 	 1
    PTVSThc 	 0.0
    FAOXC101_4Em 	 3.0
    MI1346Ptn 	 1
    PYDXPP 	 0
    r1106 	 1
    PE_HSter 	 1
    FUT97g 	 -1.0
    EX_fucfucfucgalacglcgal14acglcgalgluside_hs_LPAREN_e_RPAREN_ 	 1
    GGH_10FTHF5GLUe 	 3.0
    P45017A1r 	 -1.0
    HSPASEly 	 1
    EX_C02528_LPAREN_e_RPAREN_ 	 1
    FAOXC8C6x 	 -1.0
    RE2050R 	 2.0
    FVSTETtev 	 1
    UROLACer 	 1
    7AHGLZhr 	 1
    SPHGNtr 	 1
    EX_glyc_LPAREN_e_RPAREN_ 	 3
    r1882 	 0.0
    RE0921C 	 -1.0
    OBDHc 	 1
    INStm 	 -1.0
    r2405 	 -1.0
    r0281 	 0.0
    r0494 	 2.0
    34DHPLACOX 	 3.0
    PNTKm 	 3.0
    r1988 	 0.0
    PRODt2rL 	 0.0
    RE1943C 	 1
    RE0830N 	 1.0
    HOCTDEC2CRNe 	 1
    ECOAH1x 	 3.0
    SRTN23OX 	 2.0
    r1747 	 0.0
    CYSLYSL 	 1
    r2360 	 0
    DOGULND1 	 1
    TRPO2 	 2.0
    r0821 	 1
    HPDCACRNCPT2 	 0.0
    r2257 	 -1.0
    r2282 	 -1.0
    NDPK1 	 0.0
    PHACCOAGLNAC 	 1
    RE3259C 	 0.0
    RE3269C 	 1.0
    NP1 	 3.0
    EX_am1c9cs_LPAREN_e_RPAREN_ 	 1
    NACHEX2ly 	 3.0
    DNDPt38m 	 -1.0
    RE3258C 	 0.0
    ATP2ter 	 1
    r0321 	 1
    ACCOAtn 	 1
    VITEtl 	 3.0
    MI13PP 	 1
    PVStep 	 1
    GHMT3m 	 3.0
    r0639 	 3.0
    PI3P4Kn 	 1
    G13MT_U 	 1
    PIK3er 	 3.0
    RE3448X 	 3.0
    ABUTt2r 	 0.0
    CRNCAR3tp 	 1
    EX_cytd_LPAREN_e_RPAREN_ 	 3
    C81CRNe 	 1
    r0770 	 0.0
    7AHGLZtev 	 1
    OCDEAFABP1tc 	 -1.0
    r1780 	 0.0
    r0226 	 1
    RE3050R 	 1
    RE2974R 	 1.0
    MDZGLCtev 	 1
    DNDPt63m 	 -1.0
    NDPK7n 	 0
    6AHGLZitr 	 1
    FAOXC123C102x 	 1.0
    RE0702E 	 -1.0
    r1472 	 0.0
    LEUt5m 	 1
    r1943 	 0.0
    UDPG4E 	 2.0
    MG3er 	 3.0
    FAOXC2251836m 	 3.0
    r1423 	 1
    EX_tlacfvs_LPAREN_e_RPAREN_ 	 1
    GALNACT5g 	 1.0
    FE3MTP1 	 3.0
    EX_leuktrA4_LPAREN_e_RPAREN_ 	 3
    RE0453N 	 3.0
    FCLTm 	 0.0
    NACHEXA18ly 	 3.0
    r1597 	 3.0
    C162OHe 	 1
    FAOXC102m 	 3.0
    5HXKYNDCL 	 0.0
    4BHGLZhr 	 1
    r1675 	 0.0
    LYSMTF1n 	 3.0
    TYRTAm 	 3.0
    XOL7AH2tm 	 1
    r1941 	 0.0
    CSAPASEly 	 1
    FAOXC163_4Z_7Z_10Zx 	 -1.0
    DESFVShc 	 0.0
    H7_ETer 	 0.0
    THBPT4ACAMDASE 	 1.0
    EX_retn_LPAREN_e_RPAREN_ 	 3
    r0668 	 0.0
    RE0066M 	 -1.0
    r2042 	 -1.0
    RE3566C 	 1
    BALABETAtc 	 3.0
    r1898 	 0.0
    ASCBt 	 1
    r0636 	 1
    RE3440C 	 1
    EX_smvacid_LPAREN_e_RPAREN_ 	 1
    r1632 	 0.0
    UGCG 	 3.0
    TLACFVSitr 	 1
    FRTT 	 0.0
    3HPVSTEThc 	 1
    r0931 	 1
    RE3139X 	 3.0
    r1255 	 3.0
    FAOXC123m 	 3.0
    SPMTDe 	 1
    RPE 	 1.0
    r1848 	 0.0
    r0308 	 3.0
    SPH1Pte 	 1
    F1Atg 	 1
    EX_leuktrF4_LPAREN_e_RPAREN_ 	 1
    RE0926C 	 1
    HXANtx 	 1
    RE3520E 	 1
    r1680 	 0.0
    EX_tsacmsul_LPAREN_e_RPAREN_ 	 1
    TYMSULT 	 2.0
    TRIOK 	 1
    DNDPt18m 	 -1.0
    GLUTCOAACBP 	 3.0
    HSD3B12r 	 -1.0
    C4STMO2r 	 3.0
    r2537 	 1
    ENGASE2ly 	 1
    r1580 	 3.0
    PHEyLATthc 	 0.0
    RE3513C 	 1
    FATP5t 	 -1.0
    DNDPt56m 	 -1.0
    EX_CE2250_LPAREN_e_RPAREN_ 	 1
    FAOXC246226x 	 -1.0
    RE1537X 	 -1.0
    FACOAL181i 	 3.0
    TSTSTERONESULT 	 0
    PROFVSCOAitx 	 1
    r0645 	 1
    r0511 	 3.0
    GLYCK2 	 1
    EX_tmdm3_LPAREN_e_RPAREN_ 	 1
    NDPK7 	 0.0
    PI45PLCn 	 1.0
    EX_crtstrn_LPAREN_e_RPAREN_ 	 1
    56DHPVSteb 	 -1.0
    RE0923R 	 0
    UMPK2 	 3.0
    EX_3hsmvacid_LPAREN_e_RPAREN_ 	 1
    r0673 	 0.0
    r1942 	 0.0
    r0647 	 1
    PNTK 	 3.0
    EX_caribup_R_LPAREN_e_RPAREN_ 	 1
    SPODM 	 3.0
    RE2872C 	 1
    MTHFD2 	 3.0
    RE1830C 	 1
    r2037 	 -1.0
    RE3233C 	 1
    EX_4abut_LPAREN_e_RPAREN_ 	 1
    Uritn 	 1
    THMt2m 	 1
    SFCYSe 	 1
    AKR1D2 	 -1.0
    EX_arach_LPAREN_e_RPAREN_ 	 1
    RE3511C 	 1
    EX_akg_LPAREN_e_RPAREN_ 	 3
    LIPOti 	 0.0
    B3GALT5g 	 0.0
    SPHINGStr 	 1
    DRIBt 	 1
    ALCD1 	 3.0
    r2496 	 -1.0
    ETHP 	 -1.0
    FAOXTC182TC162m 	 3.0
    5OHFVSGLUtev 	 1
    5MTHFt2le 	 1
    r1007 	 1
    r1910 	 0.0
    O2t 	 1
    RE3402M 	 3.0
    TRDR2 	 3.0
    6MELVACIDitr 	 1
    EX_HC02179_LPAREN_e_RPAREN_ 	 1
    r2211 	 -1.0
    ADPACTD 	 1
    ACHVESSEC 	 0.0
    EX_cysam_LPAREN_e_RPAREN_ 	 1
    FAOXC201181x 	 -1.0
    GALNACT4g 	 1.0
    ICDHyrm 	 3.0
    UMPK5n 	 3.0
    ARTPLM3m 	 1
    r0927 	 1
    EX_HC02180_LPAREN_e_RPAREN_ 	 1
    RN0021X 	 0.0
    r0205 	 3.0
    FORt2m 	 1
    r2334 	 1.0
    UGALGTg 	 -1.0
    PEPCK 	 -1.0
    CRNtuNa 	 0.0
    FAOXC22OHC22r 	 0.0
    RN0028C 	 0.0
    RE1817R 	 3.0
    r0441 	 1
    GNMT 	 -1.0
    GLUTCOADHm 	 -1.0
    NABTNO 	 3.0
    EX_am1ccs_LPAREN_e_RPAREN_ 	 1
    FACOAL140i 	 2.0
    SMVHYDROhep 	 0.0
    r0391 	 2.0
    RE3488X 	 2.0
    RE1240C 	 1
    ALLOP2tu 	 0.0
    r0170 	 3.0
    ACTNMO 	 0.0
    FACOAL1831 	 3.0
    THSACMPhr 	 1
    ABO2g 	 -1.0
    RE3073X 	 1
    RE3162C 	 1
    RE1135L 	 0.0
    r1059 	 -1.0
    r0339 	 -1.0
    RE1944C 	 1
    ST8SIA11 	 -1.0
    RE3401M 	 3.0
    RE1517X 	 1
    EX_4pyrdx_LPAREN_e_RPAREN_ 	 1
    CORE7GTg 	 1
    RE3437C 	 1
    EX_gum_LPAREN_e_RPAREN_ 	 1
    r2055 	 -1.0
    r1594 	 3.0
    GLUt2m 	 2.0
    P45017A3r 	 -1.0
    r0483 	 1.0
    PAFH 	 1.0
    HTAXOLte 	 1
    RE3396M 	 3.0
    RE3417C 	 1
    FUT95g 	 -1.0
    EICOSTETCPT1 	 1.0
    DDPGAm 	 1
    CRVSM31hc 	 1
    FAOXC9070m 	 3.0
    CTPS2 	 3.0
    XYLTt 	 1
    r1833 	 0.0
    CARIBUP_Rthv 	 1
    RE2919X 	 1.0
    P45027A15m 	 0.0
    r1659 	 3.0
    3OHACMPitr 	 1
    DNDPt1m 	 -1.0
    PHETA1m 	 3.0
    r1960 	 0.0
    SRTNt6_2_r 	 -1.0
    TXA2tr 	 1
    HSD3A2r 	 -1.0
    ABO9g 	 -1.0
    XYLTer 	 1.0
    r1635 	 3.0
    RE3088X 	 3.0
    LACLt 	 1
    r1733 	 0.0
    EX_6hsmvacid_LPAREN_e_RPAREN_ 	 1
    r1326 	 1
    24_25DHVITD2tm 	 1
    THRFVStev 	 1
    NADHtru 	 1
    THYMDtm 	 -1.0
    RE3170C 	 3.0
    EX_ahcys_LPAREN_e_RPAREN_ 	 1
    NADPHtxu 	 1
    EX_galgalfucfucgalacglcgalacglcgal14acglcgalgluside_hs_LPAREN_e_RPAREN_ 	 1
    XANtx 	 1
    EX_aicar_LPAREN_e_RPAREN_ 	 1
    EBP2r 	 3.0
    WHTSTSTERONEte 	 1
    EX_adprib_LPAREN_e_RPAREN_ 	 1
    RE3436R 	 2.0
    FAOXC226C225x 	 -1.0
    FAOXC101_4Zx 	 -1.0
    SEBCOACROT 	 1.0
    24_25DHVITD3tm 	 1
    RE3002X 	 1.0
    GLCAASE8ly 	 3.0
    FAOXC164C143x 	 1.0
    CSNATr 	 -1.0
    ACt2r 	 1
    APOC_LYS_BTNPm 	 1
    EX_tsul_LPAREN_e_RPAREN_ 	 1
    FAOXC141C121m 	 3.0
    AM1C9CShr 	 1
    HS2ly 	 1.0
    FAOXC13BRC11BRx 	 0.0
    OXYPtepv 	 1
    RE3562M 	 3.0
    UMPK7n 	 3.0
    METATB0tc 	 3.0
    r2038 	 -1.0
    C226COAtx 	 1
    DHPR2 	 1
    STRDNCCPT1 	 1.0
    G14T3g 	 3.0
    r0330 	 1
    TAUBETAtc 	 3.0
    FAS180COA 	 0.0
    ACOX2x 	 -1.0
    GGH_10FTHF6GLUe 	 3.0
    S4TASE5ly 	 -1.0
    56EPPVSitr 	 1
    RE2269C 	 1
    RE2563C 	 1
    RE3415C 	 1
    RE3173R 	 1
    RE0925R 	 0
    GULNDer 	 1
    FAOXC7C5m 	 1.0
    RSVLACteb 	 -1.0
    EX_ca2_LPAREN_e_RPAREN_ 	 3
    P4502E1 	 0.0
    HDCAFAPMtc 	 3.0
    r0730 	 3.0
    CAMPt 	 1.0
    DNDPt2m 	 -1.0
    r1938 	 0.0
    S23T4g 	 0.0
    RE1514M 	 3.0
    LEUGLYHYc 	 1
    VALLAT1tc 	 1
    PGPP_hs 	 -1.0
    SERLYSNaex 	 -1.0
    r2118 	 2.0
    WHDDCAte 	 1
    H2Ot 	 1.0
    RE3533M 	 3.0
    FVSTETGLUtev 	 1
    r1812 	 0.0
    GALASE5ly 	 0
    r2433 	 -1.0
    DSPVSteb 	 -1.0
    S6TASE26ly 	 3.0
    r1616 	 0.0
    r2302 	 -1.0
    HOCTDACBP 	 3.0
    H2O2tly 	 1
    HMGCOASi 	 3.0
    r2080 	 2.0
    UNK2 	 1
    FT 	 1
    RE2292C 	 1
    RE0937C 	 0.0
    r1254 	 3.0
    LEUB0AT3tc 	 -1.0
    EX_caribupglu_S_LPAREN_e_RPAREN_ 	 1
    ACGAM6PSi 	 3.0
    3HCO3_NAt 	 -1.0
    ENGASE3ly 	 1
    GLYBtm 	 1
    M4BET2er 	 2.0
    DEDOLP1_U 	 1
    C4tmc 	 -1.0
    r2196 	 -1.0
    SPHINGStl 	 1
    ATVACIDhr 	 1
    FAOXC185C164m 	 3.0
    ALAR 	 1
    GGH_6DHFe 	 3.0
    EX_lstnm2_LPAREN_e_RPAREN_ 	 1
    LIMNENte 	 1
    ACOATA 	 0.0
    AM1ACCShr 	 1
    MHISOR 	 0.0
    RE2459C 	 1
    MINOHPtn 	 1
    RE2318X 	 3.0
    PHEtec 	 0.0
    PTE3x 	 0
    HDDACBP 	 3.0
    4HATVACIDthc 	 0.0
    APAT2rm 	 1.0
    ADPMAN 	 3.0
    ETOHMO 	 0.0
    r2124 	 2.0
    r2227 	 -1.0
    MACOXO 	 3.0
    ARACHFATPtc 	 0.0
    r0541 	 -1.0
    4BHGLZitr 	 1
    EX_3mlda_LPAREN_e_RPAREN_ 	 1
    r0819 	 1
    ALAt4 	 3.0
    13DMTtu 	 1
    r2225 	 -1.0
    RE2914M 	 3.0
    PTVSTM3itr 	 1
    RE3347C 	 1
    RE3532M 	 3.0
    r0947 	 -1.0
    UMPK2n 	 3.0
    r0086 	 3.0
    RE2677N 	 1
    r1612 	 0.0
    IMPC 	 3.0
    r1566 	 3.0
    EX_dha_LPAREN_e_RPAREN_ 	 1
    RE1711C 	 1
    RE3335R 	 3.0
    PYRSMCT1 	 -1.0
    r2294 	 -1.0
    TRIPVSteb 	 -1.0
    RE3573X 	 3.0
    FPGS 	 -1.0
    biomass_DNA 	 3
    HCOUMARINte 	 1
    FAOXC183_6Z_9Z_12Zx 	 -1.0
    EX_dcsptn1_LPAREN_e_RPAREN_ 	 1
    RE1796M 	 0.0
    FOLt2 	 -1.0
    r1552 	 3.0
    r1989 	 0.0
    RE3220L 	 0.0
    12HTACRtep 	 1
    RE1539X 	 -1.0
    r1024 	 1
    r0512 	 1.0
    EX_ala_L_LPAREN_e_RPAREN_ 	 3
    CYSGLYexR 	 0.0
    r1710 	 0.0
    EX_fucgalgbside_hs_LPAREN_e_RPAREN_ 	 1
    ENGASEly 	 1
    VACCCRNt 	 -1.0
    EX_txa2_LPAREN_e_RPAREN_ 	 1
    MAN1_6B1er 	 0.0
    STS1 	 0.0
    RE3624M 	 1
    FA161ACPH 	 -1.0
    RE0928C 	 -1.0
    FAOXC205C185m 	 3.0
    r2368 	 0
    HRETNtn 	 1
    r0245 	 3.0
    r1899 	 0.0
    RE3174C 	 3.0
    EX_btn_LPAREN_e_RPAREN_ 	 3
    EX_bvite_LPAREN_e_RPAREN_ 	 1
    EX_acgam_LPAREN_e_RPAREN_ 	 3
    RE2998M 	 -1.0
    NACHEX8ly 	 3.0
    UPPDC1 	 3.0
    EX_pydx5p_LPAREN_e_RPAREN_ 	 3
    r2297 	 -1.0
    METLEUex 	 0.0
    EX_ha_LPAREN_e_RPAREN_ 	 3
    ALDD2x 	 3.0
    EX_retfa_LPAREN_e_RPAREN_ 	 3
    RE1587R 	 0.0
    RE0828X 	 0.0
    ACSRTNMT 	 -1.0
    PI4P5Kn 	 1
    ILEtec 	 -1.0
    RE1235C 	 1
    PIK4n 	 1
    r2371 	 0.0
    r1965 	 0.0
    RAI3 	 3.0
    RE2375C 	 1
    GGH_7DHFe 	 3.0
    r0879 	 -1.0
    RE3111C 	 1
    r0892 	 1
    r0755 	 3.0
    RE3106C 	 1
    r0629 	 0.0
    RE3525X 	 1
    FATP4t 	 -1.0
    AM1CCSitr 	 1
    RE2319C 	 1
    RE3231R 	 1
    EX_HC01610_LPAREN_e_RPAREN_ 	 3
    EX_15dmt_LPAREN_e_RPAREN_ 	 1
    RE3148R 	 3.0
    EX_bz_LPAREN_e_RPAREN_ 	 1
    DURIK1m 	 -1.0
    RE3236R 	 1
    CERT1gt 	 -1.0
    BTNt2m 	 2.0
    RE1050E 	 -1.0
    FAOXC180x 	 -1.0
    ASNB0AT3tc 	 -1.0
    r0761 	 0.0
    r2395 	 -1.0
    RE2718C 	 0
    RE1952C 	 0.0
    PGLYCABCte 	 0.0
    r2011 	 -1.0
    LDH_D 	 -1.0
    TRDRm 	 0
    LNLCt 	 -1.0
    ENO 	 3.0
    r2293 	 -1.0
    r1685 	 0.0
    ST6GALNAC21 	 3.0
    RE3033C 	 1
    4OHMDZhr 	 1.0
    13DMTtep 	 1
    IPDDIx 	 3.0
    EX_fucgalfucgalacglcgalgluside_hs_LPAREN_e_RPAREN_ 	 1
    FAOXC225C226m 	 3.0
    PRAGSr 	 3.0
    RE2318C 	 1
    EPOXTACtev 	 1
    7THFtm 	 1
    EX_fucgal14acglcgalgluside_hs_LPAREN_e_RPAREN_ 	 1
    RE3421M 	 1
    RE1942C 	 1
    r0394 	 2.0
    D_3AIBt 	 1
    CLS_hs 	 3.0
    r0547 	 3.0
    r0570 	 3.0
    NACHEXA21ly 	 3.0
    r2437 	 -1.0
    OROTGLUt 	 -1.0
    AG13T14g 	 3.0
    r2353 	 0
    DASPO1p 	 -1.0
    DPROOp 	 -1.0
    EX_adn_LPAREN_e_RPAREN_ 	 3
    TSACMGLUCitr 	 1
    GLXO2p 	 -1.0
    P4507B11r 	 -1.0
    SUCCOAPET 	 0.0
    RE1233M 	 3.0
    FAOXC123_3Z_6Z_9Zm 	 3.0
    XYLTD_Dr 	 1
    r2361 	 0
    RE1447M 	 1
    LTC4Sr 	 3.0
    PI5P3K 	 3.0
    RETNtr 	 1
    FRUt1r 	 0
    C10CRNe 	 1
    NCNt 	 -1.0
    LVSTOXD3Hhep 	 1.0
    EX_cholate_LPAREN_e_RPAREN_ 	 1
    RE1942R 	 1.0
    RE3435C 	 1
    dcmpt 	 1
    DOCO13EFATP 	 0.0
    FAOXC163C164Gm 	 3.0
    GLUDC 	 1.0
    BAAT2x 	 0.0
    r2110 	 2.0
    GLCNACT5g 	 2.0
    ARTCOAL3 	 1
    DNDPt6m 	 -1.0
    RE3338X 	 3.0
    ACGBGBSIDEtl 	 1
    MI1456PKn 	 0.0
    RE3552X 	 1
    SR5AR2r 	 3.0
    EX_am1accs_LPAREN_e_RPAREN_ 	 1
    31DMThr 	 1
    PTE4x 	 0
    PROSTGE2t2 	 0.0
    EX_leuktrE4_LPAREN_e_RPAREN_ 	 3
    C16OHc 	 1.0
    CSPG_Ct 	 1
    CSNATp 	 -1.0
    SLDxm 	 3.0
    SERDGLNexR 	 -1.0
    GHMT2rm 	 3.0
    FAOXC183C163m 	 3.0
    HPACtr 	 1
    S2TASE5ly 	 1
    UMPtr 	 1
    PHYHx 	 1.0
    r1591 	 3.0
    RE3464C 	 1
    3HPVSCOAitm 	 1
    sink_tetdece1coa_LPAREN_c_RPAREN_ 	 1
    TMLYSOX 	 0.0
    GALASE8ly 	 0
    RE1907C 	 1
    EX_hpdca_LPAREN_e_RPAREN_ 	 1
    r2178 	 -1.0
    RE3533R 	 3.0
    r1392 	 3.0
    FORtr 	 1
    r0196 	 1.0
    HXPRT 	 3.0
    ATVLACitr 	 1
    RE3334M 	 1
    EX_cdpea_LPAREN_e_RPAREN_ 	 1
    ABO3g 	 -1.0
    RSVtu 	 0.0
    LPCOXp 	 -1.0
    6HMSMVACIDteb 	 1
    LTC4CP 	 1
    FAOXC81_5Zm 	 3.0
    RE3218R 	 0.0
    FACOAL245_1 	 3.0
    SEASMETtn 	 1
    r1532 	 1.0
    r1182 	 3.0
    ELAIDCPT1 	 1.0
    ASPTAm 	 3.0
    RE2994X 	 0.0
    RE1709C 	 1
    ACYP 	 0.0
    GALOR 	 3.0
    r2090 	 2.0
    EX_aflatoxin_LPAREN_e_RPAREN_ 	 1
    GPIDA2er 	 2.0
    C204CPT1 	 1.0
    LGNCCPT2 	 0.0
    RE0915E 	 0
    RE2318M 	 3.0
    RE3498N 	 2.0
    MI1PS 	 0.0
    INOSTO 	 -1.0
    LNSTLSr 	 3.0
    GP1Ctg 	 1
    VACCCPT1 	 1.0
    RE3151C 	 3.0
    RE3021C 	 1
    GLYC_St 	 1
    NACHEXA12ly 	 3.0
    DM_datp_m_ 	 1
    GCC2am 	 -1.0
    PI345P5P 	 1.0
    r1498 	 -1.0
    RE3344M 	 3.0
    RE3459C 	 1
    C160CPT1 	 1.0
    MI134PK 	 1
    ACGAM2E 	 -1.0
    LSTNM2itr 	 1
    DURItn 	 1
    RE3587C 	 1
    CSPG_Et 	 1
    RETFAt1 	 1
    PVD3 	 1
    NTD7e 	 1.0
    GLPASE1 	 3.0
    RE3225R 	 3.0
    FVSTETGLUitr 	 1
    HEX4 	 3.0
    GCALDD 	 3.0
    r1259 	 3.0
    r0093 	 -1.0
    ATVACIDhc 	 0.0
    CLOHtex2 	 -1.0
    r0319 	 0.0
    FADDP 	 -1.0
    r2306 	 -1.0
    r2115 	 2.0
    XOLTRIOLtr 	 1
    r2208 	 -1.0
    CYSAMOe 	 1
    EX_sprm_LPAREN_e_RPAREN_ 	 1
    FAOXC184x 	 1.0
    EX_ach_LPAREN_e_RPAREN_ 	 1
    RN0028X 	 0.0
    EX_hxan_LPAREN_e_RPAREN_ 	 3
    25HVITD3c 	 1.0
    RE1135C 	 0.0
    r1668 	 0.0
    EX_c6crn_ 	 1
    7AHGLZABCt 	 1
    FAOXC185m 	 3.0
    RE1807C 	 1
    r0743 	 3.0
    NTD9 	 3.0
    DRBK 	 -1.0
    r0002 	 1
    r1764 	 0.0
    SARCOXp 	 -1.0
    CSPG_Bt 	 1
    EX_HC02199_LPAREN_e_RPAREN_ 	 1
    CYTK12n 	 3.0
    H2Otg 	 1
    GLUtr 	 1
    EX_HC00229_LPAREN_e_RPAREN_ 	 1
    B3GNT36g 	 0.0
    RE3525R 	 1
    GLGNS1 	 0.0
    AM4NCStep 	 1
    RE3398M 	 3.0
    S3MEACMPhc 	 1
    r1815 	 0.0
    LEUKTRD4tr 	 1
    r2140 	 -1.0
    ASPt6 	 1.0
    RE1806R 	 0
    TTDCAtr 	 0.0
    r1825 	 0.0
    r2278 	 -1.0
    EX_dlnlcg_LPAREN_e_RPAREN_ 	 1
    r0692 	 0.0
    RE0916C 	 1
    r2472 	 -1.0
    RE1508C 	 1
    P45011A1m 	 -1.0
    RE0936C 	 0.0
    RE3017R 	 1
    UGLCNACtg 	 3.0
    SPTix 	 -1.0
    S4T2g 	 0.0
    r2355 	 0
    TCHOLAt2 	 -1.0
    RE2426C 	 1
    RE3444M 	 3.0
    TETPENT3CPT1 	 1.0
    RE3597M 	 0
    S2T4g 	 2.0
    r0456 	 1
    RE2477C 	 1
    RE3476C 	 3.0
    T4HCINNMFM 	 1
    GTHDH 	 3.0
    RE0569C 	 3.0
    r2100 	 2.0
    RE0570C 	 1
    EX_dhap_LPAREN_e_RPAREN_ 	 1
    CRVSM1hc 	 1.0
    EX_5htrp_LPAREN_e_RPAREN_ 	 1
    PHCHGSm 	 1
    RE2626M 	 0.0
    r0998 	 1
    RE3082X 	 3.0
    G14T16g 	 3.0
    P4502C18 	 0.0
    r2056 	 -1.0
    r1155 	 2.0
    GLCNACPT_U 	 1
    r1767 	 0.0
    GAO2g 	 3.0
    FAH3 	 -1.0
    PRISTtx 	 1
    ORNLEUrBATtc 	 -1.0
    CYTK4n 	 3.0
    3HLVSTAChep 	 0.0
    RE3430X 	 1
    CRVS23M24hc 	 1.0
    RE3243R 	 1
    r0120 	 3.0
    EX_hom_L_LPAREN_e_RPAREN_ 	 1
    RE3389M 	 1
    RE2605C 	 1
    r1179 	 3.0
    EX_HC02208_LPAREN_e_RPAREN_ 	 1
    EX_dopa_LPAREN_e_RPAREN_ 	 3
    ALOX52 	 2.0
    MTHFD2m 	 3.0
    BTNt2 	 2.0
    FAOXC164x 	 1.0
    DCT 	 -1.0
    EX_am1c4n9cs_LPAREN_e_RPAREN_ 	 1
    DOLP_Uter 	 1
    GLCGLUT2 	 -1.0
    r2026 	 -1.0
    GALASE9ly 	 0
    r0797 	 0.0
    r0153 	 -1.0
    FAOXC102_4Z_7Zx 	 -1.0
    r2369 	 0
    NADtx 	 1
    CRNtx 	 1
    RE3112C 	 1
    XOLTRI25te 	 1
    RE3568C 	 1
    URAt 	 -1.0
    r0121 	 3.0
    PTVSTM13hr 	 0.0
    DESAT22_1p 	 1
    HEXCCRNt 	 -1.0
    RE1582R 	 0.0
    RE2048R 	 1
    RE3270C 	 1
    RE3273C 	 -1.0
    RE1815X 	 3.0
    GCCbim 	 -1.0
    r0062 	 1
    DHCR243r 	 3.0
    35DSMVhep 	 1.0
    3HPVSitr 	 1
    SUCRe 	 0.0
    G16MT_U 	 1
    r2240 	 -1.0
    r1323 	 1
    CRVSM1SPhc 	 1
    HBZOPT10m 	 3.0
    PTVSTGLUChc 	 0
    3HPVSteb 	 -1.0
    PEFLIPm 	 1.0
    r2326 	 1.0
    3HLYTCL 	 0.0
    r2010 	 -1.0
    PRO1xm 	 1
    EX_cit_LPAREN_e_RPAREN_ 	 3
    NACHEX4ly 	 3.0
    r2324 	 1.0
    FAEL205 	 3.0
    RE3399M 	 3.0
    PYLALDOX 	 3.0
    S6TASE19ly 	 3.0
    r1026 	 -1.0
    PHEMEABCte 	 -1.0
    EX_tcynt_LPAREN_e_RPAREN_ 	 1
    r1449 	 0.0
    FAOXC5C5DCc 	 1.0
    AMETtd 	 1
    DEBRISOQUINEt 	 -1.0
    CYTK5n 	 3.0
    ST8SIA54g 	 -1.0
    TYRt 	 0.0
    NTPP9 	 3.0
    MDZtu 	 -1.0
    NDPK9n 	 0
    r1791 	 0.0
    ACMPthc 	 0.0
    RETt 	 -1.0
    r2513 	 1
    FCOAH 	 1
    AM1CSAitr 	 1
    DOLGPP_Uer 	 1
    TDPGDH 	 0.0
    r0239 	 -1.0
    FAOXC141_7Em 	 3.0
    GALT2g 	 2.0
    11DOCRTSTRNtr 	 1
    RE1516X 	 1
    IBUPGLUCtpvb 	 1
    GALASE13ly 	 0
    PGSr 	 2.0
    r1807 	 0.0
    RETNt 	 1
    EX_5mta_LPAREN_e_RPAREN_ 	 1
    r0408 	 3.0
    VITD3tm 	 1
    AM4N9CShr 	 1
    P4503A43r 	 -1.0
    GPDDA1 	 1
    S6T13g 	 0.0
    RE3337X 	 1.0
    RE3491C 	 1
    GLXtp 	 1
    r1862 	 0.0
    NAt3_1g 	 0.0
    NORANMT 	 -1.0
    r2390 	 0.0
    PAIL45P_HStn 	 1
    ARGNm 	 0.0
    TETPENT6CPT2 	 0.0
    RE1233C 	 -1.0
    r0614 	 -1.0
    S6T15g 	 0.0
    r0835 	 -1.0
    DOPAQNISO1 	 1
    VITD3t 	 1
    C162ACBP 	 3.0
    RE3220R 	 0.0
    NNDPR 	 1.0
    G14T21g 	 3.0
    r1682 	 0.0
    r2389 	 0.0
    THMATPe 	 1
    ALOX12R 	 1.0
    PECDCHe 	 1
    10FTHF6GLUtl 	 1
    RE2596C 	 1
    EX_glz_LPAREN_e_RPAREN_ 	 1
    FAOXC184_3Z_6Z_9Z_12Zm 	 3.0
    GK1 	 3.0
    It 	 -1.0
    FAOXC164GC163m 	 3.0
    LVSTtu 	 1
    GLYtm 	 1
    RE3346C 	 3.0
    GTPCIn 	 3.0
    EX_6hlvst_LPAREN_e_RPAREN_ 	 1
    3MOPte 	 1
    HMGLm 	 -1.0
    CYTK3n 	 3.0
    r2119 	 1
    EX_meracmp_LPAREN_e_RPAREN_ 	 1
    GPAM_hs 	 2.0
    RE3338C 	 3.0
    TMNDNCt 	 1
    ST6GALNAC27 	 3.0
    r1585 	 3.0
    PI4P5K 	 0.0
    r2275 	 -1.0
    r2046 	 -1.0
    EX_whddca_LPAREN_e_RPAREN_ 	 1
    RE3081X 	 1.0
    EX_12HPET_LPAREN_e_RPAREN_ 	 1
    CRVSitr 	 1
    r0388 	 -1.0
    EX_psylchol_LPAREN_e_RPAREN_ 	 1
    TAURCHAe 	 1
    r1974 	 0.0
    FADH2tx 	 1
    r2040 	 -1.0
    LYSMTF3n 	 3.0
    PROSTGI2t 	 1
    LTDCL 	 0.0
    RE3036C 	 1
    PPAt 	 1
    TACRitr 	 1
    RE3430M 	 1
    RTOTALCRNCPT2 	 0.0
    r2267 	 -1.0
    G14T8g 	 3.0
    DNDPt5m 	 -1.0
    AM4NCSitr 	 1
    HMGLx 	 -1.0
    RE3145X 	 0.0
    NDPK5 	 0.0
    r1745 	 0.0
    FOLR2 	 1.0
    RE3381C 	 0.0
    DECCRNe 	 1
    TKT1 	 3.0
    DURAD2 	 3.0
    AACOAT 	 3.0
    r2231 	 -1.0
    RE3385M 	 1
    PEPLYStn 	 1
    AM1CCSteb 	 1
    RE0579X 	 0.0
    URIt 	 -1.0
    FAOXC103C102x 	 -1.0
    RE3474C 	 1
    EX_HC02206_LPAREN_e_RPAREN_ 	 1
    7BHGLZGLCitr 	 1
    C10DCc 	 1.0
    r1765 	 0.0
    F6Tg 	 1.0
    EX_galacglcgalgbside_hs_LPAREN_e_RPAREN_ 	 1
    NACHEX12ly 	 3.0
    RE2990X 	 1.0
    UMPK7 	 3.0
    FAOXC122C101m 	 1.0
    RE2040C 	 1
    CYTDn 	 -1.0
    RE3564M 	 3.0
    r0934 	 3.0
    r1929 	 0.0
    ATPS4m 	 -1.0
    DM_atp_c_ 	 3
    r1979 	 0.0
    FAOXC163_7Z_10Z_13Zm 	 3.0
    FUT98g 	 -1.0
    r2215 	 -1.0
    r1251 	 3.0
    RE2626C 	 1
    3MEACMPhc 	 1
    EHGLAT 	 3.0
    LEUyLAThtc 	 0.0
    r1843 	 0.0
    THRCYSNaEx 	 0.0
    THP2Ctp 	 1
    PItn 	 1
    FADH2ETC 	 3.0
    Am1CSAteb 	 1
    EX_orn_LPAREN_e_RPAREN_ 	 3
    r1304 	 1
    DCATDc 	 0.0
    MI13456PK 	 1
    RE2637C 	 0.0
    r1922 	 0.0
    LCTStl 	 1
    NDPK8m 	 0.0
    PGK 	 3.0
    NACHEXA20ly 	 3.0
    NTMELYStner 	 1
    GLNTHRNaEx 	 0.0
    r0210 	 0.0
    FAOXC164_4Z_7Z_10Z_13Zx 	 -1.0
    3HPVSTETCOAhcm 	 3.0
    RE1521M 	 3.0
    3MEACMPitr 	 1
    r2047 	 -1.0
    RE3574X 	 0.0
    NT5C 	 3.0
    RE3626M 	 1
    UMPK6 	 3.0
    EX_sulpacmp_LPAREN_e_RPAREN_ 	 1
    RE3596C 	 3.0
    r1527 	 1
    r0681 	 0.0
    EX_tetpent3_LPAREN_e_RPAREN_ 	 1
    RE3408C 	 1
    CARBIBUP_SGLUthv 	 1
    FAOXC226205x 	 -1.0
    DNDPt48m 	 -1.0
    EX_pro_L_LPAREN_e_RPAREN_ 	 3
    CYSALANaEx 	 0.0
    DNAMTSEn 	 3.0
    TMDM1itr 	 1
    LCTStg 	 1
    NACHEX9ly 	 3.0
    r1586 	 3.0
    34DHPLACOX_NADP_ 	 3.0
    ABO8g 	 -1.0
    r2505 	 3.0
    1331TACRteb 	 1
    PERILLYLte 	 1
    RE3340M 	 3.0
    r2197 	 -1.0
    PDE4n 	 0.0
    RE1978C 	 1
    RN0030R 	 0.0
    EX_tre_LPAREN_e_RPAREN_ 	 3
    WHHDCAte 	 1
    SPODMx 	 3.0
    r0539 	 3.0
    MECOAS1m 	 1
    EX_pe_hs_LPAREN_e_RPAREN_ 	 3
    EX_spmd_LPAREN_e_RPAREN_ 	 3
    r1958 	 0.0
    r1589 	 3.0
    RE3273R 	 2.0
    GBA 	 0
    ALASERNaEx 	 0.0
    RE3372C 	 1
    EX_gln_L_LPAREN_e_RPAREN_ 	 3
    RE3092X 	 1.0
    DM_dttp_n_ 	 3
    DGNSKm 	 3.0
    r0678 	 0.0
    NADtn 	 1
    ASNGLNNaEx 	 0.0
    HEX10 	 3.0
    r1368 	 1
    RSVLACitr 	 1
    RE2335C 	 1
    RE3001M 	 3.0
    GTHP 	 3.0
    PETOHMr_hs 	 -1.0
    PRGNLONEtr 	 1
    PIK5 	 0.0
    RE3453C 	 1
    6CSMVACIDhep 	 0.0
    r2054 	 -1.0
    r2365 	 0
    r1887 	 0.0
    RE3184M 	 3.0
    RE1807M 	 0.0
    HS4ly 	 1.0
    RE0702M 	 1
    RE3261R 	 0
    EX_tymsf_LPAREN_e_RPAREN_ 	 1
    ORETNtn 	 1
    TTDCPT2 	 0.0
    EX_3ttetddcoacrn_ 	 1
    EX_fol_LPAREN_e_RPAREN_ 	 3
    r1846 	 0.0
    EX_asp_L_LPAREN_e_RPAREN_ 	 3
    r2123 	 1
    H2O2t 	 1
    RE0908G 	 0.0
    GALTg 	 -1.0
    PNTORDe 	 1
    r2187 	 -1.0
    CYTK11 	 3.0
    FAEL204 	 3.0
    r0714 	 1.0
    RE3241C 	 3.0
    r1378 	 3.0
    RE1846C 	 0.0
    RE1311M 	 1.0
    25HVITD3tin 	 -1.0
    r2486 	 -1.0
    r2387 	 0.0
    EX_cl_LPAREN_e_RPAREN_ 	 3
    NACHEX5ly 	 3.0
    VALTA 	 1.0
    RE3169C 	 1
    EX_for_LPAREN_e_RPAREN_ 	 3
    RN0022R 	 2.0
    DM_taur_LPAREN_c_RPAREN_ 	 1
    GUACYC 	 1.0
    r1740 	 0.0
    NTPP11 	 3.0
    EX_thm_LPAREN_e_RPAREN_ 	 3
    RE1309C 	 1
    FADH2tru 	 1
    RE3151R 	 3.0
    ADNK3 	 1
    RE3074X 	 3.0
    A4GALTg 	 -1.0
    6BHGLZABCt 	 1
    CYTK3 	 3.0
    RE3119C 	 1
    SRTNACT 	 -1.0
    PS_HSter 	 1
    RE3015C 	 1
    RE1954C 	 1
    6EPSteb 	 -1.0
    ECGISOr 	 1
    RE3125R 	 1
    RE3631C 	 1
    DHGLZhc 	 1
    EX_HC01104_LPAREN_e_RPAREN_ 	 1
    RE3018R 	 0.0
    FAOXC123_3Z_6Z_9Zx 	 3.0
    RE2249C 	 1
    r0973 	 1
    B3GALT3g 	 2.0
    NCKt 	 2.0
    r2148 	 -1.0
    r2048 	 -1.0
    THYPX 	 -1.0
    r0611 	 1
    MTRI 	 1
    r1520 	 3.0
    XYLt 	 0.0
    2HIBUPGLUC_Sitr 	 1
    r0595 	 3.0
    CSASULPhc 	 1
    31DMTitr 	 1
    ARACHCPT1 	 1.0
    NDPK10 	 0.0
    FAOXC226C225m 	 3.0
    RE1846X 	 0.0
    ACGAMK 	 3.0
    C51CPT1 	 1.0
    1531TACRtev 	 1
    EX_tmdm5_LPAREN_e_RPAREN_ 	 1
    P450LTB4r 	 1
    r2509 	 1
    GGH_6THFl 	 3.0
    ADPRDPm 	 2.0
    RETNGLCt 	 1
    FRDPtcr 	 1
    biomass_protein 	 3
    RE3629C 	 1
    r0317 	 3.0
    r0202 	 -1.0
    RE3432M 	 0
    r1646 	 3.0
    RE3086X 	 0.0
    ATPtm 	 3.0
    FAEL183 	 3.0
    EX_ptth_LPAREN_e_RPAREN_ 	 1
    AACTOOR 	 -1.0
    BAMPPALDOXm 	 3.0
    A4GALTc 	 -1.0
    RE2443M 	 3.0
    LEUKTRA4t 	 0
    FACOAL260 	 3.0
    UMPK3n 	 3.0
    EX_5mthf_LPAREN_e_RPAREN_ 	 1
    3HBCOAHLm 	 1
    TTDCAFATPtc 	 0.0
    r2174 	 -1.0
    r2301 	 -1.0
    r1164 	 1.0
    MOGAT 	 -1.0
    LSTNM5hr 	 1.0
    EX_pectintchol_LPAREN_e_RPAREN_ 	 1
    6AHGLZhr 	 1
    MI3456PK 	 1
    CRNtHa 	 0.0
    RE2029C 	 1
    EX_epoxtac_LPAREN_e_RPAREN_ 	 1
    r2247 	 -1.0
    PMI1346PH 	 0.0
    CRVS1tev 	 1
    ADNt 	 -1.0
    r1080 	 1
    FUCGALFUCGALACGLCGALGLUSIDEte 	 1
    1513DTACRhr 	 1
    LCADim 	 3.0
    PUNP5 	 3.0
    r0660 	 3.0
    DURAD 	 3.0
    DPPS 	 1
    RE0827C 	 -1.0
    ARTFR51 	 1
    EX_fad_LPAREN_e_RPAREN_ 	 1
    FAOXC2452253x 	 -1.0
    TYRDOPO 	 -1.0
    BPNT2 	 0.0
    EX_HC01609_LPAREN_e_RPAREN_ 	 3
    EX_gthox_LPAREN_e_RPAREN_ 	 3
    THSACMPitr 	 1
    r2492 	 -1.0
    MANtg 	 1
    4HATVACIDtep 	 1
    TYRDOPO3 	 -1.0
    r1797 	 0.0
    r2364 	 0
    RE3326M 	 -1.0
    RE2349M 	 0.0
    RE1533X 	 3.0
    r1831 	 0.0
    r1785 	 0.0
    FAOXC81C61m 	 1.0
    RE1860C 	 1
    r2190 	 -1.0
    BILIRUBt2 	 -1.0
    RE3272N 	 1
    RE3134C 	 0.0
    5MTAte 	 1
    S6TASE23ly 	 3.0
    PCHOLPr_hs 	 2.0
    CERT2rt 	 -1.0
    RE1812C 	 1
    S6TASE18ly 	 3.0
    MAL_Ltx 	 1
    RE3475C 	 1
    RE1514X 	 3.0
    r0752 	 3.0
    GALASE4ly 	 0
    BZt 	 1
    5HTRPVESSEC 	 0.0
    10FTHF7GLUtl 	 1
    r2251 	 -1.0
    r0433 	 1.0
    r1574 	 3.0
    CATm 	 3.0
    FACOAL226 	 3.0
    r2014 	 -1.0
    r0557 	 0
    RE3387M 	 3.0
    G1PTT 	 1
    RE2221M 	 0.0
    PI4P3K 	 3.0
    GLDBRAN 	 3.0
    RE3420C 	 1
    UDPRIBc 	 1
    EX_C04849_LPAREN_e_RPAREN_ 	 1
    AGTim 	 0.0
    PFK26 	 3.0
    EX_HC01440_LPAREN_e_RPAREN_ 	 1
    RN0021C 	 0.0
    r1521 	 3.0
    r2338 	 0.0
    THYMt 	 -1.0
    TMDM5itr 	 1
    DEOXFVSitx 	 1
    GALGALFUCFUCGALACGLCGALACGLCGAL14ACGLCGALGLUSIDEte 	 1
    EX_his_L_LPAREN_e_RPAREN_ 	 3
    CHOLtu 	 0.0
    r0737 	 0
    RE2081C 	 1
    GLACOm 	 3.0
    RE2235C 	 1
    S2T1g 	 0.0
    NRVNCCPT2 	 0.0
    XOL7AH2tr 	 1
    RE3234C 	 1
    DSREDUCr 	 3.0
    EX_lstn1gluc_LPAREN_e_RPAREN_ 	 1
    C161CPT1 	 1.0
    ANDRSTRNGLCte 	 3.0
    r2271 	 -1.0
    NMNATn 	 -1.0
    SUCOASm 	 2.0
    r1437 	 1
    CHATn 	 -1.0
    PYDAMK 	 3.0
    DEDOLP2_U 	 1
    CHOLD2m 	 -1.0
    TETTET6COAtx 	 1
    TMDPK 	 -1.0
    EX_digalsgalside_hs_LPAREN_e_RPAREN_ 	 1
    ACACT1rm 	 3.0
    FAOXC123C102m 	 1.0
    PVSGLUCteb 	 -1.0
    DRPA 	 3.0
    RE2360C 	 1
    URIt2r 	 0.0
    METtec 	 -1.0
    RETNtr2 	 1
    H8TAer 	 -1.0
    EX_HC02192_LPAREN_e_RPAREN_ 	 1
    AM1ALCSteb 	 1
    RE1827C 	 1
    1PPDCRp 	 1
    EX_hcoumarin_LPAREN_e_RPAREN_ 	 1
    SBTle 	 1
    2HATVACIDhc 	 0.0
    RE3557C 	 1
    RE2387C 	 1
    DCMPDA 	 3.0
    PHEMEe 	 1
    r2372 	 0.0
    EX_ak2lgchol_hs_LPAREN_e_RPAREN_ 	 1
    r1063 	 1
    ENMAN2g 	 1
    SCP22x 	 3.0
    r1963 	 0.0
    DHDPBMTm 	 -1.0
    GK1m 	 3.0
    DMHPTCRNte 	 1
    r1253 	 3.0
    r1411 	 0
    FACOAL205 	 3.0
    ALDD21 	 3.0
    MTHFTe 	 1
    RE3242R 	 1
    C120CRNt 	 -1.0
    CRVSM24itr 	 1
    r0984 	 -1.0
    C12DCACOT 	 -1.0
    DOLGLCP_Lter 	 1
    C8CRNe 	 1
    r2419 	 -1.0
    EX_galfuc12gal14acglcgalgluside_hs_LPAREN_e_RPAREN_ 	 1
    FAOXC101C8x 	 1.0
    r1669 	 0.0
    DHGLZABCt 	 1
    RE2897C 	 1
    ASPTA 	 3.0
    ESTRSABCtc 	 0.0
    RE3412C 	 1
    RE3567C 	 1
    SMPD4 	 -1.0
    EX_etoh_LPAREN_e_RPAREN_ 	 3
    r2402 	 -1.0
    N2M2NMASNt 	 1
    RETNGLCt2 	 1
    NADK 	 3.0
    NKCC2t 	 1.0
    BILDGLCURt 	 -1.0
    RE3227R 	 1
    r1945 	 0.0
    EX_HC02194_LPAREN_e_RPAREN_ 	 1
    OXYP2CONJ 	 1
    r1971 	 0.0
    EX_hx_LPAREN_e_RPAREN_ 	 1
    RN0032C 	 1
    FAOXC10080x 	 -1.0
    CRVStu 	 1
    OAADC 	 1
    ADPtx 	 1
    PSYTCHe 	 1
    GMPS2 	 3.0
    RE1538X 	 -1.0
    EX_HC02202_LPAREN_e_RPAREN_ 	 1
    UDPtl 	 1
    r1150 	 1
    r1633 	 3.0
    EX_fucfucgalacglcgalgluside_hs_LPAREN_e_RPAREN_ 	 1
    r1857 	 0.0
    C81CPT1 	 1.0
    RE2644C 	 -1.0
    4ABUTtcn 	 1
    FAOXC101m 	 3.0
    DCYTD 	 0.0
    r1981 	 0.0
    2HATVACIDtep 	 1
    PSP_L 	 0.0
    C142OHe 	 1
    FUCGALFUCGALACGLCGALGLUSIDEtg 	 1
    DM_dgtp_m_ 	 1
    ASPte 	 1
    STRDNCCRNt 	 -1.0
    EX_CE5798_LPAREN_e_RPAREN_ 	 1
    r0579 	 3.0
    RE2722C 	 0
    r2388 	 0.0
    SARCStp 	 1
    DNDPt53m 	 -1.0
    G5SADrm 	 1
    GLXO1 	 3.0
    r1691 	 0.0
    EX_whtststerone_LPAREN_e_RPAREN_ 	 1
    PItx 	 1
    RE2327C 	 1
    ALDD2y 	 3.0
    r2084 	 1
    r2017 	 -1.0
    RE1518X 	 1
    RE3224R 	 3.0
    r1667 	 0.0
    PSt3 	 1
    RETH2 	 1
    r1571 	 3.0
    C226CPT2 	 0.0
    r0181 	 2.0
    IMPD 	 3.0
    C161CRN2t 	 -1.0
    NACHORCTL3le 	 -1.0
    FAOXC6040m 	 3.0
    RE0830C 	 1
    r2101 	 1
    RE3500R 	 2.0
    r1303 	 1
    YVITEt 	 1
    r0594 	 3.0
    ACGALtly 	 1
    r1851 	 0.0
    RE3500X 	 0.0
    EGMESTr 	 0.0
    EX_glygn5_LPAREN_e_RPAREN_ 	 1
    DOLPGT1_Ler 	 3.0
    5HOXINDACTO2OX 	 -1.0
    r0220 	 0.0
    AVITE2t 	 1
    PYDXDH 	 -1.0
    TRIODTHYSUFt 	 1
    r0792 	 -1.0
    RE3519C 	 0.0
    PROD2m 	 0
    AVITE1t 	 -1.0
    EX_2mcit_LPAREN_e_RPAREN_ 	 1
    EAFLATOXINte 	 1
    RE3445M 	 1
    r0354 	 3.0
    DHAPtc 	 1
    GBGT1 	 -1.0
    GALGT2 	 0.0
    GNDc 	 3.0
    FK 	 -1.0
    CYTDt4 	 0.0
    DM_yvite_c_ 	 1
    DTMPKm 	 0.0
    EX_10fthf_LPAREN_e_RPAREN_ 	 1
    r0812 	 1
    NDP8 	 0.0
    q10h2tc 	 1
    TSTSTERONEtr 	 1
    DNDPt30m 	 -1.0
    MI1PP 	 3.0
    RETH 	 1.0
    RN0028R 	 2.0
    TETPENT6CRNt 	 -1.0
    ACALDtr 	 1
    RE2034C 	 1
    FAOXC185_3Z_6Z_9Z_12Z_15Zx 	 3.0
    BETALDHxm 	 1
    EX_crmp_hs_LPAREN_e_RPAREN_ 	 1
    ALAB0AT3tc 	 -1.0
    GLPASE2 	 3.0
    GTMLTe 	 -1.0
    PRFGS 	 0.0
    r1642 	 3.0
    FAOXC7050m 	 3.0
    PROt2r 	 0.0
    P5CR 	 -1.0
    4HBZFm 	 1
    PYK 	 3.0
    FAOXC162C142m 	 3.0
    MI145PK 	 1.0
    r2151 	 -1.0
    VACCCPT2 	 0.0
    SMS 	 1.0
    RE3335C 	 -1.0
    PI4P3Ker 	 3.0
    CPCTDTX 	 1.0
    FAOXC102C103m 	 3.0
    NDPK10m 	 0.0
    EX_ppi_LPAREN_e_RPAREN_ 	 1
    FAOXC181_9Zm 	 3.0
    ACSm 	 0.0
    14MDZtev 	 1
    GLYOX 	 -1.0
    RE1234C 	 1
    EX_din_LPAREN_e_RPAREN_ 	 3
    EX_dhglz_LPAREN_e_RPAREN_ 	 1
    NDPK3 	 0.0
    M16NTg 	 3.0
    ARTFR12 	 1
    S6TASE25ly 	 -1.0
    RE2992M 	 3.0
    RE2660N 	 1.0
    GTHS 	 3.0
    r1980 	 0.0
    FAOXC180 	 0.0
    RE3226R 	 3.0
    EX_ps_hs_LPAREN_e_RPAREN_ 	 3
    FAS160COA 	 0.0
    r1547 	 3.0
    r0686 	 -1.0
    BAAT4x 	 0.0
    Kt3g 	 0.0
    EX_bhb_LPAREN_e_RPAREN_ 	 3
    SERtN1 	 -1.0
    r0970 	 1
    OMEPRAZOLEte 	 1
    ORNtiDF 	 3.0
    Coqe 	 1
    RE0583N 	 1
    r1551 	 3.0
    CRVNCtr 	 1
    r2263 	 -1.0
    r1533 	 1.0
    4OHPROIMINOtc 	 -1.0
    r1367 	 1
    GLNyLATthc 	 0.0
    RE2975M 	 1.0
    ACMPGLUChr 	 0
    PYRt2p 	 2.0
    r2093 	 1
    r2518 	 3.0
    RE3345R 	 3.0
    THMPPtm 	 1
    EX_tripvs_LPAREN_e_RPAREN_ 	 1
    RE2520C 	 1
    r1735 	 0.0
    r2206 	 -1.0
    RE2386C 	 1
    1513TACRitr 	 1
    AG13T13g 	 3.0
    RE3500C 	 0.0
    PPAn 	 1
    GLYPHEHYc 	 1
    LGNCCPT1 	 1.0
    FPGS3m 	 -1.0
    r0768 	 0.0
    PYRt2r 	 2.0
    FAOXC204C205x 	 -1.0
    HSD3B13 	 0.0
    ST3GAL62g 	 0.0
    r2200 	 -1.0
    RE2916M 	 3.0
    r1861 	 0.0
    r1171 	 1.0
    SARCStex 	 1
    MGSA 	 1
    r2370 	 -1.0
    PYAM5POr 	 -1.0
    BTNTe 	 1
    FAOXC182_9E_12Em 	 3.0
    HTDCRNe 	 1
    OXYPR7tehv 	 1
    NAPQIhr 	 1
    r1262 	 3.0
    r1829 	 0.0
    7HPVShc 	 1.0
    EX_docosac_LPAREN_e_RPAREN_ 	 1
    EX_HC02214_LPAREN_e_RPAREN_ 	 1
    r1994 	 -1.0
    DESAT18_4 	 3.0
    PTVSTtep 	 1
    r1018 	 1
    4MPTNLtr 	 1
    RE3457C 	 1
    ACHEe 	 -1.0
    FE2t 	 1
    AMETtn 	 1
    DALAt2r 	 0.0
    FPGS5m 	 -1.0
    RE2718G 	 -1.0
    1MNCAMti 	 1
    DTMPK 	 0.0
    ASAH1 	 3.0
    EX_3mob_LPAREN_e_RPAREN_ 	 3
    RE3242C 	 1
    GPIDAer 	 2.0
    GLYK 	 0.0
    RE3488C 	 1
    RE3004M 	 3.0
    r2165 	 -1.0
    FBA 	 3.0
    TRIODTHYSULT 	 2.0
    r1907 	 0.0
    RE1828M 	 0.0
    RETFA 	 1
    RN0020C 	 1
    TLACFVStev 	 1
    RE0919C 	 -1.0
    DNDPt33m 	 1
    PROIMINOtc 	 -1.0
    RE3308M 	 1
    RE3515C 	 1
    RE3124C 	 3.0
    EX_mdz_LPAREN_e_RPAREN_ 	 1
    SADT 	 3.0
    RE3229C 	 1
    3SALAASPm 	 2.0
    6HMSMVitr 	 1
    EX_pan4p_LPAREN_e_RPAREN_ 	 1
    MDZitr 	 1
    FPGS9m 	 -1.0
    EX_2hibupglu_S_LPAREN_e_RPAREN_ 	 1
    CATp 	 3.0
    RE1100R 	 0.0
    RE0066R 	 -1.0
    RE3526M 	 1
    LVSTOXD6METhep 	 1.0
    EX_nfdlac_LPAREN_e_RPAREN_ 	 1
    CYSSNAT4te 	 -1.0
    N3Tg 	 0.0
    CHAT 	 -1.0
    NTD3l 	 2.0
    RE3152C 	 1
    SRTNtu 	 0.0
    EX_acmpglu_LPAREN_e_RPAREN_ 	 1
    AG13T17g 	 3.0
    PYDXNtr 	 1
    RE3520N 	 1
    MALTSULtm 	 -1.0
    SUCCACT 	 1
    r1140 	 1
    r2122 	 2.0
    DOGULNO1 	 1
    FUC13GALACGLCGAL14ACGLCGALGLUSIDEtg 	 1
    SUCCt4_2 	 -1.0
    GHMT3 	 -1.0
    S4TASE4ly 	 -1.0
    ENGASE2 	 3.0
    S2TASE1ly 	 3.0
    RADH2 	 1
    r0166 	 3.0
    RE3006M 	 3.0
    GSNtl 	 -1.0
    EX_HC01787_LPAREN_e_RPAREN_ 	 1
    Htg 	 1
    METB0AT2tc 	 -1.0
    r0389 	 -1.0
    ST6GALNAC22 	 3.0
    EBASTINEOHte 	 1
    ARAB_Lt 	 1
    r1935 	 0.0
    GGT6 	 1
    r1554 	 3.0
    FUC14GALACGLCGALGLUSIDEtg 	 1
    r2299 	 -1.0
    RE2916X 	 3.0
    FPGSm 	 -1.0
    GLUDxm 	 3.0
    r0603 	 -1.0
    C10OHc 	 1.0
    RE0579M 	 3.0
    PROB0AT2tc 	 1
    C141CPT1 	 1.0
    EX_HC02203_LPAREN_e_RPAREN_ 	 1
    GD1B2te 	 1
    6CSMVACIDteb 	 1
    C12OHc 	 1.0
    FAOXC164C163x 	 -1.0
    EX_fum_LPAREN_e_RPAREN_ 	 3
    RE2987X 	 3.0
    ESTRONESt2 	 -1.0
    CRTSTRNtr 	 1
    r2363 	 0
    5HOXINOXDA 	 3.0
    EX_gtp_LPAREN_e_RPAREN_ 	 1
    ATVACYLGLUChc 	 0
    EX_35dhpvs_LPAREN_e_RPAREN_ 	 1
    ATPtx 	 1
    TDP 	 1
    r1806 	 0.0
    NACHEX22ly 	 3.0
    UAG4E 	 2.0
    EX_onpthl_LPAREN_e_RPAREN_ 	 1
    CSNAT3x 	 -1.0
    FAOXC161C141x 	 -1.0
    THRt4 	 3.0
    FAOXC16BRx 	 0.0
    r1776 	 0.0
    LFORKYNHYD 	 3.0
    r1670 	 0.0
    GNDer 	 3.0
    ABO7g 	 -1.0
    RE2513L 	 3.0
    NTD9e 	 1.0
    EX_urate_LPAREN_e_RPAREN_ 	 1
    EX_lvstacid_LPAREN_e_RPAREN_ 	 1
    EX_cmp_LPAREN_e_RPAREN_ 	 1
    SR5ARr 	 3.0
    FAOXC184_3Z_6Z_9Z_12Zx 	 3.0
    RE3178M 	 3.0
    TETPENT3CPT2 	 0.0
    MI145P6Kn 	 0.0
    RE3225C 	 3.0
    RE2156M 	 1
    3MOX4HOXPGALDOX_NADP_ 	 3.0
    r1144 	 1
    EX_tagat_D_LPAREN_e_RPAREN_ 	 1
    CYANt 	 1
    RE1920C 	 1
    r0310 	 3.0
    r1302 	 1
    SQLEr 	 3.0
    FAOXC161_7Zm 	 3.0
    RE2680C 	 1
    CLI2tex 	 0.0
    2MB2COAc 	 3.0
    RE3388M 	 3.0
    GLUPRT 	 0.0
    r0280 	 -1.0
    AHCYStn 	 1
    4OHMDZitr 	 1
    RBK 	 -1.0
    RE3095L 	 3.0
    RE3339X 	 0.0
    r1676 	 0.0
    r2392 	 0.0
    r1828 	 0.0
    FUMSO3tm 	 -1.0
    RE2675C 	 1
    r0829 	 -1.0
    r0650 	 -1.0
    PSFLIPm 	 1.0
    B3GNT51g 	 3.0
    r1783 	 0.0
    MLTG1 	 0.0
    CRTSLt 	 1
    NDPK10n 	 0
    r2171 	 -1.0
    6THFtm 	 1
    DCTPtn 	 1
    KSIt 	 1
    GGT5r 	 -1.0
    CEPTC 	 3.0
    r2358 	 0
    FAOXC141C141OHm 	 3.0
    RE2080C 	 1
    r1695 	 0.0
    2HATVACIDitr 	 1
    B3GNT11g 	 3.0
    EX_crm_hs_LPAREN_e_RPAREN_ 	 1
    AM1ACSteb 	 1
    r0709 	 3.0
    PAN4PPe 	 1
    DM_bvite_c_ 	 1
    FAOXC140120x 	 -1.0
    SCP21x 	 3.0
    H7ET2er 	 2.0
    r2125 	 2.0
    RE1317C 	 1
    LDH_Lm 	 3.0
    NACSMCTte 	 -1.0
    r2376 	 0.0
    LRAT2 	 1
    r1014 	 1
    ARACHDt2 	 -1.0
    C60CPT2 	 0.0
    RE2028C 	 1
    DM_datp_n_ 	 3
    HEX1 	 3.0
    2HATVACIDGLUCitr 	 1
    DM_btn 	 1
    ACTLMO 	 0.0
    EX_HC01444_LPAREN_e_RPAREN_ 	 1
    EX_3thexddcoacrn_ 	 1
    ALDSTRNte 	 1
    STACMPhc 	 1
    RE1537C 	 1
    r2207 	 -1.0
    LNLNCGt 	 -1.0
    r2449 	 0.0
    EX_glyc3p_LPAREN_e_RPAREN_ 	 3
    LSTNM4hr 	 1
    r1581 	 3.0
    EX_adp_ 	 1
    EX_prpp_LPAREN_e_RPAREN_ 	 1
    EX_tdcrn_ 	 1
    EX_lnlncg_LPAREN_e_RPAREN_ 	 3
    r1701 	 0.0
    EX_met_L_LPAREN_e_RPAREN_ 	 3
    EX_am1a4ncs_LPAREN_e_RPAREN_ 	 1
    RE3239C 	 1
    RE3432C 	 3.0
    EX_lys_L_LPAREN_e_RPAREN_ 	 3
    r0287 	 3.0
    SELADT 	 3.0
    RE3198C 	 1
    PMANM 	 3.0
    r2034 	 -1.0
    TAURtcx 	 3.0
    MCOATA 	 0.0
    RE1901R 	 1
    RETH1e 	 1
    r0130 	 1.0
    RE2722G 	 -1.0
    EX_gsn_LPAREN_e_RPAREN_ 	 3
    BTNPLm 	 -1.0
    r2245 	 -1.0
    RE2972M 	 0.0
    r1445 	 0.0
    RE0569E 	 1
    EX_ebastineoh_LPAREN_e_RPAREN_ 	 1
    EX_3ohacmp_LPAREN_e_RPAREN_ 	 1
    25HVITD2tm 	 1
    3MOXTYROX 	 3.0
    CORE4GTg 	 1.0
    DNDPt36m 	 1
    r1156 	 0
    PTVSTM13itr 	 1
    RE3334X 	 0.0
    LPSe 	 1.0
    EX_gltcho_LPAREN_e_RPAREN_ 	 1
    RE1526X 	 3.0
    DM_n5m2masn_g_ 	 1
    PPA2 	 1
    HXANtl 	 -1.0
    r2223 	 -1.0
    HPDCACRNCPT1 	 1.0
    FAOXC2252053m 	 3.0
    ACGAMtly 	 1
    S3TASE2ly 	 1
    PRGNLONEtm 	 1
    EX_cgly_LPAREN_e_RPAREN_ 	 3
    15DMThr 	 1
    RE3020R 	 0.0
    NACHEXA19ly 	 3.0
    RE2250C 	 1
    TMDOATthc 	 -1.0
    H2MTer_U 	 1
    RE3556C 	 1
    r0764 	 0.0
    M1316Mg 	 3.0
    HPYRR2x 	 3.0
    RE3134R 	 0.0
    ST8SIA55g 	 -1.0
    S6TASE5ly 	 -1.0
    RTOTALCRNt 	 -1.0
    r2221 	 -1.0
    NACHEXA11ly 	 3.0
    UDPGALt2g 	 1
    RE3532R 	 3.0
    EX_fvstetglu_LPAREN_e_RPAREN_ 	 1
    GUAPRT 	 3.0
    AM1CSAhr 	 1.0
    EX_2hatvlac_LPAREN_e_RPAREN_ 	 1
    r1774 	 0.0
    GLYB0AT3tc 	 -1.0
    PECGCHLe 	 1
    RE1819X 	 0
    RE3554M 	 3.0
    RE3447M 	 3.0
    LYStm 	 -1.0
    PTVSTLAChc 	 1
    URATEtx 	 1
    LSTNM7TDhc 	 1
    DM_sprm_c_ 	 1
    r1782 	 0.0
    EX_am4n9cs_LPAREN_e_RPAREN_ 	 1
    AACTtm 	 1
    EX_ptdca_LPAREN_e_RPAREN_ 	 1
    r2237 	 -1.0
    Q10H2e 	 1
    PROGLYPRO1c 	 -1.0
    5OHFVSGLUhc 	 1
    NACHEX15ly 	 3.0
    AGPAT1 	 3.0
    RETI1 	 1
    GASNASE2ly 	 3.0
    PAPSitr 	 1
    RE3289C 	 1
    EX_lnlc_LPAREN_e_RPAREN_ 	 3
    RN0002R 	 1
    Tyr_ggnt 	 1
    GLCNACT2g 	 2.0
    DPMVDc 	 -1.0
    GLNSERNaEx 	 0.0
    r1830 	 0.0
    ALOX5 	 2.0
    RE1525X 	 3.0
    THRGLYexR 	 -1.0
    GALASE19ly 	 0
    RE2067C 	 1
    UMPK5 	 3.0
    RE3168C 	 1
    COAtr 	 1
    RN0002N 	 1
    r0766 	 0.0
    M16N6Tg 	 3.0
    RE3033N 	 2.0
    LPS2e 	 -1.0
    DM_T_antigen_g_ 	 1
    NACHEX13ly 	 3.0
    SUCCTD 	 1
    NAGA2ly 	 0.0
    EX_nad_LPAREN_e_RPAREN_ 	 1
    RE3169R 	 1
    ALAtN1 	 -1.0
    r2065 	 -1.0
    RE2993X 	 -1.0
    r1811 	 0.0
    r2335 	 1.0
    Asn_X_Ser_Thrtr 	 1
    r0734 	 1.0
    ESTRGLCABCCte 	 -1.0
    AG13T15g 	 3.0
    RE2915X 	 0.0
    r0968 	 1
    RE1917C 	 1
    S6TASE8ly 	 -1.0
    ACETONEt2 	 2.0
    PGISr 	 -1.0
    EX_lstnm1_LPAREN_e_RPAREN_ 	 1
    AHANDROSTANGLCtr 	 1
    r2535 	 -1.0
    CYTDtm 	 -1.0
    ALLOP1tu 	 0.0
    FAOXC183C163Gm 	 3.0
    HSAT4ly 	 1
    RE2849C 	 1
    ODECOAtx 	 1
    DESAT18_10 	 0.0
    GCHOLAt 	 0.0
    EX_glypro_LPAREN_e_RPAREN_ 	 1
    HSD17B8r 	 1.0
    r0440 	 1
    MAN1PT2 	 0
    PYDX5Ptm 	 1
    DHAPAx 	 3.0
    BAAT1x 	 0.0
    RE2306C 	 0.0
    1531TACRhr 	 1
    MTHFD 	 3.0
    EX_caribup_s_LPAREN_e_RPAREN_ 	 1
    RE2888E 	 -1.0
    FORMCOAtx 	 1
    RE3636C 	 1
    r1569 	 3.0
    PROSTGE1t 	 1.0
    TALA 	 3.0
    HDCEAtr 	 0.0
    ESTRAABCtc 	 0.0
    FAOXCPRIST2x 	 -1.0
    RE3235R 	 1
    PSSA1_hs 	 3.0
    FAOXC102C101x 	 -1.0
    5OXPROt 	 1
    RE3571R 	 2.0
    PPNCL3 	 2.0
    r0129 	 1.0
    C16DCc 	 1.0
    r0360 	 3.0
    RE3175C 	 1
    RE1238X 	 1
    GLCAASE1ly 	 3.0
    TTDCRNt 	 -1.0
    FACOAL2252 	 3.0
    RN0001C 	 1
    q10tm 	 1
    2HATVACIDGLUCteb 	 -1.0
    EX_hco3_LPAREN_e_RPAREN_ 	 1
    r1932 	 0.0
    ACt2m 	 1
    ACGPID 	 -1.0
    r0555 	 1.0
    EX_lipoate_LPAREN_e_RPAREN_ 	 1
    RE1819C 	 3.0
    RE3310R 	 1
    EX_ibupgluc_LPAREN_e_RPAREN_ 	 1
    FMNALKPle 	 3.0
    RNDR1 	 3.0
    NRPPHRt4_2_r 	 -1.0
    r2186 	 -1.0
    AM9CSAteb 	 1
    FAOXC2242046m 	 3.0
    PTE2x 	 0
    3HLVSTitr 	 1
    PNP 	 3.0
    RE2868C 	 1
    GAPD 	 3.0
    FUT96g 	 -1.0
    r2359 	 0
    r2511 	 1
    RE3172C 	 1
    6HTSTSTERONEte 	 1
    RE1918C 	 1
    r1877 	 0.0
    A4GNT2g 	 -1.0
    RE0452M 	 1
    PSYTDECHe 	 1
    GLCNACT3g 	 2.0
    RE1310C 	 1
    S26Tg 	 2.0
    LYSOXp 	 1
    CRGLZitr 	 1
    r2341 	 0
    4HMDZGLUChr 	 0
    SELNPS 	 3.0
    FAOXC102_4Z_7Zm 	 3.0
    ACITL 	 3.0
    PSHSABCtc 	 3.0
    ATVLACh2r 	 1
    GLAl 	 3.0
    MHGLZitr 	 1
    GMPR 	 0.0
    SIAASE4ly 	 -1.0
    S3MEACMPtev 	 1
    PHEATB0tc 	 3.0
    G14T15g 	 3.0
    6THFtl 	 1
    MI4PP 	 3.0
    FAOXC140120m 	 3.0
    PDX5PO 	 -1.0
    VD3 	 1
    LSTNM2tev 	 1
    AM9CSAtep 	 1
    DGK2m 	 1
    34DHOXPEGOX 	 3.0
    RE2265C 	 1
    NTD10 	 3.0
    RE3381L 	 3.0
    25HVITD2tin 	 1
    RE3560M 	 3.0
    GALGALGALTHCRMtg 	 1
    10FTHF7GLUtm 	 1
    EX_galgalgalthcrm_hs_LPAREN_e_RPAREN_ 	 1
    13DAMPPOX 	 0.0
    RE3344X 	 1.0
    GTACMPtev 	 1
    RE3177M 	 3.0
    GCHOLAt3 	 0.0
    13DMThr 	 1
    FAOXC200180m 	 3.0
    6HMSMVhep 	 1.0
    RE3110C 	 3.0
    FVSTETGLUhc 	 1
    M13N2Tg 	 -1.0
    4HATVLACOXDhc 	 1.0
    RN0023X 	 0.0
    RE3260R 	 0
    r1650 	 3.0
    PI45P3K 	 2.0
    FAOXC161_9Em 	 3.0
    RE3166R 	 1
    RE3526X 	 1
    EX_6bhglz_LPAREN_e_RPAREN_ 	 1
    RE2513E 	 -1.0
    FAOXC160140m 	 3.0
    r0924 	 -1.0
    DOLPMT_U 	 1
    r0735 	 -1.0
    SMPD3l 	 -1.0
    r1652 	 3.0
    PIACGT 	 -1.0
    RE2859C 	 1
    MTHFDm 	 3.0
    r2277 	 -1.0
    EX_ctp_LPAREN_e_RPAREN_ 	 1
    DM_gpi_sig_er_ 	 1
    GLUNm 	 2.0
    GALASE3ly 	 0
    EX_tacr_LPAREN_e_RPAREN_ 	 1
    PTCRTD 	 1
    GCCam 	 -1.0
    RAI2 	 1
    RE3423C 	 1
    RE2717L 	 1.0
    ADRNCOAtx 	 1
    RE1933C 	 1
    RE3434C 	 1
    r0988 	 1
    LVACLAChep 	 1
    COUCOAFm 	 1
    GALtly 	 1
    RE2068C 	 1
    L_LACtcm 	 0.0
    r1750 	 0.0
    RSVtev 	 1
    FAOXC3DC 	 1.0
    O2Stn 	 1
    HISTAtu 	 1.0
    RE0928R 	 0
    INSTt4 	 3.0
    S6TASE21ly 	 3.0
    CLOXAtex2 	 -1.0
    NTP3e 	 3.0
    DNDPt51m 	 -1.0
    EX_thf_LPAREN_e_RPAREN_ 	 1
    EX_ura_LPAREN_e_RPAREN_ 	 3
    RE3518R 	 2.0
    EX_glygn4_LPAREN_e_RPAREN_ 	 1
    ACN13ACNGALGBSIDEte 	 1
    HPYRtp 	 1
    GGH_7DHFl 	 3.0
    EX_retinol_9_cis_LPAREN_e_RPAREN_ 	 1
    RE2050C 	 1
    SERHL 	 0.0
    r2254 	 -1.0
    FAOXC12DCx 	 1.0
    GD1Cte 	 1
    DNDPt21m 	 -1.0
    DALAt2rL 	 0.0
    FOLTle 	 -1.0
    RE1063C 	 1
    sink_pre_prot_LPAREN_r_RPAREN_ 	 1
    DM_5hpet_LPAREN_r_RPAREN_ 	 1
    CPPPGO 	 0.0
    AKGDm 	 0
    C121CPT1 	 1.0
    GLYt2rL 	 0.0
    LAPCOAl 	 2.0
    ETOHtx 	 1
    PPD2CSPp 	 1
    r0142 	 1.0
    r1443 	 0.0
    DNADtn 	 1
    C102CPT1 	 1.0
    RETFAt 	 1
    FPGS9 	 -1.0
    RBFK 	 0.0
    SULFOX 	 1.0
    r1708 	 0.0
    CSNATm 	 -1.0
    NAIt 	 -1.0
    RE3444C 	 0
    FAOXC5030m 	 3.0
    STS2r 	 0.0
    PMEVKx 	 0.0
    r1906 	 0.0
    FUCFUC132GALACGLCGAL14ACGLCGALGLUSIDEtg 	 1
    S6TASE15ly 	 3.0
    UDPGNP 	 1
    EX_estriolglc_LPAREN_e_RPAREN_ 	 1
    r0438 	 1
    r1957 	 0.0
    RIBt2 	 1
    MAN2_6B1er 	 1
    ACS2 	 3.0
    GLCSGLT1le 	 0.0
    RETI2 	 1
    EX_glygly_LPAREN_e_RPAREN_ 	 1
    RE1816C 	 1
    RE1941C 	 1
    ACOAO7p 	 -1.0
    NAt5 	 -1.0
    RE2146R 	 0
    5OHFVSGLUitr 	 1
    DNDPt44m 	 -1.0
    RE2869C 	 1
    7HPVSteb 	 -1.0
    r2006 	 -1.0
    MI1P_Dtn 	 1
    DPGM 	 3.0
    r0772 	 0.0
    RE1062M 	 3.0
    r0224 	 1.0
    EX_oxa_LPAREN_e_RPAREN_ 	 3
    EX_sfcys_LPAREN_e_RPAREN_ 	 1
    GLCATg 	 -1.0
    RE1691M 	 3.0
    RE3395M 	 1
    r1673 	 0.0
    r0817 	 1
    PTDCAt 	 1
    MDHx 	 1
    r1792 	 0.0
    GLCter 	 1
    GLCt2_2 	 0.0
    G14T13g 	 3.0
    AGLPET 	 1
    RE3488N 	 2.0
    RE3559X 	 3.0
    ALCD21_L 	 3.0
    r2438 	 -1.0
    H4ET3er 	 0.0
    RE3470C 	 3.0
    r0443 	 1
    RE2921X 	 1.0
    r1454 	 2.0
    ALAGLYexR 	 -1.0
    RE2439C 	 1
    r0616 	 -1.0
    r0775 	 3.0
    H3MTer_U 	 1
    RE3190M 	 3.0
    SO4HCOtex 	 2.0
    r1665 	 0.0
    RE2476C 	 1
    SO4t4_2 	 0.0
    r1004 	 1
    EX_omeprazole_LPAREN_e_RPAREN_ 	 1
    EX_glyc_S_LPAREN_e_RPAREN_ 	 1
    PIPLC 	 1.0
    DMGDHm 	 -1.0
    r0719 	 3.0
    GLCAT2g 	 -1.0
    GLCURter 	 1
    4NPHte 	 1
    34DHPHAMT 	 3.0
    RE2522C 	 1
    r1777 	 0.0
    35DHPVSitr 	 1
    3DPHBH2 	 1
    r1926 	 0.0
    r0009 	 3.0
    SPODMe 	 -1.0
    RE0919R 	 0
    EX_xolest2_hs_LPAREN_e_RPAREN_ 	 1
    LNLCCPT1 	 1.0
    EX_udpgal_LPAREN_e_RPAREN_ 	 1
    S6T7g 	 0.0
    EX_gbside_hs_LPAREN_e_RPAREN_ 	 1
    PGMT 	 3.0
    HEXCOAACBP 	 3.0
    r1590 	 3.0
    RE3488R 	 2.0
    r2068 	 -1.0
    GCCcm 	 -1.0
    GLYtp 	 1
    RETNCOA 	 1
    r0801 	 3.0
    RE3525N 	 1
    VALtec 	 -1.0
    ST8SIA51g 	 -1.0
    RETABCtc 	 0.0
    r1365 	 1
    PTVSTM3hc 	 1.0
    PHEB0AT3tc 	 -1.0
    UDPGLDCg 	 3.0
    r1742 	 0.0
    r1052 	 1
    RE2051G 	 3.0
    r2333 	 1.0
    GLYSARPEPT1tc 	 -1.0
    r2109 	 1
    DHPM1 	 3.0
    PMTFATPtc 	 1
    CYSACMPitr 	 1
    LCADi_Dm 	 3.0
    RE3486C 	 1
    r1655 	 3.0
    ARGt4 	 3.0
    FAOXC81C61x 	 -1.0
    B3GNT314g 	 0.0
    C6DCCACT 	 -1.0
    DM_fol 	 1
    NADN 	 1
    B3GNT31g 	 0.0
    H2Otp 	 1
    r2192 	 -1.0
    r1736 	 0.0
    TSTSTERONESte 	 1
    RE3220C 	 0.0
    r1386 	 1
    r0463 	 3.0
    SIAASE2ly 	 -1.0
    DM_ethamp_r_ 	 1
    r0937 	 1
    PYNP2r 	 2.0
    FTHFDH 	 0.0
    MALTt1r 	 0.0
    DCTPtm 	 1
    ATPH1e 	 2.0
    r1067 	 1
    NACHEX1ly 	 3.0
    r1003 	 1
    FAS140COA 	 0.0
    NS26Tg 	 2.0
    FA1822ACPH 	 1
    RE3185M 	 3.0
    r0637 	 1.0
    r2019 	 -1.0
    GLYVESSEC 	 -1.0
    r1664 	 0.0
    AM19CSitr 	 1
    HAS2 	 3.0
    3DSPHR 	 2.0
    TETTET6CPT2 	 0.0
    ENGASE 	 3.0
    GLCAASE9ly 	 3.0
    GCALDDm 	 3.0
    CYTK13n 	 3.0
    PCHOL_HSter 	 1
    PAPStg 	 0
    r1944 	 0.0
    TRDR 	 3.0
    NACHEX11ly 	 3.0
    H2O2tn 	 1
    P45027A11m 	 0.0
    CO2t 	 1
    r2031 	 -1.0
    ASPNATm 	 1
    LINKDEG2ly 	 1
    METB0AT3tc 	 -1.0
    PVSGLUCtev 	 1
    r2169 	 -1.0
    r0724 	 3.0
    r1382 	 -1.0
    r2049 	 -1.0
    T2M26DCOAHLx 	 3.0
    biomass_reaction 	 3
    3MLDAt 	 1
    EX_atp_LPAREN_e_RPAREN_ 	 3
    RE3164R 	 1
    r2113 	 2.0
    RE1804C 	 1
    S6TASE12ly 	 3.0
    GULN3D 	 1
    EX_retinal_LPAREN_e_RPAREN_ 	 1
    NADPN 	 1
    ACNAM9PL 	 3.0
    LSTNM7thc 	 1
    FAOXC225_4Z_7Z_10Z_13Z_16Zx 	 -1.0
    C6DCe 	 1
    SMVitr 	 1
    DNDPt35m 	 1
    FACOAL80i 	 2.0
    r2332 	 1.0
    r2106 	 1
    SOAT12r 	 1
    r1517 	 3.0
    EX_crvs_LPAREN_e_RPAREN_ 	 1
    EX_triodthy_LPAREN_e_RPAREN_ 	 1
    EX_Tyr_ggn_LPAREN_e_RPAREN_ 	 1
    r1609 	 3.0
    EX_estroneglc_LPAREN_e_RPAREN_ 	 1
    S3T1g 	 2.0
    r2373 	 0.0
    r0575 	 3.0
    r2098 	 2.0
    r1824 	 0.0
    CHOLtg 	 1
    ACGAMPM 	 -1.0
    HSD17B7r 	 2.0
    EX_duri_LPAREN_e_RPAREN_ 	 1
    r1159 	 1
    DM_pmtcoa_LPAREN_r_RPAREN_ 	 1
    ALAPAT4te 	 -1.0
    TMNDNCCPT2 	 0.0
    NADNe 	 0.0
    3HSMVitr 	 1
    10FTHFtm 	 1
    r1325 	 1
    AHC 	 3.0
    IPDPtr 	 1
    ASPDTDe 	 1
    r2141 	 -1.0
    r1684 	 0.0
    RIBFLVt3 	 0.0
    RE3586X 	 0.0
    S2TASE4ly 	 1
    RE3434R 	 2.0
    34HPPte 	 1
    r0422 	 3.0
    3HBCDm 	 3.0
    EX_2hb_LPAREN_e_RPAREN_ 	 3
    SERGLNNaEx 	 0.0
    KSItly 	 1
    RE2452C 	 1
    DNDPt8m 	 -1.0
    DLNLCGCPT1 	 1.0
    P45017A4r 	 0.0
    FAOXC221201x 	 -1.0
    r2375 	 0.0
    EX_acetone_LPAREN_e_RPAREN_ 	 3
    CYTDt5le 	 0.0
    O16G1e 	 0.0
    RE0923C 	 -1.0
    r1976 	 0.0
    6MELVACtbc 	 -1.0
    TMDitr 	 1
    RE3409C 	 1
    r0425 	 2.0
    7BHGLZitr 	 1
    ADSELK 	 3.0
    r1017 	 1
    HSD11B2r 	 0.0
    CO2tm 	 1
    FA180ACPH 	 -1.0
    EX_eicostet_LPAREN_e_RPAREN_ 	 1
    RAtn3 	 1
    RE1826C 	 1
    r1448 	 0.0
    ESTRADIOLGLCtr 	 1
    S6TASE6ly 	 -1.0
    EX_ile_L_LPAREN_e_RPAREN_ 	 3
    ICDHxm 	 2.0
    4OHMDZtev 	 1
    VITD2tm 	 1
    DPCOAtl 	 1
    G6PDH1rer 	 0.0
    RE3307C 	 3.0
    PIPLCn 	 1.0
    GLCAE1g 	 1
    ABUTD 	 3.0
    RE3308X 	 1
    HDECEACBP 	 3.0
    CHOLATEt3 	 0.0
    C8DCe 	 1
    TSULt4_3 	 -1.0
    TYMSFt 	 1
    TDCHOLAte 	 1
    RE3113C 	 1
    r0716 	 1.0
    r0475 	 1.0
    EX_ncam_LPAREN_e_RPAREN_ 	 3
    RE1827M 	 0.0
    GLYOXm 	 -1.0
    RE1522X 	 1.0
    r1696 	 0.0
    RE2522X 	 3.0
    RE3343C 	 0
    FAOXC8060m 	 3.0
    r2030 	 -1.0
    1HMDGLUChr 	 0
    GT1Ate 	 1
    ADPACDAc 	 1
    PTHPS 	 3.0
    RE0452N 	 1
    r2158 	 -1.0
    EX_tchola_LPAREN_e_RPAREN_ 	 1
    GRTTx 	 3.0
    RE2625C 	 1
    STRDNCCPT2 	 0.0
    AG13T18g 	 3.0
    PEFLIP 	 1.0
    RE2443C 	 1
    RE2649M 	 3.0
    COAtg 	 1
    EX_hdca_LPAREN_e_RPAREN_ 	 3
    C204CRNt 	 -1.0
    DHAPA 	 3.0
    6BHGLZhr 	 0.0
    GALACGLCGALGBSIDEtg 	 1
    RE0688X 	 0.0
    RE1957G 	 0.0
    ACRNtm 	 -1.0
    GALt2_2 	 0.0
    DMATT 	 0.0
    ORNt4m 	 -1.0
    r2139 	 -1.0
    CYOOm2 	 0
    RE2407C 	 1
    LCYSTATm 	 3.0
    SERCYSNaEx 	 0.0
    RE2041C 	 1
    GBSIDEte 	 1
    RN0022C 	 0.0
    DM_ncam 	 1
    ABUTt2rL 	 0.0
    GALFUC12GAL14ACGLCGALGLUSIDEte 	 1
    EX_dchac_LPAREN_e_RPAREN_ 	 1
    r0717 	 1.0
    r2241 	 -1.0
    AKR1D 	 -1.0
    EX_tethex3_LPAREN_e_RPAREN_ 	 1
    RE2453M 	 3.0
    THMP 	 1
    NOS1 	 -1.0
    RE2155R 	 0
    RE3248M 	 3.0
    RE2387R 	 0.0
    r1672 	 0.0
    MANtly 	 1
    r2347 	 0
    CORE6GTg 	 0.0
    UREAt5 	 0.0
    CBR1 	 2.0
    r1493 	 1
    r1725 	 0.0
    NTD4l 	 2.0
    LEUKTRA4tr 	 1
    LST4EXPTDhc 	 1
    NADKm 	 3.0
    r2514 	 1
    RE2649C 	 0
    UDPACGALtl 	 1
    PGESr 	 2.0
    DNDPt4m 	 -1.0
    MTHFC 	 3.0
    PNTEHe 	 1
    TSACMGLUChr 	 1
    NACHEXA2ly 	 3.0
    PNTEH 	 0.0
    EX_thmtp_LPAREN_e_RPAREN_ 	 1
    GP1Cte 	 1
    CSPG_Dt 	 1
    r2153 	 -1.0
    EX_C02470_LPAREN_e_RPAREN_ 	 1
    B3GNT310g 	 0.0
    r0788 	 1
    EPCTX 	 -1.0
    ASCBt4 	 0.0
    RE3576X 	 0.0
    GDPtg 	 1
    RE3485N 	 2.0
    r0522 	 -1.0
    DM_13_cis_retn_n_ 	 1
    PYDXtr 	 1
    BILIRUBtr 	 1
    EX_xyl_D_LPAREN_e_RPAREN_ 	 3
    RE2799C 	 1
    r1583 	 3.0
    EX_pheacgln_LPAREN_e_RPAREN_ 	 1
    RE3335X 	 3.0
    GGH_5THFl 	 3.0
    ALASm 	 1.0
    OCDCEAFATPtc 	 -1.0
    DNDPt54m 	 -1.0
    NMNATr 	 0.0
    NDPK6n 	 0
    r0757 	 3.0
    AMANK 	 3.0
    FAOXC182C162m 	 3.0
    THYOXt2 	 0.0
    RE2069C 	 1
    P45027A12m 	 0.0
    r1788 	 0.0
    P4502F1 	 -1.0
    RE2129C 	 1
    EX_HC02205_LPAREN_e_RPAREN_ 	 1
    AGLPT 	 1
    RE2333C 	 1
    PGLYCt 	 1
    r1841 	 0.0
    P4502C94 	 0.0
    r1707 	 0.0
    r2189 	 -1.0
    GGH_10FTHF7GLUe 	 3.0
    BTNt3i 	 0.0
    GLYCTDle 	 1
    r1639 	 3.0
    MALTe 	 -1.0
    RE3577X 	 3.0
    RE3370C 	 1
    EX_leuktrB4_LPAREN_e_RPAREN_ 	 1
    LNELDCCRNt 	 -1.0
    FAOXC121C101m 	 1.0
    H6MTer_U 	 1
    EX_ala_D_LPAREN_e_RPAREN_ 	 3
    EX_csasulp_LPAREN_e_RPAREN_ 	 1
    EX_bilirub_LPAREN_e_RPAREN_ 	 1
    RE2972G 	 0.0
    EX_am1alcs_LPAREN_e_RPAREN_ 	 1
    ME2m 	 3.0
    r1892 	 0.0
    r1628 	 3.0
    ETOHt 	 1
    FAOXC101_4Zm 	 3.0
    RE3012M 	 3.0
    RE3411C 	 1
    RE3370R 	 0.0
    FACOAL240 	 3.0
    HISDC 	 0.0
    EX_c10crn_ 	 1
    C5DCe 	 1
    L_LACtm 	 2.0
    r1640 	 3.0
    r1456 	 1
    FAOXC1811601m 	 3.0
    C4CRNe 	 1
    4HATVLACitr 	 1
    RE3122R 	 1
    EX_debrisoquine_LPAREN_e_RPAREN_ 	 1
    PPItx 	 1
    RE3226C 	 3.0
    GUMDCHAe 	 1
    r0163 	 1.0
    r0301 	 3.0
    PHEMEtm 	 1
    FAS100COA 	 0.0
    RE1310M 	 1.0
    DNDPt41m 	 -1.0
    ACCOAC 	 3.0
    RE3160R 	 3.0
    r2128 	 1
    r0950 	 1
    UTPtn 	 1
    CYTK10n 	 3.0
    EX_1glyc_hs_LPAREN_e_RPAREN_ 	 1
    FAOXC12DCc 	 1.0
    RE1135G 	 0.0
    G12MT2_U 	 1
    FAOXC161_7Em 	 3.0
    r2070 	 -1.0
    r1737 	 0.0
    FAOXC143_5Z_8Z_11Zm 	 3.0
    RE3120C 	 3.0
    DNDPt40m 	 -1.0
    CLCFTRte 	 0.0
    EX_4abutn_LPAREN_e_RPAREN_ 	 1
    r2434 	 -1.0
    EX_galfucgalacglcgal14acglcgalgluside_hs_LPAREN_e_RPAREN_ 	 1
    NRVNCt 	 1
    ACHtn 	 1
    r0580 	 1
    r2088 	 2.0
    RE2898C 	 1
    r2252 	 -1.0
    GLYSNAT5tc 	 -1.0
    ONPTHLte 	 1
    6HSMVitr 	 1
    C181CPT2 	 0.0
    r0932 	 1
    GALC 	 1.0
    r1178 	 3.0
    P45011B12m 	 -1.0
    5HOMEPRAZOLEte 	 1
    C181OHc 	 1.0
    RE1809C 	 1
    TETPENT3CRNt 	 -1.0
    RE0918R 	 0.0
    RE3163R 	 1
    r1781 	 0.0
    RE3121C 	 1
    CDSm 	 3.0
    r0630 	 0.0
    GD1Ctg 	 1
    3MOPt2im 	 1
    NACHEX24ly 	 3.0
    S2T2g 	 0.0
    NDP7ex 	 3.0
    r1731 	 0.0
    RE0066C 	 1
    GALSIDEtl 	 1
    5OHFVShc 	 1.0
    C4DCe 	 1
    r2532 	 -1.0
    RE1811C 	 1
    EX_elaid_LPAREN_e_RPAREN_ 	 1
    EX_thymd_LPAREN_e_RPAREN_ 	 3
    r1744 	 0.0
    CYTD 	 0.0
    r0024 	 1
    RE3218C 	 0.0
    RE1957R 	 1.0
    ATPH2e 	 2.0
    DALAOXx 	 -1.0
    r1028 	 -1.0
    VITD3Hm 	 0.0
    DCSPTN1CPT1 	 1.0
    r0680 	 3.0
    SBTD_D2 	 3.0
    FPGS7m 	 -1.0
    LSTNM4itr 	 1
    r1167 	 1.0
    RE2860C 	 1
    HIVCACBP 	 3.0
    CRNATBtc 	 3.0
    ARABR 	 3.0
    EX_3aib_LPAREN_e_RPAREN_ 	 1
    RE3273G 	 2.0
    G3M8MASNterg 	 1
    CHTNASEe 	 -1.0
    G6PDH2r 	 0.0
    r0796 	 1
    ALLOPOXDhep 	 2.0
    DAGKn_hs 	 -1.0
    T4HCINNOX 	 -1.0
    HPYRRy 	 3.0
    PI34P3Pn 	 1
    r0708 	 3.0
    RE2624M 	 -1.0
    CO2tg 	 1
    EX_pectingchol_LPAREN_e_RPAREN_ 	 1
    AM1CCShr 	 1
    RPI 	 2.0
    LACZe 	 0
    4HATVACIDOXDhc 	 1.0
    FAOXC162_7E_10Em 	 3.0
    FERO 	 0.0
    KYNATESYN 	 1
    RN0029C 	 1
    RE2523C 	 1
    RE3340C 	 3.0
    EX_bilglcur_LPAREN_e_RPAREN_ 	 1
    GGNG 	 0.0
    RE2638X 	 0.0
    r1175 	 3.0
    UMPKm 	 1
    MI14Ptn 	 1
    4BHGLZtev 	 1
    RE3158X 	 3.0
    THMMPt4 	 -1.0
    r1896 	 0.0
    FAEL184 	 3.0
    r2285 	 -1.0
    RE1916C 	 1
    EX_ac_LPAREN_e_RPAREN_ 	 3
    AG13T16g 	 3.0
    ARGtiDF 	 3.0
    r1930 	 0.0
    CAATPS 	 3.0
    r1908 	 0.0
    DNDPt61m 	 -1.0
    ALADGLNexR 	 -1.0
    NADPNe 	 0.0
    r2510 	 1
    G6PDA 	 3.0
    r0444 	 1
    11DOCRTSLtm 	 1
    S3T3g 	 1.0
    NACHEXA9ly 	 3.0
    DNDPt27m 	 1
    HMGCOARr 	 3.0
    TAGHSTDe 	 1
    EX_fald_LPAREN_e_RPAREN_ 	 3
    TCHOLABCtc 	 0.0
    RE3066X 	 1.0
    RE2304E 	 1.0
    RE2440C 	 1
    DATPtn 	 1
    TYR3MO2 	 0.0
    PAPitr 	 1
    FAOXC2252053x 	 -1.0
    r1400 	 1
    GTHRDt 	 1
    LCADi_D 	 3.0
    RE3265C 	 -1.0
    ABO6g 	 -1.0
    r0549 	 3.0
    MMSAD1m 	 0.0
    AMPTASECGe 	 0.0
    D_LACt2 	 2.0
    HISyLATtc 	 0.0
    EX_h2o_LPAREN_e_RPAREN_ 	 3
    r1931 	 0.0
    XYLUR 	 2.0
    r0960 	 0.0
    AICARte 	 1
    GABABGTtc 	 -1.0
    DAG_HSter 	 1
    XAO2x 	 2.0
    GCC2bim 	 -1.0
    GLCAE2g 	 0.0
    r1074 	 1
    r1967 	 0.0
    AM1ACCStev 	 1
    r1897 	 0.0
    SERGLYexR 	 -1.0
    HSD3B11r 	 -1.0
    7BHGLZGLChr 	 1
    RE2677R 	 2.0
    LALDO 	 3.0
    PHYCBOXL 	 0.0
    14HMDZitr 	 1
    SERtp 	 1
    r2346 	 0
    r2242 	 -1.0
    r2286 	 -1.0
    r1895 	 0.0
    EX_2hatvlacgluc_LPAREN_e_RPAREN_ 	 1
    EX_hestratriol_LPAREN_e_RPAREN_ 	 1
    AM1ACCSitr 	 1
    7BHGLZhr 	 0.0
    EX_dca_LPAREN_e_RPAREN_ 	 1
    r2229 	 -1.0
    FUCtly 	 1
    EX_HC02197_LPAREN_e_RPAREN_ 	 1
    CHLP 	 -1.0
    PSSA2_hs 	 -1.0
    C16txc 	 -1.0
    GLYGLYPEPT1tc 	 -1.0
    r1817 	 0.0
    EX_CE2838_LPAREN_e_RPAREN_ 	 1
    2HATVLACOXDhc 	 1.0
    MM5cg 	 3.0
    r2036 	 -1.0
    DM_core7_g_ 	 1
    EX_HC00955_LPAREN_e_RPAREN_ 	 1
    FUT34g 	 1.0
    NNATn 	 -1.0
    2HATVACIDthc 	 0.0
    C161CPT22 	 0.0
    EX_adrn_LPAREN_e_RPAREN_ 	 1
    EX_dheas_LPAREN_e_RPAREN_ 	 1
    AMPtp 	 1
    r1148 	 1
    THRALANaEx 	 0.0
    HSPGtly 	 1
    C4STMO1r 	 3.0
    RE2128C 	 1
    GALGALGALTHCRMte 	 1
    FUCACNGAL14ACGLCGALGLUSIDEte 	 1
    EX_man_LPAREN_e_RPAREN_ 	 3
    EX_dpcoa_LPAREN_e_RPAREN_ 	 1
    RN0021R 	 2.0
    EX_andrstrn_LPAREN_e_RPAREN_ 	 1
    SERDGLYexR 	 -1.0
    r1000 	 1
    DMGtm 	 1
    FAOXOHC16C16DCc 	 3.0
    SMVthep 	 -1.0
    CRVSM23teb 	 -1.0
    FAOXC163x 	 1.0
    RE3039C 	 1
    r0688 	 3.0
    LPS3e 	 -1.0
    RE3167C 	 1
    RE2220C 	 1
    r0917 	 0.0
    XYLtly 	 1
    EX_dag_hs_LPAREN_e_RPAREN_ 	 3
    ARACHt 	 1
    ADPRIBt 	 0.0
    EX_tauribup_S_LPAREN_e_RPAREN_ 	 1
    LSTNM1hr 	 1.0
    PPPGOm 	 0.0
    RE3403M 	 3.0
    RE2117M 	 -1.0
    r0653 	 3.0
    MALT 	 0.0
    AM1ACShr 	 1
    RE3337M 	 3.0
    3HPVShc 	 1.0
    FAOXOHC22C22DCc 	 3.0
    FUCACGALFUCGALACGLCGALGLUSIDEte 	 1
    6HLVSTAChep 	 0.0
    RE3165R 	 3.0
    EX_6ohfvs_LPAREN_e_RPAREN_ 	 1
    RE3562C 	 -1.0
    FAOXC165C164m 	 3.0
    RE3307M 	 0
    NDPK6m 	 0.0
    NACHEX10ly 	 3.0
    r2320 	 1.0
    NACHEX3ly 	 3.0
    r2298 	 -1.0
    NACHEXA15ly 	 3.0
    r0340 	 -1.0
    RE3470M 	 0
    r1611 	 3.0
    RSVATPtu 	 -1.0
    r1886 	 0.0
    RE1077C 	 1
    LEUTAm 	 0.0
    r1553 	 3.0
    R_group_phosphotase_3 	 1
    EX_phyt_LPAREN_e_RPAREN_ 	 1
    CYTK11n 	 3.0
    r2009 	 -1.0
    LSTNitr 	 1
    NTD8 	 3.0
    CSAitr 	 1
    NRVNCCRNt 	 -1.0
    RE1904R 	 0.0
    HKYNH 	 3.0
    FPGS8 	 -1.0
    RE3367E 	 0.0
    INStl 	 -1.0
    ABO1g 	 -1.0
    NMNtn 	 1
    3HPVSTETteb 	 -1.0
    RE3485C 	 1
    GALGLUSIDEtg 	 1
    ST6GALNAC28 	 3.0
    RE2429M 	 -1.0
    EX_paf_hs_LPAREN_e_RPAREN_ 	 1
    FAOXC61C4x 	 1.0
    KSII_CORE2t 	 1
    r0634 	 3.0
    DCK1m 	 1
    RE1538C 	 1
    NDPK3n 	 0
    ORNDC 	 3.0
    XOLTRI25tc 	 1
    FUMm 	 3.0
    EX_atvlac_LPAREN_e_RPAREN_ 	 1
    r1801 	 0.0
    r1168 	 1.0
    RE2240C 	 1
    r2408 	 -1.0
    EX_strch1_LPAREN_e_RPAREN_ 	 3
    r1044 	 0.0
    DGTPtm 	 1
    GLCAT6g 	 2.0
    CAt7r 	 -1.0
    EX_6epvs_LPAREN_e_RPAREN_ 	 1
    TCYNTt 	 1
    r0173 	 3.0
    CYTDtl 	 -1.0
    EX_dmantipyrine_LPAREN_e_RPAREN_ 	 1
    AHANDROSTANGLCte 	 3.0
    56EPPVShc 	 1.0
    GLYCLTtp 	 1
    NDERSVitr 	 1
    r2471 	 0.0
    r2538 	 1
    P45027A13m 	 0.0
    5HOXINDACTOXm 	 3.0
    RE3407C 	 1
    6OHFVSteb 	 1
    EX_mhglz_LPAREN_e_RPAREN_ 	 1
    HOMt4 	 3.0
    FAOXOHMC10DC10c 	 3.0
    UDPGLCAtg 	 1
    CSPG_Ctly 	 1
    NACHEXA5ly 	 3.0
    EX_3hlvstacid_LPAREN_e_RPAREN_ 	 1
    AM4NC9CSteb 	 1
    6OHFVSGLUtev 	 1
    NRPPHRSFt 	 1
    2HATVLACtep 	 1
    OCT11EFATP 	 0.0
    S6TASE2ly 	 3.0
    GSNtm 	 -1.0
    r2264 	 -1.0
    RE3224C 	 3.0
    EX_pyr_LPAREN_e_RPAREN_ 	 3
    RE3140X 	 1.0
    FBP 	 0.0
    FTCD 	 -1.0
    RE3288R 	 0.0
    r2052 	 -1.0
    RE3521C 	 1
    ATVLACThc 	 0.0
    S2L2N2M2MASNtly 	 1
    ALCD2yf 	 3.0
    EX_3aib_D_LPAREN_e_RPAREN_ 	 1
    GALASE20ly 	 0
    r1441 	 1
    PRO_Dtde 	 1
    RE1303C 	 1
    EX_arg_L_LPAREN_e_RPAREN_ 	 3
    ARACHCOAtx 	 1
    B3GNT33g 	 0.0
    RE3019C 	 1
    r2300 	 -1.0
    FAOXC6C4x 	 -1.0
    PEt 	 1
    GD1B2tg 	 1
    EX_HC02213_LPAREN_e_RPAREN_ 	 1
    SGALSIDEtg 	 1
    EX_malthp_LPAREN_e_RPAREN_ 	 1
    ILEB0AT2tc 	 1
    ACMPShc 	 2.0
    RE3045C 	 1
    RE1815R 	 3.0
    NH4tn 	 1
    RE3448M 	 3.0
    ST6GALNAC26 	 3.0
    S2L2FN2M2MASNtly 	 1
    RE2768R 	 0.0
    NDPK4n 	 0
    CORE2GTg 	 1.0
    RE2867C 	 1
    BDMT_U 	 1
    EX_thr_L_LPAREN_e_RPAREN_ 	 3
    r2173 	 -1.0
    S6TASE9ly 	 -1.0
    ALAASNNaEx 	 0.0
    r0784 	 1
    S6T4g 	 0.0
    RE1796C 	 1
    RE3287C 	 1
    RE2272C 	 1
    FE2DMT1 	 3.0
    P4503A5 	 0
    r1924 	 0.0
    r1565 	 3.0
    r1562 	 3.0
    EX_succ_LPAREN_e_RPAREN_ 	 1
    AG13T2g 	 3.0
    13DMTitr 	 1
    RE0938C 	 0.0
    FOAXC122C101x 	 1.0
    RE0690E 	 0.0
    H2ETer 	 0.0
    FAOXC121_3Zm 	 3.0
    r0558 	 1
    EX_5oxpro_LPAREN_e_RPAREN_ 	 1
    r1627 	 3.0
    EX_3tetd7ecoacrn_ 	 1
    r1436 	 1
    RE3578X 	 1.0
    FAOXC164m 	 3.0
    3ISPVSitr 	 1
    r1257 	 3.0
    BVITEt 	 1
    CITtbm 	 0.0



```python
Recon2.optimize()
Recon2.summary(fva=True)
```

    IN FLUXES                                      OUT FLUXES                                    OBJECTIVES
    ---------------------------------------------  --------------------------------------------  -----------------------
    id                   Flux  Range               id                   Flux  Range              biomass_reac...  0.0384
    ---------------  --------  ------------------  ---------------  --------  -----------------
    o2_e             0.348     [-0.401, 1]         pi_e             0.28      [-0.04, 2.69]
    fe2_e            0.186     [-0.012, 0.545]     fe3_e            0.186     [-0.012, 0.545]
    atp_e            0.103     [-0.11, 0.8]        nh4_e            0.134     [-0.012, 3.24]
    ocdca_e          0.07      [0, 0.07]           pyr_e            0.122     [-0.012, 6.75]
    so4_e            0.05      [0.05, -0.178]      lac_L_e          0.11      [-0.09, 8.57]
    ala_D_e          0.04      [0.04, -3.21]       ura_e            0.0761    [-0.012, 1.53]
    lys_L_e          0.03      [0.0227, 0.03]      co2_e            0.0703    [0, 9.12]
    Lcystin_c        0.024     [0.024, -0.0652]    vacc_e           0.07      [0, 1.23]
    utp_e            0.02      [0.02, -0.89]       xmp_e            0.0608    [-0.012, 1.63]
    val_L_e          0.0135    [0.00153, 0.03]     slfcys_e         0.0601    [0, 0.114]
    gln_L_e          0.0125    [0.09, -1.54]       acald_e          0.0277    [-0.012, 10]
    thr_L_e          0.012     [0.012, 0.012]      strdnc_e         0.0252    [0, 0.488]
    leu_L_e          0.012     [0.00894, 0.012]    gua_e            0.0222    [-0.012, 1.39]
    crn_e            0.012     [0, 0.012]          aicar_e          0.0207    [0, 1.75]
    o2s_e            0.012     [0, 0.012]          datp_m           0.0191    [0, 0.91]
    sbt_D_e          0.012     [0, 0.012]          HC01609_e        0.012     [-0.012, 0.012]
    2hb_e            0.012     [-0.00613, 0.012]   3bcrn_e          0.00728   [0, 0.012]
    tyr_L_e          0.012     [-0.00792, 0.012]   hdcea_e          0.00649   [-0.012, 1.88]
    HC01610_e        0.012     [-0.012, 0.012]     ac_e             0.00621   [-0.012, 10]
    HC00250_e        0.012     [0.012, -0.166]     gthox_e          0.006     [0, 0.0892]
    cgly_e           0.012     [0.012, -0.166]     34hpp_e          0.00587   [0.00792, -0.012]
    cys_L_e          0.012     [0.012, -0.166]     glyb_e           0.0054    [-0.012, 1.16]
    lnlnca_e         0.012     [0.012, -0.684]     3ivcrn_e         0.00472   [0, 0.012]
    arg_L_e          0.012     [0.012, -0.807]     for_e            0.000783  [-0.012, 1.37]
    dctp_n           0.012     [0.012, -0.898]     dag_hs_e         0.000671  [-0.012, 0.871]
    dgtp_n           0.012     [0.012, -0.898]     h_e              0         [-1, 10]
    cytd_e           0.012     [0.012, -1.07]      hco3_e           0         [0, 9.12]
    dcmp_e           0.012     [0.012, -1.07]      mthgxl_e         0         [0, 8.66]
    dcyt_e           0.012     [0.012, -1.07]      glyc_e           0         [-0.21, 8.06]
    chol_e           0.012     [0.012, -1.24]      glyc_S_e         0         [0, 6.76]
    gsn_e            0.012     [0.012, -1.39]      bhb_e            0         [-0.012, 6.1]
    uri_e            0.012     [0.012, -1.45]      succ_e           0         [0, 5.79]
    ptrc_e           0.012     [0.012, -1.62]      acetone_e        0         [-0.012, 5.05]
    ocdcea_e         0.012     [0.012, -1.66]      drib_e           0         [-0.012, 4.93]
    ddca_e           0.012     [0.012, -1.84]      acac_e           0         [-0.012, 4.92]
    man_e            0.012     [0.012, -2.79]      fum_e            0         [-0.012, 4.73]
    cit_e            0.012     [0.012, -2.97]      gal_e            0         [-0.012, 3.82]
    ala_B_e          0.012     [0.012, -3.24]      abt_e            0         [0, 3.44]
    pro_L_e          0.012     [0.012, -3.24]      xylt_e           0         [-0.012, 3.43]
    fru_e            0.012     [0.012, -4.26]      akg_e            0         [-0.012, 3.39]
    ppa_e            0.012     [0.012, -8.04]      glc_D_e          0         [-1, 3.27]
    mal_L_e          0.012     [0.012, -4.73]      3aib_D_e         0         [0, 3.25]
    ile_L_e          0.011     [0.011, 0.012]      3aib_e           0         [0, 3.25]
    asn_L_e          0.0107    [0.012, -1.62]      4abut_e          0         [0, 3.25]
    phe_L_e          0.00996   [-0.00472, 0.012]   4abut_n          0         [0, 3.25]
    4mop_e           0.00894   [0.00894, 0.012]    5oxpro_e         0         [0, 3.25]
    glu_L_e          0.00821   [0.012, -3.24]      asp_L_e          0         [-0.012, 3.24]
    ps_hs_e          0.00817   [0.012, -0.822]     gly_e            0         [-0.03, 3.22]
    ser_L_e          0.00712   [0.012, -3.24]      ala_L_e          0         [-0.04, 3.21]
    met_L_e          0.00587   [0.00587, 0.012]    2mcit_e          0         [0, 3.2]
    his_L_e          0.00485   [-0.00715, 0.012]   HC00342_e        0         [0, 2.98]
    pe_hs_r          0.00213   [0.012, -0.82]      HC01444_e        0         [0, 2.75]
    biomass_other_c  0.00207   [0.00207, 0.00207]  glyc3p_e         0         [-0.012, 2.72]
    orn_e            0.00179   [0.012, -1.62]      oxa_e            0         [0, 2.49]
    strch1_e         0.00164   [0, 0.012]          4abutn_e         0         [0, 2.31]
    inost_e          0.000895  [0.012, -4.15]      ttdca_e          0         [-0.012, 2.16]
    sph1p_e          0.000671  [0.012, -1.65]      HC00229_e        0         [0, 2.14]
    pglyc_hs_e       0.000559  [0.012, -0.813]     acmana_e         0         [0, 2.13]
    trp_L_e          0.000511  [0.000511, 0.012]   HC01441_e        0         [0, 2.05]
    dttp_n           0.000502  [0.012, -0.898]     ethamp_r         0         [0, 1.98]
    h2o_e            0         [-9.79, 10]         1glyc_hs_e       0         [0, 1.89]
    chsterol_e       0         [-0.16, 1.82]       Rtotal_e         0         [0, 1.89]
    pydam_e          0         [0, 1]              xolest2_hs_e     0         [0, 1.82]
    btn_e            0         [0, 0.2]            urate_e          0         [0, 1.75]
    gthrd_e          0         [-0.0883, 0.09]     hxan_e           0         [-0.012, 1.74]
    3hpvs_e          0         [0, 0.012]          din_e            0         [-0.012, 1.72]
    4hpro_LT_m       0         [0, 0.012]          ins_e            0         [0, 1.69]
    arab_L_e         0         [0, 0.012]          urea_e           0         [0, 1.64]
    carn_e           0         [0, 0.012]          cbasp_e          0         [0, 1.63]
    caro_e           0         [0, 0.012]          gluala_e         0         [-0.012, 1.62]
    csa_e            0         [0, 0.012]          mag_hs_e         0         [-0.012, 1.61]
    csn_e            0         [0, 0.012]          hdca_e           0         [-0.3, 1.59]
    dhdascb_e        0         [0, 0.012]          lcts_e           0         [-0.5, 1.55]
    etoh_e           0         [0, 0.012]          thym_e           0         [0, 1.52]
    fol_e            0         [0, 0.012]          duri_e           0         [0, 1.48]
    gam_e            0         [0, 0.012]          orot_e           0         [0, 1.48]
    glygn2_e         0         [0, 0.012]          thymd_e          0         [-0.012, 1.44]
    n2m2nmasn_e      0         [0, 0.012]          h2o2_e           0         [0, 1.43]
    octa_e           0         [0, 0.012]          sphs1p_e         0         [0, 1.43]
    pnto_R_e         0         [0, 0.012]          ha_pre1_e        0         [0, 1.4]
    ptvstlac_e       0         [0, 0.012]          ade_e            0         [-0.012, 1.39]
    pydxn_e          0         [0, 0.012]          adn_e            0         [-0.012, 1.39]
    rib_D_e          0         [0, 0.012]          dad_2_e          0         [-0.012, 1.39]
    sucr_e           0         [0, 0.012]          dgsn_e           0         [-0.012, 1.39]
    thm_e            0         [0, 0.012]          amp_e            0         [-0.012, 1.38]
    tre_e            0         [0, 0.012]          ppi_e            0         [0, 1.36]
    ribflv_e         0         [-0.012, 0.012]     camp_e           0         [0, 1.35]
    crvnc_e          0         [0.012, -0.31]      35cgmp_e         0         [0, 1.31]
                                                   udp_e            0         [0, 1.27]
                                                   ump_e            0         [0, 1.27]
                                                   nrvnc_e          0         [0, 1.25]
                                                   elaid_e          0         [0, 1.23]
                                                   ach_e            0         [0, 1.19]
                                                   whhdca_e         0         [0, 1.19]
                                                   lnlc_e           0         [-0.06, 1.17]
                                                   whddca_e         0         [0, 1.15]
                                                   fald_e           0         [-0.012, 1.14]
                                                   meoh_e           0         [-0.012, 1.14]
                                                   citr_L_c         0         [0, 1.09]
                                                   citr_L_e         0         [0, 1.09]
                                                   creat_e          0         [-0.012, 1.08]
                                                   4pyrdx_e         0         [0, 1.04]
                                                   pydx5p_e         0         [-0.012, 1.02]
                                                   pydx_e           0         [-0.012, 1.02]
                                                   dctp_m           0         [0, 0.91]
                                                   dgtp_m           0         [0, 0.91]
                                                   dttp_m           0         [0, 0.91]
                                                   prpp_e           0         [0, 0.91]
                                                   spmd_e           0         [-0.012, 0.899]
                                                   datp_n           0         [-0.012, 0.898]
                                                   crmp_hs_e        0         [0, 0.885]
                                                   CE1940_e         0         [0, 0.833]
                                                   pe_hs_e          0         [-0.012, 0.82]
                                                   lpchol_hs_e      0         [-0.012, 0.804]
                                                   spc_hs_e         0         [0, 0.802]
                                                   no_e             0         [0, 0.772]
                                                   lnlncg_e         0         [-0.012, 0.704]
                                                   ha_e             0         [-0.012, 0.693]
                                                   pchol_hs_e       0         [0, 0.663]
                                                   tag_hs_e         0         [-0.012, 0.59]
                                                   sprm_c           0         [0, 0.562]
                                                   sprm_e           0         [0, 0.562]
                                                   dlnlcg_e         0         [0, 0.494]
                                                   eicostet_e       0         [0, 0.483]
                                                   gbside_hs_e      0         [0, 0.476]
                                                   fuc14galacgl...  0         [0, 0.466]
                                                   fucgal14acgl...  0         [0, 0.466]
                                                   fucfuc12gal1...  0         [0, 0.435]
                                                   fucfucgalacg...  0         [0, 0.435]
                                                   pheme_e          0         [0, 0.43]
                                                   fuc_L_e          0         [0, 0.422]
                                                   fucacngal14a...  0         [0, 0.42]
                                                   fucgalgbside...  0         [0, 0.416]
                                                   galfuc12gal1...  0         [0, 0.416]
                                                   galgalgalthc...  0         [0, 0.404]
                                                   tmndnc_e         0         [0, 0.397]
                                                   fucgalfucgal...  0         [0, 0.39]
                                                   clpnd_e          0         [-0.012, 0.384]
                                                   fucacgalfucg...  0         [0, 0.377]
                                                   arachd_e         0         [-0.012, 0.371]
                                                   fuc13galacgl...  0         [0, 0.361]
                                                   gt1a_hs_e        0         [0, 0.355]
                                                   galacglcgalg...  0         [0, 0.347]
                                                   fucfuc132gal...  0         [0, 0.343]
                                                   galfucgalacg...  0         [0, 0.33]
                                                   fucfucfucgal...  0         [0, 0.325]
                                                   leuktrB4_e       0         [0, 0.32]
                                                   adrn_e           0         [0, 0.311]
                                                   5HPET_c          0         [0, 0.308]
                                                   5HPET_r          0         [0, 0.308]
                                                   leuktrA4_e       0         [-0.012, 0.308]
                                                   acgalfucgala...  0         [0, 0.302]
                                                   prostgf2_e       0         [0, 0.271]
                                                   prostgd2_e       0         [0, 0.27]
                                                   fucfucfucgal...  0         [0, 0.262]
                                                   dcsptn1_e        0         [0, 0.26]
                                                   prostgi2_e       0         [0, 0.258]
                                                   txa2_e           0         [0, 0.258]
                                                   prostgh2_e       0         [-0.012, 0.258]
                                                   acgalfucgala...  0         [0, 0.245]
                                                   galgalfucfuc...  0         [0, 0.234]
                                                   digalsgalsid...  0         [0, 0.228]
                                                   sl_L_e           0         [0, 0.228]
                                                   taur_c           0         [0, 0.228]
                                                   taur_e           0         [0, 0.228]
                                                   bilirub_e        0         [0, 0.228]
                                                   co_e             0         [0, 0.228]
                                                   bilglcur_e       0         [0, 0.22]
                                                   bildglcur_e      0         [0, 0.212]
                                                   leuktrC4_e       0         [0, 0.202]
                                                   leuktrF4_e       0         [0, 0.202]
                                                   btn_c            0         [0, 0.2]
                                                   leuktrD4_e       0         [-0.012, 0.19]
                                                   leuktrE4_e       0         [-0.012, 0.19]
                                                   vitd3_e          0         [0, 0.159]
                                                   xoltri24_e       0         [0, 0.13]
                                                   xoltri25_e       0         [0, 0.13]
                                                   xoltri27_e       0         [0, 0.13]
                                                   1a25dhvitd3_n    0         [0, 0.13]
                                                   2425dhvitd3_e    0         [0, 0.13]
                                                   4mptnl_e         0         [0, 0.13]
                                                   aprgstrn_e       0         [0, 0.13]
                                                   prgstrn_e        0         [0, 0.13]
                                                   cholate_e        0         [0, 0.126]
                                                   C02528_e         0         [0, 0.121]
                                                   dgchol_e         0         [0, 0.121]
                                                   tchola_e         0         [0, 0.119]
                                                   tdchola_e        0         [0, 0.115]
                                                   5adtststerone_e  0         [0, 0.115]
                                                   andrstrn_e       0         [0, 0.115]
                                                   tststerone_e     0         [0, 0.115]
                                                   gchola_e         0         [-0.012, 0.114]
                                                   5adtststeron...  0         [0, 0.112]
                                                   andrstrnglc_e    0         [0, 0.112]
                                                   tststeroneglc_e  0         [0, 0.112]
                                                   crtstrn_e        0         [0, 0.11]
                                                   dheas_e          0         [0, 0.109]
                                                   6htststerone_e   0         [0, 0.106]
                                                   crtsl_e          0         [0, 0.102]
                                                   aldstrn_e        0         [0, 0.102]
                                                   estradiol_e      0         [0, 0.0987]
                                                   estradiolglc_e   0         [0, 0.0972]
                                                   estroneglc_e     0         [0, 0.0969]
                                                   estrones_e       0         [0, 0.0946]
                                                   hestratriol_e    0         [0, 0.0923]
                                                   13_cis_oretn_n   0         [0, 0.06]
                                                   13_cis_retn_n    0         [0, 0.06]
                                                   13_cis_retng...  0         [0, 0.06]
                                                   hretn_n          0         [0, 0.06]
                                                   oretn_n          0         [0, 0.06]
                                                   retnglc_e        0         [0, 0.06]
                                                   11_cis_retfa_e   0         [0, 0.048]
                                                   9_cis_retfa_e    0         [0, 0.048]
                                                   retinol_9_cis_e  0         [0, 0.048]
                                                   retinol_cis_...  0         [0, 0.048]
                                                   retn_e           0         [-0.012, 0.048]
                                                   malt_e           0         [-0.012, 0.036]
                                                   retfa_e          0         [-0.012, 0.036]
                                                   retinol_e        0         [-0.012, 0.036]
                                                   34dhoxpeg_e      0         [0, 0.0319]
                                                   adrnl_e          0         [0, 0.0319]
                                                   dopasf_e         0         [0, 0.0319]
                                                   mepi_e           0         [0, 0.0319]
                                                   nrpphr_e         0         [0, 0.0319]
                                                   nrpphrsf_e       0         [0, 0.0319]
                                                   q10_e            0         [0, 0.0319]
                                                   am19cs_e         0         [0, 0.024]
                                                   am4n9cs_e        0         [0, 0.024]
                                                   fad_e            0         [0, 0.024]
                                                   1mncam_e         0         [0, 0.0235]
                                                   ncam_c           0         [0, 0.0235]
                                                   ncam_e           0         [0, 0.0235]
                                                   34dhphe_e        0         [0, 0.0199]
                                                   4hphac_e         0         [0, 0.0199]
                                                   tymsf_e          0         [0, 0.0199]
                                                   dopa_e           0         [-0.012, 0.0199]
                                                   q10h2_e          0         [-0.012, 0.0199]
                                                   3mlda_e          0         [0, 0.0191]
                                                   hista_e          0         [0, 0.0191]
                                                   pheacgln_e       0         [0, 0.0167]
                                                   3mob_e           0         [-0.012, 0.0165]
                                                   3ddcrn_e         0         [0, 0.012]
                                                   3deccrn_e        0         [0, 0.012]
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
                                                   4hpro_LT_e       0         [0, 0.012]
                                                   5mthf_e          0         [0, 0.012]
                                                   Asn_X_Ser_Thr_l  0         [0, 0.012]
                                                   Tyr_ggn_e        0         [0, 0.012]
                                                   am1a4ncs_e       0         [0, 0.012]
                                                   am1accs_e        0         [0, 0.012]
                                                   am1acs_e         0         [0, 0.012]
                                                   am1alcs_e        0         [0, 0.012]
                                                   am1c4n9cs_e      0         [0, 0.012]
                                                   am1c9cs_e        0         [0, 0.012]
                                                   am1ccs_e         0         [0, 0.012]
                                                   am1cglc_e        0         [0, 0.012]
                                                   am1csa_e         0         [0, 0.012]
                                                   am4ncs_e         0         [0, 0.012]
                                                   ascb_L_e         0         [0, 0.012]
                                                   c10crn_e         0         [0, 0.012]
                                                   c4crn_e          0         [0, 0.012]
                                                   c51crn_e         0         [0, 0.012]
                                                   c5dc_e           0         [0, 0.012]
                                                   c6crn_e          0         [0, 0.012]
                                                   c8crn_e          0         [0, 0.012]
                                                   csasulp_e        0         [0, 0.012]
                                                   ddeccrn_e        0         [0, 0.012]
                                                   fol_c            0         [0, 0.012]
                                                   glygn4_e         0         [0, 0.012]
                                                   glygn5_e         0         [0, 0.012]
                                                   n5m2masn_g       0         [0, 0.012]
                                                   octdececoa_c     0         [0, 0.012]
                                                   pmtcoa_r         0         [0, 0.012]
                                                   pnto_R_c         0         [0, 0.012]
                                                   ptvst_e          0         [0, 0.012]
                                                   ptvstm13_e       0         [0, 0.012]
                                                   ptvstm3_e        0         [0, 0.012]
                                                   rbt_e            0         [0, 0.012]
                                                   thf_e            0         [0, 0.012]
                                                   thmmp_e          0         [0, 0.012]
                                                   thmtp_e          0         [0, 0.012]
                                                   ttdcrn_e         0         [0, 0.012]
                                                   am9csa_e         0         [-0.012, 0.012]
                                                   fmn_e            0         [-0.012, 0.012]
                                                   strch2_e         0         [-0.012, 0.012]
                                                   5htrp_e          0         [0, 0.0115]
                                                   C02470_e         0         [0, 0.0115]
                                                   anth_c           0         [0, 0.0115]
                                                   anth_e           0         [0, 0.0115]
                                                   srtn_e           0         [0, 0.0115]
                                                   5mta_e           0         [0, 0.00613]
                                                   ahcys_e          0         [0, 0.00613]
                                                   ivcrn_e          0         [0, 0.00306]
                                                   3mop_e           0         [0, 0.00102]
                                                   lac_D_e          0         [-0.001, 0.001]
                                                   nac_e            0         [0.0115, -0.012]



```python
%%time

from corda import CORDA


opt_MeanNormalBiopsy = CORDA(model=Recon2, confidence=conf_MeanNormalBiopsy, n=5, met_prod=metas,  penalty_factor=1000) 
opt_MeanNormalBiopsy.build()
print(opt_MeanNormalBiopsy)
model_MeanNormalBiopsy=opt_MeanNormalBiopsy.cobra_model(name="MeanNormalBiopsy")
print(model_MeanNormalBiopsy.optimize())
```

    build status: reconstruction complete
    Inc. reactions: 2336/7891
     - unclear: 167/640
     - exclude: 250/2368
     - low and medium: 722/3508
     - high: 1197/1375
    
    <Solution 0.022 at 0x7f93aa861630>
    CPU times: user 27min 44s, sys: 780 ms, total: 27min 45s
    Wall time: 27min 45s



```python
cp=model_MeanNormalBiopsy.copy()
cp.optimize()
cp.summary(fva=True)
```

    IN FLUXES                                        OUT FLUXES                                      OBJECTIVES
    -----------------------------------------------  ----------------------------------------------  ----------------------
    id                   Flux  Range                 id                   Flux  Range                biomass_reac...  0.022
    ---------------  --------  --------------------  ---------------  --------  -------------------
    o2_e             0.183     [0.138, 1]            h_e              0.294     [-0.96, 3.48]
    glc_D_e          0.0924    [0.00605, 1]          glyc_S_e         0.0916    [0, 1.08]
    h2o_e            0.0852    [-1.73, 3.13]         xmp_e            0.0814    [-0.012, 0.832]
    atp_e            0.0834    [-0.0325, 0.8]        ala_L_e          0.0805    [-0.04, 0.359]
    gthrd_e          0.0378    [-0.023, 0.09]        nh4_e            0.06      [-0.012, 1.29]
    lys_L_e          0.013     [0.013, 0.013]        cit_e            0.0541    [-0.012, 0.311]
    leu_L_e          0.012     [0.012, 0.012]        pi_e             0.0472    [0, 1.76]
    ps_hs_e          0.012     [0.00391, 0.012]      ppi_e            0.0472    [0, 0.586]
    4hpro_LT_m       0.012     [0, 0.012]            cgly_e           0.0258    [-0.012, 0.101]
    crn_e            0.012     [0, 0.012]            oxa_e            0.024     [0, 0.511]
    leuktrA4_e       0.012     [0, 0.012]            leuktrD4_e       0.012     [0, 0.012]
    o2s_e            0.012     [0, 0.012]            retn_e           0.012     [0, 0.012]
    retinol_e        0.012     [0, 0.012]            pyr_e            0.0116    [-0.012, 0.138]
    bhb_e            0.012     [-0.00931, 0.012]     utp_e            0.0108    [-0.02, 0.0458]
    glyc3p_e         0.012     [-0.0119, 0.012]      tmndnc_e         0.008     [0, 0.096]
    uri_e            0.012     [0.012, -0.0538]      lac_L_e          0.008     [0.06, -0.09]
    glu_L_e          0.012     [0.012, -0.101]       dag_hs_e         0.00771   [-0.000385, 0.0235]
    akg_e            0.012     [0.012, -0.138]       3bcrn_e          0.00669   [0, 0.012]
    ser_L_e          0.012     [0.012, -0.387]       3ivcrn_e         0.00531   [0, 0.012]
    xylt_e           0.012     [0.012, -0.625]       acac_e           0.004     [0.00931, -0.012]
    gly_e            0.0119    [0.03, -0.369]        Rtotal_e         0.00378   [0, 0.25]
    pro_L_e          0.00907   [0.00907, 0.00907]    lac_D_e          0.001     [0, 0.001]
    crvnc_e          0.008     [0, 0.012]            for_e            0.000449  [-0.012, 0.04]
    gln_L_e          0.00796   [0.09, -0.203]        co2_e            0         [0, 1.35]
    arg_L_e          0.0079    [0.0079, 0.0079]      amp_e            0         [-0.012, 0.833]
    asp_L_e          0.00776   [0, 0.012]            abt_e            0         [0, 0.637]
    val_L_e          0.00776   [0.00776, 0.00776]    fe3_e            0         [-0.012, 0.545]
    thr_L_e          0.00688   [0.00688, 0.00688]    dgtp_m           0         [0, 0.523]
    ile_L_e          0.00629   [0.00629, 0.012]      man_e            0         [-0.012, 0.504]
    asn_L_e          0.00615   [0.00615, 0.00615]    hdcea_e          0         [0, 0.428]
    phe_L_e          0.00571   [0.00571, 0.012]      vacc_e           0         [0, 0.424]
    lpchol_hs_e      0.00378   [0, 0.0119]           aicar_e          0         [0, 0.392]
    tyr_L_e          0.00351   [0.00351, 0.012]      ade_e            0         [0, 0.356]
    met_L_e          0.00337   [0.00337, 0.012]      prpp_e           0         [0, 0.356]
    his_L_e          0.00278   [0, 0.012]            elaid_e          0         [0, 0.356]
    pe_hs_e          0.00122   [0.012, -0.0189]      5oxpro_e         0         [0, 0.31]
    biomass_other_c  0.00119   [0.00119, 0.00119]    fum_e            0         [-0.012, 0.289]
    cys_L_e          0.00102   [0.00102, 0.012]      mal_L_e          0         [-0.012, 0.289]
    cytd_e           0.000859  [0.012, -0.0538]      ocdcea_e         0         [0, 0.27]
    inost_e          0.000513  [0.000513, 0.000513]  gal_e            0         [0, 0.204]
    hdca_e           0.000385  [-0.0116, 0.3]        ac_e             0         [0, 0.15]
    sph1p_e          0.000385  [0.000385, 0.012]     nrvnc_e          0         [0, 0.142]
    pglyc_hs_e       0.000321  [0.000321, 0.012]     clpnd_e          0         [0, 0.096]
    trp_L_e          0.000293  [0.000293, 0.000293]  eicostet_e       0         [0, 0.084]
    datp_n           0.00029   [0.012, -0.457]       dlnlcg_e         0         [0, 0.072]
    dttp_n           0.000288  [0.000288, -0.0117]   lnlncg_e         0         [0, 0.072]
    dgtp_n           0.000218  [0.012, -0.529]       lnlnca_e         0         [-0.012, 0.072]
    dctp_n           0.000208  [0.012, -0.0538]      dctp_m           0         [0, 0.0658]
    fe2_e            0         [-0.012, 0.545]       gthox_e          0         [0, 0.0565]
    glyc_e           0         [-0.012, 0.21]        dcmp_e           0         [-0.012, 0.0538]
    ocdca_e          0         [0, 0.07]             pheme_e          0         [0, 0.0325]
    lnlc_e           0         [-0.012, 0.06]        bilglcur_e       0         [0, 0.0311]
    3hpvs_e          0         [0, 0.012]            co_e             0         [0, 0.0311]
    3mob_e           0         [0, 0.012]            bildglcur_e      0         [0, 0.0307]
    4mop_e           0         [0, 0.012]            fald_e           0         [-0.012, 0.024]
    arab_L_e         0         [0, 0.012]            meoh_e           0         [-0.012, 0.024]
    arachd_e         0         [0, 0.012]            pe_hs_r          0         [-0.012, 0.0189]
    carn_e           0         [0, 0.012]            cholate_e        0         [0, 0.0156]
    csa_e            0         [0, 0.012]            3mlda_e          0         [0, 0.0126]
    dopa_e           0         [0, 0.012]            34dhoxpeg_e      0         [0, 0.012]
    etoh_e           0         [0, 0.012]            3ddcrn_e         0         [0, 0.012]
    fol_e            0         [0, 0.012]            3deccrn_e        0         [0, 0.012]
    gam_e            0         [0, 0.012]            3hdececrn_e      0         [0, 0.012]
    gluala_e         0         [0, 0.012]            3hexdcrn_e       0         [0, 0.012]
    ha_e             0         [0, 0.012]            3hpvstet_e       0         [0, 0.012]
    mag_hs_e         0         [0, 0.012]            3octdec2crn_e    0         [0, 0.012]
    n2m2nmasn_e      0         [0, 0.012]            3octdeccrn_e     0         [0, 0.012]
    octa_e           0         [0, 0.012]            3octdece1crn_e   0         [0, 0.012]
    ptvstlac_e       0         [0, 0.012]            3tdcrn_e         0         [0, 0.012]
    thymd_e          0         [0, 0.012]            3tetd7ecoacrn_e  0         [0, 0.012]
    gchola_e         0         [-0.00355, 0.012]     3thexddcoacrn_e  0         [0, 0.012]
    HC01609_e        0         [-0.012, 0.012]       3ttetddcoacrn_e  0         [0, 0.012]
    pnto_R_e         0         [0, 0.011]            adrn_e           0         [0, 0.012]
                                                     am19cs_e         0         [0, 0.012]
                                                     am1csa_e         0         [0, 0.012]
                                                     am9csa_e         0         [0, 0.012]
                                                     c6crn_e          0         [0, 0.012]
                                                     c8crn_e          0         [0, 0.012]
                                                     fol_c            0         [0, 0.012]
                                                     ivcrn_e          0         [0, 0.012]
                                                     leuktrB4_e       0         [0, 0.012]
                                                     n5m2masn_g       0         [0, 0.012]
                                                     ptvstm3_e        0         [0, 0.012]
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
                                                     HC01610_e        0         [0.012, -0.012]



```python
%%time

from corda import CORDA


opt_MeanCancerBiopsy = CORDA(model=Recon2, confidence=conf_MeanCancerBiopsy, n=5, met_prod=metas,  penalty_factor=1000) 
opt_MeanCancerBiopsy.build()
print(opt_MeanCancerBiopsy)
model_MeanCancerBiopsy=opt_MeanCancerBiopsy.cobra_model(name="MeanCancerBiopsy")
print(model_MeanCancerBiopsy.optimize())
```

    build status: reconstruction complete
    Inc. reactions: 2714/7891
     - unclear: 281/1381
     - exclude: 153/1275
     - low and medium: 732/3500
     - high: 1548/1735
    
    <Solution 0.033 at 0x7f9350be3080>
    CPU times: user 30min 26s, sys: 1.83 s, total: 30min 28s
    Wall time: 30min 27s



```python
cp=model_MeanCancerBiopsy.copy()
cp.optimize()
cp.summary(fva=True)
```

    IN FLUXES                                        OUT FLUXES                                    OBJECTIVES
    -----------------------------------------------  --------------------------------------------  -----------------------
    id                   Flux  Range                 id                   Flux  Range              biomass_reac...  0.0334
    ---------------  --------  --------------------  ---------------  --------  -----------------
    glc_D_e          0.381     [0.22, 1]             h_e              0.724     [-0.372, 3.23]
    o2_e             0.0592    [-0.0058, 1]          lac_L_e          0.623     [-0.09, 2.69]
    val_L_e          0.03      [-0.000222, 0.03]     pyr_e            0.0723    [-0.012, 1.37]
    gln_L_e          0.0256    [0.09, -0.142]        h2o_e            0.0458    [-0.37, 1.69]
    Lcystin_c        0.024     [-0.00522, 0.024]     cys_L_e          0.0409    [-0.012, 0.0464]
    pi_e             0.0218    [0.04, -0.0836]       glyc_e           0.0313    [-0.21, 0.299]
    lys_L_e          0.0198    [0.0198, 0.03]        ala_L_e          0.0149    [-0.0169, 0.294]
    gly_e            0.018     [0.03, -0.269]        nh4_e            0.0147    [-0.012, 0.451]
    arg_L_e          0.012     [0.012, 0.012]        co2_e            0.0145    [0, 1.72]
    thr_L_e          0.012     [0.0104, 0.012]       34dhoxpeg_e      0.012     [0, 0.012]
    leu_L_e          0.012     [0.00622, 0.012]      3mob_e           0.00849   [-0.012, 0.0182]
    dopa_e           0.012     [0, 0.012]            dag_hs_e         0.00571   [-0.012, 0.183]
    ps_hs_e          0.012     [0.012, -0.112]       so4_e            0.00557   [0, 0.0584]
    pro_L_e          0.012     [0.012, -0.287]       2hb_e            0.00156   [0.00156, -0.012]
    ser_L_e          0.012     [0.012, -0.287]       lac_D_e          0.001     [0, 0.001]
    ile_L_e          0.00956   [0.00956, 0.012]      for_e            0.000681  [-0.012, 0.0992]
    asn_L_e          0.00933   [0.00933, 0.012]      fe3_e            0         [-0.012, 0.545]
    phe_L_e          0.00867   [0.00867, 0.012]      ac_e             0         [-0.012, 0.479]
    4mop_e           0.00622   [0.00622, 0.012]      vacc_e           0         [0, 0.452]
    chol_e           0.00552   [-0.0116, 0.012]      oxa_e            0         [0, 0.449]
    tyr_L_e          0.00533   [0.00533, 0.012]      elaid_e          0         [0, 0.449]
    met_L_e          0.00511   [0.00511, 0.012]      no_e             0         [0, 0.411]
    his_L_e          0.00422   [0.00422, 0.012]      gal_e            0         [0, 0.39]
    cytd_e           0.00309   [0.012, -0.0517]      cit_e            0         [-0.012, 0.348]
    hdca_e           0.00263   [-0.0291, 0.3]        ocdcea_e         0         [0, 0.34]
    pe_hs_e          0.00208   [0.012, -0.112]       akg_e            0         [-0.012, 0.331]
    biomass_other_c  0.0018    [0.0018, 0.0018]      hdcea_e          0         [0, 0.327]
    adn_e            0.00179   [0.00868, -0.021]     xylt_e           0         [-0.012, 0.319]
    gsn_e            0.00121   [-0.0108, 0.012]      glyc_S_e         0         [0, 0.299]
    inost_e          0.000779  [0.000779, 0.000779]  acmana_e         0         [0, 0.287]
    sph1p_e          0.000584  [0.012, -0.112]       ha_pre1_e        0         [0, 0.187]
    pglyc_hs_e       0.000487  [0.000487, 0.000487]  nrvnc_e          0         [0, 0.17]
    trp_L_e          0.000444  [0.000444, 0.012]     ethamp_r         0         [0, 0.124]
    datp_n           0.00044   [0.00044, 0.00044]    pe_hs_r          0         [-0.012, 0.112]
    dttp_n           0.000437  [0.000437, -0.0408]   ha_e             0         [-0.012, 0.0868]
    dgtp_n           0.000331  [0.000331, 0.000331]  thymd_e          0         [-0.012, 0.0757]
    dctp_n           0.000315  [0.000315, 0.000315]  urate_e          0         [0, 0.0721]
    fe2_e            0         [-0.012, 0.545]       ppa_e            0         [-0.012, 0.0617]
    ocdca_e          0         [0, 0.07]             clpnd_e          0         [0, 0.06]
    lnlc_e           0         [0, 0.06]             dlnlcg_e         0         [0, 0.06]
    utp_e            0         [0, 0.02]             eicostet_e       0         [0, 0.06]
    3hpvs_e          0         [0, 0.012]            lnlnca_e         0         [0, 0.06]
    4hpro_LT_m       0         [0, 0.012]            lnlncg_e         0         [0, 0.06]
    amp_e            0         [0, 0.012]            taur_c           0         [0, 0.0584]
    arab_L_e         0         [0, 0.012]            taur_e           0         [0, 0.0584]
    arachd_e         0         [0, 0.012]            tchola_e         0         [0, 0.0537]
    asp_L_e          0         [0, 0.012]            uri_e            0         [-0.012, 0.0517]
    crn_e            0         [0, 0.012]            chsterol_e       0         [0, 0.0514]
    crvnc_e          0         [0, 0.012]            4mptnl_e         0         [0, 0.0501]
    dcyt_e           0         [0, 0.012]            aprgstrn_e       0         [0, 0.0501]
    din_e            0         [0, 0.012]            5adtststerone_e  0         [0, 0.0487]
    etoh_e           0         [0, 0.012]            andrstrn_e       0         [0, 0.0487]
    fmn_e            0         [0, 0.012]            fuc13galacgl...  0         [0, 0.048]
    gam_e            0         [0, 0.012]            fuc14galacgl...  0         [0, 0.048]
    gluala_e         0         [0, 0.012]            gchola_e         0         [-0.012, 0.0457]
    leuktrA4_e       0         [0, 0.012]            galfucgalacg...  0         [0, 0.0445]
    n2m2nmasn_e      0         [0, 0.012]            5adtststeron...  0         [0, 0.0432]
    o2s_e            0         [0, 0.012]            andrstrnglc_e    0         [0, 0.0432]
    orn_e            0         [0, 0.012]            tststeroneglc_e  0         [0, 0.043]
    ptrc_e           0         [0, 0.012]            man_e            0         [-0.012, 0.036]
    ptvstlac_e       0         [0, 0.012]            fuc_L_e          0         [0, 0.0338]
    retinol_e        0         [0, 0.012]            13_cis_retng...  0         [0, 0.024]
    spmd_e           0         [0, 0.012]            gthox_e          0         [0, 0.024]
    HC01609_e        0         [-0.012, 0.012]       retnglc_e        0         [0, 0.024]
    HC01610_e        0         [-0.012, 0.012]       fald_e           0         [-0.012, 0.024]
                                                     meoh_e           0         [-0.012, 0.024]
                                                     1mncam_e         0         [0, 0.0236]
                                                     glyb_e           0         [0, 0.0235]
                                                     Rtotal_e         0         [0, 0.0232]
                                                     lpchol_hs_e      0         [0, 0.0232]
                                                     pchol_hs_e       0         [0, 0.0232]
                                                     ppi_e            0         [0, 0.02]
                                                     udp_e            0         [0, 0.02]
                                                     ump_e            0         [0, 0.02]
                                                     fucfucfucgal...  0         [0, 0.016]
                                                     3bcrn_e          0         [0, 0.012]
                                                     3ddcrn_e         0         [0, 0.012]
                                                     3deccrn_e        0         [0, 0.012]
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
                                                     leuktrD4_e       0         [0, 0.012]
                                                     n5m2masn_g       0         [0, 0.012]
                                                     prpp_e           0         [0, 0.012]
                                                     ptvst_e          0         [0, 0.012]
                                                     ribflv_e         0         [0, 0.012]
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
                                                     nac_e            0         [0.0116, -0.012]
                                                     gthrd_e          0         [0.0584, -0.06]



```python
model_MeanCancerBiopsy.metabolites.atp_c.summary(fva=0.95)

```

    PRODUCING REACTIONS -- ATP(4-) (atp_c)
    --------------------------------------
    %        FLUX  RANGE      RXN ID      REACTION
    ---  --------  ---------  ----------  ----------------------------------------
    51%  0.772     [0, -1.6]  PGK         3pg_c + atp_c <=> 13dpg_c + adp_c
    48%  0.723     [0, -1.6]  PYK         adp_c + h_c + pep_c --> atp_c + pyr_c
    2%   0.0269    [0, -1.6]  DAGK_hs     atp_c + dag_hs_c <=> adp_c + h_c + pa...
    0%   0         [0, -1.6]  r0377       atp_c + dcyt_c --> adp_c + dcmp_c + h_c
    
    CONSUMING REACTIONS -- ATP(4-) (atp_c)
    --------------------------------------
    %        FLUX  RANGE      RXN ID      REACTION
    ---  --------  ---------  ----------  ----------------------------------------
    45%  0.69      [0, 1.6]   biomass...  0.716189801699717 ala_L_c + 0.5088668...
    26%  0.389     [0, 1.6]   HEX1        atp_c + glc_D_c --> adp_c + g6p_c + h_c
    25%  0.38      [0, 1.6]   PFK         atp_c + f6p_c --> adp_c + fdp_c + h_c
    1%   0.022     [0, 1.6]   ETHAK       atp_c + etha_c --> adp_c + ethamp_c +...
    0%   0.00442   [0, 1.6]   ADK1        amp_c + atp_c <=> 2.0 adp_c
    0%   0.00409   [0, 1.6]   DPMVDc      5dpmev_c + atp_c --> adp_c + co2_c + ...
    0%   0.00409   [0, 1.6]   MEVK1c      atp_c + mev_R_c --> 5pmev_c + adp_c +...
    0%   0.00409   [0, 1.6]   PMEVKc      5pmev_c + atp_c --> 5dpmev_c + adp_c
    0%   0.004     [0, 1.6]   CHOLK       atp_c + chol_c --> adp_c + cholp_c + h_c
    0%   0.00357   [0, 1.6]   UMPK        atp_c + ump_c <=> adp_c + udp_c
    0%   0.00247   [0, 1.6]   CYTK1       atp_c + cmp_c <=> adp_c + cdp_c
    0%   0.00247   [0, 1.6]   CYTK2       atp_c + dcmp_c <=> adp_c + dcdp_c
    0%   0.00204   [0, 1.6]   FACOAL204   arachd_c + atp_c + coa_c --> amp_c + ...
    0%   0.00179   [0, 1.6]   ADNK1       adn_c + atp_c --> adp_c + amp_c + h_c
    0%   0.00179   [0, 1.6]   biomass...  0.925862068965503 atp_c + 0.673034482...
    0%   0.00178   [0, 1.6]   URIK1       atp_c + uri_c --> adp_c + h_c + ump_c
    0%   0.0013    [0, 1.6]   CYTDK1      atp_c + cytd_c --> adp_c + cmp_c + h_c
    0%   0.00121   [0, 1.6]   GK1         atp_c + gmp_c <=> adp_c + gdp_c
    0%   0.00121   [0, 1.6]   NICRNS      atp_c + nicrns_c --> adp_c + h_c + ni...
    0%   0.000584  [0, 1.6]   FACOAL160i  atp_c + coa_c + hdca_c --> amp_c + pm...
    0%   0         [0, 1.6]   1MNCAMti    1mncam_c + atp_c + h2o_c --> 1mncam_e...
    0%   0         [0, 1.6]   3HPVSCOAhc  3hpvs_c + atp_c + coa_c --> 3hpvscoa_...
    0%   0         [0, 1.6]   5ADTSTS...  5adtststeroneglc_c + atp_c + h2o_c --...
    0%   0         [0, 1.6]   AACOAT      acac_c + atp_c + coa_c --> aacoa_c + ...
    0%   0         [0, 1.6]   ACCOAC      accoa_c + atp_c + hco3_c --> adp_c + ...
    0%   0         [0, 1.6]   ACGAMK      acgam_c + atp_c --> acgam6p_c + adp_c...
    0%   0         [0, 1.6]   ACITL       atp_c + cit_c + coa_c --> accoa_c + a...
    0%   0         [0, 1.6]   ACS         ac_c + atp_c + coa_c --> accoa_c + am...
    0%   0         [0, 1.6]   ACS2        atp_c + coa_c + ppa_c --> amp_c + ppc...
    0%   0         [0, 1.6]   ADSK        aps_c + atp_c --> adp_c + h_c + paps_c
    0%   0         [0, 1.6]   AIRCr_P...  air_c + asp_L_c + atp_c + co2_c <=> 2...
    0%   0         [0, 1.6]   AMANK       acmana_c + atp_c --> acmanap_c + adp_...
    0%   0         [0, 1.6]   ANDRSTR...  andrstrnglc_c + atp_c + h2o_c --> adp...
    0%   0         [0, 1.6]   ARGSS       asp_L_c + atp_c + citr_L_c --> amp_c ...
    0%   0         [0, 1.6]   ATP2ter     amp_r + atp_c <=> amp_c + atp_r
    0%   0         [0, 1.6]   ATPasel     atp_c + h2o_c --> adp_c + h_l + pi_c
    0%   0         [0, 1.6]   ATPtm       adp_c + atp_m --> adp_m + atp_c
    0%   0         [0, 1.6]   BILDGLC...  atp_c + bildglcur_c + h2o_c --> adp_c...
    0%   0         [0, 1.6]   BILGLCURte  atp_c + bilglcur_c + h2o_c --> adp_c ...
    0%   0         [0, 1.6]   CAATPS      atp_c + ca2_c + h2o_c --> adp_c + ca2...
    0%   0         [0, 1.6]   CBPS        2.0 atp_c + gln_L_c + h2o_c + hco3_c ...
    0%   0         [0, 1.6]   CHSTEROLt   atp_c + chsterol_c + h2o_c --> adp_c ...
    0%   0         [0, 1.6]   CKc         atp_c + creat_c <=> adp_c + h_c + pcr...
    0%   0         [0, 1.6]   CTPS1       atp_c + nh4_c + utp_c --> adp_c + ctp...
    0%   0         [0, 1.6]   CTPS2       atp_c + gln_L_c + h2o_c + utp_c --> a...
    0%   0         [0, 1.6]   DGK1        atp_c + dgmp_c <=> adp_c + dgdp_c
    0%   0         [0, 1.6]   DM_atp_c_   atp_c + h2o_c --> adp_c + h_c + pi_c
    0%   0         [0, 1.6]   DOPAVESSEC  atp_c + dopa_c + h2o_c --> adp_c + do...
    0%   0         [0, 1.6]   DPCOAK      atp_c + dpcoa_c --> adp_c + coa_c + h_c
    0%   0         [0, 1.6]   DTMPK       atp_c + dtmp_c --> adp_c + dtdp_c
    0%   0         [0, 1.6]   DURIK1      atp_c + duri_c --> adp_c + dump_c + h_c
    0%   0         [0, 1.6]   FACOAL150   atp_c + coa_c + ptdca_c <=> amp_c + p...
    0%   0         [0, 1.6]   FACOAL161   atp_c + coa_c + hdcea_c --> amp_c + h...
    0%   0         [0, 1.6]   FACOAL170   atp_c + coa_c + hpdca_c <=> amp_c + h...
    0%   0         [0, 1.6]   FACOAL180i  atp_c + coa_c + ocdca_c --> amp_c + p...
    0%   0         [0, 1.6]   FACOAL1812  atp_c + coa_c + vacc_c <=> amp_c + oc...
    0%   0         [0, 1.6]   FACOAL1813  atp_c + coa_c + elaid_c <=> amp_c + o...
    0%   0         [0, 1.6]   FACOAL181i  atp_c + coa_c + ocdcea_c --> amp_c + ...
    0%   0         [0, 1.6]   FACOAL1821  atp_c + coa_c + lnlc_c <=> amp_c + ln...
    0%   0         [0, 1.6]   FACOAL1822  atp_c + coa_c + lneldc_c <=> amp_c + ...
    0%   0         [0, 1.6]   FACOAL1831  atp_c + coa_c + lnlncg_c <=> amp_c + ...
    0%   0         [0, 1.6]   FACOAL1832  atp_c + coa_c + lnlnca_c <=> amp_c + ...
    0%   0         [0, 1.6]   FACOAL184   atp_c + coa_c + strdnc_c <=> amp_c + ...
    0%   0         [0, 1.6]   FACOAL191   atp_c + coa_c + prist_c --> amp_c + p...
    0%   0         [0, 1.6]   FACOAL2...  atp_c + coa_c + tettet6_c <=> amp_c +...
    0%   0         [0, 1.6]   FACOAL2...  atp_c + coa_c + tetpent6_c <=> amp_c ...
    0%   0         [0, 1.6]   FACOAL2...  atp_c + coa_c + tetpent3_c <=> amp_c ...
    0%   0         [0, 1.6]   FACOAL2...  atp_c + coa_c + tethex3_c <=> amp_c +...
    0%   0         [0, 1.6]   FACOAL203   atp_c + coa_c + dlnlcg_c <=> amp_c + ...
    0%   0         [0, 1.6]   FACOAL2042  atp_c + coa_c + eicostet_c <=> amp_c ...
    0%   0         [0, 1.6]   FACOAL205   atp_c + coa_c + tmndnc_c <=> amp_c + ...
    0%   0         [0, 1.6]   FACOAL224   adrn_c + atp_c + coa_c <=> adrncoa_c ...
    0%   0         [0, 1.6]   FACOAL2251  atp_c + coa_c + dcsptn1_c <=> amp_c +...
    0%   0         [0, 1.6]   FACOAL2252  atp_c + clpnd_c + coa_c <=> amp_c + c...
    0%   0         [0, 1.6]   FACOAL226   atp_c + coa_c + crvnc_c <=> amp_c + c...
    0%   0         [0, 1.6]   FACOAL240   atp_c + coa_c + lgnc_c --> amp_c + pp...
    0%   0         [0, 1.6]   FACOAL241   atp_c + coa_c + nrvnc_c --> amp_c + n...
    0%   0         [0, 1.6]   FACOAL260   atp_c + coa_c + hexc_c <=> amp_c + he...
    0%   0         [0, 1.6]   FK          atp_c + fuc_L_c --> adp_c + fuc1p_L_c...
    0%   0         [0, 1.6]   FOLABCCte   atp_c + fol_c + h2o_c --> adp_c + fol...
    0%   0         [0, 1.6]   FPGS        atp_c + 4.0 glu_L_c + thf_c --> 5thf_...
    0%   0         [0, 1.6]   FPGS2       5thf_c + atp_c + glu_L_c --> 6thf_c +...
    0%   0         [0, 1.6]   FPGS3       6thf_c + atp_c + glu_L_c --> 7thf_c +...
    0%   0         [0, 1.6]   FPGS4       atp_c + dhf_c + 4.0 glu_L_c --> 5dhf_...
    0%   0         [0, 1.6]   FPGS5       5dhf_c + atp_c + glu_L_c --> 6dhf_c +...
    0%   0         [0, 1.6]   FPGS6       6dhf_c + atp_c + glu_L_c --> 7dhf_c +...
    0%   0         [0, 1.6]   FPGS7       10fthf_c + atp_c + 4.0 glu_L_c --> 10...
    0%   0         [0, 1.6]   FPGS8       10fthf5glu_c + atp_c + glu_L_c --> 10...
    0%   0         [0, 1.6]   FPGS9       10fthf6glu_c + atp_c + glu_L_c --> 10...
    0%   0         [0, 1.6]   FTHFL       atp_c + for_c + thf_c --> 10fthf_c + ...
    0%   0         [0, 1.6]   GLNS        atp_c + glu_L_c + nh4_c --> adp_c + g...
    0%   0         [0, 1.6]   GLUCYS      atp_c + cys_L_c + glu_L_c --> adp_c +...
    0%   0         [0, 1.6]   GLYK        atp_c + glyc_c --> adp_c + glyc3p_c +...
    0%   0         [0, 1.6]   GMPS2       atp_c + gln_L_c + h2o_c + xmp_c --> a...
    0%   0         [0, 1.6]   GTHS        atp_c + glucys_c + gly_c --> adp_c + ...
    0%   0         [0, 1.6]   HEX10       atp_c + gam_c --> adp_c + gam6p_c + h_c
    0%   0         [0, 1.6]   HEX4        atp_c + man_c --> adp_c + h_c + man6p_c
    0%   0         [0, 1.6]   HEX7        atp_c + fru_c --> adp_c + f6p_c + h_c
    0%   0         [0, 1.6]   KHK2        atp_c + xylu_D_c --> adp_c + h_c + xu...
    0%   0         [0, 1.6]   LINOFATPtc  atp_c + coa_c + lnlc_e --> amp_c + ln...
    0%   0         [0, 1.6]   METAT       atp_c + h2o_c + met_L_c --> amet_c + ...
    0%   0         [0, 1.6]   MI145PK     atp_c + mi145p_c --> adp_c + h_c + mi...
    0%   0         [0, 1.6]   NADK        atp_c + nad_c --> adp_c + h_c + nadp_c
    0%   0         [0, 1.6]   NADS2       atp_c + dnad_c + gln_L_c + h2o_c --> ...
    0%   0         [0, 1.6]   NDPK10      atp_c + didp_c <=> adp_c + ditp_c
    0%   0         [0, 1.6]   NDPK6       atp_c + dudp_c <=> adp_c + dutp_c
    0%   0         [0, 1.6]   NDPK9       atp_c + idp_c <=> adp_c + itp_c
    0%   0         [0, 1.6]   NMNATr      atp_c + h_c + nmn_c <=> nad_c + ppi_c
    0%   0         [0, 1.6]   NNATr       atp_c + h_c + nicrnt_c --> dnad_c + p...
    0%   0         [0, 1.6]   NaKt        atp_c + h2o_c + k_e + na1_c --> adp_c...
    0%   0         [0, 1.6]   OPAHir      5oxpro_c + atp_c + 2.0 h2o_c --> adp_...
    0%   0         [0, 1.6]   PACCOAL     atp_c + coa_c + pac_c --> amp_c + pha...
    0%   0         [0, 1.6]   PCHOLABCtc  atp_c + h2o_c + pchol_hs_c --> adp_c ...
    0%   0         [0, 1.6]   PEFLIP      atp_c + h2o_c + pe_hs_e --> adp_c + h...
    0%   0         [0, 1.6]   PEFLIPm     atp_c + h2o_c + pe_hs_c --> adp_c + h...
    0%   0         [0, 1.6]   PFK26       atp_c + f6p_c --> adp_c + f26bp_c + h_c
    0%   0         [0, 1.6]   PGLYCABCte  atp_c + h2o_c + pglyc_hs_c --> adp_c ...
    0%   0         [0, 1.6]   PI34P5K     atp_c + pail34p_hs_c --> adp_c + h_c ...
    0%   0         [0, 1.6]   PI45P3K     atp_c + pail45p_hs_c --> adp_c + 6.0 ...
    0%   0         [0, 1.6]   PI4P3K      atp_c + pail4p_hs_c --> adp_c + h_c +...
    0%   0         [0, 1.6]   PI5P3K      atp_c + pail5p_hs_c --> adp_c + h_c +...
    0%   0         [0, 1.6]   PIK3        atp_c + pail_hs_c --> adp_c + h_c + p...
    0%   0         [0, 1.6]   PNTK        atp_c + pnto_R_c --> 4ppan_c + adp_c ...
    0%   0         [0, 1.6]   PNTOt5      atp_c + h2o_c + na1_e + pnto_R_e --> ...
    0%   0         [0, 1.6]   PPNCL3      4ppan_c + atp_c + cys_L_c --> 4ppcys_...
    0%   0         [0, 1.6]   PRAGSr      atp_c + gly_c + pram_c <=> adp_c + ga...
    0%   0         [0, 1.6]   PRFGS       atp_c + fgam_c + gln_L_c + h2o_c --> ...
    0%   0         [0, 1.6]   PRISTCO...  atp_c + h2o_c + pristcoa_c --> adp_c ...
    0%   0         [0, 1.6]   PRPPS       atp_c + r5p_c --> amp_c + h_c + prpp_c
    0%   0         [0, 1.6]   PSFLIP      atp_c + h2o_c + ps_hs_e --> adp_c + h...
    0%   0         [0, 1.6]   PSFLIPm     atp_c + h2o_c + ps_hs_c --> adp_c + h...
    0%   0         [0, 1.6]   PSHSABCtc   atp_c + h2o_c + ps_hs_c --> adp_c + h...
    0%   0         [0, 1.6]   PTPAT       atp_c + h_c + pan4p_c --> dpcoa_c + p...
    0%   0         [0, 1.6]   PYDAMK      atp_c + pydam_c --> adp_c + h_c + pya...
    0%   0         [0, 1.6]   PYDXK       atp_c + pydx_c --> adp_c + h_c + pydx...
    0%   0         [0, 1.6]   PYDXNK      atp_c + pydxn_c --> adp_c + h_c + pdx...
    0%   0         [0, 1.6]   RBFK        atp_c + ribflv_c --> adp_c + fmn_c + h_c
    0%   0         [0, 1.6]   RTOTALF...  Rtotal_e + atp_c + coa_c + 4.0 h_c --...
    0%   0         [0, 1.6]   R_group...  Rtotal_c + atp_c + coa_c + 4.0 h_c --...
    0%   0         [0, 1.6]   SADT        atp_c + h_c + so4_c --> aps_c + ppi_c
    0%   0         [0, 1.6]   SCP21cx     atp_c + h2o_c + phytcoa_c --> adp_c +...
    0%   0         [0, 1.6]   SLCBK1      atp_c + sphgn_c --> adp_c + h_c + sph...
    0%   0         [0, 1.6]   SPHK21c     atp_c + sphings_c --> adp_c + h_c + s...
    0%   0         [0, 1.6]   TCHOLAt3    atp_c + h2o_c + tchola_c --> adp_c + ...
    0%   0         [0, 1.6]   TMDK1       atp_c + thymd_c --> adp_c + dtmp_c + h_c
    0%   0         [0, 1.6]   TSTSTER...  atp_c + h2o_c + tststeroneglc_c --> a...
    0%   0         [0, 1.6]   VITEtl      atp_c + avite1_c + h2o_c --> adp_c + ...
    0%   0         [0, 1.6]   r0283       ala_B_c + atp_c + his_L_c --> amp_c +...
    0%   0         [0, 1.6]   r0301       atp_c + nh4_c + xmp_c --> amp_c + gmp...
    0%   0         [0, 1.6]   r0345       atp_c + damp_c <=> adp_c + dadp_c
    0%   0         [0, 1.6]   r0408       atp_c + s7p_c <=> HC00361_c + adp_c +...
    0%   0         [0, 1.6]   r0578       atp_c + ptth_c <=> adp_c + h_c + pan4p_c
    0%   0         [0, 1.6]   r0666       atp_c + fpram_c <=> adp_c + air_c + 2...
    0%   0         [0, 1.6]   r0679       HC01231_c + atp_c <=> 4ppcys_c + adp_...
    0%   0         [0, 1.6]   r0812       amp_x + atp_c <=> amp_c + atp_x
    0%   0         [0, 1.6]   r0853       atp_c + pi_m --> atp_m + pi_c
    0%   0         [0, 1.6]   r0860       adp_c + atp_x <=> adp_x + atp_c
    0%   0         [0, 1.6]   r1254       atp_c + coa_c + 4.0 fadh2_c + strdnc_...
    0%   0         [0, 1.6]   r1382       6.0 atp_c + 6.0 glu_L_c + 6.0 h2o_c +...
    0%   0         [0, 1.6]   r1488       atp_c + coa_c + hdcea_c <=> CE0852_c ...
    0%   0         [0, 1.6]   r1514       arachd_c + atp_c + h2o_c --> adp_c + ...
    0%   0         [0, 1.6]   r1515       atp_c + h2o_c + hdca_c --> adp_c + h_...
    0%   0         [0, 1.6]   r1516       atp_c + h2o_c + ocdcea_c --> adp_c + ...
    0%   0         [0, 1.6]   r1517       atp_c + h2o_c + strdnc_c --> adp_c + ...
    0%   0         [0, 1.6]   r1518       atp_c + h2o_c + lnlc_c --> adp_c + h_...
    0%   0         [0, 1.6]   r1519       atp_c + elaid_c + h2o_c --> adp_c + e...
    0%   0         [0, 1.6]   r1520       atp_c + h2o_c + lnlncg_c --> adp_c + ...
    0%   0         [0, 1.6]   r1521       atp_c + h2o_c + lnlnca_c --> adp_c + ...
    0%   0         [0, 1.6]   r1522       atp_c + h2o_c + lgnc_c --> adp_c + h_...
    0%   0         [0, 1.6]   r1523       atp_c + h2o_c + hdcea_c --> adp_c + h...
    0%   0         [0, 1.6]   r1525       atp_c + dlnlcg_c + h2o_c --> adp_c + ...
    0%   0         [0, 1.6]   r1529       atp_c + h2o_c + ttdca_c --> adp_c + h...
    0%   0         [0, 1.6]   r2505       HC02198_c + atp_c + h2o_c --> HC02198...
    0%   0         [0, 1.6]   r2517       atp_c + h2o_c + thcholstoic_c --> adp...
    0%   0         [0, 1.6]   r2518       atp_c + dhcholestanate_c + h2o_c --> ...



```python
for rx in model_MeanNormalBiopsy.reactions:
    print(rx.id)
```

    13DAMPPOX
    2HBO
    34DHOXPEGOX
    3HBCDm
    3SALACBOXL
    3SALATAim
    3SPYRSPm
    41R1H2MAE12BOOX
    41R2A1H12BOOX
    4HOXPACDOX_NADP_
    A_MANASEly
    AACOAT
    ABTD
    ABUTD
    ACACT10m
    ACACT1rm
    ACACT4p
    ACACT5p
    ACACT6p
    ACACT7p
    ACCOAC
    ACGAM6PSi
    ACGAMK
    ACGAMPM
    ACITL
    ACOAD10m
    ACOAD9m
    ACOAO7p
    ACONT
    ACONTm
    ADK1
    ADK1m
    ADK3m
    ADMDC
    ADNK1
    ADNK1m
    ADPT
    ADRNCPT1
    ADRNCPT2
    ADSK
    ADSL1
    ADSL2
    ADSS
    AG13T10g
    AG13T11g
    AG13T12g
    AG13T13g
    AG13T14g
    AG13T15g
    AG13T4g
    AG13T5g
    AG13T6g
    AG13T7g
    AG13T8g
    AG13T9g
    AGDC
    AGPAT1
    AGTim
    AGTix
    AHEXASE2ly
    AICART
    AKGDm
    AKR1C1
    AKR1C42
    AKR1D2
    ALASm
    ALCD1
    ALCD21_L
    ALCD22_L
    ALCD2if
    ALCD2yf
    ALDD21
    ALDD2x
    ALDD2xm
    ALDD2y
    ALR2
    ALR3
    AMACR2p
    AMPDA
    APAT2rm
    ARABR
    ARACHCPT1
    ARACHCPT2
    ARGN
    ARGNm
    R_group_phosphotase_3
    ARTPLM3
    ASAH1
    ASPTA
    ASPTAm
    ATPH1e
    B_MANNASEly
    B3GALTg
    B3GNT51g
    BAAT1x
    BAMPPALDOX
    BDHm
    BDMT_U
    BILIRED
    BPNT
    BPNT2
    C14STRr
    C161CPT1
    C180CPT1
    C180CPT2
    C181CPT1
    C181CPT2
    C204CPT1
    C204CPT2
    C226CPT1
    C226CPT2
    C2M26DCOAHLm
    C3STDH1Pr
    C3STDH1r
    C3STKR2r
    C4STMO1r
    C4STMO2Pr
    C4STMO2r
    CAT2p
    CATm
    CATp
    CDIPTr
    CDS
    CDSm
    CHSTEROLSULT
    CITL
    CKc
    CLS_hs
    CMPSAS
    CMPSASn
    CORE6GTg
    CPPPGO
    CRTNsyn
    CSm
    CSNAT2x
    CSNATp
    CSNATr
    CTPS1
    CTPS2
    CYSO
    CYSTA
    CYSTAm
    CYTDK1
    CYTDK2m
    CYTK1
    CYTK10
    CYTK10n
    CYTK11
    CYTK11n
    CYTK12
    CYTK12n
    CYTK13
    CYTK13n
    CYTK14
    CYTK14n
    CYTK1m
    CYTK1n
    CYTK2
    CYTK2n
    CYTK3
    CYTK3n
    CYTK4
    CYTK4n
    CYTK5
    CYTK5n
    CYTK6
    CYTK6n
    CYTK7
    CYTK7n
    CYTK8
    CYTK8n
    CYTK9
    CYTK9n
    DADNK
    DAGK_hs
    DASCBR
    DCSPTN1CPT1
    DCSPTN1CPT2
    DCYTD
    DDPGAm
    DESAT16_2
    DESAT18_10
    DESAT18_3
    DESAT18_4
    DESAT18_5
    DESAT18_9
    DESAT20_2
    DESAT22_1p
    DGAT
    DGK1
    DGK2m
    DGNSKm
    DGULND
    DHCR241r
    DHCR242r
    DHCR243r
    DHCR71r
    DHCR72r
    DHCRD1
    DHCRD2
    DLNLCGCPT1
    DLNLCGCPT2
    DM_atp_c_
    DMATTx
    DOLASNT_Uer
    DOLDPP_Uer
    DOLPGT1_Uer
    DOLPGT2_Uer
    DOLPGT3_Uer
    DOLPMT_U
    DOLPMT1_Uer
    DOLPMT2_Uer
    DOLPMT3_Uer
    DOLPMT4_Uer
    DOPABMO
    DPCOAK
    DPGase
    DPGM
    DPMVDx
    DRPA
    DSAT
    DURIK1
    DURIPP
    DUTPDPm
    DUTPDPn
    EBP1r
    EBP2r
    ECOAH12m
    ECOAH1m
    ECOAH9m
    EHGLAT2m
    EHGLATm
    EICOSTETCPT1
    EICOSTETCPT2
    ELAIDCPT1
    ELAIDCPT2
    ENGASE2ly
    ENGASE3ly
    ENO
    F1PGT
    F6Tg
    FACOAL160i
    FACOAL180i
    FACOAL1812
    FACOAL1813
    FACOAL1821
    FACOAL1831
    FACOAL1832
    FACOAL203
    FACOAL204
    FACOAL2042
    FACOAL205
    FACOAL224
    FACOAL2252
    FACOAL226
    FACOAL80i
    FAEL183
    FAEL184
    FAEL204
    FAEL205
    FAH3
    FALDH
    FAOXC180x
    FAOXC183806m
    FAOXC200180m
    FAOXC200180x
    FAOXC2031836m
    FAOXC2051843m
    FAOXC2242046m
    FAOXC2242046x
    FAOXC226205m
    FAOXC226205x
    FBA
    FBA2
    FBA4
    FBP26
    FCLTm
    FE3R2e
    FK
    FTHFL
    FUCASE2ly
    FUM
    FUMm
    FUT31g
    G12MT1_U
    G12MT2_U
    G13MT_U
    G14T10g
    G14T11g
    G14T12g
    G14T13g
    G14T14g
    G14T15g
    G14T16g
    G14T17g
    G14T6g
    G14T7g
    G14T8g
    G14T9g
    G14Tg
    G16MT_U
    G6PDA
    G6PDH2r
    GALASE10ly
    GALASE11ly
    GALASE12ly
    GALASE13ly
    GALASE14ly
    GALASE15ly
    GALASE19ly
    GALASE3ly
    GALASE4ly
    GALASE5ly
    GALASE6ly
    GALASE7ly
    GALASE8ly
    GALASE9ly
    GALNTg
    GALOR
    GALU
    GAPD
    GARFT
    GASNASE2ly
    GASNASE3ly
    GBA
    GBAl
    GCALDD
    GCC2cm
    GFUCS
    GGNG
    GGT5r
    GHMT2rm
    GK1
    GLB1
    GLBRAN
    GLCAASE8ly
    GLCAASE9ly
    GLCNACPT_U
    GLCNACT_U
    GLDBRAN
    GLGNS1
    GLNS
    GLPASE1
    GLPASE2
    GLUCYS
    GLUDxm
    GLUDym
    GLUPRT
    GLXO1
    GLYAMDTRc
    GLYCTO1p
    GLYOp
    GMAND
    GMPS2
    GPAM_hs
    GPDDA1
    GRTTx
    GTHO
    GTHOm
    GTHP
    GTHPe
    GTHPm
    GTHS
    GULN3D
    GULNDer
    H2CO3D
    H2O2syn
    HACD1m
    HACD1x
    HACD9m
    HEX1
    HEX10
    HEX4
    HEX7
    HIBDm
    HISDC
    HMBS
    HMGCOASi
    HMGLm
    HMGLx
    HOXG
    HPYRDC
    HPYRR2x
    HSD17B2r
    HSD17B42x
    HSD17B4x
    HSD3A1r
    HSD3A2r
    HSD3B11
    HSD3B11r
    HSD3B12r
    HSD3B13r
    HYPOE
    HYPTROX
    ICDHy
    ICDHyrm
    ILETA
    ILETAm
    IMPC
    IMPD
    IPDDIx
    KHK2
    LALDO
    LALDO2x
    LCADi
    LCADi_D
    LCYSTATm
    LDH_L
    LDH_Lm
    LGTHL
    LNELDCCPT2
    LNLCCPT1
    LNLCCPT2
    LNLNCACPT1
    LNLNCACPT2
    LNLNCGCPT1
    LNLNCGCPT2
    LNSTLSr
    LPASE
    LPS3
    LPS4e
    LSTO1r
    LSTO2r
    LTA4H
    LTC4Sr
    M1316Mg
    M13N2Tg
    M13N4Tg
    M14NTg
    M16N4Tg
    M16N6Tg
    M16NTg
    MACOXO
    MAN1PT2
    MAN6PI
    MAOX
    MCCCrm
    MCDp
    MCLOR
    MDH
    MDHm
    ME2
    METAT
    MEVK1x
    MG1er
    MG2er
    MG3er
    MGCHrm
    MGSA
    MGSA2
    MHISOR
    MI1345PP
    MI134P4P
    MI134PP
    MI13PP
    MI145PK
    MI145PP
    MI14PP
    MI1PP
    MI34PP
    MI3PP
    MI4PP
    MM5ag
    MM5bg
    MM5cg
    MM6B1ag
    MM6B1bg
    MM6B2g
    MM6bg
    MM7B1g
    MM7B2g
    MM7Cag
    MM7Cbg
    MM8Ag
    MM8Ber
    MM8Cg
    MMMm
    MMSAD3m
    MMTSADm
    MTHFCm
    MTHFD2
    MTHFD2m
    MTHFDm
    N4Tg
    NACHEX10ly
    NACHEX11ly
    NACHEX12ly
    NACHEX13ly
    NACHEX14ly
    NACHEX15ly
    NACHEX16ly
    NACHEX17ly
    NACHEX18ly
    NACHEX19ly
    NACHEX20ly
    NACHEX21ly
    NACHEX22ly
    NACHEX27ly
    NADK
    NADS2
    NAGA2ly
    NAGLCAly
    NBAHH_ir
    NDP6
    NDP7er
    NDP7ex
    NDP7g
    NDP8ex
    NDPK10
    NDPK10n
    NDPK3m
    NDPK4n
    NDPK6
    NDPK6n
    NDPK9
    NDPK9n
    NICRNS
    NORANMT
    NOS1
    NP1
    NT5C
    NTD1
    NTD10
    NTD11
    NTD2
    NTD2e
    NTD3
    NTD4
    NTD5
    NTD6
    NTD7
    NTD8
    NTD9
    OIVD1m
    OIVD2m
    ORNDC
    ORNTArm
    P45011A1m
    P45017A1r
    P45017A2r
    P45017A3r
    P45017A4r
    P45027A13m
    P45027A1m
    P4508B11r
    P4508B13r
    P450LTB4r
    P450SCC1m
    PACCOAL
    PCHOLP_hs
    PCHOLPg_hs
    PCHOLPm_hs
    PCHOLPr_hs
    PDXPP
    PEAMNO
    PETOHMm_hs
    PFK
    PFK26
    PGCD
    PGK
    PGL
    PGM
    PGMT
    PHACCOAGLNAC
    PHCDm
    PHETA1m
    PHYCBOXL
    PI345P3P
    PI345P3Pn
    PI34P5K
    PI34P5Kn
    PI3P4K
    PI3P4Kn
    PI45P3K
    PI45P3Kn
    PI45P5P
    PI45PLC
    PI4P3K
    PI4P3Kn
    PI4P5Kn
    PI4PLC
    PI5P3K
    PI5P4Kn
    PIK3
    PIK3n
    PIK4
    PIK4n
    PIK5n
    PIPLC
    PMANM
    PMEVKx
    PNTEH
    PNTK
    PPA
    PPAm
    PPAn
    PPAP
    PPBNGS
    PPCOAOm
    PPM
    PPNCL3
    PPPGOm
    PRAGSr
    PRDX
    PRFGS
    PRGNLONESULT
    PRPNCOAHYDm
    PRPPS
    PSDm_hs
    PSERT
    PSP_L
    PSSA1_hs
    PTRCOX1
    PUNP1
    PUNP2
    PUNP3
    PUNP4
    PUNP6
    PYDAMK
    PYDXK
    PYDXNK
    PYDXPP
    PYK
    PYNP2r
    RADH2
    RAI1
    RAI2
    RAI3
    RDH1
    RDH1a
    RDH2
    RDH2a
    RDH3
    RDH3a
    RETI1
    RETI2
    RNDR3
    RNDR4
    S23T3g
    S6T10g
    S6T11g
    S6T12g
    S6T13g
    S6T14g
    S6T15g
    S6T4g
    S6T5g
    S6T6g
    S6T7g
    S6T8g
    S6T9g
    S6TASE10ly
    S6TASE11ly
    S6TASE12ly
    S6TASE13ly
    S6TASE14ly
    S6TASE15ly
    S6TASE16ly
    S6TASE17ly
    S6TASE18ly
    S6TASE19ly
    S6TASE1ly
    S6TASE20ly
    S6TASE21ly
    S6TASE23ly
    S6TASE24ly
    S6TASE26ly
    S6TASE2ly
    S6TASE3ly
    SADT
    SAMHISTA
    SBTD_D2
    SBTR
    SCP2x
    SCPx
    SFGTH
    SGPL12r
    SIAASE2ly
    SLDx
    SLDxm
    SMS
    SPHK21c
    SPMDOX
    SPODM
    SPODMm
    SPODMn
    SPODMx
    SPRMS
    SPTix
    SQLEr
    SQLSr
    SR5AR2r
    SR5ARr
    SUCD1m
    T2M26DCOAHLm
    TALA
    TKT1
    TKT2
    TMDK1
    TPI
    TRDR
    TRDR2
    TRDR3
    TYRCBOX
    TYROXDAc
    TYRTAm
    UAGDP
    UDPDOLPT_U
    UDPG4E
    UDPGD
    UDPGNP
    UGCG
    UGT1A10r
    UGT1A2r
    UGT1A3r
    UMPK
    UMPK2
    UMPK2n
    UMPK3
    UMPK3n
    UMPK4
    UMPK4n
    UMPK5
    UMPK5n
    UMPK6
    UMPK6n
    UMPK7
    UMPK7n
    UMPKn
    UPP3S
    UPPDC1
    URIK1
    VLCSp
    VLCSr
    XYLTD_Dr
    XYLUR
    r0009
    r0010
    r0021
    r0022
    r0024
    r0028
    r0047
    r0051
    r0055
    r0119
    r0129
    r0130
    r0139
    r0142
    r0149
    r0160
    r0166
    r0170
    r0173
    r0181
    r0208
    r0245
    r0246
    r0249
    r0268
    r0283
    r0287
    r0301
    r0309
    r0310
    r0317
    r0340
    r0345
    r0354
    r0355
    r0357
    r0358
    r0360
    r0361
    r0363
    r0364
    r0365
    r0381
    r0383
    r0386
    r0391
    r0392
    r0393
    r0400
    r0407
    r0408
    r0410
    r0426
    r0430
    r0431
    r0433
    r0438
    r0440
    r0444
    r0445
    r0456
    r0463
    r0464
    r0494
    r0497
    r0510
    r0511
    r0537
    r0539
    r0545
    r0547
    r0549
    r0555
    r0568
    r0570
    r0575
    r0578
    r0579
    r0580
    r0584
    r0591
    r0595
    r0596
    r0604
    r0614
    r0618
    r0629
    r0633
    r0634
    r0637
    r0639
    r0641
    r0642
    r0643
    r0648
    r0649
    r0652
    r0653
    r0656
    r0660
    r0661
    r0666
    r0668
    r0669
    r0670
    r0671
    r0679
    r0680
    r0686
    r0688
    r0698
    r0706
    r0714
    r0715
    r0716
    r0717
    r0718
    r0719
    r0720
    r0721
    r0722
    r0723
    r0724
    r0726
    r0727
    r0729
    r0730
    r0731
    r0732
    r0733
    r0734
    r0739
    r0741
    r0743
    r0744
    r0750
    r0758
    r0774
    r0779
    r0781
    r0782
    r0783
    r0786
    r0787
    r0795
    r1109
    r1135
    r1146
    r1167
    r1169
    r1177
    r1181
    r1378
    r1380
    r1384
    r1418
    r1443
    r1444
    r1445
    r1446
    r1447
    r1448
    r1449
    r1450
    r1451
    r1457
    r1481
    r1488
    RE0383C
    RE0452M
    RE0512M
    RE0512X
    RE0565C
    RE0566C
    RE0567C
    RE0568C
    RE0583C
    RE0688C
    RE0827C
    RE0908C
    RE0912C
    RE1062M
    RE1100C
    RE1135C
    RE1447N
    RE1448N
    RE1516M
    RE1516X
    RE1517M
    RE1517X
    RE1520M
    RE1521M
    RE1522M
    RE1523M
    RE1523X
    RE1525M
    RE1525X
    RE1526M
    RE1526X
    RE1527M
    RE1527X
    RE1530C
    RE1531M
    RE1531X
    RE1532M
    RE1532X
    RE1533M
    RE1533X
    RE1534M
    RE1534X
    RE1573M
    RE1573X
    RE1635R
    RE1691M
    RE1796M
    RE1815R
    RE1816R
    RE1817R
    RE1818R
    RE1834C
    RE1860E
    RE2051C
    RE2051G
    RE2051R
    RE2272L
    RE2318R
    RE2319R
    RE2453M
    RE2454M
    RE2514C
    RE2514L
    RE2541L
    RE2624X
    RE2625C
    RE2626M
    RE2680C
    RE2814M
    RE2814R
    RE2908M
    RE2908X
    RE2909M
    RE2909X
    RE2910M
    RE2910X
    RE2913M
    RE2913X
    RE2914M
    RE2914X
    RE2915M
    RE2915X
    RE2916M
    RE2916X
    RE2917M
    RE2917X
    RE2919M
    RE2920M
    RE2921M
    RE2954C
    RE2985X
    RE2986X
    RE2987X
    RE2988X
    RE2989X
    RE2990X
    RE2991X
    RE2992X
    RE2993X
    RE2994X
    RE2996X
    RE2997X
    RE2998M
    RE3001M
    RE3002X
    RE3004M
    RE3005M
    RE3006M
    RE3011R
    RE3012C
    RE3012M
    RE3012R
    RE3076X
    RE3082X
    RE3088X
    RE3090X
    RE3093X
    RE3097X
    RE3106C
    RE3110C
    RE3112C
    RE3113C
    RE3119C
    RE3120C
    RE3121C
    RE3122C
    RE3139X
    RE3141X
    RE3156X
    RE3158X
    RE3177M
    RE3178M
    RE3179M
    RE3189M
    RE3190M
    RE3191M
    RE3192M
    RE3193M
    RE3194M
    RE3224C
    RE3225C
    RE3226C
    RE3227C
    RE3228C
    RE3229C
    RE3230C
    RE3231C
    RE3232C
    RE3234C
    RE3235C
    RE3236C
    RE3237C
    RE3241C
    RE3241R
    RE3242C
    RE3242R
    RE3243C
    RE3243R
    RE3244C
    RE3244R
    RE3245C
    RE3247X
    RE3248M
    RE3248X
    RE3250M
    RE3250X
    RE3269C
    RE3270C
    RE3273G
    RE3273R
    RE3301G
    RE3301R
    RE3336M
    RE3336X
    RE3337M
    RE3338C
    RE3338M
    RE3338X
    RE3340C
    RE3340M
    RE3340X
    RE3342M
    RE3342X
    RE3344M
    RE3346C
    RE3346M
    RE3346R
    RE3347C
    RE3381L
    RE3383M
    RE3384M
    RE3386M
    RE3387M
    RE3388M
    RE3390M
    RE3391M
    RE3392M
    RE3393M
    RE3394M
    RE3396M
    RE3398M
    RE3399M
    RE3400M
    RE3401M
    RE3402M
    RE3403M
    RE3404M
    RE3443M
    RE3443X
    RE3447M
    RE3448C
    RE3448M
    RE3448X
    RE3496N
    RE3521R
    RE3532R
    RE3533R
    RE3534R
    RE3554R
    RE3557R
    RE3559M
    RE3559X
    RE3563M
    RE3564C
    RE3564M
    RE3564X
    RE3572X
    RE3573X
    RE3575X
    RE3587N
    RE3626M
    HSD17B3r
    ADPCOACROT
    C10OHc
    C12OHc
    C141OHc
    C142OHc
    C14OHc
    C161OHc
    C162OHc
    C16OHc
    C181OHc
    C182OHc
    C18OHc
    C4OHc
    C4x
    C50CPT1
    C60CPT1
    C6COAt
    FAOXC101C8m
    FAOXC101C8x
    FAOXC101m
    FAOXC102C101m
    FAOXC102C103m
    FAOXC102C81m
    FAOXC102m
    FAOXC103C102m
    FAOXC10DCC8DCx
    FAOXC121C101m
    FAOXC122C101m
    FAOXC122m
    FAOXC123C102m
    FAOXC123m
    FAOXC12DCC10DCx
    FAOXC12DCx
    FAOXC141C121m
    FAOXC141C141OHm
    FAOXC142C142OHm
    FAOXC143C123m
    FAOXC14C14OHm
    FAOXC14DCC12DCx
    FAOXC161C141m
    FAOXC161C161OHm
    FAOXC162C142m
    FAOXC162C162OHm
    FAOXC163C143m
    FAOXC163C164Gm
    FAOXC163GC142m
    FAOXC163Gm
    FAOXC164C143m
    FAOXC164C165m
    FAOXC164GC163m
    FAOXC164m
    FAOXC165C164m
    FAOXC16C16OHm
    FAOXC16DCC14DCx
    FAOXC16DCr
    FAOXC16OHC16r
    FAOXC181C161m
    FAOXC181C181OHm
    FAOXC182C162m
    FAOXC182C182OHm
    FAOXC183C163Gm
    FAOXC183C163m
    FAOXC184C163m
    FAOXC184C164m
    FAOXC184m
    FAOXC185C164m
    FAOXC185m
    FAOXC18C18OHm
    FAOXC204C184m
    FAOXC205C185m
    FAOXC225C204m
    FAOXC225C226m
    FAOXC225m
    FAOXC226C205m
    FAOXC226C225m
    FAOXC226C227m
    FAOXC226m
    FAOXC227C226m
    FAOXC5C5OHm
    FAOXC5OHc
    FAOXC61m
    FAOXC6DCC4DCx
    FAOXC81C61m
    FAOXC8DCC6DCx
    FAOXOHC16C16DCc
    FAOXTC101TC102m
    FAOXTC102C101m
    FAOXTC122TC101m
    FAOXTC122m
    FAOXTC142TC122m
    FAOXTC162TC142m
    FAOXTC182TC162m
    OCD11COACPT1
    OCD11CRNCPT2
    OCTDECCPT1
    OCTDECCPT2
    SEBCOACROT
    SUBERCROT
    SUCCCROT
    SUCCOAPET
    ALAALACNc
    GLYGLYCNc
    GLYLEUHYDROc
    GLYSARCNc
    ADNK3
    ADNK4
    DPMVDc
    DSREDUCr
    FERO
    GNDc
    HMGCOARc
    INSK
    IPDDI
    MEVK1c
    PMEVKc
    NADKm
    34DHOXPEGt
    3MLDAt
    3MOBt2im
    3MOPt2im
    3SALAASPm
    4MOPt2im
    4MPTNLte
    4MPTNLtm
    5ADTSTSTERONEte
    5ADTSTSTERONEtr
    5AOPtm
    ABTti
    ACACt2m
    ACACtx
    ACALDtm
    ACALDtr
    ACALDtx
    ACGAMtly
    ACNAMlt
    ACNAMtn
    ACt2m
    ADNtm
    ADRNCOAtx
    ADRNCRNt
    ADRNt
    AHCYStr
    ALACYSNaEx
    ALASERNaEx
    ALAt4
    ALATHRNaEx
    AMETr
    ANDRSTRNGLCte
    ANDRSTRNGLCtr
    ANDRSTRNte
    ANDRSTRNtr
    APRGSTRNte
    ARAB_Lt
    ARACHCOAtx
    ARACHCRNt
    ARACHDCOAtx
    ARGtiDF
    ARTFR12
    ASCBt
    ASNt4
    Asn_X_Ser_Thrtr
    ATP1ter
    ATP2ter
    ATPasel
    ATPS4m
    ATPtm
    AVITE1t
    BALAtmr
    BHBtm
    BILDGLCURtr
    BILGLCURtr
    BILIRUBtr
    C180CRNt
    C204CRNt
    C226COAtx
    C226CRNt
    CAATPS
    CDPDAGtm
    CERT1rt
    CERT2rt
    CHLtm
    CHOLATEt
    CHOLATEt3
    CHOLtg
    CHOLtr
    CHSTEROLt
    CHSTEROLt2
    CITt4_2
    CLOXAtex2
    CLPNDt
    CMPACNAtg
    CMPACNAtn
    CO2ter
    CO2tm
    CO2tp
    COAtm
    COAtp
    COAtr
    COt
    CREATt4_2_r
    CRNCARtp
    CRVNCtr
    CTPtn
    CYOOm2
    CYOR_u10m
    CYSALANaEx
    CYSSERNaEx
    CYSTHRNaEx
    CYTDt4
    CYTDtm
    DAGt
    DATPtn
    DCSPTN1COAtx
    DCSPTN1CRNt
    DCTPtn
    DGSNtm
    DGTPtn
    DHAAt1r
    DIDPtn
    DITPtn
    D_LACt2
    DLNLCGCRNt
    DNDPt27m
    DNDPt29m
    DNDPt34m
    DNDPt35m
    DNDPt37m
    DNDPt39m
    DNDPt40m
    DNDPt41m
    DNDPt42m
    DNDPt43m
    DNDPt55m
    DNDPt57m
    DOLGLCP_Lter
    DOLMANP_Lter
    DOLP_Uter
    DOPAtu
    DOPAVESSEC
    DTDPtn
    DTTPtn
    DUDPtn
    DUMPtn
    EICOSTETCRNt
    EICOSTETt
    ELAIDCRNt
    ETOHt
    ETOHtx
    F1Atg
    FADH2tru
    FADH2tx
    FADtru
    FADtx
    FE2tm
    FORtr
    FRDPtr
    FUC14GALACGLCGALGLUSIDEte
    FUC14GALACGLCGALGLUSIDEtg
    FUCtly
    FUMSO3tm
    FUMSO4tm
    FUMtm
    G3PD2m
    GALGLUSIDEtg
    GALGLUSIDEtl
    GALt1r
    GALtly
    GAMt1r
    GCHOLAt3
    GCHOLAte
    GCHOLAtx
    GDPFUCtg
    GDPtg
    GGLUCT
    GLCMter
    GLCt4
    GLCter
    GLCtg
    GLCtly
    GLCURter
    GLCURtly
    GLNt4
    GLNtm
    GLUtr
    GLXtm
    GLXtp
    GLYBt4_2_r
    GLYCLTtp
    GLYC_St
    GLYCtm
    GLYt4
    GLYtm
    GLYtp
    GMPtg
    GTHRDtr
    GULNter
    H2O2t
    H2O2tm
    H2O2tn
    H2O2tp
    H2Ot
    H2Oter
    H2Otg
    H2Otly
    H2Otn
    H2Otp
    HAtly
    HDCAter
    HDD2COAtx
    HISt4
    HMGCOAtm
    HMGCOAtx
    HOMt4
    HPACtr
    HPYRtp
    Htg
    Htm
    Htr
    Htx
    IDPtn
    ILEt5m
    INSTt4
    IPDPtx
    ITPtn
    KCC2t
    KCCt
    KSIt
    KSItly
    LEUKTRA4tr
    LEUKTRD4tr
    LNELDCCRNt
    LNLCCRNt
    LNLNCACRNt
    LNLNCAt
    LNLNCGCRNt
    LPCHOLt
    LYStiDF
    M4MPDOL_Uter
    M7MASNBterg
    M8MASNterg
    MAGt
    MALSO3tm
    MALSO4tm
    MALtm
    MANt1r
    MANter
    MANtg
    MANtly
    MEOHt2
    N2M2NMASNt
    N2M2NMASNtly
    NADPHtru
    NADPHtxu
    NADPtru
    NADPtxu
    NaKt
    NAt
    NCKt
    NH4tp
    NRVNCt
    O2St
    O2Stm
    O2Stn
    O2Stx
    O2t
    O2ter
    O2tm
    O2tn
    O2tp
    OCCOAtm
    OCCOAtx
    OCTAt
    ORNtiDF
    PA_HSter
    PA_HStg
    PAIL_HStn
    PAIL45P_HStn
    PAIL4P_HStn
    PAPStg
    PAPtg
    PCHOL_HSter
    PCHOL_HStg
    PE_HSter
    PE_HStm
    PEFLIP
    PEt
    PGLYCt
    PHEACGLNt
    PHEMEtm
    PHEt4
    PIt2m
    PItg
    PItn
    PItx
    PPItr
    PPItx
    PPPG9tm
    PRGNLONEtm
    PRGNLONEtr
    PRISTANALtx
    PRISTtx
    PROSTGD2t
    PROSTGE1t
    PROSTGE2t
    PSFLIP
    PSFLIPm
    PSt3
    PYRt2p
    RETNt
    RETt
    RTOTALt
    Rtotaltl
    SCP21x
    SCP22x
    Ser_Thrtg
    SERALANaEx
    SERCYSNaEx
    SERt4
    SERTHRNaEx
    SERtp
    SO4OXAtex2
    SO4tl
    SPH1Pte
    SPHINGStl
    SPHS1Ptr
    STRDNCt
    SUCCt2m
    SULFOX
    TAURt4_2_r
    THMDt4
    THMt3
    THRALANaEx
    THRCYSNaEx
    THRSERNaEx
    THRt4
    TMNDNCCOAtx
    TMNDNCt
    TTDCAtr
    UDPACGALtl
    UDPGALtg
    UDPGLCAter
    UDPGLCter
    UDPGLCtg
    UDPtl
    UGALNACtg
    UGLCNACtg
    UREAt
    UREAt5
    UREAtm
    VACCt
    XOL27OHtm
    XOLDIOLONEt
    XYLTt
    r0002
    r0205
    r0801
    r0812
    r0818
    r0819
    r0821
    r0822
    r0830
    r0838
    r0840
    r0841
    r0853
    r0860
    r0892
    r0926
    r0940
    r0941
    r0942
    r0954
    r0963
    r0973
    r0974
    r0990
    r0997
    r1000
    r1002
    r1003
    r1004
    r1006
    r1010
    r1013
    r1014
    r1116
    r1117
    r1162
    r1291
    r1292
    r1318
    r1400
    r1421
    r1434
    r1459
    r1464
    r1514
    r1515
    r1516
    r1517
    r1518
    r1519
    r1520
    r1521
    r1522
    r1523
    r1525
    r1529
    r1544
    r1546
    r1547
    r1548
    r1549
    r1551
    r1552
    r1553
    r1554
    r1556
    r1557
    r1559
    r1560
    r1561
    r1562
    r1563
    r1564
    r1565
    r1566
    r1567
    r1568
    r1569
    r1570
    r1571
    r1573
    r1574
    r1575
    r1576
    r1578
    r1579
    r1580
    r1581
    r1583
    r1584
    r1585
    r1586
    r1587
    r1588
    r1589
    r1590
    r1591
    r1592
    r1593
    r1594
    r1595
    r1596
    r1597
    r1598
    r1599
    r1600
    r1602
    r1603
    r1604
    r1605
    r1606
    r1607
    r1608
    r1609
    r1610
    r1611
    r1627
    r1628
    r1629
    r1630
    r1631
    r1633
    r1634
    r1635
    r1636
    r1637
    r1638
    r1639
    r1640
    r1641
    r1642
    r1643
    r1645
    r1646
    r1648
    r1649
    r1650
    r1651
    r1652
    r1653
    r1654
    r1655
    r1656
    r1657
    r1658
    r1659
    r1660
    r1661
    r1662
    r2079
    r2080
    r2081
    r2082
    r2083
    r2085
    r2087
    r2088
    r2090
    r2093
    r2100
    r2101
    r2103
    r2104
    r2105
    r2109
    r2115
    r2120
    r2124
    r2126
    r2129
    r2133
    r2139
    r2313
    r2314
    r2315
    r2316
    r2317
    r2318
    r2319
    r2320
    r2321
    r2322
    r2323
    r2324
    r2325
    r2326
    r2327
    r2328
    r2329
    r2330
    r2331
    r2332
    r2333
    r2334
    r2335
    r2344
    r2353
    r2371
    r2372
    r2377
    r2378
    r2379
    r2404
    r2419
    r2449
    r2472
    r2508
    r2509
    r2513
    r2516
    r2517
    r2521
    r2537
    2MB2COAc
    3HBCOARc
    BCRNe
    C141ACBP
    C141OHe
    C142ACBP
    C142OHe
    C162ACBP
    C162OHe
    C6CRNe
    C8CRNe
    DDCRNe
    DECCRNe
    GLUTCOAACBP
    HDCACBP
    HDDACBP
    HDECAACBP
    HDECEACBP
    HEDCECRNe
    HEXCOAACBP
    HEXDCRNe
    HEXDIACTD
    HEXDICOAACBP
    HEXDICOAACBPx
    HIVCACBP
    HIVCRNe
    HOCDACBP
    HOCTDACBP
    HOCTDEC2CRNe
    HOCTDECCRNe
    HTDCACBP
    HTDCRNe
    IVCOAACBP
    IVCRNe
    OCD11CRNCACT
    OCTDEC2ACBP
    OCTDECCACT
    OCTDECE1CRNe
    OMHPALTD
    SUCCTD
    4OHPROIMINOtc
    CARPEPT1tc
    LINOFATPtc
    NACHORCTL3le
    VITEtl
    4HPROLTASCT1
    AHCYStd
    AMETtd
    ASPPROASCT1
    CHSTEROLtrc
    CYTDt5le
    DATPtm
    DHAPtc
    FE2DMT1
    FE3MTP1
    GLUPROASCT1
    GLYC3Ptmc
    LACLt
    NADtm
    NODe
    PPItm
    RTOTALFATPc
    SPMTDe
    THMDt5le
    q10h2tc
    q10tm
    3MOBte
    4HPRO_LTte
    4MOPte
    5MTAte
    5OXPROt
    AHCYSte
    AICARte
    MAL_Lte
    fumt
    xmpt
    dcmpt
    glyc3pte
    DM_datp_n_
    DM_dctp_m_
    DM_dctp_n_
    DM_dgtp_m_
    DM_dgtp_n_
    DM_dttp_n_
    DM_ethamp_r_
    DM_n5m2masn_g_
    sink_octdececoa_LPAREN_c_RPAREN_
    DM_taur_LPAREN_c_RPAREN_
    DM_pe_hs_LPAREN_r_RPAREN_
    DM_4hrpo
    DM_Lcystin
    DM_fol
    EX_1glyc_hs_LPAREN_e_RPAREN_
    EX_2hb_LPAREN_e_RPAREN_
    EX_34dhoxpeg_LPAREN_e_RPAREN_
    EX_3mlda_LPAREN_e_RPAREN_
    EX_4hphac_LPAREN_e_RPAREN_
    EX_4mptnl_LPAREN_e_RPAREN_
    EX_5adtststerone_LPAREN_e_RPAREN_
    EX_abt_LPAREN_e_RPAREN_
    EX_ac_LPAREN_e_RPAREN_
    EX_acac_LPAREN_e_RPAREN_
    EX_acald_LPAREN_e_RPAREN_
    EX_acetone_LPAREN_e_RPAREN_
    EX_acgam_LPAREN_e_RPAREN_
    EX_ade_LPAREN_e_RPAREN_
    EX_adn_LPAREN_e_RPAREN_
    EX_adrn_LPAREN_e_RPAREN_
    EX_akg_LPAREN_e_RPAREN_
    EX_ala_B_LPAREN_e_RPAREN_
    EX_ala_D_LPAREN_e_RPAREN_
    EX_ala_L_LPAREN_e_RPAREN_
    EX_amp_LPAREN_e_RPAREN_
    EX_andrstrn_LPAREN_e_RPAREN_
    EX_andrstrnglc_LPAREN_e_RPAREN_
    EX_aprgstrn_LPAREN_e_RPAREN_
    EX_arab_L_LPAREN_e_RPAREN_
    EX_arachd_LPAREN_e_RPAREN_
    EX_arg_L_LPAREN_e_RPAREN_
    EX_asn_L_LPAREN_e_RPAREN_
    EX_asp_L_LPAREN_e_RPAREN_
    EX_atp_LPAREN_e_RPAREN_
    EX_bhb_LPAREN_e_RPAREN_
    EX_bildglcur_LPAREN_e_RPAREN_
    EX_bilglcur_LPAREN_e_RPAREN_
    EX_btn_LPAREN_e_RPAREN_
    EX_ca2_LPAREN_e_RPAREN_
    EX_caro_LPAREN_e_RPAREN_
    EX_cgly_LPAREN_e_RPAREN_
    EX_chol_LPAREN_e_RPAREN_
    EX_cholate_LPAREN_e_RPAREN_
    EX_chsterol_LPAREN_e_RPAREN_
    EX_cit_LPAREN_e_RPAREN_
    EX_cl_LPAREN_e_RPAREN_
    EX_CLPND_LPAREN_e_RPAREN_
    EX_co_LPAREN_e_RPAREN_
    EX_co2_LPAREN_e_RPAREN_
    EX_creat_LPAREN_e_RPAREN_
    EX_crn_LPAREN_e_RPAREN_
    EX_crvnc_LPAREN_e_RPAREN_
    EX_csn_LPAREN_e_RPAREN_
    EX_cys_L_LPAREN_e_RPAREN_
    EX_cytd_LPAREN_e_RPAREN_
    EX_dad_2_LPAREN_e_RPAREN_
    EX_dag_hs_LPAREN_e_RPAREN_
    EX_dcyt_LPAREN_e_RPAREN_
    EX_dgsn_LPAREN_e_RPAREN_
    EX_dhdascb_LPAREN_e_RPAREN_
    EX_din_LPAREN_e_RPAREN_
    EX_dlnlcg_LPAREN_e_RPAREN_
    EX_dopa_LPAREN_e_RPAREN_
    EX_drib_LPAREN_e_RPAREN_
    EX_eicostet_LPAREN_e_RPAREN_
    EX_elaid_LPAREN_e_RPAREN_
    EX_etoh_LPAREN_e_RPAREN_
    EX_fe2_LPAREN_e_RPAREN_
    EX_fe3_LPAREN_e_RPAREN_
    EX_fol_LPAREN_e_RPAREN_
    EX_for_LPAREN_e_RPAREN_
    EX_fru_LPAREN_e_RPAREN_
    EX_fuc14galacglcgalgluside_hs_LPAREN_e_RPAREN_
    EX_fuc_L_LPAREN_e_RPAREN_
    EX_gal_LPAREN_e_RPAREN_
    EX_gam_LPAREN_e_RPAREN_
    EX_gchola_LPAREN_e_RPAREN_
    EX_glc_LPAREN_e_RPAREN_
    EX_gln_L_LPAREN_e_RPAREN_
    EX_gluala_LPAREN_e_RPAREN_
    EX_glu_L_LPAREN_e_RPAREN_
    EX_gly_LPAREN_e_RPAREN_
    EX_glyb_LPAREN_e_RPAREN_
    EX_glyc_LPAREN_e_RPAREN_
    EX_glyc_S_LPAREN_e_RPAREN_
    EX_glygn2_LPAREN_e_RPAREN_
    EX_gsn_LPAREN_e_RPAREN_
    EX_gthox_LPAREN_e_RPAREN_
    EX_gthrd_LPAREN_e_RPAREN_
    EX_gua_LPAREN_e_RPAREN_
    EX_h_LPAREN_e_RPAREN_
    EX_h2o_LPAREN_e_RPAREN_
    EX_ha_LPAREN_e_RPAREN_
    EX_hdca_LPAREN_e_RPAREN_
    EX_hdcea_LPAREN_e_RPAREN_
    EX_his_L_LPAREN_e_RPAREN_
    EX_hxan_LPAREN_e_RPAREN_
    EX_ile_L_LPAREN_e_RPAREN_
    EX_inost_LPAREN_e_RPAREN_
    EX_k_LPAREN_e_RPAREN_
    EX_lac_D_LPAREN_e_RPAREN_
    EX_lac_L_LPAREN_e_RPAREN_
    EX_lcts_LPAREN_e_RPAREN_
    EX_leuktrA4_LPAREN_e_RPAREN_
    EX_leuktrB4_LPAREN_e_RPAREN_
    EX_leuktrD4_LPAREN_e_RPAREN_
    EX_leuktrE4_LPAREN_e_RPAREN_
    EX_leu_L_LPAREN_e_RPAREN_
    EX_lnlc_LPAREN_e_RPAREN_
    EX_lnlnca_LPAREN_e_RPAREN_
    EX_lnlncg_LPAREN_e_RPAREN_
    EX_lpchol_hs_LPAREN_e_RPAREN_
    EX_lys_L_LPAREN_e_RPAREN_
    EX_mag_hs_LPAREN_e_RPAREN_
    EX_malt_LPAREN_e_RPAREN_
    EX_man_LPAREN_e_RPAREN_
    EX_meoh_LPAREN_e_RPAREN_
    EX_met_L_LPAREN_e_RPAREN_
    EX_n2m2nmasn_LPAREN_e_RPAREN_
    EX_na1_LPAREN_e_RPAREN_
    EX_nac_LPAREN_e_RPAREN_
    EX_ncam_LPAREN_e_RPAREN_
    EX_nh4_LPAREN_e_RPAREN_
    EX_no_LPAREN_e_RPAREN_
    EX_nrvnc_LPAREN_e_RPAREN_
    EX_o2_LPAREN_e_RPAREN_
    EX_o2s_LPAREN_e_RPAREN_
    EX_ocdca_LPAREN_e_RPAREN_
    EX_ocdcea_LPAREN_e_RPAREN_
    EX_octa_LPAREN_e_RPAREN_
    EX_orn_LPAREN_e_RPAREN_
    EX_oxa_LPAREN_e_RPAREN_
    EX_pchol_hs_LPAREN_e_RPAREN_
    EX_pe_hs_LPAREN_e_RPAREN_
    EX_pglyc_hs_LPAREN_e_RPAREN_
    EX_pheacgln_LPAREN_e_RPAREN_
    EX_phe_L_LPAREN_e_RPAREN_
    EX_pheme_LPAREN_e_RPAREN_
    EX_pi_LPAREN_e_RPAREN_
    EX_pnto_R_LPAREN_e_RPAREN_
    EX_ppa_LPAREN_e_RPAREN_
    EX_pro_L_LPAREN_e_RPAREN_
    EX_ps_hs_LPAREN_e_RPAREN_
    EX_pydam_LPAREN_e_RPAREN_
    EX_pydx_LPAREN_e_RPAREN_
    EX_pydxn_LPAREN_e_RPAREN_
    EX_pyr_LPAREN_e_RPAREN_
    EX_retfa_LPAREN_e_RPAREN_
    EX_retinol_LPAREN_e_RPAREN_
    EX_retn_LPAREN_e_RPAREN_
    EX_rib_D_LPAREN_e_RPAREN_
    EX_ribflv_LPAREN_e_RPAREN_
    EX_Rtotal_LPAREN_e_RPAREN_
    EX_sel_LPAREN_e_RPAREN_
    EX_ser_L_LPAREN_e_RPAREN_
    EX_so4_LPAREN_e_RPAREN_
    EX_sph1p_LPAREN_e_RPAREN_
    EX_strch1_LPAREN_e_RPAREN_
    EX_strch2_LPAREN_e_RPAREN_
    EX_sucr_LPAREN_e_RPAREN_
    EX_tag_hs_LPAREN_e_RPAREN_
    EX_taur_LPAREN_e_RPAREN_
    EX_thm_LPAREN_e_RPAREN_
    EX_thr_L_LPAREN_e_RPAREN_
    EX_thymd_LPAREN_e_RPAREN_
    EX_tmndnc_LPAREN_e_RPAREN_
    EX_tre_LPAREN_e_RPAREN_
    EX_trp_L_LPAREN_e_RPAREN_
    EX_ttdca_LPAREN_e_RPAREN_
    EX_tyr_L_LPAREN_e_RPAREN_
    EX_ura_LPAREN_e_RPAREN_
    EX_urea_LPAREN_e_RPAREN_
    EX_uri_LPAREN_e_RPAREN_
    EX_utp_LPAREN_e_RPAREN_
    EX_vacc_LPAREN_e_RPAREN_
    EX_val_L_LPAREN_e_RPAREN_
    EX_xyl_D_LPAREN_e_RPAREN_
    EX_xylt_LPAREN_e_RPAREN_
    EX_acmana_LPAREN_e_RPAREN_
    EX_fald_LPAREN_e_RPAREN_
    EX_HC00250_LPAREN_e_RPAREN_
    EX_HC01609_LPAREN_e_RPAREN_
    EX_HC01610_LPAREN_e_RPAREN_
    EX_prpp_LPAREN_e_RPAREN_
    EX_ptrc_LPAREN_e_RPAREN_
    EX_pydx5p_LPAREN_e_RPAREN_
    EX_spmd_LPAREN_e_RPAREN_
    EX_no2_LPAREN_e_RPAREN_
    EX_prostgh2_LPAREN_e_RPAREN_
    EX_ppi_LPAREN_e_RPAREN_
    EX_ddca_LPAREN_e_RPAREN_
    EX_glcur_LPAREN_e_RPAREN_
    EX_3bcrn_
    EX_3ddcrn_
    EX_3deccrn_
    EX_3hdececrn_
    EX_3hexdcrn_
    EX_3ivcrn_
    EX_3octdec2crn_
    EX_3octdeccrn_
    EX_3octdece1crn_
    EX_3tdcrn_
    EX_3tetd7ecoacrn_
    EX_3thexddcoacrn_
    EX_3ttetddcoacrn_
    EX_c6crn_
    EX_c8crn_
    EX_ivcrn_
    EX_adpcbl_LPAREN_e_RPAREN_
    EX_carn_LPAREN_e_RPAREN_
    EX_sbt_DASH_d_LPAREN_e_RPAREN_
    EX_fmn_LPAREN_e_RPAREN_
    EX_q10h2_LPAREN_e_RPAREN_
    EX_34hpp_
    EX_3mob_LPAREN_e_RPAREN_
    EX_4mop_LPAREN_e_RPAREN_
    EX_5mta_LPAREN_e_RPAREN_
    EX_5oxpro_LPAREN_e_RPAREN_
    EX_ahcys_LPAREN_e_RPAREN_
    EX_aicar_LPAREN_e_RPAREN_
    EX_mal_L_LPAREN_e_RPAREN_
    EX_fum_LPAREN_e_RPAREN_
    EX_xmp_LPAREN_e_RPAREN_
    EX_dcmp_LPAREN_e_RPAREN_
    EX_glyc3p_LPAREN_e_RPAREN_
    biomass_reaction
    biomass_protein
    biomass_DNA
    biomass_RNA
    biomass_carbohydrate
    biomass_lipid
    biomass_other
    BETBGTtc
    CRNATBtc
    CYTDt2r
    FOLABCCte
    H2OGLYAQPt
    HDCAFAPMtc
    LGNCFATtc
    NACSMCTte
    OCDCAFAPMtc
    OCDCAFATPtc
    PCHOLABCtc
    PHEMEABCte
    PMTFATPtc
    PSHSABCtc
    THMDt2r
    URIt2r
    3HPVSTETCOAhcm
    3HPVSTETteb
    3HPVSteb
    3HPVStep
    AM19CSALThr
    AM19CSteb
    AM1CSAhr
    AM1CSAtep
    AM9CSAhr
    AM9CSAteb
    AM9CSAtep
    Am19CStev
    Am1CSAteb
    CRVSM1hr
    CRVSM23hr
    CSAtd
    CVM1GLUChc
    CVM23GLUChc
    EX_3hpvs_LPAREN_e_RPAREN_
    EX_3hpvstet_LPAREN_e_RPAREN_
    EX_am19cs_LPAREN_e_RPAREN_
    EX_am1csa_LPAREN_e_RPAREN_
    EX_am9csa_LPAREN_e_RPAREN_
    EX_csa_LPAREN_e_RPAREN_
    EX_ptvstlac_LPAREN_e_RPAREN_
    EX_ptvstm3_LPAREN_e_RPAREN_
    PTVSTLACtev
    PTVSTM3eb
    PTVSTM3hc
    PTVSThc
    SMVACIDhep
    SMVGLUChep
    3HPVSCOAitm
    3HPVSTETCOAitm
    AM19CSitr
    CSAitr
    PTVSTLACitr
    PTVSTM3itr
    3HPVSCOAhc
    3HPVSTEThc
    AM1CSAitr
    AM9CSAitr
    AIRCr_PRASCS
    FAOXC10080m
    FAOXC101_3Em
    FAOXC101_4Em
    FAOXC101_4Zm
    FAOXC101_4Zx
    FAOXC102_4Z_7Zm
    FAOXC120100m
    FAOXC121_3Em
    FAOXC121_3Zm
    FAOXC121_5Em
    FAOXC122_3E_6Em
    FAOXC122_3Z_6Zm
    FAOXC122_3Z_6Zx
    FAOXC123_3Z_6Z_9Zm
    FAOXC140120m
    FAOXC141_5Em
    FAOXC141_5Zm
    FAOXC141_7Em
    FAOXC142_5E_8Em
    FAOXC142_5Z_8Zm
    FAOXC142_5Z_8Zx
    FAOXC143_5Z_8Z_11Zm
    FAOXC160140m
    FAOXC161_7Em
    FAOXC161_7Zm
    FAOXC161_9Em
    FAOXC161140m
    FAOXC162_7E_10Em
    FAOXC162_7Z_10Zm
    FAOXC163_4Z_7Z_10Zm
    FAOXC163_4Z_7Z_10Zx
    FAOXC163_7Z_10Z_13Zm
    FAOXC164_4Z_7Z_10Z_13Zm
    FAOXC181_11Em
    FAOXC181_9Em
    FAOXC181_9Zm
    FAOXC182_9E_12Em
    FAOXC182_9Z_12Zm
    FAOXC183_6Z_9Z_12Zm
    FAOXC183_6Z_9Z_12Zx
    FAOXC183_9Z_12Z_15Zm
    FAOXC184_3Z_6Z_9Z_12Zm
    FAOXC184_3Z_6Z_9Z_12Zx
    FAOXC184_6Z_9Z_12Z_15Zm
    FAOXC185_3Z_6Z_9Z_12Z_15Zm
    FAOXC204_5Z_8Z_11Z_14Zm
    FAOXC204_5Z_8Z_11Z_14Zx
    FAOXC204184m2
    FAOXC205_5Z_8Z_11Z_14Z_17Zm
    FAOXC225_4Z_7Z_10Z_13Z_16Zx
    FAOXC4020m
    FAOXC6040m
    FAOXC61_3Zm
    FAOXC8060m
    FAOXC81_5Zm



```python
for rx in model_MeanCancerBiopsy.reactions:
    print(rx.id)
```

    13DAMPPOX
    1PPDCRp
    2HBO
    2OXOADOXm
    34DHOXPEGOX
    34HPLFM
    3DPHBH1
    3DPHBH2
    3DSPHR
    3HAO
    3HBCDm
    3HBCOAHLm
    3SALATAi
    3SALATAim
    3SPYRSP
    3SPYRSPm
    41R1H2MAE12BOOX
    41R2A1H12BOOX
    4HBZCOAFm
    4HBZFm
    4HOXPACDOX_NADP_
    5ADTSTSTERONESULT
    A_MANASEly
    AACOAT
    AATAi
    ABO8g
    ABTArm
    ABUTD
    ACACT10m
    ACACT1r
    ACACT1rm
    ACACT1x
    ACACT4p
    ACACT5p
    ACACT6p
    ACACT7p
    ACCOAC
    ACCOALm
    ACGAM6PSi
    ACGAMK
    ACGAMPM
    ACITL
    ACNAM9PL
    ACNAMPH
    ACNMLr
    ACOAD10m
    ACOAD1fm
    ACOAD9m
    ACONTm
    ACP1_FMN_
    ACS
    ACS2
    ACSm
    ADK1
    ADK1m
    ADK3m
    ADMDC
    ADNK1
    ADNK1m
    ADPRDP
    ADPT
    ADRNCPT1
    ADRNCPT2
    ADSK
    ADSL1
    ADSL2
    ADSS
    AG13T10g
    AG13T11g
    AG13T12g
    AG13T13g
    AG13T14g
    AG13T15g
    AG13T16g
    AG13T17g
    AG13T18g
    AG13T1g
    AG13T2g
    AG13T3g
    AG13T4g
    AG13T5g
    AG13T6g
    AG13T7g
    AG13T8g
    AG13T9g
    AGPAT1
    AGTim
    AHC
    AHEXASE2ly
    AHEXASEly
    AICART
    AKR1C1
    AKR1C42
    AKR1D2
    ALASm
    ALCD1
    ALCD21_L
    ALCD22_L
    ALCD2if
    ALCD2yf
    ALDD21
    ALDD2x
    ALDD2xm
    ALDD2y
    ALR2
    ALR3
    AMACR2p
    AMANK
    ARABR
    ARACHCPT1
    ARACHCPT2
    ARGSL
    ARGSS
    R_group_phosphotase_3
    ARTPLM3
    ARTPLM3m
    ASAH1
    ASNNm
    ASPCTr
    ASPTA
    ASPTAm
    B_MANNASEly
    B3GALTg
    B3GNT11g
    B3GNT12g
    B3GNT310g
    B3GNT312g
    B3GNT315g
    B3GNT51g
    BAAT2x
    BAMPPALDOX
    BDHm
    BDMT_U
    BETALDHxm
    BILIRED
    BPNT
    BPNT2
    C14STRr
    C160CPT1
    C160CPT2
    C161CPT1
    C161CPT2
    C180CPT1
    C180CPT2
    C181CPT1
    C181CPT2
    C204CPT1
    C204CPT2
    C226CPT1
    C226CPT2
    C2M26DCOAHLm
    C2M26DCOAHLx
    C3STDH1Pr
    C3STDH1r
    C3STKR2r
    C4STMO1r
    C4STMO2r
    CAT2p
    CATm
    CATp
    CBPS
    CDIPTr
    CDS
    CDSm
    CEPTC
    CH25H
    CHLPCTD
    CHOLD2m
    CHOLK
    CKc
    CLS_hs
    CMPSAS
    COQ3m
    COQ5m
    COQ6m
    COQ7m
    CORE2GTg
    CORE3GTg
    CORE4GTg
    CORE6GTg
    COUCOAFm
    CPPPGO
    CRTNsyn
    CSm
    CTPS1
    CTPS2
    CYSGLTH
    CYSO
    CYSTA
    CYSTAm
    CYTD
    CYTDK1
    CYTDK2m
    CYTK1
    CYTK10
    CYTK10n
    CYTK11
    CYTK11n
    CYTK12
    CYTK12n
    CYTK13
    CYTK13n
    CYTK14
    CYTK14n
    CYTK1m
    CYTK1n
    CYTK2
    CYTK2n
    CYTK3
    CYTK3n
    CYTK4
    CYTK4n
    CYTK5
    CYTK5n
    CYTK6
    CYTK6n
    CYTK7
    CYTK7n
    CYTK8
    CYTK8n
    CYTK9
    CYTK9n
    DAGK_hs
    DASCBR
    DCIm
    DCMPDA
    DCSPTN1CPT1
    DCSPTN1CPT2
    DDPGAm
    DESAT16_2
    DESAT18_10
    DESAT18_3
    DESAT18_4
    DESAT18_5
    DESAT18_6
    DESAT18_7
    DESAT18_8
    DESAT18_9
    DESAT20_2
    DESAT22_1p
    DGAT
    DGK1
    DGNSKm
    DGULND
    DHCR241r
    DHCR242r
    DHCR243r
    DHCR71r
    DHCR72r
    DHCRD1
    DHCRD2
    DHDPBMTm
    DHEASULT
    DHFR
    DHORTS
    DHPM1
    DLNLCGCPT1
    DLNLCGCPT2
    DM_atp_c_
    DMATT
    DMATTx
    DOLASNT_Uer
    DOLDPP_Uer
    DOLPGT1_Uer
    DOLPGT2_Uer
    DOLPGT3_Uer
    DOLPMT_U
    DOLPMT1_Uer
    DOLPMT2_Uer
    DOLPMT3_Uer
    DOLPMT4_Uer
    DOPABMO
    DPCOAK
    DPGase
    DPGM
    DPHMBDCm
    DPMVDx
    DPPS
    DRPA
    DSAT
    DTMPK
    DURAD
    DURAD2
    DURIK1
    DURIPP
    DUTPDPm
    DUTPDPn
    EBP1r
    EBP2r
    ECOAH12m
    ECOAH1m
    ECOAH1x
    ECOAH9m
    EHGLAT2m
    EHGLATm
    EICOSTETCPT1
    EICOSTETCPT2
    ELAIDCPT1
    ELAIDCPT2
    ENGASE2ly
    ENGASE3ly
    ENGASEly
    ENO
    ETF
    ETFQO
    ETHAK
    ETHP
    F1PGT
    F6Tg
    FACOAL150
    FACOAL160i
    FACOAL161
    FACOAL170
    FACOAL180i
    FACOAL1812
    FACOAL1813
    FACOAL181i
    FACOAL1821
    FACOAL1822
    FACOAL1831
    FACOAL1832
    FACOAL184
    FACOAL191
    FACOAL203
    FACOAL204
    FACOAL2042
    FACOAL205
    FACOAL224
    FACOAL2251
    FACOAL2252
    FACOAL226
    FACOAL240
    FACOAL241
    FACOAL244_1
    FACOAL245_1
    FACOAL245_2
    FACOAL246_1
    FACOAL260
    FAEL183
    FAEL184
    FAEL204
    FAEL205
    FALDH
    FAOXC1811601m
    FAOXC182806m
    FAOXC183806m
    FAOXC200180m
    FAOXC2031836m
    FAOXC2051843m
    FAOXC2242046m
    FAOXC2251836m
    FAOXC2252053m
    FAOXC226205m
    FAS80COA_L
    FBA
    FBA2
    FBA4
    FBP26
    FCLTm
    FE3R2e
    FK
    FKYNH
    FPGS
    FPGS2
    FPGS2m
    FPGS3
    FPGS3m
    FPGS4
    FPGS4m
    FPGS5
    FPGS5m
    FPGS6
    FPGS6m
    FPGS7
    FPGS7m
    FPGS8
    FPGS8m
    FPGS9
    FPGS9m
    FPGSm
    FTHFDH
    FTHFL
    FTHFLm
    FUCASE2ly
    FUCASEe
    FUCASEly
    FUM
    FUMm
    FUT16g
    FUT31g
    FUT910g
    FUT911g
    FUT98g
    FUT99g
    G12MT1_U
    G12MT2_U
    G13MT_U
    G14T10g
    G14T11g
    G14T12g
    G14T13g
    G14T14g
    G14T15g
    G14T16g
    G14T17g
    G14T18g
    G14T19g
    G14T20g
    G14T21g
    G14T2g
    G14T3g
    G14T4g
    G14T5g
    G14T6g
    G14T7g
    G14T8g
    G14T9g
    G14Tg
    G16MT_U
    G5SADrm
    G5SDym
    G6PDA
    G6PDH1rer
    G6PDH2r
    G6PDH2rer
    GALASE10ly
    GALASE11ly
    GALASE12ly
    GALASE13ly
    GALASE14ly
    GALASE15ly
    GALASE16ly
    GALASE17ly
    GALASE18ly
    GALASE19ly
    GALASE1ly
    GALASE20ly
    GALASE3ly
    GALASE4ly
    GALASE5ly
    GALASE6ly
    GALASE7ly
    GALASE8ly
    GALASE9ly
    GALC
    GALNTg
    GALOR
    GALU
    GAPD
    GARFT
    GASNASE2ly
    GASNASE3ly
    GASNASEly
    GBA
    GCALDD
    GCC2cm
    GF6PTA
    GFUCS
    GGH_10FTHF5GLUl
    GGH_10FTHF6GLUl
    GGH_10FTHF7GLUl
    GGH_5DHFl
    GGH_5THFl
    GGH_6DHFl
    GGH_6THFl
    GGH_7DHFl
    GGH_7THFl
    GGNG
    GGT5r
    GHMT2rm
    GK1
    GK1m
    GLAl
    GLBRAN
    GLCAASE8ly
    GLCAASE9ly
    GLCNACPT_U
    GLCNACT_U
    GLDBRAN
    GLGNS1
    GLNS
    GLPASE1
    GLPASE2
    GLU5Km
    GLUCYS
    GLUDxm
    GLUDym
    GLUPRT
    GLXO1
    GLYAMDTRc
    GLYCLTDy
    GLYCTO1p
    GLYK
    GLYOp
    GMAND
    GMPR
    GMPS2
    GNDer
    GPAM_hs
    GRTT
    GRTTx
    GTHDH
    GTHO
    GTHOm
    GTHP
    GTHPe
    GTHPm
    GTHS
    GTPCI
    GTPCIn
    GUAPRT
    GULN3D
    GULNDer
    H2CO3D
    H2O2syn
    HACD1m
    HACD1x
    HACD9m
    HBZOPT10m
    HEX1
    HEX10
    HEX4
    HEX7
    HIBDm
    HISDC
    HKYNH
    HMBS
    HMGCOASi
    HOXG
    HPYRDC
    HPYRR2x
    HPYRRy
    HSD17B2r
    HSD17B42x
    HSD17B4x
    HSD3A1r
    HSD3A2r
    HSD3B11
    HSD3B11r
    HSD3B12r
    HSD3B13r
    HXPRT
    HYPOE
    HYPTROX
    ICDHy
    ICDHyp
    ICDHyrm
    ILETA
    IMPC
    IMPD
    IPDDIx
    ITCOAL1m
    KHK2
    KYN
    KYN3OX
    LALDO
    LALDO2x
    LCADi
    LCADi_D
    LCAT1e
    LCYSTAT
    LCYSTATm
    LDH_L
    LDH_Lm
    LEUTA
    LFORKYNHYD
    LGTHL
    LNELDCCPT1
    LNELDCCPT2
    LNLCCPT1
    LNLCCPT2
    LNLNCACPT1
    LNLNCACPT2
    LNLNCGCPT1
    LNLNCGCPT2
    LNSTLSr
    LPCOXp
    LSTO1r
    LSTO2r
    LTA4H
    LTC4Sr
    LYSOXp
    M1316Mg
    M13N2Tg
    M13N4Tg
    M14NTg
    M16N4Tg
    M16N6Tg
    M16NTg
    MACOXO
    MAN1PT2
    MAOX
    MCCCrm
    MCLOR
    MDH
    MDHm
    ME1m
    ME2m
    METAT
    MEVK1x
    MG1er
    MG2er
    MG3er
    MGCHrm
    MGSA
    MGSA2
    MHISOR
    MI1345PP
    MI134PP
    MI145PK
    MI14PP
    MI1PP
    MI34PP
    MI3PP
    MI4PP
    MM5ag
    MM5bg
    MM5cg
    MM6B1ag
    MM6B1bg
    MM6B2g
    MM6bg
    MM7B1g
    MM7B2g
    MM7Cag
    MM7Cbg
    MM8Ag
    MM8Ber
    MM8Cg
    MMEm
    MMMm
    MMSAD3m
    MMTSADm
    MTHFC
    MTHFCm
    MTHFD
    MTHFD2
    MTHFD2m
    MTHFDm
    N3Tg
    N4Tg
    NACHEX10ly
    NACHEX11ly
    NACHEX12ly
    NACHEX13ly
    NACHEX14ly
    NACHEX15ly
    NACHEX16ly
    NACHEX17ly
    NACHEX18ly
    NACHEX19ly
    NACHEX20ly
    NACHEX21ly
    NACHEX22ly
    NACHEX23ly
    NACHEX24ly
    NACHEX25ly
    NACHEX26ly
    NACHEX27ly
    NACHEXA10ly
    NACHEXA11ly
    NACHEXA12ly
    NACHEXA13ly
    NACHEXA14ly
    NACHEXA15ly
    NACHEXA16ly
    NACHEXA17ly
    NACHEXA18ly
    NACHEXA19ly
    NACHEXA20ly
    NACHEXA21ly
    NACHEXA22ly
    NACHEXA9ly
    NADK
    NADN
    NADS2
    NAGA2ly
    NAGLCAly
    NBAHH_ir
    NDP6
    NDP7er
    NDP7ex
    NDP7g
    NDP8ex
    NDPK10
    NDPK10n
    NDPK3m
    NDPK6
    NDPK6m
    NDPK6n
    NDPK9
    NDPK9n
    NICRNS
    NMNATr
    NMNS
    NNATr
    NNMT
    NORANMT
    NOS1
    NP1
    NT5C
    NTD1
    NTD10
    NTD11
    NTD2
    NTD2e
    NTD3
    NTD4
    NTD5
    NTD6
    NTD7
    NTD8
    NTD9
    OBDHc
    OCOAT1m
    OIVD1m
    OIVD2m
    OIVD3m
    OMPDC
    OPAHir
    ORNDC
    ORNTArm
    ORPT
    P45011A1m
    P45017A1r
    P45017A2r
    P45017A3r
    P45017A4r
    P45027A13m
    P45027A1m
    P4508B11r
    P4508B13r
    P5CDm
    P5CRm
    PACCOAL
    PCHOLP_hs
    PCHOLPm_hs
    PCm
    PDXPP
    PEAMNO
    PETOHMm_hs
    PFK
    PFK26
    PGCD
    PGI
    PGK
    PGL
    PGLer
    PGM
    PGMT
    PHACCOAGLNAC
    PHCDm
    PHETA1
    PHETA1m
    PHETHPTOX2
    PHYCBOXL
    PI345P3P
    PI345P3Pn
    PI34P5K
    PI45P3K
    PI45P3Kn
    PI45P5P
    PI45P5Pn
    PI45PLC
    PI4P3K
    PI4PLC
    PI4PP
    PI5P3K
    PIK3
    PIPLC
    PLA2_2e
    PMANM
    PMEVKx
    PNP
    PNTEH
    PNTK
    PPA
    PPAm
    PPAP
    PPBNGS
    PPCOACm
    PPCOAOm
    PPD2CSPp
    PPM
    PPNCL3
    PPPGOm
    PRAGSr
    PRDX
    PRFGS
    PROD2m
    PRPNCOAHYDm
    PRPPS
    PSDm_hs
    PSERT
    PSP_L
    PSSA1_hs
    PSSA2_hs
    PTE3x
    PTE4x
    PTPAT
    PTRCOX1
    PUNP1
    PUNP2
    PUNP3
    PUNP4
    PUNP5
    PUNP6
    PUNP7
    PYDAMK
    PYDXK
    PYDXNK
    PYDXPP
    PYK
    PYNP2r
    QUILSYN
    RADH2
    RAI1
    RAI2
    RAI3
    RBFK
    RDH1
    RDH1a
    RDH2
    RDH2a
    RDH3
    RDH3a
    RETI1
    RETI2
    RNDR1
    RNDR2
    RNDR3
    RNDR4
    RTOTALCRNCPT1
    RTOTALCRNCPT2
    S23T2g
    S23T3g
    S23T4g
    S26Tg
    S6T10g
    S6T11g
    S6T12g
    S6T13g
    S6T14g
    S6T15g
    S6T16g
    S6T17g
    S6T18g
    S6T1g
    S6T2g
    S6T3g
    S6T4g
    S6T5g
    S6T6g
    S6T7g
    S6T8g
    S6T9g
    S6TASE10ly
    S6TASE11ly
    S6TASE12ly
    S6TASE13ly
    S6TASE14ly
    S6TASE15ly
    S6TASE16ly
    S6TASE17ly
    S6TASE18ly
    S6TASE19ly
    S6TASE1ly
    S6TASE20ly
    S6TASE21ly
    S6TASE22ly
    S6TASE23ly
    S6TASE24ly
    S6TASE25ly
    S6TASE26ly
    S6TASE2ly
    S6TASE3ly
    SADT
    SAMHISTA
    SBTD_D2
    SBTR
    SCP2x
    SCPx
    SERPT
    SFGTH
    SGPL12r
    SIAASE2ly
    SIAASE3ly
    SIAASE4ly
    SIAASEly
    SLCBK1
    SLDx
    SLDxm
    SMS
    SPHK21c
    SPMDOX
    SPODM
    SPODMm
    SPODMn
    SPODMx
    SPRMS
    SQLEr
    SQLSr
    SR5AR2r
    SR5ARr
    STS1
    SUCD1m
    SUCOAS1m
    T2M26DCOAHLm
    T2M26DCOAHLx
    T4HCINNMFM
    TALA
    THBPT4ACAMDASE
    THRD_L
    TKT1
    TKT2
    TMDK1
    TMDPP
    TMDS
    TPI
    TRDR
    TRDR2
    TRDR3
    TRPHYDRO2
    TRPO2
    TYRCBOX
    TYROXDAc
    TYRTA
    TYRTAm
    UAGDP
    UDPDOLPT_U
    UDPG4E
    UDPGD
    UGCG
    UGT1A10r
    UGT1A3r
    UGT1A4r
    UGT1A5r
    UGT1A5r2
    UGT1A9r
    UMPK
    UMPK2
    UMPK2n
    UMPK3
    UMPK3n
    UMPK4
    UMPK4n
    UMPK5
    UMPK5n
    UMPK6
    UMPK6n
    UMPK7
    UMPK7n
    UMPKn
    UPP3S
    UPPDC1
    UPPN
    URIDK2m
    URIK1
    VALTA
    VLCSp
    XAO2x
    XAOx
    XYLTD_Dr
    XYLUR
    r0009
    r0010
    r0021
    r0022
    r0024
    r0027
    r0033
    r0047
    r0051
    r0055
    r0074
    r0081
    r0082
    r0083
    r0084
    r0085
    r0086
    r0113
    r0119
    r0120
    r0121
    r0142
    r0156
    r0157
    r0160
    r0166
    r0170
    r0173
    r0178
    r0186
    r0193
    r0208
    r0210
    r0224
    r0239
    r0245
    r0246
    r0281
    r0283
    r0287
    r0301
    r0309
    r0310
    r0311
    r0317
    r0340
    r0345
    r0354
    r0355
    r0357
    r0358
    r0360
    r0361
    r0363
    r0364
    r0365
    r0377
    r0381
    r0383
    r0386
    r0391
    r0392
    r0393
    r0399
    r0403
    r0407
    r0408
    r0410
    r0422
    r0423
    r0424
    r0426
    r0430
    r0431
    r0441
    r0446
    r0450
    r0463
    r0464
    r0474
    r0475
    r0480
    r0510
    r0511
    r0527
    r0537
    r0539
    r0545
    r0547
    r0549
    r0552
    r0553
    r0555
    r0568
    r0570
    r0575
    r0578
    r0579
    r0580
    r0590
    r0591
    r0594
    r0595
    r0596
    r0604
    r0614
    r0618
    r0627
    r0629
    r0630
    r0634
    r0639
    r0642
    r0643
    r0653
    r0656
    r0660
    r0661
    r0666
    r0669
    r0670
    r0671
    r0672
    r0679
    r0680
    r0686
    r0688
    r0698
    r0706
    r0708
    r0709
    r0714
    r0715
    r0716
    r0717
    r0718
    r0719
    r0721
    r0722
    r0723
    r0724
    r0726
    r0727
    r0728
    r0729
    r0730
    r0731
    r0732
    r0733
    r0734
    r0739
    r0741
    r0743
    r0744
    r0750
    r0752
    r0753
    r0754
    r0755
    r0756
    r0757
    r0758
    r0774
    r0775
    r0776
    r0777
    r0778
    r0779
    r0781
    r0782
    r0783
    r0784
    r0786
    r0787
    r0788
    r0789
    r0795
    r1109
    r1135
    r1146
    r1164
    r1165
    r1166
    r1167
    r1168
    r1169
    r1170
    r1171
    r1172
    r1174
    r1175
    r1177
    r1179
    r1181
    r1182
    r1183
    r1251
    r1253
    r1254
    r1255
    r1257
    r1259
    r1260
    r1262
    r1378
    r1380
    r1382
    r1383
    r1384
    r1391
    r1392
    r1418
    r1443
    r1444
    r1445
    r1448
    r1457
    r1472
    r1474
    r1477
    r1481
    r1487
    r1488
    RE0344M
    RE0383C
    RE0453N
    RE0512M
    RE0512X
    RE0565C
    RE0566C
    RE0567C
    RE0568C
    RE0577M
    RE0578M
    RE0579C
    RE0579M
    RE0583C
    RE0688C
    RE0689E
    RE0827C
    RE1062M
    RE1233M
    RE1254C
    RE1447N
    RE1448N
    RE1514M
    RE1514X
    RE1516M
    RE1517M
    RE1518M
    RE1520M
    RE1521M
    RE1522M
    RE1523M
    RE1525M
    RE1525X
    RE1526M
    RE1526X
    RE1527M
    RE1527X
    RE1531M
    RE1531X
    RE1532M
    RE1532X
    RE1533M
    RE1533X
    RE1534M
    RE1534X
    RE1573M
    RE1573X
    RE1635M
    RE1635R
    RE1635X
    RE1691M
    RE1804M
    RE1807M
    RE1815M
    RE1815R
    RE1815X
    RE1816M
    RE1816R
    RE1816X
    RE1817M
    RE1817R
    RE1817X
    RE1818M
    RE1818R
    RE1818X
    RE1819C
    RE1834M
    RE1835M
    RE1860E
    RE1916X
    RE2051C
    RE2051G
    RE2051R
    RE2272L
    RE2296X
    RE2318M
    RE2318R
    RE2318X
    RE2319M
    RE2319R
    RE2319X
    RE2443M
    RE2453M
    RE2454M
    RE2514C
    RE2514L
    RE2522X
    RE2523X
    RE2524X
    RE2525X
    RE2541L
    RE2625C
    RE2626M
    RE2649C
    RE2649M
    RE2666C
    RE2666G
    RE2677C
    RE2680C
    RE2814R
    RE2908M
    RE2908X
    RE2909M
    RE2909X
    RE2910M
    RE2910X
    RE2912M
    RE2912X
    RE2913M
    RE2914M
    RE2915M
    RE2916M
    RE2916X
    RE2917M
    RE2917X
    RE2919M
    RE2920M
    RE2921M
    RE2954C
    RE2987X
    RE2988X
    RE2989X
    RE2991X
    RE2992M
    RE2995M
    RE2995X
    RE2999M
    RE3000M
    RE3001M
    RE3003M
    RE3004M
    RE3005M
    RE3006M
    RE3009C
    RE3010C
    RE3011R
    RE3012C
    RE3012M
    RE3012R
    RE3014R
    RE3015R
    RE3016R
    RE3074X
    RE3076X
    RE3082X
    RE3088X
    RE3090X
    RE3093X
    RE3097X
    RE3106C
    RE3110C
    RE3112C
    RE3113C
    RE3114R
    RE3119C
    RE3120C
    RE3121C
    RE3122C
    RE3132R
    RE3136C
    RE3139X
    RE3141X
    RE3146R
    RE3156X
    RE3158X
    RE3160R
    RE3177M
    RE3178M
    RE3179M
    RE3184M
    RE3185M
    RE3186M
    RE3189M
    RE3190M
    RE3191M
    RE3192M
    RE3193M
    RE3194M
    RE3218C
    RE3224C
    RE3225C
    RE3226C
    RE3227C
    RE3228C
    RE3229C
    RE3230C
    RE3231C
    RE3232C
    RE3234C
    RE3235C
    RE3236C
    RE3237C
    RE3241C
    RE3241R
    RE3242C
    RE3242R
    RE3243C
    RE3243R
    RE3244C
    RE3244R
    RE3245C
    RE3248M
    RE3248X
    RE3250M
    RE3250X
    RE3269C
    RE3270C
    RE3272N
    RE3307C
    RE3335M
    RE3335R
    RE3335X
    RE3336M
    RE3336X
    RE3337M
    RE3338C
    RE3338M
    RE3338X
    RE3339M
    RE3340C
    RE3340M
    RE3340X
    RE3342M
    RE3342X
    RE3343M
    RE3344M
    RE3345M
    RE3345R
    RE3345X
    RE3346C
    RE3346M
    RE3346R
    RE3347C
    RE3381L
    RE3383M
    RE3384M
    RE3386M
    RE3387M
    RE3388M
    RE3390M
    RE3391M
    RE3392M
    RE3393M
    RE3394M
    RE3396M
    RE3398M
    RE3399M
    RE3400M
    RE3401M
    RE3402M
    RE3403M
    RE3404M
    RE3432C
    RE3443M
    RE3443X
    RE3444M
    RE3446M
    RE3446R
    RE3446X
    RE3447M
    RE3448C
    RE3448M
    RE3448X
    RE3470C
    RE3476C
    RE3496N
    RE3502C
    RE3521M
    RE3521R
    RE3521X
    RE3532M
    RE3532R
    RE3533M
    RE3533R
    RE3534M
    RE3534R
    RE3554M
    RE3554R
    RE3557M
    RE3557R
    RE3559M
    RE3559X
    RE3560M
    RE3562M
    RE3562R
    RE3562X
    RE3563M
    RE3564C
    RE3564M
    RE3564X
    RE3572X
    RE3573X
    RE3575X
    RE3577X
    RE3587N
    RE3596C
    RE3597C
    RE3624M
    C10OHc
    C12OHc
    C141OHc
    C142OHc
    C14OHc
    C161OHc
    C162OHc
    C16OHc
    C181OHc
    C182OHc
    C18OHc
    C4OHc
    C50CPT1
    C60CPT1
    FAOXC101C102m
    FAOXC101C8m
    FAOXC101C8x
    FAOXC101m
    FAOXC102C101m
    FAOXC102C103m
    FAOXC102C81m
    FAOXC102m
    FAOXC103C102m
    FAOXC10DCC8DCx
    FAOXC121C101m
    FAOXC121C10x
    FAOXC122C101m
    FAOXC122m
    FAOXC123C102m
    FAOXC123m
    FAOXC12C12OHm
    FAOXC12DCC10DCx
    FAOXC141C121m
    FAOXC141C141OHm
    FAOXC142C142OHm
    FAOXC143C123m
    FAOXC14C14OHm
    FAOXC14DCC12DCx
    FAOXC15ATPx
    FAOXC161C141m
    FAOXC161C161OHm
    FAOXC162C142m
    FAOXC162C162OHm
    FAOXC163C143m
    FAOXC163C164Gm
    FAOXC163GC142m
    FAOXC163Gm
    FAOXC164C143m
    FAOXC164C165m
    FAOXC164GC163m
    FAOXC164m
    FAOXC165C164m
    FAOXC16BRx
    FAOXC16C16OHm
    FAOXC16DCC14DCx
    FAOXC16DCr
    FAOXC16OHC16r
    FAOXC181C161m
    FAOXC181C181OHm
    FAOXC182C162m
    FAOXC182C182OHm
    FAOXC183C163Gm
    FAOXC183C163m
    FAOXC184C163m
    FAOXC184C164m
    FAOXC184m
    FAOXC185C164m
    FAOXC185m
    FAOXC18C18OHm
    FAOXC204C184m
    FAOXC205C185m
    FAOXC225C204m
    FAOXC225C226m
    FAOXC225m
    FAOXC226C205m
    FAOXC226C225m
    FAOXC226C227m
    FAOXC226m
    FAOXC227C226m
    FAOXC5C5DCc
    FAOXC5C5OHm
    FAOXC5OHc
    FAOXC61m
    FAOXC6DCC4DCx
    FAOXC81C61m
    FAOXC8DCC6DCx
    FAOXOHC16C16DCc
    FAOXTC101TC102m
    FAOXTC102C101m
    FAOXTC122TC101m
    FAOXTC122m
    FAOXTC142TC122m
    FAOXTC162TC142m
    FAOXTC182TC162m
    OCD11COACPT1
    OCD11CRNCPT2
    SUCCOAPET
    ALAALACNc
    GLYGLYCNc
    GLYLEUHYDROc
    GLYSARCNc
    CHOLESTle
    DPMVDc
    DSREDUCr
    FADH2ETC
    FERO
    FMNALKPle
    GNDc
    HMGCOARc
    IPDDI
    MEVK1c
    PMEVKc
    NNDPR
    NADKm
    10FTHF5GLUtl
    10FTHF5GLUtm
    10FTHF6GLUtl
    10FTHF6GLUtm
    10FTHF7GLUtl
    10FTHF7GLUtm
    10FTHFtl
    10FTHFtm
    1MNCAMti
    2AMADPTm
    2HBt2
    2OXOADPTm
    34DHOXPEGt
    3MLDAt
    3MOBt2im
    3MOPt2im
    3SALAASPm
    4ABUTtm
    4MOPt2im
    4MPTNLte
    4MPTNLtm
    5ADTSTSTERONEGLCte
    5ADTSTSTERONEGLCtr
    5ADTSTSTERONEte
    5ADTSTSTERONEtr
    5AOPtm
    5DHFtl
    5THFtl
    5THFtm
    6DHFtl
    6DHFtm
    6THFtl
    6THFtm
    7DHFtl
    7DHFtm
    7THFtl
    7THFtm
    ABTti
    ACACt2
    ACACt2m
    ACALDtm
    ACALDtr
    ACALDtx
    ACCOAgt
    ACCOAtr
    ACGAMtly
    ACNAMlt
    ADNt4
    ADRNCOAtx
    ADRNCRNt
    AHCYStr
    ALAt4
    AMETr
    ANDRSTRNGLCte
    ANDRSTRNGLCtr
    ANDRSTRNte
    ANDRSTRNtr
    APRGSTRNte
    ARAB_Lt
    ARACHCRNt
    ARACHDCOAtx
    ARACHDtr
    ARGt4
    ARGtiDF
    ARTFR12
    ASCBt4
    ASNt4
    ASNtm
    Asn_X_Ser_Thrtr
    ATP2ter
    ATPasel
    ATPtm
    AVITE1t
    BHBt
    BHBtm
    BILDGLCURte
    BILGLCURt
    BILGLCURte
    BILGLCURtr
    BILIRUBtr
    C160CRNt
    C180CRNt
    C204CRNt
    C226CRNt
    CAATPS
    CDPDAGtm
    CERT1rt
    CERT2rt
    CHLtm
    CHOLtu
    CHSTEROLt
    CHSTEROLt2
    CITtam
    CLOXAtex2
    CLPNDt
    CMPACNAtg
    CO2ter
    CO2tm
    CO2tp
    COAtm
    COAtp
    COAtr
    COt
    CREATt4_2_r
    CYOOm2
    CYOR_u10m
    CYTDt4
    CYTDtm
    DAGt
    DCSPTN1COAtx
    DCSPTN1CRNt
    DECDPtm
    DGSNtm
    DHAAt1r
    DHCHOLESTANATEtm
    DHFtl
    DHFtm
    DHORD9
    DIDPtn
    DIGALSIDEtg
    DIGALSIDEtl
    DITPtn
    D_LACt2
    DLNLCGCRNt
    DOLGLCP_Lter
    DOLMANP_Lter
    DOLP_Uter
    DOPAtu
    DOPAVESSEC
    DTDPtn
    DTTPtn
    DUDPtn
    DUMPtn
    EICOSTETCRNt
    EICOSTETt
    ELAIDCRNt
    ETOHt
    ETOHtx
    F1Atg
    FADH2tru
    FADtru
    FATP6t
    FATP7t
    FATP9t
    FE2tm
    FORt2m
    FORtr
    FRDPtr
    FUC13GALACGLCGAL14ACGLCGALGLUSIDEte
    FUC13GALACGLCGAL14ACGLCGALGLUSIDEtg
    FUC14GALACGLCGALGLUSIDEte
    FUC14GALACGLCGALGLUSIDEtg
    FUCFUCFUCGALACGLC13GALACGLCGAL14ACGLCGALGLUSIDEte
    FUCFUCFUCGALACGLC13GALACGLCGAL14ACGLCGALGLUSIDEtg
    FUCtly
    FUMSO3tm
    FUMSO4tm
    G3PD2m
    G6Pter
    GALFUCGALACGLCGAL14ACGLCGALGLUSIDEte
    GALFUCGALACGLCGAL14ACGLCGALGLUSIDEtg
    GALGLUSIDEtg
    GALSIDEtg
    GALSIDEtl
    GALt1r
    GALtly
    GAMt1r
    GCHOLAt
    GDPFUCtg
    GDPtg
    GGLUCT
    GLCMter
    GLCt4
    GLCter
    GLCtg
    GLCURtly
    GLNt4
    GLNtm
    GLUt7l
    GLUtr
    GLXtm
    GLXtp
    GLYBt4_2_r
    GLYBtm
    GLYCLTtp
    GLYC_St
    GLYCtm
    GLYt4
    GLYtm
    GLYtp
    GMPtg
    GTHRDtr
    GULNter
    H2O2tm
    H2O2tn
    H2O2tp
    H2Oter
    H2Otg
    H2Otly
    H2Otn
    H2Otp
    HAS1
    HAS2
    HAtly
    HDCAter
    HDD2COAtx
    HISt4
    HMGCOAtm
    HMGCOAtx
    HOMt4
    HPACtr
    Htg
    Htm
    Htr
    Htx
    HXANtx
    IDPtn
    INSTt4
    IPDPtx
    ITPtn
    KCC2t
    KCCt
    KSII_CORE2t
    KSII_CORE2tly
    KSII_CORE4t
    KSII_CORE4tly
    KSIt
    KSItly
    LEUKTRA4t
    LEUKTRA4tr
    LEUKTRB4t
    LEUKTRD4tr
    L_LACtm
    LNELDCCRNt
    LNLCCRNt
    LNLNCACRNt
    LNLNCGCRNt
    LYSt4
    LYStiDF
    LYStip
    M4MPDOL_Uter
    M7MASNBterg
    M8MASNterg
    MAGt
    MALSO3tm
    MALSO4tm
    MANt1r
    MANter
    MANtg
    MANtly
    MEOHt2
    N2M2NMASNt
    N2M2NMASNtly
    NADH2_u10m
    NADHtpu
    NADHtru
    NADPHtru
    NADPHtxu
    NADPtru
    NADPtxu
    NADtru
    NaKt
    NAt
    NAtx
    NCKt
    NH4tp
    O2St
    O2Stm
    O2Stn
    O2Stx
    O2t
    O2ter
    O2tm
    O2tn
    O2tp
    OCCOAtm
    OCCOAtx
    ORNt3m
    ORNtiDF
    PAIL_HStn
    PAIL45P_HStn
    PAIL4P_HStn
    PAPStg
    PAPtg
    PE_HSter
    PE_HStm
    PEFLIP
    PEFLIPm
    PEt
    PGLYCt
    PHEACGLNt
    PHEMEtm
    PHEt4
    PIt2m
    PItg
    PItn
    PItx
    PNTOt5
    PPAt
    PPAtm
    PPItr
    PPItx
    PPPG9tm
    PRGNLONEtm
    PRGNLONEtr
    PRISTANALtx
    PRISTCOAtx
    PRISTtx
    PROtm
    PSFLIP
    PSFLIPm
    PSt3
    PYRt2m
    PYRt2p
    RETNGLCt
    RETNGLCt2
    RETNGLCt2r
    RETNGLCtr
    RETNt
    RETNtr
    RETNtr2
    RETt
    RTOTALCRNt
    Rtotaltl
    S2L2FN2M2MASNt
    S2L2FN2M2MASNtly
    S2L2N2M2MASNtly
    SCP21x
    SCP22x
    Ser_Thrtg
    SERt4
    SO4OXAtex2
    SO4t4_2
    SO4tl
    SPH1Pte
    SPHINGStl
    SPHS1Ptr
    STRDNCt
    TAURt4_2_r
    TAURtcx
    TCHOLAt3
    TCHOLAtx
    THFtl
    THFtm
    THMt3
    THP2Ctp
    THRt4
    TRPt4
    TSTSTERONEGLCte
    TSTSTERONEGLCtr
    TTDCAtr
    TYRt4
    UDPACGALtl
    UDPGALtg
    UDPGLCAter
    UDPGLCter
    UDPGLCtg
    UDPtl
    UGALNACtg
    UGLCNACtg
    URATEt
    URATEtx
    URIt4
    VACCt
    XOL27OHtm
    XOLDIOLONEt
    XYLTt
    r0002
    r0205
    r0801
    r0812
    r0829
    r0830
    r0838
    r0841
    r0845
    r0853
    r0860
    r0911
    r0913
    r0915
    r0926
    r0934
    r0936
    r0937
    r0941
    r0942
    r0950
    r0954
    r0960
    r0963
    r0968
    r0973
    r1001
    r1004
    r1007
    r1010
    r1020
    r1044
    r1078
    r1147
    r1148
    r1291
    r1292
    r1400
    r1421
    r1459
    r1464
    r1514
    r1515
    r1516
    r1517
    r1518
    r1519
    r1520
    r1521
    r1522
    r1523
    r1525
    r1529
    r1544
    r1546
    r1547
    r1548
    r1549
    r1551
    r1552
    r1553
    r1554
    r1556
    r1557
    r1559
    r1560
    r1561
    r1562
    r1563
    r1564
    r1565
    r1566
    r1567
    r1568
    r1569
    r1570
    r1571
    r1573
    r1574
    r1575
    r1576
    r1578
    r1579
    r1580
    r1581
    r1583
    r1584
    r1585
    r1586
    r1587
    r1588
    r1589
    r1590
    r1591
    r1592
    r1593
    r1594
    r1595
    r1596
    r1597
    r1598
    r1599
    r1600
    r1602
    r1603
    r1604
    r1605
    r1606
    r1607
    r1608
    r1609
    r1610
    r1611
    r1627
    r1628
    r1629
    r1630
    r1631
    r1633
    r1634
    r1635
    r1636
    r1637
    r1638
    r1639
    r1640
    r1641
    r1642
    r1643
    r1645
    r1646
    r1648
    r1649
    r1650
    r1651
    r1652
    r1653
    r1654
    r1655
    r1656
    r1657
    r1658
    r1659
    r1660
    r1661
    r1662
    r2073
    r2079
    r2080
    r2082
    r2083
    r2084
    r2085
    r2086
    r2087
    r2088
    r2089
    r2090
    r2092
    r2095
    r2096
    r2099
    r2102
    r2107
    r2110
    r2111
    r2114
    r2116
    r2117
    r2119
    r2121
    r2124
    r2128
    r2150
    r2313
    r2314
    r2316
    r2317
    r2320
    r2328
    r2344
    r2353
    r2355
    r2363
    r2373
    r2375
    r2377
    r2378
    r2379
    r2382
    r2386
    r2390
    r2449
    r2465
    r2473
    r2505
    r2506
    r2508
    r2509
    r2513
    r2514
    r2516
    r2517
    r2518
    r2521
    r2537
    CITt4_4
    PIt8
    2MB2COAc
    3HBCOARc
    BCRNe
    C141ACBP
    C141OHe
    C142ACBP
    C142OHe
    C162ACBP
    C162OHe
    C5DCe
    C6CRNe
    DDCRNe
    DECCRNe
    GLUTCOAACBP
    HDCACBP
    HDDACBP
    HDECAACBP
    HDECEACBP
    HEDCECRNe
    HEXCOAACBP
    HEXDCRNe
    HEXDIACTD
    HEXDICOAACBP
    HEXDICOAACBPx
    HIVCACBP
    HIVCRNe
    HOCDACBP
    HOCTDACBP
    HOCTDEC2CRNe
    HOCTDECCRNe
    HTDCACBP
    HTDCRNe
    IVCOAACBP
    IVCRNe
    OCD11CRNCACT
    OCTDEC2ACBP
    OCTDECE1CRNe
    OMHPALTD
    PRISTCOAtcx
    SCP21cx
    SUCCTD
    4OHPROIMINOtc
    ALAATB0tc
    ARGATB0tc
    ASNATB0tc
    BALABETAtc
    CARPEPT1tc
    CYSATB0tc
    GLNATB0tc
    ILEATB0tc
    LEUATB0tc
    LINOFATPtc
    LYSATB0tc
    METATB0tc
    NACHORCTL3le
    PHEATB0tc
    SERATB0tc
    TAUBETAtc
    THRATB0tc
    TRPATB0tc
    TYRATB0tc
    VALATB0tc
    VITEtl
    ADNt5le
    AHCYStd
    AMETtd
    ASPPROASCT1
    CHSTEROLtrc
    CYTDt5le
    DHAPtc
    FE2DMT1
    FE3MTP1
    FRDPtcr
    GLYC3Ptmc
    GSNt5le
    LACLt
    NADtm
    NADtx
    NODe
    PHEMEe
    Q10H2e
    RTOTALFATPc
    SPMTDe
    SPRMTDe
    URIt5le
    q10h2tc
    q10tm
    3MOBte
    4HPRO_LTte
    4MOPte
    5MTAte
    DM_Asn_X_Ser_Thr_ly_
    DM_datp_n_
    DM_dctp_n_
    DM_dgtp_n_
    DM_dttp_n_
    DM_ethamp_r_
    DM_n5m2masn_g_
    DM_sprm_c_
    DM_taur_LPAREN_c_RPAREN_
    DM_pe_hs_LPAREN_r_RPAREN_
    DM_4hrpo
    DM_Lcystin
    DM_anth
    EX_13_cis_retnglc_LPAREN_e_RPAREN_
    EX_1mncam_LPAREN_e_RPAREN_
    EX_2hb_LPAREN_e_RPAREN_
    EX_34dhoxpeg_LPAREN_e_RPAREN_
    EX_3mlda_LPAREN_e_RPAREN_
    EX_4hphac_LPAREN_e_RPAREN_
    EX_4mptnl_LPAREN_e_RPAREN_
    EX_5adtststerone_LPAREN_e_RPAREN_
    EX_5adtststeroneglc_LPAREN_e_RPAREN_
    EX_abt_LPAREN_e_RPAREN_
    EX_ac_LPAREN_e_RPAREN_
    EX_acac_LPAREN_e_RPAREN_
    EX_acald_LPAREN_e_RPAREN_
    EX_acetone_LPAREN_e_RPAREN_
    EX_acgam_LPAREN_e_RPAREN_
    EX_ade_LPAREN_e_RPAREN_
    EX_adn_LPAREN_e_RPAREN_
    EX_adrn_LPAREN_e_RPAREN_
    EX_akg_LPAREN_e_RPAREN_
    EX_ala_B_LPAREN_e_RPAREN_
    EX_ala_D_LPAREN_e_RPAREN_
    EX_ala_L_LPAREN_e_RPAREN_
    EX_amp_LPAREN_e_RPAREN_
    EX_andrstrn_LPAREN_e_RPAREN_
    EX_andrstrnglc_LPAREN_e_RPAREN_
    EX_aprgstrn_LPAREN_e_RPAREN_
    EX_arab_L_LPAREN_e_RPAREN_
    EX_arachd_LPAREN_e_RPAREN_
    EX_arg_L_LPAREN_e_RPAREN_
    EX_asn_L_LPAREN_e_RPAREN_
    EX_asp_L_LPAREN_e_RPAREN_
    EX_atp_LPAREN_e_RPAREN_
    EX_bhb_LPAREN_e_RPAREN_
    EX_bilglcur_LPAREN_e_RPAREN_
    EX_btn_LPAREN_e_RPAREN_
    EX_ca2_LPAREN_e_RPAREN_
    EX_caro_LPAREN_e_RPAREN_
    EX_cgly_LPAREN_e_RPAREN_
    EX_chol_LPAREN_e_RPAREN_
    EX_chsterol_LPAREN_e_RPAREN_
    EX_cit_LPAREN_e_RPAREN_
    EX_cl_LPAREN_e_RPAREN_
    EX_CLPND_LPAREN_e_RPAREN_
    EX_co_LPAREN_e_RPAREN_
    EX_co2_LPAREN_e_RPAREN_
    EX_creat_LPAREN_e_RPAREN_
    EX_crn_LPAREN_e_RPAREN_
    EX_crvnc_LPAREN_e_RPAREN_
    EX_csn_LPAREN_e_RPAREN_
    EX_cys_L_LPAREN_e_RPAREN_
    EX_cytd_LPAREN_e_RPAREN_
    EX_dad_2_LPAREN_e_RPAREN_
    EX_dag_hs_LPAREN_e_RPAREN_
    EX_dcyt_LPAREN_e_RPAREN_
    EX_dgsn_LPAREN_e_RPAREN_
    EX_dhdascb_LPAREN_e_RPAREN_
    EX_din_LPAREN_e_RPAREN_
    EX_dlnlcg_LPAREN_e_RPAREN_
    EX_dopa_LPAREN_e_RPAREN_
    EX_drib_LPAREN_e_RPAREN_
    EX_eicostet_LPAREN_e_RPAREN_
    EX_elaid_LPAREN_e_RPAREN_
    EX_etoh_LPAREN_e_RPAREN_
    EX_fe2_LPAREN_e_RPAREN_
    EX_fe3_LPAREN_e_RPAREN_
    EX_fol_LPAREN_e_RPAREN_
    EX_for_LPAREN_e_RPAREN_
    EX_fru_LPAREN_e_RPAREN_
    EX_fuc13galacglcgal14acglcgalgluside_hs_LPAREN_e_RPAREN_
    EX_fuc14galacglcgalgluside_hs_LPAREN_e_RPAREN_
    EX_fucfucfucgalacglc13galacglcgal14acglcgalgluside_hs_LPAREN_e_RPAREN_
    EX_fuc_L_LPAREN_e_RPAREN_
    EX_gal_LPAREN_e_RPAREN_
    EX_galfucgalacglcgal14acglcgalgluside_hs_LPAREN_e_RPAREN_
    EX_gam_LPAREN_e_RPAREN_
    EX_gchola_LPAREN_e_RPAREN_
    EX_glc_LPAREN_e_RPAREN_
    EX_gln_L_LPAREN_e_RPAREN_
    EX_gluala_LPAREN_e_RPAREN_
    EX_glu_L_LPAREN_e_RPAREN_
    EX_gly_LPAREN_e_RPAREN_
    EX_glyb_LPAREN_e_RPAREN_
    EX_glyc_LPAREN_e_RPAREN_
    EX_glyc_S_LPAREN_e_RPAREN_
    EX_glygn2_LPAREN_e_RPAREN_
    EX_gsn_LPAREN_e_RPAREN_
    EX_gthox_LPAREN_e_RPAREN_
    EX_gthrd_LPAREN_e_RPAREN_
    EX_gua_LPAREN_e_RPAREN_
    EX_h_LPAREN_e_RPAREN_
    EX_h2o_LPAREN_e_RPAREN_
    EX_ha_LPAREN_e_RPAREN_
    EX_ha_pre1_LPAREN_e_RPAREN_
    EX_hdca_LPAREN_e_RPAREN_
    EX_hdcea_LPAREN_e_RPAREN_
    EX_his_L_LPAREN_e_RPAREN_
    EX_hxan_LPAREN_e_RPAREN_
    EX_ile_L_LPAREN_e_RPAREN_
    EX_inost_LPAREN_e_RPAREN_
    EX_k_LPAREN_e_RPAREN_
    EX_lac_D_LPAREN_e_RPAREN_
    EX_lac_L_LPAREN_e_RPAREN_
    EX_lcts_LPAREN_e_RPAREN_
    EX_leuktrA4_LPAREN_e_RPAREN_
    EX_leuktrB4_LPAREN_e_RPAREN_
    EX_leuktrD4_LPAREN_e_RPAREN_
    EX_leuktrE4_LPAREN_e_RPAREN_
    EX_leu_L_LPAREN_e_RPAREN_
    EX_lnlc_LPAREN_e_RPAREN_
    EX_lnlnca_LPAREN_e_RPAREN_
    EX_lnlncg_LPAREN_e_RPAREN_
    EX_lpchol_hs_LPAREN_e_RPAREN_
    EX_lys_L_LPAREN_e_RPAREN_
    EX_mag_hs_LPAREN_e_RPAREN_
    EX_malt_LPAREN_e_RPAREN_
    EX_man_LPAREN_e_RPAREN_
    EX_meoh_LPAREN_e_RPAREN_
    EX_met_L_LPAREN_e_RPAREN_
    EX_n2m2nmasn_LPAREN_e_RPAREN_
    EX_na1_LPAREN_e_RPAREN_
    EX_nac_LPAREN_e_RPAREN_
    EX_ncam_LPAREN_e_RPAREN_
    EX_nh4_LPAREN_e_RPAREN_
    EX_no_LPAREN_e_RPAREN_
    EX_nrvnc_LPAREN_e_RPAREN_
    EX_o2_LPAREN_e_RPAREN_
    EX_o2s_LPAREN_e_RPAREN_
    EX_ocdca_LPAREN_e_RPAREN_
    EX_ocdcea_LPAREN_e_RPAREN_
    EX_octa_LPAREN_e_RPAREN_
    EX_orn_LPAREN_e_RPAREN_
    EX_oxa_LPAREN_e_RPAREN_
    EX_pchol_hs_LPAREN_e_RPAREN_
    EX_pe_hs_LPAREN_e_RPAREN_
    EX_pglyc_hs_LPAREN_e_RPAREN_
    EX_pheacgln_LPAREN_e_RPAREN_
    EX_phe_L_LPAREN_e_RPAREN_
    EX_pheme_LPAREN_e_RPAREN_
    EX_pi_LPAREN_e_RPAREN_
    EX_pnto_R_LPAREN_e_RPAREN_
    EX_ppa_LPAREN_e_RPAREN_
    EX_pro_L_LPAREN_e_RPAREN_
    EX_ps_hs_LPAREN_e_RPAREN_
    EX_pydam_LPAREN_e_RPAREN_
    EX_pydx_LPAREN_e_RPAREN_
    EX_pydxn_LPAREN_e_RPAREN_
    EX_pyr_LPAREN_e_RPAREN_
    EX_retfa_LPAREN_e_RPAREN_
    EX_retinol_LPAREN_e_RPAREN_
    EX_retn_LPAREN_e_RPAREN_
    EX_retnglc_LPAREN_e_RPAREN_
    EX_rib_D_LPAREN_e_RPAREN_
    EX_ribflv_LPAREN_e_RPAREN_
    EX_Rtotal_LPAREN_e_RPAREN_
    EX_sel_LPAREN_e_RPAREN_
    EX_ser_L_LPAREN_e_RPAREN_
    EX_so4_LPAREN_e_RPAREN_
    EX_sph1p_LPAREN_e_RPAREN_
    EX_strch1_LPAREN_e_RPAREN_
    EX_strch2_LPAREN_e_RPAREN_
    EX_sucr_LPAREN_e_RPAREN_
    EX_tag_hs_LPAREN_e_RPAREN_
    EX_taur_LPAREN_e_RPAREN_
    EX_tchola_LPAREN_e_RPAREN_
    EX_thm_LPAREN_e_RPAREN_
    EX_thr_L_LPAREN_e_RPAREN_
    EX_thymd_LPAREN_e_RPAREN_
    EX_tre_LPAREN_e_RPAREN_
    EX_trp_L_LPAREN_e_RPAREN_
    EX_tststeroneglc_LPAREN_e_RPAREN_
    EX_ttdca_LPAREN_e_RPAREN_
    EX_tyr_L_LPAREN_e_RPAREN_
    EX_udp_LPAREN_e_RPAREN_
    EX_ump_LPAREN_e_RPAREN_
    EX_ura_LPAREN_e_RPAREN_
    EX_urate_LPAREN_e_RPAREN_
    EX_urea_LPAREN_e_RPAREN_
    EX_uri_LPAREN_e_RPAREN_
    EX_utp_LPAREN_e_RPAREN_
    EX_vacc_LPAREN_e_RPAREN_
    EX_val_L_LPAREN_e_RPAREN_
    EX_xyl_D_LPAREN_e_RPAREN_
    EX_xylt_LPAREN_e_RPAREN_
    EX_4abutn_LPAREN_e_RPAREN_
    EX_acmana_LPAREN_e_RPAREN_
    EX_fald_LPAREN_e_RPAREN_
    EX_HC00250_LPAREN_e_RPAREN_
    EX_HC01609_LPAREN_e_RPAREN_
    EX_HC01610_LPAREN_e_RPAREN_
    EX_prpp_LPAREN_e_RPAREN_
    EX_ptrc_LPAREN_e_RPAREN_
    EX_pydx5p_LPAREN_e_RPAREN_
    EX_spmd_LPAREN_e_RPAREN_
    EX_no2_LPAREN_e_RPAREN_
    EX_sprm_LPAREN_e_RPAREN_
    EX_prostgh2_LPAREN_e_RPAREN_
    EX_ppi_LPAREN_e_RPAREN_
    EX_ddca_LPAREN_e_RPAREN_
    EX_CE1940_LPAREN_e_RPAREN_
    EX_glcur_LPAREN_e_RPAREN_
    EX_3bcrn_
    EX_3ddcrn_
    EX_3deccrn_
    EX_3hdececrn_
    EX_3hexdcrn_
    EX_3ivcrn_
    EX_3octdec2crn_
    EX_3octdeccrn_
    EX_3octdece1crn_
    EX_3tdcrn_
    EX_3tetd7ecoacrn_
    EX_3thexddcoacrn_
    EX_3ttetddcoacrn_
    EX_c5dc_
    EX_c6crn_
    EX_ivcrn_
    EX_adpcbl_LPAREN_e_RPAREN_
    EX_carn_LPAREN_e_RPAREN_
    EX_sbt_DASH_d_LPAREN_e_RPAREN_
    EX_fmn_LPAREN_e_RPAREN_
    EX_q10h2_LPAREN_e_RPAREN_
    EX_34hpp_
    EX_3mob_LPAREN_e_RPAREN_
    EX_4mop_LPAREN_e_RPAREN_
    EX_5mta_LPAREN_e_RPAREN_
    EX_mal_L_LPAREN_e_RPAREN_
    EX_fum_LPAREN_e_RPAREN_
    EX_xmp_LPAREN_e_RPAREN_
    EX_dcmp_LPAREN_e_RPAREN_
    EX_glyc3p_LPAREN_e_RPAREN_
    biomass_reaction
    biomass_protein
    biomass_DNA
    biomass_RNA
    biomass_carbohydrate
    biomass_lipid
    biomass_other
    ADNCNT3tc
    BETBGTtc
    CRNATBtc
    CYTDt2r
    FOLABCCte
    GSNt2r
    H2OGLYAQPt
    HDCAFAPMtc
    LGNCFATtc
    NACSMCTte
    OCDCAFAPMtc
    PCHOLABCtc
    PGLYCABCte
    PSHSABCtc
    THMDt2r
    URIt2r
    3HPVSTETCOAhcm
    3HPVSTETtev
    3HPVSthc
    CRVSM1hr
    CRVSM23hr
    CVM1GLUChc
    CVM23GLUChc
    EX_3hpvs_LPAREN_e_RPAREN_
    EX_3hpvstet_LPAREN_e_RPAREN_
    EX_am9csa_LPAREN_e_RPAREN_
    EX_csa_LPAREN_e_RPAREN_
    EX_ptvst_LPAREN_e_RPAREN_
    EX_ptvstlac_LPAREN_e_RPAREN_
    PTVSTLACtev
    PTVSThc
    PTVSTtu
    SMVACIDhep
    SMVGLUChep
    3HPVSCOAitm
    3HPVSTETCOAitm
    PTVSTLACitr
    3HPVSCOAhc
    3HPVSTEThc
    PTVSTitr
    AIRCr_PRASCS
    FAOXC10080m
    FAOXC101_3Em
    FAOXC101_4Em
    FAOXC101_4Zm
    FAOXC101_4Zx
    FAOXC102_4Z_7Zm
    FAOXC120100m
    FAOXC121_3Em
    FAOXC121_3Zm
    FAOXC121_5Em
    FAOXC122_3E_6Em
    FAOXC122_3Z_6Zm
    FAOXC122_3Z_6Zx
    FAOXC123_3Z_6Z_9Zm
    FAOXC140120m
    FAOXC141_5Em
    FAOXC141_5Zm
    FAOXC141_7Em
    FAOXC142_5E_8Em
    FAOXC142_5Z_8Zm
    FAOXC142_5Z_8Zx
    FAOXC143_5Z_8Z_11Zm
    FAOXC160140m
    FAOXC161_7Em
    FAOXC161_7Zm
    FAOXC161_9Em
    FAOXC161140m
    FAOXC162_7E_10Em
    FAOXC162_7Z_10Zm
    FAOXC163_4Z_7Z_10Zm
    FAOXC163_4Z_7Z_10Zx
    FAOXC163_7Z_10Z_13Zm
    FAOXC164_4Z_7Z_10Z_13Zm
    FAOXC181_11Em
    FAOXC181_9Em
    FAOXC181_9Zm
    FAOXC182_9E_12Em
    FAOXC182_9Z_12Zm
    FAOXC183_6Z_9Z_12Zm
    FAOXC183_9Z_12Z_15Zm
    FAOXC184_3Z_6Z_9Z_12Zm
    FAOXC184_3Z_6Z_9Z_12Zx
    FAOXC184_6Z_9Z_12Z_15Zm
    FAOXC185_3Z_6Z_9Z_12Z_15Zm
    FAOXC204_5Z_8Z_11Z_14Zm
    FAOXC204_5Z_8Z_11Z_14Zx
    FAOXC204184m2
    FAOXC205_5Z_8Z_11Z_14Z_17Zm
    FAOXC225_4Z_7Z_10Z_13Z_16Zm
    FAOXC4020m
    FAOXC6040m
    FAOXC61_3Zm
    FAOXC8060m
    FAOXC81_5Zm



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
R_Biopsies=get_drugable_targets(
    normal_Model=model_MeanNormalBiopsy, 
    disease_Model=model_MeanCancerBiopsy, model_name="Biopsies")
R_Biopsies["Dise_ratio"]=R_Biopsies["del_dise_flux"]/(R_Biopsies["dise_flux"])
R_Biopsies["Norm_ratio"]=R_Biopsies["del_norm_flux"]/(R_Biopsies["norm_flux"])
R_Biopsies[R_Biopsies["Dise_ratio"]< 0.999].sort_values(by="Dise_ratio", ascending=True)
```

    Common reactions size 1900
    Normal unique reactions size 336
    Disease unique reactions size 710


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
      <td>0.0219966</td>
      <td>0.033402</td>
      <td>HGNC:3350 or HGNC:3354 or HGNC:3354 or HGNC:3353</td>
      <td>Biopsies</td>
      <td>0.0219966</td>
      <td>0.239507</td>
      <td>1</td>
    </tr>
    <tr>
      <th>GAPD</th>
      <td>0.008</td>
      <td>0.0219966</td>
      <td>0.033402</td>
      <td>HGNC:24864 or HGNC:4141</td>
      <td>Biopsies</td>
      <td>0.0219966</td>
      <td>0.239507</td>
      <td>1</td>
    </tr>
    <tr>
      <th>PGK</th>
      <td>0.008</td>
      <td>0.0219966</td>
      <td>0.033402</td>
      <td>HGNC:8896 or HGNC:8898</td>
      <td>Biopsies</td>
      <td>0.0219966</td>
      <td>0.239507</td>
      <td>1</td>
    </tr>
    <tr>
      <th>PGM</th>
      <td>0.008</td>
      <td>0.0219966</td>
      <td>0.033402</td>
      <td>HGNC:1093 or HGNC:8888 or HGNC:8889</td>
      <td>Biopsies</td>
      <td>0.0219966</td>
      <td>0.239507</td>
      <td>1</td>
    </tr>
    <tr>
      <th>TPI</th>
      <td>0.008</td>
      <td>0.0219966</td>
      <td>0.033402</td>
      <td>HGNC:12009</td>
      <td>Biopsies</td>
      <td>0.0219966</td>
      <td>0.239507</td>
      <td>1</td>
    </tr>
    <tr>
      <th>biomass_reaction</th>
      <td>0.01</td>
      <td>0.01</td>
      <td>0.033402</td>
      <td></td>
      <td>Biopsies</td>
      <td>0.0219966</td>
      <td>0.299383</td>
      <td>0.454617</td>
    </tr>
    <tr>
      <th>biomass_protein</th>
      <td>0.0141643</td>
      <td>0.0141643</td>
      <td>0.033402</td>
      <td></td>
      <td>Biopsies</td>
      <td>0.0219966</td>
      <td>0.424056</td>
      <td>0.643933</td>
    </tr>
    <tr>
      <th>EX_glc_LPAREN_e_RPAREN_</th>
      <td>0.0156099</td>
      <td>0.0219966</td>
      <td>0.033402</td>
      <td></td>
      <td>Biopsies</td>
      <td>0.0219966</td>
      <td>0.467336</td>
      <td>1</td>
    </tr>
    <tr>
      <th>GLCt4</th>
      <td>0.0156099</td>
      <td>0.0219966</td>
      <td>0.033402</td>
      <td>HGNC:11037 or HGNC:11038 or HGNC:22146 or HGNC...</td>
      <td>Biopsies</td>
      <td>0.0219966</td>
      <td>0.467336</td>
      <td>1</td>
    </tr>
    <tr>
      <th>EX_lys_L_LPAREN_e_RPAREN_</th>
      <td>0.0168888</td>
      <td>0.0168888</td>
      <td>0.033402</td>
      <td></td>
      <td>Biopsies</td>
      <td>0.0219966</td>
      <td>0.505621</td>
      <td>0.767791</td>
    </tr>
    <tr>
      <th>EX_arg_L_LPAREN_e_RPAREN_</th>
      <td>0.027835</td>
      <td>0.0219966</td>
      <td>0.033402</td>
      <td></td>
      <td>Biopsies</td>
      <td>0.0219966</td>
      <td>0.833333</td>
      <td>1</td>
    </tr>
    <tr>
      <th>EX_thr_L_LPAREN_e_RPAREN_</th>
      <td>0.0319806</td>
      <td>0.0219966</td>
      <td>0.033402</td>
      <td></td>
      <td>Biopsies</td>
      <td>0.0219966</td>
      <td>0.957445</td>
      <td>1</td>
    </tr>
  </tbody>
</table>
</div>




```python
cobra.io.write_sbml_model(model_MeanCancerBiopsy, "Full_SCC_model_0419_n5.sbml")
cobra.io.write_sbml_model(model_MeanNormalBiopsy, "Full_Normal_model_0419_n5.sbml")
```


```python
FBA=model_MeanCancerBiopsy.optimize()
FBA.fluxes.to_csv("CancerFluxesFBA.csv")
FBA=model_MeanNormalBiopsy.optimize()
FBA.fluxes.to_csv("NormalFluxesFBA.csv")
```


```python
FBA=model_MeanCancerBiopsy.optimize()
FBA.metabolites.ser_L_e.summary()
```

    PRODUCING REACTIONS -- L-serine (ser_L_e)
    -----------------------------------------
    %       FLUX  RXN ID      REACTION
    ----  ------  ----------  ---------------------------------------
    100%   0.012  EX_ser_...  ser_L_e <=>
    
    CONSUMING REACTIONS -- L-serine (ser_L_e)
    -----------------------------------------
    %       FLUX  RXN ID      REACTION
    ----  ------  ----------  ---------------------------------------
    100%   0.012  r1559       ala_L_e + ser_L_c <=> ala_L_c + ser_L_e



```python
FBA.metabolites.ser_L_c.summary()
```

    PRODUCING REACTIONS -- L-serine (ser_L_c)
    -----------------------------------------
    %      FLUX  RXN ID      REACTION
    ---  ------  ----------  --------------------------------------------------
    51%  0.0486  PSP_L       h2o_c + pser_L_c --> pi_c + ser_L_c
    36%  0.0338  PSSA1_hs    pchol_hs_c + ser_L_c <=> chol_c + ps_hs_c
    13%  0.012   r1559       ala_L_e + ser_L_c <=> ala_L_c + ser_L_e
    
    CONSUMING REACTIONS -- L-serine (ser_L_c)
    -----------------------------------------
    %      FLUX  RXN ID      REACTION
    ---  ------  ----------  --------------------------------------------------
    63%  0.0593  r0160       pyr_c + ser_L_c --> ala_L_c + hpyr_c
    23%  0.022   PSSA2_hs    pe_hs_c + ser_L_c <=> etha_c + ps_hs_c
    14%  0.0131  biomass...  0.716189801699717 ala_L_c + 0.508866855524079 a...



```python
for rx in Recon2.metabolites.get_by_id("dhap_c").reactions:
    print(rx.id,rx.reaction)
```

    DHAPA Rtotalcoa_c + dhap_c --> adhap_hs_c + coa_c + 4.0 h_c
    MGSA dhap_c --> mthgxl_c + pi_c
    ALKP dhap_c + h2o_c --> dha_c + pi_c
    FBA5 tag1p_D_c <=> dhap_c + glyald_c
    TPI dhap_c <=> g3p_c
    G3PD1 glyc3p_c + nad_c <=> dhap_c + h_c + nadh_c
    FBA4 xu1p_D_c <=> dhap_c + gcald_c
    G3PD2m fad_m + glyc3p_c --> dhap_c + fadh2_m
    DHAPtc dhap_m --> dhap_c
    FBA2 f1p_c <=> dhap_c + glyald_c
    r0407 HC00361_c <=> dhap_c + e4p_c
    FBA fdp_c <=> dhap_c + g3p_c



```python
for rx in Recon2.metabolites.get_by_id("dhap_c").reactions:
    print(rx.id,rx.reaction)
```


```python
Recon2.reactions.r0160.name
```




    'L-Serine:pyruvate aminotransferase Glycine, serine and threonine metabolism EC:2.6.1.51'




```python
FBA.metabolites.gly_c.summary()
```

    PRODUCING REACTIONS -- glycine (gly_c)
    --------------------------------------
    %       FLUX  RXN ID      REACTION
    ----  ------  ----------  --------------------------------------------------
    100%  0.0119  r1546       gly_e + met_L_c <=> gly_c + met_L_e
    
    CONSUMING REACTIONS -- glycine (gly_c)
    --------------------------------------
    %       FLUX  RXN ID      REACTION
    ----  ------  ----------  --------------------------------------------------
    151%  0.018   biomass...  0.716189801699717 ala_L_c + 0.508866855524079 a...



```python

# Production of ATP from glucose in anaerobic conditions for Cancer Biopsies
# test for max ATP hydrolysis flux from only glucose

closedModel=model_MeanNormalBiopsy.copy()
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


```


```python
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
closedModel.reactions.EX_o2_LPAREN_e_RPAREN_.bounds=(-10,0)
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


```python
Recon2.metabolites.gln_L
```



Con la finalidad de probar la capacidad fisiologica del modelo se verificaron las posibles fuentes de obtención de energia. 
Incluir más contexto fisiologico

Incluir mayor contexto de normoxia e hipoxia

## AtpYieldCheck
FBA=cobra.flux_analysis.pfba(model_MeanCancerBiopsy)

ATP=0

for rx in model_MeanCancerBiopsy.metabolites.atp_c.reactions:
    if FBA[rx.id]>0:
        for prod in rx.products:
            if prod.id == 'atp_c' and len(rx.compartments)==1:
                print("product",rx.id,rx.reaction,FBA[rx.id])
                ATP=ATP+abs(FBA[rx.id])
    if FBA[rx.id]<0:
        for reac in rx.reactants:
            if reac.id == 'atp_c' and len(rx.compartments)==1:
                print("reactants",rx.id,rx.reaction,FBA[rx.id])
                ATP=ATP+abs(FBA[rx.id])

                
                

for rx in model_MeanCancerBiopsy.metabolites.atp_g.reactions:
    if FBA[rx.id]>0:
        for prod in rx.products:
            if prod.id == 'atp_g' and len(rx.compartments)==1:
                print("product",rx.id,rx.reaction,FBA[rx.id])
                ATP=ATP+abs(FBA[rx.id])
    if FBA[rx.id]<0:
        for reac in rx.reactants:
            if reac.id == 'atp_g' and len(rx.compartments)==1:
                print("reactants",rx.id,rx.reaction,FBA[rx.id])
                ATP=ATP+abs(FBA[rx.id])
                
                

                
for rx in model_MeanCancerBiopsy.metabolites.atp_l.reactions:
    if FBA[rx.id]>0:
        for prod in rx.products:
            if prod.id == 'atp_l' and len(rx.compartments)==1:
                print("product",rx.id,rx.reaction,FBA[rx.id])
                ATP=ATP+abs(FBA[rx.id])

    if FBA[rx.id]<0:
        for reac in rx.reactants:
            if reac.id == 'atp_l' and len(rx.compartments)==1:
                print("reactants",rx.id,rx.reaction,FBA[rx.id])
                ATP=ATP+abs(FBA[rx.id])

                
                
                
for rx in model_MeanCancerBiopsy.metabolites.atp_m.reactions:
    if FBA[rx.id]>0:
        for prod in rx.products:
            if prod.id == 'atp_m' and len(rx.compartments)==1:
                print("product",rx.id,rx.reaction,FBA[rx.id])
                ATP=ATP+abs(FBA[rx.id])

    if FBA[rx.id]<0:
        for reac in rx.reactants:
            if reac.id == 'atp_m' and len(rx.compartments)==1:
                print("reactants",rx.id,rx.reaction,FBA[rx.id])
                ATP=ATP+abs(FBA[rx.id])

                
                
                
for rx in model_MeanCancerBiopsy.metabolites.atp_n.reactions:
    if FBA[rx.id]>0:
        for prod in rx.products:
            if prod.id == 'atp_n' and len(rx.compartments)==1:
                print("product",rx.id,rx.reaction,FBA[rx.id])
                ATP=ATP+abs(FBA[rx.id])

    if FBA[rx.id]<0:
        for reac in rx.reactants:
            if reac.id == 'atp_n' and len(rx.compartments)==1:
                print("reactants",rx.id,rx.reaction,FBA[rx.id])
                ATP=ATP+abs(FBA[rx.id])

                
                
                
for rx in model_MeanCancerBiopsy.metabolites.atp_r.reactions:
    if FBA[rx.id]>0:
        for prod in rx.products:
            if prod.id == 'atp_r' and len(rx.compartments)==1:
                print("product",rx.id,rx.reaction,FBA[rx.id])
                ATP=ATP+abs(FBA[rx.id])

    if FBA[rx.id]<0:
        for reac in rx.reactants:
            if reac.id == 'atp_r' and len(rx.compartments)==1:
                print("reactants",rx.id,rx.reaction,FBA[rx.id])
                ATP=ATP+abs(FBA[rx.id])

                
                

for rx in model_MeanCancerBiopsy.metabolites.atp_x.reactions:
    if FBA[rx.id]>0:
        for prod in rx.products:
            if prod.id == 'atp_x' and len(rx.compartments)==1:
                print("product",rx.id,rx.reaction,FBA[rx.id])
                ATP=ATP+abs(FBA[rx.id])

    if FBA[rx.id]<0:
        for reac in rx.reactants:
            if reac.id == 'atp_x' and len(rx.compartments)==1:
                print("reactants",rx.id,rx.reaction,FBA[rx.id])
                ATP=ATP+abs(FBA[rx.id])




print(ATP)
FBA=cobra.flux_analysis.pfba(model_MeanNormalBiopsy)

ATP=0

for rx in model_MeanNormalBiopsy.metabolites.atp_c.reactions:
    if FBA[rx.id]>0:
        for prod in rx.products:
            if prod.id == 'atp_c' and len(rx.compartments)==1:
                print("product",rx.id,rx.reaction,FBA[rx.id])
                ATP=ATP+abs(FBA[rx.id])
    if FBA[rx.id]<0:
        for reac in rx.reactants:
            if reac.id == 'atp_c' and len(rx.compartments)==1:
                print("reactants",rx.id,rx.reaction,FBA[rx.id])
                ATP=ATP+abs(FBA[rx.id])

                
                

for rx in model_MeanNormalBiopsy.metabolites.atp_g.reactions:
    if FBA[rx.id]>0:
        for prod in rx.products:
            if prod.id == 'atp_g' and len(rx.compartments)==1:
                print("product",rx.id,rx.reaction,FBA[rx.id])
                ATP=ATP+abs(FBA[rx.id])
    if FBA[rx.id]<0:
        for reac in rx.reactants:
            if reac.id == 'atp_g' and len(rx.compartments)==1:
                print("reactants",rx.id,rx.reaction,FBA[rx.id])
                ATP=ATP+abs(FBA[rx.id])
                
for rx in model_MeanNormalBiopsy.metabolites.atp_l.reactions:
    if FBA[rx.id]>0:
        for prod in rx.products:
            if prod.id == 'atp_l' and len(rx.compartments)==1:
                print("product",rx.id,rx.reaction,FBA[rx.id])
                ATP=ATP+abs(FBA[rx.id])

    if FBA[rx.id]<0:
        for reac in rx.reactants:
            if reac.id == 'atp_l' and len(rx.compartments)==1:
                print("reactants",rx.id,rx.reaction,FBA[rx.id])
                ATP=ATP+abs(FBA[rx.id])
                
for rx in model_MeanNormalBiopsy.metabolites.atp_m.reactions:
    if FBA[rx.id]>0:
        for prod in rx.products:
            if prod.id == 'atp_m' and len(rx.compartments)==1:
                print("product",rx.id,rx.reaction,FBA[rx.id])
                ATP=ATP+abs(FBA[rx.id])

    if FBA[rx.id]<0:
        for reac in rx.reactants:
            if reac.id == 'atp_m' and len(rx.compartments)==1:
                print("reactants",rx.id,rx.reaction,FBA[rx.id])
                ATP=ATP+abs(FBA[rx.id])                
                
for rx in model_MeanNormalBiopsy.metabolites.atp_n.reactions:
    if FBA[rx.id]>0:
        for prod in rx.products:
            if prod.id == 'atp_n' and len(rx.compartments)==1:
                print("product",rx.id,rx.reaction,FBA[rx.id])
                ATP=ATP+abs(FBA[rx.id])

    if FBA[rx.id]<0:
        for reac in rx.reactants:
            if reac.id == 'atp_n' and len(rx.compartments)==1:
                print("reactants",rx.id,rx.reaction,FBA[rx.id])
                ATP=ATP+abs(FBA[rx.id])              
                
                
for rx in model_MeanNormalBiopsy.metabolites.atp_r.reactions:
    if FBA[rx.id]>0:
        for prod in rx.products:
            if prod.id == 'atp_r' and len(rx.compartments)==1:
                print("product",rx.id,rx.reaction,FBA[rx.id])
                ATP=ATP+abs(FBA[rx.id])

    if FBA[rx.id]<0:
        for reac in rx.reactants:
            if reac.id == 'atp_r' and len(rx.compartments)==1:
                print("reactants",rx.id,rx.reaction,FBA[rx.id])
                ATP=ATP+abs(FBA[rx.id])                

for rx in model_MeanNormalBiopsy.metabolites.atp_x.reactions:
    if FBA[rx.id]>0:
        for prod in rx.products:
            if prod.id == 'atp_x' and len(rx.compartments)==1:
                print("product",rx.id,rx.reaction,FBA[rx.id])
                ATP=ATP+abs(FBA[rx.id])

    if FBA[rx.id]<0:
        for reac in rx.reactants:
            if reac.id == 'atp_x' and len(rx.compartments)==1:
                print("reactants",rx.id,rx.reaction,FBA[rx.id])
                ATP=ATP+abs(FBA[rx.id])


print(ATP)


```python
pFBA1=model_MeanCancerBiopsy.optimize()
pFBA1.fluxes.to_csv("FBA_meanCancerBiopsyBiomass.csv")
```


```python
pFBA2=model_MeanNormalBiopsy.optimize()
pFBA2.fluxes.to_csv("FBA_meanNormalBiopsyBiomass.csv")
```


```python
fluxes=pd.concat([pFBA1.fluxes, pFBA2.fluxes], axis=1)
fluxes.columns=["Cancer", "Normal"]
fluxes.head()
fluxes.to_csv("FBA_meanBiopsyBiomass.csv")
```


```python
import cobra.test
from cobra.flux_analysis import production_envelope

prod_env = production_envelope(
model_MeanNormalBiopsy, ["biomass_reaction"], objective="biomass_reaction", c_source='EX_glc_LPAREN_e_RPAREN_' )
```


```python
%matplotlib inline
prod_env.plot(kind='line', x='EX_glc_LPAREN_e_RPAREN_', y='flux');
```


```python

```

## Production_envelope


```python
import cobra.test
from cobra.flux_analysis import production_envelope

prod_env = production_envelope(model_MeanCancerBiopsy, ["HEX1", "EX_o2_LPAREN_e_RPAREN_"])
```


```python
FBA=model_MeanCancerBiopsy.optimize()
for rx in model_MeanCancerBiopsy.metabolites.glc_D_c.reactions:
    print(rx)
    print(FBA[rx.id])
    print(rx.name)
    print(rx.reaction)
```


```python
model_MeanCancerBiopsy.reactions.GLCter.gene_reaction_rule

```


```python
%matplotlib inline 
import cobra.test
from cobra.flux_analysis import production_envelope

prod_env = production_envelope(
model_MeanCancerBiopsy, ["HEX1", "EX_o2_LPAREN_e_RPAREN_"], objective="biomass_reaction") #, carbon_sources="EX_glc__D_e"
prod_env.head()

prod_env[prod_env["flux"]>0].plot(kind='line', x='HEX1', y='flux');
```


```python
prod_env2 = prod_env[prod_env["flux"] > 0]
prod_env2

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
%matplotlib inline

prod_env = production_envelope(
    model_MeanCancerBiopsy, ["EX_co2_LPAREN_e_RPAREN_","EX_glc_LPAREN_e_RPAREN_","flux",])
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

## Flux Summary dataframe


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
rxs=model_MaxCancerBiopsy.metabolites.get_by_id("pchol_hs_c").reactions

for rx in rxs:
    metbols=model_MaxCancerBiopsy.reactions.get_by_id(rx.id).metabolites
    for met in metbols:
        if(met.formula_weight>100):
            print(met.id, met.formula_weight)
```

# Sampling penalty Factor


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

res1=pool.apply_async(model_create, [Recon2, conf_MeanCancerBiopsy, 5, metas ,100] )
res2=pool.apply_async(model_create, [Recon2, conf_MeanCancerBiopsy, 5, metas ,500] )
res3=pool.apply_async(model_create, [Recon2, conf_MeanCancerBiopsy, 5, metas ,1000] )
res4=pool.apply_async(model_create, [Recon2, conf_MeanCancerBiopsy, 5, metas ,2500] )
res5=pool.apply_async(model_create, [Recon2, conf_MeanCancerBiopsy, 5, metas ,5000] )

res6=pool.apply_async(model_create, [Recon2, conf_MeanNormalBiopsy, 5, metas ,100] )
res7=pool.apply_async(model_create, [Recon2, conf_MeanNormalBiopsy, 5, metas ,500] )
res8=pool.apply_async(model_create, [Recon2, conf_MeanNormalBiopsy, 5, metas ,1000] )
res9=pool.apply_async(model_create, [Recon2, conf_MeanNormalBiopsy, 5, metas ,2500] )
res0=pool.apply_async(model_create, [Recon2, conf_MeanNormalBiopsy, 5, metas ,5000] )

opt_MeanCancerBiopsy1=res1.get() 
opt_MeanCancerBiopsy2=res2.get() 
opt_MeanCancerBiopsy3=res3.get() 
opt_MeanCancerBiopsy4=res4.get() 
opt_MeanCancerBiopsy5=res5.get() 


opt_MeanNormalBiopsy1=res6.get() 
opt_MeanNormalBiopsy2=res7.get() 
opt_MeanNormalBiopsy3=res8.get() 
opt_MeanNormalBiopsy4=res9.get() 
opt_MeanNormalBiopsy5=res0.get() 
```

    Process ForkPoolWorker-1:
    KeyboardInterrupt
    Process ForkPoolWorker-2:
    Traceback (most recent call last):
      File "/opt/conda/lib/python3.5/multiprocessing/process.py", line 249, in _bootstrap
        self.run()
      File "/opt/conda/lib/python3.5/multiprocessing/process.py", line 93, in run
        self._target(*self._args, **self._kwargs)
      File "/opt/conda/lib/python3.5/multiprocessing/pool.py", line 108, in worker
        task = get()
      File "/opt/conda/lib/python3.5/multiprocessing/queues.py", line 342, in get
        with self._rlock:
      File "/opt/conda/lib/python3.5/multiprocessing/synchronize.py", line 96, in __enter__
        return self._semlock.__enter__()
    KeyboardInterrupt
    Traceback (most recent call last):
      File "/opt/conda/lib/python3.5/multiprocessing/process.py", line 249, in _bootstrap
        self.run()
      File "/opt/conda/lib/python3.5/multiprocessing/process.py", line 93, in run
        self._target(*self._args, **self._kwargs)
      File "/opt/conda/lib/python3.5/multiprocessing/pool.py", line 108, in worker
        task = get()
      File "/opt/conda/lib/python3.5/multiprocessing/queues.py", line 342, in get
        with self._rlock:
      File "/opt/conda/lib/python3.5/multiprocessing/synchronize.py", line 96, in __enter__
        return self._semlock.__enter__()
    KeyboardInterrupt
    Process ForkPoolWorker-9:
    Process ForkPoolWorker-8:
    Process ForkPoolWorker-7:
    Process ForkPoolWorker-21:
    Process ForkPoolWorker-31:
    Process ForkPoolWorker-13:
    Process ForkPoolWorker-15:
    Process ForkPoolWorker-23:
    Process ForkPoolWorker-17:
    Process ForkPoolWorker-24:
    Process ForkPoolWorker-18:
    Process ForkPoolWorker-30:
    Process ForkPoolWorker-22:
    Process ForkPoolWorker-19:
    Process ForkPoolWorker-16:
    Process ForkPoolWorker-11:
    Process ForkPoolWorker-20:
    Process ForkPoolWorker-14:
    Process ForkPoolWorker-12:
    Process ForkPoolWorker-10:
    Process ForkPoolWorker-4:
    Process ForkPoolWorker-27:
    Process ForkPoolWorker-29:
    Process ForkPoolWorker-25:
    Process ForkPoolWorker-26:
    Process ForkPoolWorker-28:
    Process ForkPoolWorker-32:
    Process ForkPoolWorker-6:
    Process ForkPoolWorker-3:
    Traceback (most recent call last):
    Traceback (most recent call last):
      File "/opt/conda/lib/python3.5/multiprocessing/process.py", line 249, in _bootstrap
        self.run()
      File "/opt/conda/lib/python3.5/multiprocessing/process.py", line 249, in _bootstrap
        self.run()
      File "/opt/conda/lib/python3.5/multiprocessing/process.py", line 93, in run
        self._target(*self._args, **self._kwargs)
      File "/opt/conda/lib/python3.5/multiprocessing/process.py", line 93, in run
        self._target(*self._args, **self._kwargs)
    Traceback (most recent call last):
      File "/opt/conda/lib/python3.5/multiprocessing/pool.py", line 108, in worker
        task = get()
    Traceback (most recent call last):
      File "/opt/conda/lib/python3.5/multiprocessing/pool.py", line 108, in worker
        task = get()
      File "/opt/conda/lib/python3.5/multiprocessing/queues.py", line 342, in get
        with self._rlock:
      File "/opt/conda/lib/python3.5/multiprocessing/process.py", line 249, in _bootstrap
        self.run()
      File "/opt/conda/lib/python3.5/multiprocessing/process.py", line 249, in _bootstrap
        self.run()
    Traceback (most recent call last):
      File "/opt/conda/lib/python3.5/multiprocessing/queues.py", line 342, in get
        with self._rlock:
    Traceback (most recent call last):
      File "/opt/conda/lib/python3.5/multiprocessing/synchronize.py", line 96, in __enter__
        return self._semlock.__enter__()
      File "/opt/conda/lib/python3.5/multiprocessing/process.py", line 93, in run
        self._target(*self._args, **self._kwargs)
      File "/opt/conda/lib/python3.5/multiprocessing/process.py", line 93, in run
        self._target(*self._args, **self._kwargs)
      File "/opt/conda/lib/python3.5/multiprocessing/synchronize.py", line 96, in __enter__
        return self._semlock.__enter__()
      File "/opt/conda/lib/python3.5/multiprocessing/process.py", line 249, in _bootstrap
        self.run()
      File "/opt/conda/lib/python3.5/multiprocessing/pool.py", line 108, in worker
        task = get()
      File "/opt/conda/lib/python3.5/multiprocessing/process.py", line 249, in _bootstrap
        self.run()
    KeyboardInterrupt
      File "/opt/conda/lib/python3.5/multiprocessing/queues.py", line 342, in get
        with self._rlock:
      File "/opt/conda/lib/python3.5/multiprocessing/pool.py", line 108, in worker
        task = get()
    KeyboardInterrupt
      File "/opt/conda/lib/python3.5/multiprocessing/process.py", line 93, in run
        self._target(*self._args, **self._kwargs)
    Traceback (most recent call last):
      File "/opt/conda/lib/python3.5/multiprocessing/process.py", line 93, in run
        self._target(*self._args, **self._kwargs)
      File "/opt/conda/lib/python3.5/multiprocessing/synchronize.py", line 96, in __enter__
        return self._semlock.__enter__()
      File "/opt/conda/lib/python3.5/multiprocessing/queues.py", line 342, in get
        with self._rlock:
      File "/opt/conda/lib/python3.5/multiprocessing/pool.py", line 108, in worker
        task = get()
    KeyboardInterrupt
      File "/opt/conda/lib/python3.5/multiprocessing/synchronize.py", line 96, in __enter__
        return self._semlock.__enter__()
      File "/opt/conda/lib/python3.5/multiprocessing/process.py", line 249, in _bootstrap
        self.run()
      File "/opt/conda/lib/python3.5/multiprocessing/pool.py", line 108, in worker
        task = get()
      File "/opt/conda/lib/python3.5/multiprocessing/queues.py", line 342, in get
        with self._rlock:
    KeyboardInterrupt
      File "/opt/conda/lib/python3.5/multiprocessing/process.py", line 93, in run
        self._target(*self._args, **self._kwargs)
      File "/opt/conda/lib/python3.5/multiprocessing/queues.py", line 342, in get
        with self._rlock:
      File "/opt/conda/lib/python3.5/multiprocessing/synchronize.py", line 96, in __enter__
        return self._semlock.__enter__()
      File "/opt/conda/lib/python3.5/multiprocessing/pool.py", line 108, in worker
        task = get()
    KeyboardInterrupt
      File "/opt/conda/lib/python3.5/multiprocessing/synchronize.py", line 96, in __enter__
        return self._semlock.__enter__()
    KeyboardInterrupt
      File "/opt/conda/lib/python3.5/multiprocessing/queues.py", line 342, in get
        with self._rlock:
      File "/opt/conda/lib/python3.5/multiprocessing/synchronize.py", line 96, in __enter__
        return self._semlock.__enter__()
    KeyboardInterrupt
    Traceback (most recent call last):
    Traceback (most recent call last):
    Traceback (most recent call last):
    Traceback (most recent call last):
    Traceback (most recent call last):
    Traceback (most recent call last):
    Traceback (most recent call last):
    Traceback (most recent call last):
    Traceback (most recent call last):
    Traceback (most recent call last):
    Traceback (most recent call last):
    Traceback (most recent call last):
    Traceback (most recent call last):
    Traceback (most recent call last):
    Traceback (most recent call last):
    Traceback (most recent call last):
    Traceback (most recent call last):
    Traceback (most recent call last):
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
      File "/opt/conda/lib/python3.5/multiprocessing/process.py", line 249, in _bootstrap
        self.run()
      File "/opt/conda/lib/python3.5/multiprocessing/process.py", line 249, in _bootstrap
        self.run()
      File "/opt/conda/lib/python3.5/multiprocessing/process.py", line 249, in _bootstrap
        self.run()
      File "/opt/conda/lib/python3.5/multiprocessing/process.py", line 249, in _bootstrap
        self.run()
      File "/opt/conda/lib/python3.5/multiprocessing/process.py", line 249, in _bootstrap
        self.run()
      File "/opt/conda/lib/python3.5/multiprocessing/process.py", line 249, in _bootstrap
        self.run()
      File "/opt/conda/lib/python3.5/multiprocessing/process.py", line 249, in _bootstrap
        self.run()
      File "/opt/conda/lib/python3.5/multiprocessing/process.py", line 249, in _bootstrap
        self.run()
      File "/opt/conda/lib/python3.5/multiprocessing/process.py", line 249, in _bootstrap
        self.run()
      File "/opt/conda/lib/python3.5/multiprocessing/process.py", line 249, in _bootstrap
        self.run()
      File "/opt/conda/lib/python3.5/multiprocessing/process.py", line 249, in _bootstrap
        self.run()
      File "/opt/conda/lib/python3.5/multiprocessing/process.py", line 249, in _bootstrap
        self.run()
      File "/opt/conda/lib/python3.5/multiprocessing/process.py", line 249, in _bootstrap
        self.run()
      File "/opt/conda/lib/python3.5/multiprocessing/process.py", line 249, in _bootstrap
        self.run()
      File "/opt/conda/lib/python3.5/multiprocessing/process.py", line 249, in _bootstrap
        self.run()
      File "/opt/conda/lib/python3.5/multiprocessing/process.py", line 249, in _bootstrap
        self.run()
      File "/opt/conda/lib/python3.5/multiprocessing/process.py", line 249, in _bootstrap
        self.run()



    ---------------------------------------------------------------------------

    KeyboardInterrupt                         Traceback (most recent call last)

    <timed exec> in <module>()


    /opt/conda/lib/python3.5/multiprocessing/pool.py in get(self, timeout)
        600 
        601     def get(self, timeout=None):
    --> 602         self.wait(timeout)
        603         if not self.ready():
        604             raise TimeoutError


    /opt/conda/lib/python3.5/multiprocessing/pool.py in wait(self, timeout)
        597 
        598     def wait(self, timeout=None):
    --> 599         self._event.wait(timeout)
        600 
        601     def get(self, timeout=None):


    /opt/conda/lib/python3.5/threading.py in wait(self, timeout)
        547             signaled = self._flag
        548             if not signaled:
    --> 549                 signaled = self._cond.wait(timeout)
        550             return signaled
        551 


    /opt/conda/lib/python3.5/threading.py in wait(self, timeout)
        291         try:    # restore state no matter what (e.g., KeyboardInterrupt)
        292             if timeout is None:
    --> 293                 waiter.acquire()
        294                 gotit = True
        295             else:


    KeyboardInterrupt: 


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
      File "/opt/conda/lib/python3.5/multiprocessing/process.py", line 93, in run
        self._target(*self._args, **self._kwargs)
      File "/opt/conda/lib/python3.5/multiprocessing/process.py", line 93, in run
        self._target(*self._args, **self._kwargs)
      File "/opt/conda/lib/python3.5/multiprocessing/process.py", line 93, in run
        self._target(*self._args, **self._kwargs)
      File "/opt/conda/lib/python3.5/multiprocessing/process.py", line 93, in run
        self._target(*self._args, **self._kwargs)
      File "/opt/conda/lib/python3.5/multiprocessing/process.py", line 93, in run
        self._target(*self._args, **self._kwargs)
      File "/opt/conda/lib/python3.5/multiprocessing/process.py", line 93, in run
        self._target(*self._args, **self._kwargs)
      File "/opt/conda/lib/python3.5/multiprocessing/process.py", line 93, in run
        self._target(*self._args, **self._kwargs)
      File "/opt/conda/lib/python3.5/multiprocessing/process.py", line 93, in run
        self._target(*self._args, **self._kwargs)
      File "/opt/conda/lib/python3.5/multiprocessing/process.py", line 93, in run
        self._target(*self._args, **self._kwargs)
      File "/opt/conda/lib/python3.5/multiprocessing/process.py", line 93, in run
        self._target(*self._args, **self._kwargs)
      File "/opt/conda/lib/python3.5/multiprocessing/process.py", line 93, in run
        self._target(*self._args, **self._kwargs)
      File "/opt/conda/lib/python3.5/multiprocessing/process.py", line 93, in run
        self._target(*self._args, **self._kwargs)
      File "/opt/conda/lib/python3.5/multiprocessing/process.py", line 93, in run
        self._target(*self._args, **self._kwargs)
      File "/opt/conda/lib/python3.5/multiprocessing/process.py", line 93, in run
        self._target(*self._args, **self._kwargs)
      File "/opt/conda/lib/python3.5/multiprocessing/process.py", line 93, in run
        self._target(*self._args, **self._kwargs)
      File "/opt/conda/lib/python3.5/multiprocessing/process.py", line 93, in run
        self._target(*self._args, **self._kwargs)
      File "/opt/conda/lib/python3.5/multiprocessing/process.py", line 93, in run
        self._target(*self._args, **self._kwargs)
      File "/opt/conda/lib/python3.5/multiprocessing/process.py", line 93, in run
        self._target(*self._args, **self._kwargs)
      File "/opt/conda/lib/python3.5/multiprocessing/pool.py", line 108, in worker
        task = get()
      File "/opt/conda/lib/python3.5/multiprocessing/pool.py", line 108, in worker
        task = get()
      File "/opt/conda/lib/python3.5/multiprocessing/pool.py", line 108, in worker
        task = get()
      File "/opt/conda/lib/python3.5/multiprocessing/pool.py", line 108, in worker
        task = get()
      File "/opt/conda/lib/python3.5/multiprocessing/pool.py", line 108, in worker
        task = get()
      File "/opt/conda/lib/python3.5/multiprocessing/pool.py", line 108, in worker
        task = get()
      File "/opt/conda/lib/python3.5/multiprocessing/pool.py", line 108, in worker
        task = get()
      File "/opt/conda/lib/python3.5/multiprocessing/pool.py", line 108, in worker
        task = get()
      File "/opt/conda/lib/python3.5/multiprocessing/pool.py", line 108, in worker
        task = get()
      File "/opt/conda/lib/python3.5/multiprocessing/pool.py", line 108, in worker
        task = get()
      File "/opt/conda/lib/python3.5/multiprocessing/pool.py", line 108, in worker
        task = get()
      File "/opt/conda/lib/python3.5/multiprocessing/pool.py", line 108, in worker
        task = get()
      File "/opt/conda/lib/python3.5/multiprocessing/pool.py", line 108, in worker
        task = get()
      File "/opt/conda/lib/python3.5/multiprocessing/pool.py", line 108, in worker
        task = get()
      File "/opt/conda/lib/python3.5/multiprocessing/pool.py", line 108, in worker
        task = get()
      File "/opt/conda/lib/python3.5/multiprocessing/pool.py", line 108, in worker
        task = get()
      File "/opt/conda/lib/python3.5/multiprocessing/pool.py", line 108, in worker
        task = get()
      File "/opt/conda/lib/python3.5/multiprocessing/pool.py", line 108, in worker
        task = get()
      File "/opt/conda/lib/python3.5/multiprocessing/pool.py", line 108, in worker
        task = get()
      File "/opt/conda/lib/python3.5/multiprocessing/pool.py", line 108, in worker
        task = get()
      File "/opt/conda/lib/python3.5/multiprocessing/pool.py", line 108, in worker
        task = get()
      File "/opt/conda/lib/python3.5/multiprocessing/pool.py", line 108, in worker
        task = get()
      File "/opt/conda/lib/python3.5/multiprocessing/queues.py", line 342, in get
        with self._rlock:
      File "/opt/conda/lib/python3.5/multiprocessing/queues.py", line 342, in get
        with self._rlock:
      File "/opt/conda/lib/python3.5/multiprocessing/queues.py", line 342, in get
        with self._rlock:
      File "/opt/conda/lib/python3.5/multiprocessing/queues.py", line 342, in get
        with self._rlock:
      File "/opt/conda/lib/python3.5/multiprocessing/queues.py", line 342, in get
        with self._rlock:
      File "/opt/conda/lib/python3.5/multiprocessing/queues.py", line 342, in get
        with self._rlock:
      File "/opt/conda/lib/python3.5/multiprocessing/queues.py", line 342, in get
        with self._rlock:
      File "/opt/conda/lib/python3.5/multiprocessing/queues.py", line 342, in get
        with self._rlock:
      File "/opt/conda/lib/python3.5/multiprocessing/queues.py", line 342, in get
        with self._rlock:
      File "/opt/conda/lib/python3.5/multiprocessing/queues.py", line 342, in get
        with self._rlock:
      File "/opt/conda/lib/python3.5/multiprocessing/queues.py", line 342, in get
        with self._rlock:
      File "/opt/conda/lib/python3.5/multiprocessing/queues.py", line 342, in get
        with self._rlock:
      File "/opt/conda/lib/python3.5/multiprocessing/queues.py", line 342, in get
        with self._rlock:
      File "/opt/conda/lib/python3.5/multiprocessing/queues.py", line 343, in get
        res = self._reader.recv_bytes()
      File "/opt/conda/lib/python3.5/multiprocessing/queues.py", line 342, in get
        with self._rlock:
      File "/opt/conda/lib/python3.5/multiprocessing/queues.py", line 342, in get
        with self._rlock:
      File "/opt/conda/lib/python3.5/multiprocessing/queues.py", line 342, in get
        with self._rlock:
      File "/opt/conda/lib/python3.5/multiprocessing/queues.py", line 342, in get
        with self._rlock:
      File "/opt/conda/lib/python3.5/multiprocessing/queues.py", line 342, in get
        with self._rlock:
      File "/opt/conda/lib/python3.5/multiprocessing/queues.py", line 342, in get
        with self._rlock:
      File "/opt/conda/lib/python3.5/multiprocessing/queues.py", line 342, in get
        with self._rlock:
      File "/opt/conda/lib/python3.5/multiprocessing/queues.py", line 342, in get
        with self._rlock:
      File "/opt/conda/lib/python3.5/multiprocessing/synchronize.py", line 96, in __enter__
        return self._semlock.__enter__()
      File "/opt/conda/lib/python3.5/multiprocessing/synchronize.py", line 96, in __enter__
        return self._semlock.__enter__()
      File "/opt/conda/lib/python3.5/multiprocessing/connection.py", line 216, in recv_bytes
        buf = self._recv_bytes(maxlength)
      File "/opt/conda/lib/python3.5/multiprocessing/synchronize.py", line 96, in __enter__
        return self._semlock.__enter__()
      File "/opt/conda/lib/python3.5/multiprocessing/synchronize.py", line 96, in __enter__
        return self._semlock.__enter__()
      File "/opt/conda/lib/python3.5/multiprocessing/synchronize.py", line 96, in __enter__
        return self._semlock.__enter__()
      File "/opt/conda/lib/python3.5/multiprocessing/synchronize.py", line 96, in __enter__
        return self._semlock.__enter__()
      File "/opt/conda/lib/python3.5/multiprocessing/synchronize.py", line 96, in __enter__
        return self._semlock.__enter__()
      File "/opt/conda/lib/python3.5/multiprocessing/synchronize.py", line 96, in __enter__
        return self._semlock.__enter__()
      File "/opt/conda/lib/python3.5/multiprocessing/synchronize.py", line 96, in __enter__
        return self._semlock.__enter__()
      File "/opt/conda/lib/python3.5/multiprocessing/synchronize.py", line 96, in __enter__
        return self._semlock.__enter__()
      File "/opt/conda/lib/python3.5/multiprocessing/synchronize.py", line 96, in __enter__
        return self._semlock.__enter__()
      File "/opt/conda/lib/python3.5/multiprocessing/synchronize.py", line 96, in __enter__
        return self._semlock.__enter__()
      File "/opt/conda/lib/python3.5/multiprocessing/synchronize.py", line 96, in __enter__
        return self._semlock.__enter__()
      File "/opt/conda/lib/python3.5/multiprocessing/synchronize.py", line 96, in __enter__
        return self._semlock.__enter__()
      File "/opt/conda/lib/python3.5/multiprocessing/synchronize.py", line 96, in __enter__
        return self._semlock.__enter__()
      File "/opt/conda/lib/python3.5/multiprocessing/synchronize.py", line 96, in __enter__
        return self._semlock.__enter__()
      File "/opt/conda/lib/python3.5/multiprocessing/synchronize.py", line 96, in __enter__
        return self._semlock.__enter__()
      File "/opt/conda/lib/python3.5/multiprocessing/synchronize.py", line 96, in __enter__
        return self._semlock.__enter__()
      File "/opt/conda/lib/python3.5/multiprocessing/synchronize.py", line 96, in __enter__
        return self._semlock.__enter__()
      File "/opt/conda/lib/python3.5/multiprocessing/synchronize.py", line 96, in __enter__
        return self._semlock.__enter__()
      File "/opt/conda/lib/python3.5/multiprocessing/synchronize.py", line 96, in __enter__
        return self._semlock.__enter__()
    KeyboardInterrupt
    KeyboardInterrupt
    KeyboardInterrupt
    KeyboardInterrupt
    KeyboardInterrupt
    KeyboardInterrupt
    KeyboardInterrupt
    KeyboardInterrupt
      File "/opt/conda/lib/python3.5/multiprocessing/connection.py", line 407, in _recv_bytes
        buf = self._recv(4)
    KeyboardInterrupt
    KeyboardInterrupt
    KeyboardInterrupt
    KeyboardInterrupt
    KeyboardInterrupt
    KeyboardInterrupt
    KeyboardInterrupt
    KeyboardInterrupt
    KeyboardInterrupt
    KeyboardInterrupt
    KeyboardInterrupt
    KeyboardInterrupt
    KeyboardInterrupt
      File "/opt/conda/lib/python3.5/multiprocessing/connection.py", line 379, in _recv
        chunk = read(handle, remaining)



```python
model_MeanCancerBiopsy1=opt_MeanCancerBiopsy1.cobra_model(name="Cancer1")
print(model_MeanCancerBiopsy1.optimize())
model_MeanCancerBiopsy2=opt_MeanCancerBiopsy2.cobra_model(name="Cancer2")
print(model_MeanCancerBiopsy2.optimize())
model_MeanCancerBiopsy3=opt_MeanCancerBiopsy3.cobra_model(name="Cancer3")
print(model_MeanCancerBiopsy3.optimize())
model_MeanCancerBiopsy4=opt_MeanCancerBiopsy4.cobra_model(name="Cancer4")
print(model_MeanCancerBiopsy4.optimize())
model_MeanCancerBiopsy5=opt_MeanCancerBiopsy5.cobra_model(name="Cancer5")
print(model_MeanCancerBiopsy5.optimize())
```

    <Solution 0.029 at 0x7f2a1e7be748>
    <Solution 0.033 at 0x7f2a1f1ed978>
    <Solution 0.033 at 0x7f2a34213fd0>
    <Solution 0.033 at 0x7f2a25faef60>



    ---------------------------------------------------------------------------

    NameError                                 Traceback (most recent call last)

    <ipython-input-63-b97041049df6> in <module>()
          7 model_MeanCancerBiopsy4=opt_MeanCancerBiopsy4.cobra_model(name="Cancer4")
          8 print(model_MeanCancerBiopsy4.optimize())
    ----> 9 model_MeanCancerBiopsy5=opt_MeanCancerBiopsy5.cobra_model(name="Cancer5")
         10 print(model_MeanCancerBiopsy5.optimize())


    NameError: name 'opt_MeanCancerBiopsy5' is not defined



```python
model_MeanNormalBiopsy1=opt_MeanNormalBiopsy1.cobra_model(name="Norm1")
print(model_MeanNormalBiopsy1.optimize())
model_MeanNormalBiopsy2=opt_MeanNormalBiopsy2.cobra_model(name="Norm2")
print(model_MeanNormalBiopsy2.optimize())
model_MeanNormalBiopsy3=opt_MeanNormalBiopsy3.cobra_model(name="Norm3")
print(model_MeanNormalBiopsy3.optimize())
model_MeanNormalBiopsy4=opt_MeanNormalBiopsy4.cobra_model(name="Norm4")
print(model_MeanNormalBiopsy4.optimize())
model_MeanNormalBiopsy5=opt_MeanNormalBiopsy5.cobra_model(name="Norm5")
print(model_MeanNormalBiopsy5.optimize())

```


```python
model_MeanCancerBiopsy1.summary(fva=True)
```


```python
model_MeanCancerBiopsy2.summary(fva=True)
```


```python
model_MeanNormalBiopsy2.summary(fva=True)
```


```python
model_MeanCancerBiopsy3.summary(fva=True)
```


```python
model_MeanNormalBiopsy3.summary(fva=True)
```


```python
model_MeanCancerBiopsy4.summary(fva=True)
```


```python
model_MeanNormalBiopsy4.summary(fva=True)
```


```python
model_MeanCancerBiopsy5.summary(fva=True)
```


```python
model_MeanNormalBiopsy5.summary(fva=True)
```


```python
R_Biopsies=get_drugable_targets(
    normal_Model=model_MeanNormalBiopsy1, 
    disease_Model=model_MeanCancerBiopsy1, model_name="Biopsies")
R_Biopsies["Dise_ratio"]=R_Biopsies["del_dise_flux"]/R_Biopsies["dise_flux"]
R_Biopsies["Norm_ratio"]=R_Biopsies["del_norm_flux"]/R_Biopsies["norm_flux"]
R_Biopsies[R_Biopsies["Dise_ratio"]< 0.999].sort_values(by="Dise_ratio", ascending=True)
```


```python
R_Biopsies=get_drugable_targets(
    normal_Model=model_MeanNormalBiopsy2, 
    disease_Model=model_MeanCancerBiopsy2, model_name="Biopsies")
R_Biopsies["Dise_ratio"]=R_Biopsies["del_dise_flux"]/R_Biopsies["dise_flux"]
R_Biopsies["Norm_ratio"]=R_Biopsies["del_norm_flux"]/R_Biopsies["norm_flux"]
R_Biopsies[R_Biopsies["Dise_ratio"]< 0.999].sort_values(by="Dise_ratio", ascending=True)
```


```python
R_Biopsies=get_drugable_targets(
    normal_Model=model_MeanNormalBiopsy3, 
    disease_Model=model_MeanCancerBiopsy3, model_name="Biopsies")
R_Biopsies["Dise_ratio"]=R_Biopsies["del_dise_flux"]/R_Biopsies["dise_flux"]
R_Biopsies["Norm_ratio"]=R_Biopsies["del_norm_flux"]/R_Biopsies["norm_flux"]
R_Biopsies[R_Biopsies["Dise_ratio"]< 0.999].sort_values(by="Dise_ratio", ascending=True)
```


```python
R_Biopsies=get_drugable_targets(
    normal_Model=model_MeanNormalBiopsy4, 
    disease_Model=model_MeanCancerBiopsy4, model_name="Biopsies")
R_Biopsies["Dise_ratio"]=R_Biopsies["del_dise_flux"]/R_Biopsies["dise_flux"]
R_Biopsies["Norm_ratio"]=R_Biopsies["del_norm_flux"]/R_Biopsies["norm_flux"]
R_Biopsies[R_Biopsies["Dise_ratio"]< 0.999].sort_values(by="Dise_ratio", ascending=True)
```


```python
R_Biopsies=get_drugable_targets(
    normal_Model=model_MeanNormalBiopsy5, 
    disease_Model=model_MeanCancerBiopsy5, model_name="Biopsies")
R_Biopsies["Dise_ratio"]=R_Biopsies["del_dise_flux"]/R_Biopsies["dise_flux"]
R_Biopsies["Norm_ratio"]=R_Biopsies["del_norm_flux"]/R_Biopsies["norm_flux"]
R_Biopsies[R_Biopsies["Dise_ratio"]< 0.999].sort_values(by="Dise_ratio", ascending=True)
```


```python
cobra.io.write_sbml_model(model_MeanCancerBiopsy1, "Full_SCC_model1_0519_n5pf100.sbml")
cobra.io.write_sbml_model(model_MeanNormalBiopsy1, "Full_Normal_model1_0519_n5pf100.sbml")

cobra.io.write_sbml_model(model_MeanCancerBiopsy2, "Full_SCC_model2_0519_n5pf500.sbml")
cobra.io.write_sbml_model(model_MeanNormalBiopsy2, "Full_Normal_model2_0519_n5pf500.sbml")

cobra.io.write_sbml_model(model_MeanCancerBiopsy3, "Full_SCC_model3_0519_n5pf1000.sbml")
cobra.io.write_sbml_model(model_MeanNormalBiopsy3, "Full_Normal_model3_0519_n5pf1000.sbml")

cobra.io.write_sbml_model(model_MeanCancerBiopsy4, "Full_SCC_model4_0519_n5pf2500.sbml")
cobra.io.write_sbml_model(model_MeanNormalBiopsy4, "Full_Normal_model4_0519_n5pf2500.sbml")


cobra.io.write_sbml_model(model_MeanCancerBiopsy5, "Full_SCC_model5_0519_n5pf5000.sbml")
cobra.io.write_sbml_model(model_MeanNormalBiopsy5, "Full_Normal_model5_0519_n5pf5000.sbml")
```


```python
FBA1=model_MeanCancerBiopsy1.optimize()
FBA2=model_MeanNormalBiopsy1.optimize() 

result = pd.concat([FBA1.fluxes, FBA2.fluxes], axis=1)
result.columns = ["CancerBiopsy","NormalBiopsy"]
result.to_csv("BiopsiesFBA1.csv")

FBA1.fluxes.to_csv("CancerBiopsyFBA1.csv")
FBA2.fluxes.to_csv("NormalBiopsyFBA1.csv")
```


```python
FBA1=model_MeanCancerBiopsy2.optimize()
FBA2=model_MeanNormalBiopsy2.optimize() 

result = pd.concat([FBA1.fluxes, FBA2.fluxes], axis=1)
result.columns = ["CancerBiopsy","NormalBiopsy"]
result.to_csv("BiopsiesFBA2.csv")

FBA1.fluxes.to_csv("CancerBiopsyFBA2.csv")
FBA2.fluxes.to_csv("NormalBiopsyFBA2.csv")
```


```python
FBA1=model_MeanCancerBiopsy3.optimize()
FBA2=model_MeanNormalBiopsy3.optimize() 

result = pd.concat([FBA1.fluxes, FBA2.fluxes], axis=1)
result.columns = ["CancerBiopsy","NormalBiopsy"]
result.to_csv("BiopsiesFBA3.csv")

FBA1.fluxes.to_csv("CancerBiopsyFBA3.csv")
FBA2.fluxes.to_csv("NormalBiopsyFBA3.csv")
```


```python
FBA1=model_MeanCancerBiopsy4.optimize()
FBA2=model_MeanNormalBiopsy4.optimize() 

result = pd.concat([FBA1.fluxes, FBA2.fluxes], axis=1)
result.columns = ["CancerBiopsy","NormalBiopsy"]
result.to_csv("BiopsiesFBA4.csv")

FBA1.fluxes.to_csv("CancerBiopsyFBA4.csv")
FBA2.fluxes.to_csv("NormalBiopsyFBA4.csv")
```


```python
FBA1=model_MeanCancerBiopsy5.optimize()
FBA2=model_MeanNormalBiopsy5.optimize() 

result = pd.concat([FBA1.fluxes, FBA2.fluxes], axis=1)
result.columns = ["CancerBiopsy","NormalBiopsy"]
result.to_csv("BiopsiesFBA5.csv")

FBA1.fluxes.to_csv("CancerBiopsyFBA5.csv")
FBA2.fluxes.to_csv("NormalBiopsyFBA5.csv")
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
def multiple_randomized_analysis(modelC, modelN, Recon, pos, metabolites):
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
    
    opt_n= CORDA(model=Recon, confidence=nConf_scores, n=10, met_prod=metabolites,  penalty_factor=1000 ) 
    opt_c.build()
    opt_n.build()
    model_c=opt_c.cobra_model(name="C")
    model_n=opt_n.cobra_model(name="N")
     
    targets=get_drugable_targets(model_n, model_c, "biopsys" )
    f="Targets_pf_1000_Byopsies_" + str(pos) + ".tsv"                                                                  
    targets.to_csv(f)

# Setup a list of processes that we want to run
pool = mp.Pool(processes=28)

for i in range(1,1001):
    pool.apply_async(multiple_randomized_analysis, args=(conf_MeanCancerBiopsy,conf_MeanNormalBiopsy,Recon2,i,metas))

pool.close()
pool.join()


```

closedModel=model_MeanCancerBiopsy.copy()
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
closedModel.reactions.EX_co2_LPAREN_e_RPAREN_.bounds=(1,1000)
closedModel.reactions.EX_glc_LPAREN_e_RPAREN_.bounds=(-1,-1)
closedModel.reactions.EX_gln_L_LPAREN_e_RPAREN_.bounds=(0,0)
closedModel.reactions.EX_hco3_LPAREN_e_RPAREN_.bounds=(0,0)

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

FBA.fluxes.to_csv("pFBA_Glucose_aerobic_cancerBiopsy.csv")

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

FBA.fluxes.to_csv("pFBA_Glucose_anaerobic_cancerBiopsy.csv")

######################################################################
## Glutamine aerobic
######################################################################
closedModel.reactions.EX_o2_LPAREN_e_RPAREN_.bounds=(-1000,-1)
closedModel.reactions.EX_h2o_LPAREN_e_RPAREN_.bounds=(-1000,1000)
closedModel.reactions.EX_co2_LPAREN_e_RPAREN_.bounds=(0,1000)
closedModel.reactions.EX_glc_LPAREN_e_RPAREN_.bounds=(0,0)
closedModel.reactions.EX_gln_L_LPAREN_e_RPAREN_.bounds=(-1,-1)


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

FBA.fluxes.to_csv("pFBA_Glutamine_aerobic_cancerBiopsy.csv")


######################################################################
## Glutamine anaerobic
######################################################################
closedModel.reactions.EX_o2_LPAREN_e_RPAREN_.bounds=(0,0)
closedModel.reactions.EX_h2o_LPAREN_e_RPAREN_.bounds=(-1000,1000)
closedModel.reactions.EX_co2_LPAREN_e_RPAREN_.bounds=(0,1000)
closedModel.reactions.EX_glc_LPAREN_e_RPAREN_.bounds=(0,0)
closedModel.reactions.EX_gln_L_LPAREN_e_RPAREN_.bounds=(-1,-1)


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

FBA.fluxes.to_csv("pFBA_Glutamine_anaerobic_cancerBiopsy.csv")



```python
for rx in model_MeanNormalBiopsy4.reactions:
    print(rx.id)
```


```python

```


```python

```
