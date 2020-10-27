

```python
import cobra
hgu133APlus2_hela=cobra.io.read_sbml_model("Models/hgu133APlus2_hela_model_0418_DEMEM6429_n5.sbml")
hgu133APlus2_hacat=cobra.io.read_sbml_model("Models/hgu133APlus2_hacat_model_0418_DEMEM6429_n5.sbml")

hgu133APlus2_scc=cobra.io.read_sbml_model("Models/hgu133APlus2_cancer_model_0418_westDiet_n5.sbml")
hgu133APlus2_normal=cobra.io.read_sbml_model("Models/hgu133APlus2_normal_model_0418_westDiet_n5.sbml")


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
        
        results[rx]["genes"]=disease_Model.reactions.get_by_id(rx).gene_name_reaction_rule

        dmodel.reactions.get_by_id(rx).bounds=dbounds


    return(pd.DataFrame(results).transpose())

```


```python
R_clines=get_drugable_targets(hgu133APlus2_hacat, hgu133APlus2_hela, "hgu133APlus2_clines")

```

    Common reactions size 2620
    Normal unique reactions size 442
    Disease unique reactions size 582



```python
R_biopsys=get_drugable_targets(hgu133APlus2_normal, hgu133APlus2_scc, "hgu133APlus2_biopsys" )

```

    Common reactions size 1963
    Normal unique reactions size 524
    Disease unique reactions size 614



```python
R_biopsys
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
      <th>10FTHF5GLUtl</th>
      <td>0.216258</td>
      <td>1</td>
      <td></td>
      <td>hgu133APlus2_biopsys</td>
      <td>1</td>
      <td>0.128673</td>
      <td>0.128673</td>
    </tr>
    <tr>
      <th>10FTHF6GLUtl</th>
      <td>0.216258</td>
      <td>1</td>
      <td></td>
      <td>hgu133APlus2_biopsys</td>
      <td>1</td>
      <td>0.128673</td>
      <td>0.128673</td>
    </tr>
    <tr>
      <th>10FTHF7GLUtl</th>
      <td>0.216258</td>
      <td>1</td>
      <td></td>
      <td>hgu133APlus2_biopsys</td>
      <td>1</td>
      <td>0.128673</td>
      <td>0.128673</td>
    </tr>
    <tr>
      <th>10FTHFtl</th>
      <td>0.216258</td>
      <td>1</td>
      <td></td>
      <td>hgu133APlus2_biopsys</td>
      <td>1</td>
      <td>0.128673</td>
      <td>0.128673</td>
    </tr>
    <tr>
      <th>10FTHFtm</th>
      <td>0.216258</td>
      <td>1</td>
      <td></td>
      <td>hgu133APlus2_biopsys</td>
      <td>1</td>
      <td>0.128673</td>
      <td>0.128673</td>
    </tr>
    <tr>
      <th>12HTACRhr</th>
      <td>0.216258</td>
      <td>0.216258</td>
      <td>HGNC:2637 or HGNC:2638</td>
      <td>hgu133APlus2_biopsys</td>
      <td>1</td>
      <td>0.128673</td>
      <td>1</td>
    </tr>
    <tr>
      <th>12HTACRitr</th>
      <td>0.216258</td>
      <td>0.216258</td>
      <td></td>
      <td>hgu133APlus2_biopsys</td>
      <td>1</td>
      <td>0.128673</td>
      <td>1</td>
    </tr>
    <tr>
      <th>12HTACRtep</th>
      <td>0.216258</td>
      <td>0.216258</td>
      <td></td>
      <td>hgu133APlus2_biopsys</td>
      <td>1</td>
      <td>0.128673</td>
      <td>1</td>
    </tr>
    <tr>
      <th>12HTACRtu</th>
      <td>0.216258</td>
      <td>0.216258</td>
      <td></td>
      <td>hgu133APlus2_biopsys</td>
      <td>1</td>
      <td>0.128673</td>
      <td>1</td>
    </tr>
    <tr>
      <th>13DAMPPOX</th>
      <td>0.216258</td>
      <td>1</td>
      <td>HGNC:549 or HGNC:550 or HGNC:80</td>
      <td>hgu133APlus2_biopsys</td>
      <td>1</td>
      <td>0.128673</td>
      <td>1</td>
    </tr>
    <tr>
      <th>14HMDZALThr</th>
      <td>0.216258</td>
      <td>0.216258</td>
      <td>HGNC:2637 or HGNC:2638</td>
      <td>hgu133APlus2_biopsys</td>
      <td>1</td>
      <td>0.128673</td>
      <td>1</td>
    </tr>
    <tr>
      <th>14HMDZhr</th>
      <td>0.216258</td>
      <td>0.216258</td>
      <td>HGNC:2637 or HGNC:2638</td>
      <td>hgu133APlus2_biopsys</td>
      <td>1</td>
      <td>0.128673</td>
      <td>1</td>
    </tr>
    <tr>
      <th>14HMDZitr</th>
      <td>0.216258</td>
      <td>0.216258</td>
      <td></td>
      <td>hgu133APlus2_biopsys</td>
      <td>1</td>
      <td>0.128673</td>
      <td>1</td>
    </tr>
    <tr>
      <th>14MDZtev</th>
      <td>0.216258</td>
      <td>0.216258</td>
      <td></td>
      <td>hgu133APlus2_biopsys</td>
      <td>1</td>
      <td>0.128673</td>
      <td>1</td>
    </tr>
    <tr>
      <th>1OHMDZhr</th>
      <td>0.216258</td>
      <td>0.216258</td>
      <td>HGNC:2637 or HGNC:2638</td>
      <td>hgu133APlus2_biopsys</td>
      <td>1</td>
      <td>0.128673</td>
      <td>1</td>
    </tr>
    <tr>
      <th>1OHMDZitr</th>
      <td>0.216258</td>
      <td>0.216258</td>
      <td></td>
      <td>hgu133APlus2_biopsys</td>
      <td>1</td>
      <td>0.128673</td>
      <td>1</td>
    </tr>
    <tr>
      <th>1OHMDZtep</th>
      <td>0.216258</td>
      <td>0.216258</td>
      <td></td>
      <td>hgu133APlus2_biopsys</td>
      <td>1</td>
      <td>0.128673</td>
      <td>1</td>
    </tr>
    <tr>
      <th>1PPDCRp</th>
      <td>0.216258</td>
      <td>1</td>
      <td></td>
      <td>hgu133APlus2_biopsys</td>
      <td>1</td>
      <td>0.128673</td>
      <td>0.128673</td>
    </tr>
    <tr>
      <th>2AMADPTm</th>
      <td>0.216258</td>
      <td>1</td>
      <td>HGNC:14411</td>
      <td>hgu133APlus2_biopsys</td>
      <td>1</td>
      <td>0.128673</td>
      <td>0.128673</td>
    </tr>
    <tr>
      <th>2HATVACIDOXDhc</th>
      <td>0.216258</td>
      <td>0.216258</td>
      <td>HGNC:2622 or HGNC:2637 or HGNC:2638</td>
      <td>hgu133APlus2_biopsys</td>
      <td>1</td>
      <td>0.128673</td>
      <td>1</td>
    </tr>
    <tr>
      <th>2HATVACIDitr</th>
      <td>0.216258</td>
      <td>0.216258</td>
      <td></td>
      <td>hgu133APlus2_biopsys</td>
      <td>1</td>
      <td>0.128673</td>
      <td>1</td>
    </tr>
    <tr>
      <th>2HATVACIDtep</th>
      <td>0.216258</td>
      <td>0.216258</td>
      <td></td>
      <td>hgu133APlus2_biopsys</td>
      <td>1</td>
      <td>0.128673</td>
      <td>1</td>
    </tr>
    <tr>
      <th>2HATVLACGLUChr</th>
      <td>0.216258</td>
      <td>0.216258</td>
      <td>HGNC:12535 or HGNC:12536</td>
      <td>hgu133APlus2_biopsys</td>
      <td>1</td>
      <td>0.128673</td>
      <td>1</td>
    </tr>
    <tr>
      <th>2HATVLACGLUCitr</th>
      <td>0.216258</td>
      <td>0.216258</td>
      <td></td>
      <td>hgu133APlus2_biopsys</td>
      <td>1</td>
      <td>0.128673</td>
      <td>1</td>
    </tr>
    <tr>
      <th>2HATVLACGLUCteb</th>
      <td>0.216258</td>
      <td>0.216258</td>
      <td>HGNC:40 or HGNC:53</td>
      <td>hgu133APlus2_biopsys</td>
      <td>1</td>
      <td>0.128673</td>
      <td>1</td>
    </tr>
    <tr>
      <th>2HATVLACOXDhc</th>
      <td>0.216258</td>
      <td>0.216258</td>
      <td>HGNC:2637 or HGNC:2638</td>
      <td>hgu133APlus2_biopsys</td>
      <td>1</td>
      <td>0.128673</td>
      <td>1</td>
    </tr>
    <tr>
      <th>2HBO</th>
      <td>0.216258</td>
      <td>1</td>
      <td>HGNC:21481 or HGNC:28335 or HGNC:30866 or (HGN...</td>
      <td>hgu133APlus2_biopsys</td>
      <td>1</td>
      <td>0.128673</td>
      <td>1</td>
    </tr>
    <tr>
      <th>2HBt2</th>
      <td>0.216258</td>
      <td>1</td>
      <td>HGNC:10922 or HGNC:10924 or HGNC:10928</td>
      <td>hgu133APlus2_biopsys</td>
      <td>1</td>
      <td>0.128673</td>
      <td>0.128673</td>
    </tr>
    <tr>
      <th>2MB2COAc</th>
      <td>0.216258</td>
      <td>1</td>
      <td>HGNC:2690</td>
      <td>hgu133APlus2_biopsys</td>
      <td>1</td>
      <td>0.128673</td>
      <td>1</td>
    </tr>
    <tr>
      <th>2OXOADOXm</th>
      <td>0.216258</td>
      <td>1</td>
      <td>HGNC:21350 and HGNC:2898 and HGNC:2911 and HGN...</td>
      <td>hgu133APlus2_biopsys</td>
      <td>1</td>
      <td>0.128673</td>
      <td>0.128673</td>
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
      <th>r2344</th>
      <td>0.216258</td>
      <td>1</td>
      <td>HGNC:18057</td>
      <td>hgu133APlus2_biopsys</td>
      <td>1</td>
      <td>0.128673</td>
      <td>0.128673</td>
    </tr>
    <tr>
      <th>r2353</th>
      <td>0.216258</td>
      <td>1</td>
      <td>HGNC:18057</td>
      <td>hgu133APlus2_biopsys</td>
      <td>1</td>
      <td>0.128673</td>
      <td>1</td>
    </tr>
    <tr>
      <th>r2361</th>
      <td>0.216258</td>
      <td>0.216258</td>
      <td>HGNC:18057</td>
      <td>hgu133APlus2_biopsys</td>
      <td>1</td>
      <td>0.128673</td>
      <td>1</td>
    </tr>
    <tr>
      <th>r2363</th>
      <td>0.216258</td>
      <td>0.216258</td>
      <td>HGNC:18057</td>
      <td>hgu133APlus2_biopsys</td>
      <td>1</td>
      <td>0.128673</td>
      <td>1</td>
    </tr>
    <tr>
      <th>r2368</th>
      <td>0.216258</td>
      <td>1</td>
      <td>HGNC:18057</td>
      <td>hgu133APlus2_biopsys</td>
      <td>1</td>
      <td>0.128673</td>
      <td>0.128673</td>
    </tr>
    <tr>
      <th>r2369</th>
      <td>0.216258</td>
      <td>1</td>
      <td>HGNC:18057</td>
      <td>hgu133APlus2_biopsys</td>
      <td>1</td>
      <td>0.128673</td>
      <td>0.128673</td>
    </tr>
    <tr>
      <th>r2372</th>
      <td>0.216258</td>
      <td>0.216258</td>
      <td>HGNC:10979</td>
      <td>hgu133APlus2_biopsys</td>
      <td>1</td>
      <td>0.128673</td>
      <td>1</td>
    </tr>
    <tr>
      <th>r2404</th>
      <td>0.216258</td>
      <td>0.216258</td>
      <td>HGNC:22921</td>
      <td>hgu133APlus2_biopsys</td>
      <td>1</td>
      <td>0.128673</td>
      <td>1</td>
    </tr>
    <tr>
      <th>r2419</th>
      <td>0.216258</td>
      <td>1</td>
      <td>HGNC:10980</td>
      <td>hgu133APlus2_biopsys</td>
      <td>1</td>
      <td>0.128673</td>
      <td>0.128673</td>
    </tr>
    <tr>
      <th>r2472</th>
      <td>0.216258</td>
      <td>1</td>
      <td>HGNC:11024</td>
      <td>hgu133APlus2_biopsys</td>
      <td>1</td>
      <td>0.128673</td>
      <td>1</td>
    </tr>
    <tr>
      <th>r2473</th>
      <td>0.216258</td>
      <td>1</td>
      <td>HGNC:4061</td>
      <td>hgu133APlus2_biopsys</td>
      <td>1</td>
      <td>0.128673</td>
      <td>0.128673</td>
    </tr>
    <tr>
      <th>r2505</th>
      <td>0.216258</td>
      <td>1</td>
      <td>HGNC:51</td>
      <td>hgu133APlus2_biopsys</td>
      <td>1</td>
      <td>0.128673</td>
      <td>0.128673</td>
    </tr>
    <tr>
      <th>r2508</th>
      <td>0.216258</td>
      <td>0.216258</td>
      <td></td>
      <td>hgu133APlus2_biopsys</td>
      <td>1</td>
      <td>0.128673</td>
      <td>1</td>
    </tr>
    <tr>
      <th>r2509</th>
      <td>0.216258</td>
      <td>1</td>
      <td></td>
      <td>hgu133APlus2_biopsys</td>
      <td>1</td>
      <td>0.128673</td>
      <td>1</td>
    </tr>
    <tr>
      <th>r2511</th>
      <td>0.216258</td>
      <td>1</td>
      <td></td>
      <td>hgu133APlus2_biopsys</td>
      <td>1</td>
      <td>0.128673</td>
      <td>1</td>
    </tr>
    <tr>
      <th>r2513</th>
      <td>0.216258</td>
      <td>1</td>
      <td></td>
      <td>hgu133APlus2_biopsys</td>
      <td>1</td>
      <td>0.128673</td>
      <td>0.128673</td>
    </tr>
    <tr>
      <th>r2514</th>
      <td>0.216258</td>
      <td>1</td>
      <td></td>
      <td>hgu133APlus2_biopsys</td>
      <td>1</td>
      <td>0.128673</td>
      <td>0.128673</td>
    </tr>
    <tr>
      <th>r2516</th>
      <td>0.216258</td>
      <td>1</td>
      <td>HGNC:10922</td>
      <td>hgu133APlus2_biopsys</td>
      <td>1</td>
      <td>0.128673</td>
      <td>0.128673</td>
    </tr>
    <tr>
      <th>r2517</th>
      <td>0.216258</td>
      <td>1</td>
      <td>HGNC:67</td>
      <td>hgu133APlus2_biopsys</td>
      <td>1</td>
      <td>0.128673</td>
      <td>0.128673</td>
    </tr>
    <tr>
      <th>r2518</th>
      <td>0.216258</td>
      <td>1</td>
      <td>HGNC:67</td>
      <td>hgu133APlus2_biopsys</td>
      <td>1</td>
      <td>0.128673</td>
      <td>0.128673</td>
    </tr>
    <tr>
      <th>r2519</th>
      <td>0.216258</td>
      <td>1</td>
      <td></td>
      <td>hgu133APlus2_biopsys</td>
      <td>1</td>
      <td>0.128673</td>
      <td>1</td>
    </tr>
    <tr>
      <th>r2521</th>
      <td>0.216258</td>
      <td>1</td>
      <td></td>
      <td>hgu133APlus2_biopsys</td>
      <td>1</td>
      <td>0.128673</td>
      <td>1</td>
    </tr>
    <tr>
      <th>r2537</th>
      <td>0.216258</td>
      <td>1</td>
      <td></td>
      <td>hgu133APlus2_biopsys</td>
      <td>1</td>
      <td>0.128673</td>
      <td>1</td>
    </tr>
    <tr>
      <th>r2538</th>
      <td>0.216258</td>
      <td>1</td>
      <td></td>
      <td>hgu133APlus2_biopsys</td>
      <td>1</td>
      <td>0.128673</td>
      <td>1</td>
    </tr>
    <tr>
      <th>sink_citr_LPAREN_c_RPAREN_</th>
      <td>0.216258</td>
      <td>1</td>
      <td></td>
      <td>hgu133APlus2_biopsys</td>
      <td>1</td>
      <td>0.128673</td>
      <td>1</td>
    </tr>
    <tr>
      <th>sink_dd2coa_LPAREN_c_RPAREN_</th>
      <td>0.216258</td>
      <td>0.216258</td>
      <td></td>
      <td>hgu133APlus2_biopsys</td>
      <td>1</td>
      <td>0.128673</td>
      <td>1</td>
    </tr>
    <tr>
      <th>sink_octdececoa_LPAREN_c_RPAREN_</th>
      <td>0.216258</td>
      <td>0.216258</td>
      <td></td>
      <td>hgu133APlus2_biopsys</td>
      <td>1</td>
      <td>0.128673</td>
      <td>1</td>
    </tr>
    <tr>
      <th>sink_pre_prot_LPAREN_r_RPAREN_</th>
      <td>0.216258</td>
      <td>1</td>
      <td></td>
      <td>hgu133APlus2_biopsys</td>
      <td>1</td>
      <td>0.128673</td>
      <td>0.128673</td>
    </tr>
    <tr>
      <th>sink_tetdece1coa_LPAREN_c_RPAREN_</th>
      <td>0.216258</td>
      <td>0.216258</td>
      <td></td>
      <td>hgu133APlus2_biopsys</td>
      <td>1</td>
      <td>0.128673</td>
      <td>1</td>
    </tr>
    <tr>
      <th>xmpt</th>
      <td>0.216258</td>
      <td>1</td>
      <td></td>
      <td>hgu133APlus2_biopsys</td>
      <td>1</td>
      <td>0.128673</td>
      <td>1</td>
    </tr>
  </tbody>
</table>
<p>3101 rows × 7 columns</p>
</div>




```python
R_biopsys.to_csv("Drugable_targets_hgu133APlus2_biopsys.tsv",sep="\t")
R_clines.to_csv("Drugable_targets_hgu133APlus2_clines.tsv",sep="\t")
```


```python
hgu133APlus2_biopsys =  R_biopsys
hgu133APlus2_clines  =  R_clines
```


```python
hgu133APlus2_biopsys_drugable = hgu133APlus2_biopsys[ (hgu133APlus2_biopsys["norm_dise_ratio"] > 1.0000001) ].sort_values(by="norm_dise_ratio", ascending=False)
hgu133APlus2_clines_drugable = hgu133APlus2_clines[ (hgu133APlus2_clines["norm_dise_ratio"] > 1.0000001) ].sort_values(by="norm_dise_ratio", ascending=False)
```


```python
hgu133APlus2_biopsys_drugable
```


```python
Recon2.reactions.PCHOLPg_hs.name


```


```python
Recon2.metabolites.pa_hs_g.name
```


```python
hgu133APlus2_clines_drugable
```

reactionlist=["ACACT1x", "ADK1", "AGPAT1", "C14STRr", "C3STKR2r", "C4STMO1r", "CDIPTr", "CEPTC", "CHLPCTD", "CHOLK", "CLS_hs", "DMATTx", "ENO", "GAPD", "GK1", "GLCt4", "GPAM_hs", "GRTTx", "HISt4", "HMGCOARc","INSTt4", "IPDDIx", "LNSTLSr", "LYStiDF","OMPDC", "ORPT", "PGK", "PGM","r0463", "r0781", "r0787", "RPE", "SMS", "SQLEr", "TKT1", "TKT2"]


```python
Recon2.reactions.GAPD.gene_reaction_rule
```


```python
Recon2.reactions.GAPD.name
```


```python
Recon2.reactions.PCHOLPg_hs.genes
```


```python
Recon2.reactions.PCHOLPr_hs.genes
```


```python
 Recon2.reactions.PCHOLP_hs.genes
```


```python
Recon2.reactions.ACACT1x.gene_reaction_rule
```


```python
# results={}
for rx in R_biopsys.index:
    rule=hgu133APlus2_scc.reactions.get_by_id(rx).gene_reaction_rule
    results[rx]={}
    results[rx]["TypeData"]="Biopsys_HGU133APlus2"
    results[rx]["norm_flux"]=R_biopsys.loc[rx]["norm_flux"]
    results[rx]["dise_flux"]=R_biopsys.loc[rx]["dise_flux"]
    results[rx]["norm_prolif_ratio"]=R_biopsys.loc[rx]["norm_prolif_ratio"]
    results[rx]["dise_prolif_ratio"]=R_biopsys.loc[rx]["dise_prolif_ratio"]
    results[rx]["norm_dise_ratio"]=R_biopsys.loc[rx]["norm_dise_ratio"]

results=pd.DataFrame(results).transpose()
```


```python
results={}
for rx in R_clines[ R_clines["norm_dise_ratio"] > 1.2 ].index:
    rule=hgu133APlus2_hela.reactions.get_by_id(rx).gene_reaction_rule
    if(rule!=''):
        results[rx]={}
        results[rx]["Reaction"]=rx
        results[rx]["TypeData"]="cLines_HGU133APlus2"
        results[rx]["geneRule"]=rule
        results[rx]["norm_flux"]=R_clines.loc[rx]["norm_flux"]
        results[rx]["dise_flux"]=R_clines.loc[rx]["dise_flux"]
        results[rx]["norm_prolif_ratio"]=R_clines.loc[rx]["norm_prolif_ratio"]*100
        results[rx]["dise_prolif_ratio"]=R_clines.loc[rx]["dise_prolif_ratio"]*100
        results[rx]["norm_dise_ratio"]=R_clines.loc[rx]["norm_dise_ratio"]

pd.DataFrame(results).transpose()

#R_biopsys[ R_biopsys["norm_dise_ratio"] >1.1 ].sort_values(by="norm_dise_ratio", ascending=False)
```


```python
FBA_hela = hgu133APlus2_hela.optimize()
FBA_hacat =hgu133APlus2_hacat.optimize()

FBA_scc = hgu133APlus2_scc.optimize()
FBA_normal = hgu133APlus2_normal.optimize()

print("Hela proliferation:",FBA_hela.f)
print("HaCaT proliferation:",FBA_hacat.f)
print("Scc proliferation:",FBA_scc.f)
print("Normal proliferation:",FBA_normal.f)



print("Hela TPI:",FBA_hela["TPI"])
print("HaCaT TPI:",FBA_hacat["TPI"])
print("Scc TPI:",FBA_scc["TPI"])
print("Normal TPI:",FBA_normal["TPI"])

```
hgu133APlus2_hela.reactions.TPI.bounds=(-.1,.1)
hgu133APlus2_hacat.reactions.TPI.bounds=(-.1,.1)
hgu133APlus2_scc.reactions.TPI.bounds=(-.1,.1)
hgu133APlus2_normal.reactions.TPI.bounds=(-.1,.1)

FBA_hela = hgu133APlus2_hela.optimize()
FBA_hacat =hgu133APlus2_hacat.optimize()

FBA_scc = hgu133APlus2_scc.optimize()
FBA_normal = hgu133APlus2_normal.optimize()

print("Hela proliferation:",FBA_hela.f)
print("HaCaT proliferation:",FBA_hacat.f)
print("Scc proliferation:",FBA_scc.f)
print("Normal proliferation:",FBA_normal.f)



print("Hela TPI:",FBA_hela["TPI"])
print("HaCaT TPI:",FBA_hacat["TPI"])
print("Scc TPI:",FBA_scc["TPI"])
print("Normal TPI:",FBA_normal["TPI"])


```python
hgu133APlus2_hela.reactions.PGM.reaction
```


```python
import cobra
Recon2 = cobra.io.read_sbml_model("Models/recon2.2.xml")

hgu133A_hela=cobra.io.read_sbml_model("hgu133A_hela_model_0418_DEMEM6429_n5.sbml")
hgu133A_hacat=cobra.io.read_sbml_model("hgu133A_hacat_model_0418_DEMEM6429_n5.sbml")

hgu133A_scc=cobra.io.read_sbml_model("hgu133A_scc_model_0418_western_diet_n5.sbml")
hgu133A_normal=cobra.io.read_sbml_model("hgu133A_nc_model_0418_western_diet_n5.sbml")

```


```python
R_clines=get_drugable_targets(hgu133A_hacat, hgu133A_hela, "hgu133A_clines",0.001)
R_biopsys=get_drugable_targets(hgu133A_normal, hgu133A_scc, "hgu133A_biopsys" ,0.001)

R_biopsys.to_csv("Drugable_targets_hgu133A_biopsys.tsv",sep="\t")
R_clines.to_csv("Drugable_targets_hgu133A_clines.tsv",sep="\t")
```


```python
hgu133APlus2_biopsys =  pd.read_csv("Drugable_targets_hgu133APlus2_biopsys.tsv",sep="\t", index_col=0)
hgu133APlus2_clines  =  pd.read_csv("Drugable_targets_hgu133APlus2_clines.tsv",sep="\t", index_col=0)
```


```python
hgu133A_biopsys      =  pd.read_csv("Drugable_targets_hgu133A_biopsys.tsv",sep="\t", index_col=0)
hgu133A_clines       =  pd.read_csv("Drugable_targets_hgu133A_clines.tsv",sep="\t", index_col=0)
```


```python

```


```python
hgu133A_clines_drugable = hgu133A_clines[ (hgu133A_clines["norm_dise_ratio"] > 1.0000001) & (hgu133A_clines["dise_flux"] < hgu133A_clines["norm_flux"])].sort_values(by="norm_dise_ratio", ascending=False)
```


```python
hgu133A_biopsys_drugable = hgu133A_biopsys[ (hgu133A_biopsys["norm_dise_ratio"] > 1.0000001) & (hgu133A_biopsys["dise_flux"] < hgu133A_biopsys["norm_flux"])].sort_values(by="norm_dise_ratio", ascending=False)
```


```python
hgu133A_clines_drugable
```


```python
import itertools
import numpy as np

a = np.array(list(itertools.chain(hgu133APlus2_biopsys_drugable.index, hgu133APlus2_clines_drugable.index, hgu133A_clines_drugable.index, hgu133A_biopsys_drugable.index)))
u  = np.unique(a, return_index=False)
```


```python
for rx in u:
    sets = Recon2.reactions.get_by_id(rx).genes
    if(list(sets)):
        print(Recon2.reactions.get_by_id(rx).id)
        print(list(sets))

```

ACACT1x	HGNC:94
ADK1	HGNC:365,HGNC:361,HGNC:20091,HGNC:362
AGPAT1	HGNC:326,HGNC:324,HGNC:20886,HGNC:20885,HGNC:25193,HGNC:325,HGNC:20880
C14STRr	HGNC:11863
C3STKR2r	HGNC:5213
C4STMO1r	HGNC:10545
CATm	HGNC:1516
CLS_hs	HGNC:16148
DMATTx	HGNC:3631
ENO	HGNC:3354,HGNC:3353,HGNC:3350
GAPD	HGNC:24864,HGNC:4141
GLCt4	HGNC:23155,HGNC:11038,HGNC:22146,HGNC:23091,HGNC:28750,HGNC:11037
GPAM_hs	HGNC:24865,HGNC:25193
GRTTx	HGNC:3631
HISt4	HGNC:13448,HGNC:11047
HMGCOARc	HGNC:5006
INSTt4	HGNC:11038
IPDDIx	HGNC:5387,HGNC:23487
LNSTLSr	HGNC:6708
LPASE	HGNC:15449,HGNC:9035,HGNC:2014
LYStiDF	HGNC:11057,HGNC:14679,HGNC:11060,HGNC:11061
OMPDC	HGNC:12563
ORPT	HGNC:12563
PGK	HGNC:8896,HGNC:8898
PGM	HGNC:1093,HGNC:8888,HGNC:8889
RPE	HGNC:10293
SMS	HGNC:29799
SQLEr	HGNC:11279
TKT1	HGNC:11835,HGNC:11834,HGNC:25313
TKT2	HGNC:11835,HGNC:11834,HGNC:25313
r0463	HGNC:5007
r0781	HGNC:2649


```python
Recon2.reactions.EX_h2o2_LPAREN_e_RPAREN_.reaction
```


```python
hgu133APlus2_hela=cobra.io.read_sbml_model("Models/hgu133APlus2_hela_model_0418_DEMEM6429_n5.sbml")
hgu133APlus2_hacat=cobra.io.read_sbml_model("Models/hgu133APlus2_hacat_model_0418_DEMEM6429_n5.sbml")

hgu133APlus2_scc=cobra.io.read_sbml_model("Models/hgu133APlus2_cancer_model_0418_westDiet_n5.sbml")
hgu133APlus2_normal=cobra.io.read_sbml_model("Models/hgu133APlus2_normal_model_0418_westDiet_n5.sbml")


```


```python
import pandas
from time import time

import cobra.test

from cobra.flux_analysis import (
    single_gene_deletion, single_reaction_deletion, double_gene_deletion,
    double_reaction_deletion)

deletion_results=single_reaction_deletion(hgu133APlus2_scc, hgu133APlus2_scc.reactions)

hgu133APlus2_scc_del = deletion_results[deletion_results.status != 'optimal']
hgu133APlus2_scc_del[hgu133APlus2_scc_del.status != 'optimal']

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

deletion_results=single_reaction_deletion(hgu133APlus2_hacat, hgu133APlus2_hacat.reactions)

hgu133APlus2_hacat_del = deletion_results[deletion_results.status != 'optimal']
hgu133APlus2_hacat_del[hgu133APlus2_hacat_del.status != 'optimal']

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
    <tr>
      <th>C14STRr</th>
      <td>0.0</td>
      <td>infeasible</td>
    </tr>
    <tr>
      <th>C4STMO1r</th>
      <td>0.0</td>
      <td>infeasible</td>
    </tr>
    <tr>
      <th>CDIPTr</th>
      <td>0.0</td>
      <td>infeasible</td>
    </tr>
    <tr>
      <th>CHSTEROLtrc</th>
      <td>0.0</td>
      <td>infeasible</td>
    </tr>
    <tr>
      <th>CLS_hs</th>
      <td>0.0</td>
      <td>infeasible</td>
    </tr>
    <tr>
      <th>CO2ter</th>
      <td>0.0</td>
      <td>infeasible</td>
    </tr>
    <tr>
      <th>DMATTx</th>
      <td>0.0</td>
      <td>infeasible</td>
    </tr>
    <tr>
      <th>DSAT</th>
      <td>0.0</td>
      <td>infeasible</td>
    </tr>
    <tr>
      <th>EX_lys_L_LPAREN_e_RPAREN_</th>
      <td>0.0</td>
      <td>infeasible</td>
    </tr>
    <tr>
      <th>EX_met_L_LPAREN_e_RPAREN_</th>
      <td>0.0</td>
      <td>infeasible</td>
    </tr>
    <tr>
      <th>EX_pglyc_hs_LPAREN_e_RPAREN_</th>
      <td>0.0</td>
      <td>infeasible</td>
    </tr>
    <tr>
      <th>EX_phe_L_LPAREN_e_RPAREN_</th>
      <td>0.0</td>
      <td>infeasible</td>
    </tr>
    <tr>
      <th>EX_pro_L_LPAREN_e_RPAREN_</th>
      <td>0.0</td>
      <td>infeasible</td>
    </tr>
    <tr>
      <th>EX_thr_L_LPAREN_e_RPAREN_</th>
      <td>0.0</td>
      <td>infeasible</td>
    </tr>
    <tr>
      <th>EX_trp_L_LPAREN_e_RPAREN_</th>
      <td>0.0</td>
      <td>infeasible</td>
    </tr>
    <tr>
      <th>FORtr</th>
      <td>0.0</td>
      <td>infeasible</td>
    </tr>
    <tr>
      <th>FRDPtr</th>
      <td>0.0</td>
      <td>infeasible</td>
    </tr>
    <tr>
      <th>GK1</th>
      <td>0.0</td>
      <td>infeasible</td>
    </tr>
    <tr>
      <th>GRTTx</th>
      <td>0.0</td>
      <td>infeasible</td>
    </tr>
    <tr>
      <th>HMGCOARc</th>
      <td>0.0</td>
      <td>infeasible</td>
    </tr>
    <tr>
      <th>IPDDIx</th>
      <td>0.0</td>
      <td>infeasible</td>
    </tr>
    <tr>
      <th>LNSTLSr</th>
      <td>0.0</td>
      <td>infeasible</td>
    </tr>
    <tr>
      <th>LYStiDF</th>
      <td>0.0</td>
      <td>infeasible</td>
    </tr>
    <tr>
      <th>O2ter</th>
      <td>0.0</td>
      <td>infeasible</td>
    </tr>
    <tr>
      <th>PGI</th>
      <td>0.0</td>
      <td>infeasible</td>
    </tr>
    <tr>
      <th>PGLYCt</th>
      <td>0.0</td>
      <td>infeasible</td>
    </tr>
    <tr>
      <th>RE2675C</th>
      <td>0.0</td>
      <td>infeasible</td>
    </tr>
    <tr>
      <th>SMS</th>
      <td>0.0</td>
      <td>infeasible</td>
    </tr>
    <tr>
      <th>SQLEr</th>
      <td>0.0</td>
      <td>infeasible</td>
    </tr>
    <tr>
      <th>biomass_DNA</th>
      <td>0.0</td>
      <td>infeasible</td>
    </tr>
    <tr>
      <th>biomass_RNA</th>
      <td>0.0</td>
      <td>infeasible</td>
    </tr>
    <tr>
      <th>biomass_carbohydrate</th>
      <td>0.0</td>
      <td>infeasible</td>
    </tr>
    <tr>
      <th>biomass_lipid</th>
      <td>0.0</td>
      <td>infeasible</td>
    </tr>
    <tr>
      <th>biomass_other</th>
      <td>0.0</td>
      <td>infeasible</td>
    </tr>
    <tr>
      <th>biomass_protein</th>
      <td>0.0</td>
      <td>infeasible</td>
    </tr>
    <tr>
      <th>r0781</th>
      <td>0.0</td>
      <td>infeasible</td>
    </tr>
  </tbody>
</table>
</div>




```python
hgu133APlus2_scc
```


```python
import numpy as np


main_list = np.setdiff1d(hgu133APlus2_hela_del.index,hgu133APlus2_hacat_del.index)

main_list

```




    array(['AGPAT1', 'C3STKR2r', 'GPAM_hs'], dtype=object)


