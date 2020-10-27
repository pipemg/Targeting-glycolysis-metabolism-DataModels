

```python
import cobra
Recon2 = cobra.io.read_sbml_model("Models/recon2.2.xml")

hgu133APlus2_hela=cobra.io.read_sbml_model("hgu133APlus2_hela_model_0418_DEMEM6429_n5.sbml")
hgu133APlus2_hacat=cobra.io.read_sbml_model("hgu133APlus2_hacat_model_0418_DEMEM6429_n5.sbml")

hgu133APlus2_scc=cobra.io.read_sbml_model("hgu133APlus2_cancer_model_0418_westDiet_n5.sbml")
hgu133APlus2_normal=cobra.io.read_sbml_model("hgu133APlus2_normal_model_0418_westDiet_n5.sbml")


```

    cobra/io/sbml.py:235 [1;31mUserWarning[0m: M_h_x appears as a reactant and product FAOXC220200x



```python
hgu133APlus2_scc.optimize().f
```




    0.21625829891222031




```python
hgu133APlus2_normal.optimize().f
```




    0.12867308030486935




```python
import pandas as pd

def get_antimeteabolites_targets(normal_Model, disease_Model, model_name,  eps=0.1):
 
    metNids = [m.id for m in normal_Model.metabolites]
    metDids = [m.id for m in disease_Model.metabolites]

    rxNids = [m.id for m in normal_Model.reactions]
    rxDids = [m.id for m in disease_Model.reactions]
    
    nmodel=normal_Model.copy()
    dmodel=disease_Model.copy()

    common_mets = list(set(metNids) & set(metDids))
    print("Common metabolites size",len(common_mets))    
    common_rxs = list(set(rxNids) & set(rxDids))
    print("Common reactions size",len(common_rxs))
    
    
    unique_Nmets = list(set(metNids) - set(metDids))
    print("Normal unique metabolites size",len(unique_Nmets))
    unique_Nrxs = list(set(rxNids) - set(rxDids))
    print("Normal unique reactions size",len(unique_Nrxs))
    
    unique_Dmets = list(set(metDids) - set(metNids))
    print("Disease unique metabolites size",len(unique_Dmets))
    unique_Drxs = list(set(rxDids) - set(rxNids))
    print("Disease unique reactions size",len(unique_Drxs))
    
    nflx0=normal_Model.optimize().f
    dflx0=disease_Model.optimize().f

    
    results={}

    bounds={}

    for met in common_mets:
        #print(met)
        drxs=disease_Model.metabolites.get_by_id(met).reactions
        nrxs=normal_Model.metabolites.get_by_id(met).reactions
        for rx in drxs:
            if(met in rx.products):
                dmodel.reactions.get_by_id(rx.id).upper_bound=0
            if(met in rx.reactants):
                dmodel.reactions.get_by_id(rx.id).lower_bound=0           
            
        for rx in nrxs:
            if(met in rx.products):
                nmodel.reactions.get_by_id(rx.id).upper_bound=0
            if(met in rx.reactants):
                nmodel.reactions.get_by_id(rx.id).lower_bound=0   
   
        nfba=nmodel.optimize()    
        dfba=dmodel.optimize()
        
        nflx1=nfba.f
        dflx1=dfba.f
        
        
        results[met]={}
        
        results[met]["model"]=model_name
        
        results[met]["norm_flux"]=nflx1
        results[met]["dise_flux"]=dflx1
 
        results[met]["norm_prolif_ratio"]=nflx1/nflx0
        results[met]["dise_prolif_ratio"]=dflx1/dflx0
        
        results[met]["norm_dise_ratio"]=(nflx1/nflx0)/(dflx1/dflx0)

        for rx in nrxs:
            nmodel.reactions.get_by_id(rx.id).bounds=normal_Model.reactions.get_by_id(rx.id).bounds
        for rx in drxs:
            dmodel.reactions.get_by_id(rx.id).bounds=disease_Model.reactions.get_by_id(rx.id).bounds


    return(pd.DataFrame(results).transpose())
```


```python
antimets = get_antimeteabolites_targets(hgu133APlus2_normal, hgu133APlus2_scc, "hgu133APlus2_biopsys" )
```

    Common metabolites size 5324
    Common reactions size 1963
    Normal unique metabolites size 0
    Normal unique reactions size 524
    Disease unique metabolites size 0
    Disease unique reactions size 614



```python
antimets[antimets["norm_dise_ratio"]>1]
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
  </tbody>
</table>
</div>


