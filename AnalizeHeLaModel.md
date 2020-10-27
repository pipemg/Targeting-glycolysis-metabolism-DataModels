

```python
import cobra
#!wget http://www.ebi.ac.uk/biomodels-main/download?mid=MODEL1603150001 -O recon2.2.xml
Recon2 = cobra.io.read_sbml_model("Models/recon2.2.xml")
```

    cobra/io/sbml.py:235 [1;31mUserWarning[0m: M_h_x appears as a reactant and product FAOXC220200x



```python
hela_model=cobra.io.read_sbml_model(filename='hela_model_221117_DEMEM6429_n5.sbml')
```


    ---------------------------------------------------------------------------

    OSError                                   Traceback (most recent call last)

    <ipython-input-2-ef35c4c81b2b> in <module>()
    ----> 1 hela_model=cobra.io.read_sbml_model(filename='hela_model_221117_DEMEM6429_n5.sbml')
    

    /opt/conda/lib/python3.5/site-packages/cobra/io/sbml3.py in read_sbml_model(filename, number, **kwargs)
        566     if not _with_lxml:
        567         warn("Install lxml for faster SBML I/O", ImportWarning)
    --> 568     xmlfile = parse_stream(filename)
        569     xml = xmlfile.getroot()
        570     # use libsbml if not l3v1 with fbc v2


    /opt/conda/lib/python3.5/site-packages/cobra/io/sbml3.py in parse_stream(filename)
        150                 return parse(infile)
        151         else:
    --> 152             return parse(filename)
        153     except ParseError as e:
        154         raise CobraSBMLError("Malformed XML file: " + str(e))


    src/lxml/lxml.etree.pyx in lxml.etree.parse (src/lxml/lxml.etree.c:81716)()


    src/lxml/parser.pxi in lxml.etree._parseDocument (src/lxml/lxml.etree.c:118635)()


    src/lxml/parser.pxi in lxml.etree._parseDocumentFromURL (src/lxml/lxml.etree.c:118982)()


    src/lxml/parser.pxi in lxml.etree._parseDocFromFile (src/lxml/lxml.etree.c:117894)()


    src/lxml/parser.pxi in lxml.etree._BaseParser._parseDocFromFile (src/lxml/lxml.etree.c:112440)()


    src/lxml/parser.pxi in lxml.etree._ParserContext._handleParseResultDoc (src/lxml/lxml.etree.c:105896)()


    src/lxml/parser.pxi in lxml.etree._handleParseResult (src/lxml/lxml.etree.c:107604)()


    src/lxml/parser.pxi in lxml.etree._raiseParseError (src/lxml/lxml.etree.c:106415)()


    OSError: Error reading file 'hela_model_221117_DEMEM6429_n5.sbml': failed to load external entity "hela_model_221117_DEMEM6429_n5.sbml"



```python
import math
growth_rate = math.log(2)/16.2
print(growth_rate*1.20)
import math
print(math.log(2)/0.008)
```


```python
#Biomass bounds
for rxex in hela_model.exchanges:
    rxex.bounds=(-0.001,1)
    
#### AminoAcids

hela_model.reactions.EX_ala_L_LPAREN_e_RPAREN_.lower_bound=-1
hela_model.reactions.EX_arg_L_LPAREN_e_RPAREN_.lower_bound=-1
hela_model.reactions.EX_asn_L_LPAREN_e_RPAREN_.lower_bound=-1
hela_model.reactions.EX_asp_L_LPAREN_e_RPAREN_.lower_bound=-1
#hela_model.reactions.EX_cys_L_LPAREN_e_RPAREN_.lower_bound=-1
hela_model.reactions.EX_gln_L_LPAREN_e_RPAREN_.lower_bound=-1
hela_model.reactions.EX_glu_L_LPAREN_e_RPAREN_.lower_bound=-1
hela_model.reactions.EX_gly_LPAREN_e_RPAREN_.lower_bound=-1
hela_model.reactions.EX_his_L_LPAREN_e_RPAREN_.lower_bound=-1
hela_model.reactions.EX_ile_L_LPAREN_e_RPAREN_.lower_bound=-1
hela_model.reactions.EX_leu_L_LPAREN_e_RPAREN_.lower_bound=-1
hela_model.reactions.EX_lys_L_LPAREN_e_RPAREN_.lower_bound=-1 
hela_model.reactions.EX_met_L_LPAREN_e_RPAREN_.lower_bound=-1
hela_model.reactions.EX_phe_L_LPAREN_e_RPAREN_.lower_bound=-1
hela_model.reactions.EX_pro_L_LPAREN_e_RPAREN_.lower_bound=-1
hela_model.reactions.EX_ser_L_LPAREN_e_RPAREN_.lower_bound=-1
hela_model.reactions.EX_thr_L_LPAREN_e_RPAREN_.lower_bound=-1
hela_model.reactions.EX_trp_L_LPAREN_e_RPAREN_.lower_bound=-1
hela_model.reactions.EX_tyr_L_LPAREN_e_RPAREN_.lower_bound=-1
hela_model.reactions.EX_val_L_LPAREN_e_RPAREN_.lower_bound=-1

### DMEM 6429 medium
#Carbon Sources
hela_model.reactions.EX_glc_LPAREN_e_RPAREN_.lower_bound=-4.5
hela_model.reactions.EX_pyr_LPAREN_e_RPAREN_.lower_bound= -1

### DMEM 6429 medium
#Minerals, vitamins, other
#hela_model.reactions.EX_ca2_LPAREN_e_RPAREN_.lower_bound=-1
hela_model.reactions.EX_chol_LPAREN_e_RPAREN_.lower_bound=-1
#hela_model.reactions.EX_cl_LPAREN_e_RPAREN_.lower_bound=-1
#cobalt2
#cu2
hela_model.reactions.EX_fe2_LPAREN_e_RPAREN_.lower_bound=-1
#hela_model.reactions.EX_fe3_LPAREN_e_RPAREN_.lower_bound=-1
hela_model.reactions.EX_fol_LPAREN_e_RPAREN_.lower_single_gene_deletionbound=-1
#hela_model.reactions.EX_k_LPAREN_e_RPAREN_.lower_bound=-1
hela_model.reactions.EX_h2o_LPAREN_e_RPAREN_.lower_bound=-10
#h2s
hela_model.reactions.EX_inost_LPAREN_e_RPAREN_.lower_bound=-1
#magnesium
#manganase
#hela_model.reactions.EX_ncam_LPAREN_e_RPAREN_.lower_bound=-1
hela_model.reactions.EX_pi_LPAREN_e_RPAREN_.lower_bound=-10
#hela_model.reactions.EX_pnto_R_LPAREN_e_RPAREN_.lower_bound=-1
hela_model.reactions.EX_pydxn_LPAREN_e_RPAREN_.lower_bound=-1
hela_model.reactions.EX_ribflv_LPAREN_e_RPAREN_.lower_bound=-1
hela_model.reactions.EX_so4_LPAREN_e_RPAREN_.lower_bound=-1
#hela_model.reactions.EX_thm_LPAREN_e_RPAREN_.lower_bound=-1
#zn2

#EXTRA
hela_model.reactions.EX_o2_LPAREN_e_RPAREN_.bounds=(-1000,0)
#hela_model.reactions.EX_h2o2_LPAREN_e_RPAREN_.bounds=(-10,100)
#hela_model.reactions.EX_o2s_LPAREN_e_RPAREN_.bounds=(-10,100)


#hela_model.reactions.biomass_reaction.lower_bound=0.008
hela_model.reactions.biomass_reaction.bounds=(0.008,0.052)

fba_solution = hela_model.optimize()
hela_model.optimize()
```


```python
import cobra
fba_solution = hela_model.optimize()
fba_solution
```


```python
hela_model.reactions.biomass_RNA.reaction
```


```python
fba_solution.fluxes.to_csv('hela_fba_fluxes_221117_DEMEM6429_n5.csv')
```


```python
pfba_solution = cobra.flux_analysis.pfba(hela_model)
```


```python
pfba_solution.fluxes.to_csv('hela_pfba_fluxes_221117_DEMEM6429_n5.csv')
```

Doubling time of human cancerous cell lines
16.2 Hrs 
http://bionumbers.hms.harvard.edu/bionumber.aspx?&id=100685&ver=21


```python
Recon2.reactions.DM_atp
```


```python
hela_model.metabolites.get_by_id("g6p_c").summary()
```


```python
pfba_solution.fluxes.to_csv('hela_pfba_fluxes_221117_DEMEM6429_n5_nobiomassbound.csv')
```

This is to analize the different behavior of the model


```python
closedModel=hela_model.copy()
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

closedModel.reactions.r0191.bounds=(0,0)
closedModel.reactions.r0822.bounds=(0,0)
closedModel.reactions.r0821.bounds=(0,0)


closedModel.reactions.EX_o2_LPAREN_e_RPAREN_.bounds=(-1000,0.0)
closedModel.reactions.EX_h2o_LPAREN_e_RPAREN_.bounds=(-1000,1000)
closedModel.reactions.EX_co2_LPAREN_e_RPAREN_.bounds=(-0.0,1000)
closedModel.reactions.EX_glc_LPAREN_e_RPAREN_.bounds=(-1,-1)
closedModel.optimize()
```

ICDHyrm y ICDHxm tienen la misma reaccion

Se realiza glucolisis normalmente, siempre y cuando se bloqueen las reacciones r0191, r0821 y r0822


```python
cobra.flux_analysis.pfba(closedModel).fluxes.to_csv('hela_pfba_fluxes_221117_DEMEM6429_n5_glucose_aerobic.csv')
```


```python
closedModel=hela_model.copy()
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

closedModel.reactions.r0191.bounds=(0,0)
closedModel.reactions.r0822.bounds=(0,0)
closedModel.reactions.r0821.bounds=(0,0)

# test for max ATP hydrolysis flux from only o2 and the defined carbon
# source
# glucose anaerobic

closedModel.reactions.EX_o2_LPAREN_e_RPAREN_.bounds=(0,0)
closedModel.reactions.EX_h2o_LPAREN_e_RPAREN_.bounds=(-1000,1000)
closedModel.reactions.EX_co2_LPAREN_e_RPAREN_.bounds=(-0.0,1000)
closedModel.reactions.EX_glc_LPAREN_e_RPAREN_.bounds=(-1,-1)
closedModel.optimize()
```




    <Solution 2.000 at 0x7fe90d7a5860>




```python
cobra.flux_analysis.pfba(closedModel).fluxes.to_csv('hela_pfba_fluxes_221117_DEMEM6429_n5_glucose_anaerobic.csv')
```


```python
closedModel=hela_model.copy()
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




    <Solution 1.000 at 0x7fe90adf6a20>




```python
cobra.flux_analysis.pfba(closedModel).fluxes.to_csv('hela_pfba_fluxes_221117_DEMEM6429_n5_glutamine_aerobic.csv')
```


```python
closedModel=hela_model.copy()
closedModel.add_reaction(Recon2.reactions.DM_atp_c_)
closedModel.objective="DM_atp_c_"

for rx in closedModel.exchanges:
    rx.upper_bound= 100
    rx.lower_bound= -0.00000001 #-0.00474 #this is the minimum value
    

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
pfba=cobra.flux_analysis.pfba(closedModel)

```


```python
pfba.reactions.SUCOAS1m.reaction
```


```python
pfba.reactions.SUCCt2m.reaction
```


```python
pfba.metabolites.coa_m.summary()
```

La glutamina entra al ciclo de krebs por medio de alfa-ceto-glutarato y luego sigue el ciclo hasta succinato, de alli se exporta mientras que la Coenzima-A regresa a AKGDm


```python
cobra.flux_analysis.pfba(closedModel).fluxes.to_csv('hela_pfba_fluxes_221117_DEMEM6429_n5_glutamine_anaerobic.csv')
```


```python
closedModel=hela_model.copy()
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


## fru aerobic
closedModel.reactions.EX_o2_LPAREN_e_RPAREN_.bounds=(-1000,0)
closedModel.reactions.EX_h2o_LPAREN_e_RPAREN_.bounds=(-1000,1000)
closedModel.reactions.EX_co2_LPAREN_e_RPAREN_.bounds=(0,1000)
closedModel.reactions.EX_fru_LPAREN_e_RPAREN_.bounds=(-1,-1)

closedModel.optimize()
```


```python
cobra.flux_analysis.pfba(closedModel).fluxes.to_csv('hela_pfba_fluxes_221117_DEMEM6429_n5_fructose_aerobic.csv')
```


```python
closedModel=hela_model.copy()
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


## fru anaerobic
closedModel.reactions.EX_o2_LPAREN_e_RPAREN_.bounds=(0,0)
closedModel.reactions.EX_h2o_LPAREN_e_RPAREN_.bounds=(-1000,1000)
closedModel.reactions.EX_co2_LPAREN_e_RPAREN_.bounds=(0,1000)
closedModel.reactions.EX_fru_LPAREN_e_RPAREN_.bounds=(-1,-1)

closedModel.optimize()
```


```python
cobra.flux_analysis.pfba(closedModel).fluxes.to_csv('hela_pfba_fluxes_221117_DEMEM6429_n5_fructose_anaerobic.csv')
```


```python
import pandas
from time import time

import cobra.test

from cobra.flux_analysis import (
    single_gene_deletion, single_reaction_deletion, double_gene_deletion,
    double_reaction_deletion)

deletion_results=single_reaction_deletion(hela_model, hela_model.reactions)
```


```python
deletion_results[deletion_results["status"]!="optimal"]
```


```python
for rx in closedModel.reactions:
    closedModel=hela_model.copy()
    rx.bounds=(-0.001,0.001)
    fba=closedModel.optimize()
    fba2=hela_model.optimize()
    if(fba.f < fba2.f/2 ):
        print(rx.id)
```


```python
import pandas
from time import time

import cobra.test

from cobra.flux_analysis import (
    single_gene_deletion, single_reaction_deletion, double_gene_deletion,
    double_reaction_deletion)

deletion_results=single_reaction_deletion(closedModel, closedModel.reactions)

deletion_results[deletion_results.status != 'optimal']
```

double_deletation=double_reaction_deletion(
    closedModel, closedModel.exchanges, return_frame=True).round(4)

double_deletation=double_reaction_deletion(
    hela_model, hela_model.exchanges, return_frame=True).round(4)

deletion_results[deletion_results.status != 'optimal']


import pandas
from time import time

import cobra.test

from cobra.flux_analysis import (
    single_gene_deletion, single_reaction_deletion, double_gene_deletion,
    double_reaction_deletion)

deletion_results=single_reaction_deletion(hela_model, hela_model.reactions)
deletion_results[deletion_results.status != 'optimal']
%matplotlib inline
from cobra.flux_analysis import production_envelope

prod_env=production_envelope(closedModel, ["EX_glc_LPAREN_e_RPAREN_", "EX_co2_LPAREN_e_RPAREN_"], objective='ATPS4m')
prod_env[prod_env.direction == 'maximum'].plot(kind='scatter', x="EX_glc_LPAREN_e_RPAREN_", y='EX_co2_LPAREN_e_RPAREN_')
import matplotlib.backends.backend_pdf
from cobra.flux_analysis import production_envelope

pdf = matplotlib.backends.backend_pdf.PdfPages("output.pdf")
for rxex in hela_model.exchanges:
    prod_env = production_envelope(hela_model, [rxex.id, "biomass_reaction"], objective='biomass_reaction')
    plot=prod_env[prod_env.direction == 'maximum'].plot(kind='line', x='biomass_reaction', y=rxex.id)

pdf.close()        


```python

```


```python

```
