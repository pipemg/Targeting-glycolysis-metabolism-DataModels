

```python
import cobra
Recon2 = cobra.io.read_sbml_model("Models/recon2.2.xml")
```

    cobra/io/sbml.py:235 [1;31mUserWarning[0m: M_h_x appears as a reactant and product FAOXC220200x



```python
Recon2.metabolites.gln_L_c
```




    <Metabolite gln_L_c at 0x7f1c42420358>




```python
model.add_reactions(Recon2.reactions.DM_atp_c_)
model.objective="DM_atp_c_"
FBA=model.optimize()

```


```python
Recon2.reactions.DM_atp_c_.reaction
```




    'atp_c + h2o_c --> adp_c + h_c + pi_c'




```python
print(FBA["DM_atp_c_"],FBA["EX_o2_LPAREN_e_RPAREN_"],FBA["EX_h2o_LPAREN_e_RPAREN_"],FBA["EX_co2_LPAREN_e_RPAREN_"], FBA["EX_glc_LPAREN_e_RPAREN_"], FBA["EX_gln_L_LPAREN_e_RPAREN_"])
```

    10.7898785227 -0.01 -0.087748864 0.090179064 -4.5 -0.01637224

