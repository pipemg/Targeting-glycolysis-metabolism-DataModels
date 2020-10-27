

```python
import cobra
#!wget http://www.ebi.ac.uk/biomodels-main/download?mid=MODEL1603150001 -O recon2.2.xml
Recon2 = cobra.io.read_sbml_model("Models/recon2.2.xml")
```

    cobra/io/sbml.py:235 [1;31mUserWarning[0m: M_h_x appears as a reactant and product FAOXC220200x



```python
hela_model=cobra.io.read_sbml_model(filename='Full_HeLa_model_1218_DEMEM6429_n5_pf2500.sbml')
kerat_model=cobra.io.read_sbml_model(filename='Full_Kerat_model_1218_DEMEM6429_n5_pf2500.sbml')
```


```python
hela_model=cobra.io.read_sbml_model(filename='Full_HeLa_model_1218_DEMEM6429_n5_pf2500.sbml')
kerat_model=cobra.io.read_sbml_model(filename='Full_Kerat_model_1218_DEMEM6429_n5_pf2500.sbml')

scc_model=cobra.io.read_sbml_model(filename='Full_SCC_model_0119_western_diet_n5_pf2500.sbml')
nor_model=cobra.io.read_sbml_model(filename='Full_Normal_model_0119_western_diet_n5_pf2500.sbml')

```


```python
for rx in scc_model.reactions:
    print(rx.id, "\t", rx.name)
```

    13DAMPPOX 	 1,3-Diaminopropane:oxygen oxidoreductase (deaminating)
    2HBO 	 2-Hydroxybutyrate:NAD+ oxidoreductase
    34DHOXPEGOX 	 3,4-Dihydroxyphenylethyleneglycol:NAD+ oxidoreductase
    3HAO 	 3-hydroxyanthranilate 3,4-dioxygenase
    3SALACBOXL 	 3-Sulfino-L-alanine carboxy-lyase
    3SALAOX 	 cysteinesulfinic acid oxidase
    3SALATAim 	 3-sulfino-alanine transaminase (irreversible), mitochondrial
    3SPYRSPm 	 3-sulfinopyruvate hydrolase (spotaneous reaction), mitochondrial
    A_MANASEly 	 alpha-mannosidase, lysosomal
    AACOAT 	 Acetoacetyl-CoA:acetate CoA-transferase
    ABTD 	 L-arabinitol 4-dehydrogenase
    ABUTD 	 Aminobutyraldehyde dehydrogenase
    ACACT10m 	 acetyl-CoA C-acetyltransferase, mitochondrial
    ACACT1x 	 acetyl-CoA C-acetyltransferase, mitochondrial
    ACACT5p 	 acetyl-CoA C-acyltransferase (decanoyl-CoA), peroxisomal
    ACACT6p 	 acetyl-CoA C-acetyltransferase (dodecanoyl), peroxisomal
    ACACT7p 	 acetyl-CoA acyltransferase (tetradecanoyl-CoA), peroxisomal
    ACCOAC 	 acetyl-CoA carboxylase
    ACGAM6PSi 	 N-acetylglucosamine-6-phosphate synthase
    ACGAMK 	 N-acetylglucosamine kinase
    ACGAMPM 	 phosphoacetylglucosamine mutase
    ACITL 	 ATP-Citrate lyase
    ACNAM9PL 	 N-Acetylneuraminate 9-phosphate pyruvate-lyase (pyruvate-phosphorylating)
    ACNAMPH 	 N-Acetylneuraminate 9-phosphate phosphohydrolase
    ACNMLr 	 N-Acetylneuraminate lyase (reversible)
    ACOAD10m 	 acyl-CoA dehydrogenase (2-methylbutanoyl-CoA), mitochondrial
    ACONT 	 aconitase
    ACONTm 	 Aconitate hydratase
    ACP1_FMN_ 	 acid phosphatase (FMN)
    ADA 	 Adenosine deaminase
    ADAe 	 Adenosine deaminase, extracellular
    ADCim 	 Acetoacetate decarboxylation (irreversible), mitochondrial
    ADK3m 	 adentylate kinase (GTP)
    ADNK1m 	 adenosine kinase, mitochondrial
    ADPRDP 	 ADPribose diphosphatase
    ADSK 	 adenylyl-sulfate kinase
    ADSL1 	 adenylosuccinate lyase
    ADSS 	 adenylosuccinate synthase
    AG13T10g 	 N-acetyllactosaminide beta-1,3-N-acetylglucosaminyltransferase, Golgi apparatus
    AG13T11g 	 N-acetyllactosaminide beta-1,3-N-acetylglucosaminyltransferase, Golgi apparatus
    AG13T12g 	 N-acetyllactosaminide beta-1,3-N-acetylglucosaminyltransferase, Golgi apparatus
    AG13T13g 	 N-acetyllactosaminide beta-1,3-N-acetylglucosaminyltransferase, Golgi apparatus
    AG13T14g 	 N-acetyllactosaminide beta-1,3-N-acetylglucosaminyltransferase, Golgi apparatus
    AG13T15g 	 N-acetyllactosaminide beta-1,3-N-acetylglucosaminyltransferase, Golgi apparatus
    AG13T16g 	 N-acetyllactosaminide beta-1,3-N-acetylglucosaminyltransferase, Golgi apparatus
    AG13T17g 	 N-acetyllactosaminide beta-1,3-N-acetylglucosaminyltransferase, Golgi apparatus
    AG13T18g 	 N-acetyllactosaminide beta-1,3-N-acetylglucosaminyltransferase, Golgi apparatus
    AG13T1g 	 N-acetyllactosaminide beta-1,3-N-acetylglucosaminyltransferase, Golgi apparatus
    AG13T2g 	 N-acetyllactosaminide beta-1,3-N-acetylglucosaminyltransferase, Golgi apparatus
    AG13T3g 	 N-acetyllactosaminide beta-1,3-N-acetylglucosaminyltransferase, Golgi apparatus
    AG13T4g 	 N-acetyllactosaminide beta-1,3-N-acetylglucosaminyltransferase, Golgi apparatus
    AG13T5g 	 N-acetyllactosaminide beta-1,3-N-acetylglucosaminyltransferase, Golgi apparatus
    AG13T6g 	 N-acetyllactosaminide beta-1,3-N-acetylglucosaminyltransferase, Golgi apparatus
    AG13T7g 	 N-acetyllactosaminide beta-1,3-N-acetylglucosaminyltransferase, Golgi apparatus
    AG13T8g 	 N-acetyllactosaminide beta-1,3-N-acetylglucosaminyltransferase, Golgi apparatus
    AG13T9g 	 N-acetyllactosaminide beta-1,3-N-acetylglucosaminyltransferase, Golgi apparatus
    AGDC 	 N-acetylglucosamine-6-phosphate deacetylase
    AGTim 	 alanine-glyoxylate transaminase (irreversible), mitochondrial
    AGTix 	 alanine-glyoxylate transaminase (irreversible), (peroxisomal)
    AHC 	 adenosylhomocysteinase
    AHEXASE2ly 	 beta-N-acetylhexosaminidase, lysosomal
    AHEXASEly 	 beta-N-acetylhexosaminidase, lysosomal
    AKR1C42 	 aldo-keto reductase family 1, member C4 (chlordecone reductase; 3-alpha hydroxysteroid dehydrogenase, type I; dihydrodiol dehydrogenase 4)
    AKR1D2 	 aldo-keto reductase family 1, member D1 (delta 4-3-ketosteroid-5-beta-reductase)
    ALASm 	 5-aminolevulinate synthase
    ALCD1 	 alcohol dehydrogenase (methanol)
    ALCD21_L 	 alcohol dehydrogenase (L-1,2-propanediol)
    ALCD22_L 	 alcohol dehydrogenase (L-lactaldehyde)
    ALCD2if 	 alcohol dehydrogenase, forward rxn (ethanol -> acetaldehyde)
    ALCD2yf 	 alcohol dehydrogenase (ethanol, NADP), forward reaction
    ALDD21 	 aldehyde dehydrogenase (pristanal, NAD)
    ALDD2x 	 aldehyde dehydrogenase (acetaldehyde, NAD)
    ALDD2xm 	 aldehyde dehydrogenase (acetylaldehyde, NAD), mitochondrial
    ALDD2y 	 aldehyde dehydrogenase (acetaldehyde, NADP)
    ALR2 	 aldose reductase (methylglyoxal)
    ALR3 	 aldose reductase (acetol)
    AMACR2p 	 alpha-methylacyl-CoA racemase (reductase)
    AMANK 	 N-acetyl-D-mannosamine kinase
    AMPDA 	 Adenosine monophosphate deaminase
    AMY1e 	 alpha-amylase, extracellular (strch1 -> strch2)
    AMY2e 	 alpha-amylase, extracellular (glygn2 -> glygn4)
    APAT2rm 	 3-Aminopropanoate:2-oxoglutarate aminotransferase (m)
    ARABR 	 arabinose reductase
    ARACHCPT1 	 carnitine acyltransferase I
    ARACHCPT2 	 carnitine acyltransferase II
    ARGSL 	 argininosuccinate lyase
    ARGSS 	 argininosuccinate synthase
    R_group_phosphotase_3 	 R group phosphotase 3
    ARTPLM3 	 R group to palmitate conversion
    ASAH1 	 N-acylsphingosine amidohydrolase
    ASCBOX 	 ascorbic acid oxidase
    ASP1DC 	 aspartate 1-decarboxylase
    ASPTAm 	 aspartate transaminase
    B_MANNASEly 	 beta-mannosidase, lysosomal
    B3GALTg 	 Beta galactosyltransferase
    B3GNT51g 	 UDP-GlcNAc:betaGal beta-1,3-N-acetylglucosaminyltransferase 5
    BAMPPALDOX 	 beta-Aminopropion aldehyde:NAD+ oxidoreductase
    BDMT_U 	 GDPmannose:chitobiosyldiphosphodolichol beta-D-mannosyltransferase (uterus)
    BETALDHxm 	 betaine-aldehyde dehydrogenase, mitochondrial
    BILIRED 	 Nad(p)h biliverdin reductase
    BPNT 	 3,5-bisphosphate nucleotidase
    BPNT2 	 3,5-bisphosphate nucleotidase (paps)
    C14STRr 	 C-14 sterol reductase
    C161CPT1 	 production of palmitoleoylcarnitine
    C204CPT1 	 carnitine C20:4 transferase
    C226CPT1 	 carnitine C22:6 transferase
    C226CPT2 	 C226 transport into the mitochondria
    C3STDH1Pr 	 C-3 sterol dehydrogenase (4-methylzymosterol)
    C3STKR2r 	 C-3 sterol keto reductase (zymosterol)
    C4STMO1r 	 C-4 sterol methyl oxidase (4,4-dimethylzymosterol)
    C4STMO2r 	 C-4 methyl sterol oxidase
    CAT2p 	 catalase A, peroxisomal (ethanol)
    CATm 	 catalase
    CATp 	 catalase A, peroxisomal
    CDIPTr 	 phosphatidylinositol synthase (Homo sapiens)
    CDS 	 phosphatidate cytidylyltransferase
    CDSm 	 phosphatidate cytidylyltransferase
    CHOLD2m 	 choline dehydrogenase (FAD acceptor), mitochondrial
    CK 	 ATP Creatine kinase
    CLS_hs 	 cardiolipin synthase (homo sapiens)
    CMPSAS 	 CMP sialic acid synthase
    CMPSASn 	 CMP sialic acid synthase, nuclear
    CORE2GTg 	 Core 2 acetylglucosaminyltransferase, Golgi apparatus
    CORE3GTg 	 Core 3 beta-GlcNAc-transferase, Golgi apparatus
    CORE4GTg 	 Core 4 beta6-GalNAc-transferase, Golgi apparatus
    CORE5GTg 	 Core 5 alpha-GalNAc-transferase, Golgi apparatus
    CORE6GTg 	 Core 6 beta-GlcNAc-transferase A, Golgi apparatus
    CORE7GTg 	 Core 7 alpha-GalNAc-transferase, Golgi apparatus
    CPPPGO 	 coproporphyrinogen oxidase (O2 required)
    CRTNsyn 	 Creatinine synthase
    CSm 	 citrate synthase
    CYSGLTH 	 Glutathione:cystine oxidoreductase
    CYSO 	 cysteine oxidase
    CYSTA 	 cysteine transaminase
    CYSTAm 	 cysteine transaminase (mitochondrial)
    CYTDK1 	 ATPcytidine 5-phosphotransferase
    CYTDK2m 	 cytidine kinase (ATP), mitochondrial
    CYTK1 	 cytidylate kinase (CMP)
    CYTK10 	 cytidylate kinase (CMP,dGTP)
    CYTK10n 	 cytidylate kinase (CMP,dGTP),nuclear
    CYTK11 	 cytidylate kinase (dCMP,dGTP)
    CYTK11n 	 cytidylate kinase (dCMP,dGTP),nuclear
    CYTK12 	 cytidylate kinase (dCMP,dCTP)
    CYTK12n 	 cytidylate kinase (dCMP,dCTP),nuclear
    CYTK13 	 cytidylate kinase (dCMP,dATP)
    CYTK13n 	 cytidylate kinase (dCMP,dATP),nuclear
    CYTK14 	 cytidylate kinase (dCMP,UTP)
    CYTK14n 	 cytidylate kinase (dCMP,CTP),nuclear
    CYTK1m 	 cytidylate kinase (CMP),mitochondrial
    CYTK1n 	 cytidylate kinase (CMP),nuclear
    CYTK2 	 cytidylate kinase (dCMP)
    CYTK2n 	 cytidylate kinase (dCMP),nuclear
    CYTK3 	 cytidylate kinase (CMP)(GTP)
    CYTK3n 	 cytidylate kinase (dCMP,CTP),nuclear
    CYTK4 	 cytidylate kinase (dCMP)(GTP)
    CYTK4n 	 cytidylate kinase (dCMP,GTP),nuclear
    CYTK5 	 cytidylate kinase (dCMP)
    CYTK5n 	 cytidylate kinase (CMP),nuclear
    CYTK6 	 cytidylate kinase (CMP,CTP)
    CYTK6n 	 cytidylate kinase (CMP,CTP),nuclear
    CYTK7 	 cytidylate kinase (CMP,UTP)
    CYTK7n 	 cytidylate kinase (CMP,UTP),nuclear
    CYTK8 	 cytidylate kinase (CMP,dATP)
    CYTK8n 	 cytidylate kinase (CMP,dATP),nuclear
    CYTK9 	 cytidylate kinase (CMP,dCTP)
    CYTK9n 	 cytidylate kinase (CMP,dCTP),nuclear
    DADA 	 Deoxyadenosine deaminase
    DADNK 	 deoxyadenosine kinase
    DAGK_hs 	 Diacylglycerol phosphate kinase (homo sapiens)
    DASCBR 	 dehydroascorbate reductase
    DCSPTN1CPT1 	 carnitine O-palmitoyltransferase
    DCSPTN1CPT2 	 carnitine transferase
    DDPGAm 	 2-dehydro-3-deoxy-phosphogluconate aldolase, mitochondrial
    DESAT16_2 	 palmitoyl-CoA desaturase (n-C16:0CoA -> n-C16:1CoA)
    DESAT18_3 	 stearoyl-CoA desaturase (n-C18:0CoA -> n-C18:1CoA)
    DESAT18_4 	 stearoyl-CoA desaturase (n-C18:0CoA -> n-C18:1CoA)
    DESAT18_5 	 stearoyl-CoA desaturase (n-C18:0CoA -> n-C18:1CoA)
    DESAT18_9 	 fatty acyl-CoA desaturase (n-C18:2CoA -> n-C18:3CoA)
    DESAT20_1 	 fatty acyl-CoA desaturase (n-C20:3CoA -> n-C20:4CoA)
    DESAT20_2 	 fatty acyl-CoA desaturase (n-C20:4CoA -> n-C20:5CoA)
    DGAT 	 diacylglycerol acyltransferase
    DGK1 	 deoxyguanylate kinase (dGMP:ATP)
    DGK2m 	 deoxyguanylate kinase (dGMP:dATP) (mitochondrial)
    DGULND 	 dehydro-L-gulonate decarboxylase
    DHCR241r 	 24-dehydrocholesterol reductase [Precursor]
    DHCR242r 	 24-dehydrocholesterol reductase [Precursor]
    DHCR243r 	 24-dehydrocholesterol reductase [Precursor]
    DHCR71r 	 7-dehydrocholesterol reductase
    DHCR72r 	 7-dehydrocholesterol reductase
    DLNLCGCPT1 	 carnitine O-palmitoyltransferase
    DLNLCGCPT2 	 carnitine transferase
    DM_atp_c_ 	 DM atp(c)
    DMATTx 	 dimethylallyltranstransferase
    DOLASNT_Uer 	 Dolichyl-diphosphooligosaccharide:protein-L-asparagine oligopolysaccharidotransferase (uterus)
    DOLDPP_Uer 	 Dolichyl-diphosphate phosphohydrolase, human (uterus)
    DOLGPP_Uer 	 Dolichyl-beta-D-glucosyl-phosphate dolichylphosphohydrolase (uterus)
    DOLPGT1_Uer 	 dolichyl-phosphate-glucose-glycolipid alpha-glucosyltransferase (uterus)
    DOLPGT2_Uer 	 dolichyl-phosphate-glucose-glycolipid alpha-glucosyltransferase (uterus)
    DOLPGT3_Uer 	 dolichyl-phosphate-glucose-glycolipid alpha-glucosyltransferase (uterus)
    DOLPMT_U 	 Dolichyl-phosphate D-mannosyltransferase (uterus)
    DOLPMT1_Uer 	 dolichyl-phosphate-mannose-glycolipid alpha-mannosyltransferase (uterus)
    DOLPMT2_Uer 	 dolichyl-phosphate-mannose-glycolipid alpha-mannosyltransferase (uterus)
    DOLPMT3_Uer 	 dolichyl-phosphate-mannose-glycolipid alpha-mannosyltransferase (uterus)
    DOLPMT4_Uer 	 dolichyl-phosphate-mannose-glycolipid alpha-mannosyltransferase (uterus)
    DPGase 	 Diphosphoglycerate phosphatase
    DPGM 	 Diphosphoglyceromutase
    DPMVDx 	 diphosphomevalonate decarboxylase
    DRPA 	 deoxyribose-phosphate aldolase
    DSAT 	 dihydrosphingosine N-acyltransferase
    DURIK1 	 deoxyuridine kinase (ATP:Deoxyuridine)
    DURIPP 	 deoxyuridine phosphorylase
    DUTPDPm 	 dUTP diphosphatase
    DUTPDPn 	 dUTP diphosphatase, nuclear
    EBP1r 	 3-beta-hydroxysteroid-delta(8),delta(7)-isomerase
    EBP2r 	 3-beta-hydroxysteroid-delta(8),delta(7)-isomerase
    ECOAH12m 	 3-hydroxyacyl-CoA dehydratase (3-hydroxyisobutyryl-CoA) (mitochondria)
    ECOAH1m 	 3-hydroxyacyl-CoA dehydratase (3-hydroxybutanoyl-CoA) (mitochondria)
    ECOAH9m 	 2-Methylprop-2-enoyl-CoA (2-Methylbut-2-enoyl-CoA), mitochondrial
    EHGLAT2m 	 L-erythro-4-Hydroxyglutamate:2-oxoglutarate aminotransferase 2, mitochondrial
    EHGLATm 	 L-erythro-4-Hydroxyglutamate:2-oxoglutarate aminotransferase, mitochondrial
    ENGASE2ly 	 endo-beta-N-acetylglucosaminidase, lysosomal
    ENGASE3ly 	 endo-beta-N-acetylglucosaminidase, lysosomal
    ENGASEly 	 endo-beta-N-acetylglucosaminidase, lysosomal
    ENO 	 enolase
    F1PGT 	 fucose-1-phosphate guanylyltransferase
    F6Tg 	 glycoprotein 6-alpha-L-fucosyltransferase
    FACOAL160i 	 fatty-acid--CoA ligase (palmitate, n-C16:0)
    FACOAL161 	 fatty-acid--CoA ligase (hexadecenoate)
    FACOAL180i 	 fatty-acid--CoA ligase (stearate, n-C18:0)
    FACOAL1812 	 fatty-acid--CoA ligase
    FACOAL184 	 fatty-acid--CoA ligase
    FACOAL2042 	 fatty-acid--CoA ligase
    FACOAL205 	 fatty-acid--CoA ligase
    FACOAL224 	 fatty-acid--CoA ligase
    FACOAL2251 	 fatty-acid--CoA ligase
    FACOAL2252 	 fatty-acid--CoA ligase
    FACOAL226 	 fatty-acid--CoA ligase
    FACOAL80i 	 fatty-acid--CoA ligase (octanoate, n-C8:0)
    FADDP 	 FAD diphosphatase
    FAEL183 	 fatty-acyl-CoA elongation (n-C18:3CoA)
    FAEL184 	 fatty-acyl-CoA elongation (n-C20:4CoA)
    FAEL204 	 fatty-acyl-CoA elongation (n-C20:4CoA)
    FAEL205 	 fatty-acyl-CoA elongation (n-C20:5CoA)
    FALDH 	 formaldehyde dehydrogenase
    FAOXC200180m 	 fatty acid beta oxidation(C20-->C18)m
    FAOXC2031836m 	 Beta oxidation of fatty acid
    FAOXC2051843m 	 Beta oxidation of fatty acid
    FAOXC2251836m 	 Beta oxidation of fatty acid
    FAOXC226205m 	 Beta oxidation of long chain fatty acid
    FBA 	 fructose-bisphosphate aldolase
    FBA2 	 D-Fructose 1-phosphate D-glyceraldehyde-3-phosphate-lyase
    FBA4 	 D-Fructose 1-phosphate D-glyceraldehyde-3-phosphate-lyase
    FBP26 	 Fructose-2,6-bisphosphate 2-phosphatase
    FCLTm 	 Ferrochelatase, mitochondrial
    FE3R2e 	 Fe(III) reduction (ascorbate)
    FK 	 Fucokinase
    FKYNH 	 N-Formyl-L-kynurenine amidohydrolase
    FPGS 	 folylpolyglutamate synthetase
    FPGS2 	 folylpolyglutamate synthetase
    FPGS3 	 folylpolyglutamate synthetase
    FPGS4 	 folylpolyglutamate synthetase (DHF)
    FPGS5 	 folylpolyglutamate  synthetase (DHF)
    FPGS6 	 folylpolyglutamate synthetase (DHF)
    FPGS7 	 folylpolyglutamate synthetase (10fthf)
    FPGS8 	 folylpolyglutamate synthetase (10fthf)
    FPGS9 	 folylpolyglutamate synthetase (10fthf)
    FTHFL 	 Formate-tetrahydrofolate ligase
    FUCASE2ly 	 alpha-fucosidase, lysosomal
    FUCASEe 	 alpha-fucosidase, extracellular
    FUCASEly 	 alpha-fucosidase, lysosomal
    FUM 	 fumarase
    FUMm 	 fumarase, mitochondrial
    FUT31g 	 fucosyltransferase 3 (galactoside 3(4)-L-fucosyltransferase, Lewis blood group included) 1
    G12MT1_U 	 Glycolipid 1,2-alpha-D-mannosyltransferase (uterus)
    G12MT2_U 	 Glycolipid 1,2-alpha-D-mannosyltransferase (uterus)
    G13MT_U 	 Glycolipid 1,3-alpha-D-mannosyltransferase (uterus)
    G14T10g 	 beta-N-acetylglucosaminylglycopeptide beta-1,4-galactosyltransferase, Golgi
    G14T11g 	 beta-N-acetylglucosaminylglycopeptide beta-1,4-galactosyltransferase, Golgi
    G14T12g 	 beta-N-acetylglucosaminylglycopeptide beta-1,4-galactosyltransferase, Golgi
    G14T13g 	 beta-N-acetylglucosaminylglycopeptide beta-1,4-galactosyltransferase, Golgi
    G14T14g 	 beta-N-acetylglucosaminylglycopeptide beta-1,4-galactosyltransferase, Golgi
    G14T15g 	 beta-N-acetylglucosaminylglycopeptide beta-1,4-galactosyltransferase, Golgi
    G14T16g 	 beta-N-acetylglucosaminylglycopeptide beta-1,4-galactosyltransferase, Golgi
    G14T17g 	 beta-N-acetylglucosaminylglycopeptide beta-1,4-galactosyltransferase, Golgi
    G14T18g 	 beta-N-acetylglucosaminylglycopeptide beta-1,4-galactosyltransferase, Golgi
    G14T19g 	 beta-N-acetylglucosaminylglycopeptide beta-1,4-galactosyltransferase, Golgi
    G14T20g 	 beta-N-acetylglucosaminylglycopeptide beta-1,4-galactosyltransferase, Golgi
    G14T21g 	 beta-N-acetylglucosaminylglycopeptide beta-1,4-galactosyltransferase, Golgi
    G14T2g 	 beta-N-acetylglucosaminylglycopeptide beta-1,4-galactosyltransferase, Golgi
    G14T3g 	 beta-N-acetylglucosaminylglycopeptide beta-1,4-galactosyltransferase, Golgi
    G14T4g 	 beta-N-acetylglucosaminylglycopeptide beta-1,4-galactosyltransferase, Golgi
    G14T5g 	 beta-N-acetylglucosaminylglycopeptide beta-1,4-galactosyltransferase, Golgi
    G14T6g 	 beta-N-acetylglucosaminylglycopeptide beta-1,4-galactosyltransferase, Golgi
    G14T7g 	 beta-N-acetylglucosaminylglycopeptide beta-1,4-galactosyltransferase, Golgi
    G14T8g 	 beta-N-acetylglucosaminylglycopeptide beta-1,4-galactosyltransferase, Golgi
    G14T9g 	 beta-N-acetylglucosaminylglycopeptide beta-1,4-galactosyltransferase, Golgi
    G14Tg 	 beta-N-acetylglucosaminylglycopeptide beta-1,4-galactosyltransferase
    G16MT_U 	 Glycolipid 1,6-alpha-D-mannosyltransferase (uterus)
    G6PDA 	 glucosamine-6-phosphate deaminase
    G6PDH2r 	 glucose 6-phosphate dehydrogenase
    GALASE10ly 	 beta-galactosidase, lysosomal
    GALASE11ly 	 beta-galactosidase, lysosomal
    GALASE12ly 	 beta-galactosidase, lysosomal
    GALASE13ly 	 beta-galactosidase, lysosomal
    GALASE14ly 	 beta-galactosidase, lysosomal
    GALASE15ly 	 beta-galactosidase, lysosomal
    GALASE16ly 	 beta-galactosidase, lysosomal
    GALASE17ly 	 beta-galactosidase, lysosomal
    GALASE18ly 	 beta-galactosidase, lysosomal
    GALASE19ly 	 beta-galactosidase, lysosomal
    GALASE1ly 	 beta-galactosidase, lysosomal
    GALASE20ly 	 beta-galactosidase, lysosomal
    GALASE3ly 	 beta-galactosidase, lysosomal
    GALASE4ly 	 beta-galactosidase, lysosomal
    GALASE5ly 	 beta-galactosidase, lysosomal
    GALASE6ly 	 beta-galactosidase, lysosomal
    GALASE7ly 	 beta-galactosidase, lysosomal
    GALASE8ly 	 beta-galactosidase, lysosomal
    GALASE9ly 	 beta-galactosidase, lysosomal
    GALC 	 Galactocerebrosidase
    GALNTg 	 GalNAc transferase, Golgi apparatus
    GALOR 	 D-Galactose:NADP+ 1-oxidoreductase
    GALU 	 UTP-glucose-1-phosphate uridylyltransferase
    GAPD 	 glyceraldehyde-3-phosphate dehydrogenase
    GASNASE2ly 	 glycosylasparaginase, lysosomal
    GASNASE3ly 	 glycosylasparaginase, lysosomal
    GASNASEly 	 glycosylasparaginase, lysosomal
    GBA 	 Glucosylceramidase
    GCALDD 	 Glycolaldehyde dehydrogenase
    GCC2cm 	 glycine-cleavage complex (lipoamide), mitochondrial
    GFUCS 	 GDP-L-fucose synthase
    GGH_10FTHF5GLUl 	 Gamma-glutamyl hydrolase (10FTHF5GLU), lysosomal
    GGH_10FTHF6GLUl 	 Gamma-glutamyl hydrolase (10FTHF6GLU), lysosomal
    GGH_10FTHF7GLUl 	 Gamma-glutamyl hydrolase (10FTHF7GLU), lysosomal
    GGH_5DHFl 	 Gamma-glutamyl hydrolase (5DHF), lysosomal
    GGH_5THFl 	 Gamma-glutamyl hydrolase (5THF), lysosomal
    GGH_6DHFl 	 Gamma-glutamyl hydrolase (6DHF), lysosomal
    GGH_6THFl 	 Gamma-glutamyl hydrolase (6THF), lysosomal
    GGH_7DHFl 	 Gamma-glutamyl hydrolase (7DHF), lysosomal
    GGH_7THFl 	 Gamma-glutamyl hydrolase (7THF), lysosomal
    GGT5r 	 Gamma-glutamyltransferase 5
    GHMT2rm 	 glycine hydroxymethyltransferase, reversible, mitochondrial
    GLCAASE8ly 	 beta-glucuronidase, lysosomal
    GLCAASE9ly 	 beta-glucuronidase, lysosomal
    GLCNACPT_U 	 UDP-GlcNAc:dolichol-phosphate GlcNAc phosphotransferase (uterus)
    GLCNACT_U 	 UDP-GlcNAc:N-acetyl-D-glucosaminyl diphosphodolichol N-acetyl-D-glucosaminyltransferase (uterus)
    GLNS 	 glutamine synthetase
    GLUCYS 	 Gamma-glutamylcysteine synthetase
    GLUDC 	 Glutamate Decarboxylase
    GLUDxm 	 glutamate dehydrogenase (NAD) (mitochondrial)
    GLUDym 	 glutamate dehydrogenase (NADP), mitochondrial
    GLUNm 	 glutaminase (mitochondrial)
    GLXO1 	 glyoxylate oxidase
    GLYCTO1p 	 Glycolate oxidase, peroxisome
    GLYOp 	 glycine oxidase, perixosomal
    GMAND 	 GDP-D-mannose dehydratase
    GMPS2 	 GMP synthase
    GPDDA1 	 Glycerophosphodiester phosphodiesterase (Glycerophosphocholine)
    GRTTx 	 geranyltranstransferase
    GTHO 	 glutathione oxidoreductase
    GTHOm 	 glutathione oxidoreductase
    GTHP 	 glutathione peroxidase
    GTHPe 	 glutathione peroxidase (e)
    GTHPm 	 glutathione peroxidase, mitochondria
    GTHS 	 glutathione synthetase
    GTPCI 	 GTP cyclohydrolase I
    GTPCIn 	 GTP cyclohydrolase I, nuclear
    GUAD 	 guanine deaminase
    GUAPRT 	 guanine phosphoribosyltransferase
    GULN3D 	 L-gulonate 3-dehydrogenase
    GULNDer 	 gulonate dehydrogenase, endoplasmic reticulum
    H2CO3D 	 carboxylic acid dissociation
    HACD1m 	 3-hydroxyacyl-CoA dehydrogenase (acetoacetyl-CoA) (mitochondria)
    HACD1x 	 3-hydroxyacyl-CoA dehydrogenase (acetoacetyl-CoA) (peroxisome)
    HACD9m 	 3-hydroxyacyl-CoA dehydrogenase (2-Methylacetoacetyl-CoA), mitochondrial
    HEX1 	 hexokinase (D-glucose:ATP)
    HEX10 	 hexokinase (D-glucosamine:ATP)
    HEX4 	 hexokinase (D-mannose:ATP)
    HEX7 	 hexokinase (D-fructose:ATP)
    HIBDm 	 3-hydroxyisobutyrate dehydrogenase, mitochondrial
    HKYNH 	 3-Hydroxy-L-kynurenine hydrolase
    HMBS 	 hydroxymethylbilane synthase
    HMGCOASi 	 Hydroxymethylglutaryl CoA synthase (ir)
    HMGLm 	 hydroxymethylglutaryl-CoA lyase
    HOXG 	 Heme oxygenase 1
    HPYRDC 	 hydroxypyruvate decarboxylase
    HPYRR2x 	 hydroxypyruvate reductase (NADH)
    HSD17B42x 	 hydroxysteroid (17-beta) dehydrogenase 4
    HSD17B4x 	 hydroxysteroid (17-beta) dehydrogenase 4
    HXPRT 	 hypoxanthine phosphoribosyltransferase (Hypoxanthine)
    HYPOE 	 hypothetical enyme
    HYPTROX 	 Hypotaurine oxidase
    ICDHxm 	 Isocitrate dehydrogenase (NAD+)
    ICDHy 	 isocitrate dehydrogenase (NADP)
    ICDHyp 	 Isocitrate dehydrogenase (NADP+)
    ICDHyrm 	 Isocitrate dehydrogenase (NADP+)
    ILETA 	 isoleucine transaminase
    IMPD 	 IMP dehydrogenase
    IPDDIx 	 isopentenyl-diphosphate D-isomerase
    KHK2 	 ketohexokinase (D-xylulose)
    KYN 	 kynureninase
    LALDO 	 D-Lactaldehyde:NAD+ oxidoreductase (glutathione-formylating)
    LCADi 	 lactaldehyde dehydrogenase
    LCADi_D 	 lactaldehyde dehydrogenase
    LCYSTATm 	 L-Cysteate:2-oxoglutarate aminotransferase, mitochondrial
    LCYSTCBOXL 	 3-Sulfoalanine carboxy-lyase
    LDH_L 	 L-lactate dehydrogenase
    LDH_Lm 	 L-lactate dehydrogenase
    LEUTA 	 leucine transaminase
    LFORKYNHYD 	 L-Formylkynurenine hydrolase
    LGTHL 	 lactoylglutathione lyase
    LNLCCPT1 	 carnitine O-palmitoyltransferase
    LNLCCPT2 	 carnitine transferase
    LNLNCGCPT1 	 carnitine O-palmitoyltransferase
    LNLNCGCPT2 	 carnitine transferase
    LNSTLSr 	 lanosterol synthase
    LPASE 	 lysophospholipase
    LSTO1r 	 Lathosterol oxidase
    LSTO2r 	 Lathosterol oxidase
    LTA4H 	 Leukotriene A-4 hydrolase
    LTC4Sr 	 Leukotriene C4 synthase
    LTD4DP 	 Leukotriene D4 dipeptidase
    M1316Mg 	 mannosyl-oligosaccharide 1,3-1,6-alpha-mannosidase
    M13N2Tg 	 alpha-1,3-mannosyl-glycoprotein 2-beta-N-acetylglucosaminyltransferase
    M14NTg 	 beta-1,4-mannosyl-glycoprotein 4-beta-N-acetylglucosaminyltransferase
    M16NTg 	 alpha-1,6-mannosyl-glycoprotein 2-beta-N-acetylglucosaminyltransferase
    MAN1PT2 	 mannose-1-phosphate guanylyltransferase (GDP)
    MCCCrm 	 methylcrotonoyl-CoA carboxylase, mitochondrial
    MCLOR 	 3-Mercaptolactate:NAD+ oxidoreductase
    MDH 	 malate dehydrogenase
    MDHm 	 malate dehydrogenase, mitochondrial
    MECOALm 	 mesaconate--CoA ligase (ADP-forming), mitochondrial
    MECOAS1m 	 mesaconate--CoA ligase (GDP-forming)
    METAT 	 methionine adenosyltransferase
    MEVK1x 	 mevalonate kinase (atp)
    MG1er 	 mannosyl-oligosaccharide glucosidase, endoplasmic reticulum
    MG2er 	 mannosyl-oligosaccharide glucosidase, endoplasmic reticulum
    MG3er 	 mannosyl-oligosaccharide glucosidase, endoplasmic reticulum
    MGCHrm 	 methylglutaconyl-CoA hydratase (reversible), mitochondrial
    MGSA 	 methylglyoxal synthase
    MGSA2 	 methyglyoxylate synthase 2 (from g3p)
    MI1345PP 	 inositol-1,3,4,5-trisphosphate 5-phosphatase
    MI134PP 	 inositol-1,3,4-trisphosphate 1-phosphatase
    MI145PK 	 inositol-1,4,5-trisphosphate 3-kinase
    MI14PP 	 inositol-1,4-bisphosphate 1-phosphatase
    MI1PP 	 myo-inositol 1-phosphatase
    MI34PP 	 inositol-3,4-bisphosphate 4-phosphatase
    MI3PP 	 myo-inositol 3-phosphatase
    MI4PP 	 myo-inositol 4-phosphatase
    MM5ag 	 mannosyl-oligosaccharide 1,2-alpha-mannosidase, Golgi apparatus
    MM6bg 	 mannosyl-oligosaccharide 1,2-alpha-mannosidase, Golgi apparatus
    MM7Cag 	 mannosyl-oligosaccharide 1,2-alpha-mannosidase, Golgi apparatus
    MM7Cbg 	 mannosyl-oligosaccharide 1,2-alpha-mannosidase, Golgi apparatus
    MM8Ag 	 mannosyl-oligosaccharide 1,2-alpha-mannosidase, Golgi apparatus
    MM8Cg 	 mannosyl-oligosaccharide 1,2-alpha-mannosidase, Golgi apparatus
    MMSAD3m 	 methylmalonate-semialdehyde dehydrogenase (malonic semialdehyde), mitochondrial
    MTHFC 	 methenyltetrahydrofolate cyclohydrolase
    MTHFCm 	 methenyltetrahydrifikate cyclohydrolase, mitochondrial
    MTHFD2 	 methylenetetrahydrofolate dehydrogenase (NAD)
    MTHFD2m 	 methylenetetrahydrofolate dehydrogenase (NAD), mitochondrial
    MTHFDm 	 methylenetetrahydrofolate dehydrogenase (NADP), mitochondrial
    N4Tg 	 N-acetylgalactosamine 4-beta-galactosyltransferase, Golgi apparatus
    NACHEX10ly 	 beta-N-acetylhexosaminidase, lysosomal
    NACHEX11ly 	 beta-N-acetylhexosaminidase, lysosomal
    NACHEX12ly 	 beta-N-acetylhexosaminidase, lysosomal
    NACHEX13ly 	 beta-N-acetylhexosaminidase, lysosomal
    NACHEX14ly 	 beta-N-acetylhexosaminidase, lysosomal
    NACHEX15ly 	 beta-N-acetylhexosaminidase, lysosomal
    NACHEX16ly 	 beta-N-acetylhexosaminidase, lysosomal
    NACHEX17ly 	 beta-N-acetylhexosaminidase, lysosomal
    NACHEX18ly 	 beta-N-acetylhexosaminidase, lysosomal
    NACHEX19ly 	 beta-N-acetylhexosaminidase, lysosomal
    NACHEX20ly 	 beta-N-acetylhexosaminidase, lysosomal
    NACHEX21ly 	 beta-N-acetylhexosaminidase, lysosomal
    NACHEX22ly 	 beta-N-acetylhexosaminidase, lysosomal
    NACHEX23ly 	 beta-N-acetylhexosaminidase, lysosomal
    NACHEX24ly 	 beta-N-acetylhexosaminidase, lysosomal
    NACHEX25ly 	 beta-N-acetylhexosaminidase, lysosomal
    NACHEX26ly 	 beta-N-acetylhexosaminidase, lysosomal
    NACHEX27ly 	 beta-N-acetylhexosaminidase, lysosomal
    NADK 	 NAD kinase
    NADN 	 NAD nucleosidase
    NADS2 	 NAD synthase (glutamine-hydrolysing)
    NAGA2ly 	 N-acetylgalactosaminidase, alpha-
    NAGLCAly 	 N-acetylglucosaminidase, lysosomal
    NBAHH_ir 	 Nalpha-(beta-alanyl)-L-histidine hydrolase IR
    NDP7er 	 nucleoside-diphosphatase (UDP), endoplasmic reticulum
    NDP7g 	 nucleoside-diphosphatase (UDP), Golgi apparatus
    NDPK3m 	 nucleoside-diphosphate kinase (ATP:CDP), mitochondrial
    NDPK4 	 nucleoside-diphosphate kinase (ATP:dTDP)
    NDPK6 	 Nucleoside-diphosphate kinase (ATP:dUDP)
    NDPK6n 	 nucleoside-diphosphate kinase (ATP:dUDP), nuclear
    NDPK9 	 nucleoside-diphosphate kinase (ATP:IDP)
    NICRNS 	 NICRNS
    NNMT 	 Nicotinamide N-methyltransferase
    NT5C 	 Nicotinate D-ribonucleotide phosphohydrolase
    NTD1 	 5-nucleotidase (dUMP)
    NTD10 	 5-nucleotidase (XMP)
    NTD11 	 5-nucleotidase (IMP)
    NTD2 	 5-nucleotidase (UMP)
    NTD3 	 5-nucleotidase (dCMP)
    NTD4 	 5-nucleotidase (CMP)
    NTD5 	 5-nucleotidase (dTMP)
    NTD6 	 5-nucleotidase (dAMP)
    NTD7 	 5-nucleotidase (AMP)
    NTD8 	 5-nucleotidase (dGMP)
    NTD9 	 5-nucleotidase (GMP)
    OIVD1m 	 2-oxoisovalerate dehydrogenase (acylating; 4-methyl-2-oxopentaoate), mitochondrial
    ORNTArm 	 ornithine transaminase reversible (m)
    P45027A13m 	 5-beta-cholestane-3-alpha,7-alpha,12-alpha-triol 27-hydroxylase
    P45027A1m 	 Cytochrome P450 27
    P4508B11r 	 sterol 12-alpha-hydroxylase
    P4508B13r 	 sterol 12-alpha-hydroxylase (nadh)
    P450LTB4r 	 cytochrome p450 leukotriene B4
    P450SCC1m 	 cholesterol monooxygenase
    PCHOLP_hs 	 choline phosphatase
    PCHOLPm_hs 	 choline phosphatase
    PDXPP 	 Pyridoxine 5-phosphate phosphatase
    PETOHMm_hs 	 phosphatidylethanolamine N-methyltransferase
    PFK 	 phosphofructokinase
    PFK26 	 6-phosphofructo-2-kinase
    PGI 	 glucose-6-phosphate isomerase
    PGK 	 phosphoglycerate kinase
    PGL 	 6-phosphogluconolactonase
    PGM 	 phosphoglycerate mutase
    PGMT 	 phosphoglucomutase
    PGS 	 Prostaglandin G/H synthase
    PHCDm 	 L-1-pyrroline-3-hydroxy-5-carboxylate dehydrogenase
    PHETA1m 	 phenylalanine transaminase (m)
    PI345P5Pn 	 phosphatidylinositol-3,4,5-trisphosphate 5-phosphatase, nuclear
    PI34P3Pn 	 phosphatidylinositol-3,4-bisphosphate 3-phosphatase, nuclear
    PI34P4Pn 	 phosphatidylinositol-3,4-bisphosphate 4-phosphatase, nuclear
    PI34P5Kn 	 phosphatidylinositol 3,4-bisphosphate 5-kinase, nuclear
    PI3P3Pn 	 phosphatidylinositol-3-phosphate 3-phosphatase, nuclear
    PI3P4Kn 	 phosphatidylinositol 3-phosphate 4-kinase, nuclear
    PI45P5P 	 phosphatidylinositol-4,5-bisphosphate 5-phosphatase
    PI45P5Pn 	 phosphatidylinositol-4,5-bisphosphate 5-phosphatase, nuclear
    PI45PLC 	 phosphatidylinositol 4,5-bisphosphate phospholipase C
    PI4P3Kn 	 phosphatidylinositol 4-phosphate 3-kinase, nuclear
    PI4P5Kn 	 phosphatidylinositol 4-phosphate 5-kinase, nuclear
    PI4PLC 	 phosphatidylinositol 4-phosphate phospholipase C
    PI4PP 	 phosphatidylinositol-4-phosphate 4-phosphatase
    PI5P4Kn 	 phosphatidylinositol-5-phosphate 4-kinase, nuclear
    PIK3n 	 phosphatidylinositol 3-kinase, nuclear
    PIK4n 	 phosphatidylinositol 4-kinase, nuclear
    PIK5n 	 phosphatidylinositol 5-kinase, nuclear
    PIPLC 	 phosphatidylinositol phospholipase C
    PMANM 	 phosphomannomutase
    PMEVKx 	 phosphomevalonate kinase
    PNP 	 purine-nucleoside phosphorylase
    PNTEH 	 Hydrolase Class ( RXN R02973)
    PNTK 	 pantothenate kinase
    PPA 	 inorganic diphosphatase
    PPAn 	 inorganic diphosphatase, nuclear
    PPAP 	 phosphatidic acid phosphatase
    PPBNGS 	 porphobilinogen synthase
    PPM 	 phosphopentomutase
    PPPGOm 	 protoporphyrinogen oxidase, mitochondrial
    PRDX 	 Peroxidase (multiple substrates)
    PRPNCOAHYDm 	 Propenoyl-CoA hydrolase (m)
    PRPPS 	 phosphoribosylpyrophosphate synthetase
    PSSA1_hs 	 Phosphatidylserine synthase homo sapiens
    PTRCOX1 	 Putrescine:oxygen oxidoreductase (deaminating)
    PUNP1 	 purine-nucleoside phosphorylase (Adenosine)
    PUNP2 	 purine-nucleoside phosphorylase (Deoxyadenosine)
    PUNP3 	 purine-nucleoside phosphorylase (Guanosine)
    PUNP4 	 purine-nucleoside phosphorylase (Deoxyguanosine)
    PUNP5 	 purine-nucleoside phosphorylase (Inosine)
    PUNP6 	 purine-nucleoside phosphorylase (Deoxyinosine)
    PYDAMK 	 pyridoxamine kinase
    PYDXK 	 pyridoxal kinase
    PYDXNK 	 pyridoxine kinase
    PYDXPP 	 Pyridoxal 5-phosphate phosphatase
    PYK 	 pyruvate kinase
    PYNP2r 	 pyrimidine-nucleoside phosphorylase (uracil)
    QUILSYN 	 Quinolinate Synthase (Eukaryotic)
    RADH2 	 retinal dehydrogenase (NADPH)
    RAI3 	 13-cis-retinoic acid isomerase
    RBFK 	 riboflavin kinase
    RDH1 	 retinol dehydrogenase (all-trans)
    RDH1a 	 retinol dehydrogenase (all-trans,NADPH)
    RDH2 	 retinol dehydrogenase (9-cis,NADH)
    RDH2a 	 retinol dehydrogenase (9-cis,NADPH)
    RDH3 	 retinol dehydrogenase (11-cis,NADH)
    RDH3a 	 retinol dehydrogenase (11-cis,NADPH)
    RETFA 	 retinol acyltransferase
    RNDR1 	 ribonucleoside-diphosphate reductase (ADP)
    RNDR2 	 ribonucleoside-diphosphate reductase (GDP)
    RNDR3 	 ribonucleoside-diphosphate reductase (CDP)
    RNDR4 	 ribonucleoside-diphosphate reductase (UDP)
    RNMK 	 ribosylnicotinamide kinase
    RPE 	 ribulose 5-phosphate 3-epimerase
    RPI 	 ribose-5-phosphate isomerase
    S23T2g 	 beta-galactoside alpha-2,3-sialyltransferase (core 2)
    S23T3g 	 beta-galactoside alpha-2,3-sialyltransferase (complex N-glycan)
    S23T4g 	 beta-galactoside alpha-2,3-sialyltransferase
    S26Tg 	 beta-galactoside alpha-2,6-sialyltransferase
    S6T10g 	 galactose/N-acetylglucosamine 6-O-sulfotransferase, Golgi apparatus
    S6T11g 	 galactose/N-acetylglucosamine 6-O-sulfotransferase, Golgi apparatus
    S6T12g 	 galactose/N-acetylglucosamine 6-O-sulfotransferase, Golgi apparatus
    S6T13g 	 galactose/N-acetylglucosamine 6-O-sulfotransferase, Golgi apparatus
    S6T14g 	 galactose/N-acetylglucosamine 6-O-sulfotransferase, Golgi apparatus
    S6T15g 	 galactose/N-acetylglucosamine 6-O-sulfotransferase, Golgi apparatus
    S6T16g 	 galactose/N-acetylglucosamine 6-O-sulfotransferase, Golgi apparatus
    S6T17g 	 galactose/N-acetylglucosamine 6-O-sulfotransferase, Golgi apparatus
    S6T18g 	 galactose/N-acetylglucosamine 6-O-sulfotransferase, Golgi apparatus
    S6T1g 	 galactose/N-acetylglucosamine 6-O-sulfotransferase, Golgi apparatus
    S6T2g 	 galactose/N-acetylglucosamine 6-O-sulfotransferase, Golgi apparatus
    S6T3g 	 galactose/N-acetylglucosamine 6-O-sulfotransferase, Golgi apparatus
    S6T4g 	 galactose/N-acetylglucosamine 6-O-sulfotransferase, Golgi apparatus
    S6T5g 	 galactose/N-acetylglucosamine 6-O-sulfotransferase, Golgi apparatus
    S6T6g 	 galactose/N-acetylglucosamine 6-O-sulfotransferase, Golgi apparatus
    S6T7g 	 galactose/N-acetylglucosamine 6-O-sulfotransferase, Golgi apparatus
    S6T8g 	 galactose/N-acetylglucosamine 6-O-sulfotransferase, Golgi apparatus
    S6T9g 	 galactose/N-acetylglucosamine 6-O-sulfotransferase, Golgi apparatus
    S6TASE10ly 	 galactose-6-sulfate sulfatase, lysosomal
    S6TASE11ly 	 N-acetylglucosamine-6-sulfatase, lysosomal
    S6TASE12ly 	 N-acetylglucosamine-6-sulfatase, lysosomal
    S6TASE13ly 	 N-acetylglucosamine-6-sulfatase, lysosomal
    S6TASE14ly 	 N-acetylglucosamine-6-sulfatase, lysosomal
    S6TASE15ly 	 N-acetylglucosamine-6-sulfatase, lysosomal
    S6TASE16ly 	 N-acetylglucosamine-6-sulfatase, lysosomal
    S6TASE17ly 	 N-acetylglucosamine-6-sulfatase, lysosomal
    S6TASE18ly 	 N-acetylglucosamine-6-sulfatase, lysosomal
    S6TASE19ly 	 N-acetylglucosamine-6-sulfatase, lysosomal
    S6TASE1ly 	 N-acetylglucosamine-6-sulfatase, lysosomal
    S6TASE20ly 	 N-acetylglucosamine-6-sulfatase, lysosomal
    S6TASE21ly 	 N-acetylglucosamine-6-sulfatase, lysosomal
    S6TASE22ly 	 galactose-6-sulfate sulfatase, lysosomal
    S6TASE23ly 	 N-acetylglucosamine-6-sulfatase, lysosomal
    S6TASE24ly 	 N-acetylglucosamine-6-sulfatase, lysosomal
    S6TASE25ly 	 galactose-6-sulfate sulfatase, lysosomal
    S6TASE26ly 	 N-acetylglucosamine-6-sulfatase, lysosomal
    S6TASE2ly 	 N-acetylglucosamine-6-sulfatase, lysosomal
    S6TASE3ly 	 N-acetylglucosamine-6-sulfatase, lysosomal
    SADT 	 sulfate adenylyltransferase
    SBTD_D2 	 sorbitol dehydrogenase
    SBTR 	 aldose reductase
    SCP2x 	 peroxisomal thiolase 2
    SCPx 	 peroxisomal thiolase 2
    SFGTH 	 S-Formylglutathione hydralase
    SIAASE2ly 	 sialidase, lysosomal
    SIAASE3ly 	 sialidase, lysosomal
    SIAASE4ly 	 sialidase, lysosomal
    SIAASEly 	 sialidase, lysosomal
    SLCBK1 	 sphingolipid long chain base kinase (sphinganine)
    SLDx 	 L-sulfolactate dehydrogenase (NAD+)
    SLDxm 	 L-sulfolactate dehydrogenase (NAD+), mitochondrial
    SMS 	 Sphingomyelin synthase (homo sapiens)
    SPHK21c 	 sphingosine kinase 2
    SPHMDAc 	 sphingomyelin deacylase
    SPMDOX 	 Spermidine:(acceptor) oxidoreductase
    SPODM 	 superoxide dismutase
    SPODMm 	 superoxide dismutase
    SPODMn 	 superoxide dismutase, nuclear
    SPODMx 	 superoxide dismutase, peroxisome
    SPTix 	 serine-pyruvate aminotransferase (irreversible), peroxisomal
    SQLEr 	 Squalene epoxidase, endoplasmic reticular (NADP)
    SQLSr 	 Squalene synthase
    SUCD1m 	 succinate dehydrogenase
    SUCOASm 	 Succinate--CoA ligase (ADP-forming)
    TALA 	 transaldolase
    TKT1 	 transketolase
    TKT2 	 transketolase
    TMDK1 	 thymidine kinase (ATP:thymidine)
    TMDS 	 thymidylate synthase
    TPI 	 triose-phosphate isomerase
    TRIOK 	 triokinase
    TRPO2 	 L-Tryptophan:oxygen 2,3-oxidoreductase (decyclizing)
    TYRCBOX 	 L-Tyrosine carboxy-lyase
    TYRTAm 	 tyrosine transaminase, mitochondrial
    UAG4E 	 UDP-N-acetylglucosamine 4-epimerase
    UAGDP 	 UDP-N-acetylglucosamine diphosphorylase
    UDPDOLPT_U 	 UDPglucose:dolichyl-phosphate beta-D-glucosyltransferase (uterus)
    UDPG4E 	 UDPglucose 4-epimerase
    UDPGD 	 UDPglucose 6-dehydrogenase
    UGCG 	 Ceramide glucosyltransferase
    UGT1A10r 	 UDP-glucuronosyltransferase 1-10 precursor, microsomal
    UMPK 	 UMP kinase
    UMPK2 	 UMP kinase (CTP)
    UMPK2n 	 UMP kinase (CTP),nuclear
    UMPK3 	 UMP kinase (UTP)
    UMPK3n 	 UMP kinase (UTP),nuclear
    UMPK4 	 UMP kinase (GTP)
    UMPK4n 	 UMP kinase (GTP),nuclear
    UMPK5 	 UMP kinase (dATP)
    UMPK5n 	 UMP kinase (dATP),nuclear
    UMPK6 	 UMP kinase (dCTP)
    UMPK6n 	 UMP kinase (dCTP),nuclear
    UMPK7 	 UMP kinase (dGTP)
    UMPK7n 	 UMP kinase (dGTP),nuclear
    UMPKn 	 UMP kinase, nuclear
    UPP3S 	 uroporphyrinogen-III synthase
    UPPDC1 	 uroporphyrinogen decarboxylase (uroporphyrinogen III)
    URIK1 	 uridine kinase (ATP:Uridine)
    VALTA 	 valine transaminase
    VLCSp 	 Very-long-chain-fatty-acid-CoA ligase
    VLCSr 	 Very-long-chain-fatty-acid-CoA ligase
    XANDp 	 xanthine dehydrogenase, peroxisomal
    XAO2x 	 xanthine oxidase
    XAOx 	 xanthine  oxidase,peroxisomal
    XYLTD_Dr 	 xylitol dehydrogenase (D-xyulose-forming)
    r0009 	 Pyrophosphate phosphohydrolase EC:3.6.1.1
    r0010 	 hydrogen-peroxide:hydrogen-peroxide oxidoreductase EC:1.11.1.6
    r0013 	 beta-N-acetylhexosaminidase Aminosugars metabolism EC:3.2.1.52
    r0016 	 Fe(II):oxygen oxidoreductase Porphyrin and chlorophyll metabolism EC:1.16.3.1
    r0021 	 glutathione:NAD+ oxidoreductase Glutamate metabolism EC:1.8.1.7
    r0022 	 glutathione:NAD+ oxidoreductase Glutamate metabolism EC:1.8.1.7
    r0047 	 Adenosine 5-monophosphate phosphohydrolase Purine metabolism EC:3.1.3.5
    r0055 	 2-Oxopropanal:NADP+ oxidoreductase Pyruvate metabolism EC:1.2.1.49
    r0082 	 Oxalosuccinate:NADP+ oxidoreductase (decarboxylating) Citrate cycle (TCA cycle) EC:1.1.1.42
    r0084 	 Oxalosuccinate:NADP+ oxidoreductase (decarboxylating) Citrate cycle (TCA cycle) EC:1.1.1.42
    r0086 	 2-Oxoglutaramate amidohydrolase Glutamate metabolism EC:3.5.1.3
    r0120 	 GTP 7,8-8,9-dihydrolase Folate biosynthesis EC:3.5.4.16
    r0121 	 r0121
    r0139 	 CDP diphosphohydrolase Pyrimidine metabolism EC:3.6.1.5
    r0149 	 CTP diphosphohydrolase Pyrimidine metabolism EC:3.6.1.5
    r0157 	 L-Glutamine:pyruvate aminotransferase Glutamate metabolism EC:2.6.1.15
    r0160 	 L-Serine:pyruvate aminotransferase Glycine, serine and threonine metabolism EC:2.6.1.51
    r0170 	 Farnesyl-diphosphate:farnesyl-diphosphate farnesyltransferase Biosynthesis of steroids EC:2.5.1.21
    r0173 	 (S)-Lactate:NAD+ oxidoreductase Glycolysis / Gluconeogenesis / Pyruvate metabolism EC:1.1.1.27
    r0226 	 5,6,7,8-Tetrahydrofolate:NADP+ oxidoreductase One carbon pool by folate / Folate biosynthesis EC:1.5.1.3
    r0239 	 N-Formylanthranilate amidohydrolase Tryptophan metabolism EC:3.5.1.9
    r0245 	 Glycerol:NADP+ oxidoreductase Glycerolipid metabolism EC:1.1.1.2 EC:1.1.1.72
    r0246 	 Glycerol:NADP+ oxidoreductase Glycerolipid metabolism EC:1.1.1.72 EC:1.1.1.2
    r0267 	 CMP-N-acetylneuraminate,ferrocytochrome-b5:oxygen oxidoreductase (N-acetyl-hydroxylating) Aminosugars metabolism EC:1.14.18.2
    r0268 	 cytidine monophospho-N-acetylneuraminic acid hydroxylase  EC:1.14.18.2
    r0283 	 L-Histidine:beta-alanine ligase (AMP-forming) Alanine and aspartate metabolism / Histidine metabolism / beta-Alanine metabolism EC:6.3.2.11
    r0287 	 Acetyl-CoA:acetyl-CoA C-acetyltransferase Fatty acid elongation in mitochondria / Fatty acid metabolism EC:2.3.1.16
    r0301 	 Xanthosine-5-phosphate:ammonia ligase (AMP-forming) Purine metabolism EC:6.3.4.1
    r0340 	 ATP:(R)-glycerate 3-phosphotransferase Glycine, serine and threonine metabolism EC:2.7.1.31
    r0345 	 ATP:AMP phosphotransferase Purine metabolism EC:2.7.4.11
    r0354 	 Glycolysis / Glyconeogenesis EC:2.7.1.1
    r0355 	 Glycolysis / Glyconeogenesis EC:2.7.1.1
    r0357 	 Glycolysis / Glyconeogenesis EC:2.7.1.1
    r0358 	 Glycolysis / Glyconeogenesis EC:2.7.1.1
    r0360 	 Glycolysis / Glyconeogenesis EC:2.7.1.1
    r0361 	 Glycolysis / Glyconeogenesis EC:2.7.1.1
    r0363 	 Glycolysis / Glyconeogenesis EC:2.7.1.1
    r0364 	 Glycolysis / Glyconeogenesis EC:2.7.1.1
    r0365 	 3-hydroxypropanoate:NAD+ oxidoreductase beta-Alanine metabolism / Propanoate metabolism EC:1.1.1.59
    r0377 	 ATP:deoxycitidine 5-phosphotransferase Pyrimidine metabolism EC:2.7.1.74
    r0381 	 hypotaurine:NAD+ oxidoreductase Taurine and hypotaurine metabolism EC:1.8.1.3
    r0383 	 pyruvate:[dihydrolipoyllysine-residue acetyltransferase]-lipoyllysine 2-oxidoreductase (decarboxylating, acceptor-acetylating) EC:1.2.4.1
    r0386 	 4-methyl-2-oxopentanoate:[dihydrolipoyllysine-residue (2-methylpropanoyl)transferase] lipoyllysine 2-oxidoreductase (decarboxylating, acceptor-2-methylpropanoylating) EC:1.2.4.4
    r0391 	 Nicotinate D-ribonucleotide:pyrophosphate phosphoribosyltransferase Nicotinate and nicotinamide metabolism EC:2.4.2.11
    r0392 	 D-Glyceraldehyde:NAD+ oxidoreductase Glycerolipid metabolism EC:1.2.1.3
    r0393 	 D-Glyceraldehyde:NAD+ oxidoreductase Glycerolipid metabolism EC:1.2.1.3
    r0394 	 Xanthine:NAD+ oxidoreductase Purine metabolism EC:1.17.1.4
    r0395 	 Hypoxanthine:oxygen oxidoreductase Purine metabolism EC:1.17.3.2
    r0400 	 N-acetylneuraminate,ferrocytochrome-b5:oxygen oxidoreductase (N-acetyl-hydroxylating) Aminosugars metabolism EC:1.14.18.2
    r0407 	 Sedoheptulose 1,7-bisphosphate D-glyceraldehyde-3-phosphate-lyase Carbon fixation EC:4.1.2.13
    r0408 	 ATP:Sedoheptulose 7-phosphate 1-phosphotransferase EC:2.7.1.11
    r0422 	 Isocitrate:NADP+ oxidoreductase (decarboxylating) Citrate cycle (TCA cycle) EC:1.1.1.42
    r0424 	 Isocitrate:NADP+ oxidoreductase (decarboxylating) Citrate cycle (TCA cycle) EC:1.1.1.42
    r0426 	 isocitrate hydro-lyase Citrate cycle (TCA cycle) / Glyoxylate and dicarboxylate metabolism EC:4.2.1.3
    r0430 	 Palmitoyl-CoA:L-carnitine O-palmitoyltransferase Fatty acid metabolism EC:2.3.1.21
    r0431 	 Palmitoyl-CoA:L-carnitine O-palmitoyltransferase Fatty acid metabolism EC:2.3.1.21
    r0441 	 Palmitoyl-CoA:L-carnitine O-palmitoyltransferase Fatty acid metabolism EC:2.3.1.21
    r0444 	 Palmitoyl-CoA:L-carnitine O-palmitoyltransferase Fatty acid metabolism EC:2.3.1.21
    r0445 	 Palmitoyl-CoA:L-carnitine O-palmitoyltransferase Fatty acid metabolism EC:2.3.1.21
    r0456 	 ATP:deoxyguanosine 5-phosphotransferase Purine metabolism EC:2.7.1.113
    r0463 	 (S)-3-Hydroxy-3-methylglutaryl-CoA acetoacetyl-CoA-lyase (CoA-acetylating) Synthesis and degradation of ketone bodies / Valine, leucine and isoleucine degradation / Butanoate metabolism EC:2.3.3.10
    r0464 	 4-Aminobutyraldehyde:NAD+ oxidoreductase Urea cycle and metabolism of amino groups EC:1.2.1.3
    r0470 	 2-Deoxyadenosine 5-diphosphate:oxidized-thioredoxin 2-oxidoreductase Purine metabolism EC:1.17.4.1
    r0472 	 2-Deoxyadenosine 5-diphosphate:oxidized-thioredoxin 2-oxidoreductase Purine metabolism EC:1.17.4.1
    r0474 	 2-Deoxyadenosine 5-diphosphate:oxidized-thioredoxin 2-oxidoreductase Purine metabolism EC:1.17.4.1
    r0475 	 2-Deoxyadenosine 5-diphosphate:oxidized-thioredoxin 2-oxidoreductase Purine metabolism EC:1.17.4.1
    r0494 	 dTDP diphosphohydrolase Pyrimidine metabolism EC:3.6.1.5
    r0497 	 dTTP nucleotidohydrolase Pyrimidine metabolism EC:3.6.1.39
    r0502 	 xanthine:NAD+ oxidoreductase Purine metabolism EC:1.17.1.4
    r0504 	 Xanthine:oxygen oxidoreductase Purine metabolism EC:1.17.3.2
    r0510 	 steroyl-CoA,hydrogen-donor:oxygen oxidoreductase Polyunsaturated fatty acid biosynthesis EC:1.14.19.1
    r0511 	 steroyl-CoA,hydrogen-donor:oxygen oxidoreductase Polyunsaturated fatty acid biosynthesis EC:1.14.19.1
    r0527 	 Nicotinamide-D-ribonucleotide amidohydrolase Nicotinate and nicotinamide metabolism EC:3.5.1.42
    r0531 	 dUTP:cytidine 5-phosphotransferase Pyrimidine metabolism EC:2.7.1.48
    r0539 	 Cysteamine:oxygen oxidoreductase Taurine and hypotaurine metabolism EC:1.13.11.19
    r0549 	 4-aminobutanal:NAD+ 1-oxidoreductase; 4-aminobutyraldehyde:NAD+ oxidoreductase Urea cycle and metabolism of amino groups / beta-Alanine metabolism EC:1.2.1.3 EC:1.2.1.19
    r0555 	 acetyl-CoA:enzyme N6-(dihydrolipoyl)lysine S-acetyltransferase Glycolysis / Gluconeogenesis / Alanine and aspartate metabolism / Pyruvate metabolism EC:2.3.1.12
    r0568 	 (5-L-Glutamyl)-L-amino-acid 5-glutamyltransferase (cyclizing) Glutathione metabolism EC:2.3.2.4
    r0570 	 2-Deoxy-D-ribose 1-phosphate 1,5-phosphomutase Pentose phosphate pathway EC:5.4.2.7
    r0575 	 Presqualene diphosphate:farnesyl-diphosphate farnesyltransferase Biosynthesis of steroids EC:2.5.1.21
    r0578 	 ATP:pantothenate 4-phosphotransferase Pantothenate and CoA biosynthesis EC:2.7.1.33
    r0579 	 ATP:pantothenate 4-phosphotransferase Pantothenate and CoA biosynthesis EC:2.7.1.33
    r0580 	 N-((R)-Pantothenoyl)-L-cysteine carboxy-lyase Pantothenate and CoA biosynthesis EC:4.1.1.30
    r0584 	 Deamino-NAD+ nucleotidohydrolase Nicotinate and nicotinamide metabolism EC:3.6.1.9
    r0595 	 3-Mercaptopyruvate:cyanide sulfurtransferase Cysteine metabolism EC:2.8.1.2
    r0596 	 3-hydroxyisobutyryl-CoA hydrolase beta-Alanine metabolism / Propanoate metabolism EC:3.1.2.4
    r0604 	 (S)-2-methylbutanoyl-CoA:enzyme N6-(dihydrolipoyl)lysine S-(2-methylbutanoyl)transferase Valine, leucine and isoleucine degradation EC:2.3.1.168
    r0614 	 N-[(R)-4-Phosphopantothenoyl]-L-cysteine carboxy-lyase Pantothenate and CoA biosynthesis EC:4.1.1.36
    r0618 	 trans-4-Hydroxy-L-proline:NADP+ 5-oxidoreductase Arginine and proline metabolism EC:1.5.1.2
    r0634 	 Octanoyl-CoA:acetyl-CoA C-acyltransferase Fatty acid elongation in mitochondria / Fatty acid metabolism EC:2.3.1.16
    r0639 	 Lauroyl-CoA:acetyl-CoA C-acyltransferase Fatty acid elongation in mitochondria / Fatty acid metabolism EC:2.3.1.16
    r0642 	 (S)-Methylmalonate semialdehyde:NAD+ oxidoreductase Valine, leucine and isoleucine degradation EC:1.2.1.3
    r0643 	 (S)-Methylmalonate semialdehyde:NAD+ oxidoreductase Valine, leucine and isoleucine degradation EC:1.2.1.3
    r0656 	 3-methylbutanoyl-CoA:enzyme N6-(dihydrolipoyl)lysine S-(3-methylbutanoyl)transferase Valine, leucine and isoleucine degradation EC:2.3.1.168
    r0660 	 (S)-3-Hydroxydodecanoyl-CoA hydro-lyase Fatty acid elongation in mitochondria / Fatty acid metabolism EC:4.2.1.17
    r0661 	 (S)-3-Hydroxydodecanoyl-CoA hydro-lyase Fatty acid elongation in mitochondria / Fatty acid metabolism EC:4.2.1.17
    r0668 	 CTP:N-acylneuraminate cytidylyltransferase Aminosugars metabolism EC:2.7.7.43
    r0670 	 (S)-3-Methyl-2-oxopentanoate:[dihydrolipoyllysine-residue (2-methylpropanoyl)transferase] lipoyllysine 2-oxidoreductase (decarboxylating, acceptor-2-methylpropanoylating) EC:1.2.4.4
    r0671 	 (R)-4-Phosphopantothenate:L-cysteine ligase EC:6.3.2.5
    r0679 	 ATP:pantothenate 4-phosphotransferase Pantothenate and CoA biosynthesis EC:2.7.1.33
    r0680 	 ATP:pantothenate 4-phosphotransferase Pantothenate and CoA biosynthesis EC:2.7.1.33
    r0683 	 Propanoyl-CoA:(acceptor) 2,3-oxidoreductase beta-Alanine metabolism EC:1.3.3.6 EC:1.3.99.3
    r0686 	 L-1-Pyrroline-3-hydroxy-5-carboxylate:NADP+ oxidoreductase Arginine and proline metabolism EC:1.5.1.12
    r0688 	 3alpha,7alpha-Dihydroxy-5beta-cholestan-26-al:NAD+ oxidoreductase Bile acid biosynthesis EC:1.2.1.3
    r0698 	 Propanoyl-CoA:acetyl-CoA C-acyltransferase Bile acid biosynthesis EC:2.3.1.16
    r0706 	 acyl-Coenzyme A dehydrogenase family, member 9  Bile acid biosynthesis EC:1.3.3.6
    r0708 	 2-Amino-4-hydroxy-6-(erythro-1,2,3-trihydroxypropyl) dihydropteridine triphosphate 7,8-8,9-dihydrolase Folate biosynthesis EC:3.5.4.16
    r0709 	 r0709
    r0715 	 (S)-3-Hydroxyhexadecanoyl-CoA:NAD+ oxidoreductase Fatty acid elongation in mitochondria / Fatty acid metabolism EC:1.1.1.35 EC:1.1.1.211
    r0717 	 (S)-3-Hydroxyhexadecanoyl-CoA hydro-lyase Fatty acid elongation in mitochondria / Fatty acid metabolism EC:4.2.1.17
    r0718 	 (S)-3-Hydroxytetradecanoyl-CoA:NAD+ oxidoreductase Fatty acid elongation in mitochondria / Fatty acid metabolism EC:1.1.1.35 EC:1.1.1.211
    r0719 	 (S)-3-Hydroxytetradecanoyl-CoA:NAD+ oxidoreductase Fatty acid elongation in mitochondria / Fatty acid metabolism EC:1.1.1.35 EC:1.1.1.211
    r0721 	 (S)-3-Hydroxytetradecanoyl-CoA hydro-lyase Fatty acid elongation in mitochondria / Fatty acid metabolism EC:4.2.1.17
    r0722 	 (S)-3-Hydroxydodecanoyl-CoA:NAD+ oxidoreductase Fatty acid elongation in mitochondria / Fatty acid metabolism EC:1.1.1.35 EC:1.1.1.211
    r0723 	 (S)-3-Hydroxydodecanoyl-CoA:NAD+ oxidoreductase Fatty acid elongation in mitochondria / Fatty acid metabolism EC:1.1.1.211 EC:1.1.1.35
    r0726 	 (S)-Hydroxydecanoyl-CoA:NAD+ oxidoreductase Fatty acid elongation in mitochondria / Fatty acid metabolism EC:1.1.1.35 EC:1.1.1.211
    r0730 	 (S)-Hydroxyoctanoyl-CoA:NAD+ oxidoreductase Fatty acid elongation in mitochondria / Fatty acid metabolism EC:1.1.1.211
    r0731 	 (S)-Hydroxyoctanoyl-CoA hydro-lyase Fatty acid elongation in mitochondria / Fatty acid metabolism EC:4.2.1.17
    r0732 	 Hexanoyl-CoA:acetyl-CoA C-acyltransferase Fatty acid elongation in mitochondria / Fatty acid metabolism EC:2.3.1.16
    r0733 	 (S)-Hydroxyhexanoyl-CoA:NAD+ oxidoreductase Fatty acid elongation in mitochondria / Fatty acid metabolism EC:1.1.1.211 EC:1.1.1.35
    r0734 	 (S)-Hydroxyhexanoyl-CoA hydro-lyase Fatty acid elongation in mitochondria / Fatty acid metabolism EC:4.2.1.17
    r0739 	 alcohol dehydrogenase Bile acid biosynthesis EC:1.1.1.1
    r0741 	 5beta-Cholestane-3alpha,7alpha,12alpha-triol,NADPH:oxygen oxidoreductase (26-hydroxylating) Bile acid biosynthesis EC:1.14.13.15
    r0743 	 hydroxysteroid (17-beta) dehydrogenase 4  Bile acid biosynthesis EC:1.1.1.35
    r0744 	 (24R,25R)-3alpha,7alpha,12alpha,24-tetrahydroxy-5beta-cholestanoyl- CoA hydro-lyase Bile acid biosynthesis EC:4.2.1.107
    r0750 	 3alpha,7alpha,12alpha-Trihydroxy-5beta-cholestane:NADP+ oxidoreductase (B-specific); 3alpha,7alpha,12alpha-Trihydroxy-5beta-cholestane:NADP+ oxidoreductase Bile acid biosynthesis EC:1.1.1.50
    r0758 	 5-Hydroxyindoleacetaldehyde:NAD+ oxidoreductase Tryptophan metabolism EC:1.2.1.3
    r0774 	 Uroporphyrinogen I carboxy-lyase Porphyrin and chlorophyll metabolism EC:4.1.1.37
    r0775 	 Formamidopyrimidine nucleoside triphosphate 7,8-8,9-dihydrolase Folate biosynthesis EC:3.5.4.16
    r0776 	 r0776
    r0777 	 GTP 7,8-8,9-dihydrolase Folate biosynthesis EC:3.5.4.16
    r0778 	 r0778
    r0781 	 Lanosterol,NADPH:oxygen oxidoreductase (14-methyl cleaving) Biosynthesis of steroids EC:1.14.13.70
    r0782 	 GDP-L-fucose:NADP+ 4-oxidoreductase (3,5-epimerizing) Fructose and mannose metabolism EC:1.1.1.271
    r0783 	 lanosterol D24-reductase Biosynthesis of steroids EC:1.3.1.72
    r0784 	 xylitol:NAD oxidoreductase Pentose and glucuronate interconversions EC:1.1.1.15
    r0787 	 3-sn-phosphatidate phosphohydrolase Sphingolipid metabolism EC:3.1.3.4
    r0795 	 Spontaneous reaction
    r1109 	 Citrate oxaloacetate-lyase ((pro-3S)-CH2COO- -->acetate) Citrate cycle (TCA cycle) EC:4.1.3.6
    r1135 	 hydroxysteroid (17-beta) dehydrogenase 7  Biosynthesis of steroids EC:1.1.1.270
    r1146 	 Biosynthesis of steroids Enzyme catalyzed
    r1167 	 EC:2.3.1.26
    r1171 	 EC:2.3.1.26
    r1177 	 EC:3.1.1.13
    r1183 	 EC:3.1.1.13
    r1378 	 S-(hydroxymethyl)glutathione dehydrogenase Methane metabolism EC:1.1.1.284
    r1380 	 delta24-sterol reductase Biosynthesis of steroids EC:1.3.1.72
    r1382 	 folylpolyglutamyl synthetase  EC:6.3.2.17
    r1383 	 gamma-glutamyl hydrolase  EC:3.4.19.9
    r1384 	 Guanosine aminohydrolase EC:3.5.4.15
    r1418 	 Carbonic acid hydro-lyase Nitrogen metabolism EC:4.2.1.1
    r1443 	 EC:1.3.3.6
    r1444 	 EC:1.3.3.6
    r1448 	 EC:1.3.3.6
    r1450 	 EC:1.3.3.6
    r1466 	 EC:2.3.1.86
    r1488 	 EC:6.2.1.3
    RE0124C 	 RE0124
    RE0383C 	 RE0383
    RE0452M 	 RE0452
    RE0453C 	 RE0453
    RE0453N 	 RE0453
    RE0512M 	 RE0512
    RE0512X 	 RE0512
    RE0565C 	 RE0565
    RE0566C 	 RE0566
    RE0567C 	 RE0567
    RE0568C 	 RE0568
    RE0583C 	 RE0583
    RE0688C 	 RE0688
    RE0691C 	 RE0691
    RE0827C 	 RE0827
    RE0908C 	 RE0908
    RE0912C 	 RE0912
    RE1062M 	 RE1062
    RE1447N 	 RE1447
    RE1448N 	 RE1448
    RE1516M 	 RE1516
    RE1517M 	 RE1517
    RE1518M 	 RE1518
    RE1520M 	 RE1520
    RE1521M 	 RE1521
    RE1522M 	 RE1522
    RE1523M 	 RE1523
    RE1525M 	 RE1525
    RE1525X 	 RE1525
    RE1526M 	 RE1526
    RE1526X 	 RE1526
    RE1527M 	 RE1527
    RE1527X 	 RE1527
    RE1530C 	 RE1530
    RE1531M 	 RE1531
    RE1532M 	 RE1532
    RE1533M 	 RE1533
    RE1534M 	 RE1534
    RE1573M 	 RE1573
    RE1573X 	 RE1573
    RE1635M 	 RE1635
    RE1635R 	 RE1635
    RE1796M 	 RE1796
    RE1807C 	 RE1807
    RE1815M 	 RE1815
    RE1815R 	 RE1815
    RE1816M 	 RE1816
    RE1816R 	 RE1816
    RE1817M 	 RE1817
    RE1817R 	 RE1817
    RE1818M 	 RE1818
    RE1818R 	 RE1818
    RE1834C 	 RE1834
    RE1897C 	 RE1897
    RE1898C 	 RE1898
    RE2051C 	 RE2051
    RE2051G 	 RE2051
    RE2051R 	 RE2051
    RE2272L 	 RE2272
    RE2318M 	 RE2318
    RE2318R 	 RE2318
    RE2319M 	 RE2319
    RE2319R 	 RE2319
    RE2349C 	 RE2349
    RE2443M 	 RE2443
    RE2514C 	 RE2514
    RE2514L 	 RE2514
    RE2541L 	 RE2541
    RE2625C 	 RE2625
    RE2626C 	 RE2626
    RE2626M 	 RE2626
    RE2649C 	 RE2649
    RE2675C 	 RE2675
    RE2677C 	 RE2677
    RE2680C 	 RE2680
    RE2814M 	 RE2814
    RE2814R 	 RE2814
    RE2908M 	 RE2908
    RE2908X 	 RE2908
    RE2909M 	 RE2909
    RE2909X 	 RE2909
    RE2910M 	 RE2910
    RE2910X 	 RE2910
    RE2916M 	 RE2916
    RE2917M 	 RE2917
    RE2919M 	 RE2919
    RE2920M 	 RE2920
    RE2921M 	 RE2921
    RE2954C 	 RE2954
    RE2987X 	 RE2987
    RE2991X 	 RE2991
    RE2998M 	 RE2998
    RE3001M 	 RE3001
    RE3004M 	 RE3004
    RE3005M 	 RE3005
    RE3006M 	 RE3006
    RE3009C 	 RE3009
    RE3011R 	 RE3011
    RE3012C 	 RE3012
    RE3012M 	 RE3012
    RE3012R 	 RE3012
    RE3076X 	 RE3076
    RE3082X 	 RE3082
    RE3088X 	 RE3088
    RE3093X 	 RE3093
    RE3106C 	 RE3106
    RE3110C 	 RE3110
    RE3112C 	 RE3112
    RE3113C 	 RE3113
    RE3119C 	 RE3119
    RE3120C 	 RE3120
    RE3121C 	 RE3121
    RE3122C 	 RE3122
    RE3136C 	 RE3136
    RE3141X 	 RE3141
    RE3156X 	 RE3156
    RE3177M 	 RE3177
    RE3178M 	 RE3178
    RE3179M 	 RE3179
    RE3189M 	 RE3189
    RE3190M 	 RE3190
    RE3191M 	 RE3191
    RE3192M 	 RE3192
    RE3193M 	 RE3193
    RE3194M 	 RE3194
    RE3224C 	 RE3224
    RE3225C 	 RE3225
    RE3226C 	 RE3226
    RE3227C 	 RE3227
    RE3228C 	 RE3228
    RE3229C 	 RE3229
    RE3230C 	 RE3230
    RE3231C 	 RE3231
    RE3232C 	 RE3232
    RE3234C 	 RE3234
    RE3235C 	 RE3235
    RE3236C 	 RE3236
    RE3237C 	 RE3237
    RE3241C 	 RE3241
    RE3241R 	 RE3241
    RE3242C 	 RE3242
    RE3242R 	 RE3242
    RE3243C 	 RE3243
    RE3243R 	 RE3243
    RE3244C 	 RE3244
    RE3244R 	 RE3244
    RE3245C 	 RE3245
    RE3248M 	 RE3248
    RE3248X 	 RE3248
    RE3250M 	 RE3250
    RE3250X 	 RE3250
    RE3272N 	 RE3272
    RE3336M 	 RE3336
    RE3337M 	 RE3337
    RE3338M 	 RE3338
    RE3338X 	 RE3338
    RE3340M 	 RE3340
    RE3340X 	 RE3340
    RE3342M 	 RE3342
    RE3344M 	 RE3344
    RE3346C 	 RE3346
    RE3346M 	 RE3346
    RE3347C 	 RE3347
    RE3381L 	 RE3381
    RE3383M 	 RE3383
    RE3384M 	 RE3384
    RE3386M 	 RE3386
    RE3387M 	 RE3387
    RE3388M 	 RE3388
    RE3390M 	 RE3390
    RE3391M 	 RE3391
    RE3392M 	 RE3392
    RE3393M 	 RE3393
    RE3394M 	 RE3394
    RE3396M 	 RE3396
    RE3398M 	 RE3398
    RE3399M 	 RE3399
    RE3400M 	 RE3400
    RE3401M 	 RE3401
    RE3402M 	 RE3402
    RE3403M 	 RE3403
    RE3404M 	 RE3404
    RE3443M 	 RE3443
    RE3447M 	 RE3447
    RE3448M 	 RE3448
    RE3448X 	 RE3448
    RE3496N 	 RE3496
    RE3521M 	 RE3521
    RE3521R 	 RE3521
    RE3532M 	 RE3532
    RE3532R 	 RE3532
    RE3533M 	 RE3533
    RE3533R 	 RE3533
    RE3534M 	 RE3534
    RE3534R 	 RE3534
    RE3554M 	 RE3554
    RE3554R 	 RE3554
    RE3557M 	 RE3557
    RE3557R 	 RE3557
    RE3559M 	 RE3559
    RE3563M 	 RE3563
    RE3564M 	 RE3564
    RE3564X 	 RE3564
    RE3572X 	 RE3572
    RE3573X 	 RE3573
    RE3587N 	 RE3587
    RE3626M 	 RE3626
    C10OHc 	 production of 3-OHdecanoylcarnitinec
    C120CPT1 	 production of dodecanoylcarnitine
    C12OHc 	 production of 3-OHdodecanoylcarnitinec
    C142OHc 	 production of 3-hydroxytetradeca dienoyl carnitine
    C14OHc 	 production of 3-hydroxytetradecanoylcarnitine
    C18OHc 	 production of 3-hydroxyoctadecanoylcarnitine
    C50CPT1 	 production of isovalerylcarnitine
    C51CPT1 	 production of tiglylcarnitine
    C60CPT1 	 production of hexanoylcarnitine
    C80CPT1 	 production of octanoylcarnitine
    FAOXC101C8m 	 fatty acid beta oxidation(C10:1-->C8)m
    FAOXC101m 	 isomerization(C10:1)m
    FAOXC102C101m 	 reduction(C10:2-->C10:1)m
    FAOXC102C103m 	 fatty acid beta oxidation(C10:2-->C10:3)m
    FAOXC102C81m 	 fatty acid beta oxidation(C10:2-->C8:1)m
    FAOXC102m 	 isomerization(C10:2)m
    FAOXC103C102m 	 fatty acid beta oxidation(C10:3-->C10:2)m
    FAOXC10C10OHm 	 fatty acid beta oxidation(C10-->OHC10)m
    FAOXC10DCC8DCx 	 fatty acid beta oxidation(C10DC-->C8DC)
    FAOXC122C101m 	 fatty acid beta oxidation(C12:2-->C10:1)m
    FAOXC122m 	 isomerization(C12:2)m
    FAOXC123C102m 	 fatty acid beta oxidation(C12:3-->C10:2)m
    FAOXC123m 	 isomerization(C12:3)m
    FAOXC12C12OHm 	 fatty acid beta oxidation(C12-->OHC12)m
    FAOXC12DCC10DCx 	 fatty acid beta oxidation(C12DC-->C10DC)x
    FAOXC142C142OHm 	 fatty acid beta oxidation(C14:2-->C14:2OH)m
    FAOXC143C123m 	 fatty acid beta oxidation(C14:3-->C12:3)m
    FAOXC14C14OHm 	 fatty acid beta oxidation(C14-->C14OH)m
    FAOXC14DCC12DCx 	 fatty acid beta oxidation(C14DC-->C12DC)x
    FAOXC162C142m 	 fatty acid beta oxidation(C16:2-->C14:2)m
    FAOXC163C164Gm 	 fatty acid beta oxidation(C16:3-->C16:4)gm
    FAOXC163GC142m 	 fatty acid beta oxidation(C16:3-->C14:2)gm
    FAOXC163Gm 	 isomerization(C16:3)gm
    FAOXC164C143m 	 fatty acid beta oxidation(C16:4-->C14:3)m
    FAOXC164C165m 	 fatty acid beta oxidation(C16:4-->C16:5)m
    FAOXC164GC163m 	 fatty acid beta oxidation(C16:4-->C16:3)gm
    FAOXC164m 	 isomerization(C16:4)m
    FAOXC165C164m 	 fatty acid beta oxidation (C16:5-->C16:4)m
    FAOXC16DCC14DCx 	 fatty acid beta oxidation(C16DC-->C14DC)x
    FAOXC16DCr 	 fatty acid activation(C16DC)er
    FAOXC181C161m 	 fatty acid beta oxidation(C18:1-->C16:1)m
    FAOXC182C162m 	 fatty acid beta oxidation(C18:2-->C16:2)m
    FAOXC183C163Gm 	 fatty acid beta oxidation(C18:3-->C16:3)gm
    FAOXC184C163m 	 fatty acid beta oxidation (C18:4-->C16:3)m
    FAOXC184C164m 	 fatty acid beta oxidation(C18:4-->C16:4)m
    FAOXC184m 	 isomerization of (C18:4)m
    FAOXC18C18OHm 	 fatty acid beta oxidation(C18-->C18OH)m
    FAOXC204C184m 	 fatty acid beta oxidation(C20:4-->C18:4)m
    FAOXC225C204m 	 fatty acid beta oxidation(C22:5-->C20:4)m
    FAOXC225C226m 	 fatty acid beta oxidation(C22:5-->C22:6)m
    FAOXC225m 	 isomerization(C22:5)m
    FAOXC226C205m 	 fatty acid beta oxidation(C22:6-->C20:5)m
    FAOXC226C225m 	 fatty acid beta oxidation(C22:6-->C22:5)m
    FAOXC226C227m 	 fatty acid beta oxidation(C22:6-->C22:7)m
    FAOXC226m 	 isomerization(C22:6)m
    FAOXC227C226m 	 fatty acid beta oxidation(C22:7-->C22:6)m
    FAOXC5C5OHm 	 fatty acid beta oxidation(C5-->C5OH)m
    FAOXC5OHc 	 (3-hydroxyisovalerylcoa-->3-hydroxyisovalerylcarnitine)c
    FAOXC61m 	 isomerization(C6:1)m
    FAOXC6DCC4DCx 	 fatty acid beta oxidation(C6DC-->C4DC)x
    FAOXC81C61m 	 fatty acid beta oxidation(C8:1-->C6:1)m
    FAOXC8DCC6DCx 	 fatty acid beta oxidation(C8DC-->C6DC)x
    FAOXOHC16C16DCc 	 fatty acid omega oxidation(w-OHC16-->C16DC)c
    FAOXTC101TC102m 	 fatty acid beta oxidation trans(C10:1-->C10:2)m
    FAOXTC102C101m 	 fatty acid beta oxidation trans(C10:2-->C10:1)m
    FAOXTC122TC101m 	 fatty acid beta oxidation trans(C12:2-->C10:1)m
    FAOXTC122m 	 isomerization trans(C12:2)m
    FAOXTC142TC122m 	 fatty acid beta oxidation trans(C14:2-->C12:2)m
    OCD11COACPT1 	 transport of 11-octadecenoyl carnitine into the mitochondrial matrix
    OCD11CRNCPT2 	 transport of 11-octadecenoyl carnitine into the mitochondrial matrix
    SUCCOAPET 	 thioesterification of succinyl coa for release into cytosol
    ALAALACNc 	 ALAALACNc
    GLYGLYCNc 	 hydrolysis of glycylycine for uptake
    GLYLEUHYDROc 	 hydrolysis of Glycylleucine in the small intestine for cellular uptake
    GLYSARCNc 	 Hydrolysis of glycylsarcosine for uptake
    OAADC 	 oxaloacetate decarboxylase
    ADNK3 	 adenosine kinase
    ADNK4 	 adenosine kinase
    DPMVDc 	 diphosphomevalonate decarboxylase, cytosol
    DSREDUCr 	 Desmosterol reductase
    DUTPDP 	 dUTP diphosphatase
    FERO 	 ferroxidase
    FMNALKPle 	 FMNALKPle
    GLYC3PFADm 	 glycerophosphate shuttle for trasnfer of reducing equivalents
    HMGCOARc 	 Hydroxymethylglutaryl CoA reductase (ir) in cytosol
    INSK 	 insosine kinase
    MEVK1c 	 mevalonate kinase (atp) cytosol
    PMEVKc 	 phosphomevalonate kinase, cytosol
    SFCYSc 	 Formation of sulfocysteine
    NNDPR 	 NNDPR
    NADKm 	 NADKm
    10FTHF5GLUtl 	 5-glutamyl-10FTHF transport, lysosomal
    10FTHF6GLUtl 	 6-glutamyl-10FTHF transport, lysosomal
    10FTHF7GLUtl 	 7-glutamyl-10FTHF transport, lysosomal
    10FTHFtl 	 10-Formyltetrahydrofolate lysosomal transport via diffusion
    1MNCAMti 	 N1-Methylnicotinamide transport
    3MOPt2im 	 3-Methyl-2-oxopentanoate mitochondrial transport via proton symport
    3SALAASPm 	 cysteinesulfinate-aspartate mitochondrial shuttle
    4MOPt2im 	 4-methyl-2-oxopentanoate mitochondrial transport  via proton symport
    5AOPtm 	 5-Aminolevulinate mitochondrial transport
    5DHFtl 	 5-glutamyl-DHF transport, lysosomal
    5THFtl 	 5-glutamyl-THF transport, lysosomal
    6DHFtl 	 6-glutamyl-DHF transport, lysosomal
    6THFtl 	 6-glutamyl-THF transport, lysosomal
    7DHFtl 	 7-glutamyl-DHF transport, lysosomal
    7THFtl 	 7-glutamyl-THF transport, lysosomal
    ABTti 	 L-arabinitol transport via passive diffusion
    ACACt2m 	 Acetoacetate mitochondrial transport via H+ symport
    ACALDt 	 acetaldehyde reversible transport
    ACALDtm 	 acetaldehyde mitochondrial diffusion
    ACALDtx 	 acetaldehyde peroxisomal diffusion
    ACCOAgt 	 acetyl-coa transport
    ACCOAtr 	 acetyl-coa transport
    ACGAMtly 	 N-acetyl-glucosamine lysosomal efflux
    ACNAMlt 	 N-acetylneuraminate transport into lysososme
    ACNAMtn 	 N-acetylneuraminate nuclear import
    ACt2m 	 acetate mitochondrial transport via proton symport
    ACt2r 	 acetate reversible transport via proton symport
    ADNtm 	 adenosine facilated transport in mitochondria
    ADRNt 	 fatty acid transport via diffusion
    AHCYStr 	 S-Adenosyl-L-homocysteine intracellular diffusion
    ALAt2r 	 L-alanine reversible transport via proton symport
    ALAt4 	 Alanine-Sodium symporter
    AMETr 	 S-Adenosyl-L-methionine intracellular diffusion
    ARAB_Lt 	 L-arabinoase extracellular transport
    ARACHCRNt 	 carnitine/acylcarnitine translocase
    ARACHDCOAtx 	 fatty acid intracellular transport
    ARACHDtr 	 intracellular transport
    ARGtiDF 	 L-arginine transport via diffusion (extracellular to cytosol)
    ARTFR12 	 R group artificial flux (C16:1)
    ASCBt 	 L-ascorbate transport via facilitated diffusion
    ASNt4 	 L-asparagine transport in via sodium symport
    Asn_X_Ser_Thrtr 	 Asn-X-Ser/Thr transport (from ER to lysosome)
    ASPGLUm 	 aspartate-glutamate mitochondrial shuttle
    ASPt6 	 L-aspartate transport via Na, H symport and K antiport
    ATP1ter 	 ADP/ATP transporter, endoplasmic reticulum
    ATP2ter 	 AMP/ATP transporter, endoplasmic reticulum
    ATPasel 	 V-type ATPase, H+ transporting, lysosomal
    ATPtm 	 ADP/ATP transporter, mitochondrial
    BALAtmr 	 Beta-alanine reversible mitochondrial transport (diffusion)
    BHBt 	 (R)-3-Hydroxybutanoate transport via H+ symport
    BILGLCURtr 	 glucuronidated compound transport
    BILIRUBtr 	 lipid, flip-flop intracellular transport
    C226CRNt 	 C226 transport into the mitochondria
    CDPDAGtm 	 intracellular transport
    CERT1rt 	 ceramide transport protein
    CERT2rt 	 ceramide transport protein
    CHLtm 	 choline transport via diffusion (cytosol to mitochondria)
    CHOLATEt 	 cholate transport via bicarbonate countertransport
    CHOLATEt3 	 ABC bile acid transporter
    CHSTEROLt2 	 cholesterol intracellular transport
    CLOXAtex2 	 chloride transport via oxalate countertransport (2:1)
    CLPNDt 	 fatty acid transport via diffusion
    CMPACNAtg 	 CMP-Sia Golgi transport via CMP antiport
    CMPACNAtn 	 CMP-Sia nuclear export
    CO2t 	 CO2 transporter via diffusion
    CO2ter 	 CO2 endoplasmic reticular transport via diffusion
    CO2tm 	 CO2 transport (diffusion), mitochondrial
    CO2tp 	 CO2 peroxisomal transport
    COAtm 	 CoA transporter
    COAtp 	 coenzyme A transport, peroxisomal
    COAtr 	 COA transporter, endoplasmic reticulum
    COt 	 CO transporter via diffusion
    CREATt4_2_r 	 Creatine transport (sodium symport) (2:1)
    CREATtmdiffir 	 Creatine transport to/from mitochondria via diffusion
    CRVNCtr 	 fatty acid transport via diffusion
    CTPtn 	 CTP diffusion in nucleus
    CYOOm2 	 cytochrome c oxidase, mitochondrial Complex IV
    CYOR_u10m 	 ubiquinol-6 cytochrome c reductase, Complex III
    CYTDtm 	 cytidine facilated transport in mitochondria
    DAGt 	 diacylglycerol transport
    DATPtn 	 dATP diffusion in nucleus
    DCSPTN1CRNt 	 transport into the mitochondria (carnitine)
    DCSPTN1t 	 fatty acid transport via diffusion
    DCTPtn 	 dCTP diffusion in nucleus
    DGTPtn 	 dGTP diffusion in nucleus
    DHAAt1r 	 dehydroascorbate transport (uniport)
    DHFtl 	 dihydrofolate reversible lysosomal transport
    DHFtm 	 dihydrofolate reversible mitochondrial transport
    D_LACt2 	 D-lactate transport via proton symport
    DLNLCGCRNt 	 transport into the mitochondria (carnitine)
    DNDPt1m 	 dATP transport via ADP antiport
    DNDPt27m 	 dCDP transport via dTDP antiport
    DNDPt29m 	 dCDP transport via dADP antiport
    DNDPt2m 	 dATP transport via ATP antiport
    DNDPt34m 	 dGDP transport via dTDP antiport
    DNDPt35m 	 dGDP transport via dADP antiport
    DNDPt37m 	 dUTP transport via dTDP antiport
    DNDPt39m 	 dUTP transport via dGDP antiport
    DNDPt40m 	 dUTP transport via dADP antiport
    DNDPt41m 	 dUTP transport via dCDP antiport
    DNDPt42m 	 dUTP transport via ADP antiport
    DNDPt43m 	 dUTP transport via ATP antiport
    DNDPt55m 	 dCTP transport via ADP antiport
    DNDPt56m 	 dCTP transport via ATP antiport
    DNDPt57m 	 dGTP transport via ATP antiport
    DNDPt58m 	 dGTP transport via ADP antiport
    DOLGLCP_Lter 	 Dolichyl beta-D-glucosyl phosphate flippase (liver)
    DOLMANP_Lter 	 dolichol-phosphate mannose flippase (liver)
    DOLP_Uter 	 dolichol phosphate flippase (uterus)
    DTDPtn 	 dTDP nuclear transport
    DTTPtn 	 dTTP diffusion in nucleus
    DUDPtn 	 dUDP nuclear transport
    DUMPtn 	 dUMP nuclear transport
    EICOSTETt 	 fatty acid transport via diffusion
    ETOHt 	 ethanol reversible transport
    ETOHtx 	 ethanol reversible peroxisomal transport
    F1Atg 	 F1alpha transport (from golgi to lysosome)
    FADH2tru 	 FADH2 transporter, endoplasmic reticulum
    FADtru 	 FAD transporter, endoplasmic reticulum
    FE2t 	 iron (II) transport
    FE2tm 	 iron (II) transport
    FORtr 	 FOR transporter, endoplasmic reticulum
    FRDPtr 	 lipid, flip-flop intracellular transport
    FUC14GALACGLCGALGLUSIDEte 	 blood group intracellular transport
    FUC14GALACGLCGALGLUSIDEtg 	 blood group intracellular transport
    FUCtly 	 L-fucose efflux from lysosome
    FUMSO3tm 	 Fumarate:sulfite antiport, mitochondrial
    FUMtm 	 fumarate transport, mitochondrial
    GALGLUSIDEtg 	 galgluside hs intracellular transport
    GALSIDEtl 	 galactocerebroside intracellular transport
    GALt1r 	 galactose transport (uniport)
    GALtly 	 galactose efflux from lysosome
    GAMt1r 	 glucosamine transport (uniport)
    GDPFUCtg 	 GDPFuc Golgi transport via CMP antiport
    GDPtg 	 GDP intracellular transport
    GGLUCT 	 gamma-glutamylcyclotransferase
    GLCMter 	 glucose transport via membrane vesicle
    GLCt4 	 glucose transport via sodium symport
    GLCter 	 glucose transport, endoplasmic reticulum
    GLCtg 	 glucose transport, Golgi apparatus
    GLCURter 	 glucuronate endoplasmic reticular transport
    GLCURtly 	 glucuronate transport into lysososme
    GLNt4 	 L-glutamine reversible transport via sodium symport
    GLNtm 	 L-glutamine transport via electroneutral transporter
    GLUt2m 	 L-glutamate reversible transport via proton symport, mitochondrial
    GLUt6 	 Glutamate transport via Na, H symport and K antiport
    GLUt7l 	 Glutamate transport, lysosomal
    GLUtr 	 intracellular transport
    GLXtm 	 glyoxylate transport, mitochondrial
    GLXtp 	 glyoxylate transport, peroxisomal
    GLYBt4_2_r 	 Betaine transport (sodium symport) (2:1)
    GLYBtm 	 Glycine betaine transport via diffusion (mitochondria to cytosol)
    GLYCLTtp 	 glycolate transport into peroxisome
    GLYC_St 	 L-glycerate export
    GLYCtm 	 glycerol transport
    GLYt2r 	 glycine reversible transport via proton symport
    GLYt4 	 glycine transport via sodium symport
    GLYtm 	 glycine passive transport to mitochondria
    GLYtp 	 glycine passive transport to peroxisome
    GMPtg 	 GMP transport (Golgi)
    GTHRDtr 	 glutathione transport via diffusion
    GULNter 	 L-gulonate endoplasmic reticular export
    H2O2t 	 hydrogen peroxide transport via diffusion
    H2O2tm 	 hydrogen peroxide mitochondrial transport
    H2O2tn 	 hydrogen peroxide nuclear transport
    H2O2tp 	 hydrogen peroxide peroxisomal transport via diffusion
    H2Ot 	 H2O transport via diffusion
    H2Oter 	 H2O endoplasmic reticulum transport
    H2Otg 	 H2O transport, Golgi apparatus
    H2Otly 	 H2O transport, lysosomal
    H2Otn 	 H2O transport, nuclear
    H2Otp 	 H2O transport, peroxisomal
    HAS1 	 hyaluronan synthase
    HAS2 	 hyaluronan synthase
    HAtly 	 hyaluronan transport, extracellular to lysosome
    HDD2COAtx 	 hdd2coa intracellular transport
    HISt4 	 L-histidine transport in via sodium symport
    HMGCOAtm 	 Hydroxymethylglutaryl-CoA reversible mitochondrial transport
    HMGCOAtx 	 Hydroxymethylglutaryl-CoA reversible peroxisomal transport
    HOMt4 	 L-homoserine via sodium symport
    HPYRtp 	 hydroxypyruate transport, peroxisomal
    Htg 	 proton diffusion (Golgi)
    Htm 	 Uncoupling protein
    Htr 	 H transporter, endoplasmic reticulum
    Htx 	 H transporter, peroxisome
    HXANtx 	 hypoxanthine diffusion in peroxisome
    INSTt4 	 inositol transport via sodium symport
    IPDPtx 	 Isopentenyl diphosphate transport (peroxisome)
    KSII_CORE2t 	 keratan sulfate II (core2) transport, golgi to extracellular
    KSII_CORE2tly 	 keratan sulfate II (core 2) transport, extracellular to lysosome
    KSII_CORE4t 	 keratan sulfate II (core4) transport, golgi to extracellular
    KSII_CORE4tly 	 keratan sulfate II (core 4) transport, extracellular to lysosome
    KSIt 	 keratan sulfate I transport, golgi to extracellular
    KSItly 	 keratan sulfate I transport, extracellular to lysosome
    LEUKTRA4t 	 leukotriene A4 transport
    LEUKTRA4tr 	 leukotriene intracellular transport
    LEUKTRB4t 	 leukotriene B4 transport
    LEUKTRD4tr 	 leukotriene intracellular transport
    LEUKTRE4t 	 leukotriene E4 transport
    L_LACtm 	 L-lactate transport, mitochondrial
    LNLCCRNt 	 transport into the mitochondria (carnitine)
    LNLNCGCRNt 	 transport into the mitochondria (carnitine)
    LPCHOLt 	 lysophosphatidylcholine transport
    LYStiDF 	 L-lysine transport via diffusion (extracellular to cytosol)
    M4MPDOL_Uter 	 m4mpdol flippase
    M8MASNterg 	 m8masn transport from ER to Golgi apparatus
    MAGt 	 monoacylglycerol 2 transport
    MALSO3tm 	 Malate:sulfite antiport, mitochondrial
    MALtm 	 malate transport, mitochondrial
    MANt1r 	 mannose transport (uniport)
    MANtg 	 mannose efflux from Golgi apparatus
    MANtly 	 mannose efflux from lysosome
    MEOHt2 	 Methanol diffusion
    N2M2NMASNt 	 n2m2nmasn transport, Golgi to extracellular
    N2M2NMASNtly 	 n2m2nmasn transport, extracellular to lysosome
    NADHtpu 	 NADH transporter, peroxisome
    NADHtru 	 NADH transporter, endoplasmic reticulum
    NADPHtru 	 NADPH transporter, endoplasmic reticulum
    NADPHtxu 	 NADPH transporter, peroxisome
    NADPtru 	 NADP transporter, endoplasmic reticulum
    NADPtxu 	 NADP transporter, peroxisome
    NADtru 	 NAD transporter, endoplasmic reticulum
    NaKt 	 Na+/K+ exchanging ATPase
    NAt 	 sodium transport (uniport)
    NAt3_1 	 sodium proton antiporter (H:NA is 1:1)
    NH4tp 	 ammonia peroxisomal transport
    NKCC2t 	 Na+-K+-Cl- cotransport (NH4+)
    NKCCt 	 Na+-K+-Cl- cotransport
    NRVNCt 	 fatty acid transport via diffusion
    O2St 	 superoxide anion transport via diffusion (extracellular)
    O2Stm 	 superoxide anion transport via diffusion (mitochondria)
    O2Stn 	 superoxide anion transport via diffusion (nucleus)
    O2Stx 	 superoxide anion transport via diffusion (peroxisome)
    O2t 	 o2 transport (diffusion)
    O2ter 	 O2  transport, endoplasmic reticulum
    O2tm 	 O2 transport (diffusion)
    O2tn 	 O2 nuclear transport
    O2tp 	 O2  transport, peroxisomal
    OCCOAtm 	 octanoyl(n-C8:0)-CoA transport
    OCCOAtx 	 fatty acid intracellular transport
    OCTAt 	 Octanoate transport via diffusion
    ORNtiDF 	 ornithine transport via diffusion (extracellular to cytosol)
    OXAHCOtex 	 oxalate transport via bicarbonate countertransport
    PAIL_HStn 	 phosphatidylinositol nuclear transport (diffusion)
    PAIL45P_HStn 	 phosphatidylinositol 4,5-bisphosphate nuclear transport (diffusion)
    PAIL4P_HStn 	 phosphatidylinositol 4-phosphate nuclear transport (diffusion)
    PAPStg 	 3-Phosphoadenylyl sufate Golgi transport
    PAPtg 	 adenosine 3,5-bisphosphate Golgi transport
    PCREATtmdiffir 	 Phosphocreatine transport to/from mitochondria via diffusion
    PE_HSter 	 phosphatidylethanolamine scramblase
    PE_HStm 	 phosphatidylethanolamine scramblase
    PEt 	 phosphatidylethanolamine transport
    PGLYCt 	 phosphatidylglycerol transport
    PHEMEtm 	 Heme transport to cytosol
    PHEt4 	 L-phenylalanine transport in via sodium symport
    PIt2m 	 phosphate transporter, mitochondrial
    PItg 	 phosphate transport, Golgi apparatus
    PItn 	 phosphate transport, nuclear
    PItx 	 Phosphate transporter, peroxisome
    PPAt 	 Propionate transport, diffusion
    PPItr 	 Diphosphate transporter, endoplasmic reticulum
    PPItx 	 Diphosphate transporter, peroxisome
    PPPG9tm 	 protoporphyrinogen IX mitochondrial transport
    PRISTANALtx 	 pristanal peroxisomal transport
    PRISTtx 	 prist peroxisomal transport
    PROSTGH2t 	 Prostaglandin H2 transport
    PROSTGI2t 	 Prostaglandin I2 transport
    PROt2r 	 L-proline reversible transport via proton symport
    PROtm 	 L-proline transport, mitochondrial
    PSt3 	 phosphatidylserine transport
    PYRt2m 	 pyruvate mitochondrial transport via proton symport
    PYRt2p 	 pyruvate peroxisomal transport via proton symport
    RETFAt 	 fatty acid retinol efflux
    RETNt 	 retinoic acid transport
    RETt 	 retinol trasnport by STRA6
    RTOTALt 	 RTOTAL transport
    Rtotaltl 	 fatty acid intracellular transport
    S2L2FN2M2MASNt 	 s2l2fn2m2masn transport, Golgi to extracellular
    S2L2FN2M2MASNtly 	 s2l2fn2m2masn transport, extracellular to lysosome
    S2L2N2M2MASNtly 	 s2l2n2m2masn transport, extracellular to lysosome
    SCP21x 	 Sterol carrier protein 2
    SCP22x 	 Sterol carrier protein 2
    Ser_Thrtg 	 Ser/Thr transport (from golgi to lysosome)
    SERt4 	 L-serine via sodium symport
    SERtp 	 L-serine transport, peroxisomal
    SO4OXAtex2 	 sulfate transport via oxalate countertransport (2:1)
    SO4tl 	 Sulfate transport (lysosome)
    SPC_HSt 	 sphingosylphosphorylcholine transport (diffusion)
    SPH1Pte 	 sph1p transport
    SPHINGStl 	 sphingosine intracellular transport
    SPHS1Pte 	 sphingosine-1-phosphate transport
    STRDNCt 	 fatty acid transport via diffusion
    SUCCt2m 	 succinate transport, mitochondrial
    TCHOLAt 	 taurocholate transport via bicarbonate countertransport
    TCHOLAte 	 bile acid intracellular transport
    THFtl 	 5,6,7,8-Tetrahydrofolate transport, diffusion, lysosomal
    THFtm 	 5,6,7,8-Tetrahydrofolate transport, diffusion, mitochondrial
    THRt4 	 L-threonine  via sodium symport
    TMNDNCt 	 fatty acid transport via diffusion
    TRPt4 	 L-tryptophan transport in via sodium symport
    UDPACGALtl 	 udpacgal intracellular transport
    UDPGALtg 	 UDP-Gal Golgi transport via CMP antiport
    UDPGLCAter 	 UDPGlcA endoplasmic reticulum transport via UMP antiport
    UDPGLCter 	 UDP-Glc endoplasmic reticulum transport via CMP antiport
    UDPtl 	 udp intracellular transport
    UGALNACtg 	 UDP-GalNAc Golgi transport via CMP antiport
    UGLCNACtg 	 UDP-GlcNAc Golgi transport via CMP antiport
    URATEt 	 urate export from cytosol
    URATEtx 	 urate export from peroxisome
    VACCt 	 fatty acid transport via diffusion
    WHHDCAte 	 xenobiotic transport
    XANtx 	 xanthine diffusion in peroxisome
    XOL27OHtm 	 27 hydroxy cholesterol transport
    XOLDIOLONEt 	 lipid, flip-flop intracellular transport
    XYLTt 	 Xylitol transport via passive diffusion
    r0002 	 Active transport
    r0801 	 Mitochondrial Carrier (MC) TCDB:2.A.29.21.1
    r0809 	 Facilitated diffusion
    r0812 	 Mitochondrial Carrier (MC) TCDB:2.A.29.20.1
    r0818 	 Transport reaction
    r0819 	 Mitochondrial Carrier (MC) TCDB:2.A.29.2.2
    r0821 	 Mitochondrial Carrier (MC) TCDB:2.A.29.2.2
    r0822 	 Mitochondrial Carrier (MC) TCDB:2.A.29.2.2
    r0830 	 Mitochondrial Carrier (MC) TCDB:2.A.29.2.2
    r0838 	 Free diffusion
    r0853 	 Facilitated diffusion
    r0860 	 Utilized transport
    r0870 	 Facilitated diffusion
    r0892 	 Facilitated diffusion
    r0899 	 Amino acid transporter ATB0+ Facilitated diffusion
    r0911 	 Facilitated diffusion
    r0913 	 Mitochondrial Carrier (MC) TCDB:2.A.29.7.2
    r0917 	 Mitochondrial Carrier (MC) TCDB:2.A.29.7.2
    r0926 	 Postulated transport reaction
    r0927 	 Free diffusion
    r0940 	 Free diffusion
    r0941 	 Free diffusion
    r0942 	 Neurotransmitter:Sodium Symporter (NSS) TCDB:2.A.22.3.4
    r0950 	 Facilitated diffusion
    r0954 	 MCT 1 Transport reaction
    r0973 	 Facilitated diffusion
    r0990 	 Postulated transport reaction
    r1000 	 Facilitated diffusion
    r1004 	 Facilitated diffusion
    r1007 	 Facilitated diffusion
    r1010 	 Postulated transport reaction
    r1013 	 Transport reaction
    r1014 	 Postulated transport reaction
    r1030 	 Transport reaction
    r1088 	 Active transport
    r1106 	 Postulated transport reaction
    r1116 	 ATP Exporter (ATP-E) TCDB:9.A.6.1.1
    r1117 	 Facilitated diffusion
    r1143 	 Major Facilitator(MFS) TCDB:2.A.18.6.3
    r1144 	 Major Facilitator(MFS) TCDB:2.A.18.6.3
    r1147 	 Mitochondrial Carrier (MC) TCDB:2.A.29.7.2
    r1148 	 Active transport
    r1162 	 Active transport
    r1291 	 Postulated transport reaction
    r1292 	 Postulated transport reaction
    r1375 	 Transport reaction
    r1400 	 Active transport
    r1421 	 Free diffusion
    r1434 	 Transport reaction
    r1455 	 Transport reaction
    r1456 	 Transport reaction
    r1459 	 Transport reaction
    r1464 	 Active transport
    r1493 	 Utilized transport
    r1515 	 ATP-binding Cassette (ABC) TCDB:3.A.1.211.1
    r1516 	 ATP-binding Cassette (ABC) TCDB:3.A.1.211.1
    r1522 	 ATP-binding Cassette (ABC) TCDB:3.A.1.211.1
    r1523 	 ATP-binding Cassette (ABC) TCDB:3.A.1.211.1
    r1531 	 ATP-binding Cassette (ABC) TCDB:3.A.1.208.15
    r1532 	 ATP-binding Cassette (ABC) TCDB:3.A.1.208.15
    r1546 	 Amino Acid-Polyamine-Organocation (APC) TCDB:2.A.3.8.1
    r1547 	 Amino Acid-Polyamine-Organocation (APC) TCDB:2.A.3.8.1
    r1556 	 Amino Acid-Polyamine-Organocation (APC) TCDB:2.A.3.8.1
    r1560 	 Amino Acid-Polyamine-Organocation (APC) TCDB:2.A.3.8.1
    r1561 	 Amino Acid-Polyamine-Organocation (APC) TCDB:2.A.3.8.1
    r1570 	 Amino Acid-Polyamine-Organocation (APC) TCDB:2.A.3.8.1
    r1573 	 Amino Acid-Polyamine-Organocation (APC) TCDB:2.A.3.8.1
    r1574 	 Amino Acid-Polyamine-Organocation (APC) TCDB:2.A.3.8.1
    r1583 	 Amino Acid-Polyamine-Organocation (APC) TCDB:2.A.3.8.1
    r1585 	 Amino Acid-Polyamine-Organocation (APC) TCDB:2.A.3.8.1
    r1586 	 Amino Acid-Polyamine-Organocation (APC) TCDB:2.A.3.8.1
    r1595 	 Amino Acid-Polyamine-Organocation (APC) TCDB:2.A.3.8.1
    r1598 	 Amino Acid-Polyamine-Organocation (APC) TCDB:2.A.3.8.1
    r1602 	 Amino Acid-Polyamine-Organocation (APC) TCDB:2.A.3.8.1
    r1603 	 Amino Acid-Polyamine-Organocation (APC) TCDB:2.A.3.8.1
    r1606 	 Amino Acid-Polyamine-Organocation (APC) TCDB:2.A.3.8.1
    r1633 	 Amino Acid-Polyamine-Organocation (APC) TCDB:2.A.3.8.1
    r1640 	 Amino Acid-Polyamine-Organocation (APC) TCDB:2.A.3.8.1
    r1646 	 Amino Acid-Polyamine-Organocation (APC) TCDB:2.A.3.8.1
    r1651 	 Amino Acid-Polyamine-Organocation (APC) TCDB:2.A.3.8.1
    r1655 	 Amino Acid-Polyamine-Organocation (APC) TCDB:2.A.3.8.1
    r1658 	 Amino Acid-Polyamine-Organocation (APC) TCDB:2.A.3.8.1
    r1660 	 Amino Acid-Polyamine-Organocation (APC) TCDB:2.A.3.8.1
    r1662 	 Amino Acid-Polyamine-Organocation (APC) TCDB:2.A.3.8.1
    r2084 	 Major Facilitator(MFS) TCDB:2.A.1.13.1
    r2093 	 Major Facilitator(MFS) TCDB:2.A.1.13.1
    r2099 	 Major Facilitator(MFS) TCDB:2.A.1.13.1
    r2101 	 Major Facilitator(MFS) TCDB:2.A.1.13.1
    r2102 	 Major Facilitator(MFS) TCDB:2.A.1.13.1
    r2103 	 Major Facilitator(MFS) TCDB:2.A.1.13.1
    r2104 	 Major Facilitator(MFS) TCDB:2.A.1.13.1
    r2105 	 Major Facilitator(MFS) TCDB:2.A.1.13.1
    r2107 	 Major Facilitator(MFS) TCDB:2.A.1.13.1
    r2108 	 Major Facilitator(MFS) TCDB:2.A.1.13.1
    r2114 	 Major Facilitator(MFS) TCDB:2.A.1.13.1
    r2119 	 Major Facilitator(MFS) TCDB:2.A.1.13.1
    r2123 	 Major Facilitator(MFS) TCDB:2.A.1.13.1
    r2128 	 Major Facilitator(MFS) TCDB:2.A.1.13.1
    r2132 	 Major Facilitator(MFS) TCDB:2.A.1.13.1
    r2139 	 Resistance-Nodulation-Cell Division (RND) TCDB:2.A.60.1.14
    r2151 	 Resistance-Nodulation-Cell Division (RND) TCDB:2.A.60.1.14
    r2152 	 Resistance-Nodulation-Cell Division (RND) TCDB:2.A.60.1.14
    r2353 	 Organic anion transporter 5 Utilized transport
    r2368 	 Organic anion transporter 5 Utilized transport
    r2369 	 Organic anion transporter 5 Utilized transport
    r2372 	 Mitochondrial Carrier (MC) TCDB:2.A.29.7.2
    r2375 	 Mitochondrial Carrier (MC) TCDB:2.A.29.7.2
    r2419 	 Mitochondrial Carrier (MC) TCDB:2.A.29.2.7
    r2444 	 Proposed Fatty Acid Transporter (FAT) TCDB:4.C.1.1.5
    r2472 	 Major Facilitator(MFS) TCDB:2.A.1.4.7
    r2509 	 Utilized transport
    r2517 	 ATP-binding Cassette (ABC) TCDB:3.A.1.203.1
    r2521 	 Utilized transport
    r2537 	 Utilized transport
    r2538 	 Utilized transport
    PIt8 	 phosphate transport in/out via Na+ symporter
    2MB2COAc 	 transport of 2-methylcrotonoyl-CoA into cytosol
    3HBCOARc 	 transport of (R)-3-hydroxybutanoyl-CoA into cytosol
    C141ACBP 	 transport of 3-hydroxytetradecenoyl coa from mitochondria into cytosol
    C142ACBP 	 transport of 3-hydroxytetradeca dienoyl coa from mitochondria into cytosol
    C142OHe 	 excretion of C14:2-OH
    C162ACBP 	 transport of 3-hydroxy trans7,10-hexadecadienoylcoa from mitochondria into cytosol
    C6CRNe 	 transport of hexanoyl carnitine into the extra cellular fluid
    C8CRNe 	 transport of octanoyl carnitine into the extra cellular space
    DDCRNe 	 transport of 3-hydroxydodecanoyl carnitine into extra cellular space
    DDECCRNe 	 transport of lauroyl carnitine into extra cellular space
    DECCRNe 	 transport of 3-hydroxydecanoyl carnitine into extra cellular space
    GLUTCOAACBP 	 transport of glutaryl-CoA(5-) from mitochondria into cytosol
    HDCACBP 	 transport of (S)-3-Hydroxydecanoyl-CoA from mitochondria into the cytosol
    HDDACBP 	 transport of (S)-3-Hydroxydodecanoyl-CoA from mitochondria into the cytosol
    HDECAACBP 	 transport of 3-hydroxyhexadecanoylcoa from mitochondria into the cytosol
    HDECEACBP 	 transport of 3-hydroxyhexadecenoylcoa from mitochondria into the cytosol
    HEXCOAACBP 	 hexanoyl(n-C6:0)-CoA transport
    HEXDIACTD 	 transport of hexadecanedioc acid by diffusion
    HEXDICOAACBP 	 transport of hexadecanedioyl coa
    HEXDICOAACBPx 	 transport of hexadecanedioyl coa
    HIVCACBP 	 transport of 3-hydroxyisovalerylcoa from mitochondria into the cytosol
    HIVCRNe 	 transport of 3-hydroxy-isovaleryl carnitine into extra cellular space
    HOCDACBP 	 transport of (S)-3-Hydroxyoctadecanoyl-CoA from mitochondria into the cytosol
    HOCTDACBP 	 transport of 3-hydroxyoctadecenoylcoa from mitochondria into the cytosol
    HOCTDECCRNe 	 transport of 3-hydroxyoctadecanoyl carnitine into extra cellular space
    HTDCACBP 	 transport of (S)-3-Hydroxytetradecanoyl-CoA from mitochondria into the cytosol
    HTDCRNe 	 transport of 3-hydroxy-tetradecanoyl carnitine into extra cellular space
    IVCOAACBP 	 transport of Isovaleryl-CoA from mitochondria into cytosol
    IVCRNe 	 transport of isovaleryl carnitine into the extra cellular space
    OCD11CRNCACT 	 transport of 11-octadecenoyl carnitine into the mitochondrial matrix
    OCTDEC2ACBP 	 transport of 3-hydroxyoctadecadienoyl coa from mitochondria into cytosol
    SUCCTD 	 transport of succinate by diffusion
    TIGCRNe 	 transport of tiglyl carnitine into the extra cellular fluid
    4OHPROIMINOtc 	 transport of L-OH-Proline by the apical IMINO amino acid transporters in kidney and intestine
    CARPEPT1tc 	 transport of L-Carnosine by the apical PEPT1 amino acid transporters across the brush border cells of the enterocytes of the intestine and renal cells
    CYSPHELAT2tc 	 transport of L-Cysteine into the cell and efflux of L-Phenylalanine out of the cell by LAT2 on the basolateral surfaces of kidney and intestine.
    ILELAT1tc 	 transport of L-Isoleucine by LAT1 in association with 4F2hc, across the apical surface of the memebranes.
    ILEPHELAT2tc 	 transport of L-Isoleucine into the cell and efflux of L-Phenylalanine out of the cell by LAT2 on the basolateral surfaces of kidney and intestine.
    LEUPHELAT2tc 	 Transport of L-Leucine into the cell and efflux of L-Phenylalanine out of the cell by LAT2 on the basolateral surfaces of kidney and intestine.
    LINOFATPtc 	 uptake of linoleic acid by the enterocytes
    NACDe 	 release of nicotinate at the basolateral surface
    SBTle 	 diffusion of sorbitol into the enterocytes
    THRPHELAT2tc 	 transport of L-Threonine into the cell and efflux of L-Phenylalanine out of the cell by LAT2 on the basolateral surfaces of kidney and intestine.
    TRPATB0tc 	 transport of L-Tryptophan into the intestinal cells by ATB0 transporter
    TYRPHELAT2tc 	 transport of L-Tyrosine into the cell and efflux of L-Phenylalanine out of the cell by LAT2 on the basolateral surfaces of kidney and intestine.
    VALLAT1tc 	 transport of L-Valine by LAT1 in association with 4F2hc, across the apical surface of the memebranes.
    VALPHELAT2tc 	 transport of L-Valine into the cell and efflux of L-Phenylalanine out of the cell by LAT2 on the basolateral surfaces of kidney and intestine.
    AHCYStd 	 diffusion of S-Adenosyl-L-homocysteine
    AMETtd 	 diffusion of S-Adenosyl-L-methionine
    ASPDTDe 	 D-aspartate transport, extracellular
    CHSTEROLtrc 	 transport of cholesterol into the cytosol
    Coqe 	 Transport of q10
    DATPtm 	 transport of dATP into mitochondria
    DCTPtm 	 transport of dCTP into mitochondria
    DHAPtc 	 transport of DHAP into cytosol
    DTTPtm 	 transport of dTTP into mitochondria
    FE3MTP1 	 transport of ferrous iron into blood
    GLYC3Ptmc 	 glycerol-3-phopshate transport, cytoplasm
    GLYCTDle 	 difussion of glycerol accross the brush border membrane
    KHte 	 KHte
    LACLt 	 Lactate transport
    NADtm 	 transport of NAD into mitochondria
    NADtx 	 transport of NAD into peroxisome
    PHEMEe 	 release of heme into the blood
    PPItm 	 diphopshate transporter, mitochondrial
    PRO_Dtde 	 D-proline transport, extracellular
    PTCRTD 	 diffusion of putriscine into the endothelial cells
    Q10H2e 	 transport of ubiquinol into lymph
    RTOTALFATPc 	 uptake of Rtotal by enterocytes
    SFCYSe 	 Exit of sulfocysteine into extra-cellular space
    SPMTDe 	 SPMTDe
    q10h2tc 	 transport of ubiquinol into cytosol
    q10tm 	 transport of ubiquinone into mitochondria
    3MOBte 	 Transport of 3-methyl-2-oxobutanoate
    3MOPte 	 Transport of 3-methyl-2-oxopentanoate
    4MOPte 	 Transport of 4-methyl-2-oxopentanoate
    5OXPROt 	 Transport of 5-oxoprolinate
    AHCYSte 	 Transport of S-adenosyl-L-homocysteine
    MAL_Lte 	 Transport of L-malate
    fumt 	 Fumarate transport
    xmpt 	 Xanthosine 5-phosphate transport
    dcmpt 	 dCMP transport
    HC00342te 	 cis-aconitate transport
    glyc3pte 	 glycerol 3-phosphate transport
    DM_Asn_X_Ser_Thr_ly_ 	 DM Asn-X-Ser/Thr(ly)
    DM_core5_g_ 	 DM core5(g)
    DM_core7_g_ 	 DM core7(g)
    DM_datp_m_ 	 dATP demand
    DM_datp_n_ 	 DM datp(n)
    DM_dctp_m_ 	 dCTP demand
    DM_dctp_n_ 	 DM dctp(n)
    DM_dgtp_m_ 	 DM dgtp(m)
    DM_dgtp_n_ 	 DM dgtp(n)
    DM_dttp_m_ 	 dTTP demand
    DM_dttp_n_ 	 DM dttp(n)
    DM_Ser_Thr_ly_ 	 DM Ser/Thr(ly)
    DM_sprm_c_ 	 DM sprm(c)
    DM_T_antigen_g_ 	 DM T antigen(g)
    sink_citr_LPAREN_c_RPAREN_ 	 citrulline sink
    DM_taur_LPAREN_c_RPAREN_ 	 demand reaction for taurine
    DM_pe_hs_LPAREN_r_RPAREN_ 	 demand reaction for pe_hs[r]
    DM_4hrpo 	 Demand for trans-4-hydroxy-L-proline
    DM_Lcystin 	 Demand for L-cystine
    DM_anth 	 Demand for anthranilate
    DM_ncam 	 Demand for nicotinamide
    EX_1mncam_LPAREN_e_RPAREN_ 	 1-Methylnicotinamide exchange
    EX_abt_LPAREN_e_RPAREN_ 	 L-Arabinitol exchange
    EX_ac_LPAREN_e_RPAREN_ 	 Acetate exchange
    EX_acac_LPAREN_e_RPAREN_ 	 Acetoacetate exchange
    EX_acald_LPAREN_e_RPAREN_ 	 Acetaldehyde exchange
    EX_acetone_LPAREN_e_RPAREN_ 	 Acetone exchange
    EX_adn_LPAREN_e_RPAREN_ 	 exchange reaction for Adenosine
    EX_adrn_LPAREN_e_RPAREN_ 	 adrenic acid exchange
    EX_akg_LPAREN_e_RPAREN_ 	 2-Oxoglutarate exchange
    EX_ala_B_LPAREN_e_RPAREN_ 	 beta-Alanine exchange
    EX_ala_L_LPAREN_e_RPAREN_ 	 exchange reaction for L-alanine
    EX_arab_L_LPAREN_e_RPAREN_ 	 L-Arabinose exchange
    EX_arg_L_LPAREN_e_RPAREN_ 	 L-Arginine exchange
    EX_ascb_L_LPAREN_e_RPAREN_ 	 L-Ascorbate exchange
    EX_asn_L_LPAREN_e_RPAREN_ 	 exchange reaction for L-asparagine
    EX_asp_D_LPAREN_e_RPAREN_ 	 D-Aspartate exchange
    EX_asp_L_LPAREN_e_RPAREN_ 	 L-Aspartate exchange
    EX_atp_LPAREN_e_RPAREN_ 	 ATP exchange
    EX_bhb_LPAREN_e_RPAREN_ 	 (R)-3-Hydroxybutanoate transport via H+ symport
    EX_bilglcur_LPAREN_e_RPAREN_ 	 exchange reaction for blirubin mono-glucuronide
    EX_cholate_LPAREN_e_RPAREN_ 	 exchange reaction for cholate
    EX_cit_LPAREN_e_RPAREN_ 	 Citrate exchange
    EX_CLPND_LPAREN_e_RPAREN_ 	 clupanodonic acid (docosapentaenoic (n-3)) exchange
    EX_cmp_LPAREN_e_RPAREN_ 	 CMP exchange
    EX_co_LPAREN_e_RPAREN_ 	 Carbon monoxide exchange
    EX_co2_LPAREN_e_RPAREN_ 	 CO2 exchange
    EX_creat_LPAREN_e_RPAREN_ 	 Creatine exchange
    EX_crn_LPAREN_e_RPAREN_ 	 exchange reaction for L-Carnitine
    EX_crvnc_LPAREN_e_RPAREN_ 	 nC22:6 exchange
    EX_cys_L_LPAREN_e_RPAREN_ 	 exchange reaction for L-cysteine
    EX_dag_hs_LPAREN_e_RPAREN_ 	 exchange reaction for diglyceride
    EX_dcsptn1_LPAREN_e_RPAREN_ 	 docosa-4,7,10,13,16-pentaenoic acid (n-6) exchange
    EX_dhdascb_LPAREN_e_RPAREN_ 	 exchange reaction for dehydroascorbide(1-)
    EX_eicostet_LPAREN_e_RPAREN_ 	 eicosatetranoic acid exchange
    EX_etoh_LPAREN_e_RPAREN_ 	 Ethanol exchange
    EX_fe2_LPAREN_e_RPAREN_ 	 EX_fe2(u)
    EX_fe3_LPAREN_e_RPAREN_ 	 Fe3+ exchange
    EX_for_LPAREN_e_RPAREN_ 	 Formate exchange
    EX_fuc14galacglcgalgluside_hs_LPAREN_e_RPAREN_ 	 Lea glycolipid exchange
    EX_fuc_L_LPAREN_e_RPAREN_ 	 L-Fucose exchange
    EX_gal_LPAREN_e_RPAREN_ 	 exchange reaction for D-Galactose
    EX_gam_LPAREN_e_RPAREN_ 	 exchange reaction for D-Glucosamine
    EX_glc_LPAREN_e_RPAREN_ 	 D-Glucose exchange
    EX_gln_L_LPAREN_e_RPAREN_ 	 exchange reaction for L-glutamine
    EX_gluala_LPAREN_e_RPAREN_ 	 5-L-Glutamyl-L-alanine exchange
    EX_glu_L_LPAREN_e_RPAREN_ 	 L-Glutamate exchange
    EX_gly_LPAREN_e_RPAREN_ 	 exchange reaction for Glycine
    EX_glyb_LPAREN_e_RPAREN_ 	 Glycine betaine exchange
    EX_glyc_LPAREN_e_RPAREN_ 	 Glycerol exchange
    EX_glyc_S_LPAREN_e_RPAREN_ 	 (S)-Glycerate exchange
    EX_glygn2_LPAREN_e_RPAREN_ 	 glycogen, structure 2 (glycogenin-1,6-{7[1,4-Glc], 4[1,4-Glc]}) exchange
    EX_glygn4_LPAREN_e_RPAREN_ 	 exchange reaction for glycogen, structure 4 (glycogenin-1,6-{2[1,4-Glc], [1,4-Glc]})
    EX_gthox_LPAREN_e_RPAREN_ 	 Oxidized glutathione exchange
    EX_gthrd_LPAREN_e_RPAREN_ 	 Reduced glutathione exchange
    EX_h_LPAREN_e_RPAREN_ 	 exchange reaction for proton
    EX_h2o_LPAREN_e_RPAREN_ 	 H2O exchange
    EX_h2o2_LPAREN_e_RPAREN_ 	 Hydrogen peroxide exchange
    EX_ha_LPAREN_e_RPAREN_ 	 hyaluronan exchange
    EX_ha_pre1_LPAREN_e_RPAREN_ 	 hyaluronan biosynthesis, precursor 1 exchange
    EX_hdca_LPAREN_e_RPAREN_ 	 Hexadecanoate (n-C16:0) exchange
    EX_hdcea_LPAREN_e_RPAREN_ 	 exchange reaction for Hexadecenoate (n-C16:1)
    EX_his_L_LPAREN_e_RPAREN_ 	 exchange reaction for L-histidine
    EX_ile_L_LPAREN_e_RPAREN_ 	 L-Isoleucine exchange
    EX_inost_LPAREN_e_RPAREN_ 	 myo-Inositol exchange
    EX_ins_LPAREN_e_RPAREN_ 	 Inosine exchange
    EX_lac_D_LPAREN_e_RPAREN_ 	 D-lactate exchange
    EX_lac_L_LPAREN_e_RPAREN_ 	 L-Lactate exchange
    EX_leuktrA4_LPAREN_e_RPAREN_ 	 leukotriene A4 exchange
    EX_leuktrB4_LPAREN_e_RPAREN_ 	 leukotriene B4 exchange
    EX_leuktrE4_LPAREN_e_RPAREN_ 	 leukotriene E4 exchange
    EX_leu_L_LPAREN_e_RPAREN_ 	 L-Leucine exchange
    EX_lnlc_LPAREN_e_RPAREN_ 	 linoleic acid (all cis C18:2) exchange
    EX_lpchol_hs_LPAREN_e_RPAREN_ 	 exchange reaction for lysophosphatidylcholine
    EX_lys_L_LPAREN_e_RPAREN_ 	 L-Lysine exchange
    EX_mag_hs_LPAREN_e_RPAREN_ 	 Monoacylglycerol 2 (homo sapiens) exchange
    EX_man_LPAREN_e_RPAREN_ 	 exchange reaction for D-Mannose
    EX_meoh_LPAREN_e_RPAREN_ 	 methanol exchange
    EX_met_L_LPAREN_e_RPAREN_ 	 L-Methionine exchange
    EX_n2m2nmasn_LPAREN_e_RPAREN_ 	 N-Acetyl-beta-D-glucosaminyl-1,2-alpha-D-mannosyl-1,3-(N-acetyl-beta-D-glucosaminyl-1,2-alpha-D-mannosyl-1,6)-(N-acetyl-beta-D-glucosaminyl-1,4)-beta-D-mannosyl-1,4-N-acetyl-beta-D-glucosaminyl-R exchange
    EX_nac_LPAREN_e_RPAREN_ 	 Nicotinate exchange
    EX_nh4_LPAREN_e_RPAREN_ 	 Ammonia exchange
    EX_nrvnc_LPAREN_e_RPAREN_ 	 nervonic acid exchange
    EX_o2_LPAREN_e_RPAREN_ 	 exchange reaction for oxugen
    EX_o2s_LPAREN_e_RPAREN_ 	 Superoxide anion exchange
    EX_ocdca_LPAREN_e_RPAREN_ 	 echange reaction for octadecanoate (n-C18:0)
    EX_ocdcea_LPAREN_e_RPAREN_ 	 octadecenoate (n-C18:1) exchange
    EX_octa_LPAREN_e_RPAREN_ 	 exchange reaction for octanoate (n-C8:0)
    EX_orn_LPAREN_e_RPAREN_ 	 Ornithine exchange
    EX_oxa_LPAREN_e_RPAREN_ 	 Oxalate exchange
    EX_pe_hs_LPAREN_e_RPAREN_ 	 exchange reaction for phosphatidylethanolamine
    EX_pglyc_hs_LPAREN_e_RPAREN_ 	 phosphatidylglycerol (homo sapiens) exchange
    EX_phe_L_LPAREN_e_RPAREN_ 	 exchange reaction for L-phenylalanine
    EX_pheme_LPAREN_e_RPAREN_ 	 exchange reaction for heme
    EX_pi_LPAREN_e_RPAREN_ 	 Phosphate exchange
    EX_ppa_LPAREN_e_RPAREN_ 	 Propionate exchange
    EX_pro_D_LPAREN_e_RPAREN_ 	 D-Proline exchange
    EX_pro_L_LPAREN_e_RPAREN_ 	 L-Proline exchange
    EX_ps_hs_LPAREN_e_RPAREN_ 	 phosphatidylserine (homo sapiens) exchange
    EX_pyr_LPAREN_e_RPAREN_ 	 Pyruvate exchange
    EX_retfa_LPAREN_e_RPAREN_ 	 fatty acid retinol exchange
    EX_retinol_LPAREN_e_RPAREN_ 	 Retinol exchange
    EX_retn_LPAREN_e_RPAREN_ 	 Retinoate exchange
    EX_ribflv_LPAREN_e_RPAREN_ 	 exchange reaction for Riboflavin
    EX_Rtotal_LPAREN_e_RPAREN_ 	 R total exchange
    EX_ser_L_LPAREN_e_RPAREN_ 	 exchange reaction for L-serine
    EX_so4_LPAREN_e_RPAREN_ 	 Sulfate exchange
    EX_spc_hs_LPAREN_e_RPAREN_ 	 sphingosylphosphorylcholine (homo sapiens) exchange
    EX_sph1p_LPAREN_e_RPAREN_ 	 Sphinganine 1-phosphate exchange
    EX_sphs1p_LPAREN_e_RPAREN_ 	 Sphingosine 1-phosphate exchange
    EX_strch1_LPAREN_e_RPAREN_ 	 starch, structure 1 (1,6-{7[1,4-Glc], 4[1,4-Glc]}) exchange
    EX_strch2_LPAREN_e_RPAREN_ 	 exchange reaction for starch, structure 2 (1,6-{2[1,4-Glc], [1,4-Glc]})
    EX_strdnc_LPAREN_e_RPAREN_ 	 stearidonic acid exchange
    EX_thr_L_LPAREN_e_RPAREN_ 	 L-Threonine exchange
    EX_tmndnc_LPAREN_e_RPAREN_ 	 timnodonic acid exchange
    EX_trp_L_LPAREN_e_RPAREN_ 	 L-Tryptophan exchange
    EX_tyr_L_LPAREN_e_RPAREN_ 	 L-Tyrosine exchange
    EX_urate_LPAREN_e_RPAREN_ 	 Urate exchange
    EX_utp_LPAREN_e_RPAREN_ 	 UTP exchange
    EX_vacc_LPAREN_e_RPAREN_ 	 vaccenic acid exchange
    EX_val_L_LPAREN_e_RPAREN_ 	 L-Valine exchange
    EX_whhdca_LPAREN_e_RPAREN_ 	 omega hydroxy hexadecanoate (n-C16:0) exchange
    EX_xylt_LPAREN_e_RPAREN_ 	 exchange reaction for xylitol
    EX_ctp_LPAREN_e_RPAREN_ 	 Exchange of CTP(4-)
    EX_dtmp_LPAREN_e_RPAREN_ 	 Exchange of dTMP(2-)
    EX_dttp_LPAREN_e_RPAREN_ 	 Exchange of dTTP(4-)
    EX_fad_LPAREN_e_RPAREN_ 	 Exchange of Flavin adenine dinucleotide oxidized
    EX_fald_LPAREN_e_RPAREN_ 	 Exchange of formaldehyde
    EX_HC00250_LPAREN_e_RPAREN_ 	 Exchange of hydrosulfide
    EX_HC01609_LPAREN_e_RPAREN_ 	 Exchange of UroporphyrinogenI
    EX_HC01610_LPAREN_e_RPAREN_ 	 Exchange of CoproporphyrinogenI
    EX_ptrc_LPAREN_e_RPAREN_ 	 Exchange of 1,4-butanediammonium
    EX_spmd_LPAREN_e_RPAREN_ 	 Exchange of spermidine(3+)
    EX_prostgh2_LPAREN_e_RPAREN_ 	 prostaglandin H2(1-) exchange
    EX_C02470_LPAREN_e_RPAREN_ 	 Xanthurenic acid exchange
    EX_HC00822_LPAREN_e_RPAREN_ 	 Chitobiose exchange
    EX_ddca_LPAREN_e_RPAREN_ 	 laurate exchange
    EX_3ddcrn_ 	 exchange reaction for 3-hydroxydodecanoyl carnitine
    EX_3deccrn_ 	 exchange reaction for 3-hydroxydecanoyl carnitine
    EX_3ivcrn_ 	 exchange reaction for 3-hydroxy-isovaleryl carnitine
    EX_3octdeccrn_ 	 exchange reaction for 3-hydroxyoctadecanoyl carnitine
    EX_3tdcrn_ 	 exchange reaction for 3-hydroxy-tetradecanoyl carnitine
    EX_3ttetddcoacrn_ 	 exchange reaction for 3-hydroxy trans5,8tetradecadienoyl carnitine
    EX_c51crn_ 	 exchange reaction for tiglyl carnitine
    EX_c6crn_ 	 exchange reaction for hexanoyl carnitine
    EX_c8crn_ 	 exchange reaction for octanoyl carnitine
    EX_ddecrn_ 	 exchange reaction for lauroyl carnitine
    EX_ivcrn_ 	 exchange reaction for isovaleryl carnitine
    EX_4hpro_LPAREN_e_RPAREN_ 	 exchange reaction for hydroxy proline
    EX_carn_LPAREN_e_RPAREN_ 	 exchange reaction for carnosine
    EX_sbt_DASH_d_LPAREN_e_RPAREN_ 	 exchange reaction for D-Sorbitol
    EX_sfcys_LPAREN_e_RPAREN_ 	 exchange reaction for sulfocysteine
    EX_fmn_LPAREN_e_RPAREN_ 	 exchange reaction for FMN
    EX_q10_LPAREN_e_RPAREN_ 	 exchange reaction for ubiquinone
    EX_q10h2_LPAREN_e_RPAREN_ 	 exchange reaction for ubiquinol
    EX_3mob_LPAREN_e_RPAREN_ 	 Exchange of 3-methyl-2-oxobutanoate
    EX_3mop_LPAREN_e_RPAREN_ 	 Exchange of 3-methyl-2-oxopentanoate
    EX_4mop_LPAREN_e_RPAREN_ 	 Exchange of 4-methyl-2-oxopentanoate
    EX_5oxpro_LPAREN_e_RPAREN_ 	 Exchange of 5-oxoprolinate
    EX_ahcys_LPAREN_e_RPAREN_ 	 Exchange of S-adenosyl-L-homocysteine
    EX_mal_L_LPAREN_e_RPAREN_ 	 Exchange of L-malate
    EX_fum_LPAREN_e_RPAREN_ 	 Exchange of Fumarate
    EX_xmp_LPAREN_e_RPAREN_ 	 Exchange of Xanthosine 5-phosphate
    EX_dcmp_LPAREN_e_RPAREN_ 	 Exchange of dCMP
    EX_HC00342_LPAREN_e_RPAREN_ 	 exhange of cis-aconitate
    EX_glyc3p_LPAREN_e_RPAREN_ 	 exhange of glycerol 3-phosphate
    biomass_reaction 	 Generic human biomass reaction
    biomass_protein 	 protein component of biomass
    biomass_DNA 	 DNA component of biomass
    biomass_RNA 	 RNA component of biomass
    biomass_carbohydrate 	 carbohydrate component of biomass
    biomass_lipid 	 lipid component of biomass
    biomass_other 	 other component of biomass
    BETBGTtc 	 betaine transport by BGT
    CRNATBtc 	 carnitine transport by ATB0
    H2OGLYAQPt 	 water and glycerol transport by AQP
    HDCAFAPMtc 	 hexadecanoate transport by FAT
    ILEB0AT2tc 	 isoleucine transport by B0AT2
    LEUB0AT2tc 	 leucine transport by B0AT2
    LGNCFATtc 	 lignocerate transport by FAT
    OCDCAFAPMtc 	 octadecanoate tranport by FAT
    PROB0AT2tc 	 proline transport by B0AT2
    VALB0AT2tc 	 valine transport by B0AT2
    3HPVSTETCOAhcm 	 beta oxidation of 3''-S-hydroxy-pravastatin-CoA to tatranor-CoA derivative in hepatocytes, mitochondria
    3HPVSTETtev 	 exit of 3''-S-hydroxy-pravastatin-tetranor into hepatic vein
    3HPVStep 	 3''-hydroxy pravastatin exit into portal blood
    AM4N9CShc 	 demethylation of AM9 to AM4N9 in hepatocytes
    AM4N9CStev 	 efflux of AM4N9-cyclosporine into hepatic vein
    AM4NCShr 	 demethylation of cyclosporine to AM4N in hepatocytes
    AM4NCStep 	 efflux of AM4N (cyclosporine) into portal blood
    AM9CSAtep 	 efflux of AM9 (cyclosporine) into portal blood
    CRVSM1hr 	 beta-glucuronidase of cerivastatin-M1-glucuronide to M1 form
    CRVSM23hr 	 beta-glucuronidase of cerivastatin-M23-glucuronide to M23 form
    CSAtd 	 uptake of cyclosporine by the enterocytes
    CVM1GLUChc 	 glucuronidation of cerivastatin-M1
    CVM23GLUChc 	 glucuronidation of cerivastatin-M23
    EPOXTAChr 	 formation of epoxy derivative of tacrolimus in hepatocytes
    EPOXTACtev 	 exit of 31-O-Desmethyl,19-Hydroxy,37, 39-Epoxy-tacrolimus into hepatic vein
    EX_3hpvs_LPAREN_e_RPAREN_ 	 3'-S-hydroxy-pravastatin exchange
    EX_3hpvstet_LPAREN_e_RPAREN_ 	 3'-S-hydroxy-pravastatin-tetranor exchange
    EX_am4n9cs_LPAREN_e_RPAREN_ 	 AM4N9 (cyclosporine) exchange
    EX_am4ncs_LPAREN_e_RPAREN_ 	 AM4N (cyclosporine) exchange
    EX_am9csa_LPAREN_e_RPAREN_ 	 AM9 (cyclosporine) exchange
    EX_csa_LPAREN_e_RPAREN_ 	 cyclosporine exchange
    EX_epoxtac_LPAREN_e_RPAREN_ 	 31-O-Desmethyl,19-Hydroxy,37, 39-Epoxy-tacrolimus exchange
    EX_tacr_LPAREN_e_RPAREN_ 	 Tacrolimus exchange
    SMVACIDhep 	 beta-glucuronidation of simvastatin-acyl-glucuronide to simvastatin dihydroxy acid form
    SMVGLUChep 	 glucuronidation of the open acid form of simvastatin
    TACRDtsc 	 uptake of tacrolimus by enterocytes
    TMDOATPtsc 	 uptake of torasemide into enterocytes via antiport
    TMDOATtev 	 efflux of Torasemide into hepatic vein
    TMDtd 	 uptake of torasemide into enterocytes via diffusion
    3HPVSCOAitm 	 3HPVSCOAitm
    3HPVSTETCOAitm 	 3HPVSTETCOAitm
    EPOXTACitr 	 EPOXTACitr
    TACRitr 	 TACRitr
    3HPVSCOAhc 	 3HPVSCOAhc
    3HPVSTEThc 	 3HPVSTEThc
    EX_hx_LPAREN_e_RPAREN_ 	 exchange reaction for hexanoate (n-C6:0)
    FACOAL120i 	 fatty-acid--CoA ligase (dodecanoate, n-C12:0)
    FACOAL60i 	 fatty-acid--CoA ligase (hexanoate, n-C6:0)
    HXt 	 hexanoate transport via diffusion
    AIRCr_PRASCS 	 phosphoribosylaminoimidazole carboxylase / phosphoribosylaminoimidazolesuccinocarboxamide synthase
    FAOXC10080x 	 R_FAOXC10080x
    FAOXC101_3Em 	 R_FAOXC101_3Em
    FAOXC101_4Zx 	 R_FAOXC101_4Zx
    FAOXC102_4Z_7Zm 	 R_FAOXC102_4Z_7Zm
    FAOXC120100m 	 R_FAOXC120100m
    FAOXC120100x 	 R_FAOXC120100x
    FAOXC122_3Z_6Zx 	 R_FAOXC122_3Z_6Zx
    FAOXC123_3Z_6Z_9Zm 	 R_FAOXC123_3Z_6Z_9Zm
    FAOXC140120x 	 R_FAOXC140120x
    FAOXC142_5Z_8Zx 	 R_FAOXC142_5Z_8Zx
    FAOXC143_5Z_8Z_11Zm 	 R_FAOXC143_5Z_8Z_11Zm
    FAOXC161140m 	 R_FAOXC161140m
    FAOXC163_4Z_7Z_10Zx 	 R_FAOXC163_4Z_7Z_10Zx
    FAOXC164_4Z_7Z_10Z_13Zm 	 R_FAOXC164_4Z_7Z_10Z_13Zm
    FAOXC184_3Z_6Z_9Z_12Zx 	 R_FAOXC184_3Z_6Z_9Z_12Zx
    FAOXC184_6Z_9Z_12Z_15Zm 	 R_FAOXC184_6Z_9Z_12Z_15Zm
    FAOXC185_3Z_6Z_9Z_12Z_15Zm 	 R_FAOXC185_3Z_6Z_9Z_12Z_15Zm
    FAOXC204_5Z_8Z_11Z_14Zx 	 R_FAOXC204_5Z_8Z_11Z_14Zx
    FAOXC205_5Z_8Z_11Z_14Z_17Zm 	 R_FAOXC205_5Z_8Z_11Z_14Z_17Zm
    FAOXC4020m 	 R_FAOXC4020m
    FAOXC61_3Zm 	 R_FAOXC61_3Zm
    FAOXC8060m 	 R_FAOXC8060m
    FAOXC81_5Zm 	 R_FAOXC81_5Zm



```python
import pandas as pd


df = 
df.to_csv(index=False)
```


    ---------------------------------------------------------------------------

    ValueError                                Traceback (most recent call last)

    <ipython-input-24-5004c5668cb6> in <module>()
          2 
          3 
    ----> 4 df = pd.DataFrame(scc_model.reactions.index)
          5 df.to_csv(index=False)


    /opt/conda/lib/python3.5/site-packages/pandas/core/frame.py in __init__(self, data, index, columns, dtype, copy)
        352                                          copy=False)
        353             else:
    --> 354                 raise ValueError('DataFrame constructor not properly called!')
        355 
        356         NDFrame.__init__(self, mgr, fastpath=True)


    ValueError: DataFrame constructor not properly called!



```python
list = ["RE2319R", "C141ACBP", "TYRTAm", "HDCACBP", "DPGM", "RE1635R", "3HBCOARc", "FAEL184", "GLXO1", "r0173", "FBA2", "r0570", "NaKt", "EHGLATm", "RE3401M", "OCDCAFAPMtc", "MCLOR", "NTD5", "SFGTH", "ARGtiDF", "SERt4", "r0801", "r0358", "r0354", "biomass_RNA", "HDECAACBP", "RE3386M", "PPM", "ORNtiDF", "NTD8", "ECOAH12m", "PGK", "NT5C", "GLNt4", "RE3587N", "HTDCACBP", "NTD4", "RE1818R", "FBA4", "RE3110C", "RE3192M", "RE3120C", "GTHPm", "FAEL204", "GTHPe", "C162ACBP", "3SALATAim", "NTD6", "ACITL", "r0047", "ALAt4", "RE3398M", "HOCDACBP", "RE1520M", "biomass_carbohydrate", "ECOAH1m", "DM_atp_c_", "ALCD21_L", "PMANM", "r0575", "IMPD", "HEXCOAACBP", "SPODM", "biomass_lipid", "biomass_protein", "r0407", "RE2514C", "RE3344M", "ADK3m", "HIVCACBP", "r0360", "SPODMx", "RE3390M", "RE3241C", "RE3563M", "PPA", "2MB2COAc", "RE3248M", "r0739", "r0170", "HEX10", "GALU", "HEX1", "RE0565C", "RE3225C", "r0361", "SPODMn", "HDCAFAPMtc", "IVCOAACBP", "LYStiDF", "LALDO", "RE1816R", "HDDACBP", "RDH1a", "DUTPDPn", "NTD7", "EHGLAT2m", "RE1521M", "RE3194M", "ENO", "NTD10", "C4STMO1r", "NTD1", "NTD3", "2HBO", "34DHOXPEGOX", "RE1062M", "NADKm", "RE3241R", "OCTDEC2ACBP", "PHETA1m", "RE2920M", "r0355", "RE3521R", "ATPtm", "PRPNCOAHYDm", "HISt4", "NTD9", "r0656", "RAI3", "r0357", "ASPTAm", "MDHm", "ASNt4", "HOMt4", "RE3532R", "HEX4", "DUTPDPm", "RE3447M", "r0364", "FBA", "SQLSr", "LDH_L", "RE3557R", "GLYt4", "ALCD1", "RE3496N", "NTD2", "RE2514L", "GLUTCOAACBP", "RE3193M", "FAEL183", "RE3533R", "PRDX", "HDECEACBP", "LTC4Sr", "PHEt4", "LCYSTATm", "LGTHL", "RE3004M", "RE3224C", "HEXDICOAACBPx", "biomass_DNA", "CYSTAm", "r1378", "RE1523M", "RE2919M", "biomass_other", "PGM", "ALCD22_L", "RE3226C", "r0604", "PGMT", "RE2318R", "NTD11", "C142ACBP", "DASCBR", "RE3554R", "r0363", "GTHP", "r0009", "HOCTDACBP", "FALDH", "PIt2m", "ECOAH9m", "RE3534R", "H2OGLYAQPt", "THRt4", "RE1522M", "HEXDICOAACBP", "SLDxm", "MDH", "RE3337M", "HPYRR2x", "LDH_Lm", "RE1815R", "RE1817R", "C4STMO2r", "LGNCFATtc", "HEX7", "biomass_reaction", "GAPD", "CSm", "ALCD2if", "TPI", "SLDx", "r0660", "FAEL205", "GALNTg", "RE2921M", "DPGase"]

```


```python
for rx in list:
    print(Recon2.reactions.get_by_id(rx).gene_reaction_rule)
```

    HGNC:7061 or HGNC:7063 or HGNC:7064
    HGNC:2690
    HGNC:11573 or HGNC:4433
    HGNC:2690
    HGNC:1093 or HGNC:8888 or HGNC:8889
    HGNC:7061 or HGNC:7063 or HGNC:7064
    HGNC:2690
    HGNC:14415 or HGNC:14416 or HGNC:15829 or HGNC:21308
    HGNC:19708 or HGNC:21481 or HGNC:28335 or HGNC:30866 or (HGNC:6535 and HGNC:6541) or HGNC:6535 or HGNC:6541 or HGNC:6544
    HGNC:6535 or HGNC:6541
    HGNC:414 or HGNC:417 or HGNC:418
    HGNC:8906
    (HGNC:14073 and HGNC:804) or (HGNC:14073 and HGNC:808) or (HGNC:799 and HGNC:804) or (HGNC:799 and HGNC:805) or (HGNC:799 and HGNC:806) or (HGNC:799 and HGNC:808) or (HGNC:800 and HGNC:804) or (HGNC:800 and HGNC:805) or (HGNC:800 and HGNC:808) or (HGNC:801 and HGNC:804) or (HGNC:801 and HGNC:805)
    HGNC:4433
    HGNC:3151 or HGNC:3247 or HGNC:4801 or HGNC:4803
    HGNC:1663 or HGNC:4433
    HGNC:21481 or HGNC:28335 or HGNC:30866 or (HGNC:6535 and HGNC:6541) or HGNC:6535 or HGNC:6541 or HGNC:6544
    HGNC:17144 or HGNC:17818 or HGNC:17820 or HGNC:8021 or HGNC:8022
    HGNC:3465
    HGNC:11057 or HGNC:11060 or HGNC:11061 or HGNC:14679
    HGNC:11047 or HGNC:13447 or HGNC:13448 or (HGNC:13557 and HGNC:27960) or HGNC:14679 or (HGNC:27960 and HGNC:29437)
    HGNC:10992
    HGNC:4922 or HGNC:4923 or HGNC:4925
    HGNC:4195 or HGNC:4922 or HGNC:4923 or HGNC:4925
    
    HGNC:2690
    HGNC:3151 or HGNC:3247 or HGNC:4801 or HGNC:4803
    HGNC:8905 or HGNC:8906
    HGNC:11057 or HGNC:11060 or HGNC:11061
    HGNC:17144 or HGNC:17818 or HGNC:17820 or HGNC:8021 or HGNC:8022
    HGNC:3151 or (HGNC:4801 and HGNC:4803)
    HGNC:8896 or HGNC:8898
    HGNC:17144 or HGNC:17818 or HGNC:17820 or HGNC:8021 or HGNC:8022
    HGNC:11047 or HGNC:13447 or HGNC:13448 or (HGNC:13557 and HGNC:27960) or HGNC:14679 or (HGNC:27960 and HGNC:29437)
    HGNC:16049
    HGNC:2690
    HGNC:17144 or HGNC:17818 or HGNC:17820 or HGNC:8021 or HGNC:8022
    HGNC:7061 or HGNC:7063 or HGNC:7064
    HGNC:414 or HGNC:417 or HGNC:418
    HGNC:14415 or HGNC:14416 or HGNC:15829 or HGNC:21308
    HGNC:3151 or HGNC:3247 or HGNC:4801 or HGNC:4803
    HGNC:14415 or HGNC:14416 or HGNC:15829 or HGNC:21308
    HGNC:4553 or HGNC:4556 or HGNC:9354
    HGNC:14415 or HGNC:14416 or HGNC:15829 or HGNC:21308
    HGNC:4553 or HGNC:4554 or HGNC:4555 or HGNC:4558 or HGNC:4559
    HGNC:2690
    HGNC:4433
    HGNC:17144 or HGNC:17818 or HGNC:17820 or HGNC:8021 or HGNC:8022
    HGNC:115
    HGNC:15769 or HGNC:17144 or HGNC:17818 or HGNC:17820 or HGNC:17819 or HGNC:38831 or HGNC:8021 or HGNC:8022
    HGNC:11047 or HGNC:11057 or HGNC:11060 or HGNC:11061 or HGNC:13447 or HGNC:13448 or (HGNC:13557 and HGNC:27960) or HGNC:14679 or (HGNC:27960 and HGNC:29437)
    HGNC:3151 or HGNC:3247 or HGNC:4801 or HGNC:4803
    HGNC:2690
    HGNC:3151 or HGNC:3247 or HGNC:4801 or HGNC:4803
    
    HGNC:3151 or (HGNC:4801 and HGNC:4803) or HGNC:890
    
    HGNC:16354 or (HGNC:249 and HGNC:250) or (HGNC:249 and HGNC:251) or HGNC:249 or (HGNC:250 and HGNC:251) or HGNC:250 or HGNC:251 or HGNC:252 or HGNC:253 or HGNC:255 or HGNC:256 or HGNC:28697
    HGNC:8906 or HGNC:9114 or HGNC:9115
    HGNC:3629
    HGNC:6052 or HGNC:6053
    HGNC:2690
    HGNC:11179
    
    
    HGNC:414 or HGNC:417 or HGNC:418
    HGNC:12015 or HGNC:16753
    HGNC:3151 or HGNC:3247 or HGNC:4801 or HGNC:4803
    HGNC:17376
    HGNC:2690
    HGNC:4922 or HGNC:4923 or HGNC:4925
    HGNC:11179
    HGNC:3151 or HGNC:3247 or HGNC:4801 or HGNC:4803
    HGNC:14415 or HGNC:14416 or HGNC:15829 or HGNC:21308
    HGNC:3151 or HGNC:3247 or HGNC:4801 or HGNC:4803
    HGNC:30042 or HGNC:9226
    HGNC:2690
    HGNC:3151 or HGNC:3247 or HGNC:4801 or HGNC:4803
    (HGNC:249 and HGNC:250 and HGNC:251) or HGNC:252 or HGNC:253 or HGNC:255 or HGNC:256
    HGNC:3629
    HGNC:23302 or HGNC:4195 or HGNC:4922 or HGNC:4923 or HGNC:4925
    HGNC:12527 or HGNC:12527
    HGNC:23302 or HGNC:4195 or HGNC:4922 or HGNC:4923 or HGNC:4925
    HGNC:14415 or HGNC:14416 or HGNC:15829 or HGNC:21308
    HGNC:14415 or HGNC:14416 or HGNC:15829 or HGNC:21308
    HGNC:4922 or HGNC:4923 or HGNC:4925
    HGNC:11179
    HGNC:1663 or HGNC:4433
    HGNC:2690
    HGNC:11057 or HGNC:11060 or HGNC:11061 or HGNC:14679
    HGNC:253
    HGNC:7061 or HGNC:7063 or HGNC:7064
    HGNC:2690
    HGNC:14423 or HGNC:17964 or HGNC:19975 or HGNC:19977 or HGNC:19978 or HGNC:19979 or HGNC:30311
    HGNC:3078
    HGNC:17144 or HGNC:17818 or HGNC:17820 or HGNC:8022
    HGNC:4433
    HGNC:3151 or HGNC:3247 or HGNC:4801 or HGNC:4803
    HGNC:3151 or HGNC:3247 or HGNC:4801 or HGNC:4803
    HGNC:3350 or HGNC:3354 or HGNC:3354 or HGNC:3353
    HGNC:17144 or HGNC:17818 or HGNC:17820 or HGNC:8021 or HGNC:8022
    HGNC:10545
    HGNC:17144 or HGNC:17820
    HGNC:17144 or HGNC:17818 or HGNC:17820 or HGNC:4678 or HGNC:8021
    HGNC:21481 or HGNC:28335 or HGNC:30866 or (HGNC:6535 and HGNC:6541) or HGNC:6535 or HGNC:6541 or HGNC:6544
    HGNC:16354 or (HGNC:249 and HGNC:250) or (HGNC:249 and HGNC:251) or (HGNC:250 and HGNC:251) or HGNC:252 or HGNC:253 or HGNC:255 or HGNC:256
    HGNC:16058
    HGNC:26404
    HGNC:14416 or HGNC:14418 or HGNC:15829 or HGNC:21308
    HGNC:2690
    HGNC:11573 or HGNC:4433
    HGNC:3151 or HGNC:3247 or HGNC:4801 or HGNC:4803
    HGNC:4195 or HGNC:4922 or HGNC:4923 or HGNC:4925
    HGNC:7061 or HGNC:7063 or HGNC:7064
    HGNC:10990 or HGNC:10991 or HGNC:10992
    HGNC:3151 or HGNC:3247 or (HGNC:4801 and HGNC:4803) or HGNC:890
    HGNC:11047 or HGNC:13448
    HGNC:17818 or HGNC:17820 or HGNC:8021 or HGNC:8022
    HGNC:2698
    HGNC:4632 or HGNC:4632 or HGNC:4638 or HGNC:4626
    HGNC:4922 or HGNC:4923 or HGNC:4925
    HGNC:4433
    HGNC:6970 or HGNC:6971
    HGNC:11047 or HGNC:13447 or HGNC:13448 or (HGNC:13557 and HGNC:27960) or HGNC:14679 or (HGNC:27960 and HGNC:29437)
    HGNC:11057 or HGNC:11060 or HGNC:11061
    HGNC:7061 or HGNC:7063 or HGNC:7064
    HGNC:23302 or HGNC:4195 or HGNC:4922 or HGNC:4923 or HGNC:4925
    HGNC:3078 or HGNC:6176
    HGNC:3151 or HGNC:3247 or HGNC:4801 or HGNC:4803
    HGNC:4922 or HGNC:4923 or HGNC:4925
    HGNC:414 or HGNC:417 or HGNC:418
    HGNC:3629
    HGNC:21481 or HGNC:28335 or HGNC:30866 or (HGNC:6535 and HGNC:6541) or HGNC:6535 or HGNC:6541 or HGNC:6544
    HGNC:7061 or HGNC:7063 or HGNC:7064
    HGNC:11047 or HGNC:13447 or HGNC:13448 or (HGNC:13557 and HGNC:27960) or HGNC:14679 or (HGNC:27960 and HGNC:29437)
    HGNC:16354 or (HGNC:249 and HGNC:250) or (HGNC:249 and HGNC:251) or HGNC:249 or (HGNC:250 and HGNC:251) or HGNC:250 or HGNC:251 or HGNC:252 or HGNC:253 or HGNC:255 or HGNC:256 or HGNC:28697
    HGNC:16049
    HGNC:17820
    HGNC:16753 or HGNC:7218
    HGNC:2690
    HGNC:3151 or HGNC:3247 or HGNC:4801 or HGNC:4803
    HGNC:14415 or HGNC:14416 or HGNC:15829 or HGNC:21308
    HGNC:7061 or HGNC:7063 or HGNC:7064
    HGNC:16753 or HGNC:3423
    HGNC:2690
    HGNC:6719 or HGNC:7063 or HGNC:7064
    HGNC:11047 or HGNC:11057 or HGNC:11060 or HGNC:11061 or HGNC:13447 or HGNC:13448 or (HGNC:13557 and HGNC:27960) or HGNC:14679 or (HGNC:27960 and HGNC:29437)
    HGNC:4433
    HGNC:4323
    HGNC:3151 or HGNC:3247 or HGNC:4801 or HGNC:4803
    HGNC:14415 or HGNC:14416 or HGNC:15829 or HGNC:21308
    HGNC:2690
    
    HGNC:4433
    HGNC:253
    HGNC:3151 or HGNC:3247 or HGNC:4801 or HGNC:4803
    HGNC:3151 or HGNC:3247 or HGNC:4801 or HGNC:4803
    
    HGNC:1093 or HGNC:8888 or HGNC:8889
    HGNC:16354 or (HGNC:249 and HGNC:250) or (HGNC:249 and HGNC:251) or HGNC:249 or (HGNC:250 and HGNC:251) or HGNC:250 or HGNC:251 or HGNC:252 or HGNC:253 or HGNC:255 or HGNC:256 or HGNC:28697
    HGNC:14415 or HGNC:14416 or HGNC:15829 or HGNC:21308
    HGNC:2698
    HGNC:8905 or HGNC:8906
    HGNC:7061 or HGNC:7063 or HGNC:7064
    HGNC:17818 or HGNC:17820 or HGNC:8021 or HGNC:8022
    HGNC:2690
    HGNC:16065 or HGNC:4330 or HGNC:13312
    HGNC:7061 or HGNC:7063 or HGNC:7064
    HGNC:4922 or HGNC:4923 or HGNC:4925
    HGNC:4553 or HGNC:4554 or HGNC:4556 or HGNC:9352 or HGNC:9353
    HGNC:28883 or HGNC:30042 or HGNC:9226
    HGNC:2690
    HGNC:253
    HGNC:10989
    HGNC:3151 or (HGNC:4801 and HGNC:4803)
    HGNC:7061 or HGNC:7063 or HGNC:7064
    HGNC:16029 or HGNC:636 or HGNC:640 or HGNC:643
    HGNC:11047 or HGNC:13447 or HGNC:13448 or (HGNC:13557 and HGNC:27960) or (HGNC:27960 and HGNC:29437)
    HGNC:3151 or HGNC:3247 or HGNC:4801 or HGNC:4803
    HGNC:2690
    HGNC:6971
    HGNC:17836 or HGNC:6970
    HGNC:3151 or HGNC:3247 or HGNC:4801 or HGNC:4803
    HGNC:19708 or HGNC:21481 or HGNC:28335 or (HGNC:6535 and HGNC:6541) or HGNC:6535 or HGNC:6541 or HGNC:6544
    HGNC:19708 or HGNC:6535 or HGNC:6541
    HGNC:7061 or HGNC:7063 or HGNC:7064
    HGNC:7061 or HGNC:7063 or HGNC:7064
    HGNC:10545 or HGNC:13398
    HGNC:1663 or HGNC:4433
    HGNC:23302 or HGNC:4195 or HGNC:4922 or HGNC:4923 or HGNC:4925
    
    HGNC:24864 or HGNC:4141
    HGNC:2422
    HGNC:16354 or (HGNC:249 and HGNC:250) or (HGNC:249 and HGNC:251) or HGNC:249 or (HGNC:250 and HGNC:251) or HGNC:250 or HGNC:251 or HGNC:252 or HGNC:253 or HGNC:255 or HGNC:256 or HGNC:28697
    HGNC:12009
    HGNC:17836 or HGNC:6970
    (HGNC:4801 and HGNC:4803) or HGNC:3151
    HGNC:14415 or HGNC:14416 or HGNC:15829 or HGNC:21308
    HGNC:19873 or HGNC:19875 or HGNC:19877 or HGNC:21531 or HGNC:22946 or HGNC:23242 or HGNC:4123 or HGNC:4124 or HGNC:4125 or HGNC:4126 or HGNC:4127 or HGNC:4128 or HGNC:4129 or HGNC:4130 or HGNC:4131
    HGNC:3151 or HGNC:3247 or HGNC:4801 or HGNC:4803
    HGNC:1093 or HGNC:8888 or HGNC:8889



```python
Recon2.metabolites.HC00250_e.name
```




    'hydrosulfide'




```python
Recon2.reactions.PSSA1_hs.gene_reaction_rule
```




    'HGNC:9587'




```python
#Biomass bounds
for rxex in hela_model.exchanges:
    rxex.bounds=(-0.01,1)
    
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
hela_model.reactions.EX_chol_LPAREN_e_RPAREN_.lower_bound=-1
hela_model.reactions.EX_fe2_LPAREN_e_RPAREN_.lower_bound=-1
hela_model.reactions.EX_fol_LPAREN_e_RPAREN_.lower_single_gene_deletionbound=-1
hela_model.reactions.EX_h2o_LPAREN_e_RPAREN_.lower_bound=-10
hela_model.reactions.EX_inost_LPAREN_e_RPAREN_.lower_bound=-1
hela_model.reactions.EX_pi_LPAREN_e_RPAREN_.lower_bound=-10
hela_model.reactions.EX_pydxn_LPAREN_e_RPAREN_.lower_bound=-1
hela_model.reactions.EX_ribflv_LPAREN_e_RPAREN_.lower_bound=-1
hela_model.reactions.EX_so4_LPAREN_e_RPAREN_.lower_bound=-1

#EXTRA
hela_model.reactions.EX_o2_LPAREN_e_RPAREN_.bounds=(-1000,0)

hela_model.reactions.biomass_reaction.bounds=(0.0,1000)
```


```python
#Biomass bounds
for rxex in hacat_model.exchanges:
    rxex.bounds=(-0.01,1)
    
### DMEM 6429 medium
#### AminoAcids

hacat_model.reactions.EX_ala_L_LPAREN_e_RPAREN_.lower_bound=-1
hacat_model.reactions.EX_arg_L_LPAREN_e_RPAREN_.lower_bound=-1
hacat_model.reactions.EX_asn_L_LPAREN_e_RPAREN_.lower_bound=-1
hacat_model.reactions.EX_asp_L_LPAREN_e_RPAREN_.lower_bound=-1
hacat_model.reactions.EX_cys_L_LPAREN_e_RPAREN_.lower_bound=-1
hacat_model.reactions.EX_gln_L_LPAREN_e_RPAREN_.lower_bound=-1
hacat_model.reactions.EX_glu_L_LPAREN_e_RPAREN_.lower_bound=-1
hacat_model.reactions.EX_gly_LPAREN_e_RPAREN_.lower_bound=-1
hacat_model.reactions.EX_his_L_LPAREN_e_RPAREN_.lower_bound=-1
hacat_model.reactions.EX_ile_L_LPAREN_e_RPAREN_.lower_bound=-1
hacat_model.reactions.EX_leu_L_LPAREN_e_RPAREN_.lower_bound=-1
hacat_model.reactions.EX_lys_L_LPAREN_e_RPAREN_.lower_bound=-1 
hacat_model.reactions.EX_met_L_LPAREN_e_RPAREN_.lower_bound=-1
hacat_model.reactions.EX_phe_L_LPAREN_e_RPAREN_.lower_bound=-1
hacat_model.reactions.EX_pro_L_LPAREN_e_RPAREN_.lower_bound=-1
hacat_model.reactions.EX_ser_L_LPAREN_e_RPAREN_.lower_bound=-1
hacat_model.reactions.EX_thr_L_LPAREN_e_RPAREN_.lower_bound=-1
hacat_model.reactions.EX_trp_L_LPAREN_e_RPAREN_.lower_bound=-1
hacat_model.reactions.EX_tyr_L_LPAREN_e_RPAREN_.lower_bound=-1
hacat_model.reactions.EX_val_L_LPAREN_e_RPAREN_.lower_bound=-1



### DMEM 6429 medium
#Carbon Sources
hacat_model.reactions.EX_glc_LPAREN_e_RPAREN_.lower_bound=-4.5
hacat_model.reactions.EX_pyr_LPAREN_e_RPAREN_.lower_bound= -1

### DMEM 6429 medium
#Minerals, vitamins, other
hacat_model.reactions.EX_chol_LPAREN_e_RPAREN_.lower_bound=-1
hacat_model.reactions.EX_fe2_LPAREN_e_RPAREN_.lower_bound=-1
hacat_model.reactions.EX_fol_LPAREN_e_RPAREN_.lower_single_gene_deletionbound=-1
hacat_model.reactions.EX_h2o_LPAREN_e_RPAREN_.lower_bound=-10
#h2s
hacat_model.reactions.EX_inost_LPAREN_e_RPAREN_.lower_bound=-1
hacat_model.reactions.EX_pi_LPAREN_e_RPAREN_.lower_bound=-10
hacat_model.reactions.EX_ribflv_LPAREN_e_RPAREN_.lower_bound=-1
hacat_model.reactions.EX_so4_LPAREN_e_RPAREN_.lower_bound=-1

#EXTRA
hacat_model.reactions.EX_o2_LPAREN_e_RPAREN_.bounds=(-1000,0)

hacat_model.reactions.biomass_reaction.lower_bound=0.000
hacat_model.reactions.biomass_reaction.upper_bound=1000

```


```python
hela_fba  = hela_model.optimize()
hacat_fba = hacat_model.optimize()
```

For each flux in HeLa we reduce the flux to 0.001 and check infleasible fluxes


```python
import pandas as pd

def get_drugable_targets(normal_Model, disease_Model,  eps=0.001):
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
        print(rx)

        nbounds=nmodel.reactions.get_by_id(rx).bounds
        dbounds=dmodel.reactions.get_by_id(rx).bounds
        
        nmodel.reactions.get_by_id(rx).bounds=(-eps,eps)
        dmodel.reactions.get_by_id(rx).bounds=(-eps,eps)
                
        nflx1=nmodel.optimize().f
        dflx1=dmodel.optimize().f
        
        results[rx]={}
        results[rx]["norm_flux"]=nflx1
        results[rx]["dise_flux"]=dflx1
 
        results[rx]["norm_prolif_ratio"]=nflx1/nflx0
        results[rx]["dise_prolif_ratio"]=dflx1/dflx0
        
        
        results[rx]["norm_dise_ratio"]=(nflx1/nflx0)/(dflx1/dflx0)
        nmodel.reactions.get_by_id(rx).bounds=nbounds
        dmodel.reactions.get_by_id(rx).bounds=dbounds

        
    for rx in unique_Nrx:
        print(rx)
        
        nbounds=nmodel.reactions.get_by_id(rx).bounds
        
        nmodel.reactions.get_by_id(rx).bounds=(-eps,eps)
        nflx1=nmodel.optimize().f
        results[rx]={}

        results[rx]["norm_flux"]=nflx1
        results[rx]["dise_flux"]=dflx0
 
        results[rx]["norm_prolif_ratio"]=nflx1/nflx0
        results[rx]["dise_prolif_ratio"]=dflx0
        
        results[rx]["norm_dise_ratio"]=nflx1/nflx0
        nmodel.reactions.get_by_id(rx).bounds=nbounds

    for rx in unique_Drx:
        print(rx)
        
        dbounds=dmodel.reactions.get_by_id(rx).bounds
    
        dmodel.reactions.get_by_id(rx).bounds=(-eps,eps)
        dflx1=dmodel.optimize().f
        
        results[rx]={}
        results[rx]["norm_flux"]=nflx0
        results[rx]["dise_flux"]=dflx1
 
        results[rx]["norm_prolif_ratio"]=nflx0
        results[rx]["dise_prolif_ratio"]=dflx1/dflx0
        
        results[rx]["norm_dise_ratio"]=1/(dflx1/dflx0)
        dmodel.reactions.get_by_id(rx).bounds=dbounds


    return(pd.DataFrame(results).transpose())
       
  
               

```


```python
R=get_drugable_targets(hacat_model, hela_model)
```


    ---------------------------------------------------------------------------

    NameError                                 Traceback (most recent call last)

    <ipython-input-5-272e3856561b> in <module>()
    ----> 1 R=get_drugable_targets(hacat_model, hela_model)
    

    NameError: name 'hacat_model' is not defined



```python
R[R["dise_prolif_ratio"]< 1 ].sort_values(by="norm_dise_ratio", ascending=False)
```


    ---------------------------------------------------------------------------

    NameError                                 Traceback (most recent call last)

    <ipython-input-3-78e08c28ed6e> in <module>()
    ----> 1 R[  R["dise_prolif_ratio"]< 1 ].sort_values(by="norm_dise_ratio", ascending=False)
    

    NameError: name 'R' is not defined



```python
hela_model.reactions.PMEVKc.gene_reaction_rule
```


    ---------------------------------------------------------------------------

    NameError                                 Traceback (most recent call last)

    <ipython-input-6-1554ae659e02> in <module>()
    ----> 1 hela_model.genes.get_by_id("HGNC:9068")
    

    NameError: name 'hela_model' is not defined



```python
%matplotlib inline

import seaborn as sns
import matplotlib.pyplot as plt
import numpy as np
from matplotlib.colors import ListedColormap

R["FoldChange"]=np.log2(R["norm_dise_ratio"])

#sns.boxplot(x="FoldChange", data=R, )


plt.scatter( np.log2(R["dise_prolif_ratio"]) ,np.log2(R["norm_prolif_ratio"]))
```




    <matplotlib.collections.PathCollection at 0x7f6f84dd7278>




![png](output_17_1.png)



```python
import pandas as pd

f1=cobra.flux_analysis.pfba(hela_model).fluxes
f2=cobra.flux_analysis.pfba(hacat_model).fluxes

f1=hela_model.optimize().fluxes
f2=hacat_model.optimize().fluxes


df=pd.DataFrame(dict(hela_model = f1, hacat_model = f2)).reset_index()          

```


```python
df.to_csv("fba_hela_hacat_fluxes.csv", header=True)

```


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
      <th>index</th>
      <th>hacat_model</th>
      <th>hela_model</th>
    </tr>
  </thead>
  <tbody>
    <tr>
      <th>0</th>
      <td>13DAMPPOX</td>
      <td>0.000000e+00</td>
      <td>0.000000e+00</td>
    </tr>
    <tr>
      <th>1</th>
      <td>1PPDCRp</td>
      <td>0.000000e+00</td>
      <td>0.000000e+00</td>
    </tr>
    <tr>
      <th>2</th>
      <td>2AMADPTm</td>
      <td>0.000000e+00</td>
      <td>0.000000e+00</td>
    </tr>
    <tr>
      <th>3</th>
      <td>2DR1PP</td>
      <td>NaN</td>
      <td>0.000000e+00</td>
    </tr>
    <tr>
      <th>4</th>
      <td>2HBO</td>
      <td>NaN</td>
      <td>0.000000e+00</td>
    </tr>
    <tr>
      <th>5</th>
      <td>2HBt2</td>
      <td>NaN</td>
      <td>0.000000e+00</td>
    </tr>
    <tr>
      <th>6</th>
      <td>2MB2COAc</td>
      <td>0.000000e+00</td>
      <td>4.695740e-03</td>
    </tr>
    <tr>
      <th>7</th>
      <td>2MCITt</td>
      <td>0.000000e+00</td>
      <td>NaN</td>
    </tr>
    <tr>
      <th>8</th>
      <td>2OXOADOXm</td>
      <td>0.000000e+00</td>
      <td>0.000000e+00</td>
    </tr>
    <tr>
      <th>9</th>
      <td>2OXOADPTm</td>
      <td>0.000000e+00</td>
      <td>0.000000e+00</td>
    </tr>
    <tr>
      <th>10</th>
      <td>34DHOXPEGOX</td>
      <td>0.000000e+00</td>
      <td>0.000000e+00</td>
    </tr>
    <tr>
      <th>11</th>
      <td>34HPPte</td>
      <td>0.000000e+00</td>
      <td>0.000000e+00</td>
    </tr>
    <tr>
      <th>12</th>
      <td>3AIBTm</td>
      <td>NaN</td>
      <td>0.000000e+00</td>
    </tr>
    <tr>
      <th>13</th>
      <td>3AIB_Dtm</td>
      <td>0.000000e+00</td>
      <td>NaN</td>
    </tr>
    <tr>
      <th>14</th>
      <td>3AIBt</td>
      <td>NaN</td>
      <td>0.000000e+00</td>
    </tr>
    <tr>
      <th>15</th>
      <td>3AIBtm</td>
      <td>NaN</td>
      <td>0.000000e+00</td>
    </tr>
    <tr>
      <th>16</th>
      <td>3HAO</td>
      <td>0.000000e+00</td>
      <td>NaN</td>
    </tr>
    <tr>
      <th>17</th>
      <td>3HBCDm</td>
      <td>0.000000e+00</td>
      <td>NaN</td>
    </tr>
    <tr>
      <th>18</th>
      <td>3HBCOAHLm</td>
      <td>0.000000e+00</td>
      <td>0.000000e+00</td>
    </tr>
    <tr>
      <th>19</th>
      <td>3HBCOARc</td>
      <td>0.000000e+00</td>
      <td>0.000000e+00</td>
    </tr>
    <tr>
      <th>20</th>
      <td>3HCO3_NAt</td>
      <td>0.000000e+00</td>
      <td>0.000000e+00</td>
    </tr>
    <tr>
      <th>21</th>
      <td>3MLDAt</td>
      <td>0.000000e+00</td>
      <td>0.000000e+00</td>
    </tr>
    <tr>
      <th>22</th>
      <td>3MOBt2im</td>
      <td>0.000000e+00</td>
      <td>0.000000e+00</td>
    </tr>
    <tr>
      <th>23</th>
      <td>3MOBte</td>
      <td>0.000000e+00</td>
      <td>0.000000e+00</td>
    </tr>
    <tr>
      <th>24</th>
      <td>3MOPt2im</td>
      <td>0.000000e+00</td>
      <td>0.000000e+00</td>
    </tr>
    <tr>
      <th>25</th>
      <td>3MOPte</td>
      <td>0.000000e+00</td>
      <td>0.000000e+00</td>
    </tr>
    <tr>
      <th>26</th>
      <td>4ABUTtcn</td>
      <td>0.000000e+00</td>
      <td>0.000000e+00</td>
    </tr>
    <tr>
      <th>27</th>
      <td>4HOXPACDOX_NADP_</td>
      <td>0.000000e+00</td>
      <td>0.000000e+00</td>
    </tr>
    <tr>
      <th>28</th>
      <td>4MOPt2im</td>
      <td>0.000000e+00</td>
      <td>0.000000e+00</td>
    </tr>
    <tr>
      <th>29</th>
      <td>4MOPte</td>
      <td>0.000000e+00</td>
      <td>0.000000e+00</td>
    </tr>
    <tr>
      <th>...</th>
      <td>...</td>
      <td>...</td>
      <td>...</td>
    </tr>
    <tr>
      <th>2767</th>
      <td>r2433</td>
      <td>0.000000e+00</td>
      <td>2.725616e-03</td>
    </tr>
    <tr>
      <th>2768</th>
      <td>r2438</td>
      <td>0.000000e+00</td>
      <td>0.000000e+00</td>
    </tr>
    <tr>
      <th>2769</th>
      <td>r2444</td>
      <td>NaN</td>
      <td>1.000000e-02</td>
    </tr>
    <tr>
      <th>2770</th>
      <td>r2447</td>
      <td>0.000000e+00</td>
      <td>1.000000e-02</td>
    </tr>
    <tr>
      <th>2771</th>
      <td>r2449</td>
      <td>0.000000e+00</td>
      <td>1.000000e-02</td>
    </tr>
    <tr>
      <th>2772</th>
      <td>r2471</td>
      <td>0.000000e+00</td>
      <td>0.000000e+00</td>
    </tr>
    <tr>
      <th>2773</th>
      <td>r2472</td>
      <td>0.000000e+00</td>
      <td>0.000000e+00</td>
    </tr>
    <tr>
      <th>2774</th>
      <td>r2473</td>
      <td>0.000000e+00</td>
      <td>0.000000e+00</td>
    </tr>
    <tr>
      <th>2775</th>
      <td>r2493</td>
      <td>0.000000e+00</td>
      <td>NaN</td>
    </tr>
    <tr>
      <th>2776</th>
      <td>r2495</td>
      <td>NaN</td>
      <td>0.000000e+00</td>
    </tr>
    <tr>
      <th>2777</th>
      <td>r2502</td>
      <td>0.000000e+00</td>
      <td>0.000000e+00</td>
    </tr>
    <tr>
      <th>2778</th>
      <td>r2505</td>
      <td>0.000000e+00</td>
      <td>0.000000e+00</td>
    </tr>
    <tr>
      <th>2779</th>
      <td>r2506</td>
      <td>0.000000e+00</td>
      <td>NaN</td>
    </tr>
    <tr>
      <th>2780</th>
      <td>r2508</td>
      <td>0.000000e+00</td>
      <td>NaN</td>
    </tr>
    <tr>
      <th>2781</th>
      <td>r2509</td>
      <td>0.000000e+00</td>
      <td>0.000000e+00</td>
    </tr>
    <tr>
      <th>2782</th>
      <td>r2511</td>
      <td>0.000000e+00</td>
      <td>0.000000e+00</td>
    </tr>
    <tr>
      <th>2783</th>
      <td>r2514</td>
      <td>-1.509674e-03</td>
      <td>-2.121704e-03</td>
    </tr>
    <tr>
      <th>2784</th>
      <td>r2516</td>
      <td>0.000000e+00</td>
      <td>0.000000e+00</td>
    </tr>
    <tr>
      <th>2785</th>
      <td>r2517</td>
      <td>0.000000e+00</td>
      <td>0.000000e+00</td>
    </tr>
    <tr>
      <th>2786</th>
      <td>r2518</td>
      <td>0.000000e+00</td>
      <td>0.000000e+00</td>
    </tr>
    <tr>
      <th>2787</th>
      <td>r2520</td>
      <td>NaN</td>
      <td>0.000000e+00</td>
    </tr>
    <tr>
      <th>2788</th>
      <td>r2521</td>
      <td>0.000000e+00</td>
      <td>0.000000e+00</td>
    </tr>
    <tr>
      <th>2789</th>
      <td>r2525</td>
      <td>0.000000e+00</td>
      <td>0.000000e+00</td>
    </tr>
    <tr>
      <th>2790</th>
      <td>r2526</td>
      <td>0.000000e+00</td>
      <td>4.071467e-02</td>
    </tr>
    <tr>
      <th>2791</th>
      <td>r2532</td>
      <td>0.000000e+00</td>
      <td>0.000000e+00</td>
    </tr>
    <tr>
      <th>2792</th>
      <td>r2534</td>
      <td>0.000000e+00</td>
      <td>-1.402281e-02</td>
    </tr>
    <tr>
      <th>2793</th>
      <td>r2537</td>
      <td>0.000000e+00</td>
      <td>0.000000e+00</td>
    </tr>
    <tr>
      <th>2794</th>
      <td>sink_citr_LPAREN_c_RPAREN_</td>
      <td>0.000000e+00</td>
      <td>0.000000e+00</td>
    </tr>
    <tr>
      <th>2795</th>
      <td>sink_octdececoa_LPAREN_c_RPAREN_</td>
      <td>-3.469447e-18</td>
      <td>7.285839e-17</td>
    </tr>
    <tr>
      <th>2796</th>
      <td>xmpt</td>
      <td>1.000000e-02</td>
      <td>1.000000e-02</td>
    </tr>
  </tbody>
</table>
<p>2797 rows × 3 columns</p>
</div>


