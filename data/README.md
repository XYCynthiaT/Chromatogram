# Change log

## Untargeted_singletp_fin.R
- *untargeted_singletp_fin1_0.05.rds*:  
    istd = istd selected by CV.  
    Features select by the pval comparing blank samples to pooled samples, p<0.05.  
- *untargeted_singletp_fin1_0.01.rds*:  
    istd = istd selected by CV.  
    Features select by the pval comparing blank samples to pooled samples, p<0.01.  
- *untargeted_singletp_fin2_0.01.rds*:  
    istd = istd selected by the minimun p-values from the the comparison between blank samples and pooled samples. Then selected by minimun CV.  
    Features select by the pval comparing blank samples to pooled samples, p<0.01.   
- *untargeted_singletp_fin2_0.05.rds*:  
    istd = istd selected by the minimun p-values from the the comparison between blank samples and pooled samples. Then selected by minimun CV.  
    Features select by the pval comparing blank samples to pooled samples, p<0.05.  
- *untargeted_singletp_fin3_0.01.rds\**:  
    quant: nfeature = 347; peak: nfeature = 461  
    Features were selected by the pval comparing blank samples to pooled samples, p<0.01.  
    cid = "data/untargeted_singletp_fin_cid.rds".   
    Fatures with InChI Key that couldn't be found in PubChem (w/o CID) were excluded.  
    filter_by_cv(): aimed at compunds with identical InChI Key. (some issues, see notes).  
    Assign lipid class from LIPID MAPS and local table.  
    istd = "data/istd(all).rds". Not selected  
    Concentration calculated by the istd with most similar unsaturated bond, then chain length.  
    Cholesterol issue.  
- *untargeted_singletp_fin3_0.05.rds*:  
    Same as *untargeted_singletp_fin3_0.01.rds*, but pval cutoff = 0.05.
- *untargeted_singletp_fin_cid.rds*:  
    - InChI Key to CID (623,2). NA excluded.
- *untargeted_singletp_fin4_all.rds*:
    - istd = istd(all).rds
    - Didn't filter by blank
- *untargeted_singletp_fin4_all(single_istd).rds*:
    - istd = istd(abundance).rds
    - Didn't filter by blank
<br>   

## Targeted_singletp_fin.R
- targeted_singletp_fin3_all.rds:  
    Fatures with InChI Key that couldn't be found in PubChem (w/o CID) were excluded.
    cid = "data/targeted_singletp_fin_cid.rds".   
    filter_by_cv(): aimed at compunds with identical InChI Key. (some issues, see notes).  
    Assign lipid class from LIPID MAPS and local table.    
    Cholesterol added, CE, DG fixed.  
    EtherX exluded.
    nfeatures = 468

## istd.R
- *istd.rds*:
    One istd per lipid class. Selected by lowest CV.
- *istd(pval).rds*:
    One istd per lipid class.  
    Criteria: 1. The mean of blank samples \< pooled smaples (p\<0.05); 2. Lowest CV in pooled samples.
- *istd(all).rds*:
    Multiple istd per lipid class.
    Criteria: no exclusion.
- *istd(abundant).rds*
    the most abundant istd for each class.
<br>

## istd2.R
istd peak heights from mx733759_Agus_Lipids_Single-point quant_Submit_CS_08-14-23 (1).xlsx. It's different from the peak heights in mx 733759_Agus_lipidomics_human plasma lipoproteins_07-2023 submit.xlsx
- istd2(all).rds: istd filtered by pval of blank vs. pool and maximum sample average peak height.
