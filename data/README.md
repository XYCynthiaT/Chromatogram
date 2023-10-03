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
    istd = "data/istd(all).rds". Selected by the minimun p-values from the the comparison between blank samples and pooled samples.  
    Concentration calculated by the istd with most similar chain length.  
    Cholesterol issue.  
- *untargeted_singletp_fin3_0.05.rds*:  
    Same as *untargeted_singletp_fin3_0.01.rds*, but pval cutoff = 0.05.
- untargeted_singletp_fin_cid.rds:  
    - InChI Key to CID (623,2). NA excluded.
<br>   

## istd.R
- *istd.rds*:
    One istd per lipid class. Selected by lowest CV.
- *istd(pval).rds*:
    One istd per lipid class.  
    Criteria: 1. The mean of blank samples \< pooled smaples (p\<0.05); 2. Lowest CV in pooled samples.
- *istd(all).rds*:
    Multiple istd per lipid class.
    Criteria: 1. The mean of blank samples \< pooled samples (p\<0.05).
