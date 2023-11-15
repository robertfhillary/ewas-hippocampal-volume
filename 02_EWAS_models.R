###########################################
# CMD Line: EWAS models  ##################
########################################### 

screen -S ewas 

### Model 1 ####

## Males and females 
cd /Cluster_Filespace/Marioni_Group/Rob/Hippocampal_EWAS_LBC/Phenotypes/ 
  for i in *.phen
do  
# extracts phenotype name 
A=$(echo "$i"|rev|cut -d. -f1  --complement|rev)
osca_Linux \
--linear \
--befile ../output_filename \
--pheno $i \
--qcovar ../quant_Model1.cov \
--covar  ../fact_Model1.cov \
--fast-linear \
--out ../Outputs/LBC1936_450K_EA_${A}_M1_Both_WholeBlood_14062022 \
--methylation-beta 
done 

## Males 
cd /Cluster_Filespace/Marioni_Group/Rob/Hippocampal_EWAS_LBC/Phenotypes/ 
  for i in *.phen
do  
# extracts phenotype name 
A=$(echo "$i"|rev|cut -d. -f1  --complement|rev)
osca_Linux \
--linear \
--befile ../output_filename \
--pheno $i \
--keep ../males_keep.txt \
--qcovar ../quant_Model1.cov \
--covar  ../fact_no_sex_Model1.cov \
--fast-linear \
--out ../Outputs/LBC1936_450K_EA_${A}_M1_Males_WholeBlood_14062022 \
--methylation-beta 
done 

## Females  
cd /Cluster_Filespace/Marioni_Group/Rob/Hippocampal_EWAS_LBC/Phenotypes/ 
  for i in *.phen
do  
# extracts phenotype name 
A=$(echo "$i"|rev|cut -d. -f1  --complement|rev)
osca_Linux \
--linear \
--befile ../output_filename \
--pheno $i \
--keep ../females_keep.txt \
--qcovar ../quant_Model1.cov \
--covar  ../fact_no_sex_Model1.cov \
--fast-linear \
--out ../Outputs/LBC1936_450K_EA_${A}_M1_Females_WholeBlood_14062022 \
--methylation-beta 
done 



### Model 2 ####
## Males and females 
cd /Cluster_Filespace/Marioni_Group/Rob/Hippocampal_EWAS_LBC/Phenotypes/ 
  for i in *.phen
do  
# extracts phenotype name 
A=$(echo "$i"|rev|cut -d. -f1  --complement|rev)
osca_Linux \
--linear \
--befile ../output_filename \
--pheno $i \
--qcovar ../quant_Model2.cov \
--covar  ../fact_Model2.cov \
--fast-linear \
--out ../Outputs/LBC1936_450K_EA_${A}_M2_Both_WholeBlood_14062022 \
--methylation-beta 
done 

## Males 
cd /Cluster_Filespace/Marioni_Group/Rob/Hippocampal_EWAS_LBC/Phenotypes/ 
  for i in *.phen
do  
# extracts phenotype name 
A=$(echo "$i"|rev|cut -d. -f1  --complement|rev)
osca_Linux \
--linear \
--befile ../output_filename \
--pheno $i \
--keep ../males_keep.txt \
--qcovar ../quant_Model2.cov \
--covar  ../fact_no_sex_Model2.cov \
--fast-linear \
--out ../Outputs/LBC1936_450K_EA_${A}_M2_Males_WholeBlood_14062022 \
--methylation-beta 
done 

## Females 
cd /Cluster_Filespace/Marioni_Group/Rob/Hippocampal_EWAS_LBC/Phenotypes/ 
  for i in *.phen
do  
# extracts phenotype name 
A=$(echo "$i"|rev|cut -d. -f1  --complement|rev)
osca_Linux \
--linear \
--befile ../output_filename \
--pheno $i \
--keep ../females_keep.txt \
--qcovar ../quant_Model2.cov \
--covar  ../fact_no_sex_Model2.cov \
--fast-linear \
--out ../Outputs/LBC1936_450K_EA_${A}_M2_Females_WholeBlood_14062022 \
--methylation-beta 
done 



### Model 3 ####
## Males and females 
cd /Cluster_Filespace/Marioni_Group/Rob/Hippocampal_EWAS_LBC/Phenotypes/ 
  for i in *.phen
do  
# extracts phenotype name 
A=$(echo "$i"|rev|cut -d. -f1  --complement|rev)
osca_Linux \
--linear \
--befile ../output_filename \
--pheno $i \
--qcovar ../quant_Model3.cov \
--covar  ../fact_Model3.cov \
--fast-linear \
--out ../Outputs/LBC1936_450K_EA_${A}_M3_Both_WholeBlood_14062022 \
--methylation-beta 
done 

## Males 
cd /Cluster_Filespace/Marioni_Group/Rob/Hippocampal_EWAS_LBC/Phenotypes/ 
  for i in *.phen
do  
# extracts phenotype name 
A=$(echo "$i"|rev|cut -d. -f1  --complement|rev)
osca_Linux \
--linear \
--befile ../output_filename \
--pheno $i \
--keep ../males_keep.txt \
--qcovar ../quant_Model3.cov \
--covar  ../fact_no_sex_Model3.cov \
--fast-linear \
--out ../Outputs/LBC1936_450K_EA_${A}_M3_Males_WholeBlood_14062022 \
--methylation-beta 
done 

## Females 
cd /Cluster_Filespace/Marioni_Group/Rob/Hippocampal_EWAS_LBC/Phenotypes/ 
  for i in *.phen
do  
# extracts phenotype name 
A=$(echo "$i"|rev|cut -d. -f1  --complement|rev)
osca_Linux \
--linear \
--befile ../output_filename \
--pheno $i \
--keep ../females_keep.txt \
--qcovar ../quant_Model3.cov \
--covar  ../fact_no_sex_Model3.cov \
--fast-linear \
--out ../Outputs/LBC1936_450K_EA_${A}_M3_Females_WholeBlood_14062022 \
--methylation-beta 
done 

