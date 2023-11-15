########################
# Merge in covariates and adjust for them #
########################################### 

## Model 1: age, sex, array, pos (already in file), need: WBC props, PCs, DNAm PCs, ICV
## Model 2: above and smoking and education
## Model 3: above and handedness 

## WBC props 
wbcs=read.csv("../GrimAge_LBC36_allwaves.csv")
## subset to wave 2
wbcs=wbcs[wbcs$WAVE ==2,]
## merge in with pheno1 file 
names(wbcs)[1] <- "Basename"
pheno2=merge(wbcs[,c("Basename", "NK", "Mono", "Bcell", "Gran", "CD4T")], pheno1, by = "Basename", all.y= T)

## PCs 
pc1=read.csv("LBC1936_PCA_from_Gail.csv")
names(pc1)[1] <- "ID"
pheno2=merge(pc1[,c("ID", "C1", "C2", "C3", "C4")], pheno2, by = "ID", all.y= T)

## ICV 
icv=read.spss("../GrimAge/LBC1936_GrimAge_RH_13MAY2019.sav",to.data.frame=T)
icv1=icv[,c("lbc36no","ICV_mm3_wX", "ICV_wave_wXn")]
names(icv1)[1] <- "ID"
icv1=icv1[icv1$ICV_wave_wXn ==2,]
pheno2=merge(icv1[,c("ID", "ICV_mm3_wX")], pheno2, by = "ID", all.y= T)

## Smoking and education 
cov=read.spss("../GrimAge/LBC1936_DNA_MethylationBasedEstimatorOfTeloLength_RM_28FEB2019.sav", to.data.frame=T)
cov1=cov[,c("lbc36no", "smokcat_w1", "yrsedu_w1")]
names(cov1)[1] <- "ID"
cov1$smokcat_w1 <- as.character(cov1$smokcat_w1)
cov1[which(cov1$smokcat_w1 %in%"never smoked"),"smokcat_w1"] <- "never"
cov1[which(cov1$smokcat_w1 %in%"ex-smoker"),"smokcat_w1"] <- "ex"
cov1[which(cov1$smokcat_w1 %in%"current smoker"),"smokcat_w1"] <- "current"
pheno3=merge(cov1[,c("ID", "smokcat_w1", "yrsedu_w1")], pheno2, by = "ID", all.y= T)

## Handedness 
hand=read.spss("LBC1936_handnes_w2_RH_13JUN2022.sav",to.data.frame=T)
pheno4=merge(pheno3, hand,by.x="ID", by.y="lbc36no", all.x=T)


###########################################
# Write out phenotype and covariate files #
########################################### 

## set up loop to run lm and save out phenotype file 
var=c("Hippo_Left_mm3_w2","Hippo_Right_mm3_w2", "Mean_Hippo","Asym_Hippo", "L.GM.Vol", "R.GM.Vol", "GM.Asymmetry")
var.names=c("L.Hippo.Vol", "R.Hippo.Vol", "M.Hippo.Vol", "Hippo.Asymmetry","L.GM.Vol", "R.GM.Vol", "GM.Asymmetry")

## Write out Phenotypes 
for(i in 1:length(var)){ 
  tmp <- data.frame(FID = pheno2$Basename, 
                    IID = pheno2$Basename,
                    phen = pheno2[,var[i]])
  write.table(tmp, file=paste0("Phenotypes/", var.names[i], ".phen"), row.names=F, sep=' ')
} 

## Model1
quant_cov_m1 <- data.frame(FID = pheno2$Basename, 
                           IID = pheno2$Basename,
                           Age = pheno2$ageMRI_w2,
                           CD4T = pheno2$CD4T,
                           Mono = pheno2$Mono,
                           Bcell = pheno2$Bcell,
                           Gran = pheno2$Gran,
                           NK = pheno2$NK,
                           C1 = pheno2$C1,
                           C2 = pheno2$C2,
                           C3 = pheno2$C3, 
                           C4 = pheno2$C4,
                           ICV = pheno2$ICV_mm3_wX)
write.table(quant_cov_m1, file="quant_Model1.cov", row.names=F, sep=' ')

fact_cov_m1 <- data.frame(FID = pheno2$Basename, 
                          IID = pheno2$Basename,
                          sex = pheno2$sex.x,
                          set=pheno2$set) 
write.table(fact_cov_m1, file="fact_Model1.cov", row.names=F, sep=' ')

fact_cov_no_sex_m1 <- data.frame(FID = pheno2$Basename, 
                                 IID = pheno2$Basename,
                                 smok = pheno2$set)
write.table(fact_cov_no_sex_m1, file="fact_no_sex_Model1.cov", row.names=F, sep=' ')



## Model2 
quant_cov_m2 <- data.frame(FID = pheno3$Basename, 
                           IID = pheno3$Basename,
                           Age = pheno3$ageMRI_w2,
                           CD4T = pheno3$CD4T,
                           Mono = pheno3$Mono,
                           Bcell = pheno3$Bcell,
                           Gran = pheno3$Gran,
                           NK = pheno3$NK,
                           C1 = pheno3$C1,
                           C2 = pheno3$C2,
                           C3 = pheno3$C3, 
                           C4 = pheno3$C4,
                           edu = pheno3$yrsedu_w1,
                           ICV = pheno3$ICV_mm3_wX) 
write.table(quant_cov_m2, file="quant_Model2.cov", row.names=F, sep=' ')

fact_cov_m2 <- data.frame(FID = pheno3$Basename, 
                          IID = pheno3$Basename,
                          sex = pheno3$sex.x,
                          set=pheno3$set,
                          smok = pheno3$smokcat_w1)
write.table(fact_cov_m2, file="fact_Model2.cov", row.names=F, sep=' ', quote =T)

fact_cov_no_sex_m2 <- data.frame(FID = pheno3$Basename, 
                                 IID = pheno3$Basename,
                                 smok = pheno3$smokcat_w1)
write.table(fact_cov_no_sex_m2, file="fact_no_sex_Model2.cov", row.names=F, sep=' ')



## Model3
quant_cov_m3 <- data.frame(FID = pheno4$Basename, 
                           IID = pheno4$Basename,
                           Age = pheno4$ageMRI_w2,
                           CD4T = pheno4$CD4T,
                           Mono = pheno4$Mono,
                           Bcell = pheno4$Bcell,
                           Gran = pheno4$Gran,
                           NK = pheno4$NK,
                           C1 = pheno4$C1,
                           C2 = pheno4$C2,
                           C3 = pheno4$C3, 
                           C4 = pheno4$C4,
                           edu = pheno4$yrsedu_w1,
                           ICV = pheno4$ICV_mm3_wX) 
write.table(quant_cov_m3, file="quant_Model3.cov", row.names=F, sep=' ')

fact_cov_m3 <- data.frame(FID = pheno4$Basename, 
                          IID = pheno4$Basename,
                          sex = pheno4$sex.x,
                          set = pheno4$set,
                          smok = pheno4$smokcat_w1,
                          hand = pheno4$handnes_w2)
write.table(fact_cov_m3, file="fact_Model3.cov", row.names=F, sep=' ')

fact_cov_no_sex_m3 <- data.frame(FID = pheno4$Basename, 
                                 IID = pheno4$Basename,
                                 set = pheno4$set,
                                 smok = pheno4$smokcat_w1,
                                 hand = pheno4$handnes_w2)
write.table(fact_cov_no_sex_m3, file="fact_no_sex_Model3.cov", row.names=F, sep=' ')


## Make files to indicate male or female subsets in analyses 
#males 
males=as.data.frame(pheno2[pheno2$sex.x %in% "Male","Basename"])
names(males)[1] <- "FID"
males$IID <- males$FID
write.table(males, "males_keep.txt",row.names =F, col.names=F,  sep=' ', quote = T)
#females 
females=as.data.frame(pheno2[pheno2$sex.x %in% "Female","Basename"])
names(females)[1] <- "FID"
females$IID <- females$FID
write.table(females, "females_keep.txt",row.names =F, col.names=F,  sep=' ', quote = T)



###########################################
# Make DNAm txt file for OSCA  ############
########################################### 

## Subset DNAm to those with complete phenotypic data
dat1=dat[,which(colnames(dat) %in% tar1$Basename)]
## ensure IDs match 
ids = tar1$Basename 
dat1=dat1[,match(ids,colnames(dat1))]
table(colnames(dat1)==tar1$Basename)
## Save out file 
osca_dat <- data.frame(IID=colnames(dat1), FID=colnames(dat1))
osca_dat <- cbind(osca_dat, t(dat1))
write.table(osca_dat, file="myprofile.txt", row.names=F, sep=' ')



###########################################
# CMD Line: Make binary DNAm file  ########
########################################### 
cd /Cluster_Filespace/Marioni_Group/Rob/Hippocampal_EWAS_LBC/ 
  osca_Linux --efile myprofile.txt --methylation-beta --make-bod --out output_filename


###########################################
# R: Make .opi annotation file  ###########
########################################### 

anno=readRDS("/Cluster_Filespace/Marioni_Group/Daniel/450k_annotation.rds")
opi <- anno[colnames(osca_dat)[3:ncol(osca_dat)],c("chr","Name", "pos","UCSC_RefGene_Name", "strand")]
opi$chr <- gsub("chr", "", opi$chr)
opi$chr <- gsub("X", "23", opi$chr)
opi$chr <- gsub("Y", "24", opi$chr)
opi$chr <- as.numeric(opi$chr)
opi$UCSC_RefGene_Name	 <- as.factor(opi$UCSC_RefGene_Name	)
opi$strand <- as.factor(opi$strand)
opi[which(opi$UCSC_RefGene_Name	==""), "UCSC_RefGene_Name"] <- NA
write.table(opi, file="output_filename.opi",
            col.names=F,
            row.names=F,
            quote=F, sep='\t')

