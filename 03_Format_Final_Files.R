
###########################################
# R: Format final results files  ##########
########################################### 

## Set working directory 
setwd("/Cluster_Filespace/Marioni_Group/Rob/Hippocampal_EWAS_LBC/")
library(data.table)

## Get mean beta values for CpGs across the 543 individuals used in Model1 + Model2 
meth = load("/Cluster_Filespace/Marioni_Group/LBC/LBC_methylation/LBC_betas_3489_bloodonly.rds")
meth <- dat 

## Read in target file to subset to LBC36 - Wave 2 only (when brain scans were taken)
tar=readRDS("/Cluster_Filespace/Marioni_Group/LBC/LBC_methylation/targets_3489_bloodonly.rds")
tar1=tar[tar$cohort=="LBC36",]
tar1=tar1[tar1$WAVE==2,]

## Subset methylation data 
meth1=meth[,which(colnames(meth) %in% tar1$Basename)]

## Read in list of males and females - will use this in loop step to find correct means and sds of cpgs in analysis 
males=read.table("males_keep.txt",header=F)
females=read.table("females_keep.txt",header=F)


## Extract files to tidy up 
loop1=list.files("Outputs/", ".linear")

## Loop through files within a given folder 
for(file in loop1){
  
  ## Read in file 
  tmp = as.data.frame(fread(paste0("Outputs/", file))) 
  
  ## Extract filename 
  filename = gsub(".linear", "", file)
  
  ## Extract t-value 
  tmp$t_val = tmp$b/tmp$se
  
  ## Extract required columns 
  tmp1 <- tmp[,c("Probe","b","se","t_val","p","NMISS")]
  
  ## Rearrange columns as requested 
  names(tmp1) <- c("cpg","beta","se", "t_val", "p_val","N")
  
  ## extract pheno name 
  pheno.name=strsplit(file, "_")[[1]][4]
  
  ## extract model 
  model.name=strsplit(file, "_")[[1]][5]
  ## reformat model name so that we can pull it from covariate information
  model.name=gsub("M", "Model", model.name)
  
  ## extract whether its sex-specific analyses or not 
  type.name=strsplit(file, "_")[[1]][6]
  
  ## Subset methylation dataframe to appropriate people 
  ## Read in model 1 IDs
  m1=read.table(paste0("quant_", model.name, ".cov"),header=T)
  m1.c=read.table(paste0("fact_", model.name, ".cov"), header = T)
  m1.c=m1.c[complete.cases(m1.c),]
  m1=m1[complete.cases(m1),]
  m1= merge(m1,m1.c[,c("IID", "sex")],by="IID")
  m1.p=read.table(paste0("Phenotypes/",pheno.name, ".phen"),header=T)
  m1.p=m1.p[complete.cases(m1.p),]
  m1=merge(m1,m1.p[,c("IID","phen")],by="IID")
  meth.tmp=meth1[,which(colnames(meth1) %in% m1$IID)]
  
  if(type.name %in% "Both"){ 
    
    means=apply(meth.tmp, 1, function(x) mean(x,na.rm = T))
    #get sds
    sds=apply(meth.tmp, 1, function(x) sd(x,na.rm = T))
    #combine information
    cpgs=as.data.frame(cbind(means,sds))
    cpgs$cpg <- row.names(cpgs)
    names(cpgs) <- c("mean","sd","cpg") 
    ## Merge in CpG mean and SD 
    tmp1=merge(tmp1, cpgs, by = "cpg")
    print(tmp1$N[1] == ncol(meth.tmp))
    print(model.name)
    print(type.name)
    
  } else if (type.name %in% "Males") { 
    
    meth.tmp <- meth.tmp[,which(colnames(meth.tmp) %in% males$V1)]
    means=apply(meth.tmp, 1, function(x) mean(x,na.rm = T))
    #get sds
    sds=apply(meth.tmp, 1, function(x) sd(x,na.rm = T))
    #combine information
    cpgs=as.data.frame(cbind(means,sds))
    cpgs$cpg <- row.names(cpgs)
    names(cpgs) <- c("mean","sd","cpg") 
    ## Merge in CpG mean and SD 
    tmp1=merge(tmp1, cpgs, by = "cpg")
    print(tmp1$N[1] == ncol(meth.tmp))
    print(model.name)
    print(type.name)
    
  } else { 
    meth.tmp <- meth.tmp[,which(colnames(meth.tmp) %in% females$V1)]
    means=apply(meth.tmp, 1, function(x) mean(x,na.rm = T))
    #get sds
    sds=apply(meth.tmp, 1, function(x) sd(x,na.rm = T))
    #combine information
    cpgs=as.data.frame(cbind(means,sds))
    cpgs$cpg <- row.names(cpgs)
    names(cpgs) <- c("mean","sd","cpg") 
    ## Merge in CpG mean and SD 
    tmp1=merge(tmp1, cpgs, by = "cpg")
    print(tmp1$N[1] == ncol(meth.tmp))
    print(model.name)
    print(type.name)
    
    
  }
  
  ## Save out file as requested 
  fwrite(tmp1, paste0("Final_Outputs/",filename,".txt"), sep="\t",row.names = F)
} 



###########################################
# CMD Line: Gzip final files  #############
########################################### 

screen -S gzip 

cd /Cluster_Filespace/Marioni_Group/Rob/Hippocampal_EWAS_LBC/Final_Outputs/
  gzip *
  
  
  ## Subset methylation dataframe to appropriate people 
  ## Read in model 1 IDs
  m1=read.table(paste0("quant_Model3.cov"),header=T)
m1.c=read.table(paste0("fact_Model3.cov"), header = T)
m1.c=m1.c[complete.cases(m1.c),]
m1=m1[complete.cases(m1),]
m1= merge(m1,m1.c,by="IID")
m1.p=read.table(paste0("Phenotypes/L.Hippo.Vol.phen"),header=T)
m1.p=m1.p[complete.cases(m1.p),]
m1=merge(m1,m1.p[,c("IID","phen")],by="IID")

