####################################################################
####################################################################
####################################################################
###                                                              ###
###                  FIGURE 1 of Jankowski et al 2022            ###
###                         Daphnia-Bioassay                     ###
###                                                              ###
###                                                              ###
####################################################################
####################################################################
####################################################################

library(drc)
library(dplyr)

#Load Data
dat2=read.csv("short ppf data for EC10 analysis.csv"); attach(dat2); View(dat2)

#Check data structure
str(dat2) 
#Convert data from integers to numeric
dat2$Males=as.numeric(dat2$Males)
dat2$Females=as.numeric(dat2$Females)
dat2$Total=as.numeric(dat2$Total)
dat2$Conc_nom=as.integer(dat2$Conc_nom)

#Selected nominal dose-response model and selected percentiles for effective dose (ED) levels
mod2=drm(Males/Total~Conc_nom,weights=Total,data=dat1,fct=LL.2(),type="binomial") 
modelFit(mod2)
summary(mod2)
ED(mod2,c(5,10,20,25,50))
plot(mod2)

#Plotting Figure 1 for publication
tiff("Figure 1.tiff",width=2000,height=2000,res = 300)
plot(mod2,type="confidence",xlim = c(.01,1000), xlab="Pyriproxyfen (ng/L)",
     ylab="Fraction male neonates",
     confidence.level = 0.95)
points(x=jitter(Conc_nom,.5),y=jitter(Males/Total,5),pch=21)
xtick=seq(0,500, by=50)
axis(side=1,at=xtick,labels=FALSE)
dev.off()

####################################################################
####################################################################
####################################################################
###                                                              ###
###         Gene Expression Analysis in Jankowski et al 2022     ###
###                         Daphnia-Bioassay                     ###
###                                                              ###
###                                                              ###
####################################################################
####################################################################
####################################################################

#source("https://bioconductor.org/biocLite.R")
#biocLite("DESeq")
#biocLite("limma")
library(DESeq2)
library(glmnet)
library(ComplexHeatmap)
library(ggfortify)
library(limma)

rm(list=ls())

set.seed(1001)

setwd("~/MPCA_dat/")

pull_cor_block = function(data_set, gene_name){
  return(cor(data_set[,gene_name], data_set))
}

##Stress code is only for informational purposes. These analyses were not the focus of the paper.
load_stress = function(){
  ##Load Stress Datasets
  #Read in Stress experiments
  stress_exp = read.table("StressExpression.csv",header=TRUE, sep=",")
  rownames(stress_exp) = stress_exp[,1]
  stress_exp = t(stress_exp[,-1])
  stress_exp = data.frame(sample_name = rownames(stress_exp), stress_exp)
  stress_ids = read.table("StressIds1.csv",header=TRUE, sep=",")[,c(1,3,5)]
  stress_dat = merge(stress_ids,stress_exp, by.x="sample.name",by.y="sample_name") # was joint_dat
  rownames(stress_dat) = paste(stress_dat[,3], stress_dat[,2],sep="")
  return(stress_dat)
}

load_wetland = function(){
  mpca_dat = read.table(file = "kallisto.counts", header=TRUE)
  row.names(mpca_dat) = mpca_dat[,1]
  w_sample_info = read.csv("Water_sample_ref.csv") # Select water Samples
  wetland_sample_info = read.csv("Wetland information_MJ_092817.csv")
  wetland_sample_info[,2] = paste("X", wetland_sample_info[,2], sep="")
  wetland_dat = mpca_dat[,-which(colnames(mpca_dat) %in% as.character(w_sample_info[,2]))][,-1]
  wetland_dat = data.frame(Sample_ID = colnames(wetland_dat), t(wetland_dat))
  wetland_all = merge(wetland_sample_info, wetland_dat, by.x="Sample_Num", by.y="Sample_ID")[,-1]
  return(wetland_all)
}

load_stormwater = function(){
  mpca_dat = read.table(file = "kallisto.counts", header=TRUE)
  row.names(mpca_dat) = mpca_dat[,1]
  w_sample_info = read.csv("Water_sample_ref.csv")
  wetland_sample_info = read.csv("wetland_info.csv")
  wetland_sample_info[,2] = paste("X", wetland_sample_info[,2], sep="")
  water_dat = mpca_dat[,as.character(w_sample_info[,2])]
  water_dat = as.data.frame(t(water_dat))
  water_dat[,"Sample_ID"] = rownames(water_dat)
  water_all = merge(w_sample_info,water_dat, by="Sample_ID")
  return(water_all)
}

scale_type = "sum"
#Load and scale data
stress_dat = load_stress()
stress_scale = stress_dat
stress_bad=names(which(apply(stress_scale[,-1:-3],2,var)==0))

water_all = load_stormwater()
water_scale = water_all
water_bad=names(which(apply(water_scale[,-1:-5],2,var)==0))

wetland_all = load_wetland()
wetland_all[,3] = as.character(wetland_all[,3])
wetland_all[,4] = as.character(wetland_all[,4])
wetland_all[which(is.na(wetland_all[,4])),4] = "Control"
wetland_all[,5] = as.numeric(wetland_all[,5])
wetland_all[which(is.na(wetland_all[,5])),5] = 0
wetland_scale = wetland_all
wetland_bad=names(which(apply(wetland_scale[,-1:-5],2,var)==0))

#Remove non-varying genes in any dataset from all datasets, also only keeps genes occuring in all datasets
all_bad = union(union(stress_bad,water_bad),wetland_bad)
good_list = setdiff(intersect(intersect(names(wetland_scale),names(stress_scale)),names(water_scale)),all_bad)

stress_scale = cbind(stress_scale[,1:3], stress_scale[,good_list])
water_scale = cbind(water_scale[,1:5], water_scale[,good_list])
wetland_scale = cbind(wetland_scale[,1:5], wetland_scale[,good_list])

#Rescale each dataset
if(scale_type == "var"){
  stress_scale[,-1:-3] = t(scale(t(stress_scale[,-1:-3])))
}else if(scale_type == "sum"){
  stress_scale[,-1:-3] = stress_scale[,-1:-3]/rowSums(stress_scale[,-1:-3])
}

if(scale_type == "var"){
  water_scale[,-1:-5] = t(scale(t(water_scale[,-1:-5])))
}else if(scale_type == "sum"){
  water_scale[,-1:-5] = water_scale[,-1:-5]/rowSums(water_scale[,-1:-5])
}
rownames(water_scale) = paste(water_scale[,"Location"], water_scale[,"Season"],water_scale[,"Replicate"],sep="-")

if(scale_type == "var"){
  wetland_scale[,-1:-5] = t(scale(t(wetland_scale[,-1:-5])))
}else if(scale_type == "sum"){
  wetland_scale[,-1:-5] = wetland_scale[,-1:-5]/rowSums(wetland_scale[,-1:-5])
}
rownames(wetland_scale) = wetland_scale[,"ShortName"]


merge_scale = rbind(stress_scale[,-1:-3],wetland_scale[,-1:-5],water_scale[,-1:-5])
water_scale_drop = water_scale[-which(water_scale[,2] == "P Cr Lake" & water_scale[,3] == "Spring" & water_scale[,4] == 3),]
merge_scale_drop = rbind(stress_scale[,-1:-3],wetland_scale[,-1:-5],water_scale_drop[,-1:-5])
merge_scale_labels = data.frame(data_source=c(rep("stress", nrow(stress_scale)),rep("wetland", nrow(wetland_scale)),rep("stormwater", nrow(water_scale_drop))))


#rep_mark = c(stress_scale[,3], wetland_scale[,2], paste(water_scale[,2],water_scale[,3]))
#design <- model.matrix(~rep_mark)
batch <- c(sapply(strsplit(as.character(stress_scale[,2])," "), function(x) return(x[1])), rep("Wetland",nrow(wetland_scale)), rep("Water", nrow(water_scale)))
merge_scale <- t(removeBatchEffect(t(as.matrix(merge_scale)), batch = batch))

stress_scale[,-1:-3] = merge_scale[which(batch=="I" | batch=="X" | batch=="XI"),]
wetland_scale[,-1:-5] = merge_scale[which(batch=="Wetland"),]
water_scale[,-1:-5] = merge_scale[which(batch=="Water"),]

#rep_mark = c(stress_scale[,3], wetland_scale[,2], paste(water_scale_drop[,2],water_scale_drop[,3]))
#design <- model.matrix(~rep_mark)
batch <- c(sapply(strsplit(as.character(stress_scale[,2])," "), function(x) return(x[1])), rep("Wetland",nrow(wetland_scale)), rep("Water", nrow(water_scale_drop)))
merge_scale_drop <- t(removeBatchEffect(t(as.matrix(merge_scale_drop)), batch = batch))

#generate correlation tables
heat_margins=c(10,10)

corr_stress = cor(t(stress_scale[,-1:-3])) 
pdf("Corr for Stress samples.pdf", width=8, height=8)
heatmap(corr_stress,symm = TRUE, margins=heat_margins)
dev.off()

corr_water = cor(t(water_scale[,-1:-5])) 
pdf("Corr for Stormwater samples.pdf", width=16, height=16)
heatmap(corr_water,symm = TRUE, margins=heat_margins)
dev.off()

corr_wetland = cor(t(wetland_scale[,-1:-5])) 
pdf("Corr for Wetland samples.pdf", width=8, height=8)
heatmap(corr_wetland,symm = TRUE, margins=heat_margins)
dev.off()

corr_merge = cor(t(merge_scale)) 
pdf("Corr for Merged samples.pdf", width=24, height=24)
heatmap(corr_merge,symm = TRUE, margins=heat_margins)
dev.off()

#stress_table = joint_dat

#PCA for Stress data (Not shown in paper)
pc_stress = prcomp(stress_scale[,-1:-3],scale.=TRUE)
stress_labels = stress_scale[,1:3]
stress_labels[,2] = sapply(strsplit(as.character(stress_labels[,2]),split=" "), function(x) x[1])
plot(pc_stress)
summary(pc_stress)

pdf("PCA_Stress_Strain_Perturbation.pdf", width=8, height=8)
autoplot(pc_stress, data = stress_labels, shape = 'strain',colour='environmental.perturbation')
dev.off()

#Wetland PCA (Figure 3A)
pc_wetland = prcomp(wetland_scale[,-1:-5],scale.=TRUE)
wetland_labels = wetland_scale[,1:5]
plot(pc_wetland)
summary(pc_wetland)

pdf("PCA_Wetland_IBI.pdf", width=8, height=8)
autoplot(pc_wetland, data = wetland_labels,colour='Category')
dev.off()

#Figure 3A
pdf("PCA_Wetland_Location_PPF.pdf", width=8, height=8)
autoplot(pc_wetland, data = wetland_labels,colour='Location', shape="PPF")
dev.off()

#Stormwater PCA (Figure 4A)
pc_stormwater_outlier = prcomp(water_scale[,-1:-5],scale.=TRUE)
water_labels_outlier = water_scale[,1:5]
plot(pc_stormwater_outlier)
summary(pc_stormwater_outlier)

water_bad_scale=names(which(apply(water_scale_drop[,-1:-5],2,var)==0))
water_scale_drop = water_scale_drop[,setdiff(names(water_scale_drop), water_bad_scale)]
water_labels = water_scale_drop[,1:5]
pc_stormwater = prcomp(water_scale_drop[,-1:-5],scale.=TRUE)
water_labels = water_scale_drop[,1:5]
#plot(pc_stormwater)
#summary(pc_stormwater)

pdf("PCA_Stormwater_Season_Type_outlier.pdf", width=8, height=8)
autoplot(pc_stormwater_outlier, data = water_labels_outlier,shape='Season',colour='Location_type')
dev.off()

#Figure 4A
pdf("PCA_Stormwater_Season_Type.pdf", width=8, height=8)
autoplot(pc_stormwater, data = water_labels,shape='Season',colour='Location_type')
dev.off()

#PCA of all samples
pc_combined = prcomp(merge_scale_drop,scale.=TRUE)
water_labels_outlier = water_scale[,1:5]
plot(pc_combined)
summary(pc_combined)

pdf("PCA_Combined_datasource.pdf", width=8, height=8)
autoplot(pc_combined, data=merge_scale_labels, colour='data_source')
dev.off()

#Model for Wetland (Figure 3B and Table 1)
wetland_scale[which(wetland_scale[,"PPF"]=="0"),"PPF"] = "Low"
wetland_scale[which(wetland_scale[,"PPF"]=="119"),"PPF"] = "High"
wetland_scale[which(wetland_scale[,"PPF"]=="239"),"PPF"] = "High"

cv_PPF_model = cv.glmnet(as.matrix(wetland_scale[,-1:-5]),as.factor(wetland_scale[,"PPF"]),family="binomial", type.measure = "class")
cv_w_location_model = cv.glmnet(as.matrix(wetland_scale[,-1:-5]),as.factor(wetland_scale[,"Location"]),family="multinomial", type.measure = "class",nfolds=3)
cv_quality_model = cv.glmnet(as.matrix(wetland_scale[which(wetland_scale[,"Category"] != "Control"),-1:-5]),as.factor(wetland_scale[which(wetland_scale[,"Category"] != "Control"),"Category"]),family="binomial", type.measure = "class",nfolds=6)

#Model for Stormwater (Figure 4B and Table 2)
cv_s_season_model = cv.glmnet(as.matrix(water_scale_drop[,-1:-5]),as.factor(water_scale_drop[,"Season"]),family="multinomial", type.measure = "class")
cv_s_location_type_model = cv.glmnet(as.matrix(water_scale_drop[,-1:-5]),as.factor(water_scale_drop[,"Location_type"]),family="multinomial", type.measure = "class")

#Pval-generation Wetland
rep_lab_cv = function(){
  cv_rand_model = cv.glmnet(as.matrix(wetland_scale[,-1:-5]),as.factor(sample(wetland_scale[,"PPF"])),family="binomial", type.measure = "class")
  return(cv_rand_model$cvm[which(cv_rand_model$lambda == cv_rand_model$lambda.1se)])
}
cv_PPF_misclass_rand = replicate(20, rep_lab_cv())

rep_lab_cv = function(){
  cv_rand_model = cv.glmnet(as.matrix(water_scale_drop[,-1:-5]),as.factor(sample(water_scale_drop[,"Season"])),family="multinomial", type.measure = "class")
  return(cv_rand_model$cvm[which(cv_rand_model$lambda == cv_rand_model$lambda.1se)])
}
cv_w_location_misclass_rand = replicate(20, rep_lab_cv())

rep_lab_cv = function(){
  cv_rand_model = cv.glmnet(as.matrix(wetland_scale[which(wetland_scale[,"Category"] != "Control"),-1:-5]),as.factor(sample(wetland_scale[which(wetland_scale[,"Category"] != "Control"),"Category"])),family="binomial", type.measure = "class",nfolds=6)
  return(cv_rand_model$cvm[which(cv_rand_model$lambda == cv_rand_model$lambda.1se)])
}
cv_w_quality_misclass_rand = replicate(20, rep_lab_cv())

#Pval-generation Stormwater
rep_lab_cv = function(){
  cv_rand_model = cv.glmnet(as.matrix(wetland_scale[,-1:-5]),as.factor(sample(wetland_scale[,"PPF"])),family="binomial", type.measure = "class")
  return(cv_rand_model$cvm[which(cv_rand_model$lambda == cv_rand_model$lambda.1se)])
}
cv_s_season_misclass_rand = replicate(20, rep_lab_cv())

rep_lab_cv = function(){
  cv_rand_model = cv.glmnet(as.matrix(water_scale_drop[,-1:-5]),as.factor(sample(water_scale_drop[,"Location_type"])),family="multinomial", type.measure = "class")
  return(cv_rand_model$cvm[which(cv_rand_model$lambda == cv_rand_model$lambda.1se)])
}
cv_s_location_type_misclass_rand = replicate(20, rep_lab_cv())

#Figure 3B
pdf("Binomial_PPF_Wetlands.pdf",height=8,width=8)
plot(cv_PPF_model)
dev.off()

pdf("Multinomial_Location_Wetlands.pdf",height=8,width=8)
plot(cv_w_location_model)
dev.off()

pdf("Binomial_Quality_Wetlands.pdf",height=8,width=8)
plot(cv_quality_model)
dev.off()

#Figure 4B
pdf("Multinomial_Season_Stormwater.pdf",height=8,width=8)
plot(cv_s_season_model)
dev.off()

pdf("Multinomial_Location_type_Stormwater.pdf",height=8,width=8)
plot(cv_s_location_type_model)
dev.off()

##Select Lambda Value and select genes of interest: By PPF
extract_param = function(model_in,data_set,name_array){
  #Load Annot, name_array is the set of column names to pass on. 
  gene_annot_2 = read.table('cddrps-dapmaevg14.gotab2', header=FALSE, sep="\t", quote = "")
  gene_annot_4 = read.table('cddrps-dapmaevg14.gotab4', header=FALSE, sep="\t")
  gene_annot_4[,1] = substring(gene_annot_4[,1],12)
  
  lambda_loc = which(model_in$lambda == model_in$lambda.1se)
  
  if(class(model_in$glmnet.fit$beta) == "dgCMatrix"){
    gene_set = names(which(0 != model_in$glmnet.fit$beta[,lambda_loc ]))
  }else{
    gene_set = c()
    for(i in model_in$glmnet.fit$beta){
      gene_set = c(gene_set,names(which(0 != i[,lambda_loc ])))
    }
    gene_set = unique(gene_set)
  }
  
  PPF_annot_2 = gene_annot_2[which(gene_annot_2[,1] %in% gene_set),]
  PPF_annot_4 = gene_annot_4[which(gene_annot_4[,1] %in% gene_set),]
  
  selected_genes =data_set[,c(name_array,gene_set)]
  stress_genes = stress_scale[,c("sample.name", "strain", "environmental.perturbation",gene_set)]
  stress_genes = data.frame( stress_genes[,"sample.name"], "Stress Condition", "Stress Condition", stress_genes[,-1:-3,drop=FALSE])
  
  colnames(stress_genes)[1:3] = name_array[1:3]
  
  
  all_genes = rbind(selected_genes, stress_genes)
  all_genes[,-1:-3] = scale(all_genes[,-1:-3])
  return(list(gene_set, PPF_annot_2,PPF_annot_4, all_genes))
}

PPF_genes = extract_param(cv_PPF_model,wetland_scale,c("ShortName", "Location", "PPF"))
w_location_genes = extract_param(cv_w_location_model,wetland_scale,c("ShortName", "Location", "PPF"))
quality_genes = extract_param(cv_quality_model,wetland_scale,c("ShortName", "Category", "PPF"))

s_season_genes = extract_param(cv_s_season_model,water_scale_drop,c("Sample_ID", "Location_type", "Season"))
s_location_genes = extract_param(cv_s_location_type_model,water_scale_drop,c("Sample_ID", "Location_type", "Season"))



id_dat = data.frame(Sample=as.factor(PPF_genes[[4]][,"ShortName"]),Location=as.factor(PPF_genes[[4]][,"Location"]),PPF=as.factor(PPF_genes[[4]][,3]))
ha = HeatmapAnnotation(PPF_genes[[4]][,2:3],col = list(Location = c("Lab Water"="lightgreen","Kandiyohi Wetland" = "red", "Kerk Wetland" = "blue", "LeSouer Wetland"="orange", "Breen Wetland"="yellow","Woodland Wetland"="steelblue", "Douglas Wetland"="purple", "Stress Condition"="green"),
                                                       PPF = c("High"="red", "Low"="blue", "Stress Condition"="green")))

#Figure 3C
pdf("Heatmap of PPF Selected Genes.pdf", width=20, height=8)
Heatmap(as.matrix(t(PPF_genes[[4]][,-1:-3])), top_annotation=ha,clustering_distance_rows = "pearson",column_names_gp = gpar(fontsize=10))
dev.off()

id_dat = data.frame(Sample=as.factor(w_location_genes[[4]][,"ShortName"]),Location=as.factor(w_location_genes[[4]][,"Location"]),PPF=as.factor(w_location_genes[[4]][,3]))
ha = HeatmapAnnotation(w_location_genes[[4]][,2:3],col = list(Location = c("Lab Water"="lightgreen","Kandiyohi Wetland" = "red", "Kerk Wetland" = "blue", "LeSouer Wetland"="orange", "Breen Wetland"="yellow","Woodland Wetland"="steelblue", "Douglas Wetland"="purple", "Stress Condition"="green"),
                                                              PPF = c("High"="red", "Low"="blue", "Stress Condition"="green")))
pdf("Heatmap of Wetland Location Selected Genes.pdf", width=20, height=8)
Heatmap(as.matrix(t(w_location_genes[[4]][,-1:-5])), top_annotation=ha,clustering_distance_rows = "pearson",column_names_gp = gpar(fontsize=10))
dev.off()

#Quality Analysis
id_dat = data.frame(Sample=as.factor(quality_genes[[4]][,"ShortName"]),Category=as.factor(quality_genes[[4]][,"Category"]),PPF=as.factor(quality_genes[[4]][,3]))
ha = HeatmapAnnotation(quality_genes[[4]][,2:3],col = list(Category = c("Control"="yellow", "Low"="blue", "High"="red","Stress Condition"="green"),
                                                           PPF = c("High"="red", "Low"="blue", "Stress Condition"="green")))
pdf("Heatmap of Quality Selected Genes.pdf", width=20, height=8)
Heatmap(as.matrix(t(quality_genes[[4]][,-1:-5])), top_annotation=ha,clustering_distance_rows = "pearson",column_names_gp = gpar(fontsize=10))
dev.off()

mean(replicate(100,sum(wetland_scale[which(wetland_scale[,"PPF"] != "High"),2] != sample(wetland_scale[which(wetland_scale[,"PPF"] != "High"),2]))/nrow(wetland_scale)))

#Stormwater
id_dat = data.frame(Sample=as.factor(s_season_genes[[4]][,"Sample_ID"]),Location=as.factor(s_season_genes[[4]][,"Location_type"]),Season=as.factor(s_season_genes[[4]][,"Season"]))
ha = HeatmapAnnotation(s_season_genes[[4]][,2:3],col = list(Location_type = c("Out"="lightgreen","In" = "red", "Outflow" = "blue", "Control"="orange","Stress Condition"="green"),
                                                            Season = c("Spring"="yellow", "Summer1"="red", "Summer2"="blue", "Stress Condition"="green")))

#Figure 4C
pdf("Heatmap of Stormwater Season Selected Genes.pdf", width=20, height=8)
Heatmap(as.matrix(t(s_season_genes[[4]][,-1:-3])), top_annotation=ha,clustering_distance_rows = "pearson",column_names_gp = gpar(fontsize=10))
dev.off()

id_dat = data.frame(Sample=as.factor(s_location_genes[[4]][,"Sample_ID"]),Location=as.factor(s_location_genes[[4]][,"Location_type"]),Season=as.factor(s_location_genes[[4]][,"Season"]))
ha = HeatmapAnnotation(s_location_genes[[4]][,2:3],col = list(Location_type = c("Out"="lightgreen","In" = "red", "Outflow" = "blue", "Control"="orange","Stress Condition"="green"),
                                                              Season = c("Spring"="yellow", "Summer1"="red", "Summer2"="blue", "Stress Condition"="green")))
pdf("Heatmap of Stormwater Location Selected Genes.pdf", width=20, height=8)
Heatmap(as.matrix(t(s_location_genes[[4]][,-1:-3])), top_annotation=ha,clustering_distance_rows = "pearson",column_names_gp = gpar(fontsize=10))
dev.off()

print_go_annot = function(gene_annot, file_root){
  write.table(gene_annot[[2]], paste("go2_",file_root,".tab",sep=""), col.names=FALSE,row.names=FALSE)
  write.table(gene_annot[[3]], paste("go4_",file_root,".tab",sep=""), col.names=FALSE,row.names=FALSE)
}

print_go_annot(PPF_genes, "PPF_wetland")
print_go_annot(w_location_genes, "Location_wetland")
print_go_annot(quality_genes, "Quality_wetland")
print_go_annot(s_season_genes, "Season_stormwater")
print_go_annot(s_location_genes, "Location_stormwater")

####################################################################
####################################################################
####################################################################
###                                                              ###
###               Creation of DEG Sets 1-4 and Figure S2         ###
###                    For Jankowski et al 2022                  ###
###                       Daphnia-Bioassay                       ###
###                                                              ###
###                                                              ###
####################################################################
####################################################################
####################################################################

##############################
#Determine DEGs for each file#
##############################
dat1_w_names=read.csv("dat1_w_names.csv");attach(dat1_w_names);View(dat1_w_names)
dat2_w_names=read.csv("dat2_w_names.csv");attach(dat2_w_names);View(dat2_w_names)
dat3_w_names=read.csv("dat3_w_names.csv");attach(dat3_w_names);View(dat3_w_names)
dat4_w_names=read.csv("dat4_w_names.csv");attach(dat4_w_names);View(dat4_w_names)

#Subset the data by Padj < 0.05 and fold change >=2.0
DEGs1=subset(dat1_w_names,padj.Lab0vsWet0<0.05 & (log2FC.Lab0vsWet0>=2|log2FC.Lab0vsWet0<=-2));DEGs1=DEGs1[c(-12,-13)]
DEGs2=subset(dat2_w_names,padj.Lab0vs119<0.05 & (log2FC.Lab0vs119>=2.0|log2FC.Lab0vs119<=-2))
DEGs3=subset(dat3_w_names,padj.Wet0vs119<0.05 & (log2FC.Wet0vs119>=2|log2FC.Wet0vs119<=-2))
DEGs4=subset(dat4_w_names,padj.119LabvsWet<0.05 & (log2FC.119LabvsWet>=2|log2FC.119LabvsWet<=-2))

#Count number of Up and Down DEGs for each gene set
#Up
length(which(DEGs1$log2FC.Lab0vsWet0>=2))
length(which(DEGs2$log2FC.Lab0vs119>=2))
length(which(DEGs3$log2FC.Wet0vs119>=2))
length(which(DEGs4$log2FC.119LabvsWet>=2))

#Down 
length(which(DEGs1$log2FC.Lab0vsWet0<=-2))
length(which(DEGs2$log2FC.Lab0vs119<=-2))
length(which(DEGs3$log2FC.Wet0vs119<=-2))
length(which(DEGs4$log2FC.119LabvsWet<=-2))

#Change transcriptID column names relevant to each set of DEGs and then cbind those 4 columns into one file
DEGs1=DEGs1%>%rename(DEG1tID=transcriptID)
DEGs2=DEGs2%>%rename(DEG2tID=transcriptID)
DEGs3=DEGs3%>%rename(DEG3tID=transcriptID)
DEGs4=DEGs4%>%rename(DEG4tID=transcriptID)

#Creating one file with sig DEGs for each condition
##Next step was to bring only trascriptID columns from all 4 DEGs files into one VennDEG file##
#Choice is to do this in Excel, so exporting files as below#
#Exporting files to use in Excel#
write.csv(DEGs1,"DEGS1.csv")
write.csv(DEGs2,"DEGS2.csv")
write.csv(DEGs3,"DEGS3.csv")
write.csv(DEGs4,"DEGS4.csv")

#Bring 4 column file back into R#
VennDEG=read.csv("VennDEG.csv");attach(VennDEG)
#Replace empty cells with NA because columns are of different length
VennDEG=VennDEG%>%mutate_all(na_if,"")
vdat=VennDEG


####Creation of Figure S2 (Venn Diagram) for Genes across water conditions####
library(VennDiagram)
library(readxl)

genevenn=venn.diagram(x=list(vdat$DEG1ID,vdat$DEG2ID,vdat$DEG3ID,vdat$DEG4ID),
                      category.names=c("Lab.v.Wet(no PPF)","PPF.in.LW","PPF.in.Wet","Lab.v.Wet(PPF)"),
                      filename = "allgenes.tiff",output=TRUE,na="none"
                      ,imagetype = "png",height=1500,width = 1500,resolution = 300,compression = "lzw",
                      lwd=2,lty=1,fill=c(3,2,1,4),
                      cex=.6,fontface="plain",fontfamily="sans", cat.cex=0.6,
                      cat.fontface="bold",cat.default.pos="outer",cat.pos=c(-10,10,-20,0),
                      cat.dist=c(0.2,0.2,0.1,0.08),cat.fontfamily="sans")

#Display saved image                  
options(repr.plot.height=12,plot.width=12)
library(png) 
grid.newpage()
plot=readPNG("prac.png")
plot.new()
rasterImage(plot,0,1.1,0,1.1)


####################################################################
####################################################################
####################################################################
###                                                              ###
###        Figure 5 Redundancy Analysis in Jankowski et al 2022  ###
###                         Daphnia-Bioassay                     ###
###                                                              ###
###                                                              ###
####################################################################
####################################################################
####################################################################

##Code for "full RDA model" as illustrated in Figure 5A
data<- read.csv("RDA.full.csv", header=T, stringsAsFactors=F)
colnames(data)
responses<-data[,33:135]
chemsitu<-data[,4:32]
datarcs<-data.frame(chemsitu,responses)
rdarcs.all<-rda(responses~.,chemsitu)
rdarcs.0<-rda(responses~1,chemsitu)
anova(rdarcs.all)
R2rcs<-RsquareAdj(rdarcs.all)
R2rcs
screeplot(rdarcs.all)
FWrcs.full<-ordistep (rdarcs.0, scope = formula (rdarcs.all), R2scope=R2rcs, direction = 'forward', permutations=999)
summary(FWrcs.full)
Plotrcs<-ordiplot(FWrcs.full)
text(FWrcs, "species", col="red", cex=0.9)


##Code for "reduced RDA model", as illustrated in Figure 5B
data<- read.csv("RDA.reduced.csv", header=T, stringsAsFactors=F)
colnames(data)
responses<-data[,33:72]
chemsitu<-data[,4:32]
datarcs.reduced<-data.frame(chemsitu,responses)
rdarcs.reduced<-rda(responses~.,chemsitu)
rdarcs.0<-rda(responses~1,chemsitu)
anova(rdarcs.reduced)
R2rcs<-RsquareAdj(rdarcs.reduced)
R2rcs
screeplot(rdarcs.reduced)
FWrcs.reduced<-ordistep (rdarcs.0, scope = formula (rdarcs.reduced), R2scope=R2rcs, direction = 'forward', permutations=999)
screeplot(rdarcs.reduced)
summary(FWrcs.reduced)
Plotrcs<-ordiplot(FWrcs.reduced)
text(FWrcs.reduced, "species", col="red", cex=0.9)

