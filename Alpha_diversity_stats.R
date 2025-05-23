library(vegan)
library(ggplot2)
library(GUniFrac)
library(patchwork)
library(emmeans)
library(lme4)

setwd("~/Documents/Bioinformatics_scripts/R_scripts/Kendall_optimization/")

asv=read.table("outputs/FractionationOptimization_QCfiltered.txt", header = TRUE, row.names = 1, sep = "\t")[,1:48]
tax=read.table("outputs/FractionationOptimization_QCfiltered.txt", header = TRUE, row.names = 1, sep = "\t")[,49:56]
met=read.table("Inputs/FractionationOptimization_metafile.txt", header = TRUE,row.names = 1, sep = "\t")
met$Group=paste(met$SampleType, met$NucleicAcid, sep="-")
met$Group=gsub("FRE", "Fresh", met$Group) 
met$Group=gsub("FRO", "Frozen", met$Group) 
names(asv)
# removing Z.F2.D.P2.C2 due to low number of reads 

###rarefying
cnts=t(asv[,!colnames(asv) == "Z.F2.D.P2.C2"])
min(rowSums(cnts)) # determine sample with lowest counts
asv.rar=Rarefy(cnts, min(rowSums(cnts)))$otu.tab.rff

############################################################
##################### Alpha-diversity ######################
############################################################

alpha=as.data.frame(t(estimateR(asv.rar)))
alpha$Shannon=diversity(asv.rar, index = "shannon")
alpha$SampleType=met$SampleType[match(rownames(alpha), rownames(met))]
alpha$Fraction=met$Fraction[match(rownames(alpha), rownames(met))]
alpha$NucleicAcid=met$NucleicAcid[match(rownames(alpha), rownames(met))]
alpha$Protocol=met$Protocol[match(rownames(alpha), rownames(met))]
alpha$Group=paste(alpha$NucleicAcid, alpha$SampleType, alpha$Protocol, sep = "-")
alpha$Group2=paste(alpha$NucleicAcid, alpha$SampleType, sep = "-")
#colnames(alpha)[6]="Shannon"
##################################################
##################### Stats ######################
##################################################
#Shannon
# sha_model=lm(Shannon~ Fraction +Group, data = alpha)
# anova(sha_model) 
# summary(sha_model) 
# 
# #pairwise 
# sha_pairs = emmeans(sha_model, pairwise ~ Fraction|Group, weights = "proportional", adjust="none")
# rbind(sha_pairs$contrasts, adjust="fdr")
# sha_pairs_df=as.data.frame(sha_pairs[["contrasts"]])

### only unfrationated samples
alpha.f1=subset(alpha, Fraction == "F1")
cha_model=lm(S.chao1~ NucleicAcid+SampleType+Protocol, data = alpha.f1)
anova(cha_model) 

#All samples

cha_model=lm(S.chao1~ Fraction +NucleicAcid+SampleType+Protocol, data = alpha)
anova(cha_model) 
summary(cha_model) 

#pairwise F1 vs F2
#alpha$Group2=paste(alpha$NucleicAcid, alpha$SampleType, sep = "-")
cha_model=lm(S.chao1~ Fraction +NucleicAcid+SampleType+Protocol, data = alpha)
cha_pairs = emmeans(cha_model, pairwise ~ Fraction|NucleicAcid+SampleType+Protocol, weights = "proportional", adjust="none")
rbind(cha_pairs$contrasts, adjust="fdr")
cha_pairs_df=as.data.frame(cha_pairs[["contrasts"]])

#### per group
alpha_sub=subset(alpha,  NucleicAcid == "DNA" & SampleType == "Frozen") #& !Protocol == "P1")
cha_model=lm(S.chao1~ Protocol, data = alpha_sub)
anova(cha_model) 

#overall model


## QC models

res_sha <- resid(sha_model)
plot(fitted(sha_model), res_sha)
qqnorm(res_sha)
qqline(res_sha) 
plot(density(res_sha))

res_cha <- resid(cha_model)
plot(fitted(cha_model), res_cha)
qqnorm(res_cha)
qqline(res_cha) 
plot(density(res_cha))


