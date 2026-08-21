
#### Input data files
###############################################
#### spatial signatures (a list of vectors)
gsets.file = "path/to/spatial.clust_gene.sets_A.rds"

#### TCGA clinical data from cBio
dat.clin.file = "path/to/data_clinical_patient.txt"

#### TCGA expression profiles from cBio
dat.rna.file = "path/to/data_mrna_seq_v2_rsem.txt"



####### Read in and process the data 
###############################################
############ spatial signatures
gsets = readRDS(gsets.file)

###############################################
############ clinical data
datclin = read.table(dat.clin.file, sep = "\t", header = T)
dim(datclin)
rownames(datclin) = gsub("-",".",datclin$PATIENT_ID)

##### survival
## use PFS for survival data
x = datclin$PFS_STATUS
y = sapply(strsplit(x,":"),`[`,1)
datclin$surv_status = as.numeric(y)
datclin$surv_time = datclin$PFS_MONTHS

## truncate survival to 120 months
datclin$surv_status[which(datclin$surv_time > 120)] = 0
datclin$surv_time[which(datclin$surv_time > 120)] = 120


##### stage
table(datclin$AJCC_PATHOLOGIC_TUMOR_STAGE, useNA = "a")

stagetransf = vector()
stagetransf["STAGE X"] = "stX"
stagetransf[""] = "stX"
stagetransf["STAGE IIB"] = "st2"
stagetransf["STAGE IA"] = "st1"
stagetransf["STAGE IIIA"] = "st3"
stagetransf["STAGE IIA"] = "st2"
stagetransf["STAGE I"] = "st1"
stagetransf["STAGE IIIC"] = "st3"
stagetransf["STAGE IV"] = "st4"
stagetransf["STAGE IB"] = "st1"
stagetransf["STAGE IIIB"] = "st3"
stagetransf["STAGE III"] = "st3"
stagetransf["STAGE II"] = "st2"

datclin$stage = stagetransf[datclin$AJCC_PATHOLOGIC_TUMOR_STAGE]
table(datclin$stage, useNA = "a")

## subtype
table(datclin$SUBTYPE, useNA = "a")
x = datclin$SUBTYPE
x = gsub("BRCA_","",x)
x[which(x == "")] = "Other"
x[which(x == "Normal")] = "Other"
table(x,useNA="a")



###############################################
### RNAseq

datr1 = read.table(dat.rna.file, sep = "\t", header = T)
sum(datr1$Hugo_Symbol == "")  ## 13
datx = datr1[which(datr1$Hugo_Symbol != ""),]
datx = datx[,-2]  ## remove the entrez column

## for duplicated gene symbols, average the profiles
sum(duplicated(datx$Hugo_Symbol))  # 7
dupes = datx$Hugo_Symbol[duplicated(datx$Hugo_Symbol)]
xdup = datx[datx$Hugo_Symbol %in% dupes,]

x1 = datx[which(!datx$Hugo_Symbol %in% dupes),]
rownames(x1) = x1$Hugo_Symbol
x1 = x1[,-1]

x2 = data.frame()
sam = colnames(x1)
for (g in dupes){
   x2[g,sam] = apply(xdup[xdup$Hugo_Symbol == g,sam],2,mean)
}

dim(x1)
dim(x2)
sum(colnames(x1) == colnames(x2))
dat2 = rbind(x1,x2)
dat2 = as.matrix(dat2)
is.numeric(dat2)
hist(dat2)
dat2 = log2(dat2+1)
colnames(dat2) = gsub(".01$", "", colnames(dat2), perl=T)


dim(datclin)
dim(dat2)
length(intersect(colnames(dat2), rownames(datclin)))  ## 1082
datclin = datclin[colnames(dat2),]


###############################################
### ssGSEA scores

library(GSVA)

dgsva = ssgseaParam(dat2, gsets, minSize=3, normalize = T)
sco1 = gsva(dgsva)
hist(sco1, breaks = 50)

## add scores to clinical data frame
x = t(sco1)
datclin[rownames(x),colnames(x)] = x
dim(datclin)

#############
### subtypes boxplots

library(ggplot2)
library(ggpubr)

dd = as.data.frame(as.table(sco1))  ## long df
dd[1:10,]
colnames(dd) = c("gene.set","sample","enr.score")
dd$gene.set = factor(dd$gene.set, levels = names(gsets))
dd$subtype = datclin[dd$sample,"subtype"]
dd$subtype = factor(dd$subtype, levels = c("Basal","Her2","LumA","LumB","Other"))

ggplot(dd, aes(x = subtype, y = enr.score, fill = subtype))+
   geom_violin(width = 0.8, alpha = 0.43, scale = "width", color = "white")+
   geom_boxplot(outliers=F, linewidth=1, width = 0.4)+
   stat_compare_means(size = 2.5)+
   theme_bw() +
   theme(strip.text = element_text(size=13),
      axis.text.x = element_text(angle=0, hjust = 0.5,face = "bold",size = 9),
      axis.title.x = element_text(size = 13,face="bold"),
      axis.title.y = element_text(size = 13,face="bold") )+
   facet_wrap(~gene.set, ncol = 4, scale = "free_y")
ggsave("spat.signatures_vs_subtypes_v1.pdf", height = 3.5, width = 11)


###############################################
###############################################
## survival analysis

library(survival)

## filter the clinical data
df1 = datclin[which(!is.na(datclin$surv_status) & !is.na(datclin$stage) & datclin$stage != "stX"),]
dim(df1)


dfx = data.frame()  ## records p-values and Cox regression coefficients

## v1 = signature only
for (gs in names(gsets)){
   nm = paste0("v1.",gs)
   form = as.formula( paste0("Surv(surv_time, surv_status) ~ ", gs) )
   coxph.fit <- coxph(form, data = df1)

   x = summary(coxph.fit)$coefficient
   x = x[,c("Pr(>|z|)","coef"), drop = F]
   colnames(x) = c(paste0("Pval_",nm),paste0("coef_",nm))
   
   dfx[rownames(x),colnames(x)] = x
}

## v2 = signature + stage
for (gs in names(gsets)){
   nm = paste0("v2.",gs)
   form = as.formula( paste0("Surv(surv_time, surv_status) ~ ", gs, " + stage") )
   coxph.fit <- coxph(form, data = df1)
   
   x = summary(coxph.fit)$coefficient
   x = x[,c("Pr(>|z|)","coef"), drop = F]
   colnames(x) = c(paste0("Pval_",nm),paste0("coef_",nm))
   
   dfx[rownames(x),colnames(x)] = x
}


## v3 = signature + stage + subtype
for (gs in names(gsets)){
   nm = paste0("v3.",gs)
   form = as.formula( paste0("Surv(surv_time, surv_status) ~ ", gs, " + stage + subtype") )
   coxph.fit <- coxph(form, data = df1)
   
   x = summary(coxph.fit)$coefficient
   x = x[,c("Pr(>|z|)","coef"), drop = F]
   colnames(x) = c(paste0("Pval_",nm),paste0("coef_",nm))
   
   dfx[rownames(x),colnames(x)] = x
}


## v4 = signature + stage + subtypes (the 4 subtypes together)
nm = paste0("v4")
form = as.formula( paste0("Surv(surv_time, surv_status) ~ m1 + m2 + m3 + m4 + stage + subtype") )
coxph.fit <- coxph(form, data = df1)

x = summary(coxph.fit)$coefficient
x = x[,c("Pr(>|z|)","coef"), drop = F]
colnames(x) = c(paste0("Pval_",nm),paste0("coef_",nm))

dfx[rownames(x),colnames(x)] = x


## v5 = stage + subtypes (no signatures)
nm = paste0("v5")
form = as.formula( paste0("Surv(surv_time, surv_status) ~ stage + subtype") )
coxph.fit <- coxph(form, data = df1)

x = summary(coxph.fit)$coefficient
x = x[,c("Pr(>|z|)","coef"), drop = F]
colnames(x) = c(paste0("Pval_",nm),paste0("coef_",nm))

dfx[rownames(x),colnames(x)] = x

library(openxlsx)
write.xlsx(dfx,"surv_analysis.xlsx", rowNames = T)


#######################
## Kaplan Mayer (KM) plots
library(ggsurvfit)
library(survminer)

## plotting m3 only here
varx = "m3"
delim = quantile(df1[,varx],0.5)
df1$bin = ifelse(df1[,varx] >= delim, "High", "Low")
table(df1$bin)

## unadjusted KM plot:
p <- survfit2(Surv(surv_time,surv_status) ~ bin, data=df1) 
ggsurvplot(p, type = "survival",pval=TRUE) 
ggsave(paste0("KM_m3_median_unadjusted.pdf"), height = 5.5, width = 8)

## adjusted plot (covar = stage + subtype)
fit <- coxph(Surv(surv_time, surv_status) ~ bin + stage + subtype, data = df1)
pval = signif(summary(fit)$coefficients[1,5],3)
ggadjustedcurves(fit, data = df1, variable = "bin", method = "average", 
                 xlab = "Months", ylab = "Adjusted PFS probability",
                 title = paste0("PFS curves adjusted for stage & subtype; Signature = ", varx, 
                                "\nStratification = median; Pval = ",pval) )
ggsave(paste0("KM_adj_m3_median.pdf"), height = 6, width = 8.7)
