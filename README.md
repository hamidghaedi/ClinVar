# ClinVar variant class prediction

To classify clinvar variants into two diffrent groups ("conflict", and "no-conflict") based on the variant annotations (original clinvar annotation + VEP annotation)

### Summary on the variants review status in ClinVar (*clinvar_20210619.vcf.gz*)



|     stars       | Review status  | number of variant   |Description |
| :-------------- |:--------------|:--------------------|:----------|
|four	 | practice guideline	                                     |51  |    practice guideline |
|three	 | reviewed by expert panel	                               |11660 |  reviewed by expert panel |
|two	   |criteria provided, multiple submitter, no conflicts	  |138617  |  Two or more submitters with assertion criteria and evidence (or a public contact) provided the same interpretation. |
|one	   | criteria provided, conflicting interpretations	         |49026  |  Multiple submitters provided assertion criteria and evidence (or a public contact) but there are conflicting interpretations. The independent values are enumerated for clinical significance. |
|one	   |criteria provided, single submitter	                  |696333 |  One submitter provided an interpretation with assertion criteria and evidence (or a public contact). |
|none	 | no assertion for the individual variant	              |   -  |   The allele was not interpreted directly in any submission; it was submitted to ClinVar only as a component of a haplotype or a genotype. |
|none	 | no assertion criteria provided	                         |53236 |  The allele was included in a submission with an interpretation but without assertion criteria and evidence (or a public contact). |
|none	 | no assertion provided	                                 |10653  | The allele was included in a submission that did not provide an interpretation. |




How many variants in each class:

Review status pie chart                 | 
:----------------------------------------------------------------------------------------------------------------:|
<img src="https://raw.githubusercontent.com/hamidghaedi/clinvar/main/figs/class_dist.PNG?token=AQUBCE7EVI2C2T4CU23RZRTBBLVKS" width="800" height="400"> |

Receiving more submissions from the community, variants in the class "criteria provided, single submitter" will turn into one of the "criteria provided, multiple submitter, no conflicts" or "criteria provided, conflicting interpretations" classes. To find out which features  of variants may contibute to this conversion, I am doing this project.


--------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
#### file download and preprocessing
ClinVar variant file were download from ClinVar FTP under the name "clinvar_20210619.vcf.gz". In order to convert format , file were converted to a tab delimited file by
```bash
#ComputeCanada instance
module load StdEnv/2020 vcflib
vcf2tsv clinvar_20210619.vcf.gz > clinvar_tsv.txt
```
To add annotation to the ClinVar file,  online version of VEP was used. To return one line per variant, filter setting was set to "show most severe consequence per variant".
This code block shows reading converted ClinVar , defining a "class" variable, reading VEP annotation, making a subste ready for the rest of analysis.

```R
df <- data.table::fread("~/clinvar/clinvar_tsv.txt")

# definig a class variable

df$class = ifelse(df$CLNREVSTAT == "criteria_provided,_conflicting_interpretations", 1,
                   ifelse(df$CLNREVSTAT == "criteria_provided,_multiple_submitters,_no_conflicts", 0,
                          ifelse(df$CLNREVSTAT == "criteria_provided,_single_submitter", 2,3)))

subDf <- df[df$CLNREVSTAT %in% c("criteria_provided,_conflicting_interpretations",
                                 "criteria_provided,_multiple_submitters,_no_conflicts"), ]

#write.table(subDf, file = "~/clinvar/filteredClinVar.txt", sep = "\t", row.names = F)
# reading VEP output
vep <- data.table::fread("~/clinvar/vep_2021-07-07.txt", check.names = F)
vep = data.frame(vep)

names(vep)[1] <- "ID"

subVep <- vep[vep$ID %in% subDf$ID,]
#rm(vep)

# joining datasets
subClin = dplyr::left_join(subDf, subVep, by = "ID")
#saveRDS(subClin, file = "~/clinvar/subClin.rds")

subClin <- readRDS("~/clinvar/subClin.rds")
subClin <- data.frame(subClin)
# :), changing my mind! so new encoding:
subClin$class = ifelse(subClin$class == 0, "noConflict", "conflict")

# how many of each group
table(subClin$class)
# 0       1 
# 138617  49026 
```

#### Analysis: allele frequecny features
```R

# allele frequency data in dataset
#### data source clinvar
# AF_EXAC: the ExAC Consortium
# AF_ESP : GO-ESP [https://esp.gs.washington.edu/drupal/]
# AF_TGP : the 1000 Genomes Project

##### data source added by vep
# col 76 to 92
# col 121 to 133
# col 164 , 166
# col 176 to 223
# col 292 to 295
# col 321 to 434
# 


#### operate on numerical columns
# checking mode for columns
sapply(subClin, mode)


#
library(dplyr)
col = c()
pVal  = c()
mean0 = c()
mean1 = c()
median0 = c()
median1 = c()


for(i in 1:ncol(subClin)){
  if(all(is.na((as.numeric(subClin[,i])))) == FALSE){
    col[i] <- colnames(subClin)[i]
    val = as.numeric(subClin[,i])
    class = as.factor(subClin$class)
    df <- data.frame(val = val, class = class)
    pVal[i] = wilcox.test(val ~ class, data = df, exact = FALSE)[["p.value"]]
    mean0[i] = mean(df$val[df$class== 0], na.rm = T)
    mean1[i] = mean(df$val[df$class== 1], na.rm = T)
    median0[i] = median(df$val[df$class== 0], na.rm = T)
    median1[i] = median(df$val[df$class== 1], na.rm = T)
  }
}

resultTable <- data.frame(col = col,
                          pVal  = pVal,
                          mean0 = mean0,
                          mean1 = mean1,
                          median0 = median0,
                          median1 = median1)



# selecting allele frequency columns

af <- subClin[, c(30,8:10,76:92,121:133,164:166,176:223,292:295,321:434)]
# convert . and - to NA
names(af)
af$AF_ESP[af$AF_ESP == "."] <- NA
af$AF_EXAC[af$AF_EXAC == "."] <- NA
af$AF_TGP[af$AF_TGP == "."] <- NA

for (i in 5:ncol(af)){
  af[,i][af[,i]== "-"] = NA
}

# 
sapply(af, mode)
# convert to numbers
for (i in 2:ncol(af)){
  af[,i] = as.numeric(af[,i])
}

# set class to factor
af$class <- factor(af$class, levels = c("0", "1"), labels = c("no-conflict", "conflict"))


# test
# to see if any diffrences between NA cases and those that have a value
subAfMelted$binAf <- ifelse(is.na(subAfMelted$value), "notAvail", "Avail")

table(subAfMelted$class, subAfMelted$binAf)

#               Avail   notAvail
# no-conflict   161655   254196
# conflict       71326    75752

# so there is no genetically meaningfull inference from the above table; half of conflict
# are avail and another half is notAvail! 

# numerical testing:
t.test(AF_ESP ~ class, data = af, alternative = "two.sided", var.equal = FALSE)
# data:  AF_ESP by class
# t = 74.582, df = 47857, p-value < 2.2e-16
# alternative hypothesis: true difference in means is not equal to 0
# 95 percent    confidence interval:
# 0.05431033    0.05724193
# sample estimates:
#   mean in group no-conflict     mean in group conflict 
#   0.057109283                   0.001333154 
t.test(AF_EXAC ~ class, data = af, alternative = "two.sided", var.equal = FALSE)
# Welch Two Sample t-test
# data:  AF_EXAC by class
# t = 67.847, df = 72589, p-value < 2.2e-16
# alternative hypothesis: true difference in means is not equal to 0
# 95  percent confidence interval:
#     0.03261408 0.03455449
# sample estimates:
#   mean in group no-conflict    mean in group conflict 
#   0.034426714                   0.000842429 

t.test(AF_TGP ~ class, data = af, alternative = "two.sided", var.equal = FALSE)
# Welch Two Sample t-test
# 
# data:  AF_TGP by class
# t = 78.863, df = 45140, p-value < 2.2e-16
# alternative hypothesis: true difference in means is not equal to 0
# 95 percent confidence interval:
#   0.06362815 0.06687151
# sample estimates:
#   mean in group no-conflict       mean in group conflict 
#   0.066950827                    0.001700998 

t.test(gnomAD_AF ~ class, data = af, alternative = "two.sided", var.equal = FALSE)
# Welch Two Sample t-test
# 
# data:  gnomAD_AF by class
# t = 71.339, df = 107892, p-value < 2.2e-16
# alternative hypothesis: true difference in means is not equal to 0
# 95 percent confidence interval:
#   0.02695596 0.02847898
# sample estimates:
#   mean in group no-conflict    mean in group conflict 
# 0.028838145               0.001120677 

t.test(ExAC_Adj_AF ~ class, data = af, alternative = "two.sided", var.equal = FALSE)
# Welch Two Sample t-test
# 
# data:  ExAC_Adj_AF by class
# t = 35.911, df = 42375, p-value < 2.2e-16
# alternative hypothesis: true difference in means is not equal to 0
# 95 percent confidence interval:
#   0.01622607 0.01809958
# sample estimates:
#   mean in group no-conflict    mean in group conflict 
# 0.0180856399                    0.0009228158 

t.test(ExAC_nonTCGA_AF ~ class, data = af, alternative = "two.sided", var.equal = FALSE)
# Welch Two Sample t-test
# 
# data:  ExAC_nonTCGA_AF by class
# t = 35.457, df = 40302, p-value < 2.2e-16
# alternative hypothesis: true difference in means is not equal to 0
# 95 percent confidence interval:
#   0.01641607 0.01833721
# sample estimates:
#   mean in group no-conflict    mean in group conflict 
# 0.0182245308              0.0008478884 
t.test(ExAC_nonTCGA_Adj_AF ~ class, data = af, alternative = "two.sided", var.equal = FALSE)
# Welch Two Sample t-test
# 
# data:  ExAC_nonTCGA_Adj_AF by class
# t = 36.034, df = 40387, p-value < 2.2e-16
# alternative hypothesis: true difference in means is not equal to 0
# 95 percent confidence interval:
#   0.01704728 0.01900847
# sample estimates:
#   mean in group no-conflict    mean in group conflict 
# 0.0189630638                   0.0009351883


# subAf
subAf <- af[, c(1:4,13,45,55,61)]

subAfMelted <- reshape2::melt(subAf)
# saving as a csv
#write.table(subAfMelted, file = "~/clinvar/alleleFreqClinvarData.csv", sep = ",", row.names = F, quote = F)

tiff(filename = "freq_boxplot.tiff", width = 8, height = 5, units = "in", res = 300)
ggplot(subAfMelted, aes(x = variable, y = log10(value))) + 
  geom_boxplot(aes(fill = class), position = position_dodge(0.9)) +
  scale_fill_manual(values = c("#999999", "#E69F00")) + 
  labs(x="", y= "-log10(AF)", fill = "class") +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 45, vjust = 0.5, hjust=0.5))
dev.off()
```
<img src="https://raw.githubusercontent.com/hamidghaedi/clinvar/main/figs/frq_boxplot.PNG?token=AQUBCE6GOCAU6CL4WSOCQYDBBLWK2" width="800" height="400">
```R
#### looking at AUC performance
library(pROC)
auc <-c()
varName <- c()
auc_ci_ll <- c()
auc_ci_ul <- c()


for (i in 2:ncol(subAf)){
  df <- data.frame(class = subAf$class, 
                   val = subAf[,i])
  df <- df[!is.na(df$val), ]
  d <- roc(class ~ val, data=df, percent = TRUE, ci = TRUE, auc= TRUE)
  varName[i] <- colnames(subAf)[i]
  auc[i] <- d$auc
  auc_ci_ll[i] <- round(data.frame(d$ci)[1,],2)
  auc_ci_ul[i] <- round(data.frame(d$ci)[3,],2)
}

res_auc <- data.frame(variable = varName,
                              AUC = auc,
                              LL_CI =auc_ci_ll ,
                              UL_CI = auc_ci_ul)

# using AF_TGP , AUC [95% CI]: 73.23 [72.86 - 73.61]
## calculating summary statistics fro groups 

varName <- c()
mean_conflict = c()
mean_noConflict = c()
n_conflict= c()
n_noConflict= c()
median_conflict = c()
median_noConflict = c()
sd_Conflict = c()
sd_noConflict = c()

for (i in 2:ncol(af)){
  varName[i] = colnames(af)[i]
  d = data.frame(cls = af$class, val = af[,i])
  d = d[!is.na(d$val),]
  mean_conflict[i] = mean(d$val[d$cls == "conflict"])
  mean_noConflict[i] = mean(d$val[d$cls == "no-conflict"])
  n_conflict[i] = length(d$val[d$cls == "conflict"])
  n_noConflict[i] = length(d$val[d$cls == "no-conflict"])
  median_conflict[i] = median(d$val[d$cls == "conflict"])
  median_noConflict[i] = median(d$val[d$cls == "no-conflict"])
  sd_Conflict[i] = sd(d$val[d$cls == "conflict"])
  sd_noConflict[i] = sd(d$val[d$cls == "no-conflict"])
}

resDF <- data.frame(population = varName,
                    n_conflict = n_conflict,
                    n_noConflict= n_noConflict,
                    mean_conflict = mean_conflict,
                    mean_noConflict = mean_noConflict,
                    median_conflict = median_conflict,
                    median_noConflict = median_noConflict,
                    sd_Conflict = sd_Conflict,
                    sd_noConflict = sd_noConflict)


# looking at sub-populations NON gnomAD_EXAC
afSubPop = af[, c(1,6:12, 22:37, 87,89)]
afSubPop = afSubPop[, c(-10, -12, -14, -16, -18, -20, -23)]
afSubPop = afSubPop[, c(-9, -10, -11, -12, -13, -14)]

afSubPopMelted <- reshape2::melt(afSubPop)
afSubPopMelted = afSubPopMelted[!is.na(afSubPopMelted$value),]

tiff(filename = "~/clinvar/subPopAF_non-gnomAD_boxplot.tiff", width = 8, height = 5, units = "in", res = 300)
ggplot(afSubPopMelted, aes(x = variable, y = log10(value))) + 
  geom_boxplot(aes(fill = class), position = position_dodge(0.9)) +
  scale_fill_manual(values = c("#999999", "#E69F00")) + 
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=0.5))
  theme_bw()
  dev.off()
  ```
<img src="https://raw.githubusercontent.com/hamidghaedi/clinvar/main/figs/sub_pop_allele-freq.PNG?token=AQUBCE5KFTPLRQTIWIFWGGLBBLWS4" width="800" height="400">

#### Analysis: variant consequnce annotation
#### Pathogenicity prediction scores
#### Conservation scores
