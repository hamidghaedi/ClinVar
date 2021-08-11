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
```R
# consequences
unique(subClin$Consequence)

# [1] "missense_variant"                                                                                       
# [2] "synonymous_variant"                                                                                     
# [3] "intron_variant"                                                                                         
# [4] "splice_region_variant,synonymous_variant"                                                               
# [5] "splice_region_variant,intron_variant"                                                                   
# [6] "3_prime_UTR_variant"                                                                                    
# [7] "inframe_insertion"                                                                                      
# [8] "frameshift_variant"                                                                                     
# [9] "downstream_gene_variant"                                                                                
# [10] "5_prime_UTR_variant"                                                                                    
# [11] "inframe_deletion"                                                                                       
# [12] "start_lost"                                                                                             
# [13] "missense_variant,splice_region_variant"                                                                 
# [14] NA                                                                                                       
# [15] "splice_donor_variant"                                                                                   
# [16] "stop_gained"                                                                                            
# [17] "stop_gained,splice_region_variant"                                                                      
# [18] "upstream_gene_variant"                                                                                  
# [19] "splice_acceptor_variant"                                                                                
# [20] "splice_region_variant,5_prime_UTR_variant"                                                              
# [21] "frameshift_variant,splice_region_variant"                                                               
# [22] "frameshift_variant,stop_lost"                                                                           
# [23] "start_lost,5_prime_UTR_variant"                                                                         
# [24] "frameshift_variant,start_lost"                                                                          
# [25] "splice_region_variant,intron_variant,NMD_transcript_variant"                                            
# [26] "intron_variant,NMD_transcript_variant"                                                                  
# [27] "missense_variant,NMD_transcript_variant"                                                                
# [28] "synonymous_variant,NMD_transcript_variant"                                                              
# [29] "3_prime_UTR_variant,NMD_transcript_variant"                                                             
# [30] "splice_region_variant,3_prime_UTR_variant,NMD_transcript_variant"                                       
# [31] "missense_variant,splice_region_variant,NMD_transcript_variant"                                          
# [32] "frameshift_variant,NMD_transcript_variant"                                                              
# [33] "stop_gained,NMD_transcript_variant"                                                                     
# [34] "splice_region_variant,synonymous_variant,NMD_transcript_variant"                                        
# [35] "splice_acceptor_variant,NMD_transcript_variant"                                                         
# [36] "splice_donor_variant,intron_variant"                                                                    
# [37] "splice_acceptor_variant,intron_variant"                                                                 
# [38] "stop_gained,frameshift_variant"                                                                         
# [39] "protein_altering_variant"                                                                               
# [40] "splice_donor_variant,coding_sequence_variant,intron_variant"                                            
# [41] "stop_lost"                                                                                              
# [42] "inframe_deletion,splice_region_variant"                                                                 
# [43] "splice_donor_variant,coding_sequence_variant"                                                           
# [44] "splice_acceptor_variant,coding_sequence_variant"                                                        
# [45] "splice_acceptor_variant,coding_sequence_variant,intron_variant"                                         
# [46] "splice_donor_variant,NMD_transcript_variant"                                                            
# [47] "stop_retained_variant"                                                                                  
# [48] "stop_gained,inframe_deletion"                                                                           
# [49] "stop_gained,frameshift_variant,splice_region_variant"                                                   
# [50] "stop_lost,3_prime_UTR_variant"                                                                          
# [51] "frameshift_variant,stop_retained_variant"                                                               
# [52] "inframe_insertion,NMD_transcript_variant"                                                               
# [53] "inframe_deletion,NMD_transcript_variant"                                                                
# [54] "intron_variant,non_coding_transcript_variant"                                                           
# [55] "inframe_insertion,splice_region_variant"                                                                
# [56] "start_lost,splice_region_variant"                                                                       
# [57] "5_prime_UTR_variant,NMD_transcript_variant"                                                             
# [58] "non_coding_transcript_exon_variant"                                                                     
# [59] "splice_donor_variant,5_prime_UTR_variant,intron_variant"                                                
# [60] "stop_gained,protein_altering_variant"                                                                   
# [61] "coding_sequence_variant,intron_variant"                                                                 
# [62] "protein_altering_variant,splice_region_variant"                                                         
# [63] "splice_region_variant,coding_sequence_variant,intron_variant"                                           
# [64] "stop_retained_variant,3_prime_UTR_variant"                                                              
# [65] "splice_acceptor_variant,intron_variant,NMD_transcript_variant"                                          
# [66] "splice_donor_variant,3_prime_UTR_variant,NMD_transcript_variant"                                        
# [67] "splice_donor_variant,splice_acceptor_variant,coding_sequence_variant,intron_variant"                    
# [68] "splice_donor_variant,coding_sequence_variant,5_prime_UTR_variant,intron_variant,NMD_transcript_variant" 
# [69] "protein_altering_variant,NMD_transcript_variant"                                                        
# [70] "start_lost,NMD_transcript_variant"                                                                      
# [71] "splice_region_variant,non_coding_transcript_exon_variant"                                               
# [72] "stop_gained,inframe_insertion"                                                                          
# [73] "frameshift_variant,splice_region_variant,NMD_transcript_variant"                                        
# [74] "splice_region_variant,3_prime_UTR_variant"                                                              
# [75] "splice_donor_variant,splice_acceptor_variant,coding_sequence_variant,5_prime_UTR_variant,intron_variant"
# [76] "transcript_ablation" 

tmp = data.frame(unclass(table(subClin$class, subClin$Consequence)))
rownames(tmp) <- c("noConflict", "conflict")
tmp = data.frame(t(tmp))

tmp$conflict_precent = round(tmp$conflict/sum(tmp$conflict)*100,2)
tmp$noConflict_precent = round(tmp$noConflict/sum(tmp$noConflict)*100,2)

filtTmp = tmp[tmp$conflict_precent >= 1 | tmp$noConflict_precent >= 1, ]
filtTmp = filtTmp[, c(-1,-2)]
rownames(filtTmp)[1] <- "3_prime_UTR_variant"

chisq <- chisq.test(filtTmp)
chisq

library(corrplot)
png(filename = "~/clinvar/conseq.png", width = 7, height = 10, units = "in", res = 300)
corrplot::corrplot(chisq$residuals, is.cor = FALSE, cl.pos = 'n')
dev.off()
```
<img src="https://raw.githubusercontent.com/hamidghaedi/clinvar/main/figs/conseq.png?token=AQUBCE64NO4HCEAQ24OEB5DBBLXFI" width="500" height="800">

```R
# IMPACT

table(subClin$class, subClin$IMPACT)

tmp = data.frame(unclass(table(subClin$class, subClin$IMPACT)))
tmp = data.frame(t(tmp))

tmp$conflict_precent = round(tmp$conflict/sum(tmp$conflict)*100,2)
tmp$noConflict_precent = round(tmp$noConflict/sum(tmp$noConflict)*100,2)


chisq <- chisq.test(tmp[, c(1:2)])
chisq

png(filename = "~/clinvar/impact.png", width = 7, height = 10, units = "in", res = 300)
corrplot::corrplot(chisq$residuals, is.cor = FALSE, cl.pos = 'n')
dev.off()
```
<img src="https://raw.githubusercontent.com/hamidghaedi/clinvar/main/figs/impact.png?token=AQUBCE4KBL6DJBTX4SYADVTBBLXNE" width="500" height="600">

```R
# transcript length

subClin <- readRDS("~/clinvar/subClin.rds")
subClin <- data.frame(subClin)
subClin$class = ifelse(subClin$class == 0, "noConflict", "conflict")


trs = unique(subClin$Feature)

library(biomaRt)

mart <- useDataset("hsapiens_gene_ensembl", useMart("ensembl"))
grep("length", listAttributes(mart = mart)[,1])
atr = c("ensembl_transcript_id_version", "transcript_start", "transcript_end", "transcript_length")
flts ="ensembl_transcript_id_version" 

transMart = getBM(filters= flts, 
                  attributes= atr,
                  values=trs, 
                  mart= mart)
names(transMart)[1] <- "Feature"

tmp = subClin[, c(30, 38,40,41,44,45,46)]
tmp = dplyr::left_join(tmp, transMart)
# tmp = tmp[,c(1,3,4,7)]
tmp = tmp[!is.na(tmp$transcript_length),]

## Calculate bin breaks for numeric variables with respect to their relationships with the outcome variable class
bins =scorecard::woebin(tmp[, c('transcript_length', 'class')], y = 'class', positive = 'conflict')

# visualization
scorecard::woebin_plot(bins$transcript_length)$transcript_length

tmp = tmp[tmp$transcript_length < 40000, ]

t.test(transcript_length ~ class, data = tmp, alternative = "two.sided", var.equal = FALSE)


library(ggplot2)
png(filename = "~/clinvar/transcript_length.png", width = 7, height = 10, units = "in", res = 300)
ggplot(tmp, aes(x = as.factor(class), y = log10(transcript_length))) + 
  geom_boxplot(aes(fill = class), position = position_dodge(0.9)) +
  scale_fill_manual(values = c("#999999", "#E69F00")) + 
  theme_bw()
dev.off()
```
<img src="https://raw.githubusercontent.com/hamidghaedi/clinvar/main/figs/transcript_length.png?token=AQUBCE4I2PWFUIDDFO6UIEDBBLX4G" width="400" height="500">

```R
# variant location in exons

table(subClin$EXON == "-")

# FALSE   TRUE 
# 153639  33854

#tmp = subClin[, c(30, 40, 41)]
tmp = tmp[!is.na(tmp$EXON),]


tmp$ratio = tmp$EXON
tmp$ratio = ifelse(tmp$ratio != "-", tmp$ratio, tmp$INTRON)
tmp = cbind(tmp, stringr::str_split_fixed(tmp$ratio, "/", 2))
names(tmp)[c(12:13)] <- c("nExon", "tExon")
tmp$nratio = as.numeric(tmp$nExon)/as.numeric(tmp$tExon)
# those variants with "-" as Exon and Intron will be removed since they generate NA
# those variants encompass more than one Exon/Intron will be removed since they generate NA
colSums(is.na(tmp))
# class   EXON INTRON  ratio  nExon  tExon nratio 
# 0       0      0      0      0      0     7698
tmp = tmp[!is.na(tmp$nratio),]

ggplot(tmp, aes(x = as.factor(class), y = nratio)) + 
  geom_boxplot(aes(fill = class), position = position_dodge(0.9)) +
  scale_fill_manual(values = c("#999999", "#E69F00")) + 
  theme_bw()

# discertization visualization
## Import libraries
library(scorecard)
library(ggplot2)

# read more on discertization : https://nextjournal.com/eda/discretize-cont-var
## Calculate bin breaks for numeric variables with respect to their relationships with the outcome variable class
bins = scorecard::woebin(tmp[, c('nratio', 'class')], y = 'class', positive = 'conflict')

# visualization
scorecard::woebin_plot(bins$nratio)$nratio
# interpertaion for the graph: highest probability for conflicting group is found in bin 0.96 and greater


#############################other position based annotation#################################################
tmp$cDNA_position = as.numeric(tmp$cDNA_position)
tmp$CDS_position = as.numeric(tmp$CDS_position)
tmp$Protein_position = as.numeric(tmp$Protein_position)
tmp = tmp[!is.na(tmp$cDNA_position),]

bins =scorecard::woebin(tmp[, c('cDNA_position','CDS_position','Protein_position','class')], y = 'class', positive = 'conflict')
# visualization
scorecard::woebin_plot(bins$cDNA_position)
scorecard::woebin_plot(bins$CDS_position)
scorecard::woebin_plot(bins$Protein_position)

## calculate correlation between position based metrics
summary(tmp)

selc = c("cDNA_position","CDS_position","Protein_position","transcript_length", "nratio")

tmp = tmp[, selc]
for (i in 1:5){
  tmp[,i] <- as.numeric(tmp[,i])
}

tmp = na.omit(tmp)

# calculate correlation matrix
res <- cor(tmp)

png(filename = "~/clinvar/postion_metrics_corrplot.png", width = 7, height = 10, units = "in", res = 300)
corrplot::corrplot(res, type = "upper", order = "hclust", addCoef.col = 'black', 
         tl.col = "black", tl.srt = 45)
dev.off()
```
<img src="https://raw.githubusercontent.com/hamidghaedi/clinvar/main/figs/postion_metrics_corrplot.png?token=AQUBCEY6TJQZ34E5DTDGBR3BBLYJK" width="400" height="500">


#### Pathogenicity prediction scores

The version of dbNSFP that VEP  used to annotate ClinVar file (version 4)  compiles prediction scores from 37 prediction algorithms (SIFT, SIFT4G, Polyphen2-HDIV, Polyphen2-HVAR, LRT, MutationTaster2, MutationAssessor, FATHMM, MetaSVM, MetaLR, CADD, CADD_hg19, VEST4, PROVEAN, FATHMM-MKL coding, FATHMM-XF coding, fitCons x 4, LINSIGHT, DANN, GenoCanyon, Eigen, Eigen-PC, M-CAP, REVEL, MutPred, MVP, MPC, PrimateAI, GEOGEN2, BayesDel_addAF, BayesDel_noAF, ClinPred, LIST-S2, ALoFT). Since different score has a different scaling system, the dbNSFP developers have created a rank score for each score so that it is comparable between scores. This can help with easier interpertaion of the result. The rank score has a scale 0 to 1 and shows the percentage of scores that are less damaging in dbNSFP, e.g., a rank score of 0.9 means the top 10% most damaging.

Here I first tried to identify columns in the dataset that provide data on pathogenecity prediction score, visualized scores between groups , calculated correlation matrix to remove highly correlated features(r > 0.7) and again visualized those selected feature. 

```R
# pathogenicity  score 

subClin <- readRDS("~/clinvar/subClin.rds")
subClin <- data.frame(subClin)
subClin$class = ifelse(subClin$class == 0, "noConflict", "conflict")

#knowing pathogeneicty prediction columns are providing a kind of "score", so :
pathPredict = c(pathPredict,colnames(subClin)[grepl("score", colnames(subClin))] )

patDat = subClin[, which(names(subClin) %in% pathPredict)]

# converting sift and polyphen to numbers
df = cbind(t(data.frame(strsplit(patDat$SIFT, "\\("))), t(data.frame(strsplit(patDat$PolyPhen, "\\("))))
colnames(df) <- c("SIFT_Pred", "SIFT_score", "PolyPhen_Pred", "PolyPhen_score")
df = data.frame(df)
df$SIFT_score = as.numeric(sub("\\)", "", df$SIFT_score))
df$PolyPhen_score = as.numeric(sub("\\)", "", df$PolyPhen_score))
rownames(df) <- NULL

# removing original SIFT and polyPhen columns
patDat = cbind(patDat, df)
patDat = patDat[, -c(2,3)]


# replacing "-" with NA
patDat[patDat == "-"] <- NA

# counting  NAs and reducing data frame size by NA numbers
patDat = patDat[, colSums(is.na(patDat)) <= 150000]
patDat = patDat[rowSums(is.na(patDat)) <= 50, ]


# droping off columns that have more than one number/chr
patDat = patDat[, -c(12,17,28,36,44,48,55,69,71)]

# knowing type of data
sapply(patDat, mode)

# setting datatype as numeric
patDat[, -1] = sapply(patDat[, -1], as.numeric)

# identifying conservation score columns
consScore = c("SiPhy_29way_logOdds_rankscore", "bStatistic_converted_rankscore", "integrated_fitCons_score",
              "phastCons100way_vertebrate_rankscore", "phastCons17way_primate_rankscore",
              "phastCons30way_mammalian_rankscore","phyloP100way_vertebrate_rankscore",
              "phyloP17way_primate_rankscore", "phyloP30way_mammalian_rankscore",
              "GERP.._RS_rankscore", "GM12878_fitCons_rankscore", "GM12878_fitCons_score",
              "H1.hESC_fitCons_rankscore", "H1.hESC_fitCons_score", "HUVEC_fitCons_rankscore",
              "HUVEC_fitCons_score", "integrated_fitCons_rankscore")

# removing conservation scores for now
patDat = patDat[, -which(colnames(patDat) %in% consScore)]

# visualization
library(corrplot)
library(RColorBrewer)
library(ggplot2)

# distribution

plotDf <- reshape2::melt(patDat, id.vars = 'class')

png(filename = "~/clinvar/path_predict_dist.png", width = 17, height = 15, units = "in", res = 300)
ggplot(plotDf, aes(x=variable, y=value, fill=class)) + 
  geom_boxplot(outlier.colour=NA) +
  scale_fill_manual(values = c("#999999", "#E69F00")) +
  facet_wrap(~variable, scale="free")
dev.off()
```
<img src="https://raw.githubusercontent.com/hamidghaedi/clinvar/main/figs/path_predict_dist.png?token=AQUBCE7XEALJ7EMOXC4RHKDBCM56O" width="500" height="500">

```R
#correlation
M <-cor(patDat[, -1], use = "complete.obs")
png(filename = "~/clinvar/path_predict_corrplot.png", width = 15, height = 15, units = "in", res = 300)
corrplot(M, type="upper", order="hclust",
         col=brewer.pal(n=8, name="RdYlBu"))
dev.off()
```
<img src="https://raw.githubusercontent.com/hamidghaedi/clinvar/main/figs/path_predict_corrplot.png?token=AQUBCE7T3DRZ6HUNSPLTRYDBCM5VG" width="500" height="500">

```R
# Keeping only one variable out of pairs which show cor > 0.7
highlyCorrelated <- caret::findCorrelation(M, cutoff=(0.7),verbose = FALSE)
important_var=colnames(patDat[,-highlyCorrelated])
tmp = patDat[, c("class", important_var)]
plotDf <- reshape2::melt(tmp, id.vars = 'class')

png(filename = "~/clinvar/select_path_predict_dist.png", width = 9, height = 7, units = "in", res = 300)
ggplot(plotDf, aes(x=variable, y=value, fill=class)) + 
  geom_boxplot(outlier.colour=NA) +
  scale_fill_manual(values = c("#999999", "#E69F00")) +
  facet_wrap(~variable, scale="free")
dev.off()
```
<img src="https://raw.githubusercontent.com/hamidghaedi/clinvar/main/figs/select_path_predict_dist.png?token=AQUBCE6QJUTZ5HJAPA3TZ5DBCM6G4" width="500" height="500">


#### Conservation scores
The dbNSFP v4, provides data on 9 diffrent conservation score algorithms. The same approach to what I undertook for pathogenecity score prediction was applied here as well. 
