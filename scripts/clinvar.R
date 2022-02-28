# clinvar processing
setwd("~/clinvar")


# data coming from this command on CC:
# module load StdEnv/2020 vcflib
# vcf2tsv clinvar_20210619.vcf.gz > clinvar_tsv.txt

df <- data.table::fread("~/clinvar/clinvar_tsv.txt")

# selecting those with multiple submitter
table(df$CLNREVSTAT)


# #stars	Review status	                                        #No     Description
#--------|---------------------------------------------------|--------|----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------|
# four	 practice guideline	                                    51     practice guideline
# three	 reviewed by expert panel	                              11660  reviewed by expert panel
# two	   criteria provided, multiple submitter, no conflicts	 138617  Two or more submitters with assertion criteria and evidence (or a public contact) provided the same interpretation.
# one	   criteria provided, conflicting interpretations	        49026  Multiple submitters provided assertion criteria and evidence (or a public contact) but there are conflicting interpretations. The independent values are enumerated for clinical significance.
# one	   criteria provided, single submitter	                 696333  One submitter provided an interpretation with assertion criteria and evidence (or a public contact).
# none	 no assertion for the individual variant	                -    The allele was not interpreted directly in any submission; it was submitted to ClinVar only as a component of a haplotype or a genotype.
# none	 no assertion criteria provided	                        53236  The allele was included in a submission with an interpretation but without assertion criteria and evidence (or a public contact).
# none	 no assertion provided	                                10653  The allele was included in a submission that did not provide an interpretation.

# making a chart 
d <- data.frame(table(df$CLNREVSTAT))
d$abstract <- c("criteria_provided,_conflicting_interpretations",
                "criteria_provided,_multiple_submitters,_no_conflicts",
                "criteria_provided,_single_submitter", ".",".",".",".","." )
# visualization
pie(d$Freq, d$abstract, col = rainbow(length(d$Freq)), radius = 1)


# defining a class variable

df$class = ifelse(df$CLNREVSTAT == "criteria_provided,_conflicting_interpretations", 1,
                   ifelse(df$CLNREVSTAT == "criteria_provided,_multiple_submitters,_no_conflicts", 0,
                          ifelse(df$CLNREVSTAT == "criteria_provided,_single_submitter", 2,3)))


subDf <- df[df$CLNREVSTAT %in% c("criteria_provided,_conflicting_interpretations",
                                 "criteria_provided,_multiple_submitters,_no_conflicts"), ]

#write.table(subDf, file = "~/clinvar/filteredClinVar.txt", sep = "\t", row.names = F)
# reading vep output
# online vep was used for annotation, to return one line per variant, filter setting was set to 
# "show most severe consequence per variant"

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
subClin$class = ifelse(subClin$class == 0, "noConflict", "conflict")
subClin = na.omit(subClin)

# how many of each group
table(subClin$class)
# 0       1 
# 138617  49026 

################################################################################
########################### EXPLORATORY ANALYSIS#################################


##################### allele frequency data in dataset###################
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



# selecting allele frequency columns
popFreq = c(8L, 9L, 10L, 76L, 77L, 78L, 79L, 80L, 81L, 82L, 83L, 84L, 85L, 86L, 87L, 88L, 89L, 90L, 91L, 92L, 120L, 121L, 122L, 123L, 124L, 125L, 126L, 127L, 128L, 129L, 130L, 131L, 132L, 133L, 163L, 164L, 165L, 166L, 176L, 177L, 178L, 179L, 180L, 181L, 182L, 183L, 184L, 185L, 186L, 187L, 188L, 189L, 190L, 191L, 192L, 193L, 194L, 195L, 196L, 197L, 198L, 199L, 200L, 201L, 202L, 203L, 204L, 205L, 206L, 207L, 208L, 209L, 210L, 211L, 212L, 213L, 214L, 215L, 216L, 217L, 218L, 219L, 220L, 221L, 222L, 223L, 292L, 293L, 294L, 295L, 321L, 322L, 323L, 324L, 325L, 326L, 327L, 328L, 329L, 330L, 331L, 332L, 333L, 334L, 335L, 336L, 337L, 338L, 339L, 340L, 341L, 342L, 343L, 344L, 345L, 346L, 347L, 348L, 349L, 350L, 351L, 352L, 353L, 354L, 355L, 356L, 357L, 358L, 359L, 360L, 361L, 362L, 363L, 364L, 365L, 366L, 367L, 368L, 369L, 370L, 371L, 372L, 373L, 374L, 375L, 376L, 377L, 378L, 379L, 380L, 381L, 382L, 383L, 384L, 385L, 386L, 387L, 388L, 389L, 390L, 391L, 393L, 394L, 395L, 396L, 397L, 398L, 399L, 400L, 401L, 402L, 403L, 404L, 405L, 406L, 407L, 408L, 409L, 410L, 411L, 412L, 413L, 414L, 415L, 416L, 417L, 418L, 419L, 420L, 421L, 422L, 423L, 424L, 425L, 426L, 427L, 428L, 429L, 430L, 431L, 432L, 434L)
popfreqName = colnames(subClin)[popFreq]

AlleleFreq = c("class",popfreqName[grepl("AF", popfreqName)])
AlleleCount = c("class",popfreqName[grepl("AC", popfreqName)])
AlleleNumber = c("class",popfreqName[grepl("_AN", popfreqName)])
Allelenhomalt = c("class",popfreqName[grepl("_nhomalt", popfreqName)])

#AlleleFreq
###### 
tmp <- subClin[, colnames(subClin) %in% AlleleFreq]

# frequency of clinsig in diffrent subpopulations
clnsig = AlleleFreq = c("class","CLNSIG","CLNVC",popfreqName[grepl("AF", popfreqName)])

tmp = subClin[, colnames(subClin) %in% clnsig]
# tmp = tmp[tmp$CLNSIG %in% c('Benign',
#                               'Benign/Likely_benign',
#                               'Pathogenic',
#                               'Pathogenic/Likely_pathogenic',
#                               'Uncertain_significance',
#                               'Conflicting_interpretations_of_pathogenicity'),]

tmp[tmp == "."] <- NA
tmp[tmp == "-"] <- NA

# index
idx = !is.na(tmp$AF_TGP) & !is.na(tmp$gnomAD_AF)

stmp = tmp[idx,]

# subAf
subAf <- stmp[, c(5,3,1,14,2)]
subAf[,-1] = lapply(subAf[,-1], as.numeric)

subAfMelted <- reshape2::melt(subAf, id.vars = "class")
# saving as a csv
#write.table(subAfMelted, file = "~/clinvar/alleleFreqClinvarData.csv", sep = ",", row.names = F, quote = F)

png(filename = "~/clinvar/freq_boxplot_TGP_gnomAD_intersection.png", width = 8, height = 5, units = "in", res = 300)
ggplot(subAfMelted, aes(x = variable, y = log10(value))) + 
  geom_boxplot(aes(fill = class), position = position_dodge(0.9)) +
  scale_fill_manual(values = c("#999999", "#E69F00")) + 
  labs(x="", y= "-log10(AF)", fill = "class") +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 45, vjust = 0.5, hjust=0.5))
dev.off()


# TGP, gnomAD_TGP_variants, and gnomAD_non_TGP_variants

# index
idx = !is.na(tmp$AF_TGP) & !is.na(tmp$gnomAD_AF)
idx2 = is.na(tmp$AF_TGP) & !is.na(tmp$gnomAD_AF)
stmp = tmp[idx,]
stmp2  = tmp[idx2,]

stmp = data.frame(class = c(rep(stmp$class,2),stmp2$class),
                   variable = c(rep('AF_TGP',nrow(stmp)),
                             rep('gnomAD_in_TGP', nrow(stmp)),
                             rep('gnomAD_non_TGP', nrow(stmp2))),
                  value = as.numeric(c(stmp$AF_TGP, stmp$gnomAD_AF, stmp2$gnomAD_AF)))

# saving as a csv
#write.table(subAfMelted, file = "~/clinvar/alleleFreqClinvarData.csv", sep = ",", row.names = F, quote = F)

png(filename = "~/clinvar/freq_boxplot_TGP_gnomAD_in_and_not.png", width = 8, height = 5, units = "in", res = 300)
ggplot(stmp, aes(x = variable, y = log10(value))) + 
  geom_boxplot(aes(fill = class), position = position_dodge(0.9)) +
  scale_fill_manual(values = c("#999999", "#E69F00")) + 
  labs(x="", y= "-log10(AF)", fill = "class") +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 45, vjust = 0.5, hjust=0.5))
dev.off()

##

stmp$group = "gnomAD_TGP"
stmp2$group = "gnomAD_non_TGP"

gnomad = data.frame(rbind(stmp, stmp2))

table(gnomad$class, gnomad$group)
chisq.test(table(gnomad$class, gnomad$group))
prop.table(table(gnomad$class, gnomad$group), 2)

# tmpgnomAD
tgnomAD = gnomad[gnomad$CLNSIG %in% c('Benign',
                                      'Benign/Likely_benign',
                                      'Pathogenic',
                                      'Pathogenic/Likely_pathogenic',
                                      'Uncertain_significance'),]

prop.table(table(tgnomAD$CLNSIG, tgnomAD$group), 2)


#
tgp = subClin
tgp[tgp == "."] <- NA
tgp[tgp == "-"] <- NA
idx = !is.na(tgp$AF_TGP) & !is.na(tgp$gnomAD_AF)
idx2 = is.na(tgp$AF_TGP) & !is.na(tgp$gnomAD_AF)
#
tgp$group = ifelse(!is.na(tgp$AF_TGP) & !is.na(tgp$gnomAD_AF), "gnomAD_TGP", "gnomAD_non_TGP")

tgp = tgp[tgp$CLNSIG %in% c('Benign',
                                      'Benign/Likely_benign',
                                      'Pathogenic',
                                      'Pathogenic/Likely_pathogenic',
                                      'Uncertain_significance'),]
tmp_tgp = tgp[, c(30,18, 453)]

tmp_tgp %>% 
  group_by(group, class) %>%
  count(CLNSIG) %>%
  mutate(freq = n / sum(n))


# convert . and - to NA
tmp[tmp == "."] <- NA
tmp[tmp == "-"] <- NA

######_____________######
#tmp = tmp[, c(4,1:3,5:97)]


# count missing values
missNumber = data.frame(mis = colSums(is.na(tmp)))
missNumber$percent = missNumber$mis/dim(tmp)[1]

# drop columns with more than 80% missing values
tmp = tmp[, colnames(tmp) %in% rownames(missNumber)[missNumber$percent < 0.8]]

#reorder tmp column
tmp = tmp[,c(4,1:3,5:dim(tmp)[2])]

# convert to numeric
tmp[,-1] = sapply(tmp[,-1], as.numeric)


########___________________________##########
wc <- wilcox.test(tmp[,2] ~ class, data = tmp,
                   exact = FALSE)
re = data.frame(class =character(0), count =numeric(0), median=numeric(0),
                q25=numeric(0),q75 =numeric(0), population=character(0))
for(i in 2:length(tmp)){
t = data.frame(class = tmp[,1],
               pop = tmp[,i])
res = group_by(t, class) %>%
    summarise(
    count = n(),
    median = median(pop, na.rm = TRUE),
    q25 = quantile(pop,0.25, na.rm = TRUE),
    q75 = quantile(pop,0.75, na.rm = TRUE))
res$population = colnames(tmp)[i]
re = rbind(re, res)
}



######____________________________###########
# computing mean for missing values
for(i in 1:ncol(tmp[,-1])){
  tmp[,-1][is.na(tmp[,-1][,i]), i] <- mean(tmp[,-1][,i], na.rm = TRUE)
}


M <-cor(tmp[,-1])


#Keeping only one variable out of pairs which show cor > 0.7
highlyCorrelated <- caret::findCorrelation(M, cutoff=(0.7),verbose = FALSE)
important_var_af=colnames(tmp[,-highlyCorrelated])

important_var_af
# [1] "AF_TGP"                        "AF"                           
# [3] "gnomAD_exomes_controls_AFR_AF" "gnomAD_exomes_controls_ASJ_AF"
# [5] "gnomAD_genomes_AFR_AF"         "gnomAD_genomes_SAS_AF"   


# AlleleCount
#####
AlleleCount = c("class",popfreqName[grepl("AC", popfreqName)])

tmp <- subClin[, colnames(subClin) %in% AlleleCount]

# convert . and - to NA
tmp[tmp == "."] <- NA
tmp[tmp == "-"] <- NA

# count missing values
missNumber = data.frame(mis = colSums(is.na(tmp)))
missNumber$percent = missNumber$mis/dim(tmp)[1]

# drop columns with more than 80% missing values
tmp = tmp[, colnames(tmp) %in% rownames(missNumber)[missNumber$percent < 0.8]]

#reorder tmp column
tmp = tmp[,c(2,1,3:dim(tmp)[2])]

# convert to numeric
tmp[,-1] = sapply(tmp[,-1], as.numeric)

# computing mean for missing values
for(i in 1:ncol(tmp[,-1])){
  tmp[,-1][is.na(tmp[,-1][,i]), i] <- mean(tmp[,-1][,i], na.rm = TRUE)
}


M <-cor(tmp[,-1])


#Keeping only one variable out of pairs which show cor > 0.7
highlyCorrelated <- caret::findCorrelation(M, cutoff=(0.7),verbose = FALSE)
important_var_ac=colnames(tmp[,-highlyCorrelated])

important_var_ac
# [1] "class"                 "ExAC_nonTCGA_AF"       "gnomAD_exomes_NFE_AC" 
# [4] "gnomAD_genomes_SAS_AC"
important_var_ac = important_var_ac[-1]


# Allele number
#####
AlleleNumber = c("class",popfreqName[grepl("_AN", popfreqName)])

tmp <- subClin[, colnames(subClin) %in% AlleleNumber]

# convert . and - to NA
tmp[tmp == "."] <- NA
tmp[tmp == "-"] <- NA

# count missing values
missNumber = data.frame(mis = colSums(is.na(tmp)))
missNumber$percent = missNumber$mis/dim(tmp)[1]

# drop columns with more than 80% missing values
tmp = tmp[, colnames(tmp) %in% rownames(missNumber)[missNumber$percent < 0.8]]


# convert to numeric
tmp[,-1] = sapply(tmp[,-1], as.numeric)

# computing mean for missing values
for(i in 1:ncol(tmp[,-1])){
  tmp[,-1][is.na(tmp[,-1][,i]), i] <- mean(tmp[,-1][,i], na.rm = TRUE)
}


M <-cor(tmp[,-1])


#Keeping only one variable out of pairs which show cor > 0.7
highlyCorrelated <- caret::findCorrelation(M, cutoff=(0.7),verbose = FALSE)
important_var_an=colnames(tmp[,-highlyCorrelated])

important_var_an
# [1] "gnomAD_exomes_NFE_AN"          "gnomAD_exomes_controls_ASJ_AN"
# [3] "gnomAD_exomes_controls_NFE_AN" "gnomAD_exomes_controls_SAS_AN"
# [5] "gnomAD_genomes_NFE_AN"         "gnomAD_genomes_SAS_AN"
#####
# Allele nhomalt
#####
Allelenhomalt = c("class",popfreqName[grepl("_nhomalt", popfreqName)])

tmp <- subClin[, colnames(subClin) %in% Allelenhomalt]

# convert . and - to NA
tmp[tmp == "."] <- NA
tmp[tmp == "-"] <- NA

# count missing values
missNumber = data.frame(mis = colSums(is.na(tmp)))
missNumber$percent = missNumber$mis/dim(tmp)[1]

# drop columns with more than 80% missing values
tmp = tmp[, colnames(tmp) %in% rownames(missNumber)[missNumber$percent < 0.8]]


# convert to numeric
tmp[,-1] = sapply(tmp[,-1], as.numeric)

# computing mean for missing values
for(i in 1:ncol(tmp[,-1])){
  tmp[,-1][is.na(tmp[,-1][,i]), i] <- mean(tmp[,-1][,i], na.rm = TRUE)
}


M <-cor(tmp[,-1])


#Keeping only one variable out of pairs which show cor > 0.7
highlyCorrelated <- caret::findCorrelation(M, cutoff=(0.7),verbose = FALSE)
important_var_anhomalt=colnames(tmp[,-highlyCorrelated])

important_var_anhomalt
#####

# combine all to see if there is significant correlation 
tmp = subClin[,colnames(subClin) %in% c(important_var_ac, important_var_af, 
                                        important_var_an, important_var_anhomalt)]

# convert . and - to NA
tmp[tmp == "."] <- NA
tmp[tmp == "-"] <- NA

# convert to numeric
tmp = sapply(tmp, as.numeric)

# computing mean for missing values
for(i in 1:ncol(tmp)){
  tmp[is.na(tmp[,i]), i] <- mean(tmp[,i], na.rm = TRUE)
}


M <-cor(tmp)

library( )
png(filename = "~/clinvar/AS_AC_AN_nhomalt_corr.png", width = 15, height = 15, units = "in", res = 300)
corrplot(M, type="upper", order="hclust",
         col=brewer.pal(n=8, name="RdYlBu"))
dev.off()

# final selection
#Keeping only one variable out of pairs which show cor > 0.8
highlyCorrelated <- caret::findCorrelation(M, cutoff=(0.8),verbose = FALSE)
important_var_AS_AC_AN_nhomalt=colnames(tmp[,-highlyCorrelated])

important_var_AS_AC_AN_nhomalt

# SELECTED FEATURES

# [1] "AF_TGP"                        "AF"                           
# [3] "gnomAD_exomes_controls_NFE_AN" "gnomAD_genomes_AFR_AF"        
# [5] "gnomAD_genomes_NFE_AN" 

Fs = c(AlleleCount[-1], AlleleFreq[-1], AlleleNumber[-1], Allelenhomalt[-1])
toDropFeatures = unique(Fs[-which(Fs %in% important_var_AS_AC_AN_nhomalt)])
# dropping columns from subclin
finClin = subClin[, -which(colnames(subClin) %in% toDropFeatures)]

# More data exolorations:
# Calculating summary stats for allele frequency metrics
library(dplyr)
col = c()
pVal  = c()
mean0 = c()
mean1 = c()
median0 = c()
median1 = c()


for(i in 2:ncol(af)){
  if(all(is.na((as.numeric(af[,i])))) == FALSE){
    col[i] <- colnames(af)[i]
    val = as.numeric(af[,i])
    class = as.factor(af$class)
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

################## for paper
## calculating summary statistics fro groups 

af = tmp[,c(4,5:10, 1,11,12, 13:21, 41:48)]
# removing unwanted pop
af = af[,-21]

af[,-1] = sapply(af[,-1],as.numeric)

varName <- c()
n_conflict= c()
n_noConflict= c()
median_conflict = c()
median_noConflict = c()
q1_Conflict = c()
q3_Conflict = c()
q1_noConflict = c()
q3_noConflict = c()

for (i in 2:ncol(af)){
  varName[i] = colnames(af)[i]
  d = data.frame(cls = af$class, val = af[,i])
  d = d[!is.na(d$val),]
  n_conflict[i] = length(d$val[d$cls == "conflict"])
  n_noConflict[i] = length(d$val[d$cls == "noConflict"])
  median_conflict[i] = median(d$val[d$cls == "conflict"])
  median_noConflict[i] = median(d$val[d$cls == "noConflict"])
  q1_Conflict[i] = quantile(d$val[d$cls == "conflict"])[2]
  q3_Conflict[i] = quantile(d$val[d$cls == "conflict"])[4]
  q1_noConflict[i] = quantile(d$val[d$cls == "noConflict"])[2]
  q3_noConflict[i] = quantile(d$val[d$cls == "noConflict"])[4]
}

resDF <- data.frame(population = varName,
                    n_conflict = n_conflict,
                    n_noConflict= n_noConflict,
                    median_conflict = median_conflict,
                    median_noConflict = median_noConflict,
                    q1_Conflict = q1_Conflict,
                    q3_Conflict = q3_Conflict ,
                    q1_noConflict = q1_noConflict,
                    q3_noConflict =q3_noConflict )

resDF = na.omit(resDF)
resDF = data.frame(Population = resDF$population,
                   conflict_N = resDF$n_conflict,
                   conflictMedianQR = paste0(resDF$median_conflict, " [", resDF$q1_Conflict,"-", resDF$q3_Conflict, "]"),
                   noConflict_N = na.omit(n_noConflict),
                   noConflictMedianQR = paste0(resDF$median_noConflict, " [", resDF$q1_noConflict,"-", resDF$q3_noConflict, "]"))



library(sjPlot)

tab_df(resDF,
       alternate.rows = T, # this colors the rows
       # title = "Population allele frequency", #always give
       digits = 5,
       file = "~/clinvar/the100GP_AF.doc")


# Visualization
# 1 super populations
afSuperPop = af[, c(1,2,8,11,20)]
afSuperPop <- na.omit(reshape2::melt(afSuperPop, id.var = "class"))

tiff(filename = "~/clinvar/paper/SuperPopulationAlellefreq.png", width = 8, height = 5, units = "in", res = 300)
ggplot(afSuperPop, aes(x = variable, y = log10(value))) + 
  geom_boxplot(aes(fill = class), position = position_dodge(0.9)) +
  scale_fill_manual(values = c("#999999", "#E69F00")) + 
  xlab("Population project") + 
  ylab("log-scaled allele frequency") +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=0.5)) +
  scale_x_discrete(labels=c("AF" = "TGP", "AF_ESP" = "ESP",
                              "gnomAD_AF" = "gnomAD", "ExAC_nonTCGA_AF" = "ExAC")) + 
 theme_bw()
dev.off()
# 2 populations in super populations

afSubPop = af[, -c(2,8,11,20)]
afSubPop <- na.omit(reshape2::melt(afSubPop, id.var = "class"))


tiff(filename = "~/clinvar/paper/subpopulationAF_SUPP.png", width = 8, height = 5, units = "in", res = 300)
ggplot(afSubPop, aes(x = variable, y = log10(value))) + 
  geom_boxplot(aes(fill = class), position = position_dodge(0.9)) +
  scale_fill_manual(values = c("#999999", "#E69F00")) + 
  xlab("Population project") + 
  ylab("log-scaled allele frequency") +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=0.5))
  theme_bw()
  dev.off()

  # statitical test
  
  # for TGP
  wilcox.test(AF ~ class, data = af, exact = FALSE)
  
  af_NA_rm = na.omit(af[, c(1,2)])
  wilcox.test(AF ~ class, data = af_NA_rm, exact = FALSE)
  
  # FOR ESP
  wilcox.test(gnomAD_AF ~ class, data = af, exact = FALSE)
  
  
  
##################################################3


########################### Gene Type#########################
# replace "-" in symbol column with ENS gene id from the Gene column
subClin$SYMBOL = ifelse(subClin$SYMBOL == "-", subClin$Gene, subClin$SYMBOL)
  
tmp = data.frame(unclass(table(subClin$SYMBOL, subClin$class)))
tmp$gene = rownames(tmp)
rownames(tmp) = NULL


# visualization
library(ggplot2)
library(ggpubr)
library(dplyr)

plotDf = reshape2::melt(tmp, id.var = "gene")
# top 50 genes
tmp$total = tmp$conflict+tmp$noConflict

# top 50 most sumitted genes
tmp = tmp %>% arrange(desc(total))
top50 = tmp$gene[1:50]

# reducing plotDf
plotDf = plotDf[plotDf$gene %in% top50,]

# Arrange/sort and compute cumulative summs
plotDf <- plotDf %>%
  group_by(gene) %>%
  arrange(desc(value)) %>%
  mutate(lab_ypos = cumsum(value) - 0.5 * value) 

# Create stacked bar graphs with labels
png(filename = "~/clinvar/geneType.png", width = 16, height = 8.135, units = "in", res = 300)
 ggplot(data = plotDf, aes(x = reorder(gene, -value), y = value)) +
  geom_col(aes(fill = variable), width = 0.9)+
  scale_fill_manual(values = c("#999999", "#E69F00")) +
  theme_bw() + 
  theme(axis.text.x=element_text(angle=90)) +
  xlab("Top 50 most submitted genes") + 
  ylab("variant count in each class") +
   theme(
     axis.title.x = element_text(size = 12),
     axis.text.x = element_text(size = 10),
     axis.title.y = element_text(size = 12))
dev.off()


p1 <-  ggplot(data = plotDf, aes(x = reorder(gene, -value), y = value)) +
  geom_col(aes(fill = variable), width = 0.9)+
  scale_fill_manual(values = c("#999999", "#E69F00")) +
  theme_bw() + 
  theme(axis.text.x=element_text(angle=90)) +
  xlab("Top 50 most submitted genes") + 
  ylab("variant count") +
  theme(
    axis.title.x = element_text(size = 12),
    axis.text.x = element_text(size = 10),
    axis.title.y = element_text(size = 12))


# plotting with percent values
gene = tmp$gene
conflict = round(100*tmp$conflict/(tmp$conflict+tmp$noConflict),0)
noConflict = round(100*tmp$noConflict/(tmp$conflict+tmp$noConflict),0)
tmp2 = data.frame(gene = gene, conflict = conflict, noConflict = noConflict)

plotDf = reshape2::melt(tmp2, id.var = "gene")
plotDf = plotDf[plotDf$gene %in% top50,]
plotDf = plotDf[plotDf$variable == "conflict",]
# like what mentioned above
plotDf <- plotDf %>%
  group_by(gene) %>%
  arrange(desc(value)) %>%
  mutate(lab_ypos = cumsum(value) - 0.5 * value) 

# set order
df = tmp2[tmp2$gene %in% top50,] %>% arrange(desc(conflict))
positions <- df$gene

# Create stacked bar graphs with labels
png(filename = "~/clinvar/geneType_percent.png", width = 16, height = 8.135, units = "in", res = 300)
  ggplot(data = plotDf, aes(x = gene, y = value)) +
  geom_col(aes(fill = variable), width = 0.9)+
  geom_text(aes(y = lab_ypos, label = value, group =variable), color = "white") +
  theme_bw() +
  theme(axis.text.x=element_text(angle=90)) +
  scale_fill_manual(values = c(c("#999999", "#E69F00"))) +
  scale_x_discrete(limits = positions) +
  xlab("Top 50 most submitted genes") + 
  ylab("% of variant count in each class") +
  theme(
    axis.title.x = element_text(size = 12),
    axis.text.x = element_text(size = 10),
    axis.title.y = element_text(size = 12))
dev.off

p2 <- ggplot(data = plotDf, aes(x = gene, y = value)) +
  geom_col(aes(fill = variable), width = 0.9)+
  geom_text(aes(y = lab_ypos, label = value, group =variable), color = "white") +
  theme_bw() +
  theme(axis.text.x=element_text(angle=90)) +
  scale_fill_manual(values = c(c("#999999", "#E69F00"))) +
  scale_x_discrete(limits = positions) +
  xlab("Top 50 most submitted genes") + 
  ylab("% of variant count") +
  theme(
    axis.title.x = element_text(size = 12),
    axis.text.x = element_text(size = 10),
    axis.title.y = element_text(size = 12))

# merging two plots together
png(filename = "~/clinvar/paper/geneType.png", width = 16, height = 8.135, units = "in", res = 300)
ggarrange(p1, p2, 
          labels = c("A", "B"),
          ncol = 1, nrow = 2)
dev.off()

## pathway enrichement
hpo_dat = data.table::fread("~/clinvar/sandbox/phenotype_to_genes.txt")
subHPO = aggregate(. ~ V2, hpo_dat[, c(2,4)], FUN = function(x) 
  toString(x), na.action = NULL)

#converting comma to tab seprated value
l = strsplit(as.character(subHPO$V4), ",")
# enrichment analysis
pathway = c()
genes  = list()
pV = numeric()

d = data.frame(matrix(nrow = length(l), ncol = 3))

for(i in 1:length(l)){
  d[i,1] = subHPO$V2[i]
  group1 = top50
  group2 = gsub(" ", "", l[[i]])
  overlap = length(intersect(group1, group2))
  if(overlap > 0){
  genes[[i]] = paste(intersect(group1, group2))
  total = 21306 # human protein coding gene 
  #phyper(overlap-1, length(group2), total-length(group2), length(group1),lower.tail= FALSE)
  d[i,3] = fisher.test(matrix(c(overlap, length(group2)-overlap, 
                       length(group1)-overlap, total-length(group2)-length(group1) +overlap), 2, 2), alternative='greater')$p.value
  print(i)
  }
}

for(i in 1:nrow(d)){
d$X2[i] <- list(genes[[i]])
print(i)
}

d = d[d$X2 != "NULL",]
d$bonferroniP = d$X3

d = d[d$bonferroniP <= 0.05,]

#
pkhd1 = d[grep("PKHD1", d$X2),]



f = numeric()
for(i in 1:length(l)){
  f = c(f, length(l[[i]]))
}
# determine the df dimention
ll = list()
for(i in 1:length(l)){
  X <- vector(mode="character", length=max(f))
  X[1:length(l[[i]])] = l[[i]]
  ll[[i]] = X
  print(paste(100*i/length(l), "%"))
}

d = do.call(rbind, ll)

hpo = data.frame(cbind(subHPO[,1], d))

#removing terms greater than 500 gene
count <- apply(hpo, 1, function(x) length(which(x !="")))
#
subHPO = hpo[count < 500 & count > 3,]
subHPO[is.na(subHPO)] <- NA
write.table(subHPO, "~/clinvar/sandbox/hpo.gmt", col.names = F, row.names = F, quote = F, sep = "\t")

# saving gmt 
hpo_geneset <- fgsea::gmtPathways("~/clinvar/sandbox/hpo.gmt")

# making gene set
meltedHPO = reshape2::melt(subHPO, id.var = "X1")
meltedHPO = meltedHPO[meltedHPO$value != "",]
meltedHPO = meltedHPO[, c(1,3)]
names(meltedHPO) = c("Pathway", "SYMBOL")
# converting symbols to ENS ids
sym2ens <- AnnotationDbi::select(org.Hs.eg.db,
                                 key=meltedHPO$SYMBOL, 
                                 columns="ENSEMBL",
                                 keytype="SYMBOL")

# final dataset preparation
finClin$SYMBOL =  subClin$SYMBOL 

########################allele type#############################################

tmp = subClin[, c(30,4,5)]
tmp$allele = paste0(tmp$REF, ">", tmp$ALT)

tmp$REFlenCat = ifelse(nchar(tmp$REF) ==1, "1", ">1")
tmp$ALTlenCat = ifelse(nchar(tmp$ALT) ==1, "1", ">1")
# coding REF multi_nts and ALT multi_nts
tmp$allele = ifelse(tmp$REFlenCat != 1, "REFmulti_nts", 
                    ifelse(tmp$ALTlenCat != 1,"ALTmulti_nts", tmp$allele))

# visualization REF type
tmp = data.frame(unclass(table(tmp$allele, tmp$class)))
tmp$allele = rownames(tmp)
tmp$total = tmp$noConflict+tmp$conflict
tmp = tmp %>% arrange(desc(allele))

# adding total number to group name
tmp$allele = paste0(tmp$allele, " [", tmp$total, "]")

# plotting with percent values
tmp$conflict = tmp$conflict/tmp$total
tmp$noConflict = tmp$noConflict/tmp$total

# rounding and *100
tmp$conflict = round(tmp$conflict * 100,0)
tmp$noConflict = round(tmp$noConflict * 100,0)

df = tmp[order(tmp$conflict, decreasing = T),]
positions <- df$allele

plotDf = reshape2::melt(tmp[-4], id.var = "allele")

# like what mentioned above
plotDf <- plotDf %>%
  group_by(allele) %>%
  arrange(desc(value)) %>%
  mutate(lab_ypos = cumsum(value) - 0.5 * value) 

#vis

png(filename = "~/clinvar/alleleType_percent.png", width = 16, height = 8.135, units = "in", res = 300)
ggplot(data = plotDf, aes(x = allele, y = value)) +
  geom_col(aes(fill = variable), width = 0.9)+
  geom_text(aes(y = lab_ypos, label = value, group =variable), color = "white") +
  theme_bw() +
  theme(axis.text.x=element_text(angle=45, vjust =0.9, hjust = 1)) +
  scale_fill_manual(values = c("#999999", "#E69F00")) +
  scale_x_discrete(limits = positions) +
  xlab("Allele change [total number of events]") + 
  ylab("% of variant count in each class") +
  theme(
    axis.title.x = element_text(size = 14),
    axis.text.x = element_text(size = 12),
    axis.title.y = element_text(size = 14))
dev.off()

# create mosaic plot to see association
tmp = subClin[, c(30,4,5)]
tmp$allele = paste0(tmp$REF, ">", tmp$ALT)

tmp$REFlenCat = ifelse(nchar(tmp$REF) ==1, "1", ">1")
tmp$ALTlenCat = ifelse(nchar(tmp$ALT) ==1, "1", ">1")
# coding REF multi_nts and ALT multi_nts
tmp$allele = ifelse(tmp$REFlenCat != 1, "REFmulti_nts", 
                    ifelse(tmp$ALTlenCat != 1,"ALTmulti_nts", tmp$allele))

# creating data table
dt = as.table(as.matrix(table(tmp$allele, tmp$class)))

png(filename = "~/clinvar/mosaic_allele_change.png", width = 16, height = 8.135, units = "in", res = 300)
mosaicplot(dt, shade = T, las=2,
           main = "Allele change", cex.axis = 1.5)
dev.off()



# final dataset preparation:
# dropping column 5 and 4
names(subClin)[4:5]
# [1] "REF" "ALT"
names(finClin)[4:5]
# [1] "REF" "ALT"

finClin = finClin[, -c(4,5)]

# adding new allele column as allele column
finClin$Allele = tmp$allele

# up to now : SYMBOL, Allele
######################################## Molecular type ###########################################

names(subClin)


table(subClin$class, subClin$CLNVC)

#                  Deletion Duplication  Indel    Insertion  Inversion  Microsatellite              SNP        Variation
# noConfilict     5176        2362       541       229        50           1873                    128380         6
# conflict        758         431        108        51        13            565                     47100         0


tmp = data.frame(unclass(table(subClin$class, subClin$CLNVC)))
tmp = data.frame(t(tmp))
tmp$type = rownames(tmp)
tmp$total = tmp$noConflict+tmp$conflict
tmp = tmp %>% arrange(desc(total))
tmp$type[1] = "SNV"

# adding total number to group name
tmp$type = paste0(tmp$type, " [", tmp$total, "]")

# plotting with percent values
tmp$conflict = tmp$conflict/tmp$total
tmp$noConflict = tmp$noConflict/tmp$total

# rounding and *100
tmp$conflict = round(tmp$conflict * 100,0)
tmp$noConflict = round(tmp$noConflict * 100,0)



plotDf = reshape2::melt(tmp[-4], id.var = "type")

# like what mentioned above
plotDf <- plotDf %>%
  group_by(type) %>%
  arrange(desc(value)) %>%
  mutate(lab_ypos = cumsum(value) - 0.5 * value) 

#vis
positions = c("SNV [175340]","Deletion [5931]","Duplication [2793]","Microsatellite [2438]","Indel [648]","Insertion [280]","Inversion [63]")

png(filename = "~/clinvar/Type_percent.png", width = 16, height = 8.135, units = "in", res = 300)
ggplot(data = plotDf, aes(x = type, y = value)) +
  geom_col(aes(fill = variable), width = 0.6)+
  geom_text(aes(y = lab_ypos, label = value, group =variable), color = "white") +
  theme_bw() +
  theme(axis.text.x=element_text(angle=45, vjust =0.9, hjust = 1)) +
  scale_fill_manual(values = c(c("#999999", "#E69F00"))) +
  scale_x_discrete(limits = positions) +
  xlab("Molecular type of variant [number of mutations]") + 
  ylab("% of variant count in each class") +
  theme(
    axis.title.x = element_text(size = 14),
    axis.text.x = element_text(size = 12),
    axis.title.y = element_text(size = 14))
dev.off()


#Final dataset preparation

# up to now : "SYMBOL", "Allele", "CLNVC"

##################################Consequences############################################

# 
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
#rownames(tmp) <- c("noConflict", "conflict")
tmp = data.frame(t(tmp))
tmp$conseq = rownames(tmp)
tmp$total = tmp$noConflict+tmp$conflict
tmp = tmp %>% arrange(desc(total))


# selection this category with data larger than 1K:
tmp = tmp[tmp$total > 1000,]

# melting dataframe
plotDf = reshape2::melt(tmp[-4], id.var = "conseq")

# correcting names
plotDf$conseq[plotDf$conseq == "X3_prime_UTR_variant"] <- "3_prime_UTR_variant"
plotDf$conseq[plotDf$conseq == "X5_prime_UTR_variant"] <- "5_prime_UTR_variant"
plotDf$value2 = log10(plotDf$value)

# plotting
png(filename = "~/clinvar/conseqBarplot.png", width = 16, height = 8.135, units = "in", res = 300)
p1 = ggplot(data = plotDf, aes(x = reorder(conseq, -value), y = value)) +
  geom_col(aes(fill = variable), width = 0.9)+
  scale_fill_manual(values = c("#999999", "#E69F00")) +
  theme_bw() + 
  theme(axis.text.x=element_text(angle=90)) +
  xlab("Consequence groups with >1K variants") + 
  ylab("variant count in each class") +
  theme(
    axis.title.x = element_text(size = 14),
    axis.text.x = element_text(size = 12),
    axis.title.y = element_text(size = 14))
dev.off()

#
# plotting with percent values
tmp$conflict = tmp$conflict/tmp$total
tmp$noConflict = tmp$noConflict/tmp$total

# rounding and *100
tmp$conflict = round(tmp$conflict * 100,0)
tmp$noConflict = round(tmp$noConflict * 100,0)
# renaming variables
tmp$conseq[tmp$conseq == "X3_prime_UTR_variant"] <- "3_prime_UTR_variant"
tmp$conseq[tmp$conseq == "X5_prime_UTR_variant"] <- "5_prime_UTR_variant"


# adding total number to group name
tmp$conseq = paste0(tmp$conseq, " [", tmp$total, "]")

plotDf = reshape2::melt(tmp[-4], id.var = "conseq")

# like what mentioned above
plotDf <- plotDf %>%
  group_by(conseq) %>%
  arrange(desc(value)) %>%
  mutate(lab_ypos = cumsum(value) - 0.5 * value) 

# error in 3' prime number position
plotDf$lab_ypos[plotDf$conseq == '3_prime_UTR_variant [1559]' & plotDf$variable == 'noConflict'] = 32
plotDf$lab_ypos[plotDf$conseq == '3_prime_UTR_variant [1559]' & plotDf$variable == 'conflict'] = 82
  

# vis
positions = tmp$conseq

png(filename = "~/clinvar/Conseq_percent.png", width = 16, height = 8.135, units = "in", res = 300)
ggplot(data = plotDf, aes(x = conseq, y = value)) +
  geom_col(aes(fill = variable), width = 0.9)+
  geom_text(aes(y = lab_ypos, label = value, group =variable), color = "white") +
  theme_bw() +
  theme(axis.text.x=element_text(angle=60, hjust = 1.1)) +
  scale_fill_manual(values = c(c("#999999", "#E69F00"))) +
  scale_x_discrete(limits = positions) +
  xlab("Consequence groups with >1K variants") + 
  ylab("% of variant count in each class") +
  theme(
    axis.title.x = element_text(size = 14),
    axis.text.x = element_text(size = 12),
    axis.title.y = element_text(size = 14))
dev.off()



#Final dataset preparation

# up to now : "SYMBOL", "Allele", "CLNVC", "Consequence"

#########################################   IMPACT   ######################################
 
table(subClin$class, subClin$IMPACT)

tmp = data.frame(unclass(table(subClin$class, subClin$IMPACT)))
tmp = data.frame(t(tmp))
tmp$impact = rownames(tmp)
tmp$total = tmp$noConflict+tmp$conflict
tmp = tmp %>% arrange(desc(total))

# adding total number to group name
tmp$impact = paste0(tmp$impact, " [", tmp$total, "]")

# visualization of chi-square residuals
rownames(tmp) = tmp$impact
tmp = tmp[, c(1,2)]

chisq <- chisq.test(t(tmp))
chisq

library(corrplot)
corrplot(chisq$residuals, is.cor = FALSE,
         tl.cex = 1,
         tl.col = "black")


# plotting with percent values
tmp$conflict = tmp$conflict/tmp$total
tmp$noConflict = tmp$noConflict/tmp$total

# rounding and *100
tmp$conflict = round(tmp$conflict * 100,0)
tmp$noConflict = round(tmp$noConflict * 100,0)



plotDf = reshape2::melt(tmp[-4], id.var = "impact")

# like what mentioned above
plotDf <- plotDf %>%
  group_by(impact) %>%
  arrange(desc(value)) %>%
  mutate(lab_ypos = cumsum(value) - 0.5 * value) 

#vis

positions = c("MODIFIER [24672]", "LOW [58738]", "MODERATE [88013]", "HIGH [16070]")

png(filename = "~/clinvar/IMPACT_percent.png", width = 16, height = 8.135, units = "in", res = 300)
ggplot(data = plotDf, aes(x = impact, y = value)) +
  geom_col(aes(fill = variable), width = 0.6)+
  geom_text(aes(y = lab_ypos, label = value, group =variable), color = "white") +
  theme_bw() +
  theme(axis.text.x=element_text(angle=45, vjust =0.9, hjust = 1)) +
  scale_fill_manual(values = c(c("#999999", "#E69F00"))) +
  scale_x_discrete(limits = positions) +
  xlab("VEP impact group [number of mutations]") + 
  ylab("% of variant count in each class") +
  theme(
    axis.title.x = element_text(size = 14),
    axis.text.x = element_text(size = 12),
    axis.title.y = element_text(size = 14))
dev.off()


#

#Final dataset preparation

# up to now : "SYMBOL", "Allele", "CLNVC", "Consequence", "IMPACT"

###################################transcript length #################################################
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


library(ggplot2)
png(filename = "~/clinvar/transcript_length.png", width = 7, height = 10, units = "in", res = 300)
ggplot(tmp, aes(x = as.factor(class), y = log10(transcript_length))) + 
  geom_boxplot(aes(fill = class), position = position_dodge(0.9)) +
  scale_fill_manual(values = c("#999999", "#E69F00")) + 
  theme_bw()
dev.off()

#
finClin$transcript_length  = tmp$transcript_length
#Final dataset preparation

# up to now : "SYMBOL", "Allele", "CLNVC", "Consequence", "IMPACT", "transcript_length"
######################### variant location in exons ##########################################
View(table(subClin$EXON))

table(subClin$EXON == "-")

# FALSE   TRUE 
# 153639  33854

tmp = subClin[, c(30, 40, 41)]
tmp = tmp[!is.na(tmp$EXON),]


tmp$ratio = tmp$EXON
tmp$ratio = ifelse(tmp$ratio != "-", tmp$ratio, tmp$INTRON)
tmp = cbind(tmp, stringr::str_split_fixed(tmp$ratio, "/", 2))
names(tmp)[c(5:6)] <- c("nExon", "totalExon")
tmp$nratio = as.numeric(tmp$nExon)/as.numeric(tmp$totalExon)

#
tmp_aggregate = tmp

t.test(tmp$nratio ~ tmp$class, data = tmp, var.equal = TRUE)

# ridge plot
png(filename = "~/clinvar/paper/ExonIntron_number_ridgeplot.png", width = 16, height = 8.135, units = "in", res = 300)
p1 = ggplot(tmp, aes(x = nratio, y = class)) +
  geom_density_ridges(aes(fill = class),) +
  scale_fill_manual(values = c("#999999", "#E69F00")) +
  xlab("Relative location of variant in transcript") + 
  ylab("") +
  scale_x_continuous(n.breaks = 5) +
  ggtitle("")
dev.off()
# those variants with "-" as Exon and Intron will be removed since they generate NA
# those variants encompass more than one Exon/Intron will be removed since they generate NA
#colSums(is.na(tmp))
# class   EXON INTRON  ratio  nExon  tExon nratio 
# 0       0      0      0      0      0     7698
#tmp = tmp[!is.na(tmp$nratio),]
# coding 0 for introns
tmp$nIntron = tmp$nExon
tmp$nIntron = as.factor(ifelse(tmp$INTRON == "-", tmp$nExon, "0"))

plotDf= data.frame(unclass(table(tmp$nIntron, tmp$class)))
plotDf$ExInt = as.numeric(rownames(plotDf))
plotDf =plotDf %>% arrange(ExInt)

plotDf = plotDf[c(1:20),]

plotDf = reshape2::melt(plotDf, id.var = "ExInt")
plotDf$ExInt = as.character(plotDf$ExInt)
#visualization
#png(filename = "~/clinvar/ExonIntron_number.png", width = 16, height = 8.135, units = "in", res = 300)
ggplot(data = plotDf, aes(x = reorder(ExInt, -value), y = value)) +
  geom_col(aes(fill = variable), width = 0.9)+
  scale_fill_manual(values = c("#999999", "#E69F00")) +
  theme_bw() + 
  theme(axis.text.x=element_text(angle=90)) +
  xlab("First 20 exons") + 
  ylab("variant count in each class") +
  theme(
    axis.title.x = element_text(size = 12),
    axis.text.x = element_text(size = 10),
    axis.title.y = element_text(size = 12))
#dev.off()



# discertization visualization
## Import libraries
library(scorecard)
library(ggplot2)
library(hrbrthemes)

tt = tmp[, c('nratio', 'class')]
tt = tt[!is.na(tt$nratio),]
# read more on discertization : https://nextjournal.com/eda/discretize-cont-var
## Calculate bin breaks for numeric variables with respect to their relationships with the outcome variable class
bins = scorecard::woebin(tt,
                         y = 'class',
                         positive = 'conflict',
                         )

bins2_tree = woebin(tt, 
                    y="class",
                    x=c("nratio"), 
                    method="tree", 
                    positive = 'conflict',
                    breaks_list = list(nratio = c(0,0.2,0.4,0.6, 0.8)))
bins2_tree




# visualization
scorecard::woebin_plot(bins2_tree$nratio,
                       #line_color='#FC8D59',
                       bar_color=c("#999999", "#E69F00"),
                       show_iv = F)
# interpertaion for the graph: highest probability for conflicting group is found in bin 0.96 and greater

# density plot
ggplot(data=tmp, aes(x=nratio, group=class, fill=class)) +
  geom_density(adjust=1.5, position="fill") +
  scale_fill_manual(values = c("#999999", "#E69F00")) +
  scale_x_continuous(breaks = seq(from = 0, to = 1, by = 0.2)) +
  xlab("Relative position of a variant in the transcript") + 
  ylab("Frequency of mutations from diffrent classes") +
  theme_ipsum()

# calculate significance
tmp$relLoc = ifelse(tmp$nratio <= 0.2, "0.0-0.2", 
                    ifelse(tmp$nratio > 0.2 & tmp$nratio <=0.4, "0.2-0.4",
                           ifelse(tmp$nratio > 0.4 & tmp$nratio <=0.6, "0.4-0.6",
                                  ifelse(tmp$nratio > 0.6 & tmp$nratio <=0.8, "0.6-0.8",
                                         ifelse(tmp$nratio > 0.8 & tmp$nratio <=1.0, "0.8-1.0", NA)))))


table(tmp$relLoc)

table(tmp$relLoc, tmp$class)

chisq.posthoc.test::chisq.posthoc.test(table(tmp$relLoc, tmp$relLoc))

dt <- as.table(table(tmp$relLoc, tmp$class))


mosaicplot(dt, shade = TRUE,main = "B")
my_base <- recordPlot() 
# combining plots
library("ggplotify")
library("patchwork")
library(ggpubr)
theme_set(theme_pubr())

p2 = as.ggplot(~mosaicplot(dt, shade = TRUE,main = ""))

png(filename = "~/clinvar/density_mosaic.png", width = 16, height = 8.135, units = "in", res = 300)
ggarrange(p1, p2,
          labels = c("A", "B"),
          ncol = 2, nrow = 1,
          widths = c(2:1))
dev.off()
  
# proportion table by row
prop.table(table(tmp$relLoc, tmp$class), 1)

# proportion table by column
prop.table(table(tmp$relLoc, tmp$class), 2)

# to calculate OR

tmp$range1 = ifelse(tmp$relLoc == "0.0-0.2","0.0-0.2", "other" )
tmp$range2 = ifelse(tmp$relLoc == "0.2-0.4","0.2-0.4", "other" )
tmp$range3 = ifelse(tmp$relLoc == "0.4-0.6","0.4-0.6", "other" )
tmp$range4 = ifelse(tmp$relLoc == "0.6-0.8","0.6-0.8", "other" )
tmp$range5 = ifelse(tmp$relLoc == "0.8-1.0","0.8-1.0", "other" )

# calculate OR
abd::odds.ratio(table(tmp$class, tmp$range3))

######other position based annotation
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
# so nratio and protein position will be selected

# two variable need to be added to finClin
finClin$relativeLocationRatio = tmp$nratio
finClin$ExIntron = tmp$nIntron

#Final dataset preparation

# up to now : "SYMBOL", "Allele", "CLNVC", "Consequence", "IMPACT", "relativeLocationRatio",

#"ExIntron", "Protein_position"

# Adding data on rebviewd_by expert panel variants
df <- data.table::fread("~/clinvar/clinvar_tsv.txt")

table(df$CLNREVSTAT)

df = df[df$CLNREVSTAT == "reviewed_by_expert_panel",]

# vep output
vep = data.table::fread("~/clinvar/sandbox/review_by_panel_vep.txt")
#vep = data.table::fread("~/clinvar/sandbox/freq_var_vep.txt")

# overlap
table(df$ID %in% vep$`#Uploaded_variation`)

table(vep$`#Uploaded_variation` %in% df$ID)

#
table(duplicated(vep$`#Uploaded_variation`))


# reduce vep to only contain data present in the dataset
vep = vep[vep$`#Uploaded_variation` %in% df$ID,]
# calculating RLV

View(table(vep$EXON))

vep$class = "expert_panel_reviewed"
#vep$class = "frequent_variants"


tmp = vep[, c(44, 11, 12)]
tmp = tmp[!is.na(tmp$EXON),]


tmp$ratio = tmp$EXON
tmp$ratio = ifelse(tmp$ratio != "-", tmp$ratio, tmp$INTRON)
tmp = cbind(tmp, stringr::str_split_fixed(tmp$ratio, "/", 2))
names(tmp)[c(5:6)] <- c("nExon", "totalExon")
tmp$nratio = as.numeric(tmp$nExon)/as.numeric(tmp$totalExon)

tmp_expert = tmp
#tmp_freq = tmp
tmp = rbind(tmp_aggregate, tmp_expert)

# tmp for common variants and rare variants
v1 = data.table::fread("~/clinvar/sandbox/common_20K_vep.txt")
v2 = data.table::fread("~/clinvar/sandbox/rare_20K _vep.txt")
v1$class = "frequent_variants"
v2$class = "rare_variants"

#
t1 = v1[, c(43, 11, 12)]
t1 = t1[!is.na(t1$EXON),]

t2 = v2[, c(43, 11, 12)]
t2 = t2[!is.na(t2$EXON),]


t1$ratio = t1$EXON
t2$ratio = t2$EXON

t1$ratio = ifelse(t1$ratio != "-", t1$ratio, t1$INTRON)
t1 = cbind(t1, stringr::str_split_fixed(t1$ratio, "/", 2))
names(t1)[c(5:6)] <- c("nExon", "totalExon")
t1$nratio = as.numeric(t1$nExon)/as.numeric(t1$totalExon)
#
t2$ratio = ifelse(t2$ratio != "-", t2$ratio, t2$INTRON)
t2 = cbind(t2, stringr::str_split_fixed(t2$ratio, "/", 2))
names(t2)[c(5:6)] <- c("nExon", "totalExon")
t2$nratio = as.numeric(t2$nExon)/as.numeric(t2$totalExon)

t = rbind(t1,t2)


tmp = rbind(tmp_aggregate, tmp_expert, t)

tmp$class = factor(tmp$class, levels = c("conflict", "noConflict", 
                                         "expert_panel_reviewed", 
                                         "frequent_variants", "rare_variants"))


# ridge plot
png(filename = "~/clinvar/paper/RLV_five_group.png", width = 16, height = 8.135, units = "in", res = 300)
ggplot(tmp, aes(x = nratio, y = class)) +
  geom_density_ridges(aes(fill = class),) +
  scale_fill_manual(values = c("#999999", "#E69F00", "#8A9045FF", '#00A1D5FF', '#B24745FF'),
                    labels = c("conflict (n = 49019)", "noConflict (n = 138474)",
                               "expert_panel_reviewed (n = 11528)", "frequent_variants (n= 20994)",
                               "rare_variants (n = 20808)")) + 
  xlab("Relative location of variant in transcript") + 
  ylab("") +
  scale_x_continuous(n.breaks = 5) +
  ggtitle("")
dev.off()




# due to time out error, a combination of UCSC and Galaxy was used to retrieve
# dbsnp id and maf



# then on the CC:
#/home/ghaedi/projects/def-gooding-ab/ghaedi/clinvar/
#cat <(grep '^#' clinvar.vcf) <(grep 'reviewed_by_expert_panel' clinvar.vcf) > review_by_expert_panel.vcf



############################### codon type ############################################
d = data.frame(table(subClin$Codons))
d = d[d$Var1 != "-",]
d = d[with(d, order(-Freq)), ]
rownames(d) <- NULL
top60 = d$Var1[1:70]

#
tmp = data.frame(unclass(table(subClin$Codons, subClin$class)))
tmp$codon = rownames(tmp)
rownames(tmp) = NULL

# count missing value:
tmp$conflict[tmp$codon == "-"] + tmp$noConflict[tmp$codon == "-"]
#[1] 37378

# dropping missing value
tmp = tmp[tmp$codon != "-",]


plotDf = reshape2::melt(tmp, id.var = "codon")

# reducing plotDf
plotDf = plotDf[plotDf$codon %in% top60,]


# 
df = plotDf
df$codon1 = toupper(substr(plotDf$codon,1,3))
df$codon2 = toupper(substr(plotDf$codon,5,7))

codonVec = data.frame(codon = c("TTT", "TTC", "TTA", "TTG", "CTT", "CTC", "CTA", "CTG", "ATT", "ATC", "ATA", "ATG", "GTT", "GTC", "GTA", "GTG", "TGT", "TGC", "TGA", "TGG", "CGT", "CGC", "CGA", "CGG", "AGT", "AGC", "AGA", "AGG", "GGT", "GGC", "GGA", "GGG", "TAT", "TAC", "TAA", "TAG", "CAT", "CAC", "CAA", "CAG", "AAT", "AAC", "AAA", "AAG", "GAT", "GAC", "GAA", "GAG", "TCT", "TCC", "TCA", "TCG", "CCT", "CCC", "CCA", "CCG", "ACT", "ACC", "ACA", "ACG", "GCT", "GCC", "GCA", "GCG"),
                      aa = c("Phe", "Phe", "Leu", "Leu", "Leu", "Leu", "Leu", "Leu", "Ile", "Ile", "Ile", "Met", "Val", "Val", "Val", "Val", "Cys", "Cys", "Ter", "Trp", "Arg", "Arg", "Arg", "Arg", "Ser", "Ser", "Arg", "Arg", "Gly", "Gly", "Gly", "Gly", "Tyr", "Tyr", "Ter", "Ter", "His", "His", "Gln", "Gln", "Asn", "Asn", "Lys", "Lys", "Asp", "Asp", "Glu", "Glu", "Ser", "Ser", "Ser", "Ser", "Pro", "Pro", "Pro", "Pro", "Thr", "Thr", "Thr", "Thr", "Ala", "Ala", "Ala", "Ala")
)

vec = codonVec$aa
names(vec) = codonVec$codon

# adding aminon acid to the tables
df$aa1 = stringr::str_replace_all(df$codon1, vec)
df$aa2 = stringr::str_replace_all(df$codon2, vec)
df$silent = ifelse(df$aa1 == df$aa2, TRUE, FALSE)

# Create stacked bar graphs with labels
png(filename = "~/clinvar/codonType.png", width = 16, height = 8.135, units = "in", res = 300)
p1 = ggplot(data = plotDf, aes(x = reorder(codon, -value), y = value)) +
  geom_col(aes(fill = variable), width = 0.9)+
  scale_fill_manual(values = c("#999999", "#E69F00")) +
  theme_bw() + 
  theme(axis.text.x=element_text(angle=90)) +
  xlab("Top 70 frequent codon types") + 
  ylab("variant count in each class") +
  theme(
    axis.title.x = element_text(size = 12),
    axis.text.x = element_text(size = 10),
    axis.title.y = element_text(size = 12))
dev.off()

# plotting with percent values
codon = tmp$codon
conflict = round(100*tmp$conflict/(tmp$conflict+tmp$noConflict),0)
noConflict = round(100*tmp$noConflict/(tmp$conflict+tmp$noConflict),0)
tmp2 = data.frame(codon = codon, conflict = conflict, noConflict = noConflict)

plotDf = reshape2::melt(tmp2, id.var = "codon")
plotDf = plotDf[plotDf$codon %in% top60,]
# like what mentioned above
plotDf <- plotDf %>%
  group_by(codon) %>%
  arrange(desc(value)) %>%
  mutate(lab_ypos = cumsum(value) - 0.5 * value) 

# set order
df = tmp2[tmp2$codon %in% top60,] %>% arrange(desc(conflict))
positions <- df$codon

# Create stacked bar graphs with labels
png(filename = "~/clinvar/codonType_percent.png", width = 16, height = 8.135, units = "in", res = 300)
p2 = ggplot(data = plotDf, aes(x = codon, y = value)) +
  geom_col(aes(fill = variable), width = 0.9)+
  geom_text(aes(y = lab_ypos, label = value, group =variable), color = "white") +
  theme_bw() +
  theme(axis.text.x=element_text(angle=90)) +
  scale_fill_manual(values = c(c("#999999", "#E69F00"))) +
  scale_x_discrete(limits = positions) +
  xlab("Top 70 frequent codon types") + 
  ylab("% of variant count in each class") +
  theme(
    axis.title.x = element_text(size = 12),
    axis.text.x = element_text(size = 10),
    axis.title.y = element_text(size = 12))
dev.off()

png(filename = "~/clinvar/paper/CodonType.png", width = 16, height = 12, units = "in", res = 300)
ggarrange(p1, p2, 
          labels = c("A", "B"),
          ncol = 1, nrow = 2)
dev.off()

# codons less frequent than 100 times, recoded as other:

codon = subClin$Codons
codon100 = d$Var1[d$Freq > 100]
codonLessthan100 = d$Var1[d$Freq < 100]
codon = ifelse(codon %in% codon100, codon, 
               ifelse(codon %in% codonLessthan100, "other", NA))

finClin$Codons = codon
#Final dataset preparation

# up to now : "SYMBOL", "Allele", "CLNVC", "Consequence", "IMPACT", "relativeLocationRatio",

#"ExIntron", "Protein_position", "Codons"

############################### Transcript Support Level ############################################
# The Transcript Support Level (TSL) is a method to highlight the 
# well-supported and poorly-supported transcript models for users. 
# The method relies on the primary data that can support full-length transcript 
# structure: mRNA and EST alignments supplied by UCSC and Ensembl.
# tsl1 - all splice junctions of the transcript are supported by at least one non-suspect mRNA
# tsl2 - the best supporting mRNA is flagged as suspect or the support is from multiple ESTs
# tsl3 - the only support is from a single EST
# tsl4 - the best supporting EST is flagged as suspect
# tsl5 - no single transcript supports the model structure

table(subClin$TSL)
# -      1         2      3      4      5 
# 16101  136935   3782   2947    120  27608

d = data.frame(table(subClin$TSL))
d = d[d$Var1 != "-",]
d = d[with(d, order(-Freq)), ]
rownames(d) <- NULL

#
tmp = data.frame(unclass(table(subClin$TSL, subClin$class)))
tmp$TSL = rownames(tmp)
rownames(tmp) = NULL

tmp$pconf = round(tmp$conflict/ (tmp$conflict + tmp$noConflict),2)
tmp$pNconf = round(tmp$noConflict/ (tmp$conflict + tmp$noConflict),2)

tmp$conflict = paste0(tmp$conflict, " (", tmp$pconf, ")")
tmp$noConflict = paste0(tmp$noConflict, " (", tmp$pNconf, ")")


# writing to a word document
tab_df(tmp,
       alternate.rows = T, # this colors the rows
       # title = "Population allele frequency", #always give
       digits = 5,
       file = "~/clinvar/paper/TSL_table.doc")


# # count missing value:
# tmp$conflict[tmp$TSL == "-"] + tmp$noConflict[tmp$TSL == "-"]
# #[1] 16101
# 
# # dropping missing value
# tmp = tmp[tmp$TSL != "-",]
# 
# 
# 
# 
# plotDf = reshape2::melt(tmp, id.var = "TSL")
# plotDf = plotDf[plotDf$TSL != "-",]
# plotDf$TSL = paste0("TSL_", plotDf$TSL)
# 
# # Arrange/sort and compute cumulative summs
# plotDf <- plotDf %>%
#   group_by(TSL) %>%
#   arrange(desc(value)) %>%
#   mutate(lab_ypos = cumsum(value) - 0.5 * value) 
# 
# # calculating value bands (n = 5)
# quantile(plotDf$value, c(.20, .40, .60, .80))
# #   20%     40%      60%     80% 
# #   243.0   2006.4   4935.2  22620.6 
# 
# plotDf$band = ifelse(plotDf$value < 250, "<250",
#                      ifelse(plotDf$value >250 & plotDf$value < 2000, "250-2K",
#                             ifelse(plotDf$value >2000 & plotDf$value < 5000, "2K-5K",
#                                    ifelse(plotDf$value > 5000 & plotDf$value < 20000, "5k-20K", ">20K"))))
# 
# # made number as minus
# plotDf$log_value = log10(plotDf$value)
# plotDf$value[plotDf$variable == "conflict"] = -(plotDf$value[plotDf$variable == "conflict"])
# plotDf$log_value[plotDf$variable == "conflict"] = -(plotDf$log_value[plotDf$variable == "conflict"])
# 
# 
# #Plotting
# p1 <- plotDf %>% 
#   ggplot(aes(TSL, log_value)) + 
#   facet_wrap(~ variable, scales = "free_x") + 
#   geom_col(aes(fill = variable)) +
#   scale_fill_manual(values = c("#999999", "#E69F00")) +
#   coord_flip() +
#   theme(axis.text.y = element_text(hjust = 0.5, margin = margin(0, 0, 0, 0)),
#         axis.ticks.length = unit(0, "pt")) +
#   scale_y_continuous(labels = function(x) paste0(abs(x) * 100, "%")) +
#   labs(x = NULL)
# 
# p1_g <- ggplotGrob(p1)
# 
# p1_g$widths[7] <- p1_g$widths[4] + unit(1, "cm")
# 
# p1g_axis <- gtable_filter(p1_g, "axis-l-1-1") 
# 
# p1_g %>% 
#   gtable_add_grob(p1g_axis, l = 7, t = 8, name = "middle_axis") %>% # add the axis to the middle
#   delete_col("axis-l-1-1") %>% # delete the original axis
#   gtable_add_grob(textGrob("Label", gp = gpar(fontsize = 11)), l = 7, t = 7) %>% # add the top label
#   grid.draw() # draw the result
# 
# 
# ggplot(plotDf, aes(TSL, value)) + 
#   facet_wrap(~ variable, scales = "free_x") + 
#   geom_col(aes(fill = variable)) +
#   scale_fill_manual(values = c("#999999", "#E69F00")) +
#   coord_flip() +
#   scale_y_continuous(expand = c(0, 0)) +
#   theme(panel.spacing.x = unit(0, "mm")) +
#   scale_y_continuous(
#     expand = c(0, 0), 
#     labels = function(x) signif(abs(x), 3)
#   ) +
#   theme_bw()
# 
# 
# 
# # Create stacked bar graphs with labels
# png(filename = "~/clinvar/TSLType.png", width = 16, height = 8.135, units = "in", res = 300)
# ggplot(data = plotDf, aes(x = reorder(TSL, -value), y = value)) +
#   geom_col(aes(fill = variable), width = 0.9)+
#   scale_fill_manual(values = c("#999999", "#E69F00")) +
#   theme_bw() + 
#   theme(axis.text.x=element_text(angle=90)) +
#   xlab("TSL types") + 
#   ylab("variant count in each class") +
#   theme(
#     axis.title.x = element_text(size = 12),
#     axis.text.x = element_text(size = 10),
#     axis.title.y = element_text(size = 12))
# dev.off()
# 
# # plotting with percent values
# TSL = tmp$TSL
# conflict = round(100*tmp$conflict/(tmp$conflict+tmp$noConflict),0)
# noConflict = round(100*tmp$noConflict/(tmp$conflict+tmp$noConflict),0)
# tmp2 = data.frame(TSL = TSL, conflict = conflict, noConflict = noConflict)
# 
# plotDf = reshape2::melt(tmp2, id.var = "TSL")
# plotDf = plotDf[plotDf$TSL != "-",]
# 
# # like what mentioned above
# plotDf <- plotDf %>%
#   group_by(TSL) %>%
#   arrange(desc(value)) %>%
#   mutate(lab_ypos = cumsum(value) - 0.5 * value) 
# 
# # set order
# positions <- c("1", "5", "2", "3", "4" )
# 
# # Create stacked bar graphs with labels
# png(filename = "~/clinvar/TSLType_percent.png", width = 16, height = 8.135, units = "in", res = 300)
# ggplot(data = plotDf, aes(x = TSL, y = value)) +
#   geom_col(aes(fill = variable), width = 0.9)+
#   geom_text(aes(y = lab_ypos, label = value, group =variable), color = "white") +
#   theme_bw() +
#   theme(axis.text.x=element_text(angle=90)) +
#   scale_fill_manual(values = c(c("#999999", "#E69F00"))) +
#   scale_x_discrete(limits = positions) +
#   xlab("The TSL types") + 
#   ylab("% of variant count in each class") +
#   theme(
#     axis.title.x = element_text(size = 12),
#     axis.text.x = element_text(size = 10),
#     axis.title.y = element_text(size = 12))
# dev.off
# 

#######################Final dataset preparation

# up to now : "SYMBOL", "Allele", "CLNVC", "Consequence", "IMPACT", "relativeLocationRatio",

#"ExIntron", "Protein_position", "Codons", "TSL", 
############################### Ancestral genotype/allele ############################################
# AltaiNeandertal
# Ancestral_allele
# Denisova
# VindijiaNeandertal

table(subClin$AltaiNeandertal)
# -      ./.   A/A     A/C   A/G   A/T   C/C     C/G   C/T   G/G   T/T 
# 95587  10087 12438     4    27     1  29069     8    20   28585  11667 

table(subClin$Ancestral_allele)
# -        a     A      c     C     g     G       N     t     T 
# 97528   296   13514   878 31116   800  30419    15   260 12667 

table(subClin$Denisova)
# -       ./.   A/A      A/C   A/G   A/T   C/C     C/G   C/T   G/G     G/T   T/T 
# 95587 10104   12429     3    28     2   29065     9    32   28575     5   11654 

table(subClin$VindijiaNeandertal)
# -      ./.   A/A      A/C   A/G   A/T   C/C    C/G    C/T   G/G     G/T   T/T 
# 95587 10087  12434     3    25     3   29054    10    30   28592     1   11667

# cross-tabulation of Denisova and class:
#dropping missing values

tmp = subClin[c(30,162)]
tmp = tmp[tmp$Denisova != "-",]

# dropping low count group labels

tmp = tmp[-which(tmp$Denisova %in% c("A/C", "A/G", "A/T", "C/G","C/T", "G/T")),]

# ./.    A/A    C/C   G/G    T/T 
# 10104  12429  29065 28575  11654 


CrossTable(tmp$class, tmp$Denisova, prop.t=TRUE, prop.r=TRUE, prop.c=TRUE)


tmp = data.frame(unclass(table(tmp$Denisova, tmp$class)))
tmp$Denisova = rownames(tmp)
rownames(tmp) = NULL


plotDf = reshape2::melt(tmp, id.var = "Denisova")

# Arrange/sort and compute cumulative summs
plotDf <- plotDf %>%
  group_by(Denisova) %>%
  arrange(desc(value)) %>%
  mutate(lab_ypos = cumsum(value) - 0.5 * value) 

# Create stacked bar graphs with labels
png(filename = "~/clinvar/AncestralGenotypes.png", width = 16, height = 8.135, units = "in", res = 300)
ggplot(data = plotDf, aes(x = reorder(Denisova, -value), y = value)) +
  geom_col(aes(fill = variable), width = 0.9)+
  scale_fill_manual(values = c("#999999", "#E69F00")) +
  theme_bw() + 
  theme(axis.text.x=element_text(angle=90)) +
  xlab("TSL types") + 
  ylab("variant count in each class") +
  theme(
    axis.title.x = element_text(size = 12),
    axis.text.x = element_text(size = 10),
    axis.title.y = element_text(size = 12))
dev.off()

# plotting with percent values
Denisova = tmp$Denisova
conflict = round(100*tmp$conflict/(tmp$conflict+tmp$noConflict),0)
noConflict = round(100*tmp$noConflict/(tmp$conflict+tmp$noConflict),0)
tmp2 = data.frame(Denisova = Denisova, conflict = conflict, noConflict = noConflict)

plotDf = reshape2::melt(tmp2, id.var = "Denisova")
plotDf = plotDf[plotDf$Denisova != "-",]

# like what mentioned above
plotDf <- plotDf %>%
  group_by(Denisova) %>%
  arrange(desc(value)) %>%
  mutate(lab_ypos = cumsum(value) - 0.5 * value) 

# set order
positions <- c("C/C", "G/G", "A/A", "T/T", "./." )

# Create stacked bar graphs with labels
png(filename = "~/clinvar/DenisovaGenotype_percent.png", width = 16, height = 8.135, units = "in", res = 300)
ggplot(data = plotDf, aes(x = Denisova, y = value)) +
  geom_col(aes(fill = variable), width = 0.9)+
  geom_text(aes(y = lab_ypos, label = value, group =variable), color = "white") +
  theme_bw() +
  theme(axis.text.x=element_text(angle=90)) +
  scale_fill_manual(values = c(c("#999999", "#E69F00"))) +
  scale_x_discrete(limits = positions) +
  xlab("Denisova  Genotypes") + 
  ylab("% of variant count in each class") +
  theme(
    axis.title.x = element_text(size = 12),
    axis.text.x = element_text(size = 10),
    axis.title.y = element_text(size = 12))
dev.off

# grouping low freq genotype as other
table(subClin$Denisova)

Denisova = subClin$Denisova
lowfreq = c("A/C","A/G","A/T","C/G","C/T","G/T")
highfreq = c("./.", "A/A", "C/C", "G/G", "T/T")

Denisova = ifelse(Denisova %in% lowfreq, "other", 
                  ifelse(Denisova %in% highfreq, Denisova, NA))

finClin$Denisova = Denisova

# up to now : "SYMBOL", "Allele", "CLNVC", "Consequence", "IMPACT", "relativeLocationRatio",

#"ExIntron", "Protein_position", "Codons", "TSL", "Denisova",




############################### pathogenicity  score ############################################
### PUtting EVE result together
setwd("~/clinvar/variant_files/")

flist = list.files()

eveGene = gsub(".csv", "", flist)

uniprotID = read.table("~/clinvar/uniprotConversionTable.tab", fill = T, header = T)

ids = unique(subClin$GENEINFO)
ids = data.frame(stringr::str_split_fixed(ids, ":", 2))
names(ids) = c("symbol", "ENTREZ_ID")

# correcting ENTREZ numbers
entz = strsplit(ids$ENTREZ_ID,split='|',fixed=TRUE)
length(entz)
entz_num = c()
for (i in 1:length(entz)){
  entz_num[i] = entz[[i]][1]
}

ids$ENTREZ_ID = as.numeric(entz_num)


# joing two IDS

geneIDs = dplyr::left_join(ids, uniprotID)
geneIDs = geneIDs[!duplicated(geneIDs$symbol),]

# to see how many of genes with NAs in uniprot in genID dataset,
# Are really protein coding:

library(biomaRt)
mart <- useDataset("hsapiens_gene_ensembl", useMart("ensembl"))
genes <- getBM(attributes= c("ensembl_gene_id","entrezgene_id","hgnc_symbol", "gene_biotype"), 
               mart= mart)

protGene = genes[genes$gene_biotype == "protein_coding",]

table(geneIDs$ENTREZ_ID[is.na(geneIDs$uniprot)] %in% protGene$entrezgene_id)

# protein codng gene not in EVE:
notInEVE = geneIDs$ENTREZ_ID[is.na(geneIDs$uniprot)][geneIDs$ENTREZ_ID[is.na(geneIDs$uniprot)] %in% protGene$entrezgene_id]

notInEVE = genes[genes$entrezgene_id %in% notInEVE,]


# reading files and adding gene name to each files
vec =c("A" = "Ala",
       "G" = "Gly",
       "P" = "Pro",
       "T" = "Thr",
       "R" = "Arg",
       "N" = "Asn",
       "D" = "Asp", 
       "C" = "Cys",
       "E" = "Glu",
       "Q" = "Gln",
       "H" = "His",
       "I" = "Ile",
       "L" = "Leu",
       "K" = "Lys",
       "M" = "Met",
       "F" = "Phe",
       "S" = "Ser",
       "W" = "Trp",
       "Y" = "Tyr",
       "V" = "Val")

#
unik = subClin$GENEINFO[!duplicated(subClin$GENEINFO)]
eveList = list()

for (gene in 1:length(unik)){
  #selecting per gene
  subgen = subClin[subClin$GENEINFO == unik[gene],]
  # converting geneinfo to ids and merging with subgen
  ids = unique(subgen$GENEINFO)
  ids = data.frame(stringr::str_split_fixed(ids, ":", 2))
  names(ids) = c("symbol", "ENTREZ_ID")
  # correcting ENTREZ numbers
  entz = strsplit(ids$ENTREZ_ID,split='|',fixed=TRUE)
  entz_num = c()
  for (i in 1:length(entz)){
    entz_num[i] = entz[[i]][1]
  }
  ids$ENTREZ_ID = as.numeric(entz_num)
  geneIDs = dplyr::left_join(ids, uniprotID, by = "ENTREZ_ID")
  geneIDs = geneIDs[!duplicated(geneIDs$symbol),]
  subgen$uniprotID = geneIDs$uniprot
  subgen$ENTREZ_ID = geneIDs$ENTREZ_ID
  subgen$symbol = geneIDs$symbol
  #making mutation column in this dataset
  mut = data.frame(stringr::str_split_fixed(subgen$HGVSp, ":", 2))
  mut$symbol = subgen$symbol
  mut$mutation = paste0(mut$symbol,"_", mut$X2 )
  subgen$mutation  = mut$mutation
  if(!is.na(geneIDs$uniprot)){
    eveFrame = read.csv(paste0("~/clinvar/variant_files/", geneIDs$uniprot, ".csv"))
    eveFrame$uniprotID = geneIDs$uniprot
    eveFrame$symbol = geneIDs$symbol
    # substitution of one letter aa to three letter code
    eveFrame$wt_aaa =stringr::str_replace_all(eveFrame$wt_aa, vec)
    eveFrame$mt_aaa =stringr::str_replace_all(eveFrame$mt_aa, vec)
    eveFrame$mutation = paste0(eveFrame$symbol, "_p.", eveFrame$wt_aaa, eveFrame$position,
                               eveFrame$mt_aaa)
  }
  subgen = dplyr::left_join(subgen, eveFrame, by = "mutation")
  subgen = subgen[, c(1:456, 467,468,479)]
  subgen$unik = paste0(subgen$X.CHROM, subgen$POS, subgen$ID, subgen$REF, 
                       subgen$ALT, subgen$GENEINFO)
  subgen = subgen[, c(460, 456,457,458,459)]
  subgen = subgen[!is.na(subgen$EVE_scores_ASM),]
  # adding result to a list
  eveList[[gene]] = subgen
print(paste0(gene, " out of 4074, gene = ", eveFrame$symbol[1], " [" , round((100*gene/4047),2), "% completed]"))
}

eveFinalDataframe = do.call(rbind, eveList)


# adding eve data to subclin data

subClin$unik = paste0(subClin$X.CHROM, subClin$POS, subClin$ID, subClin$REF, 
                      subClin$ALT, subClin$GENEINFO)

subClin = dplyr::left_join(subClin, eveFinalDataframe, by = "unik")
#
#saveRDS(subClin, "~/clinvar/subClin_eve.RDS")

# adding eve to the dataset that I used for modeling inn python : finClin_2.csv
finalClin_2 = data.table::fread("~/clinvar/finalClin_2.csv")

finalClin_2 = dplyr::left_join(finalClin_2, subClin[, c(3,455)], by  = "ID")
finalClin_2 = finalClin_2[, -1]

# save 
data.table::fwrite(finalClin_2, "~/clinvar/finalClin_3.csv")


## lets do this for the jackson laboratory

jlab  = read.csv("~/perfect/JacksonLab_dataset.csv")
names(jlab)[1] <- names(geneIDs)[2]

jlab = dplyr::left_join(jlab, geneIDs, by = "ENTREZ_ID")

uniJlab = na.omit(jlab$uniprot[!duplicated(jlab$uniprot)])

uniJlabFile = paste0(uniJlab, ".csv")

list = c()

for(i in 1:length(uniJlabFile)){
  df = read.csv(paste0("~/clinvar/variant_files/", uniJlabFile[i]))
  df$uniprot_id = substr(uniJlabFile[i],1, nchar(uniJlabFile[i])-4)
  list[[i]] = df
}

eveDf = do.call(rbind, list)

eveDf$Variant = paste0(eveDf$uniprot_id, "_", eveDf$wt_aa, eveDf$position, eveDf$mt_aa)
jlab$Variant = paste0(jlab$uniprot, "_", jlab$Variant)

finDat = dplyr::inner_join(jlab, eveDf, by = "Variant")
finalDataFrame = finDat[, c(1:8, 21,22,33,54)]

finalDataFrame$EVE_classes_75_pct_retained_ASM[finalDataFrame$EVE_classes_75_pct_retained_ASM == ""] <- NA

View(unclass(table(finalDataFrame$Protein.Effect, finalDataFrame$EVE_classes_75_pct_retained_ASM)))

# read final data:
finalDataFrame = read.csv("~/perfect/Jlab_EVE_dataset.csv")


# visualization
library(ggplot2)
theme_set(
  theme_classic() +
    theme(legend.position = "top")
)
#
df = data.frame(table(finalDataFrame$Protein.Effect, finalDataFrame$EVE_classes_75_pct_retained_ASM))

ggplot(df, aes(x = Var1, y = Freq)) +
  geom_col(aes(color = Var2, fill = Var2), position = position_dodge(0.8), width = 0.7)
#


#based on the gens
df = data.frame(table(finalDataFrame$Gene_symbol, paste0(finalDataFrame$Protein.Effect, "_", finalDataFrame$EVE_classes_75_pct_retained_ASM)))
# names(df) = c("symbol", "jLabAnno", "EVE_pred", "Freq")
# df$jLabAnno[df$jLabAnno == "gain of function"] <- "GoF"

d = reshape(df, idvar = "Var1", timevar = "Var2", direction = "wide")
names(d) <- substr(names(d), 6, nchar(names(d)))
d = d[, - grep("_NA", names(d))]
rownames(d) = d[,1]
d = d[, -1]



library(pheatmap)
pheatmap(scale(d),main = "Jackson & EVE annotation")


#
#
# gene = df$symbol[!duplicated(df$symbol)]
#
# graphList = list()
# for(g in gene){
# p <- ggplot(df[df$symbol == g, ], aes(x = jLabAnno, y = Freq, group = EVE_pred, fill = EVE_pred)) +
#   geom_bar(stat="identity", position=position_dodge(width=0.9)) +
#   facet_wrap(. ~ symbol, ncol= 5)
# graphList[[g]] = p
# }
#
# library(ggpubr)
# theme_set(theme_pubr())
#
# ggarrange(
#   graphList[[1]],graphList[[2]] ,graphList[[3]],
#   common.legend = TRUE, legend = "bottom"
# )


# data on pathogenicity prediction scores are stored in metaData.xlsx:
# That data are grouped into three classes:
# 1- rank scores:
# Because different score has a different scaling system, dbNSFP developers,  created 
# a rank score for each score so that it is comparable between scores. 
#The rank score has a scale 0 to 1 and represents the percentage of scores 
#that are less damaging in dbNSFP, e.g., a rank score of 0.9 means the top 10% most damaging.

#
IdxRankScore = c(143L, 146L, 152L, 153L, 155L, 157L, 160L, 169L, 172L, 246L, 248L, 251L, 255L, 257L, 259L, 262L, 265L, 270L, 273L, 280L, 276L, 284L, 286L, 298L, 315L, 318L, 224L, 455L)

tmp = subClin[, c(30,IdxRankScore)]

# converting "-" to NA
tmp[tmp == "-"] <- NA
# removing rows with all NAs and converting dataset to numeric
tmp = tmp[rowSums(is.na(tmp)) != ncol(tmp)-1, ]

# see if columns contain more than one value
for (i in 1:ncol(tmp)){
  #print(table((stringr::str_count(tmp[,i], ","))))
  print(colnames(tmp)[i])
  print(table(grepl(",", tmp[,i])))
}

# convert to numbers

tmp[, -1] = sapply(tmp[, -1], as.numeric)

# identify columns with high number of NAs:
for(i in 1:ncol(tmp)){
  coN = colnames(tmp)[i]
  MissNum =colSums(is.na(tmp[i]))
  print(paste0("Number of missing values in ", coN, " is: ", MissNum, " ", "[", round(100*(MissNum/nrow(tmp)),2), "%]"))
}

missNumber = colSums(is.na(tmp))

# The "LINSIGHT_rankscore" and "MutPred_rankscore" have NAs more than 50% even after 
# removing rows with all NAs. These two will be dropped from final selection

tmp = tmp[, -which(colnames(tmp) %in% c("LINSIGHT_rankscore","MutPred_rankscore"))]

# correlation between ranks scores
# calculate mean for missing values
for(i in 1:ncol(tmp[,-1])){
  tmp[,-1][is.na(tmp[,-1][,i]), i] <- mean(tmp[,-1][,i], na.rm = TRUE)
}

library(corrplot)
library(RColorBrewer)

#trim column names in tmp
names(tmp) = gsub("_converted_rankscore", "", names(tmp))
names(tmp) = gsub("_rankscore", "", names(tmp))
names(tmp) = gsub("_scores_ASM", "", names(tmp))


M <-cor(tmp[,-1])

#png(filename = "~/clinvar/rankscore_corr.png", width = 15, height = 15, units = "in", res = 300)
corrplot(M, type="upper", order="hclust",
         col=brewer.pal(n=8, name="RdYlBu"))
dev.off()

# compute the pvalue for correlations

# mat : is a matrix of data
# ... : further arguments to pass to the native R cor.test function
cor.mtest <- function(mat, ...) {
  mat <- as.matrix(mat)
  n <- ncol(mat)
  p.mat<- matrix(NA, n, n)
  diag(p.mat) <- 0
  for (i in 1:(n - 1)) {
    for (j in (i + 1):n) {
      tmp <- cor.test(mat[, i], mat[, j], ...)
      p.mat[i, j] <- p.mat[j, i] <- tmp$p.value
    }
  }
  colnames(p.mat) <- rownames(p.mat) <- colnames(mat)
  p.mat
}
# matrix of the p-value of the correlation
p.mat <- cor.mtest(tmp[-1])

# visualization
col <- colorRampPalette(c("#BB4444", "#EE9988", "#FFFFFF", "#77AADD", "#4477AA"))
png(filename = "~/clinvar/rankscore_corr.png", width = 15, height = 15, units = "in", res = 300)

corrplot(M, method="color", col=col(200),  
         type="upper", order="hclust", 
         addCoef.col = "black", # Add coefficient of correlation
         tl.col="black", tl.srt=45, #Text label color and rotation
         tl.cex = 0.9, #for the size of text label (variable names)
         cl.cex = 0.9, # cex of number-label in colorlabel
         # Combine with significance
         p.mat = p.mat, sig.level = 0.01, insig = "blank", 
         # hide correlation coefficient on the principal diagonal
         diag=FALSE ,
         number.cex = 0.8
)
dev.off()

#Keeping only one variable out of pairs which show cor > 0.7
highlyCorrelated <- caret::findCorrelation(M, cutoff=(0.7),verbose = FALSE)
important_var_pathRank=colnames(tmp[,-highlyCorrelated])

important_var_pathRank


tmp = tmp[, c("class", important_var_pathRank)]

plotDf <- reshape2::melt(tmp, id.vars = 'class')

png(filename = "~/clinvar/selectedRankscore.png", width = 9, height = 7, units = "in", res = 300)
ggplot(plotDf, aes(x=variable, y=value, fill=class)) + 
  geom_boxplot(outlier.colour=NA) +
  scale_fill_manual(values = c("#999999", "#E69F00")) +
  facet_wrap(~variable, scale="free")
dev.off()


# melting
melt = reshape2::melt(tmp, id.var = "class")

ggplot(melt, aes(x = value, y = class)) +
  geom_density_ridges(aes(fill = class),) +
  scale_fill_manual(values = c("#999999", "#E69F00")) +
  ylab("") +
  scale_x_continuous(n.breaks = 5) +
  facet_wrap(~variable, scale="free")



#
chisq.test(table(subClin$class, subClin$EVE_classes_75_pct_retained_ASM))

# saving as table
tab = as.table(table(subClin$class, subClin$EVE_classes_75_pct_retained_ASM))

# nameing the axis
dimnames(tab) <- list(
  "Variant group" = c("conflict", "noConflict"),
  "EVE pathogenecity prediction" = c("Benign", "Pathogenic", "Uncertain")
)

spineplot(tab)

x <- prop.table(margin.table(tab, 1))
y <- prop.table(tab, 1)[2, ]
spineplot(tab, col = c("firebrick", "brown",  "goldenrod1"), axes = FALSE)
axis(1, at = c(0, x[1]) + x/2, labels = rownames(tab), tick = FALSE)
axis(2)
axis(4, at = c(0, y[1]) + y/2, labels = colnames(tab), tick = FALSE)

png(filename = "~/clinvar/EVE_classes_mosaic_plot.png", width = 9, height = 7, units = "in", res = 300)
mosaicplot(t(table(subClin$class, subClin$EVE_classes_75_pct_retained_ASM)), shade = TRUE,
           main = "EVE pathogeneicty classes and variant groups",
           cex.axis = 0.7)
dev.off()
# 2- original scores
IdxScore =c(72L, 73L, 108L, 109L, 110L, 111L, 113L, 114L, 115L, 116L, 117L, 118L, 119L, 144L, 147L, 148L, 149L, 150L, 151L, 156L, 158L, 161L, 167L, 168L, 170L, 171L, 245L, 249L, 250L, 253L, 256L, 258L, 260L, 263L, 266L, 271L, 274L, 279L, 282L, 285L, 287L, 299L, 316L, 319L, 450L, 452L, 226L, 455L)
 

tmp = subClin[, c(30,IdxScore)]
# dropping variables with high rate of NAs
tmp = tmp[, -which(colnames(tmp) %in% c("LINSIGHT","MutPred_score"))]


# converting SIFT/PolyPhen to numeric score and prediction
df = cbind(t(data.frame(strsplit(tmp$SIFT, "\\("))), t(data.frame(strsplit(tmp$PolyPhen, "\\("))))
colnames(df) <- c("SIFT_Pred", "SIFT_score", "PolyPhen_Pred", "PolyPhen_score")
df = data.frame(df)
df$SIFT_score = as.numeric(sub("\\)", "", df$SIFT_score))
df$PolyPhen_score = as.numeric(sub("\\)", "", df$PolyPhen_score))
rownames(df) <- NULL

tmp$SIFT = df$SIFT_score
tmp$PolyPhen = df$PolyPhen_score

# converting "-" to NA
tmp[tmp == "-"] <- NA
# removing rows with all NAs and converting dataset to numeric
tmp = tmp[rowSums(is.na(tmp)) != ncol(tmp)-1, ]

# see if columns contain more than one value
for (i in 1:ncol(tmp)){
  #print(table((stringr::str_count(tmp[,i], ","))))
  print(colnames(tmp)[i])
  print(table(grepl(",", tmp[,i])))
} 

multiScore = c("FATHMM_score", "VEST4_score", "PROVEAN_score", "MutationTaster_score","MutationAssessor_score",
               "MVP_score", "MPC_score", "LIST.S2_score","DEOGEN2_score")
# each feature should be inspected separately.
#################################################################
# FATHMM_score

  f = stringr::str_split_fixed(tmp[,"FATHMM_score"], ",", n = 100)
  f[f == ""] <- NA
  f = f[, colSums(is.na(f)) != nrow(f)]
  #convert rows with "." to have NAs
  f[f=="."] <- NA

  # loop
  u = c()
  for(i in 1:nrow(f)){
    if(length(na.omit(f[i,])) == 0){
      u[i] = NA
    }
    if(length(na.omit(f[i,])) == 1){
      u[i] = na.omit(f[i,])
    }
    if(length(unique(na.omit(f[i,]))) == 1){
      u[i] = unique(na.omit(f[i,]))
    }
    if(length(unique(na.omit(f[i,]))) > 1){ #If there's more than one scores associated with the same NS due to isoforms, the smallest score (most damaging) was used
      u[i] = min(as.numeric(na.omit(f[i,])))
    }
  }
tmp$FATHMM_score = u

# VEST4 : a random forest classifier for missense mutation that calculate probabality of pathigenicity
f = stringr::str_split_fixed(tmp[,"VEST4_score"], ",", n = 100)
f[f == ""] <- NA
f = f[, colSums(is.na(f)) != nrow(f)]
#convert rows with "." to have NAs
f[f=="."] <- NA
# check type of numbers
g = c()
for (i in 1:nrow(f)){
  mini = min(as.numeric(na.omit(f[i,])))
  maxi = max(as.numeric(na.omit(f[i,])))
  g[i] = ifelse(mini < 0 & maxi < 0, "neg", 
                ifelse(mini > 0 & maxi > 0, "pos",
                       ifelse(abs(mini) == Inf, NA,"mix")))
}
table(g)
# pos 
# 91688 > So, the largest number will be choose for analysis


# loop
u = c()
for(i in 1:nrow(f)){
  if(length(na.omit(f[i,])) == 0){
    u[i] = NA
  }
  if(length(na.omit(f[i,])) == 1){
    u[i] = na.omit(f[i,])
  }
  if(length(unique(na.omit(f[i,]))) == 1){
    u[i] = unique(na.omit(f[i,]))
  }
  if(length(unique(na.omit(f[i,]))) > 1){ 
    u[i] = max(as.numeric(na.omit(f[i,])))
  }
}
tmp$VEST4_score = u

# PROVEAN_score: 
#more negative score, more sever consequence [http://provean.jcvi.org/about.php#about_1]

f = stringr::str_split_fixed(tmp[,"PROVEAN_score"], ",", n = 100)
f[f == ""] <- NA
f = f[, colSums(is.na(f)) != nrow(f)]
#convert rows with "." to have NAs
f[f=="."] <- NA
# check type of numbers
g = c()
for (i in 1:nrow(f)){
  mini = min(as.numeric(na.omit(f[i,])))
  maxi = max(as.numeric(na.omit(f[i,])))
  g[i] = ifelse(mini < 0 & maxi < 0, "neg", 
                ifelse(mini > 0 & maxi > 0, "pos",
                       ifelse(abs(mini) == Inf, NA,"mix")))
}
table(g)
# mix   neg   pos 
# 1527 76683  6109  >>> so a mix of negative, positive score


# loop
u = c()
for(i in 1:nrow(f)){
  if(length(na.omit(f[i,])) == 0){
    u[i] = NA
  }
  if(length(na.omit(f[i,])) == 1){
    u[i] = na.omit(f[i,])
  }
  if(length(unique(na.omit(f[i,]))) == 1){
    u[i] = unique(na.omit(f[i,]))
  }
  if(length(unique(na.omit(f[i,]))) > 1){ 
    u[i] = min(as.numeric(na.omit(f[i,])))
  }
}
tmp$PROVEAN_score = u

# MutationTaster_score: 

f = stringr::str_split_fixed(tmp[,"MutationTaster_score"], ",", n = 100)
f[f == ""] <- NA
f = f[, colSums(is.na(f)) != nrow(f)]
#convert rows with "." to have NAs
f[f=="."] <- NA
# check type of numbers
g = c()
for (i in 1:nrow(f)){
  mini = min(as.numeric(na.omit(f[i,])))
  maxi = max(as.numeric(na.omit(f[i,])))
  g[i] = ifelse(mini < 0 & maxi < 0, "neg", 
                ifelse(mini > 0 & maxi > 0, "pos",
                       ifelse(abs(mini) == Inf, NA,"mix")))
}
table(g)
# mix   pos 
# 1     91491 >> mostly positive, higher more sever consequence


# loop
u = c()
for(i in 1:nrow(f)){
  if(length(na.omit(f[i,])) == 0){
    u[i] = NA
  }
  if(length(na.omit(f[i,])) == 1){
    u[i] = na.omit(f[i,])
  }
  if(length(unique(na.omit(f[i,]))) == 1){
    u[i] = unique(na.omit(f[i,]))
  }
  if(length(unique(na.omit(f[i,]))) > 1){ 
    u[i] = max(as.numeric(na.omit(f[i,])))
  }
}
tmp$MutationTaster_score = u

# MutationAssessor_score:
#MutationAssessor score (functional impact combined score) ranges from -5.545 to 5.975; 
#the larger the score the more likely it will be deleterious.

f = stringr::str_split_fixed(tmp[,"MutationAssessor_score"], ",", n = 100)
f[f == ""] <- NA
f = f[, colSums(is.na(f)) != nrow(f)]
#convert rows with "." to have NAs
f[f=="."] <- NA
# check type of numbers
g = c()
for (i in 1:nrow(f)){
  mini = min(as.numeric(na.omit(f[i,])))
  maxi = max(as.numeric(na.omit(f[i,])))
  g[i] = ifelse(mini < 0 & maxi < 0, "neg", 
                ifelse(mini > 0 & maxi > 0, "pos",
                       ifelse(abs(mini) == Inf, NA,"mix")))
}
table(g)
# mix   neg   pos 
# 1788  5428  69041 


# loop
u = c()
for(i in 1:nrow(f)){
  if(length(na.omit(f[i,])) == 0){
    u[i] = NA
  }
  if(length(na.omit(f[i,])) == 1){
    u[i] = na.omit(f[i,])
  }
  if(length(unique(na.omit(f[i,]))) == 1){
    u[i] = unique(na.omit(f[i,]))
  }
  if(length(unique(na.omit(f[i,]))) > 1){ 
    u[i] = max(as.numeric(na.omit(f[i,])))
  }
}
tmp$MutationAssessor_score = u

# MVP_score
# A pathogenicity prediction score for missense variants using deep learning approach. 
# The range of MVP score is from 0 to 1. The larger the score, the more likely the variant is pathogenic.

f = stringr::str_split_fixed(tmp[,"MVP_score"], ",", n = 100)
f[f == ""] <- NA
f = f[, colSums(is.na(f)) != nrow(f)]
#convert rows with "." to have NAs
f[f=="."] <- NA
# check type of numbers
g = c()
for (i in 1:nrow(f)){
  mini = min(as.numeric(na.omit(f[i,])))
  maxi = max(as.numeric(na.omit(f[i,])))
  g[i] = ifelse(mini < 0 & maxi < 0, "neg", 
                ifelse(mini > 0 & maxi > 0, "pos",
                       ifelse(abs(mini) == Inf, NA,"mix")))
}
table(g)
# pos 
# 82287 


# loop
u = c()
for(i in 1:nrow(f)){
  if(length(na.omit(f[i,])) == 0){
    u[i] = NA
  }
  if(length(na.omit(f[i,])) == 1){
    u[i] = na.omit(f[i,])
  }
  if(length(unique(na.omit(f[i,]))) == 1){
    u[i] = unique(na.omit(f[i,]))
  }
  if(length(unique(na.omit(f[i,]))) > 1){ 
    u[i] = max(as.numeric(na.omit(f[i,])))
  }
}
tmp$MVP_score = u

# MPC_score
#The range of MPC score is 0 to 5. 
#The larger the score, the more likely the variant is pathogenic. 
f = stringr::str_split_fixed(tmp[,"MPC_score"], ",", n = 100)
f[f == ""] <- NA
f = f[, colSums(is.na(f)) != nrow(f)]
#convert rows with "." to have NAs
f[f=="."] <- NA
# check type of numbers
g = c()
for (i in 1:nrow(f)){
  mini = min(as.numeric(na.omit(f[i,])))
  maxi = max(as.numeric(na.omit(f[i,])))
  g[i] = ifelse(mini < 0 & maxi < 0, "neg", 
                ifelse(mini > 0 & maxi > 0, "pos",
                       ifelse(abs(mini) == Inf, NA,"mix")))
}
table(g)



# loop
u = c()
for(i in 1:nrow(f)){
  if(length(na.omit(f[i,])) == 0){
    u[i] = NA
  }
  if(length(na.omit(f[i,])) == 1){
    u[i] = na.omit(f[i,])
  }
  if(length(unique(na.omit(f[i,]))) == 1){
    u[i] = unique(na.omit(f[i,]))
  }
  if(length(unique(na.omit(f[i,]))) > 1){ 
    u[i] = max(as.numeric(na.omit(f[i,])))
  }
}
tmp$MPC_score = u

# LIST.S2_score

# A deleteriousness preidction score for nonsynonymous SNVs. See https://doi.org/10.1093/nar/gkaa288.
# for details. The range of the score in dbNSFP is from 0 to 1.
# The higher the score, the more likely the variant is pathogenic.

f = stringr::str_split_fixed(tmp[,"LIST.S2_score"], ",", n = 100)
f[f == ""] <- NA
f = f[, colSums(is.na(f)) != nrow(f)]
#convert rows with "." to have NAs
f[f=="."] <- NA
# check type of numbers
g = c()
for (i in 1:nrow(f)){
  mini = min(as.numeric(na.omit(f[i,])))
  maxi = max(as.numeric(na.omit(f[i,])))
  g[i] = ifelse(mini < 0 & maxi < 0, "neg", 
                ifelse(mini > 0 & maxi > 0, "pos",
                       ifelse(abs(mini) == Inf, NA,"mix")))
}
table(g)



# loop
u = c()
for(i in 1:nrow(f)){
  if(length(na.omit(f[i,])) == 0){
    u[i] = NA
  }
  if(length(na.omit(f[i,])) == 1){
    u[i] = na.omit(f[i,])
  }
  if(length(unique(na.omit(f[i,]))) == 1){
    u[i] = unique(na.omit(f[i,]))
  }
  if(length(unique(na.omit(f[i,]))) > 1){ 
    u[i] = max(as.numeric(na.omit(f[i,])))
  }
}
tmp$LIST.S2_score = u


# DEOGEN2_score
# A deleteriousness prediction score "which incorporates heterogeneous information about
# the molecular effects of the variants, the domains involved, the relevance of the gene and the 
# interactions in which it participates". It ranges from 0 to 1. The larger the score, the more 
# likely the variant is deleterious.
f = stringr::str_split_fixed(tmp[,"DEOGEN2_score"], ",", n = 100)
f[f == ""] <- NA
f = f[, colSums(is.na(f)) != nrow(f)]
#convert rows with "." to have NAs
f[f=="."] <- NA
# check type of numbers
g = c()
for (i in 1:nrow(f)){
  mini = min(as.numeric(na.omit(f[i,])))
  maxi = max(as.numeric(na.omit(f[i,])))
  g[i] = ifelse(mini < 0 & maxi < 0, "neg", 
                ifelse(mini > 0 & maxi > 0, "pos",
                       ifelse(abs(mini) == Inf, NA,"mix")))
}
table(g)



# loop
u = c()
for(i in 1:nrow(f)){
  if(length(na.omit(f[i,])) == 0){
    u[i] = NA
  }
  if(length(na.omit(f[i,])) == 1){
    u[i] = na.omit(f[i,])
  }
  if(length(unique(na.omit(f[i,]))) == 1){
    u[i] = unique(na.omit(f[i,]))
  }
  if(length(unique(na.omit(f[i,]))) > 1){ 
    u[i] = max(as.numeric(na.omit(f[i,])))
  }
}
tmp$DEOGEN2_score = u
#################################################################
# correlation analysis and selecting features:
# identify columns with high number of NAs:

missNumber = data.frame(missnum = colSums(is.na(tmp)))
missNumber$percent = missNumber$missnum/dim(tmp)[1]
missNumber = missNumber[order(missNumber$percent),]

highPercentageNAs = rownames(missNumber)[missNumber$percent > 0.80]

tmp = tmp[, -which(colnames(tmp) %in% highPercentageNAs)]

# calculate mean for missing values
tmp[,-1] = sapply(tmp[,-1], as.numeric)

for(i in 1:ncol(tmp[,-1])){
  tmp[,-1][is.na(tmp[,-1][,i]), i] <- mean(tmp[,-1][,i], na.rm = TRUE)
}

names(tmp) = gsub("_score","", names(tmp))
# correlation plot
M <-cor(tmp[,-1])


# matrix of the p-value of the correlation
p.mat <- cor.mtest(tmp[-1])

# visualization
col <- colorRampPalette(c("#BB4444", "#EE9988", "#FFFFFF", "#77AADD", "#4477AA"))
png(filename = "~/clinvar/dbpnsfpScore_corr.png", width = 15, height = 15, units = "in", res = 300)

corrplot(M, method="color", col=col(200),  
         type="upper", order="hclust", 
         addCoef.col = "black", # Add coefficient of correlation
         tl.col="black", tl.srt=45, #Text label color and rotation
         tl.cex = 0.9, #for the size of text label (variable names)
         cl.cex = 0.9, # cex of number-label in colorlabel
         # Combine with significance
         p.mat = p.mat, sig.level = 0.01, insig = "blank", 
         # hide correlation coefficient on the principal diagonal
         diag=FALSE ,
         number.cex = 0.5
)
dev.off()



png(filename = "~/clinvar/dbpnsfpScore_corr.png", width = 15, height = 15, units = "in", res = 300)
corrplot(M, type="upper", order="hclust",
         col=brewer.pal(n=8, name="RdYlBu"))
dev.off()

#Keeping only one variable out of pairs which show cor > 0.7
highlyCorrelated <- caret::findCorrelation(M, cutoff=(0.7),verbose = FALSE)
important_var=colnames(tmp[,-highlyCorrelated])


important_var

# [1] "class"                   "SIFT"                    "PolyPhen"               
# [4] "SpliceAI_pred_DS_AG"     "SpliceAI_pred_DS_AL"     "SpliceAI_pred_DS_DG"    
# [7] "SpliceAI_pred_DS_DL"     "CADD_raw"                "CADD_raw_hg19"          
# [10] "ClinPred_score"          "DANN_score"              "DEOGEN2_score"          
# [13] "Eigen.raw_coding"        "LIST.S2_score"           "LRT_Omega"              
# [16] "LRT_score"               "M.CAP_score"             "MPC_score"              
# [19] "MetaSVM_score"           "MutationAssessor_score"  "MutationTaster_score"   
# [22] "PROVEAN_score"           "fathmm.MKL_coding_score" "fathmm.XF_coding_score" 
# [25] "LoFtool"                 "BLOSUM62"                "FATHMM_score"   



# Original dbNSFP paper provides a correlation matrix. Looking at this correlation matrix also might be helpful
# to select feature from here.
corMat = read.csv("~/clinvar/13073_2020_803_MOESM2_ESM.csv", row.names = 1)
corMat<-as.matrix(corMat)
# Getting rid of the lower triangle
# Upper triangle: correlation coefficients. Lower triangle: observation counts.
corMat[lower.tri(corMat)]<-NA

#
highlyCorrelated <- caret::findCorrelation(corMat, cutoff=(0.7),verbose = FALSE)
important_var=colnames(corMat[,-highlyCorrelated])


important_var

# [1] "SIFT"                       "SIFT4G"                     "Polyphen2_HDIV"            
# [4] "Polyphen2_HVAR"             "LRT"                        "MutationTaster"            
# [7] "MutationAssessor"           "DEOGEN2"                    "BayesDel_addAF"            
# [10] "BayesDel_noAF"              "ClinPred"                   "LIST1S2"                   
# [13] "REVEL"                      "MutPred"                    "MVP"                       
# [16] "MPC"                        "CADD"                       "CADD_hg19"                 
# [19] "Eigen1PC"                   "GenoCanyon"                 "fitCons_i6"                
# [22] "fitCons_gm"                 "GERP.."                     "phyloP100way_vertebrate"   
# [25] "phyloP30way_mammalian"      "phyloP17way_primate"        "phastCons100way_vertebrate"
# [28] "phastCons30way_mammalian"   "phastCons17way_primate"     "SiPhy"                     
# [31] "bStatistic"    



# 3- prediction class
# Numeric scores is more favorable by algorithms, so will leave this part for now.

###################
predIdx = c(142L, 145L, 154L, 159L, 247L, 252L, 254L, 261L, 264L, 272L, 278L, 281L, 283L, 314L, 317L, 225L, 457)

tmp = subClin[, c(30,predIdx)]
# dropping variables with high rate of NAs: all seems to have < 60% missing counts
tmp[tmp == "-"] <- NA
tmp[tmp == ".,."] <- NA



missNumber = data.frame(missingCount = colSums(is.na(tmp)))
missNumber$percent = missNumber$missingCount/dim(missNumber)[1]

# looking at values in different columns
for(i in 1:ncol(tmp)){
  print(colnames(tmp)[i])
  print(dim(table(tmp[,i])))
}

# [1] "class"
# [1] 2
# [1] "BayesDel_addAF_pred"
# [1] 2
# [1] "BayesDel_noAF_pred"
# [1] 2
# [1] "ClinPred_pred"
# [1] 2
# [1] "DEOGEN2_pred"
# [1] 2222
# [1] "LIST.S2_pred"
# [1] 2214
# [1] "LRT_pred"
# [1] 3
# [1] "M.CAP_pred"
# [1] 2
# [1] "MetaLR_pred"
# [1] 2
# [1] "MetaSVM_pred"
# [1] 2
# [1] "MutationAssessor_pred"
# [1] 2338
# [1] "MutationTaster_pred"
# [1] 308
# [1] "PROVEAN_pred"
# [1] 2548
# [1] "PrimateAI_pred"
# [1] 2
# [1] "fathmm.MKL_coding_pred"
# [1] 2
# [1] "fathmm.XF_coding_pred"
# [1] 2
# [1] "FATHMM_pred"
# [1] 1819


# dealing with feature with more than one value per row
multival  = c("DEOGEN2_pred", "LIST.S2_pred", "MutationAssessor_pred","MutationTaster_pred", "PROVEAN_pred", "FATHMM_pred")
u = c()

for(score in 1:length(multival)){
  f = stringr::str_split_fixed(tmp[,multival[score]], ",", n = 100)
  f[f == ""] <- NA
  f[f== "."] <- NA
  f = f[, colSums(is.na(f)) != nrow(f)]
  for(i in 1:nrow(f)){
    if(is.na(unique(as.vector((f[i,]))))){
      u[i] = NA
    }
    if(!is.na(unique(as.vector((f[i,]))))){
    u[i] = unique(as.vector(na.omit(f[i,])))    
    }
  }
tmp[, multival[score]] = u
}

# looking at values in different columns to check again
for(i in 1:ncol(tmp)){
  print(colnames(tmp)[i])
  print(table(tmp[,i]))
}


# see if columns contain more than one value
for (i in 1:ncol(tmp)){
  #print(table((stringr::str_count(tmp[,i], ","))))
  print(colnames(tmp)[i])
  print(table(grepl(",", tmp[,i])))
} 



# converting SIFT/PolyPhen to numeric score and prediction
df = cbind(t(data.frame(strsplit(subClin$SIFT, "\\("))), t(data.frame(strsplit(subClin$PolyPhen, "\\("))))
colnames(df) <- c("SIFT_Pred", "SIFT_score", "PolyPhen_Pred", "PolyPhen_score")
df = data.frame(df)
df$SIFT_score = as.numeric(sub("\\)", "", df$SIFT_score))
df$PolyPhen_score = as.numeric(sub("\\)", "", df$PolyPhen_score))
rownames(df) <- NULL


# removing original SIFT and polyPhen columns
patDat = cbind(tmp, df)
patDat = patDat[, -c(20,22)]

# replacing - with NA
patDat[patDat == "-"] <- NA
# removing rows with all NAs and converting dataset to numeric
patDat = patDat[rowSums(is.na(patDat)) != ncol(patDat)-1, ]

# how many NAs and reducing by NA numbers
patDat = patDat[, colSums(is.na(patDat)) <= 150000]
patDat = patDat[rowSums(is.na(patDat)) <= 50, ]

names(patDat) = sub("_pred", "", names(patDat))
names(patDat) = sub("_Pred", "", names(patDat))
names(patDat) = sub("_classes_75_pct_retained_ASM", "", names(patDat))

# workinh of SIFT and polyphen
patDat$PolyPhen[patDat$PolyPhen == "unknown"] <- NA
patDat$PolyPhen[patDat$PolyPhen == "possibly_damaging"] <- "possiblyD"
patDat$PolyPhen[patDat$PolyPhen == "probably_damaging"] <- "probablyD"

patDat$SIFT[patDat$SIFT == "deleterious_low_confidence"] <- "deleterious"
patDat$SIFT[patDat$SIFT == "tolerated_low_confidence"] <- "tolerated"


crosTab = data.frame(unclass(table(patDat$class,patDat$BayesDel_addAF)))
#names(crosTab) <- paste0("BayesDel_addAF_", names(crosTab))
crosTab$score = names(crosTab)
crosTab$predAlgo = "BayesDel_addAF"


for (i in 3:dim(patDat[2])){
  print(names(patDat[i]))
  print(i)
  d = data.frame(unclass(table(patDat$class,patDat[,i])))
  names(d) <- paste0(names(patDat)[i],"_", names(d))
  crosTab = cbind(crosTab, d)
}

crosTab = t(crosTab)
dt <- as.table(as.matrix(crosTab))

# melt patDat
meltpatDat = reshape2::melt(patDat, id.var = "class")




library(ggplot2)
theme_set(
  theme_classic() +
    theme(legend.position = "top")
)

pvalvector = c()

#
for (i in 2:dim(patDat)[2]){
  df = data.frame(table(patDat$class, patDat[,i]))
  c= chisq.test(table(patDat$class, patDat[,i]))
  pvalvector = c(pvalvector,c$p.value)
 }

# calculate effect size
# or
OR = list()

for(i in 2:dim(patDat)[2]){
  if(dim(table(patDat$class, patDat[,i]))[2] != 2){
    print(paste0("Other table", names(patDat)[i]))
    carmer = lsr::cramersV(table(patDat$class, patDat[,i]))# to interpret carmers V
    OR[[i]] = paste0("df = ",dim(table(patDat$class, patDat[,i]))[2]-1, ", Carmer's V: " , round(carmer,2))
  } else {
    print(paste0("2x2 table", names(patDat)[i]))
    or = abd::odds.ratio(table(patDat$class, patDat[,i]))
    OR[[i]] = paste0("OR :" , round(or$OR,2), "[", round(or$lower,2), "-", round(or$upper,2), "]") 
  }
}


adjPval = p.adjust(pvalvector, method = "BH")
adjPval = as.character(format(adjPval, scientific = T, digits = 2))

#
images = list()
for (i in 2:dim(patDat)[2]){
  df = data.frame(table(patDat$class, patDat[,i]))
  p = ggplot(df, aes(x = Var2, y =Freq)) +
    geom_col(aes(fill = Var1), position = position_dodge(0.8), width = 0.7) +
    scale_fill_manual(values = c("#999999", "#E69F00")) +
    labs(x=paste0(names(patDat)[i],"\n",OR[i], ", adjP: ", adjPval[i-1] ), y= "Frequency", fill = "class")
  images[[i]] <- p
}

library(ggpubr)
theme_set(theme_pubr())

png(filename = "~/clinvar/pathPred.png", width = 20, height = 12, units = "in", res = 300)
ggarrange(
  images[[2]], images[[3]],images[[4]], images[[5]],images[[6]],
  images[[7]], images[[8]],images[[9]], images[[10]],images[[11]],
  images[[12]], images[[13]],images[[14]], images[[15]],images[[16]],
  images[[17]], images[[18]],images[[19]], images[[20]],
  ncol = 5, nrow = 4,
  common.legend = TRUE, legend = "bottom"
)
dev.off()

################
# final dataset preparation:
df = cbind(t(data.frame(strsplit(tmp$SIFT, "\\("))), t(data.frame(strsplit(tmp$PolyPhen, "\\("))))
colnames(df) <- c("SIFT_Pred", "SIFT_score", "PolyPhen_Pred", "PolyPhen_score")
df = data.frame(df)
df$SIFT_score = as.numeric(sub("\\)", "", df$SIFT_score))
df$PolyPhen_score = as.numeric(sub("\\)", "", df$PolyPhen_score))
rownames(df) <- NULL

finClin$SIFT = df$SIFT_score
finClin$PolyPhen = df$PolyPhen_score


# up to now : "SYMBOL", "Allele", "CLNVC", "Consequence", "IMPACT", "relativeLocationRatio",

#"ExIntron", "Protein_position", "Codons", "TSL", "Denisova",SIFT, "PolyPhen
# plus c(important_var_pathRank, important_var_AS_AC_AN_nhomalt)

##############################Conservation score ##############################

# rankscore
consRankScoreIdx = c(229L, 231L, 235L, 239L, 242L, 290L, 302L, 436L, 439L, 441L, 443L, 445L, 447L, 449L)


tmp = subClin[ ,c(30,consRankScoreIdx)]

# remove rows with nearlly all NAs
table(rowSums(is.na(tmp)) > 14)

# FALSE 
# 187643 

#convert to numeric
tmp[,-1] = sapply(tmp[,-1], as.numeric)

# count missing values
missNumber = data.frame(mis = colSums(is.na(tmp)))
missNumber$percent = missNumber$mis/dim(tmp)[1]

# computing mean for missing values
for(i in 1:ncol(tmp[,-1])){
  tmp[,-1][is.na(tmp[,-1][,i]), i] <- mean(tmp[,-1][,i], na.rm = TRUE)
}

names(tmp) = sub("_rankscore", "", names(tmp))





M <-cor(tmp[,-1])

#png(filename = "~/clinvar/ConsRankscore_corr.png", width = 15, height = 15, units = "in", res = 300)
corrplot(M, type="upper", order="hclust",
         col=brewer.pal(n=8, name="RdYlBu"))
#dev.off()

# visualization:

plotDf <- reshape2::melt(tmp, id.vars = 'class')

#png(filename = "~/clinvar/selectedRankscore.png", width = 9, height = 7, units = "in", res = 300)
ggplot(plotDf, aes(x=variable, y=value, fill=class)) + 
  geom_bar(outlier.colour=NA) +
  scale_fill_manual(values = c("#999999", "#E69F00")) +
  facet_wrap(~variable, scale="free")
#dev.off()


# melting
melt = reshape2::melt(tmp, id.var = "class")

ggplot(melt, aes(x = value, y = class)) +
  geom_density_ridges(aes(fill = class),) +
  scale_fill_manual(values = c("#999999", "#E69F00")) +
  ylab("") +
  scale_x_continuous(n.breaks = 5) +
  facet_wrap(~variable, scale="free")



# visualization by pvalue

ggbarplot(plotDf, x = "class", y = "value",add = "mean_se",
               fill  = "class", palette = c("#999999", "#E69F00"),
               position = position_dodge(0.8),
               facet.by = "variable", short.panel.labs = FALSE) +
              stat_compare_means(aes(group = class), label = "p.signif")


library(dplyr)
group_by(tmp, class) %>%
  summarise(
    count = n(),
    mean = mean(GERP.._RS, na.rm = TRUE),
    sd = sd(GERP.._RS, na.rm = TRUE)
  )
res <- t.test(GERP.._RS ~ class, data = tmp, var.equal = TRUE)
res

#Keeping only one variable out of pairs which show cor > 0.7
highlyCorrelated <- caret::findCorrelation(M, cutoff=(0.7),verbose = FALSE)
important_var_ConsRank=colnames(tmp[,-highlyCorrelated])

important_var_ConsRank
# [1] "GERP.._RS_rankscore"                  "GM12878_fitCons_rankscore"           
# [3] "GenoCanyon_rankscore"                 "H1.hESC_fitCons_rankscore"           
# [5] "HUVEC_fitCons_rankscore"              "SiPhy_29way_logOdds_rankscore"       
# [7] "integrated_fitCons_rankscore"         "phastCons100way_vertebrate_rankscore"
# [9] "phyloP100way_vertebrate_rankscore"    "phyloP30way_mammalian_rankscore" 

# original scores

consScoreIdx = c(227L, 228L, 230L, 232L, 233L, 234L, 236L, 237L, 238L, 240L, 241L, 243L, 289L, 291L, 301L, 435L, 437L, 438L, 440L, 442L, 444L, 446L, 448L)

tmp = subClin[ ,c(30,consScoreIdx)]
tmp[tmp == "-"] <- NA
#convert to numeric

tmp[,-1] = sapply(tmp[,-1], as.numeric)

# count missing values
missNumber = data.frame(mis = colSums(is.na(tmp)))
missNumber$percent = missNumber$mis/dim(tmp)[1]

# dropping feature with hign rate of NAs
ex = c("GTEx_V8_gene", "GTEx_V8_tissue", "Geuvadis_eQTL_target_gene", "SiPhy_29way_pi")

tmp = tmp[, -which(colnames(tmp) %in% ex)]

conscoreRes = data.frame(score = c(), class = c(), AVG = c(), SD = c(), effectSize = c(), p_value = c())

for(i in 2:dim(tmp)[2]){
class = c("conflict", "noConflict")
AVG = c(round(median(tmp[,i][tmp$class == "conflict"], na.rm = T),2), round(median(tmp[,i][tmp$class == "noConflict"], na.rm = T),2))
SD = c(paste0(round(quantile(tmp[,i][tmp$class == "conflict"], na.rm = T)[2],2), "-",
       round(quantile(tmp[,i][tmp$class == "conflict"], na.rm = T)[4],2)),
       paste0(round(quantile(tmp[,i][tmp$class == "noConflict"], na.rm = T)[2],2), "-",
       round(quantile(tmp[,i][tmp$class == "noConflict"], na.rm = T)[4],2)))
res <- wilcox.test(tmp[,i] ~ class, data = tmp)
atest = data.frame(a = tmp[, i], class = tmp$class)
wilcoxEffsize= rstatix::wilcox_effsize(a ~ class, data = atest)
effectSize = paste0(round(wilcoxEffsize$effsize,2), " [", wilcoxEffsize$magnitude, "]")
score = names(tmp)[i]
p_value = res$p.value

d = data.frame(score, class, AVG, SD,effectSize, p_value)
conscoreRes = rbind(conscoreRes,d)
}

pvec =conscoreRes$p_value[!duplicated(conscoreRes$p_value)]
adjP = p.adjust(pvec, method = 'bonferroni')
conscoreRes$adj.p_value <- NA
vec = seq(1,38,2)

for(i in 1:length(vec)){
  conscoreRes$adj.p_value[vec[i]] = adjP[i]
}

conscoreRes$AVG = paste0(conscoreRes$AVG, "+/-", conscoreRes$SD )
conscoreRes$p_value = as.character(conscoreRes$p_value)
conscoreRes$adj.p_value = as.character(conscoreRes$adj.p_value)


library(rstatix)
wilcox_effsize(GERP.._NR ~ class, data = tmp)


# save into a word document
sjPlot::tab_df(conscoreRes, digits = 4,
               file="~/clinvar/paper/conservation score.doc") 
#


# computing mean for missing values
for(i in 1:ncol(tmp[,-1])){
  tmp[,-1][is.na(tmp[,-1][,i]), i] <- mean(tmp[,-1][,i], na.rm = TRUE)
}

# correlation
M <-cor(tmp[,-1])

#png(filename = "~/clinvar/ConsScore_corr.png", width = 15, height = 15, units = "in", res = 300)
corrplot(M, type="upper", order="hclust",
         col=brewer.pal(n=8, name="RdYlBu"))
#dev.off()

#Keeping only one variable out of pairs which show cor > 0.7
highlyCorrelated <- caret::findCorrelation(M, cutoff=(0.7),verbose = FALSE)
important_var_ConsOrgScore=colnames(tmp[,-highlyCorrelated])

important_var_ConsOrgScore = important_var_ConsOrgScore[-1]
important_var_ConsOrgScore
# correlation between rank and org score:
tmp = subClin[, c(important_var_ConsOrgScore,important_var_ConsRank)]
tmp  = sapply(tmp, as.numeric)
# computing mean for missing values
for(i in 1:ncol(tmp)){
  tmp[is.na(tmp[,i]), i] <- mean(tmp[,i], na.rm = TRUE)
}

# correlation
M <-cor(tmp)

png(filename = "~/clinvar/ConsScoreRankscore_corr.png", width = 15, height = 15, units = "in", res = 300)
corrplot(M, type="upper", order="hclust",
         col=brewer.pal(n=8, name="RdYlBu"))
dev.off()

#Keeping only one variable out of pairs which show cor > 0.7
highlyCorrelated <- caret::findCorrelation(M, cutoff=(0.7),verbose = FALSE)
important_var_ConsOrgRankScore=colnames(tmp[,-highlyCorrelated])

important_var_ConsOrgRankScore


# up to now : "SYMBOL", "Allele", "CLNVC", "Consequence", "IMPACT", "relativeLocationRatio",

#"ExIntron", "Protein_position", "Codons", "TSL", "Denisova",SIFT, "PolyPhen
# plus c(important_var_pathRank, important_var_AS_AC_AN_nhomalt, important_var_ConsOrgRankScore)

################################ submission/submitter info #######################################
# there different files (data downloaded from ClinVar FTP on 2021-08-12)

subInfo1 <- data.table::fread("~/clinvar/submission_summary.txt")
subInfo2 <- data.table::fread("~/clinvar/summary_of_conflicting_interpretations.txt")
subInfo3 <- data.table::fread("~/clinvar/variant_summary.txt")



# distribution of submitter for variants
submitterCount = data.frame(table(subInfo1$`#VariationID`))
colnames(submitterCount) = c("ID", "Subcount")
#overlap with Subclin
table(subClin$ID %in% submitterCount$ID)

# joining two datasets
subClin$ID  = as.factor(subClin$ID)
subCount = dplyr::left_join(submitterCount[which(submitterCount$ID %in% subClin$ID),], subClin[, c(3,30)])


# mean of submission count between two groups
t.test(Subcount ~ class, data = subCount, alternative = "two.sided", var.equal = FALSE)
rstatix::cohens_d(Subcount ~ class, data = subCount, var.equal = TRUE)
# Welch Two Sample t-test
# 
# data:  Subcount by class
# t = 47.8, df = 70431, p-value < 2.2e-16
# alternative hypothesis: true difference in means is not equal to 0
# 95 percent confidence interval:
#   0.6216439 0.6748036
#  sample estimates:
#   mean in group conflict mean in group noConflict 
# 3.894831                 3.246608 

# given the non-normal distribution of data 
library(dplyr)
group_by(subCount, class) %>%
  summarise(
    count = n(),
    median = median(Subcount, na.rm = TRUE),
    LowerQ = summary(s)[2],
    UpperQ = summary(subCount)[4]
    
  )

# visualization
png(filename = "~/clinvar/submitter_count.png", width = 9, height = 7, units = "in", res = 300)
ggplot(subCount, aes(x=class, y=log10(Subcount), fill=class)) + 
  geom_boxplot() +
  scale_fill_manual(values = c("#999999", "#E69F00"))
dev.off()

# making a df for each variant ID , first submitter, date and count of submission and other variables

# filtering based on variant ids

subSum = subInfo1[subInfo1$`#VariationID` %in% subClin$ID,]

ids = unique(subSum$`#VariationID`)
varID = c()
submitterNo = c()
firstSumitter = c()
firstSumitterDate = c()
CollectionMethod = c()

for (i in 1:length(ids)){
  df = subSum[subSum$`#VariationID` == ids[i], ]
  varID[i] = ids[i]
  submitterNo [i] = dim(df)[1]
  df$DateLastEvaluated = sub(",", "", df$DateLastEvaluated)
  df$DateLastEvaluated = as.Date(df$DateLastEvaluated,format = "%B %d %Y")
  df = df[order(df$DateLastEvaluated),]
  firstSumitterDate[i] = as.character(df$DateLastEvaluated[1])
  firstSumitter[i] = df$Submitter[1]
  CollectionMethod[i] = df$CollectionMethod[1]
}

submitterDf = data.frame(ID = as.factor(varID), firstSumitter = firstSumitter,firstSumitterDate = firstSumitterDate,
                         submitterNo = submitterNo, CollectionMethod = CollectionMethod)

# 
tmp = data.frame(unclass(table(submitterDf$firstSumitter)))
tmp$submitter = rownames(tmp)
topSubmitter = tmp$submitter[tmp$unclass.table.submitterDf.firstSumitter.. >= 1000]

# A total of 507 centers contributed to submission of 187643  variants.
# top centers with highest number of submission

submitterDf = dplyr::left_join(submitterDf, subClin[, c(3,30)], by = "ID")

library(dplyr)

tmp = submitterDf %>%
  group_by(firstSumitter, class) %>%
  summarise(count = n())

tmp = tmp[tmp$firstSumitter %in% topSubmitter, ]
tmp = data.frame(tmp)

# making long names, shorter
tmp$firstSumitter[tmp$firstSumitter == "ARUP Laboratories, Molecular Genetics and Genomics,ARUP Laboratories"] <- "ARUP Laboratories"
tmp$firstSumitter[tmp$firstSumitter == "Biesecker Lab/Clinical Genomics Section,National Institutes of Health"] <- "Biesecker Lab_NIH"
tmp$firstSumitter[tmp$firstSumitter == "CeGaT Praxis fuer Humangenetik Tuebingen"] <- "CeGaT Tuebingen"
tmp$firstSumitter[tmp$firstSumitter == "Center for Pediatric Genomic Medicine,Children's Mercy Hospital and Clinics"] <- "Children's Mercy Hospital"
tmp$firstSumitter[tmp$firstSumitter == "EGL Genetic Diagnostics, Eurofins Clinical Diagnostics"] <- "Eurofins Clinical Diagnostics"
tmp$firstSumitter[tmp$firstSumitter == "Genetic Services Laboratory, University of Chicago"] <- "University of Chicago"
tmp$firstSumitter[tmp$firstSumitter == "Illumina Clinical Services Laboratory,Illumina"] <- "Illumina Clinical Laboratory"
tmp$firstSumitter[tmp$firstSumitter == "Laboratory for Molecular Medicine, Partners HealthCare Personalized Medicine"] <- "Partners HealthCare Person. Med."
tmp$firstSumitter[tmp$firstSumitter == "Quest Diagnostics Nichols Institute San Juan Capistrano"] <- "Quest Diagnostics"
tmp$firstSumitter[tmp$firstSumitter == "Women's Health and Genetics/Laboratory Corporation of America, LabCorp"] <- "LabCorp"

# convert counts into percent values
num = seq(1,37, 2)
relFreq = as.numeric()

for(i in num){
  relFreq[i] = round(100 * (tmp$count[i]/(tmp$count[i] + tmp$count[i+1])),2)
  relFreq[i+1] = round(100 * (tmp$count[i+1]/(tmp$count[i] + tmp$count[i+1])),2)
}

tmp$relfreq = relFreq

# adding total number to the center name
tot = c()
num = seq(1,37, 2)
for(i in num){
  tot[i] = tmp$count[i] + tmp$count[i+1]
  tot[i+1] = tot[i]
}
tmp$tot = tot

# adding total number of variant to the center name
centerTotCount = paste0(tmp$firstSumitter," [", tmp$tot, "]")
tmp$firstSumitter = centerTotCount

# make a group value as negative
tmp$relfreq[tmp$class == "conflict"] = -tmp$relfreq[tmp$class == "conflict"]


# X Axis Breaks and Labels 
brks <- seq(-100, 100, 10)
brks
lbls <- c('100','90','80','70','60','50','40','30','20','10','0','10', '20','30','40','50','60','70','80','90','100')

png(filename = "~/clinvar/submitter_class_2.png", width = 12, height = 9, units = "in", res = 300)

ggplot(tmp, aes(x = reorder(firstSumitter, -tot), y = relfreq, fill = class)) +   # Fill column
  geom_bar(stat = "identity", width = .6) +   # draw the bars
  scale_y_continuous(breaks = brks, labels = lbls) + # Labels
  coord_flip() +  # Flip axes
  labs(title="Relative frequency for conflict/no-conflict submission") +
  xlab('Genetic center [total submissions as the first submitter]') +
  ylab('%Submission') +
  theme_tufte() +  # Tufte theme from ggthemes
  theme(plot.title = element_text(hjust = .5),
        axis.ticks = element_blank(),
        text=element_text(size=15)) +   # Centre plot title
  scale_fill_manual(values=c('#999999','#E69F00')) +  # Color palette
  guides(fill=guide_legend("Class")) +
  geom_hline(yintercept = 50, colour="blue", linetype="dashed") +
  geom_hline(yintercept = 60, colour="blue", linetype="dashed") +
  geom_hline(yintercept = 70, colour="blue", linetype="dashed") +
  geom_hline(yintercept = 80, colour="blue", linetype="dashed") +
  geom_hline(yintercept = 90, colour="blue", linetype="dashed") +
  geom_hline(yintercept = -10, colour="blue", linetype="dashed") +
  geom_hline(yintercept = -20, colour="blue", linetype="dashed") +
  geom_hline(yintercept = -30, colour="blue", linetype="dashed") +
  geom_hline(yintercept = -40, colour="blue", linetype="dashed")

dev.off()

#make a df for visualization

tmp = submitterDf[submitterDf$firstSumitter %in% topSubmitter, ]

# making long names, shorter
tmp$firstSumitter[tmp$firstSumitter == "ARUP Laboratories, Molecular Genetics and Genomics,ARUP Laboratories"] <- "ARUP Laboratories"
tmp$firstSumitter[tmp$firstSumitter == "Biesecker Lab/Clinical Genomics Section,National Institutes of Health"] <- "Biesecker Lab_NIH"
tmp$firstSumitter[tmp$firstSumitter == "CeGaT Praxis fuer Humangenetik Tuebingen"] <- "CeGaT Tuebingen"
tmp$firstSumitter[tmp$firstSumitter == "Center for Pediatric Genomic Medicine,Children's Mercy Hospital and Clinics"] <- "Children's Mercy Hospital"
tmp$firstSumitter[tmp$firstSumitter == "EGL Genetic Diagnostics, Eurofins Clinical Diagnostics"] <- "Eurofins Clinical Diagnostics"
tmp$firstSumitter[tmp$firstSumitter == "Genetic Services Laboratory, University of Chicago"] <- "University of Chicago"
tmp$firstSumitter[tmp$firstSumitter == "Illumina Clinical Services Laboratory,Illumina"] <- "Illumina Clinical Laboratory"
tmp$firstSumitter[tmp$firstSumitter == "Laboratory for Molecular Medicine, Partners HealthCare Personalized Medicine"] <- "Partners HealthCare Person. Med."
tmp$firstSumitter[tmp$firstSumitter == "Quest Diagnostics Nichols Institute San Juan Capistrano"] <- "Quest Diagnostics"
tmp$firstSumitter[tmp$firstSumitter == "Women's Health and Genetics/Laboratory Corporation of America, LabCorp"] <- "LabCorp"

# creating df
submitterList = vector(mode = "list", length = length(topSubmitter))

for (i in 1:length(topSubmitter)){
  df = tmp[tmp$firstSumitter == topSubmitter[i],]
  tmpDf = unclass(table(df$firstSumitter, df$class))
  submitterList[[i]] <- tmpDf
}

tabularDF= data.frame(do.call(rbind, submitterList))


# making long names, shorter
rownames(tabularDF)[rownames(tabularDF) == "ARUP Laboratories, Molecular Genetics and Genomics,ARUP Laboratories"] <- "ARUP Laboratories"
rownames(tabularDF)[rownames(tabularDF) == "Biesecker Lab/Clinical Genomics Section,National Institutes of Health"] <- "Biesecker Lab_NIH"
rownames(tabularDF)[rownames(tabularDF) == "CeGaT Praxis fuer Humangenetik Tuebingen"] <- "CeGaT Tuebingen"
rownames(tabularDF)[rownames(tabularDF) == "Center for Pediatric Genomic Medicine,Children's Mercy Hospital and Clinics"] <- "Children's Mercy Hospital"
rownames(tabularDF)[rownames(tabularDF) == "EGL Genetic Diagnostics, Eurofins Clinical Diagnostics"] <- "Eurofins Clinical Diagnostics"
rownames(tabularDF)[rownames(tabularDF) == "Genetic Services Laboratory, University of Chicago"] <- "University of Chicago"
rownames(tabularDF)[rownames(tabularDF) == "Illumina Clinical Services Laboratory,Illumina"] <- "Illumina Clinical Laboratory"
rownames(tabularDF)[rownames(tabularDF) == "Laboratory for Molecular Medicine, Partners HealthCare Personalized Medicine"] <- "Partners HealthCare Person. Med."
rownames(tabularDF)[rownames(tabularDF) == "Quest Diagnostics Nichols Institute San Juan Capistrano"] <- "Quest Diagnostics"
rownames(tabularDF)[rownames(tabularDF) == "Women's Health and Genetics/Laboratory Corporation of America, LabCorp"] <- "LabCorp"


# convert the data as a table
dt <- as.table(as.matrix(tabularDF))
#
library("graphics")
png(filename = "~/clinvar/submitter_class_mosaic.png", width = 12, height = 10, units = "in", res = 300)
mosaicplot(dt, shade = TRUE, las=2,
           main = "")
dev.off()



### working on date of last evaluated to make sure about the label assignment as conflicting/no-conflicting

subInfo3$LastEvaluated[subInfo3$LastEvaluated == "-"] <- NA
subInfo3$LastEvaluated = as.Date(subInfo3$LastEvaluated, format = '%b %d,%y')

before2018 <- subset(subInfo3, LastEvaluated < as.Date("2021-05-01"))



############################# methods of collection #######################################
submitterDf = dplyr::left_join(submitterDf, subClin[, c(3,30)], by = "ID")
d = data.frame(table(submitterDf$CollectionMethod))
names(d) <- c("CollectionMethod", "total")
submitterDf = dplyr::left_join(submitterDf, d, by = "CollectionMethod")

# subsetting based on the total count
df = submitterDf[submitterDf$total >= 80,]
df$CollectionMethod = paste0(df$CollectionMethod, " [", df$total, "]")
# creating viz datafarme
tmp = data.frame(unclass(table(df$CollectionMethod, df$class.x)))
tmp$colMethod = rownames(tmp)
tmp$total = tmp$conflict + tmp$noConflict
# plotting with percent values
tmp$conflict = tmp$conflict/tmp$total
tmp$noConflict = tmp$noConflict/tmp$total

# rounding and *100
tmp$conflict = round(tmp$conflict * 100,0)
tmp$noConflict = round(tmp$noConflict * 100,0)

# # making other group and adding to the dataframe
# CollectionMethod = c("others", "others")
# class = c("conflict", "noConflict")
# 
# df3 = df[-which(df$CollectionMethod %in% freqMethod),]
# conflictCount = sum(df3$count[df3$class == "conflict"])
# noConflictCount = sum(df3$count[df3$class == "noConflict"])
# count = c(conflictCount,noConflictCount)
# 
# # adding to the main dataset
# tmp = rbind(df2, data.frame(CollectionMethod = CollectionMethod, class = class, count = count))
# 
# # convert counts into percent values
# num = seq(1,15, 2)
# relFreq = as.numeric()
# 
# for(i in num){
#   relFreq[i] = round(100 * (tmp$count[i]/(tmp$count[i] + tmp$count[i+1])),2)
#   relFreq[i+1] = round(100 * (tmp$count[i+1]/(tmp$count[i] + tmp$count[i+1])),2)
# }
# 
# tmp$relfreq = relFreq

# make a group value as negative
tmp$conflict = -tmp$conflict

# making viz dataframe
plotDf = reshape2::melt(tmp[-4], id.var = "colMethod")

# X Axis Breaks and Labels 
brks <- seq(-100, 100, 10)
brks
lbls <- c('100','90','80','70','60','50','40','30','20','10','0','10', '20','30','40','50','60','70','80','90','100')

png(filename = "~/clinvar/collectionMethod_class_2.png", width = 12, height = 9, units = "in", res = 300)

ggplot(plotDf, aes(x = reorder(colMethod, -abs(value)), y = value, fill = variable)) +   # Fill column
  geom_bar(stat = "identity", width = .6) +   # draw the bars
  scale_y_continuous(breaks = brks, labels = lbls) + # Labels
  coord_flip() +  # Flip axes
  labs(title="Relative frequency for conflict/no-conflict submission") +
  xlab('CollectionMethod [total number of submissions]') +
  ylab('%Submission') +
  theme_tufte() +  # Tufte theme from ggthemes
  theme(plot.title = element_text(hjust = .5),
        axis.ticks = element_blank(),
        text=element_text(size=15)) +   # Centre plot title
  scale_fill_manual(values=c('#999999','#E69F00')) +  # Color palette
  guides(fill=guide_legend("Class")) 

dev.off()

# mosaic plot

# plot
png(filename = "~/clinvar/class_mosaic_2.png", width = 12, height = 10, units = "in", res = 300)
mosaicplot(table(df$CollectionMethod,df$class), shade = TRUE, las=2,main = "")
dev.off()


# final dataset preparation
df = submitterDf[,-6]
finClin$ID = as.factor(finClin$ID)

finClin = dplyr::left_join(finClin, df, by = "ID")

d = data.frame(table(finClin$CollectionMethod))
colmethod100 = d$Var1[d$Freq > 10]
colmethod = ifelse(finClin$CollectionMethod %in% colmethod100, finClin$CollectionMethod, "other")

finClin$CollectionMethod = colmethod

saveRDS(finClin, "~/clinvar/finClin.rds")

# up to now : "SYMBOL", "Allele", "CLNVC", "Consequence", "IMPACT", "relativeLocationRatio",

#"ExIntron", "Protein_position", "Codons", "TSL", "Denisova",SIFT, "PolyPhen", 
# "firstSumitter", "submitterNo", "CollectionMethod"
# plus c(important_var_pathRank, important_var_AS_AC_AN_nhomalt, important_var_ConsOrgRankScore)

################################################################################
###########################  Dataset Preparation #################################

catfeature = c("class", "ID", "SYMBOL", "Allele", "CLNVC", "Consequence", "IMPACT", "relativeLocationRatio",
               "ExIntron", "Protein_position", "Codons", "TSL", "Denisova","SIFT", "PolyPhen",
               "firstSumitter", "submitterNo", "CollectionMethod")

FEATURES = c(catfeature,important_var_pathRank, important_var_AS_AC_AN_nhomalt, important_var_ConsOrgRankScore)

finalClin = finClin[, FEATURES]
finalClin[finalClin == "-"] <- NA
finalClin$firstSumitter = gsub("[[:punct:]]", "", finalClin$firstSumitter)
finalClin$firstSumitter = gsub(" ", "", finalClin$firstSumitter)
finalClin$firstSumitter = gsub("__", "_", finalClin$firstSumitter)

d = data.frame(table(finalClin$firstSumitter))
topSub = d$Var1[d$Freq > 10]
finalClin$firstSumitter = ifelse(finalClin$firstSumitter %in% topSub, finalClin$firstSumitter, "other")
finalClin$firstSumitter = substr(finalClin$firstSumitter,1,35)

# correcting names
# CentredeGÃnÃtiqueMolÃculaireetChrom
finalClin$firstSumitter[finalClin$firstSumitter == "CentredeGÃnÃtiqueMolÃculaireetChrom"] <- "CentredeGAnAtiqueMolAculaireetChrom"
#GÃnÃtiquedesMaladiesduDÃveloppement
finalClin$firstSumitter[finalClin$firstSumitter == "GÃnÃtiquedesMaladiesduDÃveloppement"] <- "GAnAtiquedesMaladiesduDAveloppement"
# ThreeDMedClinicalLaboratoryInc"
finalClin$firstSumitter[finalClin$firstSumitter == "3DMedClinicalLaboratoryInc"] <- "ThreeDMedClinicalLaboratoryInc"

colnames(finalClin)[16] <- "firstSubmitter" 


# collection method
finalClin$CollectionMethod[finalClin$CollectionMethod == "case-control"] <- "CaseControl"
finalClin$CollectionMethod[finalClin$CollectionMethod == "clinical testing"] <-"ClinicalTesting"
finalClin$CollectionMethod[finalClin$CollectionMethod == "clinical testing;curation"] <-"ClinicalTesting_Curation"
finalClin$CollectionMethod[finalClin$CollectionMethod == "clinical testing;in vitro"] <-"ClinicalTesting_inVitro"
finalClin$CollectionMethod[finalClin$CollectionMethod == "clinical testing;provider interpretation"] <-"ClinicalTesting_ProviderInterpretation"
finalClin$CollectionMethod[finalClin$CollectionMethod == "curation;literature only"] <-"Curation_LiteratureOnly"
finalClin$CollectionMethod[finalClin$CollectionMethod == "in vitro"] <-"inVitro"
finalClin$CollectionMethod[finalClin$CollectionMethod == "literature only"] <-"LiteratureOnly"
finalClin$CollectionMethod[finalClin$CollectionMethod == "provider interpretation"] <-"ProviderInterpretation"
finalClin$CollectionMethod[finalClin$CollectionMethod == "reference population"] <-"ReferencePopulation"

#
finalClin$class = ifelse(finalClin$class == "noConflict", 0, 1)
finalClin$AF_TGP[finalClin$AF_TGP == "."] <- NA
#
saveRDS(finalClin, file = "~/clinvar/finalClin.RDS")
write.csv(finalClin, "~/clinvar/finalClin.csv")


### dissecting variant type for PKHD1 gene
subInfo1 <- data.table::fread("~/clinvar/submission_summary.txt")

acmgTerms = c("Benign","Likely benign","Uncertain significance","Likely pathogenic","Pathogenic")

####__________________________ PKHD1
pk = subInfo1[subInfo1$`#VariationID` %in% subClin$ID[subClin$SYMBOL == "PKHD1" &
                                                        subClin$class == "conflict"],]
# subsetting
subPK = pk[, c(1:2)]
castDF = reshape2::dcast(subPK, `#VariationID` ~ ClinicalSignificance,
                         value.var = "ClinicalSignificance")
castDF = castDF[,acmgTerms]
for(col in 1:5){
  castDF[,col] <- ifelse(castDF[,col] == 0, 0,1)
}
castDF$tot = as.character(paste0(castDF[,1],castDF[,2],castDF[,3],castDF[,4],castDF[,5]))

pkhd1=table(castDF$tot)[order(table(castDF$tot), decreasing = T)]

####__________________________ TTN
pk = subInfo1[subInfo1$`#VariationID` %in% subClin$ID[subClin$SYMBOL == "TTN" &
                                                        subClin$class == "conflict"],]
# subsetting
subPK = pk[, c(1:2)]
castDF = reshape2::dcast(subPK, `#VariationID` ~ ClinicalSignificance,
                         value.var = "ClinicalSignificance")
castDF = castDF[,acmgTerms]
for(col in 1:5){
  castDF[,col] <- ifelse(castDF[,col] == 0, 0,1)
}
castDF$tot = as.character(paste0(castDF[,1],castDF[,2],castDF[,3],castDF[,4],castDF[,5]))

ttn=table(castDF$tot)[order(table(castDF$tot), decreasing = T)]

####__________________________ COL6A3
pk = subInfo1[subInfo1$`#VariationID` %in% subClin$ID[subClin$SYMBOL == "COL6A3" &
                                                        subClin$class == "conflict"],]
# subsetting
subPK = pk[, c(1:2)]
castDF = reshape2::dcast(subPK, `#VariationID` ~ ClinicalSignificance,
                         value.var = "ClinicalSignificance")
castDF = castDF[,acmgTerms]
for(col in 1:5){
  castDF[,col] <- ifelse(castDF[,col] == 0, 0,1)
}
castDF$tot = as.character(paste0(castDF[,1],castDF[,2],castDF[,3],castDF[,4],castDF[,5]))

col6a3=table(castDF$tot)[order(table(castDF$tot), decreasing = T)]

## for loop for all

ge = c("PKHD1", "TTN", "COL6A3", "SYNE1", "DYSF", "PLEC", "ADGRV1", "APOB", "CFTR", "DMD", "NEB", "RYR1", "TSC1", "MYH7", "TSC2", "LDLR", "USH2A", "PKD1", "SCN5A", "FBN1", "KCNH2", "TP53", "RET", "MYBPC3", "RYR2", "BRCA2", "BRCA1", "PTCH1", "STK11", "DSP", "HPDL", "MLH1", "PMS2", "MSH2", "NBN", "NF1", "POLD1", "ATM", "CHEK2", "BRIP1", "CDH1", "MSH6", "PALB2", "POLE", "APC", "BARD1", "ENSG00000258864", "SMARCA4", "DICER1", "RAD50")

glist <- vector(mode = "list", length = 50)
for(g in c(1:5, 7:36, 38:46, 48)){
  pk = subInfo1[subInfo1$`#VariationID` %in% subClin$ID[subClin$SYMBOL == ge[g] &
                                                          subClin$class == "conflict"],]
  # subsetting
  subPK = pk[, c(1:2)]
  
  castDF = reshape2::dcast(subPK, `#VariationID` ~ ClinicalSignificance, value.var = "ClinicalSignificance")
  
  castDF = castDF[,c("Benign","Likely benign","Uncertain significance","Likely pathogenic","Pathogenic")]
  for(col in 1:5){
    castDF[,col] <- ifelse(castDF[,col] == 0, 0,1)
  }
  castDF$tot = as.character(paste0(castDF[,1],castDF[,2],castDF[,3],castDF[,4],castDF[,5]))
  glist[[g]] = table(castDF$tot)[order(table(castDF$tot), decreasing = T)]
  names(glist)[g] = ge[g]
  print(ge[g])
  print(g)
}

# adding missing genes to the list
g = 6
pk = subInfo1[subInfo1$`#VariationID` %in% subClin$ID[subClin$SYMBOL == ge[g] &
                                                        subClin$class == "conflict"],]
# subsetting
subPK = pk[, c(1:2)]
castDF = reshape2::dcast(subPK, `#VariationID` ~ ClinicalSignificance, value.var = "ClinicalSignificance")
missingCol = c("Benign","Likely benign","Uncertain significance","Likely pathogenic","Pathogenic")[-which(c("Benign","Likely benign","Uncertain significance","Likely pathogenic","Pathogenic") %in%
                                                                                                     colnames(castDF))]
castDF$Likely_pathogenic = 0
colnames(castDF)[6] = "Likely pathogenic"
castDF = castDF[,c("Benign","Likely benign","Uncertain significance","Likely pathogenic","Pathogenic")]
for(col in 1:5){
  castDF[,col] <- ifelse(castDF[,col] == 0, 0,1)
}
castDF$tot = as.character(paste0(castDF[,1],castDF[,2],castDF[,3],castDF[,4],castDF[,5]))
glist[[g]] = table(castDF$tot)[order(table(castDF$tot), decreasing = T)]
names(glist)[g] = ge[g]
print(ge[g])

#########
g = 37
pk = subInfo1[subInfo1$`#VariationID` %in% subClin$ID[subClin$SYMBOL == ge[g] &
                                                        subClin$class == "conflict"],]
# subsetting
subPK = pk[, c(1:2)]
castDF = reshape2::dcast(subPK, `#VariationID` ~ ClinicalSignificance, value.var = "ClinicalSignificance")
missingCol = c("Benign","Likely benign","Uncertain significance","Likely pathogenic","Pathogenic")[-which(c("Benign","Likely benign","Uncertain significance","Likely pathogenic","Pathogenic") %in%
                                                                                                            colnames(castDF))]
castDF$Likely_pathogenic = 0
castDF$Pathogenic = 0
colnames(castDF)[6] = "Likely pathogenic"
castDF = castDF[,c("Benign","Likely benign","Uncertain significance","Likely pathogenic","Pathogenic")]
for(col in 1:5){
  castDF[,col] <- ifelse(castDF[,col] == 0, 0,1)
}
castDF$tot = as.character(paste0(castDF[,1],castDF[,2],castDF[,3],castDF[,4],castDF[,5]))
glist[[g]] = table(castDF$tot)[order(table(castDF$tot), decreasing = T)]
names(glist)[g] = ge[g]
print(ge[g])

#########
g = 48
pk = subInfo1[subInfo1$`#VariationID` %in% subClin$ID[subClin$SYMBOL == ge[g] &
                                                        subClin$class == "conflict"],]
# subsetting
subPK = pk[, c(1:2)]
castDF = reshape2::dcast(subPK, `#VariationID` ~ ClinicalSignificance, value.var = "ClinicalSignificance")
missingCol = c("Benign","Likely benign","Uncertain significance","Likely pathogenic","Pathogenic")[-which(c("Benign","Likely benign","Uncertain significance","Likely pathogenic","Pathogenic") %in%
                                                                                                            colnames(castDF))]

castDF = castDF[,c("Benign","Likely benign","Uncertain significance","Likely pathogenic","Pathogenic")]
for(col in 1:5){
  castDF[,col] <- ifelse(castDF[,col] == 0, 0,1)
}
castDF$tot = as.character(paste0(castDF[,1],castDF[,2],castDF[,3],castDF[,4],castDF[,5]))
glist[[g]] = table(castDF$tot)[order(table(castDF$tot), decreasing = T)]
names(glist)[g] = ge[g]
print(ge[g])
#########
g = 49
pk = subInfo1[subInfo1$`#VariationID` %in% subClin$ID[subClin$SYMBOL == ge[g] &
                                                        subClin$class == "conflict"],]
# subsetting
subPK = pk[, c(1:2)]
castDF = reshape2::dcast(subPK, `#VariationID` ~ ClinicalSignificance, value.var = "ClinicalSignificance")
missingCol = c("Benign","Likely benign","Uncertain significance","Likely pathogenic","Pathogenic")[-which(c("Benign","Likely benign","Uncertain significance","Likely pathogenic","Pathogenic") %in%
                                                                                                            colnames(castDF))]

castDF$Likely_pathogenic = 0
colnames(castDF)[7] = "Likely pathogenic"

castDF = castDF[,c("Benign","Likely benign","Uncertain significance","Likely pathogenic","Pathogenic")]
for(col in 1:5){
  castDF[,col] <- ifelse(castDF[,col] == 0, 0,1)
}
castDF$tot = as.character(paste0(castDF[,1],castDF[,2],castDF[,3],castDF[,4],castDF[,5]))
glist[[g]] = table(castDF$tot)[order(table(castDF$tot), decreasing = T)]
names(glist)[g] = ge[g]
print(ge[g])

#########
g = 50
pk = subInfo1[subInfo1$`#VariationID` %in% subClin$ID[subClin$SYMBOL == ge[g] &
                                                        subClin$class == "conflict"],]
# subsetting
subPK = pk[, c(1:2)]
castDF = reshape2::dcast(subPK, `#VariationID` ~ ClinicalSignificance, value.var = "ClinicalSignificance")
missingCol = c("Benign","Likely benign","Uncertain significance","Likely pathogenic","Pathogenic")[-which(c("Benign","Likely benign","Uncertain significance","Likely pathogenic","Pathogenic") %in%
                                                                                                            colnames(castDF))]

castDF$Pathogenic = 0

castDF = castDF[,c("Benign","Likely benign","Uncertain significance","Likely pathogenic","Pathogenic")]
for(col in 1:5){
  castDF[,col] <- ifelse(castDF[,col] == 0, 0,1)
}
castDF$tot = as.character(paste0(castDF[,1],castDF[,2],castDF[,3],castDF[,4],castDF[,5]))
glist[[g]] = table(castDF$tot)[order(table(castDF$tot), decreasing = T)]
names(glist)[g] = ge[g]
print(ge[g])

#########
subClin$SYMBOL = ifelse(subClin$SYMBOL == "-", subClin$Gene, subClin$SYMBOL)

g = 47
pk = subInfo1[subInfo1$`#VariationID` %in% subClin$ID[subClin$SYMBOL == ge[g] &
                                                        subClin$class == "conflict"],]
# subsetting
subPK = pk[, c(1:2)]
castDF = reshape2::dcast(subPK, `#VariationID` ~ ClinicalSignificance, value.var = "ClinicalSignificance")
missingCol = c("Benign","Likely benign","Uncertain significance","Likely pathogenic","Pathogenic")[-which(c("Benign","Likely benign","Uncertain significance","Likely pathogenic","Pathogenic") %in%
                                                                                                            colnames(castDF))]


castDF = castDF[,c("Benign","Likely benign","Uncertain significance","Likely pathogenic","Pathogenic")]
for(col in 1:5){
  castDF[,col] <- ifelse(castDF[,col] == 0, 0,1)
}
castDF$tot = as.character(paste0(castDF[,1],castDF[,2],castDF[,3],castDF[,4],castDF[,5]))
glist[[g]] = table(castDF$tot)[order(table(castDF$tot), decreasing = T)]
names(glist)[g] = ge[g]
print(ge[g])

# making a dataframe
df = data.frame()

for(i in 1:length(glist)){
  d = data.frame(t(data.frame(as.list(glist[[i]]), check.names = F)))
  names(d)[1] <- "count"
  d$gene = names(glist)[i]
  d$conflict = as.character(rownames(d))
  df = rbind(df,d)
}

dfWide = reshape2::dcast(df, gene ~ conflict, value.var = "count")

dfWide[is.na(dfWide)] <- 0
rownames(dfWide) = dfWide$gene
dfWide = dfWide[,-1]

# making more abstract representation

vec = c("01100" = "Benign_VUS",
       "11100"= "Benign_VUS",
       "00110" = 'VUS_Pathogenic',
       "00111" = 'VUS_Pathogenic',
       "10100" = "Benign_VUS",
       "00101" = 'VUS_Pathogenic',
       "01110" = 'Benign_VUS_Pathogenic',
       "01111" = 'Benign_VUS_Pathogenic',
       "11110" = 'Benign_VUS_Pathogenic',
       "01101" = 'Benign_VUS_Pathogenic',
       "11101"= 'Benign_VUS_Pathogenic',
       "10110" = 'Benign_VUS_Pathogenic',
       "11000" = 'Benign_Bening',
       "01001" = 'Benign_Pathogenic',
       "01011" = 'Benign_Pathogenic',
       "10001" = 'Benign_Pathogenic',
       "10011" = 'Benign_Pathogenic',
       "11111" = 'Benign_VUS_Pathogenic',
       "10111" = 'Benign_VUS_Pathogenic',
       "11011" = 'Benign_Pathogenic',
       "01010" = 'Benign_Pathogenic',
       "11010" = 'Benign_Pathogenic',
       "11001" = 'Benign_Pathogenic',
       "10101" = 'Benign_VUS_Pathogenic',
       "01000" = "Benign"
)

df$conflict_chr = stringr::str_replace_all(df$conflict, vec)
# collapse rwos and add numer of variants based on the conflict group

df2 <- df %>% 
  select(gene,conflict_chr,count) %>% 
  na.omit() %>% 
  group_by(gene, conflict_chr) %>% 
  summarize(stat=sum(count))

df2 = df2[-which(df2$conflict_chr %in% c("Benign", "Benign_Bening")),]







df2Wide = reshape2::dcast(df2, gene ~ conflict_chr, value.var = "stat")

df2Wide[is.na(df2Wide)] <- 0
rownames(df2Wide) = df2Wide$gene
df2Wide = df2Wide[,-1]


# 
library(pheatmap)
pheatmap(scale(df2Wide))

# working on the top 10


png(filename="~/clinvar/conflict_type.png", width = 30, height = 15, units = "in", res = 300)
op <- par(mfrow=c(4,3))
          barplot(glist[[1]], main = names(glist)[1])
          barplot(glist[[2]], main = names(glist)[2])
          barplot(glist[[3]], main = names(glist)[3])
          barplot(glist[[4]], main = names(glist)[4])
          barplot(glist[[5]], main = names(glist)[5])
          barplot(glist[[6]], main = names(glist)[6])
          barplot(glist[[7]], main = names(glist)[7])
          barplot(glist[[8]], main = names(glist)[8])
          barplot(glist[[10]], main = names(glist)[9])
          barplot(glist[[10]], main = names(glist)[10])
          
par(op)
dev.off()


######## validation dataset

clinvar2022 = data.table::fread("~/clinvar/sandbox/clinvar_2022.txt")

d <- data.frame(table(clinvar2022$CLNREVSTAT))

subDf <- clinvar2022[clinvar2022$CLNREVSTAT %in% c("criteria_provided,_conflicting_interpretations",
                                 "criteria_provided,_multiple_submitters,_no_conflicts"), ]

subDf <- subDf[-which(subDf$ID %in% subClin$ID), ]

#
table(subDf$CLNREVSTAT)
# saving file and getting annotations
#write.csv(subDf, "~/clinvar/sandbox/val_df.csv")
# adding clas group to the subdf
subDf$class = ifelse(subDf$CLNREVSTAT == "criteria_provided,_conflicting_interpretations", "conflict", "noConflict")

# adding submitter info
subInfo1 <- data.table::fread("~/clinvar/sandbox/submission_summary_2022_02.txt")


# distribution of submitter for variants
submitterCount = data.frame(table(subInfo1$`#VariationID`))
colnames(submitterCount) = c("ID", "Subcount")
submitterCount$ID = as.factor(submitterCount$ID)
#overlap with Subclin
table(submitterCount$ID %in% subDf$ID)

# joining two datasets
subDf$ID  = as.factor(subDf$ID)
subCount = dplyr::left_join(submitterCount[which(submitterCount$ID %in% subDf$ID),], 
                            subDf[, c(3,30)])


# making a df for each variant ID , first submitter, date and count of submission and other variables

# filtering based on variant ids

subSum = subInfo1[subInfo1$`#VariationID` %in% subDf$ID,]

ids = unique(subSum$`#VariationID`)
varID = c()
submitterNo = c()
firstSumitter = c()
firstSumitterDate = c()
CollectionMethod = c()

for (i in 1:length(ids)){
  df = subSum[subSum$`#VariationID` == ids[i], ]
  varID[i] = ids[i]
  submitterNo [i] = dim(df)[1]
  df$DateLastEvaluated = sub(",", "", df$DateLastEvaluated)
  df$DateLastEvaluated = as.Date(df$DateLastEvaluated,format = "%B %d %Y")
  df = df[order(df$DateLastEvaluated),]
  firstSumitterDate[i] = as.character(df$DateLastEvaluated[1])
  firstSumitter[i] = df$Submitter[1]
  CollectionMethod[i] = df$CollectionMethod[1]
}

submitterDf = data.frame(ID = as.factor(varID), firstSumitter = firstSumitter,firstSumitterDate = firstSumitterDate,
                         submitterNo = submitterNo, CollectionMethod = CollectionMethod)

# 
tmp = data.frame(unclass(table(submitterDf$firstSumitter)))
tmp$submitter = rownames(tmp)
topSubmitter = tmp$submitter[tmp$unclass.table.submitterDf.firstSumitter.. >= 1000]

# A total of 507 centers contributed to submission of 187643  variants.
# top centers with highest number of submission

submitterDf = dplyr::left_join(submitterDf, subDf[, c(3,30)], by = "ID")

## adding vep annptation to the dataset
vep <- data.table::fread("~/clinvar/sandbox/vep_clinvar_2022.txt", check.names = F)
vep = data.frame(vep)

names(vep)[1] <- "ID"
vep$ID = as.factor(vep$ID)
#rm(vep)

# joining datasets
subDf = dplyr::left_join(subDf, vep, by = "ID")
subDf = dplyr::left_join(subDf, submitterDf, by = "ID")
#save object
#saveRDS(subDf, "~/clinvar/sandbox/subClin2022.rds")
########################################################################################
# preparing the dataset
subClin2022 = readRDS("~/clinvar/sandbox/subClin2022.rds")

# selecting columns present in the final data set-that used for model building 
finClin = read.csv("~/clinvar/finalClin_3.csv")

# 
table(colnames(subClin2022) %in% colnames(finClin))

colnames(finClin)[-which(colnames(finClin) %in% colnames(subClin2022))]

# correcting names
names(subClin2022)[704] <- "class"
names(subClin2022)[700] <- "firstSubmitter"
#
colnames(finClin)[-which(colnames(finClin) %in% colnames(subClin2022))]
#[1] "relativeLocationRatio" "ExIntron"              "EVE_scores_ASM"  
# these three names will be added to the dataset later

idx = colnames(subClin2022) %in% colnames(finClin)

subClin2022 = data.frame(subClin2022, check.names = F)
subClin2022 = subClin2022[, idx]
# converting . and - to NA 
subClin2022[subClin2022 == "."] <- NA
subClin2022[subClin2022 == "-"] <- NA

# Allele frequency, gene type and .. are good to go 

######################### variant location in exons ##########################################
subClin = readRDS("~/clinvar/sandbox/subClin2022.rds")

View(table(subClin$EXON))

table(subClin$EXON == "-")

# FALSE   TRUE 
# 153639  33854

tmp = subClin[, c(30, 40, 41)]
#tmp = tmp[!is.na(tmp$EXON),]


tmp$ratio = tmp$EXON
tmp$ratio = ifelse(tmp$ratio != "-", tmp$ratio, tmp$INTRON)
tmp = cbind(tmp, stringr::str_split_fixed(tmp$ratio, "/", 2))
names(tmp)[c(5:6)] <- c("nExon", "totalExon")
tmp$nratio = as.numeric(tmp$nExon)/as.numeric(tmp$totalExon)
##@@
tmp$nIntron = tmp$nExon
tmp$nIntron = as.factor(ifelse(tmp$INTRON == "-", tmp$nExon, "0"))
tmp$nIntron[tmp$nIntron == "-"] <- NA

# two variable need to be added to finClin
subClin2022$relativeLocationRatio = tmp$nratio
subClin2022$ExIntron = tmp$nIntron

############################### pathogenicity  score ############################################
### PUtting EVE result together
setwd("~/clinvar/variant_files/")

flist = list.files()

eveGene = gsub(".csv", "", flist)

uniprotID = read.table("~/clinvar/uniprotConversionTable.tab", fill = T, header = T)

ids = unique(subClin$GENEINFO)
ids = data.frame(stringr::str_split_fixed(ids, ":", 2))
names(ids) = c("symbol", "ENTREZ_ID")

# correcting ENTREZ numbers
entz = strsplit(ids$ENTREZ_ID,split='|',fixed=TRUE)
length(entz)
entz_num = c()
for (i in 1:length(entz)){
  entz_num[i] = entz[[i]][1]
}

ids$ENTREZ_ID = as.numeric(entz_num)


# joing two IDS

geneIDs = dplyr::left_join(ids, uniprotID)
geneIDs = geneIDs[!duplicated(geneIDs$symbol),]

# to see how many of genes with NAs in uniprot in genID dataset,
# Are really protein coding:

library(biomaRt)
mart <- useDataset("hsapiens_gene_ensembl", useMart("ensembl"))
genes <- getBM(attributes= c("ensembl_gene_id","entrezgene_id","hgnc_symbol", "gene_biotype"), 
               mart= mart)

protGene = genes[genes$gene_biotype == "protein_coding",]

table(geneIDs$ENTREZ_ID[is.na(geneIDs$uniprot)] %in% protGene$entrezgene_id)

# protein codng gene not in EVE:
notInEVE = geneIDs$ENTREZ_ID[is.na(geneIDs$uniprot)][geneIDs$ENTREZ_ID[is.na(geneIDs$uniprot)] %in% protGene$entrezgene_id]

notInEVE = genes[genes$entrezgene_id %in% notInEVE,]

#
ss = subClin[grep("synonymous", subClin$Consequence),]
ss$Amino_acids = paste0(ss$Amino_acids, "/", ss$Amino_acids)
ss = rbind(ss, subClin[grep("missense", subClin$Consequence)])
unik = ss$GENEINFO[!duplicated(ss$GENEINFO)]
eveList = list()

for (gene in 1:length(unik)){
  #selecting per gene
  subgen = ss[ss$GENEINFO == unik[gene],]
  # converting geneinfo to ids and merging with subgen
  ids = unique(subgen$GENEINFO)
  ids = data.frame(stringr::str_split_fixed(ids, ":", 2))
  names(ids) = c("symbol", "ENTREZ_ID")
  # correcting ENTREZ numbers
  entz = strsplit(ids$ENTREZ_ID,split='|',fixed=TRUE)
  entz_num = c()
  for (i in 1:length(entz)){
    entz_num[i] = entz[[i]][1]
  }
  ids$ENTREZ_ID = as.numeric(entz_num)
  geneIDs = dplyr::left_join(ids, uniprotID, by = "ENTREZ_ID")
  geneIDs = geneIDs[!duplicated(geneIDs$symbol),]
  subgen$uniprotID = geneIDs$uniprot
  subgen$ENTREZ_ID = geneIDs$ENTREZ_ID
  subgen$symbol = geneIDs$symbol
  #making mutation column in this dataset
  mut = data.frame(X1 = subgen$Protein_position, X2 = subgen$Amino_acids)
  mut$symbol = subgen$symbol
  mut$mutation = paste0(mut$symbol,"_", mut$X1, "_", mut$X2 )
  subgen$mutation  = mut$mutation
  if(!is.na(geneIDs$uniprot)){
    eveFrame = read.csv(paste0("~/clinvar/variant_files/", geneIDs$uniprot, ".csv"))
    eveFrame$uniprotID = geneIDs$uniprot
    eveFrame$symbol = geneIDs$symbol
    # substitution of one letter aa to three letter code
    eveFrame$mutation = paste0(eveFrame$symbol, "_",
                               eveFrame$position, "_",
                               eveFrame$wt_aa, "/",
                               eveFrame$mt_aa)
  }
  subgen = dplyr::left_join(subgen, eveFrame, by = "mutation")
  subgen = subgen[, c(3,719)]
  eveList[[gene]] = subgen
  print(paste0(gene, " out of", length(unik),", gene = ", 
               eveFrame$symbol[1], " [" , round((100*gene/length(unik)),2), "% completed]"))
}

eveFinalDataframe = do.call(rbind, eveList)


# adding eve data to subclin data
subClin2022 = dplyr::left_join(subClin2022, eveFinalDataframe, by = "ID")
#

# converting SIFT/PolyPhen to numeric score and prediction
df = cbind(t(data.frame(strsplit(subClin2022$SIFT, "\\("))), t(data.frame(strsplit(subClin2022$PolyPhen, "\\("))))
colnames(df) <- c("SIFT_Pred", "SIFT_score", "PolyPhen_Pred", "PolyPhen_score")
df = data.frame(df)
df$SIFT_score = as.numeric(sub("\\)", "", df$SIFT_score))
df$PolyPhen_score = as.numeric(sub("\\)", "", df$PolyPhen_score))
rownames(df) <- NULL

subClin2022$SIFT = df$SIFT_score
subClin2022$PolyPhen = df$PolyPhen_score

# checking the first author format
val_df = subClin2022
fc = finClin
View(table(fc$firstSubmitter))
# centers with <10 submission , encoded as "other"
# samething shold be applied for the val_df
# Before doing so, spaces have to be removed from the names

val_df$firstSubmitter = gsub(" ", "", val_df$firstSubmitter)
table(val_df$firstSubmitter %in% fc$firstSubmitter)

# recoding those centers with lower than 10 submission into other + converting centers present in the val df but not i the fc to other:
tmp = data.frame(table(val_df$firstSubmitter))
otherCent=  tmp$Var1[tmp$Freq <= 10]
val_df$firstSubmitter = ifelse(val_df$firstSubmitter %in% otherCent, "other", val_df$firstSubmitter)
# checking to see how things get improved 
table(val_df$firstSubmitter %in% fc$firstSubmitter)
#  converting centers present in the val df but not i the fc to other:
unCent = unique(val_df$firstSubmitter[-which(val_df$firstSubmitter %in% fc$firstSubmitter)])

val_df$firstSubmitter = ifelse(val_df$firstSubmitter %in% unCent, "other", val_df$firstSubmitter)

# save val_RDS
#saveRDS(val_df, "~/clinvar/sandbox/val_df.RDS")

# Consequence
val_df = val_df[which(val_df$Consequence %in% fc$Consequence),]
#
table(val_df$Consequence %in% fc$Consequence)

# SYMBOL remain unchanges however not all from val_df are in fc
val_df = val_df[val_df$SYMBOL %in% fc$SYMBOL,]
# Collection method
val_df$CollectionMethod[val_df$CollectionMethod == "case-control"] <- "CaseControl"
val_df$CollectionMethod[val_df$CollectionMethod == "clinical testing"] <- "ClinicalTesting"
val_df$CollectionMethod[val_df$CollectionMethod == "clinical testing;curation"] <- "ClinicalTesting_Curation"
val_df$CollectionMethod[val_df$CollectionMethod == "clinical testing;in vitro"] <- "ClinicalTesting_inVitro"
val_df$CollectionMethod[val_df$CollectionMethod == "clinical testing;provider interpretation"] <- "ClinicalTesting_ProviderInterpretation"
val_df$CollectionMethod[val_df$CollectionMethod == "in vitro"] <- "inVitro"
val_df$CollectionMethod[val_df$CollectionMethod == "literature only"] <- "LiteratureOnly"
val_df$CollectionMethod[val_df$CollectionMethod == "literature only;research"] <- "other"
val_df$CollectionMethod[val_df$CollectionMethod == "provider interpretation"] <- "ProviderInterpretation"
val_df$CollectionMethod[val_df$CollectionMethod == "reference population"] <- "ReferencePopulation"
val_df$CollectionMethod[val_df$CollectionMethod == "research"] <- "research"

# convet class into code
val_df$class = ifelse(val_df$class == "noConflict", "0", "1")

# Denisova
val_df$Denisova[val_df$Denisova == "C/T"] <- "other"
val_df$Denisova[val_df$Denisova == "C/G"] <- "other"
val_df$Denisova[val_df$Denisova == "A/T"] <- "other"
val_df$Denisova[val_df$Denisova == "A/G"] <- "other"
val_df$Denisova[val_df$Denisova == "A/C"] <- ""

# checking
table(val_df$Denisova)

# is the same
table(fc$Denisova)

# Allele
table(val_df$Allele %in% fc$Allele)

# defining a new column i the orginal datfraeme and thenm adding that to te final selection
tmp = subClin[, c(30,4,5)]
names(tmp)[1] <- "class"
tmp$allele = paste0(tmp$REF, ">", tmp$ALT)

tmp$REFlenCat = ifelse(nchar(tmp$REF) ==1, "1", ">1")
tmp$ALTlenCat = ifelse(nchar(tmp$ALT) ==1, "1", ">1")
# coding REF multi_nts and ALT multi_nts
tmp$allele = ifelse(tmp$REFlenCat != 1, "REFmulti_nts", 
                    ifelse(tmp$ALTlenCat != 1,"ALTmulti_nts", tmp$allele))
tmp$ID = subClin$ID
tmp = tmp[, c(7,4)]
#
val_df = dplyr::left_join(val_df, tmp, by = "ID")
val_df$Allele = val_df$allele 
val_df = val_df[, -51]
#
saveRDS(val_df, "~/clinvar/sandbox/val_df.RDS")
write.csv(val_df, "~/clinvar/val_df.csv", row.names = F)



# Preparing single sublitter data
# 
clinvar2022 = data.table::fread("~/clinvar/sandbox/clinvar_2022.txt")

table(clinvar2022$CLNREVSTAT)

ss = clinvar2022[clinvar2022$CLNREVSTAT == "criteria_provided,_single_submitter",]
ss = ss[, c(1:8)]
ss[, 6] <- ss[, 7]
ss[, 8] <- ss[, 7]

# removing variants with extra large allele
alLength = nchar(ss$REF) + nchar(ss$ALT)
alLength = alLength < 15
ss = ss[alLength,]
write.table(ss, "~/clinvar/sandbox/ss.csv", sep="\t",  col.names=FALSE, row.names = F, quote = F)


# Adding first submitter and saubmitter number
ss = clinvar2022[clinvar2022$CLNREVSTAT == "criteria_provided,_single_submitter",]
# adding submitter info
subInfo1 <- data.table::fread("~/clinvar/sandbox/submission_summary_2022_02.txt")


# distribution of submitter for variants
submitterCount = data.frame(table(subInfo1$`#VariationID`))
colnames(submitterCount) = c("ID", "Subcount")
submitterCount$ID = as.factor(submitterCount$ID)
#overlap with Subclin
table(submitterCount$ID %in% ss$ID)
submitterCount = submitterCount[submitterCount$ID %in% ss$ID,]

# joining two datasets
ss$ID  = as.factor(ss$ID)
ss = dplyr::left_join(ss, submitterCount)
ss = ss[,-30]

# making a df for each variant ID , first submitter, date and count of submission and other variables

# filtering based on variant ids

subSum = subInfo1[subInfo1$`#VariationID` %in% ss$ID,]

ids = unique(subSum$`#VariationID`)
varID = c()
submitterNo = c()
firstSumitter = c()
firstSumitterDate = c()
CollectionMethod = c()

for (i in 1:length(ids)){
  df = subSum[subSum$`#VariationID` == ids[i], ]
  varID[i] = ids[i]
  submitterNo [i] = dim(df)[1]
  df$DateLastEvaluated = sub(",", "", df$DateLastEvaluated)
  df$DateLastEvaluated = as.Date(df$DateLastEvaluated,format = "%B %d %Y")
  df = df[order(df$DateLastEvaluated),]
  firstSumitterDate[i] = as.character(df$DateLastEvaluated[1])
  firstSumitter[i] = df$Submitter[1]
  CollectionMethod[i] = df$CollectionMethod[1]
  print(paste0(i, " out of 894191" , "(", 100*(i/894191), "% )"))
}

submitterDf = data.frame(ID = as.factor(varID), firstSumitter = firstSumitter,firstSumitterDate = firstSumitterDate,
                         submitterNo = submitterNo, CollectionMethod = CollectionMethod)
submitterDf = readRDS("~/clinvar/sandbox/ss_submitterDF.RDS")
# 
tmp = data.frame(unclass(table(submitterDf$firstSumitter)))
tmp$submitter = rownames(tmp)
topSubmitter = tmp$submitter[tmp$unclass.table.submitterDf.firstSumitter.. >= 1000]

# A total of 507 centers contributed to submission of 187643  variants.
# top centers with highest number of submission

ss = dplyr::left_join(ss, submitterDf, by = "ID")

## adding vep annptation to the dataset
vep <- data.table::fread("~/clinvar/sandbox/single_submitter.vep.txt", check.names = F)
vep = data.frame(vep)

names(vep)[1] <- "ID"
vep$ID = as.factor(vep$ID)
#rm(vep)

# joining datasets
ss = dplyr::left_join(ss, vep, by = "ID")

# selecting columns present in the final data set-that used for model building 
finClin = read.csv("~/clinvar/finalClin_3.csv")

#
names(ss)[30] <- "firstSubmitter" 
table(colnames(ss) %in% colnames(finClin))

colnames(finClin)[-which(colnames(finClin) %in% colnames(ss))]

#[1] "relativeLocationRatio" "ExIntron"              "EVE_scores_ASM"  
# these three names will be added to the dataset later

idx = colnames(subClin2022) %in% colnames(finClin)

subClin2022 = data.frame(subClin2022, check.names = F)
subClin2022 = subClin2022[, idx]
# converting . and - to NA 
subClin2022[subClin2022 == "."] <- NA
subClin2022[subClin2022 == "-"] <- NA

# Allele frequency, gene type and .. are good to go 

######################### variant location in exons ##########################################
subClin = readRDS("~/clinvar/sandbox/subClin2022.rds")

View(table(subClin$EXON))

table(subClin$EXON == "-")

# FALSE   TRUE 
# 153639  33854

tmp = subClin[, c(30, 40, 41)]
#tmp = tmp[!is.na(tmp$EXON),]


tmp$ratio = tmp$EXON
tmp$ratio = ifelse(tmp$ratio != "-", tmp$ratio, tmp$INTRON)
tmp = cbind(tmp, stringr::str_split_fixed(tmp$ratio, "/", 2))
names(tmp)[c(5:6)] <- c("nExon", "totalExon")
tmp$nratio = as.numeric(tmp$nExon)/as.numeric(tmp$totalExon)
##@@
tmp$nIntron = tmp$nExon
tmp$nIntron = as.factor(ifelse(tmp$INTRON == "-", tmp$nExon, "0"))
tmp$nIntron[tmp$nIntron == "-"] <- NA

# two variable need to be added to finClin
subClin2022$relativeLocationRatio = tmp$nratio
subClin2022$ExIntron = tmp$nIntron

############################### pathogenicity  score ############################################
### PUtting EVE result together
setwd("~/clinvar/variant_files/")

flist = list.files()

eveGene = gsub(".csv", "", flist)

uniprotID = read.table("~/clinvar/uniprotConversionTable.tab", fill = T, header = T)

ids = unique(ss$GENEINFO)
ids = data.frame(stringr::str_split_fixed(ids, ":", 2))
names(ids) = c("symbol", "ENTREZ_ID")

# correcting ENTREZ numbers
entz = strsplit(ids$ENTREZ_ID,split='|',fixed=TRUE)
length(entz)
entz_num = c()
for (i in 1:length(entz)){
  entz_num[i] = entz[[i]][1]
}

ids$ENTREZ_ID = as.numeric(entz_num)


# joing two IDS

geneIDs = dplyr::left_join(ids, uniprotID)
geneIDs = geneIDs[!duplicated(geneIDs$symbol),]

# to see how many of genes with NAs in uniprot in genID dataset,
# Are really protein coding:

library(biomaRt)
mart <- useDataset("hsapiens_gene_ensembl", useMart("ensembl"))
genes <- getBM(attributes= c("ensembl_gene_id","entrezgene_id","hgnc_symbol", "gene_biotype"), 
               mart= mart)

protGene = genes[genes$gene_biotype == "protein_coding",]

table(geneIDs$ENTREZ_ID[is.na(geneIDs$uniprot)] %in% protGene$entrezgene_id)

# protein codng gene not in EVE:
notInEVE = geneIDs$ENTREZ_ID[is.na(geneIDs$uniprot)][geneIDs$ENTREZ_ID[is.na(geneIDs$uniprot)] %in% protGene$entrezgene_id]

notInEVE = genes[genes$entrezgene_id %in% notInEVE,]
subClin = ss
#
ss = subClin[grep("synonymous", subClin$Consequence),]
ss$Amino_acids = paste0(ss$Amino_acids, "/", ss$Amino_acids)
ss = rbind(ss, subClin[grep("missense", subClin$Consequence)])
unik = ss$GENEINFO[!duplicated(ss$GENEINFO)]
eveList = list()

for (gene in 1:length(unik)){
  #selecting per gene
  subgen = ss[ss$GENEINFO == unik[gene],]
  # converting geneinfo to ids and merging with subgen
  ids = unique(subgen$GENEINFO)
  ids = data.frame(stringr::str_split_fixed(ids, ":", 2))
  names(ids) = c("symbol", "ENTREZ_ID")
  # correcting ENTREZ numbers
  entz = strsplit(ids$ENTREZ_ID,split='|',fixed=TRUE)
  entz_num = c()
  for (i in 1:length(entz)){
    entz_num[i] = entz[[i]][1]
  }
  ids$ENTREZ_ID = as.numeric(entz_num)
  geneIDs = dplyr::left_join(ids, uniprotID, by = "ENTREZ_ID")
  geneIDs = geneIDs[!duplicated(geneIDs$symbol),]
  subgen$uniprotID = geneIDs$uniprot
  subgen$ENTREZ_ID = geneIDs$ENTREZ_ID
  subgen$symbol = geneIDs$symbol
  #making mutation column in this dataset
  mut = data.frame(X1 = subgen$Protein_position, X2 = subgen$Amino_acids)
  mut = cbind(mut, stringr::str_split_fixed(mut$X2, "/", 2))
  names(mut)[3:4] <- c("X3", "X4")
  mut$symbol = subgen$symbol
  mut$mutation = paste0(mut$symbol,"_", mut$X3,mut$X1, mut$X4 )
  subgen$mutation  = mut$mutation
  if(!is.na(geneIDs$uniprot)){
    eveFrame = read.csv(paste0("~/clinvar/variant_files/", geneIDs$uniprot, ".csv"))
    eveFrame$uniprotID = geneIDs$uniprot
    eveFrame$symbol = geneIDs$symbol
    # substitution of one letter aa to three letter code
    eveFrame$mutation = paste0(eveFrame$symbol, "_",
                               eveFrame$wt_aa,
                               eveFrame$position,
                               eveFrame$mt_aa)
  }
  subgen = dplyr::left_join(subgen, eveFrame, by = "mutation")
  subgen = subgen[, c(3,163)]
  eveList[[gene]] = subgen
  print(paste0(gene, " out of", length(unik),", gene = ", 
               eveFrame$symbol[1], " [" , round((100*gene/length(unik)),2), "% completed]"))
}

eveFinalDataframe = do.call(rbind, eveList)


# adding eve data to subclin data
subClin2022 = dplyr::left_join(subClin2022, eveFinalDataframe, by = "ID")
#

# converting SIFT/PolyPhen to numeric score and prediction
df = cbind(t(data.frame(strsplit(subClin2022$SIFT, "\\("))), t(data.frame(strsplit(subClin2022$PolyPhen, "\\("))))
colnames(df) <- c("SIFT_Pred", "SIFT_score", "PolyPhen_Pred", "PolyPhen_score")
df = data.frame(df)
df$SIFT_score = as.numeric(sub("\\)", "", df$SIFT_score))
df$PolyPhen_score = as.numeric(sub("\\)", "", df$PolyPhen_score))
rownames(df) <- NULL

subClin2022$SIFT = df$SIFT_score
subClin2022$PolyPhen = df$PolyPhen_score

# checking the first author format
val_df = subClin2022
fc = finClin
View(table(fc$firstSubmitter))
# centers with <10 submission , encoded as "other"
# samething shold be applied for the val_df
# Before doing so, spaces have to be removed from the names

val_df$firstSubmitter = gsub(" ", "", val_df$firstSubmitter)
table(val_df$firstSubmitter %in% fc$firstSubmitter)

# recoding those centers with lower than 10 submission into other + converting centers present in the val df but not i the fc to other:
tmp = data.frame(table(val_df$firstSubmitter))
otherCent=  tmp$Var1[tmp$Freq <= 10]
val_df$firstSubmitter = ifelse(val_df$firstSubmitter %in% otherCent, "other", val_df$firstSubmitter)
# checking to see how things get improved 
table(val_df$firstSubmitter %in% fc$firstSubmitter)
#  converting centers present in the val df but not i the fc to other:
unCent = unique(val_df$firstSubmitter[-which(val_df$firstSubmitter %in% fc$firstSubmitter)])

val_df$firstSubmitter = ifelse(val_df$firstSubmitter %in% unCent, "other", val_df$firstSubmitter)

# save val_RDS
#saveRDS(val_df, "~/clinvar/sandbox/val_df.RDS")

# Consequence
val_df = val_df[which(val_df$Consequence %in% fc$Consequence),]
#
table(val_df$Consequence %in% fc$Consequence)

# SYMBOL remain unchanges however not all from val_df are in fc
val_df = val_df[val_df$SYMBOL %in% fc$SYMBOL,]
# Collection method
val_df$CollectionMethod[val_df$CollectionMethod == "case-control"] <- "CaseControl"
val_df$CollectionMethod[val_df$CollectionMethod == "clinical testing"] <- "ClinicalTesting"
val_df$CollectionMethod[val_df$CollectionMethod == "clinical testing;curation"] <- "ClinicalTesting_Curation"
val_df$CollectionMethod[val_df$CollectionMethod == "clinical testing;in vitro"] <- "ClinicalTesting_inVitro"
val_df$CollectionMethod[val_df$CollectionMethod == "clinical testing;provider interpretation"] <- "ClinicalTesting_ProviderInterpretation"
val_df$CollectionMethod[val_df$CollectionMethod == "in vitro"] <- "inVitro"
val_df$CollectionMethod[val_df$CollectionMethod == "literature only"] <- "LiteratureOnly"
val_df$CollectionMethod[val_df$CollectionMethod == "literature only;research"] <- "other"
val_df$CollectionMethod[val_df$CollectionMethod == "provider interpretation"] <- "ProviderInterpretation"
val_df$CollectionMethod[val_df$CollectionMethod == "reference population"] <- "ReferencePopulation"
val_df$CollectionMethod[val_df$CollectionMethod == "research"] <- "research"

# convet class into code
val_df$class = ifelse(val_df$class == "noConflict", "0", "1")

# Denisova
val_df$Denisova[val_df$Denisova == "C/T"] <- "other"
val_df$Denisova[val_df$Denisova == "C/G"] <- "other"
val_df$Denisova[val_df$Denisova == "A/T"] <- "other"
val_df$Denisova[val_df$Denisova == "A/G"] <- "other"
val_df$Denisova[val_df$Denisova == "A/C"] <- ""

# checking
table(val_df$Denisova)

# is the same
table(fc$Denisova)

# Allele
table(val_df$Allele %in% fc$Allele)

# defining a new column i the orginal datfraeme and thenm adding that to te final selection
tmp = subClin[, c(30,4,5)]
names(tmp)[1] <- "class"
tmp$allele = paste0(tmp$REF, ">", tmp$ALT)

tmp$REFlenCat = ifelse(nchar(tmp$REF) ==1, "1", ">1")
tmp$ALTlenCat = ifelse(nchar(tmp$ALT) ==1, "1", ">1")
# coding REF multi_nts and ALT multi_nts
tmp$allele = ifelse(tmp$REFlenCat != 1, "REFmulti_nts", 
                    ifelse(tmp$ALTlenCat != 1,"ALTmulti_nts", tmp$allele))
tmp$ID = subClin$ID
tmp = tmp[, c(7,4)]
#
val_df = dplyr::left_join(val_df, tmp, by = "ID")
val_df$Allele = val_df$allele 
val_df = val_df[, -51]
#
saveRDS(val_df, "~/clinvar/sandbox/val_df.RDS")
write.csv(val_df, "~/clinvar/val_df.csv", row.names = F)
