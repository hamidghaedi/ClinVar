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
<img src="https://raw.githubusercontent.com/hamidghaedi/clinvar/main/figs/class_dist.PNG" width="800" height="400"> |

Receiving more submissions from the community, variants in the class "criteria provided, single submitter" will turn into one of the "criteria provided, multiple submitter, no conflicts" or "criteria provided, conflicting interpretations" classes. To find out which feature  of variants may contibute to this conversion, I am doing this project.


--------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
#### file download and preprocessing
ClinVar variant file were download from ClinVar FTP under the name "clinvar_20210619.vcf.gz". In order to convert format , file were converted to a tab delimited file by
```bash
#ComputeCanada instance
module load StdEnv/2020 vcflib
vcf2tsv clinvar_20210619.vcf.gz > clinvar_tsv.txt
```
To add annotation to the file online version of VEP was used. To return one line per variant, filter setting was set to "show most severe consequence per variant".
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
#### Analysis: variant consequnce annotation
#### Pathogenicity prediction scores
#### Conservation scores
