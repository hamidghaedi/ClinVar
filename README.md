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
