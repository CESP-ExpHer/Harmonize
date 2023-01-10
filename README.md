# Standardize and Harmonize

Created by: Seehyun Park

Creation date: 01 Jan 2023  

https://cesp.inserm.fr/en/equipe/exposome-and-heredity  


## A tutorial for performing harmonization of GWAS summary statistics 

## Contents
- [Notes](#notes)
- [Purpose](#purpose)
- [Modules and How to use](#modules)
  - [Standardize](#standardize)
  - [Harmonize](#harmonize)
        
## <a id="notes" /> Notes
The use of external Summary Statistics in genome-wide association study (GWAS) can significantly increase the size and diversity of the sample, increasing the power to detect association analysis. However, due to batch effects, genotyping errors, and the use of different genotyping platforms, the aggregation of multiple GWAS summary statistics can be quite challenging and difficult. Recent GWASs typically present the effects of a SNP with respect to alleles in the forward strand. However, the forward strands change from time to time as reference panels are updated, and GWASs from a few years ago cannot be guaranteed to use forward strand conventions.
If these GWAS summary statistics are not carefully quality controlled, the incorrect results might be derived when performing meta-analysis. 


## <a id="purpose" /> Purpose
The purpose of this package is to provide the data harmonization pipeline that allows users to check the quality of GWAS summary statistics before performing a meta-analysis. 

## <a id="modules" /> Modules and How to use
### <a id="standardize" /> Standardize
The Standardize module scans throught the GWAS summary statistics and standardize (***remove***) the SNPs according to the following criteria.
- **Null:** Null or NA values are genotyped at eigher beta, effect allele frequency or p-value.
- **Weird:** Genotypes are not based on combinations of 'ATCG'.
- **Duplicate:** The identical SNP is expressed in duplicate.
- **Palindromic:** The palindromic SNP, such that the alleles on the forward strand are the same as on the reverse strand (A/T on forward is T/A on the reverse). 


The Standardize module take 4 arguments as paramters
- **fileName:** The name of GWAS summary statistics file with a full pathway
- **columnInfoDict:** A python dictionary where user must enter the following information SNP, Chromosome, Position, non-Effect Allele, Effect Allele, EAF, Beta, SD, and p-value in the **VALUES** part of the columnInfoDict **{keys : values}**.
- **palindromicThreshold:** If **palindromicThreshold == None**, it implies that we ultimately remove all palindromic SNPs.  
However, we can sometimes infer which genotype is the forward strand by looking at the effect allele frequencies. The palindromic threshold allows palindromic SNPs when the EAF is lower than palindromic threshold or higher than 1 - palindromicThreshold. 
For instance, if user enters palindromicThreshold = 0.2, the EAF of palindromic SNPs below 0.2 and above 0.8 will remain in the GWAS summary statistics. 
- **verbose:** It prints the number of SNPs removed.

After creating an object of the Standardize class, user can use the **saveResult** method to save the result using the output directory taken as an input argument. 
- **outputDir:** The name of the output directory where user wants to save to. 

If so, you will have the following four **.txt** files in the output directory.
- **standardized.txt**
- **palindromicSNP.txt**
- **duplicatedSNP.txt**
- **weirdAllelesSNP.txt**

#### Standardize Code
```c
from standardizeModule import *

# Define the directory where input file is located
fileName = '/home2/users/park/python_projects/harmonize/test/GWAS_Thyr_Eur_CtrlsSafe_DTC_chrALL.txt'

# Define the column name in your own GWAS Summary Statistics
columnInfoDict = {'SNP': 'ID',
                  'Chromosome': '#CHROM',
                  'Position': 'POS',
                  'nonEffectAllele': 'ALT',
                  'EffectAllele': 'REF',
                  'EAF': 'A1_FREQ',
                  'Beta': 'BETA',
                  'SD': 'SE',
                  'p-value': 'P'}


# Users can set up the palindromic threshold value.
palindromicThreshold = 0.2

# Define the output directory where you want to save the result.
outputDir = '/home2/users/park/python_projects/harmonize/test/UKBB'

# Create your own object using `Standardize`
test = Standardize(fileName=fileName,
                   columnInfoDict=columnInfoDict,
                   palindromicThreshold=palindromicThreshold,
                   verbose=True)

# Save the result
test.saveResult(outputDir=outputDir)
```


### <a id="harmonize" /> Harmonize
The Harmonize module takes two GWAS summary statistics in order to make the effect of a SNP on the reference and the effect of that SNP on the target must correspond the the same allele. 

**!!! IMPORTANT !!!**

- You need to **Standardize** each GWAS summary statistics first before running **Harmonize**

Some examples are shown below:
#### Correct, unambigious
```c
reference beta = 0.13
reference effect allele = A
reference non-effect allele = G
reference EAF = 0.28

target beta = 0.22
target effect allele = A
target non-effect allele = G
target EAF = 0.26
```
Here the effect allele on the reference and the target is the same

#### Incorrect reference, unambigious
```c
reference beta = -0.48
reference effect allele = G
reference non-effect allele = T
reference EAF = 0.40

target beta = 0.056
target effect allele = T
target non-effect allele = G
target EAF = 0.61
```
Here the target allele is presendting the effect for the alternate allele on the reverse strand. We need to **flip** the target effect to **-0.056** to correspond to the same alleles as the reference on the forward strand, and also we need to inverser the EAF as **1-0.61 = 0.39**

#### Palindromic SNP, inferrable
It depends on the **palindromicThreshold** value. If user initailly enter **None**, there would be no palindromic SNPs appeared after standardized GWAS summary statistics. On the other hands, if **palindromicThreshold = 0.1** for instance, then we might have a below case. 

```c
reference beta = 0.20
reference effect allele = G
reference non-effect allele = C
reference EAF = 0.10

target beta = -0.046
target effect allele = G
target non-effect allele = C
target EAF = 0.91
```
Since we are giving the information about the EAF, we can infer that the target GWAS is presenting the effect on the reverse strand for the alternative allele. We need to flip the effect to **0.046** and also EAF to **1-0.91 = 0.09**.


The Harmonize module take 6 arguments as paramters
- **referenceFile:** The name of reference GWAS summary statistics file with a full pathway
- **targetFile:** The name of target GWAS summary statistics file with a full pathway
- **referenceDict:** A python dictionary where user must enter the following **reference GWAS summary statistics** information SNP, Chromosome, Position, non-Effect Allele, Effect Allele, and EAF in the **VALUES** part of the columnInfoDict **{keys : values}**.
- **targetDict:** A python dictionary where user must enter the following **target GWAS summary statistics** information SNP, Chromosome, Position, non-Effect Allele, Effect Allele, EAF, Beta, SD, and p-value in the **VALUES** part of the columnInfoDict **{keys : values}**.
- **intersection:** True, if user only want to harmonize the overlapped SNPs (**intersection**) 
- **verbose:** It prints the number of SNPs unharmonized.



After creating an object of the Harmonize class, user can use the **saveResult** method to save the result using the output directory taken as an input argument. 
- **outputDir:** The name of the output directory where user wants to save to. 

If so, you will have the following four **.txt** files in the output directory.
- **harmonized.txt** 
- **diffGenotyped.txt** (SNPs that are completely differently genotyped between reference and target file)
- **changed.txt** (SNPs with altered beta and EAF)
- **removed.txt** (SNPs removed due to intersection)

#### Harmonize Code
```c
# Define the full path where reference file is located.
referenceFile = '/home2/users/park/python_projects/harmonize/test/EPITHYR/standardized.txt'

# Define tge full path where target file is located.
targetFile = '/home2/users/park/python_projects/harmonize/test/UKBB/standardized.txt'

# Define the column name in your own GWAS Summary Statistics
refInfoDict = {'SNP': 'SNP',
               'Chromosome': 'chr',
               'Position': 'pos',
               'nonEffectAllele': 'A1',
               'EffectAllele': 'A2',
               'EAF': 'freq2'}   

# User needs to modify the 'value' part of the targetInfoDict
targetInfoDict = {'SNP': 'ID',
                  'Chromosome': '#CHROM',
                  'Position': 'POS',
                  'nonEffectAllele': 'ALT',
                  'EffectAllele': 'REF',
                  'EAF': 'A1_FREQ',
                  'Beta': 'BETA',
                  'SD': 'SE',
                  'p-value': 'P'}

# True if user wants to include only overlapped SNPs (intersection)
intersection = True

# Define the output directory where you want to save the result.
outputDir = '/home2/users/park/python_projects/harmonize/test/EPITHYR_UKBB'

# Create your own object using `Harmoniza`
test = Harmonize(referenceFile=referenceFile,
                 targetFile=targetFile,
                 referenceDict=refInfoDict,
                 targetDict=targetInfoDict,
                 intersection=intersection,
                 verbose=True)

# In case you want to save the result
test.saveResult(outputDir=outputDir)

```


We are now ready to do a meta-analysis

--The END--
