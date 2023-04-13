# 					CODE FOR REFORMATTING THE GWAS DATA
# ================================================================================
# Summary: The code reads the GWAS Summary Statistics data
# and writes it in a proper format that could be used for our future analyses
# ================================================================================
# Written first by: Elise Lucotte
# Modified by: Yazdan Asgari and Seehyun PARK
# Initial Creation Date: 12/2020
# Edited Date: 1/2023
# https://cesp.inserm.fr/en/equipe/exposome-and-heredity
# ================================================================================

import os
import sys
from collections import defaultdict
from standardizeModule import *

class Harmonize():
    def __init__(self, referenceFile, targetFile, referenceDict, targetDict, intersection, verbose,sampleSize):
        """
        :param referenceFile: The full path of your standardized reference file.
        :param targetFile: The full path of your standardized target file which you want to harmonize with reference
        :param referenceDict: The dictionary which contains the information of column name of reference file
        :param targetDict: The dictionary which contains the information of column name of target file
        :param intersection: to decide wheter to only include intersection part of both GWAS summary statistics or not
        :param smapleSize: the sample size number
        """

        self.referenceFile = referenceFile
        self.targetFile = targetFile
        self.referenceDict = referenceDict
        self.targetDict = targetDict
        self.intersection = intersection
        self.verbose = verbose
        self.sampleSize = sampleSize

    def getChrPosEAFInfo(self):
        ChrPosEAFDict = dict()
        try:
            with open(self.referenceFile, 'r') as inputFile:
                header = inputFile.readline().strip().split('\t')
                indexPosition = findIndex(inputList=header, infoDict=self.referenceDict)

                # Read from the second line after the header
                lines = inputFile.readlines()
                for line in lines:
                    line = line.strip().split('\t')
                    try:
                        ChrPosEAFDict[line[indexPosition['Chromosome']]][line[indexPosition['Position']]] = \
                            [line[indexPosition['SNP']], line[indexPosition['nonEffectAllele']],
                             line[indexPosition['EffectAllele']], line[indexPosition['EAF']]]
                    except KeyError:
                        ChrPosEAFDict[line[indexPosition['Chromosome']]] = {}
                        ChrPosEAFDict[line[indexPosition['Chromosome']]][line[indexPosition['Position']]] = \
                            [line[indexPosition['SNP']], line[indexPosition['nonEffectAllele']],
                             line[indexPosition['EffectAllele']], line[indexPosition['EAF']]]

                return ChrPosEAFDict
        except FileNotFoundError:
            print("Unable to open " + self.referenceFile + ":" + getLastExceptionDescription())
        except:
            print("Critical exception occured, aborting " + getLastExceptionDescription())
            raise
        finally:
            inputFile.close()

    def harmonizeData(self):
        referenceDictInfo = self.getChrPosEAFInfo()
        nb_diffGenotyped = 0
        nb_removed = 0
        nb_changed = 0

        # make result dictionary
        resultDict = defaultdict(list)
        # Define column name for each file
        resultDict['harmonized'].append([self.targetDict['SNP'],
                                         self.targetDict['Chromosome'],
                                         self.targetDict['Position'],
                                         self.targetDict['nonEffectAllele'],
                                         self.targetDict['EffectAllele'],
                                         self.targetDict['EAF'],
                                         self.targetDict['Beta'],
                                         self.targetDict['SD'],
                                         self.targetDict['p-value'],
                                         'MAF',
                                         'N',
                                         'chr_pos'])
        resultDict['diffGenotyped'].append([self.targetDict['SNP'],
                                            self.targetDict['Chromosome'],
                                            self.targetDict['Position'],
                                            'ref_nonEffectAllele',
                                            'ref_EffectAllele',
                                            self.targetDict['nonEffectAllele'],
                                            self.targetDict['EffectAllele']])

        resultDict['removed'].append([self.targetDict['SNP'],
                                      self.targetDict['Chromosome'],
                                      self.targetDict['Position'],
                                      self.targetDict['nonEffectAllele'],
                                      self.targetDict['EffectAllele'],
                                      self.targetDict['EAF'],
                                      self.targetDict['Beta'],
                                      self.targetDict['SD'],
                                      self.targetDict['p-value'],
                                      'MAF'])

        resultDict['changedBeta'].append([self.targetDict['SNP'],
                                          self.targetDict['Chromosome'],
                                          self.targetDict['Position'],
                                          'ref_nonEffectAllele',
                                          'ref_EffectAllele',
                                          self.targetDict['nonEffectAllele'],
                                          self.targetDict['EffectAllele']])

        standardizedTarget = open(self.targetFile, 'r').readlines()
        header = standardizedTarget[0].strip().split('\t')
        indexPosition = findIndex(inputList=header, infoDict=self.targetDict)

        for individual in standardizedTarget[1:]:
            individual = individual.strip().split('\t')

            #Define the target infomation
            target_SNP = individual[indexPosition['SNP']]
            target_chr = individual[indexPosition['Chromosome']]
            target_pos = individual[indexPosition['Position']]
            target_nonEffectAllele = individual[indexPosition['nonEffectAllele']]
            target_EffectAllele = individual[indexPosition['EffectAllele']]
            target_EAF = individual[indexPosition['EAF']]
            target_beta = individual[indexPosition['Beta']]
            target_sd = individual[indexPosition['SD']]
            target_pValue = individual[indexPosition['p-value']]
            target_chr_pos = target_chr + ":" + target_pos

            # Define the alternative allele
            alt_target_nonEffectAllele = findAlternativeAllele(target_nonEffectAllele)
            alt_target_EffectAllele = findAlternativeAllele(target_EffectAllele)

            # Define the minor allele frequency
            if float(target_EAF) < 0.5:
                MAF = float(target_EAF)
            else:
                MAF = 1 - float(target_EAF)

            try:
                ref_nonEffectAllele = referenceDictInfo[target_chr][target_pos][1]
                ref_EffectAllele = referenceDictInfo[target_chr][target_pos][2]
                ref_EAF = referenceDictInfo[target_chr][target_pos][3]

                # Check if both GWAS are not palindromic
                if not (isPalindromic(EA=target_EffectAllele, nEA=target_nonEffectAllele) and
                        isPalindromic(EA=ref_EffectAllele, nEA=ref_nonEffectAllele)):

                    # Leave Beta and EAF as it is
                    if (target_nonEffectAllele == ref_nonEffectAllele and target_EffectAllele == ref_EffectAllele) or \
                            (alt_target_nonEffectAllele == ref_nonEffectAllele and alt_target_EffectAllele == ref_EffectAllele):
                        harmonizedList = [target_SNP, target_chr, target_pos, target_nonEffectAllele,target_EffectAllele,
                                           str(target_EAF), str(target_beta), str(target_sd), str(target_pValue), str(MAF),
                                          str(self.sampleSize), target_chr_pos]
                        resultDict['harmonized'].append(harmonizedList)

                    # Change the sign of beta & EAF to 1-EAF
                    elif (target_nonEffectAllele == ref_EffectAllele and target_EffectAllele == ref_nonEffectAllele) or \
                            (alt_target_nonEffectAllele == ref_EffectAllele and alt_target_EffectAllele == ref_nonEffectAllele):
                        nb_changed += 1
                        harmonizedList = [target_SNP, target_chr, target_pos, target_EffectAllele,target_nonEffectAllele,
                                           str(1 - float(target_EAF)), str(-float(target_beta)), str(target_sd),str(target_pValue), str(MAF),
                                          str(self.sampleSize), target_chr_pos]
                        changedList = [target_SNP, target_chr, target_pos, ref_nonEffectAllele, ref_EffectAllele, target_nonEffectAllele, target_EffectAllele]
                        resultDict['harmonized'].append(harmonizedList)
                        resultDict['changedBeta'].append(changedList)

                    else:
                        diffGenotypedList = [target_SNP, target_chr, target_pos, ref_nonEffectAllele, ref_EffectAllele, target_nonEffectAllele, target_EffectAllele]
                        nb_diffGenotyped += 1
                        resultDict['diffGenotyped'].append(diffGenotypedList)

                # if one GWAS is palindromic and others not, vice versa
                elif ((isPalindromic(EA=ref_EffectAllele, nEA=ref_nonEffectAllele) and not
                isPalindromic(EA=target_EffectAllele, nEA=target_nonEffectAllele)) or \
                      (isPalindromic(EA=target_EffectAllele, nEA=target_nonEffectAllele) and not
                      isPalindromic(EA=ref_EffectAllele, nEA=ref_nonEffectAllele))):
                    diffGenotypedList = [target_SNP, target_chr, target_pos,ref_nonEffectAllele,
                                         ref_EffectAllele, target_nonEffectAllele, target_EffectAllele]
                    nb_diffGenotyped += 1
                    resultDict['diffGenotyped'].append(diffGenotypedList)

                # if both GWAS are palindromic
                else:
                    # Both are same palindromic allele
                    if (target_nonEffectAllele == ref_nonEffectAllele and target_EffectAllele == ref_EffectAllele) or \
                            (alt_target_nonEffectAllele == ref_nonEffectAllele and alt_target_EffectAllele == ref_EffectAllele):

                        # The defaulf value of the maximum difference between ref_EAF and target_EAF is 0.6 ;
                        # => In other words, comparison between at least 0.3 vs 0.7 is the miminum
                        if abs(float(ref_EAF) - float(target_EAF)) >= 0.4:
                            nb_changed += 1
                            harmonizedList = [target_SNP, target_chr, target_pos, target_EffectAllele, target_nonEffectAllele,
                                               str(1 - float(target_EAF)), str(-float(target_beta)), str(target_sd),
                                               str(target_pValue), str(MAF),
                                              str(self.sampleSize), target_chr_pos]
                            changedList = [target_SNP, target_chr, target_pos, ref_nonEffectAllele,
                                           ref_EffectAllele, target_nonEffectAllele, target_EffectAllele]
                            resultDict['harmonized'].append(harmonizedList)
                            resultDict['changedBeta'].append(changedList)

                        else:
                            harmonizedList = [target_SNP, target_chr, target_pos, target_nonEffectAllele, target_EffectAllele,
                                               str(target_EAF), str(target_beta), str(target_sd), str(target_pValue),str(MAF),
                                              str(self.sampleSize), target_chr_pos]
                            resultDict['harmonized'].append(harmonizedList)

                    # Both are palindromic, but different genotyped (discard)
                    else:
                        diffGenotypedList = [target_SNP, target_chr, target_pos, ref_nonEffectAllele,
                                             ref_EffectAllele,target_nonEffectAllele, target_EffectAllele]
                        nb_diffGenotyped += 1
                        resultDict['diffGenotyped'].append(diffGenotypedList)

            except KeyError:
                if self.intersection == True:
                    nb_removed += 1
                    removedList = [target_SNP, target_chr, target_pos, target_nonEffectAllele,target_EffectAllele,
                                           str(target_EAF), str(target_beta), str(target_sd), str(target_pValue), str(MAF)]
                    resultDict['removed'].append(removedList)
                else:
                    harmonizedList = [target_SNP, target_chr, target_pos, target_nonEffectAllele, target_EffectAllele,
                                       str(target_EAF), str(target_beta), str(target_sd), str(target_pValue), str(MAF),
                                      str(self.sampleSize), target_chr_pos]
                    resultDict['harmonized'].append(harmonizedList)

        if self.verbose:
            print("The number of differently genotyped SNP is ", nb_diffGenotyped)
            print("The number of removed due to intersection SNP is ", nb_removed)
            print("The number of changed beta SNP is ", nb_changed)

        return resultDict

    def saveResult(self, outputDir):
        """
        :param outputDir: The path of result output directory
        """
        isExist = os.path.exists(outputDir)
        if not isExist:
            os.makedirs(outputDir)
            print("The output directory is created! ")
        else:
            print("The output directory is already existed")

        name = os.path.basename(self.targetFile)
        name_without_ext = name[:name.rindex(".")]

        harmonized = open(f'{outputDir}/{name_without_ext}_harmonized.txt', 'w')
        diffGenotyped = open(f'{outputDir}/{name_without_ext}_diffGenotyped.txt', 'w')
        removed = open(f'{outputDir}/{name_without_ext}_removed.txt', 'w')
        changed = open(f'{outputDir}/{name_without_ext}_changed.txt', 'w')

        resultDict = self.harmonizeData()
        print("\n********* Writing & Saving results  *********\n")
        print('\n'.join(map('\t'.join, resultDict['harmonized'])), file=harmonized)
        print('\n'.join(map('\t'.join, resultDict['diffGenotyped'])), file=diffGenotyped)
        print('\n'.join(map('\t'.join, resultDict['removed'])), file=removed)
        print('\n'.join(map('\t'.join, resultDict['changedBeta'])), file=changed)

        harmonized.close()
        diffGenotyped.close()
        removed.close()
        changed.close()
        print("********* DONE *********")




if __name__ == "__main__":

    # Define the full path where reference file is located.
    referenceFile = '/home2/users/park/python_projects/harmonize/test/20230406_8STUDY/all_Epithyr_EPIC_UKBB_deCODE_0.3/all_Epithyr_EPIC_UKBB_deCODE_standardized.txt'

    # Define tge full path where target file is located.
    targetFile = '/home2/users/park/python_projects/harmonize/test/20230406_8STUDY/SummStat_deCODE_5GWAS.txt'

    # Define the column name in your own GWAS Summary Statistics
    # User needs to modify the 'value' part of the refInfoDict
    refInfoDict = {'SNP': 'SNP',
                   'Chromosome': 'chr',
                   'Position': 'pos',
                   'nonEffectAllele': 'OA',
                   'EffectAllele': 'EA',
                   'EAF': 'EAF'}

    # User needs to modify the 'value' part of the targetInfoDict
    targetInfoDict = {'SNP': 'SNP',
                      'Chromosome': 'chr',
                      'Position': 'pos',
                      'nonEffectAllele': 'OA',
                      'EffectAllele': 'EA',
                      'EAF': 'EAF',
                      'Beta': 'Beta',
                      'SD': 'SE',
                      'p-value': 'P'}

    # True if user wants to include only overlapped SNPs (intersection)
    intersection = False

    # Sample Size
    sampleSize = 290551

    # Define the output directory where you want to save the result.
    outputDir = '/home2/users/park/python_projects/harmonize/test/20230406_8STUDY/deCODE_5GWAS'

    # Create your own object using `Harmonize`
    test = Harmonize(referenceFile=referenceFile,
                     targetFile=targetFile,
                     referenceDict=refInfoDict,
                     targetDict=targetInfoDict,
                     intersection=intersection,
                     verbose=True,
                     sampleSize=sampleSize)

    # In case you want to save the result
    test.saveResult(outputDir=outputDir)




