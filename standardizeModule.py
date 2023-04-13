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


class InvalidAlleleException(Exception):
    pass


def getLastExceptionDescription():
    return str(sys.exc_info()[0]) + " / " + str(sys.exc_info()[1])


def createComplementDict():
    """
    :return: the dictionary with the complement of nucleotides
    """
    complement = defaultdict(lambda: str)
    complement['A'] = 'T'
    complement['T'] = 'A'
    complement['C'] = 'G'
    complement['G'] = 'C'
    return complement


def findIndex(inputList, infoDict):
    """
    :return: the dictionary that contains the index position of inputList from the infoDict
    """
    positionDict = dict()
    for key, value in infoDict.items():
        positionDict[key] = inputList.index(value)
    return positionDict


def findAlternativeAllele(allele, verbose=False):
    '''
    :param allele: the string format of the allele
    :param verbose: shows the invalid allele if 'verbose = True'
    :return: return the alternative allele corresponding to the 'AlleleComplement'
    '''

    alleleComplement = createComplementDict()

    try:
        allele = allele.upper()
        if set(allele.upper()) <= set('ATCG'):
            alternativeAllele = ''
            if len(allele) != 1:
                for i in range(0, len(allele)):
                    try:
                        alternativeAllele = alternativeAllele + alleleComplement[allele[i]]
                    except TypeError:
                        pass
            else:
                alternativeAllele = alleleComplement[allele]

            return alternativeAllele

        else:
            raise InvalidAlleleException

    except InvalidAlleleException:
        if verbose == True:
            print('This is invalid allele ==> ' + allele)
        else:
            pass


def isPalindromic(EA, nEA):
    """
    :param EA: Effect Allele
    :param nEA: non Effect Allele
    :return: True if palindromic
    """
    return nEA == findAlternativeAllele(EA)


class Standardize:
    def __init__(self, fileName, columnInfoDict, palindromicThreshold, verbose=False):
        """
        :param fileName: The full pathway of the input file
        :param columnInfoDict: The dictionary which contains the proper name of each column
        :param palindromicThreshold: The value of palindromic Threshold
        :param verbose: Print the number of the duplicated chromosoem:position
        """
        self.fileName = fileName
        self.columnInfoDict = columnInfoDict
        self.palindromicThreshold = palindromicThreshold
        self.verbose = verbose


    def checkChrPosDuplicated(self, chrName, posName):
        """
        :param chrName: Chromosome column name
        :param posName: Position column name
        :return: Return the list of the duplicated chromosome:position format
        """
        try:
            with open(self.fileName, 'r') as inputFile:
                # Read from the second line after the header
                header = inputFile.readline().strip().split('\t')
                chrPosition = header.index(chrName)
                posPosition = header.index(posName)
                lines = inputFile.readlines()

                setDuplicated = set()
                setChrPos = set()
                for line in lines:
                    line = line.split('\t')
                    chr = str(line[chrPosition])
                    pos = str(line[posPosition])
                    ChrPos = chr + ":" + pos

                    # Check if duplicated SNPs are existed in duplicatedSNP set
                    if ChrPos in setChrPos:
                        setDuplicated.add(ChrPos)

                    setChrPos.add(ChrPos)

                return setDuplicated

        except FileNotFoundError:
            print("Unable to open " + self.fileName + ":" + getLastExceptionDescription())
        except:
            print("Critical exception occured, aborting " + getLastExceptionDescription())
            raise
        finally:
            inputFile.close()

    def standardizeData(self):
        try:
            with open(self.fileName, 'r') as inputFile:
                header = inputFile.readline().strip().split('\t')
                indexPosition = findIndex(inputList=header, infoDict=self.columnInfoDict)

                # Read from the second line after the header
                lines = inputFile.readlines()

                # Keeping track of the different categories of SNPs:
                nb_SNP_standardized = 0  # 1-SNPs which standardized,
                nb_SNP_palindromic = 0  # 2-SNPs which are palindromic,
                nb_SNP_duplicated = 0  # 3-SNPs which are duplicated,
                nb_SNP_weird_alleles = 0  # 4-SNPs with an allele that are not composed of ATGC,

                ChrPosDuplicated = self.checkChrPosDuplicated(chrName=self.columnInfoDict['Chromosome'], posName=self.columnInfoDict['Position'])

                # make result dictionary
                resultDict = defaultdict(list)
                # Define column name for each file
                resultDict['standardized'].append([self.columnInfoDict['SNP'],
                                                   self.columnInfoDict['Chromosome'],
                                                   self.columnInfoDict['Position'],
                                                   self.columnInfoDict['nonEffectAllele'],
                                                   self.columnInfoDict['EffectAllele'],
                                                   self.columnInfoDict['EAF'],
                                                   self.columnInfoDict['Beta'],
                                                   self.columnInfoDict['SD'],
                                                   self.columnInfoDict['p-value'],
                                                   'MAF'])
                resultDict['palindromic'].append([self.columnInfoDict['SNP'],
                                                  self.columnInfoDict['Chromosome'],
                                                  self.columnInfoDict['Position'],
                                                  self.columnInfoDict['nonEffectAllele'],
                                                  self.columnInfoDict['EffectAllele'],
                                                  self.columnInfoDict['EAF'],
                                                  self.columnInfoDict['Beta'],
                                                  self.columnInfoDict['SD'],
                                                  self.columnInfoDict['p-value'],
                                                  'MAF'])
                resultDict['weird'].append([self.columnInfoDict['SNP'],
                                            self.columnInfoDict['Chromosome'],
                                            self.columnInfoDict['Position'],
                                            self.columnInfoDict['nonEffectAllele'],
                                            self.columnInfoDict['EffectAllele']])
                resultDict['duplicated'].append([self.columnInfoDict['SNP'],
                                                 self.columnInfoDict['Chromosome'],
                                                 self.columnInfoDict['Position'],
                                                 self.columnInfoDict['nonEffectAllele'],
                                                 self.columnInfoDict['EffectAllele']])

                for line in lines:
                    line = line.strip().split('\t')

                    # Define the SNP information
                    SNP = line[indexPosition['SNP']]
                    chr = line[indexPosition['Chromosome']]
                    pos = line[indexPosition['Position']]
                    ChrPos = chr + ":" + pos
                    nonEffectAllele = line[indexPosition['nonEffectAllele']]  # non effect allele
                    EffectAllele = line[indexPosition['EffectAllele']]  # effect allele
                    EAF = float(line[indexPosition['EAF']])  # effect allele frequency
                    nEAF = 1 - float(EAF)
                    beta = float(line[indexPosition['Beta']])
                    sd = float(line[indexPosition['SD']])
                    pValue = line[indexPosition['p-value']]

                    # Define alternative allele
                    altEffectAllele = findAlternativeAllele(EffectAllele)
                    altNonEffectAllele = findAlternativeAllele(nonEffectAllele)

                    # Define minor allele frequency
                    if float(EAF) < 0.5:
                        MAF = EAF
                    else:
                        MAF = 1 - float(EAF)

                    # if (eaf or beta or sd or P_value) == NULL or NA, it removes the SNP (means it skips the line)
                    if EAF in {'NULL', 'NA'} or beta in {'NULL', 'NA'} or pValue in {'NULL','NA'} or beta == 0 or sd == 0:
                        continue

                    if altEffectAllele is None or altNonEffectAllele is None:
                        nb_SNP_weird_alleles += 1
                        weirdList = [SNP, chr, pos, nonEffectAllele, EffectAllele]
                        resultDict['weird'].append(weirdList)
                        continue

                    # if the SNP is an ambiguous SNPs (Palindromic), it removes the SNP (means it skips the line) and prints the SNP in its output file
                    if len(EffectAllele) == len(nonEffectAllele):
                        if isPalindromic(EA=EffectAllele, nEA=nonEffectAllele):
                            if self.palindromicThreshold:
                                if EAF >= self.palindromicThreshold and EAF <= 1 - self.palindromicThreshold:
                                    nb_SNP_palindromic += 1
                                    palindromicList = [SNP, chr, pos, nonEffectAllele, EffectAllele, str(EAF),
                                                       str(beta), str(sd), str(pValue), str(MAF)]
                                    resultDict['palindromic'].append(palindromicList)
                                    continue
                            else:
                                nb_SNP_palindromic += 1
                                palindromicList = [SNP, chr, pos, nonEffectAllele, EffectAllele, str(EAF), str(beta),
                                                   str(sd), str(pValue), str(MAF)]
                                resultDict['palindromic'].append(palindromicList)
                                continue

                    if ChrPos in ChrPosDuplicated:
                        nb_SNP_duplicated += 1
                        duplicatedList = [SNP, chr, pos, nonEffectAllele, EffectAllele]
                        resultDict['duplicated'].append(duplicatedList)
                        continue

                    standardizedList = [SNP, chr, pos, nonEffectAllele, EffectAllele, str(EAF), str(beta), str(sd),
                                      str(pValue), str(MAF)]
                    nb_SNP_standardized += 1

                    # Save into resultDict
                    resultDict['standardized'].append(standardizedList)

                if self.verbose:
                    print("The number of palindromic SNP is ", nb_SNP_palindromic)
                    print("The number of duplicated SNP is ", nb_SNP_duplicated)
                    print("The number of weird SNP is ", nb_SNP_weird_alleles)

                return resultDict

        except FileNotFoundError:
            print("Unable to open " + self.fileName + ":" + getLastExceptionDescription())
        except:
            print("Critical exception occured, aborting " + getLastExceptionDescription())
            raise
        finally:
            inputFile.close()

    def saveResult(self, outputDir = None):
        """
        :param outputDir: The path of result output directory
        """

        isExist = os.path.exists(outputDir)
        if not isExist:
            os.makedirs(outputDir)
            print("The output directory is created! ")
        else:
            print("The output directory is already existed")

        name = os.path.basename(self.fileName)
        name_without_ext = name[:name.rindex(".")]

        standardized = open(f'{outputDir}/{name_without_ext}_standardized.txt', 'w')
        palindromic = open(f'{outputDir}/{name_without_ext}_palindromicSNP.txt', 'w')
        duplicated = open(f'{outputDir}/{name_without_ext}_duplicatedSNP.txt', 'w')
        weird = open(f'{outputDir}/{name_without_ext}_weirdAllelesSNP.txt', 'w')

        resultDict = self.standardizeData()
        print("\n********* Writing & Saving results  *********\n")
        print('\n'.join(map('\t'.join, resultDict['standardized'])), file=standardized)
        print('\n'.join(map('\t'.join, resultDict['palindromic'])), file=palindromic)
        print('\n'.join(map('\t'.join, resultDict['duplicated'])), file=duplicated)
        print('\n'.join(map('\t'.join, resultDict['weird'])), file=weird)

        standardized.close()
        palindromic.close()
        duplicated.close()
        weird.close()
        print("********* DONE *********")


if __name__ == "__main__":

    # Define the directory where input file is located
    fileName = '/home2/users/park/python_projects/harmonize/test/20230406_8STUDY/all_Epithyr_EPIC_UKBB_deCODE.txt'

    # Define the column name in your own GWAS Summary Statistics
    # User needs to modify the 'value' part of the columnInfoDict
    columnInfoDict = {'SNP': 'SNP',
                      'Chromosome': 'chr',
                      'Position': 'pos',
                      'nonEffectAllele': 'OA',
                      'EffectAllele': 'EA',
                      'EAF': 'EAF',
                      'Beta': 'Beta',
                      'SD': 'SE',
                      'p-value': 'P'}

    # Define the output directory where you want to save the result.
    outputDir = '/home2/users/park/python_projects/harmonize/test/20230406_8STUDY/all_Epithyr_EPIC_UKBB_deCODE_0.3'

    # Users can set up the palindromic threshold value.
    # Allow palindromic allele when their EAF is below than palindromic Threshold
    # For instance, if EAF is 0.29, we allow this allele even though palindromic.
    palindromicThreshold = 0.3

    # Create your own object using `Standardize`
    test = Standardize(fileName=fileName,
                       columnInfoDict=columnInfoDict,
                       palindromicThreshold=palindromicThreshold,
                       verbose=True)

    # In case you want to save the result
    test.saveResult(outputDir=outputDir)



